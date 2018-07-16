#' Determine presence of female in group
#'
#' This function uses census and biography data from the PACE database to
#' determine the periods of time females with specified minimum age
#' \code{minAgeF} were present in a \code{group} from \code{start} to
#' \code{end}.
#' The function does only work for females with known DOB
#'
#' @param census_pace Dataframe from the PACE database with monthly capuchin census data
#' @param biography_pace Dataframe from the PACE database with individual biography data
#' @param group The name of the social group for which the presence of females
#'   will be determined
#' @param start,end Define first and last day (YYYY-MM-DD) of the time period
#'   for which the presence of females will be determined
#' @param minAgeF The minimum age of females to be included (in years)
#'
#' @export
#'
#' @examples
#' female_presence(census_pace, biogaphy_pace, group = "GUAN", start =
#' "2007-01-01", end = "2011-12-31", minAgeF = 5)

determine_presence <- function (census_pace, biography_pace, group, start, end, minAgeF) {

  # Only include census data from females in specified group between the start and end month
  censusdata <- census_pace %>%
    filter(GroupCode == group & Sex == "F" &
             year(CensusDateOf) >= year(start) & month(CensusDateOf) >= month(start) &
             year(CensusDateOf) <= year(end) & month(CensusDateOf) <= month(end)) %>%
    arrange(CensusDateOf)

  # Check if last month has census data
  if(month(last(censusdata$CensusDateOf)) != month(end)) stop ("Last month without census data")

  # Determine females that fullfill age conditions for at least at one point during the analyses interval
  females <- unique((censusdata %>%
                       filter(difftime(CensusDateOf, DateOfBirth, units = "days")/365.25 >= minAgeF))$NameOf)

  # Create census table only including these females
  temp_census <- censusdata %>%
    filter(NameOf %in% females) %>%
    select(GroupCode, CensusDateOf, CensusYear, CensusMonth, NameOf, DateOfBirth, StatusCodeLong) %>%
    mutate(age_study_start = as.numeric(difftime(start, DateOfBirth, units = "days")/365.25)) %>%
    arrange(NameOf, CensusDateOf)

  # If a female was younger than "minAgeF" at "start" --> presence start date
  # (start_int) is her birthday day that she reached "minAgeF", otherwise
  # start_int = start
  temp_f_start <- temp_census %>%
    group_by(NameOf) %>%
    mutate(start_int = if_else(age_study_start < minAgeF, DateOfBirth + minAgeF*365.25, as.Date(start)))

  # Get biography data for all females to determine depart dates
  temp_biography <- biography_pace %>%
    filter (NameOf %in% females) %>%
    select(NameOf, DepartDate)

  temp_f_end <- temp_f_start %>%
    left_join(temp_biography, by = "NameOf") %>%
    group_by(NameOf) %>%
    mutate(last_census = last(CensusDateOf), last_status = last(StatusCodeLong)) %>%
    # Set end_int to "end" if last census on female was done in last month of the study period and female still alive
    # Set end_int to DepartDate if last status was "Missing" or "Dead" and within 3 months (±45 days) of DepartDate
    # (to prevent cases where females was just missing for that month but generally still around?)
    # otherwise "error
    mutate(end_int = case_when((month(last_census) == month(end) & last_status == "Alive") ~ as.character(end),
                               (last_status %in% c("Missing", "Dead") &
                                  DepartDate %within% interval(as.Date(last_census)-45, as.Date(last_census)+45)) ~ as.character(DepartDate),
           TRUE ~ "error"))

  if(any(temp_f_end$end_int == "error")) {stop(paste("Problem with last census and depart date in", group))} else{
    female_presence <- temp_f_end %>%
      ungroup() %>%
      distinct(NameOf, .keep_all = TRUE) %>%
      mutate_at(vars(start_int, end_int), funs(as.Date))  %>%
      select(GroupCode, NameOf, DateOfBirth, start_int, end_int) %>%
      mutate(female_pres_int = interval(paste(start_int, "00:00:00", sep = " "), paste(end_int, "00:00:00", sep = " "))) %>%
      rename(female_pres_start = start_int, female_pres_end = end_int) %>%
      mutate(assessed_int_start = start, assessed_int_end = end, min_age = minAgeF)
  }
  return (female_presence)
}

#' Calculate focal time for female A while female B was present
#'
#' Calculates focal effort with female A while female B was present for the
#' specifiec group and time period
#'
#' @param focaldata Dataset with focal data from capuchin monkeys
#' @inheritParams determine_presence
#'
#' @export
#'
#' @examples
#' effort_AB(focaldata = dataset, census_pace = Vcensus, biography_pace =
#' Vbiography, group = "GUAN", start = "2007-01-01", end = "2011-12-31", minAgeF
#' = 5)


calc_effort_AB <- function (focaldata, census_pace, biography_pace, group, start, end, minAgeF) {

  female_presence <- determine_presence(census_pace, biography_pace, group, start, end, minAgeF)
  females <- unique(female_presence$NameOf)

  # 1. For each year
  years <- seq(from = year(start), to = year(end), by = 1)

  for(i in years){
    # Create a dataframe with all possible dyads and set dyadic focaleffort to 0
    effort_temp <- data.frame(year_of = i,
                              ind_A = females,
                              ind_B = females,
                              effort_A_while_B_present = 0) %>%
      complete(ind_A, ind_B, fill = list(effort_A_while_B_present = 0, year_of = i)) %>%
      mutate_at(c("ind_A", "ind_B"), as.character)

    # 2. For each female create a new dataframe with one line per focal
    for(j in females){
      focals_ind_A <- focaldata %>%
        filter(name_of == j) %>%
        distinct(FocalID, .keep_all = TRUE) %>%
        filter(age_at_focal >= minAgeF) %>%
        mutate(year_of = year(FocalBegin), date_of_focal = as.Date(FocalBegin)) %>%
        filter(year_of == i) %>%
        select(FocalID, name_of, year_of, date_of_focal, foal_begin = FocalBegin, focal_duration_corrected = FocalDurationCorrected)

      # and check for all other females (potential focal partners), which focals
      # with ind_A were conducted while ind_B was around and summarise
      # focalduration
      for(k in females){
        temp <- focals_ind_A %>%
          filter(date_of_focal %within% female_presence[female_presence$name_of==k,]$female_pres_int) %>%
          summarise(effort_A_while_B_present = sum(focal_duration_corrected))

        effort_temp[effort_temp$ind_A == j & effort_temp$ind_B == k,]$effort_A_while_B_present = temp$effort_A_while_B_present
      }}

    # Create table effort_AB if doesn't exist, or add effort_temp to effort_AB
    # if exists
    if(!exists("effort_A_B")){
      final_effort_table <- effort_temp
    } else final_effort_table <- bind_rows(final_effort_table, effort_temp)

  }

  effort_A_B <- final_effort_table %>%
    mutate(GroupCode = group) %>%
    select(GroupCode, year_of, ind_A, ind_B, effort_A_while_B_present)

  return (effort_A_B)
}
