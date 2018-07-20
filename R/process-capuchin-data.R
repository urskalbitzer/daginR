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
#' @param group Name of the social group
#' @param start,end Define first and last day (YYYY-MM-DD) of the desired time period
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
#' @param focaldata Dataset with focaldata from capuchin monkeys
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
    effort_temp <- data.frame(YearOf = i,
                              IndA = females,
                              IndB = females,
                              Effort_AB = 0) %>%
      complete(IndA, IndB, fill = list(Effort_AB = 0, YearOf = i)) %>%
      mutate_at(c("IndA", "IndB"), as.character)

    # 2. For each female create a new dataframe with one line per focal
    for(j in females){
      focals_IndA <- focaldata %>%
        filter(NameOf == j) %>%
        distinct(FocalID, .keep_all = TRUE) %>%
        filter(AgeAtFocal >= minAgeF) %>%
        mutate(YearOf = year(FocalBegin), DateOfFocal = as.Date(FocalBegin)) %>%
        filter(YearOf == i) %>%
        select(FocalID, NameOf, YearOf, DateOfFocal, FocalBegin, FocalDurationCorrected)

      # and check for all other females (potential focal partners), which focals
      # with IndA were conducted while IndB was around and summarise
      # focalduration
      for(k in females){
        temp <- focals_IndA %>%
          filter(DateOfFocal %within% female_presence[female_presence$NameOf==k,]$female_pres_int) %>%
          summarise(Effort_AB = sum(FocalDurationCorrected))

        effort_temp[effort_temp$IndA == j & effort_temp$IndB == k,]$Effort_AB = temp$Effort_AB
      }
    }

    # Create table effort_AB if doesn't exist, or add effort_temp to effort_AB
    # if exists
    if(!exists("final_effort_table")){
      final_effort_table <- effort_temp
    } else final_effort_table <- bind_rows(final_effort_table, effort_temp)

  }

  effort_A_B <- final_effort_table %>%
    mutate(GroupCode = group) %>%
    select(GroupCode, YearOf, IndA, IndB, Effort_AB)

  return (effort_A_B)
}



#' Calculation of behavior rates for female dyads
#'
#' This function calculates rates of different behaviors for all female dyads
#' meeting specific criteria (age, focal effort) within a specified group and
#' time period (only works with one group and one time period at a time).
#'
#' @param behaviours_duration Behaviours with duration (e.g. grooming) to be
#'   included
#' @param beahviours_events Behaviours without duration (e.g. approach) to be
#'   included
#' @param minEffortDyad Minimum duration (in hours) a dyad has to be observed to
#'   be included (observation time of IndA while IndB was around + observation
#'   time IndB while IndA was around)
#' @param minEffortInd The minimum duration (in hours) an individual has to be
#'   observed to be included
#' @param roles_included Which roles should be included, e.g., (D)irected,
#'   (R)eceived, or (M)utual
#' @param dyads_summarized Should dyadic values be summarized? If FALSE, each
#'   behavior will be listed as IndA->IndB and as IndB->IndA separately (per
#'   dyad and per role)
#' @inheritParams determine_presence
#' @inheritParams calc_effort_AB
#'
#' @return If values are summarized, the rates IndA->IndB and IndB->IndA are
#'   averaged. This is slightly different from calculating behaviour counts (or
#'   durations) IndA->IndB + IndB->IndA divided by EffortA + EffortB (unless
#'   EffortA = EffortB). Thus, rates are not weighted by focal effort of IndA
#'   and IndB
#'
#' @export
#'
#' @examples

calc_dyadic_rates <- function (focaldata, census_pace, biography_pace,
                               group, start, end, minAgeF,
                               behaviours_duration, behaviours_events,
                               minEffortDyad = 0, minEffortInd = 0,
                               roles_included = c("D", "R", "M"),
                               dyads_summarized = TRUE) {

  # Check if min effort is defined
  if(is.na(minEffortDyad) | minEffortDyad < 0) stop("minEffortDyad NA or negative")
  if(is.na(minEffortInd) | minEffortInd < 0) stop("minEffortInd NA or negative")

  # Function only works with one group and one time period at a time, so check if there's only one of each
  if(is.na(group) | length(group) > 1) stop("Group not defined or more than one group")
  if(length(start)>1) stop("More than one start date")
  if(length(end)>1) stop("More than one end date")

  # Only use focaldata from females that 1) have minAgeF, belong to group and within study-period
  focaldata_temp <- focaldata %>%
    filter(AgeAtFocal > minAgeF) %>%
    filter (GroupCode == group) %>%
    filter (FocalBegin >= as.Date (start) & FocalBegin <= as.Date (end)) %>%
    rename(IndA = NameOf, IndB = InteractNameOf)

  # Get focalhours per individual during study-period
  effort_ind <- focaldata_temp %>%
    distinct (FocalID, .keep_all = TRUE) %>%
    group_by (IndA) %>%
    summarise (Effort_Ind = sum (FocalDurationCorrected)/60)

  # Get the time IndA was observed while IndB was around and join with individual observation time
  effort_A_while_B <- calc_effort_AB(focaldata, census_pace, biography_pace, group, start, end, minAgeF) %>%
    mutate(Effort_AB = Effort_AB/60)

  # Determine all individuals within the (focal)-dataset
  individuals <- unique(effort_A_while_B$IndA)

  #### Sum up duration/counts of behaviours for each individual and partner during study-period
  # Behaviours with duration
  b_d_dyadic <- focaldata_temp %>%
    filter(BehavName %in% behaviours_duration & Role %in% roles_included) %>%
    filter(IndB %in% individuals) %>%
    group_by(GroupCode, IndA, IndB, BehavName, Role) %>%
    summarise(Absolute = sum (BehavDuration)/60)

  # Event behaviours
  b_e_dyadic <- focaldata_temp %>%
    filter(BehavName %in% behaviours_events & Role %in% roles_included) %>%
    filter(IndB %in% individuals) %>%
    group_by(GroupCode, IndA, IndB, BehavName, Role) %>%
    summarise (Absolute = n())

  #### Combine tables with duration- and event-behaviour
  # Determine all occuring roles
  roles <- unique(rbind (b_d_dyadic, b_e_dyadic)$Role)

  # Check which of the variables (individuals, roles, behaviours_durations, behaviours_events) is the longest
  # to create a table with minimum length to create a "complete" dataframe with all combinations of ind, roles, behav
  # i.e. create dataframe with all possible combination of dyads, behaviours and roles
  maxlength <- max(length(individuals), length(roles), (length(behaviours_duration) + length(behaviours_events)))
  alldyads_template <- data.frame(GroupCode = group, IndA = rep(individuals, len = maxlength),
                                  IndB = rep(individuals, len = maxlength),
                                  BehavName = rep(c(behaviours_duration, behaviours_events), len = maxlength),
                                  Role = rep(roles, len = maxlength), Absolute = 0) %>%
    complete(IndA, IndB, BehavName, Role, fill = list(GroupCode = group, Absolute = 0)) %>%
    mutate_at(vars(IndA, IndB, BehavName, Role, GroupCode), funs(as.character(.)))

  # In order to get a table with all behavior for all dyads, bind tables with duration- and event-behaviours.
  # Then remove these cases from the alldyads-table (where "Absolute" is always 0)
  # and then bind tables with duration- and event-behaviour to the reduced alldyad-table (to have a full table again)
  b_dyadic_1 <- rbind(b_e_dyadic, b_d_dyadic) %>%
    anti_join(alldyads_template, ., by = c("IndA", "IndB", "BehavName", "Role")) %>%
    bind_rows(., b_e_dyadic, b_d_dyadic)

  # Then spread roles to calculate absolute values for all roles combined
  # Different roles are not further needed below but left in query in
  # case that only specific directions would be included into the function
  b_dyadic_2 <- b_dyadic_1 %>%
    spread (Role, Absolute, fill = 0) %>%
    mutate(AllRoles = rowSums(.[roles])) %>%
    gather_ ("Role", "Absolute", c(roles, "AllRoles")) %>%
    mutate(Role = if_else(Role == "AllRoles", paste(roles, collapse = "_"), Role)) %>%
    filter(!is.na(IndA) & !is.na(IndB))

  # Calculate rates for each behavior and dyad
  # This requires the effort for each IndA while IndB was present to correct for coresidence
  dyads_bothdirections <- b_dyadic_2 %>%
    left_join(., effort_A_while_B, by = c("IndA", "IndB", "GroupCode")) %>%
    filter(Effort_AB > 0) %>%
    mutate (Rate = Absolute/Effort_AB)

  # Also add the total focal time for both IndA and IndB and replace NA by 0
  dyads_bothdirections_2 <- dyads_bothdirections %>%
    left_join(., select(effort_ind, IndA, IndA_Effort = Effort_Ind), by = "IndA") %>%
    left_join(., select(effort_ind, IndB = IndA, IndB_Effort = Effort_Ind), by = "IndB") %>%
    replace_na(list(IndA_Effort = 0, IndB_Effort = 0))

  # Removing "self-dyads" and create a column "Dyad" in which the Ind higher up in the alphabet is always the first Ind
  output_1 <- dyads_bothdirections_2 %>%
    filter(IndA != IndB) %>%
    mutate(Dyad = ifelse (as.character(IndA) < as.character(IndB),
                          paste (IndA, IndB, sep = "_"),
                          paste (IndB, IndA, sep = "_")))

  # Check which dyads do not meet minEffortDyad and minEffortInd criteria, print a message, and filter out.
  # 1. Check
  dyads_to_filter_out <- output_1 %>%
    distinct(IndA, IndB, .keep_all = T) %>%
    group_by(Dyad) %>%
    mutate(Effort_sum_AB_BA = sum(Effort_AB)) %>%
    filter(IndA_Effort < minEffortInd | IndB_Effort < minEffortInd | Effort_sum_AB_BA < minEffortDyad)
  # 2. Message
  if(nrow(dyads_to_filter_out) > 0) print(paste("The following dyad has been sorted out (not enough data):", dyads_to_filter_out$Dyad))

  # 3. Filter
  output_2 <- filter(output_1, !(Dyad %in% dyads_to_filter_out$Dyad))

  # If "dyads_summarized = TRUE", summarize dyads so that each dyad is only listed once with an averaged rate per behavior/dyad/role
  # (Otherwise each dyad is listed twice with each individual once as IndA and once as IndB)
  if(dyads_summarized){
    output_3 <- output_2 %>%
      group_by(Dyad, BehavName, Role) %>%
      mutate(effort_sum_AB_BA = sum(Effort_AB),
             MeanRate_Dyad = mean (Rate)) %>%
      ungroup %>%
      select(Dyad, BehavName, GroupCode, Role, YearOf, MeanRate_Dyad) %>%
      # mutate(MeanRate_Dyad = round(MeanRate_Dyad, 3)) %>%
      distinct(Dyad, BehavName, Role, .keep_all = TRUE)
  }else{
    output_3 <- output_2 %>%
      select(Dyad, IndA, IndB, BehavName, GroupCode, Year, Role, Absolute,
             Effort_AB, Rate)
    # mutate(Absolute = round(Absolute, 5),
    #        Effort_AB = round(Effort_AB, 2),
    #        Rate = round(Rate, 6))
  }

  output_4 <- output_3 %>%
    mutate(StartDate = start, EndDate = end)

  return(output_4)

}
