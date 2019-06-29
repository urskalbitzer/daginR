#' Determine presence of adult individuals in the Small Camp Red Colobus group
#' for the specified period of time
#'
#' This function uses scan data in combination with data on the
#' first-adult-month of individuals to determine for each month whether an
#' individual was 1) adult and 2) present
#'
#' Default is that individuals are considered present if recorded for at least
#' \code{n_min} (default 5) times per month for 3 months in a row
#'
#'
#' @param scandata_df Data frame with red colobus scan data
#' @param first_month_adult_df Dataframe with first-month-adult dates
#' @param from_date,to_date Start and end date of analyzed period of time.
#'   Format should be "yyyy-mm-dd"
#' @param n_min Minimum number of times an individual has to be recorded
#'   (Individual + NN) to be considered present. If recorded < n_min for 3 month
#'   in row, considered as absent (including the first month considered)
#' @param details Provide details for each individual and each month (default FALSE)
#'
#' @export
#'
#' @examples
#'
#'

determine_presence_rc <- function(scandata_df, first_month_adult_df,
                         from_date, to_date, n_min = 5, details = FALSE){
  ind_per_month <- scandata_df %>%
    group_by(YearOf, MonthOf, Individual) %>%
    summarize(ind_count = n()) %>%
    ungroup()

  nn_per_month <- scandata_df %>%
    group_by(YearOf, MonthOf, NN) %>%
    summarize(nn_count = n()) %>%
    ungroup()

  name_per_month <- ind_per_month %>%
    left_join(., nn_per_month,
              by = c("YearOf", "MonthOf", "Individual" = "NN")) %>%
    mutate_at(vars(contains("count")), ~replace_na(., 0)) %>%
    mutate(count_total = ind_count + nn_count,
           Year_Month = paste0(YearOf, "-", str_pad(MonthOf, 2, pad = 0))) %>%
    select(-ind_count, -nn_count, -YearOf, -MonthOf) %>%
    filter(!is.na(Individual)) %>%
    complete(Individual, Year_Month,
             fill = list(count_total = 0)) %>%
    arrange(Individual, Year_Month)

  adults_per_month <- name_per_month %>%
    left_join(select(first_month_adult_df, Individual = Name, first_month_adult = Year_month),
              by = "Individual") %>%
    filter(!is.na(first_month_adult)) %>%
    mutate(Date_adult = as.Date(paste0(first_month_adult, "-01")),
           Date_count = as.Date(paste0(Year_Month, "-01")),
           Adult = if_else(Date_count >= Date_adult, 1, 0)) %>%
    select(-Date_adult, -Date_count, -first_month_adult)

  presence_month <- adults_per_month %>%
    group_by(Individual) %>%
    mutate(Present = if_else((count_total < n_min &
                                lead(count_total, n = 1, default = n_min) < n_min &
                                lead(count_total, n = 2, default = n_min) < n_min) |
                               (count_total < n_min &
                                  lag(count_total, n = 1, default = n_min) < n_min &
                                  lag(count_total, n = 2, default = n_min) < n_min) |
                               (count_total < n_min &
                                  lead(count_total, n = 1, default = n_min) < n_min &
                                  lag(count_total, n = 1, default = n_min) < n_min),
                             0, 1)) %>%
    ungroup %>%
    mutate(Present_Adult = if_else(Present == 1 & Adult == 1, 1, 0))

  if(!details){
    presence_month <- presence_month %>%
      select(-count_total, -Adult, -Present)
  }
  return(presence_month)
}


#' Red Colobus: Create table to link all synonyms to an individual name, and
#' also includes the sex of individuals
#'
#' This function uses the IndividualsRC table to create a longer table with
#' names, link names (synonyms) and sex of individuals. It also included unknown
#' individuals (e.g. Adult Female) with the code number as synonym, and the sex
#' (if known)
#'
#' @param individual_table IndividualsRC table should include all individuals
#'   with all possible alternative names for each individual
#'
#' @export
#'
#' @examples
#'
#'

create_link_table <- function(individual_table){

  # Preparatory steps with the individual table (should be corrected in the original copy?)
  id_table <- individual_table %>%
    filter(!(Name_New %in% c("FPF")))

  # Create a "link-table" for known individuals:
  # One column with official names, containing all names from the column Name_New
  # One column with link-names, containing all names from Name_New PLUS all alternative and old names

  link_ids_1 <- id_table %>%
    select(Name = Name_New) %>% mutate(Link_Name = Name)

  link_ids_2 <-  id_table %>%
    filter(!is.na(Name_Old)) %>% select(Name = Name_New, Link_Name = Name_Old)

  link_ids_3 <-  id_table %>%
    filter(!is.na(Alternative_Name1)) %>% select(Name = Name_New, Link_Name = Alternative_Name1)

  link_ids_4 <-  id_table %>%
    filter(!is.na(Alternative_Name2)) %>% select(Name = Name_New, Link_Name = Alternative_Name2)

  link_ids_5 <-  id_table %>%
    filter(!is.na(Alternative_Name3)) %>% select(Name = Name_New, Link_Name = Alternative_Name3)

  link_ids <- rbind(link_ids_1, link_ids_2, link_ids_3, link_ids_4, link_ids_5) %>%
    arrange(Name) %>%
    left_join(., select(individual_table, Name_New, Sex), by = c("Name" = "Name_New")) %>%
    mutate(Sex = if_else(Sex == "M", "Male", if_else(Sex == "F", "Female", NA_character_)))

  # Add codes for unknown individuals
  link_unknown <- tibble(Name = c("Adult Male", "Adult Female", "Adult Female with Infant",
                                  "Adult", "Subadult Male", "Subadult Female", "Subadult",
                                  "Juvenile", "Infant", "Unknown", "Juvenile Male",
                                  "Juvenile Female", "Infant Male", "Infant Female",
                                  "Other Species"),
                         Link_Name = as.character(c(1:14, 20)),
                         Sex = c("Male", "Female", "Female", "Unknown", "Male",
                                 "Female", "Unknown", "Unknown", "Unknown", "Unknown",
                                 "Male", "Female", "Male", "Female",
                                 "Unknown"))

  link_ids <- rbind(link_ids, link_unknown)

  return(link_ids)
}


#' Red Colobus: Create age-sex codes based on the link table
#'
#' @param link_table Table creates with \code{create_link_table}
#'
#' @export
#'
#' @examples
#'
#'

create_age_sex_codes <- function(link_table){
  age_sex_codes <- link_table %>%
    filter(Link_Name %in% as.character(c(1:14, 20))) %>%
    mutate(Age = case_when(str_detect(Name, "Adult") ~ "Adult",
                           str_detect(Name, "Subadult") ~ "Subadult",
                           str_detect(Name, "Juvenile") ~ "Juvenile",
                           str_detect(Name, "Infant") ~ "Infant",
                           TRUE ~ "Unknown"),
           Age_Sex_Code = as.numeric(Link_Name)) %>%
    select(Description = Name, Age_Sex_Code, Age, Sex)

  return(age_sex_codes)
}


#' Count number of scans for a period with Individual while it is present and adult
#' Then count for each other Individual number of times it was Nearest Neighbour (NN) or could have
#' been NN while adult and present in group
#'
#'
#' @param scandata_df Data frame with red colobus scan data
#' @param first_month_adult_df Dataframe with first-month-adult dates
#' @param from_date,to_date Start and end date of analyzed period of time.
#'   Format should be "yyyy-mm-dd". Includes both first and last day.
#' @param id_sex_df Table with sex of individuals. This is the link table creates with `create_link_table`
#' @param max_nn_dist The maximum distance for which an individual is still considered NN. Default is 5
#'
#' @export
#'
#' @examples
#'
#'

get_dyadic_association_rc <- function(scan_data_df, first_month_adult_df, id_sex_df,
                                      from_date, to_date, max_nn_dist = 5){
  # TO DO Add argument checks

  require(tidyverse)
  require(lubridate)
  require(daginR)

  # Limit scan_data_df to time period
  scan_temp <- scandata %>%
    filter(DateTime >= as.Date(from_date) &
             DateTime <= as.Date(to_date) &
             is.na(Usable_Record))

  # Only keep required columns
  scan_temp <- scan_temp %>%
    select(DateTime, Scan_ID, Individual, Ind_Age_Sex, Activity, NN, NN_Age_Sex, Dist)

  # Create additional date columns columns
  scan_temp$YearOf <- year(scan_temp$DateTime)
  scan_temp$MonthOf <- month(scan_temp$DateTime)
  scan_temp$DayOf <- day(scan_temp$DateTime)

  # Create presence table using the `determine_presence_rc` function
  # This tables indicates for each individual and month during the study period whether this individual was present and adult
  presence <- determine_presence_rc(scandata_df = scan_temp,
                                    first_month_adult_df,
                                    from_date,
                                    to_date,
                                    n_min = 5, details = FALSE)
  # Now, the actual calculations:
  # For all individuals, for each month, count lines with:
  #   1. A as Individual while A and B in the group as adults
  #   2. A as Individual with B as NN

  # Monthly counts of Ind and Ind_NN
  ind_count <- scan_temp %>%
    mutate(Year_Month = paste0(YearOf, "-", str_pad(MonthOf, 2, pad = 0))) %>%
    group_by(Year_Month, Individual) %>%
    summarize(Count_Ind = n()) %>%
    ungroup

  ind_nn_count<- scan_temp %>%
    filter(Dist <= max_nn_dist) %>%
    mutate(Year_Month = paste0(YearOf, "-", str_pad(MonthOf, 2, pad = 0))) %>%
    group_by(Year_Month, Individual, NN) %>%
    summarize(Count_Ind_NN = n()) %>%
    ungroup

  # Setup dataframe with all dyads (both directions) per month and information on adult-presence
  Individual <-  unique(presence$Individual)
  NN <-  Individual
  Year_Month = unique(presence$Year_Month)

  dyad_summaries <- expand_grid(Individual, NN, Year_Month) %>%
    filter(Individual != NN) %>%
    left_join(select(presence, Year_Month, Individual, Ind_Pres_Adult = Present_Adult),
              by = c("Individual", "Year_Month")) %>%
    left_join(select(presence, Year_Month, NN = Individual, NN_Pres_Adult = Present_Adult),
              by = c("Year_Month", "NN"))

  # Combine with counts to get the
  # 1) total counts per individual per year_month and
  # 2) counts per individual with NN per year_month

  dyad_summaries <- dyad_summaries %>%
    left_join(ind_count, by = c("Year_Month", "Individual")) %>%
    left_join(ind_nn_count, by = c("Year_Month", "Individual", "NN"))

  # Now summarize for the entire period
  # Filter month were both Individual and NN were present and adult
  # Set NA values to 0 (because both were present as adult but never observed alone/together)
  # Summarize for each Individual-NN combination for entire period

  dyad_summaries <- dyad_summaries %>%
    filter(Ind_Pres_Adult == 1 & NN_Pres_Adult) %>%
    replace_na(list(Count_Ind = 0, Count_Ind_NN = 0)) %>%
    group_by(Individual, NN) %>%
    summarize(Count_Ind_NN_Present = sum(Count_Ind), Count_Ind_NN = sum(Count_Ind_NN))

  # Add information on sex of Individual and NN
  dyad_summaries <- dyad_summaries %>%
    left_join(select(id_sex_df, Individual = Name, Ind_Sex = Sex), by = "Individual") %>%
    left_join(select(id_sex_df, NN = Name, NN_Sex = Sex), by = "NN")

  dyad_summaries <- dyad_summaries %>%
    ungroup

  return(dyad_summaries)

}
