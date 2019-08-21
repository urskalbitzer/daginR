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

#' Determine presence of adult individuals in the Small Camp Red Colobus group
#' for the specified period of time
#'
#' This function uses scan data in combination with data on the
#' first-adult-month of individuals to determine for each month whether an
#' individual was 1) adult and 2) present
#'
#' Note that estimating the presence for short periods of time can be diffent
#' from estimates from longer periods of time. Information what happens before
#' or after the specified period of time is unknown and the function assumes
#' that preceding and following month have a count of n_min. Ideally include 2
#' months before and 2 months after each observation period. (Note to myself:
#' this could also be included in the function.)
#'
#' Default is that individuals are considered absent if recorded for less than
#' \code{n_min} (default 5) times per month for 3 months in a row
#'
#'
#' @param scandata_df Data frame with red colobus scan data
#' @param first_month_adult_df Dataframe with first-month-adult dates
#' @param from_date,to_date Start and end date of analyzed period of time.
#'   Format should be "yyyy-mm-dd". First and last date will be included!
#' @param n_min Minimum number of times an individual has to be recorded
#'   (Individual + NN) to be considered present. If recorded < n_min for 3 month
#'   in row, considered as absent (including the first month considered)
#' @param details Provide details for each individual and each month (default
#'   FALSE)
#'
#' @export
#'
#' @examples
#'
#'

get_presence_rc <- function(scandata_df, first_month_adult_df,
                            from_date, to_date, n_min = 5, details = FALSE){

  # Limit scandata table to required columns, and derive year and month column
  # Use as.Date(DateTime) because otherwise from_data and to_date will be considered as datetime
  df <- scandata_df %>%
    filter(!is.na(Individual),
           as.Date(DateTime) >= as.Date(from_date) &
             as.Date(DateTime) <= as.Date(to_date)) %>%
    select(DateTime, Individual, NN) %>%
    mutate(YearOf = lubridate::year(DateTime), MonthOf = lubridate::month(DateTime))

  # Count the number of times each name was recorded as Individual per month
  ind_per_month <- df %>%
    group_by(YearOf, MonthOf, Individual) %>%
    summarize(ind_count = n()) %>%
    ungroup()

  # Count the number of times each name was recorded as NN per month
  nn_per_month <- df %>%
    group_by(YearOf, MonthOf, NN) %>%
    summarize(nn_count = n()) %>%
    ungroup()

  # Combines tables to get count how often each name was recorded in total (Individual + NN)
  name_per_month <- ind_per_month %>%
    left_join(., nn_per_month,
              by = c("YearOf", "MonthOf", "Individual" = "NN")) %>%
    mutate_at(vars(contains("count")), ~replace_na(., 0)) %>%
    mutate(count_total = ind_count + nn_count) %>%
    select(YearOf, MonthOf, Name = Individual, count_total) %>%
    complete(Name, nesting(YearOf, MonthOf),
             fill = list(count_total = 0)) %>%
    arrange(Name, YearOf, MonthOf)

  # Determine for each name in each month whether it was adult
  adults_per_month <- name_per_month %>%
    left_join(select(first_month_adult_df, Name, first_month_adult = Year_month),
              by = "Name") %>%
    filter(!is.na(first_month_adult)) %>%
    mutate(Date_adult = as.Date(paste0(first_month_adult, "-01")),
           Date_count = as.Date(paste0(YearOf, "-", MonthOf, "-01")),
           Adult = if_else(Date_count >= Date_adult, 1, 0)) %>%
    select(-Date_adult, -Date_count, -first_month_adult)

  # Now, determine the presence
  # An individual is considered as absent if observed less than n_min times in 3 consecutive months:
  # 1. Observed less than n_min times in current and two following months
  # 2. Observed less than n_min times in current and two past months
  # 3. Observed less than n_min times in past, current, and following month
  presence_month <- adults_per_month %>%
    group_by(Name) %>%
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
    mutate(Present_As_Adult = if_else(Present == 1 & Adult == 1, 1, 0))

  if(!details){
    presence_month <- presence_month %>%
      select(-count_total, -Adult, -Present)
  }
  return(presence_month)
}

#' Prepare Red Colobus scandata for determination of association index:
#' Trim the data frame from start to end, and only keep relevant columns
#'
#' @param presence_df Dataframe created with `get_presence_rc`
#' @param id_sex_df Table with sex of individuals. This is the link table
#'   creates with `create_link_table`
#' @param filter_adults Only keep lines where Individual and NN were adult?
#' @param keep_activity_interactant If FALSE, re-sort columns and sort out activity, interactant, and inter_age_sex
#' @inheritParams get_presence_rc
#'
#' @export
#'
#' @examples
#'
#'

prep_rc_scandata_for_associations <- function(scandata_df,
                                              from_date, to_date,
                                              presence_df, id_sex_df,
                                              filter_adults = TRUE,
                                              keep_activity_interactant = FALSE){
  # Only keep data within time period and remove non-usable records and rows without Individual

  df <- scandata_df %>%
    filter(DateTime >= as.Date(from_date) &
             DateTime <= as.Date(to_date) &
             !is.na(Individual) &
             is.na(Usable_Record))

  # Only keep required columns

  df <- df %>%
    select(DateTime, Scan_ID, Individual, Ind_Age_Sex, NN, NN_Age_Sex, Dist,
           Activity, Interactant, Inter_Age_Sex)

  # Create additional date columns with year, month, and day
  df <- df %>%
    mutate(YearOf = lubridate::year(DateTime),
           MonthOf = lubridate::month(DateTime),
           DayOf = lubridate::day(DateTime))

  # Combine with presence_df for 1) Individual and 2) NN
  # Replace NA by 0 because NA means Name was not adult in that month
  # (only individuals that were ever present as adults are listed in presence table)

  df <- df %>%
    left_join(rename(presence_df, Ind_Adult = Present_As_Adult),
              by = c("Individual" = "Name", "YearOf", "MonthOf")) %>%
    left_join(rename(presence_df, NN_Adult = Present_As_Adult),
              by = c("NN" = "Name", "YearOf", "MonthOf")) %>%
    mutate_at(vars(Ind_Adult, NN_Adult), ~replace_na(., 0))

  # Combine with id_sex_df to get sex for Individuals and NNs from individual table
  df <- df %>%
    left_join(rename(id_sex_df, Ind_Sex = Sex), by = c("Individual" = "Name")) %>%
    left_join(rename(id_sex_df, NN_Sex = Sex), by = c("NN" = "Name"))

  # Re-sort columns and sort out activity, interactant, and inter_age_sex if keep_activity_interactant = FALSE
  if(keep_activity_interactant){
    df <- df %>%
      select(Scan_ID, DateTime, YearOf, MonthOf, DayOf,
             Individual, Ind_Sex, Ind_Adult,
             Activity, Interactant, Inter_Age_Sex,
             NN, NN_Sex, NN_Adult, Dist)

  } else {
    df <- df %>%
      select(Scan_ID, DateTime, YearOf, MonthOf, DayOf,
             Individual, Ind_Sex, Ind_Adult,
             NN, NN_Sex, NN_Adult, Dist)
  }

  # If filter_adult is TRUE, only keep lines where Individual and NN were adult
  # Then also sort out columns Ind_Adult and NN_Adult, because they all have 1
  if(filter_adults){
    df <- df %>%
      filter(Ind_Adult == 1 & NN_Adult == 1)
  }

  return(df)
}


#' For each dyad, this function counts the number of times and Individual A was
#' recorded as Individual while Individual A and B were both adult and present
#' (`Ind_NN_Present`). Then, it calculates how often Individual B was recorded
#' as the Nearest Neighbour (NN) of Individual A for this time period (where
#' both were adult and present; `Ind_NN_together`)
#'
#' Furthermore, this function does pre-network data permutations `n_permutation`
#' times (0 for no permutations) using the function `permute_scan_data`. For
#' each permutation, `n_iter_per_permutation` NNs are swapped.
#'
#' The results are provided as 3-dimensional array:
#' Dimension 1 and 2 are the different individuals (i.e. dyad matrix)
#' Dimension 3 is data. i = 1: Ind_NN_Present; i = 2: Ind_NN_Together;
#' 3 - n_permutations: Permutations of Ind_NN_Together
#'
#' Limitations: This function cannot (yet) permute only within
#' sex classes. To achieve that, subset scandata and only include females or
#' males before running the function. Set `ind_sex_permutation` and
#' `nn_sex_permutation` to "Male"/"Female" respectively (although this should be
#' working with default arguments)
#'
#' @param scandata_df Data frame with red colobus scan data
#' @param id_sex_df Table with sex of individuals. This is the link table
#'   creates with `create_link_table`
#' @param max_nn_dist The maximum distance for which an individual is still
#'   considered NN. Default is 5 (m)
#' @param n_permutations How many permutations should be conducted?
#' @inheritParams get_presence_rc
#' @inheritParams prep_rc_scandata_for_associations
#' @inheritParams permute_scan_data
#'
#' @export
#'
#' @examples
#'

get_dyadic_association_rc <- function(scandata_df,
                                      id_sex_df,
                                      presence_df,
                                      max_nn_dist = 5,
                                      n_permutations = 0,
                                      ind_sex_permutation = c("Male", "Female"),
                                      nn_sex_permutation = c("Male", "Female"),
                                      n_iter_per_permutation = 1, ...){
  # TO DO:
  # 1. Include argument checks (including checks of provided data)

  # Setup dataframe with all dyads (both directions) per month and information on adult-presence and sex
  Names <- union(unique(scandata_df$Individual), unique(scandata_df$NN))
  YearOf = distinct(scandata_df, YearOf, MonthOf)$YearOf
  MonthOf = distinct(scandata_df, YearOf, MonthOf)$MonthOf
  dyad_table <- expand_grid(IndA = Names, IndB = Names,
                            nesting(YearOf, MonthOf)) %>%
    filter(IndA != IndB) %>%
    left_join(rename(presence_df, IndA_Adult = Present_As_Adult),
              by = c("IndA" = "Name", "YearOf", "MonthOf")) %>%
    left_join(rename(presence_df, IndB_Adult = Present_As_Adult),
              by = c("IndB" = "Name", "YearOf", "MonthOf"))

  # Setup array for entire study period
  # Dimension 1 and 2 are individuals
  # Dimension 3 is Ind_NN_Present, Ind_NN_Together, and all permutations of Ind_NN_Together
  dyad_array <- array(, dim = c(length(Names),
                                length(Names),
                                n_permutations + 2))
  d1_2_names <- sort(unique(Names))
  d3_names <- c("Ind_NN_Present",
                paste0("Ind_NN_Together_perm_",
                       stringr::str_pad(0:n_permutations,
                                        width = nchar(n_permutations),
                                        side = "left", pad = "0")))
  dimnames(dyad_array) <- list(d1_2_names, d1_2_names, d3_names)

  rm(list = c("Names", "YearOf", "MonthOf", "d1_2_names", "d3_names"))

  # Setup temporary scan data table (will be modified in permutation below)
  temp_scan_df <- scandata_df %>%
    select(YearOf, MonthOf, Individual, Ind_Adult, Ind_Sex,
           NN, NN_Adult, NN_Sex, Dist)

  # Calculate for all individuals, for each month lines where
  # A is Individual while A and B both in the group as adults

  # Monthly counts of Ind
  ind_count <- temp_scan_df %>%
    group_by(YearOf, MonthOf, Individual) %>%
    summarize(Count_Ind = n()) %>% ungroup

  # Combine with dyad summary - one value per direction (IndA-IndB and IndB-IndA)
  d_sum <- dyad_table %>%
    left_join(ind_count, by = c("YearOf", "MonthOf", "IndA" = "Individual"))

  # Only keep lines in which both indidivuals (IndA, IndB) were adult and present
  # Then, filter out NAs (which means both were around but IndA not observed. Equal to 0)

  d_sum <- d_sum %>%
    filter(IndA_Adult == 1 & IndB_Adult == 1) %>%
    filter(!is.na(Count_Ind))

  # Now, summarize values for entire period
  d_sum <- d_sum %>%
    group_by(IndA, IndB) %>%
    summarize(Count_Ind_NN_Present = sum(Count_Ind))

  # Create temporary matrix, fill with Count_Ind_NN_Present values
  # Then add to first layer of dyadic array ("Ind_NN_Present")
  temp_matrix <- as.matrix(dyad_array[,,"Ind_NN_Present"])
  temp_matrix[as.matrix(d_sum[c("IndA", "IndB")])] <-     d_sum$Count_Ind_NN_Present
  dyad_array[,, "Ind_NN_Present"] <- temp_matrix

  ## Permute scan data, for each permutation calculate number of times IndA
  ## observed with IndB as NN

  for(i in 0:n_permutations){
    cat("\rPermutation ", i, "/", n_permutations)
    # First calculation is on original data frame
    if(i == 0) temp_scan_df <- temp_scan_df
    # For each following permutation, use the temp dataset
    # Permute n_iter_per_permutation times
    if(i > 0){
      temp_scan_df <- permute_scan_data(ind_sex_permutation = ind_sex_permutation,
                                        nn_sex_permutation = nn_sex_permutation,
                                        scandata_df = temp_scan_df,
                                        presence_df = presence_df,
                                        n_iter_per_permutation = n_iter_per_permutation, ...)
    }

    # Calculate for all individuals, for each month lines where
    # A was Individual with B as NN

    # Monthly counts of Ind per NN
    ind_nn_count <- temp_scan_df %>%
      filter(Dist <= max_nn_dist) %>%
      group_by(YearOf, MonthOf, Individual, NN) %>%
      summarize(Count_Ind_NN_Together = n()) %>% ungroup

    # Combine with dyad summary - one value per direction (IndA-IndB and IndB-IndA)
    d_sum <- dyad_table %>%
      left_join(ind_nn_count, by = c("YearOf", "MonthOf", "IndA" = "Individual", "IndB" = "NN"))

    # Only keep lines in which both individuals were adult and present
    # Then, filter out NAs (which means both were around but IndB never NN of IndA. Equal to 0)
    d_sum <- d_sum %>%
      filter(IndA_Adult == 1 & IndB_Adult == 1) %>%
      filter(!is.na(Count_Ind_NN_Together))

    # Now, summarize values for entire period
    d_sum <- d_sum %>%
      group_by(IndA, IndB) %>%
      summarize(Count_Ind_NN_Together = sum(Count_Ind_NN_Together))

    # Fill dyadic arrays with Count_Ind_NN_Together values. Use (i + 2) because i
    # starts with 0 (original dataframe) but 1 is Count_Ind_NN_Present
    temp_matrix <- as.matrix(dyad_array[,,(i+2)]) # Use empty layer
    temp_matrix[as.matrix(d_sum[c("IndA", "IndB")])] <- d_sum$Count_Ind_NN_Together
    dyad_array[,, (i+2)] <- temp_matrix
  }

  # And return the array
  return(dyad_array)
}


#' Permutation of Nearest Neighbours in Red Colobus scan data: Following
#' suggestion by Farine (2017, MEE), this function permutes raw observations
#' ('pre-network) data. This method has the advantage that "account for
#' underlying structure in the generated social network, and thus can reduce
#' both type I and type II error rates." (Farine, 2017)
#'
#' It does so by choosing a random line from the scan data, and then checks if
#' Individual (1) and NN (1) are both adult and have the sex according to arguments
#' `ind_sex_permutation` and `nn_sex_permutation`.
#'
#' Then, a second (random) row is chosen, and the same checks are done for
#' Individual (2) and NN (2). Furthermore, it is controlled that NN1 was present at
#' date of second row, and NN2 present at date of first row.
#'
#' If these conditions are met, the two NNs are swapped.
#'
#' This procedure is repeated `n_iter_per_permutation` times.
#'
#'
#' @param ind_sex_permutation default `c("Male", "Female")`
#' @param nn_sex_permutation default `c("Male", "Female")`
#' @param n_iter_per_permutation Number of iterations ('NN swaps') per
#'   permutation
#' @param show_i Show progress>
#' @param set.seed If TRUE, use `set.seed(1209)`
#' @inheritParams get_presence_rc
#' @inheritParams prep_rc_scandata_for_associations
#'
#' @export
#'
#' @examples
#'

permute_scan_data <- function(ind_sex_permutation = c("Male", "Female"),
                              nn_sex_permutation = c("Male", "Female"),
                              scandata_df,
                              presence_df,
                              n_iter_per_permutation,
                              show_i = TRUE){

  # Create new dataframe
  scan_permuted <- scandata_df

  # Change NNs i n_iter_per_permutation times
  for(i in 1:n_iter_per_permutation){
    if(show_i) cat("\rChanging NN number ", i, "/", n_iter_per_permutation)
    # Pick line 1 and check if Ind + NN adult and have correct Sex class (as
    # defined by ind_sex_permutation and nn_sex_permutation)
    l1check <- FALSE
    while(isFALSE(l1check)){
      l1 <- sample(nrow(scan_permuted), 1)
      row1 <- scan_permuted[l1,]
      if(row1$Ind_Adult == 1){
        if(row1$Ind_Sex %in% ind_sex_permutation){
          if(row1$NN_Adult == 1){
            if(row1$NN_Sex %in% nn_sex_permutation){
              l1check = TRUE
            }}}}
    }

    # Pick line 2 and check if Ind + NN adult and have correct Sex class
    # Then check if NN2 present at date of row1 and NN1 present at date of row2
    l2check <- FALSE
    while(isFALSE(l2check)){
      l2 <- sample(nrow(scan_permuted), 1)
      row2 <- scan_permuted[l2,]
      if(row2$Ind_Adult == 1){
        if(row2$Ind_Sex %in% ind_sex_permutation){
          if(row2$NN_Adult == 1){
            if(row2$NN_Sex %in% nn_sex_permutation){
              # Check if NN2 was present and adult at date of row1
              if(filter(presence_df,
                        Name == row2$NN,
                        YearOf == row1$YearOf,
                        MonthOf == row1$MonthOf)$Present_As_Adult == 1){
                # Check if NN1 was present and adult at date of row2
                if(filter(presence_df,
                          Name == row1$NN,
                          YearOf == row2$YearOf,
                          MonthOf == row2$MonthOf)$Present_As_Adult == 1){
                  l2check = TRUE
                }}}}}}
    }

    # Swap NNs
    scan_permuted[l1,]$NN <- row2$NN
    scan_permuted[l1,]$NN_Sex <- row2$NN_Sex
    scan_permuted[l2,]$NN <- row1$NN
    scan_permuted[l2,]$NN_Sex <- row1$NN_Sex
  }

  # Return permuted data
  return(scan_permuted)
}

