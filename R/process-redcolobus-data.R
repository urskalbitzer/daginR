#' Determine presence of adult individuals in the Red colobus group for the specified period of time
#'
#' This function uses scan data in combination with data on the first-adult-month of individual
#' to determine for each month whether an individual was adult and present
#'
#' Default is that individuals are considered absent if recorded less than \code{n_min} (default 5) times per month for 3 months in a row
#'
#'
#' @param scandata Scan data frame
#' @param first_month_adult_df Dataframe with first-month_adult
#' @param from_date,to_date Start and end date of analyzed period of time. Format should be "yyyy-mm-dd"
#' @param n_min Min. number of times and individual has to be recorded (Individual + NN) to be considered present
#' @param details Provide details for each individual and each month?
#'
#' @export
#'
#' @examples
#'
#'

determine_rc_presence <- function(scandata, first_month_adult_df,
                         from_date, to_date, n_min = 5, details = FALSE){
  ind_per_month <- scandata %>%
    group_by(YearOf, MonthOf, Individual) %>%
    summarize(ind_count = n()) %>%
    ungroup()

  nn_per_month <- scandata %>%
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
    mutate(Present_Adult = if_else(Present == 1 & Adult == 1, 1, 0))

  if(!details){
    presence_month <- presence_month %>%
      select(-count_total, -Adult, -Present)
  }
  return(presence_month)
}
