library(tidyverse)

# get names
var_names <- read_csv("C:/Users/joshc/Downloads/biomarker-raw.csv",
                      # prevent interpreting the first row as header
                      col_names = F,
                      # read only 2 rows
                      n_max = 2,
                      # exclude first two columns
                      col_select = -(1:2)
) %>%
  # switching rows and columns
  t() %>%
  as_tibble() %>%
  # name, abbreviate are the new name. V1, V2 are the the default old neames
  rename(
    name = V1,
    abbreviation = V2
  ) %>%
  na.omit()



# read in data
biomarker_clean <- read_csv("C:/Users/joshc/Downloads/biomarker-raw.csv",
                            skip = 2,
                            col_select = -2L,
                            col_names = c(
                              "group",
                              "empty",
                              pull(var_names, abbreviation),
                              "ados"
                            ),
                            # Treat '-' and '' as NA.
                            na = c("-", "")
) %>%
  filter(!is.na(group)) %>%
  mutate(across(
    .cols = -c(group, ados),
    ~ scale(log10(.x))[, 1]
  )) %>%
  select(group, ados, everything())

outlier_threshold <- 3 
outlier_df <- biomarker_clean %>%
  mutate(across(
    .cols = -c(group, ados),  
    ~ ifelse(abs(.) > outlier_threshold, ., NA) 
  ))

outlier_summary <- outlier_df %>%
  select(group, ados, everything()) %>%
  pivot_longer(-c(group, ados), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(value)) %>%  # Keep only outliers
  group_by(group, ados, biomarker) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 2)

view(outlier_summary)
