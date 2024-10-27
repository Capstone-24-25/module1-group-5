library(tidyverse)


# get names
var_names <- read_csv("data/biomarker-raw.csv",
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

# function for trimming outliers (good idea??)
trim <- function(x, .at) {
  x[abs(x) > .at] <- sign(x[abs(x) > .at]) * .at
  # sign() is used to preserve the sign of original value
  return(x)
}

# read in data
biomarker_clean <- read_csv("data/biomarker-raw.csv",
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
  # log transform, center and scale, and trim
  # accross apply the functions to all cols at once, but exclude group and ados here.
  mutate(across(
    .cols = -c(group, ados),
    ~ trim(scale(log10(.x))[, 1], .at = 3)
  )) %>%
  # reorder columns
  select(group, ados, everything())

# export as r binary
save(
  list = "biomarker_clean",
  file = "data/biomarker-clean.RData"
)
