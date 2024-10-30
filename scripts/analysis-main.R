# Primary script that contains the code that produces our findings

library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

biomarker_raw <- read.csv('data/biomarker-raw.csv')
biomarker_raw

biomarker_num <- biomarker_raw[-1,]
biomarker_num <- biomarker_num %>% mutate_all(~ as.numeric(as.character(.)))
str(biomarker_num)

# plot histogram of each protein to see distribution of raw data
p11 <- ggplot(biomarker_num, aes(x = Gamma.enolase)) +
  geom_histogram(bins = 15)

p21 <- ggplot(biomarker_num, aes(x = E3.ubiquitin.protein.ligase.CHIP)) +
  geom_histogram(bins = 10)

p31 <- ggplot(biomarker_num, aes(x = CCAAT.enhancer.binding.protein.beta)) +
  geom_histogram(bins = 15)

p41 <- ggplot(biomarker_num, aes(x = Vitronectin)) +
  geom_histogram(bins = 20)

p51 <- ggplot(biomarker_num, aes(x = Histone.H3.1)) +
  geom_histogram(bins = 15)

p61 <- ggplot(biomarker_num, aes(x = Semaphorin.5A)) +
  geom_histogram(bins = 20)

p71 <- ggplot(biomarker_num, aes(x = Protein.S100.A6)) +
  geom_histogram(bins = 15)



# log transform the raw proteins and plot
p12 <- ggplot(biomarker_num, aes(x = log(Gamma.enolase))) +
  geom_histogram(bins = 15)

p22 <- ggplot(biomarker_num, aes(x = log(E3.ubiquitin.protein.ligase.CHIP))) +
  geom_histogram(bins = 10)

p32 <- ggplot(biomarker_num, aes(x = log(CCAAT.enhancer.binding.protein.beta))) +
  geom_histogram(bins = 15)

p42 <- ggplot(biomarker_num, aes(x = log(Vitronectin))) +
  geom_histogram(bins = 20)

p52 <- ggplot(biomarker_num, aes(x = log(Histone.H3.1))) +
  geom_histogram(bins = 15)

p62 <- ggplot(biomarker_num, aes(x = log(Semaphorin.5A))) +
  geom_histogram(bins = 20)

p72 <- ggplot(biomarker_num, aes(x = log(Protein.S100.A6))) +
  geom_histogram(bins = 15)

# plot untransformed and transformed histograms side by side
grid.arrange(p11, p12,  ncol = 2)
grid.arrange(p21, p22,  ncol = 2)
grid.arrange(p31, p32,  ncol = 2)
grid.arrange(p41, p42,  ncol = 2)
grid.arrange(p51, p52,  ncol = 2)
grid.arrange(p61, p62,  ncol = 2)
grid.arrange(p71, p72,  ncol = 2)

#################

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


##########################

load("/Users/wentaozhang/Downloads/biomarker-clean.RData")

## MULTIPLE TESTING


# function to compute tests
test_fn <- function(.df) {
  t_test(.df,
         # compute t test on level, comparing the two levels of group
         formula = level ~ group,
         # mean ASD - mean TD
         order = c("ASD", "TD"),
         alternative = "two-sided",
         var.equal = F
  )
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group,
               # the group column should be left as-is
               names_to = "protein",
               # Each unique column name in the original data
               # becomes a row in the protein column
               values_to = "level"
               # he values from each column corresponding to
               # each protein will be stored in the level column.
  ) %>%
  # nest by protein
  nest(data = c(level, group)) %>%
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(
    m = n(),
    hm = log(m) + 1 / (2 * m) - digamma(1),
    rank = row_number(),
    p.adj = m * hm * p_value / rank
  )

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)


## RANDOM FOREST

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>%
  pull(group) %>%
  factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(
  x = predictors,
  y = response,
  ntree = 1000,
  importance = T
)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)


### LOGISTIC REGRESSION


# select subset of interest 
# use union
proteins_sstar <- union(proteins_s1, proteins_s2)


biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == "ASD")) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ .,
           data = training(biomarker_split),
           family = "binomial"
)


# evaluate errors on test set
class_metrics <- metric_set(
  sensitivity,
  specificity,
  accuracy,
  roc_auc
)

testing(biomarker_split) %>%
  add_predictions(fit, type = "response") %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(
    estimate = est,
    truth = tr_c, pred,
    event_level = "second"
  )

####################






