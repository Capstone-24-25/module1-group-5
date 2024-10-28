# What do you imagine is the reason for log-transforming the protein levels in biomarker-raw.csv? 
#   (Hint: look at the distribution of raw values for a sample of proteins.)
library(ggplot2)
library(dplyr)

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


# install.packages("gridExtra")
library(gridExtra)



grid.arrange(p11, p12,  ncol = 2)
grid.arrange(p21, p22,  ncol = 2)
grid.arrange(p31, p32,  ncol = 2)
grid.arrange(p41, p42,  ncol = 2)
grid.arrange(p51, p52,  ncol = 2)
grid.arrange(p61, p62,  ncol = 2)
grid.arrange(p71, p72,  ncol = 2)


# From observing the distribution of a sample of the proteins from the raw file, it is clear that 
# many of them have skewed distributions. These skewed distributions could be due to high 
# variability in the protein levels which leads to outliers that are affecting our data and 
# predictions. To improve our model and predictions, we can log transform the proteins to help
# normalize the distributions. By doing this, it helps to improve our models performance since
# machine learning techniques like random forest will give you better predictions when the 
# input data is more normally distributed. As we can observe from the difference in the raw and
# the log transformed histograms for this sample of proteins, the transformation helps to 
# normalize the distributions of the raw proteins in our dataset.


