####
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(caret)
library(ggpubr)
library(rstatix)
library(patchwork)

data1 = read.csv("exercise_hrv.csv")
str(data1)

data1 = data1 %>%
  mutate(training_level = as.factor(training_level), sex = as.factor(sex),
         sleep_quality = as.factor(sleep_quality))
str(data1)

# Check effect of sleep quality and training level on hrv values and resting hr values

# 1 - Check for hrv values

data1 %>%
  ggplot(aes( x = training_level, y = hrv_rmssd_ms, fill = sleep_quality))+
  geom_boxplot()+
  scale_x_discrete(limits = c("None", "Low", "Moderate", "High"))+
  theme_minimal(base_size = 8)+
  theme(legend.position = "bottom")

stats = data1 %>%
  group_by(training_level, sleep_quality) %>% 
  summarise(count = n(), mean = mean(hrv_rmssd_ms), sd = sd(hrv_rmssd_ms)) %>%
  as.data.frame()
stats 

anova_two = aov(hrv_rmssd_ms ~ training_level * sleep_quality, data = data1)

# Assumption 1 - Independent samples: Yes, different individuals
# Assumption 2 - Normality of residuals

res = anova_two$residuals
par(mfrow = c(1,3))
hist(res)
boxplot(res)
qqnorm(res)
qqline(res)
r
# Shapiro test
anova_two$xlevels$training_level

h1 <- ggplot(data = anova_two, aes(y = residuals, fill = xlevels$training_level)) +
  geom_boxplot(color = "black", position = 'identity') + 
  labs(title = "Histogram", y = "Phenotype", x = "Height") +
  facet_grid(.~ anova_two$xlevels$training_level) # allows me to split by Females and Males
h1


# Create data frame for plotting
plot_data <- data.frame(
  residuals = residuals(anova_two),
  training_level = anova_two$model$training_level
)

# Plot
ggplot(plot_data, aes(x = training_level, y = residuals, fill = training_level)) +
  geom_boxplot(color = "black") +
  labs(title = "Residuals by Training Level", 
       y = "Residuals", 
       x = "Training Level")



library(ggplot2)
library(gridExtra)

# Data frame for plots
plot_data <- data.frame(
  residuals = residuals(anova_two),
  training_level = anova_two$model$training_level
)

# Boxplot by group
p1 <- ggplot(plot_data, aes(x = training_level, y = residuals, fill = training_level)) +
  geom_boxplot(color = "black") +
  labs(title = "Residuals by Training Level", x = "Training Level", y = "Residuals")

# Histogram of residuals
p2 <- ggplot(plot_data, aes(x = residuals, fill = training_level)) +
  geom_histogram(color = "black", bins = 20, alpha = 0.6, position = "identity") +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Count")

# QQ plot with reference line
p3 <- ggplot(plot_data, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Residuals")

# Arrange all plots
gridExtra::grid.arrange(p1, p2, p3, ncol = 2)








shapiro.test(res)

# Normality fails, but since we have a big sample size (>30 for each group), we can ignore this.

# Assumption 3 - Equality of variances

levenetest <- levene_test(hrv_rmssd_ms ~ training_level * sleep_quality, data = data1) %>%
  as.data.frame()
levenetest$p

# 0.66 > 0.05, we can conclude that there is not a difference between the variances of the two 
# groups

#4- Outliers: There should be no significant outliers in any group

data1 %>% 
  identify_outliers(hrv_rmssd_ms)

# We have some extreme outliers

# See how many extreme outliers we have

data1 %>%
  identify_outliers(hrv_rmssd_ms) %>%
  filter(is.extreme == TRUE)

# We have 6 extreme outliers with ids 313, 395, 439, 526, 654, 928

# Let's remove them

data1_clean = data1 %>%
  filter(!id %in% c(313, 395, 439, 526, 654, 928))
str(data1_clean)

# Checking for new outliers

data1_clean %>% 
  identify_outliers(hrv_rmssd_ms) %>%
  filter(is.extreme == TRUE)

# Removing the new outlier

data1_clean2 = data1_clean %>%
  filter(!id %in% c(134))
str(data1_clean2)

data1_clean2 %>% 
  identify_outliers(hrv_rmssd_ms) %>%
  filter(is.extreme == TRUE)

# No more extreme outliers

anova_two_clean = aov(hrv_rmssd_ms ~ training_level + sleep_quality, data = data1_clean2) # + Because sleep and train don't relate

res_clean = anova_two_clean$residuals
par(mfrow = c(1,3))
hist(res_clean)
boxplot(res_clean)
qqnorm(res_clean)
qqline(res_clean)

shapiro.test(res_clean)

# Normality fails, but since we have a big sample size (>30 for each group), we can ignore this.

levenetest_clean <- levene_test(hrv_rmssd_ms ~ training_level * sleep_quality, data = data1_clean2) %>%
  as.data.frame()
levenetest_clean$p

# 0.58 > 0.05, we can conclude that there is not a difference between the variances of the two 
# groups

### Post Hoc ###

summary(anova_two_clean)
TUKEY <- TukeyHSD(anova_two_clean)
TUKEY

# There is no significant relation between training and sleep, but there is between hrv values and both training and sleep.

with(data1_clean2, interaction.plot(sleep_quality, training_level, hrv_rmssd_ms, 
                          type = "b", col = c(2,3), leg.bty ="o", 
                          lwd = 2, pch = c(18,24), 
                          xlab = "Sleep and Training", ylab = "HRV level") )

# 1 - Check for resting HR values

str(data1)

data1 %>%
  ggplot(aes( x = training_level, y = resting_hr_bpm, fill = sleep_quality))+
  geom_boxplot()+
  scale_x_discrete(limits = c("None", "Low", "Moderate", "High"))+
  theme_minimal(base_size = 8)+
  theme(legend.position = "bottom")

stats2 = data1 %>%
  group_by(training_level, sleep_quality) %>% 
  summarise(count = n(), mean = mean(resting_hr_bpm), sd = sd(resting_hr_bpm)) %>%
  as.data.frame()
stats2 

anova_two2 = aov(resting_hr_bpm ~ training_level + sleep_quality, data = data1) # "+" For same reason as above

# Assumption 1 - Independent samples: Yes, different individuals
# Assumption 2 - Normality of residuals

res2 = anova_two2$residuals
par(mfrow = c(1,3))
hist(res2)
boxplot(res2)
qqnorm(res2)
qqline(res2)

# Shapiro test

shapiro.test(res2)

# Almost perfect normality.

# Assumption 3 - Equality of variances

levenetest2 <- levene_test(resting_hr_bpm ~ training_level * sleep_quality, data = data1) %>%
  as.data.frame()
levenetest2$p

# 0.78 > 0.05, we can conclude that there is not a difference between the variances of the two 
# groups

#4- Outliers: There should be no significant outliers in any group

data1 %>% 
  identify_outliers(resting_hr_bpm)

# We have no extreme outliers

### Post Hoc ###

summary(anova_two2)
TUKEY <- TukeyHSD(anova_two2)
TUKEY

# There is no significant relation between training and sleep.

with(data1, interaction.plot(sleep_quality, training_level, resting_hr_bpm, 
                                    type = "b", col = c(2,3), leg.bty ="o", 
                                    lwd = 2, pch = c(18,24), 
                                    xlab = "Sleep and Training", ylab = "HRV level") )

