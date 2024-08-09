# Load necessary libraries
library(readr)  # for read_tsv
library(ggplot2)
library(data.table)
library(dplyr)
# Replace 'path/to/your/data.tsv' with the actual path to your TSV file
B6vsnonB6 <- read.delim("\\\\wsl.localhost/Ubuntu/home/abed/B6vsnonB6.tsv", stringsAsFactors = FALSE)
# Plot histograms for each column
par(mfrow = c(1, 3))  # Arrange plots in 1 row, 3 columns

hist(B6vsnonB6$GT_0.0, breaks = 50, main = "GT_0/0 Distribution", xlab = "GT_0/0")
hist(B6vsnonB6$GT_0.1, breaks = 50, main = "GT_0/1 Distribution", xlab = "GT_0/1")
hist(B6vsnonB6$GT_1.1, breaks = 50,main = "GT_1/1 Distribution", xlab = "GT_1/1")

par(mfrow = c(1, 1))  # Reset plotting parameters to default
########
# Plot boxplots for each column
boxplot(B6vsnonB6[, c("GT_0.0", "GT_0.1", "GT_1.1")], 
        main = "Genotype Distributions", 
        names = c("GT_0/0", "GT_0/1", "GT_1/1"),
        ylab = "Counts")
##
# Plot density plots for each column
par(mfrow = c(1, 3))  # Arrange plots in 1 row, 3 columns

plot(density(B6vsnonB6$GT_0.0), main = "GT_0/0 Density")
plot(density(B6vsnonB6$GT_0.1), main = "GT_0/1 Density")
plot(density(B6vsnonB6$GT_1.1), main = "GT_1/1 Density")

par(mfrow = c(1, 1))  # Reset plotting parameters to default
print(is.data.table(B6vsnonB6))
B6vsnonB6 <- B6vsnonB6 %>%
  mutate(B6ratio = GT_0.0 / (GT_0.0 + GT_0.1 + GT_1.1))
hist(B6vsnonB6$B6ratio, breaks = 50,main = "B6 Ratio Across All Regions", xlab = "Ratio B6")
summary(B6vsnonB6$B6ratio)
print(summary(B6vsnonB6$B6ratio), mean)
std_dev_B6ratio <- sd(B6vsnonB6$B6ratio, na.rm = TRUE)
print(2*std_dev_B6ratio + mean )
filtered_B6vsnonB6 <- B6vsnonB6 %>%
  filter(B6ratio >= 2*std_dev_B6ratio +0.4849)
