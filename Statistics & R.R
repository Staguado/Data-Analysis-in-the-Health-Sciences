# By: Santiago Taguado 
# Statistics & R
# Harvard Course

# Packages Required
library(dplyr)
library(rafalib)
library(downloader)
library(dslabs)
library(UsingR)

# Downloading data from my computer
data <- read.csv(file.choose(),header = T)
head(data)

# Array Extraction of Bodyweight Data
data$Bodyweight[11]
length(data$Bodyweight)
mean(data[13:24])

# Random Sampling Example
set.seed(1)
i <- sample(13:24, 1)
data$Bodyweight[i]

#---------------------------------------------------------------------------------------------------
# Comparing Sleeping Rates of Rodents & Primates

# Downloading Sleeping Rate Data from Github Repository
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- basename(url)
download(url,filename)
data_2 <- read.csv(filename)
class(data_2)

# Extraction of Primate Sleeping Times
only_primate_data <- data_2 %>% filter(order == "Primates")
nrow(only_primate_data)
class(only_primate_data)

# Comparison of Rodents & Primates Sleeping Time
data_2 %>% filter(order == "Rodentia") %>% summarize(mean(sleep_total))
data_2 %>% filter(order=="Primates") %>% summarize( mean( sleep_total) )

#---------------------------------------------------------------------------------------------------

# On gg plots & Forloop Visualizations

# Downloading data from my computer
load("/Users/santiagotaguado/Desktop/R-Coding/R & Statistics Project/Data/skew.rdata")
head(dat)

# Visualizations QQ Plots: Rudimentary 
new_data <- as.vector(dat)
class(new_data)
V1 <- new_data[1:1000]
V2 <- new_data[1001:2000]
V3 <- new_data[2001:3000]
V4 <- new_data[3001:4000]
V5 <- new_data[4001:5000]
V6 <- new_data[5001:6000]
V7 <- new_data[6001:7000]
V8 <- new_data[7001:8000]
V9 <- new_data[8001:9000]
mypar(3,3)
qqnorm(V1)
qqline(V1)
qqnorm(V2)
qqline(V2)
qqnorm(V3)
qqline(V3)
qqnorm(V4)
qqline(V4)
qqnorm(V5)
qqline(V5)
qqnorm(V6)
qqline(V6)
qqnorm(V7)
qqline(V7)
qqnorm(V8)
qqline(V8)
qqnorm(V9)
qqline(V9)

# Histogram Visualization
mypar(2,1)
hist(V9)
hist(V4)

# Performing QQ-Plots: For-loop
mypar(3,3)
for (i in 1:9) {
  qqnorm(dat[,i])
  qqline(dat[,i])
}

#---------------------------------------------------------------------------------------------------

# Insect Spray Analysis

insect_boxplot_data <- InsectSprays

# Boxplot of Spray
mypar(1,1)
boxplot(InsectSprays$count ~ InsectSprays$spray)

#---------------------------------------------------------------------------------------------------

# New York Marathon Analysis

data(nym.2002, package="UsingR")
data(nym.2002)

x <- nym.2002 %>% filter(gender == "Male" & time) 
hist(x$time,breaks = 20, main = paste("Male Running Time"))

y <- nym.2002 %>% filter(gender == "Female" & time) 
hist(y$time, breaks = 20, main = paste("Female Running Time"))

nym.2002 %>% filter(gender == "Male") %>% ggplot(aes(data = time, aes(nym.2002$time))) + geom_histogram()

#---------------------------------------------------------------------------------------------------

# Female Mice Control Population Data Analysis

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

# Mean of Data
mean(x)
set.seed(1)

# Random Sample Mean
set.seed(1)
y <- sample(x,5)
abs(mean(y) - mean(x))

# Random Sample Mean 2
set.seed(5)
y <- sample(x,5)
abs(mean(y) - mean(x))

#-------------------------------------------------------------------------------------------------

# Exercise: Set the seed at 1, then using a for-loop take a random sample of 5 mice 1,000 times. 
#Save these averages.

# What proportion of these 1,000 averages are more than 1 gram away from the average of x ?

# Downloading Data
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

# Running a For-loop 1000 times of a random sample of 5 mice
set.seed(1)
nulls <- vector("numeric",n)
n <- 1000
for (i in 1:n){
  randoms <- sample(x,50)
  nulls[i] <- mean(randoms)
}
mean(x) +1

# Probability of Controls Shifting by 1 gram
(sum(nulls > mean(x)+1) + sum(nulls < mean(x)-1))/1000

# As the sample increases in size, the probability that controls are away from the 
# mean by more than 1 gram decreases.

#-------------------------------------------------------------------------------------------------

# Comparison of Random Sample Size Averages

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

# Random Sample of 5
set.seed(1)
n <- 1000
averages5 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,5)
  averages5[i] <- mean(X)
}

# Random Sample of 50
set.seed(1)
n <- 1000
averages50 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,50)
  averages50[i] <- mean(X)
}

# Histogram Comparison
mypar(2)
hist(averages5)
hist(averages50)

# Descriptive Statistics Assuming Normal Distribution: Data between 23-25 grams
mean(averages50)
sd(averages50)
pnorm(25, mean(averages50),sd(averages50)) - pnorm(23, mean(averages50),sd(averages50))

#------------------------------------------------------------------------------------------

# Chow & High Fat Diet Comparison

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 

# Using dplyr to create a vector x with the body weight of all males on the control (chow) diet.
dat <- na.omit( dat )
x <- dat %>% filter(Diet == "chow" & Sex == "M")

# Mean of Chow Diet
x_mean <- mean(x$Bodyweight)
popsd(x$Bodyweight)

# Random Sample of Chow Diet on Males
set.seed(1)
random_sample_25 <- sample(x$Bodyweight, 25)
mean(random_sample_25)

# Random Sample of High Fat on Females
y <- dat %>% filter(Sex == "F" & Diet == "hf")
y_mean <- mean(y$Bodyweight)
popsd(y$Bodyweight)
random_sample_25_y <- sample(y$Bodyweight, 25)
mean(random_sample_25_y)

#------------------------------------------------------------------------------------------

# Chow & High Fat Visualizations

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- na.omit( read.csv(filename) )

# On 95% and 99% Confidence Intervals
y <- dat %>% filter(Diet == "chow" & Sex == "M")
y_mean <- mean(y$Bodyweight)
y_sd <- popsd(y$Bodyweight)
upper_limit <- y_mean + y_sd
lower_limit <- y_mean - y_sd
upper_limit_2 <- y_mean + 2*y_sd
lower_limit_2 <- y_mean - 2*y_sd
upper_limit_3 <- y_mean + 3*y_sd
lower_limit_3 <- y_mean - 3*y_sd

1 - ((sum(y$Bodyweight > upper_limit_3) + sum(y$Bodyweight < lower_limit_3))/223)

# Creating 4- Graphs of q-q plots Side_by_Side
mypar(2,2)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)

# Creating a histogram and q-q plots Side_by_Side
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
set.seed(1)
avgs <- replicate(10000, mean( sample(y, 25)))

# QQ Plots and Hisograms
mypar(1,2)
hist(avgs)
qqnorm(avgs)
qqline(avgs)


#------------------------------------------------------------------------------------------

# Replication of Random Variables

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if(!file.exists("femaleMiceWeights.csv")) download(url,destfile=filename)
dat <- read.csv(filename)

# Replication Function Example
set.seed(1)
n <- 30
sides <- 6
p <- 0.5
zs <- replicate(10000,{
  x <- sample(1:sides,n,replace=TRUE)
  (mean(x==6) - p) / sqrt(p*(1-p)/n)
}) 
qqnorm(zs)
abline(0,1) #confirmation it's well approximated with normal distribution
mean(abs(zs) > 2)

# Running 4 Types of Distributions
ps <- c(0.5,0.5,0.01,0.01)
ns <- c(5,30,30,100)
mypar(4,2)
for(i in 1:4){
  p <- ps[i]
  sides <- 1/p
  n <- ns[i]
  zs <- replicate(10000,{
    x <- sample(1:sides,n,replace=TRUE)
    (mean(x==1) - p) / sqrt(p*(1-p)/n)
  }) 
  hist(zs,nclass=7)
  qqnorm(zs)
  abline(0,1)
}

#------------------------------------------------------------------------------------------