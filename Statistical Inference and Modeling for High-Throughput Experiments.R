#-----------------------------------------------------------------------------------------------

# By Santiago Taguado Menza
# Statistical Inference & Modeling for High-throughput Experiments
# Harvard Course

#-----------------------------------------------------------------------------------------------

# Testing P-values without Adjustment

# Using Swirl
dat[ dat[,3] > k , ]
install.packages("swirl")
library(swirl)
swirl()

# Downloading Training Data
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables

exams_on_july_27_2005<- sampleInfo %>% filter(date == "2005-06-27")
nrow(exams_on_july_27_2005)

# Alternatively, use the following:

sum(sampleInfo$date=="2005-06-27")

# How many of the genes represented in this particular technology are on chromosome Y?
# Remove features that have NA in the column of interest.

gene_Y <- geneAnnotation %>% filter(CHR == "chrY")
nrow(gene_Y)

# Alternatively, use the following:

sum(geneAnnotation$CHR=="chrY",na.rm=TRUE)

# We need the na.rm=TRUE option because some features are controls and have NA in the CHR column.
# What is the log expression value of the for gene ARPC1A on the one subject that was measured on 2005-06-10?

exams_on_june_10_2005<- sampleInfo %>% filter(date == "2005-06-10")
exams_on_june_10_2005

ARPCIA_gene <- geneAnnotation %>% filter(SYMBOL == "ARPC1A")                    
ARPCIA_gene
# GSM136727.CEL.gz
# 200950_at
geneExpression["200950_at","GSM136727.CEL.gz"]

# Alternatively, use the following:

i = which(geneAnnotation$SYMBOL=="ARPC1A")
j = which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]

median_column <- apply(geneExpression,MARGIN = 2, median)
median(median_column)

# Note that we can also use the colMedians() function from the matrixStats package.
g <- factor(sampleInfo$group)
p_values <- function(e) {t.test(e[g==1], e[g==0])$p.value}
p_values_data <- apply(geneExpression,1, p_values)
min(p_values_data)

myttest <- function(e,group){
  x <- e[group==1]
  y <- e[group==0]
  return( t.test(x,y)$p.value )
}
g <- factor(sampleInfo$group)
pvals <- apply(geneExpression,1,myttest, group=g)
min( pvals )

#-----------------------------------------------------------------------------------------------

# Testing P-values without Adjustment Part 2

# Using Swirl
dat[ dat[,3] > k , ]
install.packages("swirl")
library(swirl)
swirl()

# Downloading Training Data
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables

exams_on_july_27_2005<- sampleInfo %>% filter(date == "2005-06-27")
nrow(exams_on_july_27_2005)

# Alternatively, use the following:

sum(sampleInfo$date=="2005-06-27")

# How many of the genes represented in this particular technology are on chromosome Y?
# Remove features that have NA in the column of interest.

gene_Y <- geneAnnotation %>% filter(CHR == "chrY")
nrow(gene_Y)

# Alternatively, use the following:

sum(geneAnnotation$CHR=="chrY",na.rm=TRUE)

# We need the na.rm=TRUE option because some features are controls and have NA in the CHR column.
# What is the log expression value of the for gene ARPC1A on the one subject that was measured on 2005-06-10?

exams_on_june_10_2005<- sampleInfo %>% filter(date == "2005-06-10")
exams_on_june_10_2005

ARPCIA_gene <- geneAnnotation %>% filter(SYMBOL == "ARPC1A")                    
ARPCIA_gene
# GSM136727.CEL.gz
# 200950_at
geneExpression["200950_at","GSM136727.CEL.gz"]

# Alternatively, use the following:

i = which(geneAnnotation$SYMBOL=="ARPC1A")
j = which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]

median_column <- apply(geneExpression,MARGIN = 2, median)
median(median_column)

# Note that we can also use the colMedians() function from the matrixStats package.
g <- factor(sampleInfo$group)
p_values <- function(e) {t.test(e[g==1], e[g==0])$p.value}
p_values_data <- apply(geneExpression,1, p_values)
min(p_values_data)

myttest <- function(e,group){
  x <- e[group==1]
  y <- e[group==0]
  return( t.test(x,y)$p.value )
}
g <- factor(sampleInfo$group)
pvals <- apply(geneExpression,1,myttest, group=g)
min( pvals )

#-----------------------------------------------------------------------------------------------

# Computationally Proving P-value are Random Variables

# Note that we will later learn about the rowttests() function from the genefilter package, 
# which performs this operation.
# Inference in Practice Exercises #1

# These exercises will help clarify that p-values are random variables and some of the properties 
# of these p-values. Note that just like the sample average is a random variable because it is based 
#  on a random sample, the p-values are based on random variables (sample mean and sample standard 
# deviation for example) and thus it is also a random variable.

# To see this, let's see how p-values change when we take different samples.

set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
  control = sample(population[,1],12)
  treatment = sample(population[,1],12)
  t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)

# What proportion of the p-values is below 0.05?

mean(pvals < 0.05)

# What proportion of the p-values are below 0.01?

mean(pvals < 0.01)

# Assume you are testing the effectiveness of 20 diets on mice weight. 
# For each of the 20 diets, you run an experiment with 10 control mice and 10 treated mice. 
# Assume the null hypothesis that the diet has no effect is true for all 20 diets and 
# that mice weights follow a normal distribution with mean 30 grams and a standard deviation of 2 grams. 
# Run a Monte Carlo simulation for one of these studies:

cases = rnorm(10,30,2)
controls = rnorm(10,30,2)
t.test(cases,controls)$p.value

# Now run a Monte Carlo simulation imitating the results for the experiment for all 20 diets. 
# If you set the seed at 100, set.seed(100), and use the same code as above inside a call to replicate(), 
# how many of the p-values (number not proportion) are below 0.05?

set.seed(100)
pvals_20_diets <- replicate(20, {cases = rnorm(10,30,2)
controls = rnorm(10,30,2)
p_values = t.test(cases,controls)$p.value
sum(p_values <= 0.05)})


set.seed(100)
mousePval = function(){
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  return(t.test(cases,controls)$p.value)
}
pvals = c()
for (i in 1:20){
  pvals[i] = mousePval()
}
length(pvals[pvals < 0.05])

# Now create a simulation to learn about the distribution of the number of p-values that are less than 0.05. 
# In the previous question, we ran the 20 diet experiment once. 
# Now we will run these 20 experiments 1,000 times and each time save the number of p-values that are 
# less than 0.05.

# Set the seed at 100 again, set.seed(100), run the code from the previous question 1,000 times, 
# and save the number of times the p-value is less than 0.05 for each of the 1,000 instances.

# What is the average of these 1,000 numbers? 
# Note that this is the expected number of tests (out of the 20 we run) 
# that we will reject when the null is true.

set.seed(100)
pvals_20_diets_1000 <- replicate(1000,{
  pvals = replicate(20,{
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  t.test(cases,controls)$p.value
  })
  sum(pvals <= 0.05)
})

table(pvals_20_diets_1000)
mean(pvals_20_diets_1000)

# Inference in Practice Exercises #5
# 1 point possible (graded)
# Note that what the answer to question #4 says is that on average, 
# we expect some p-value to be 0.05 even when the null is true for all diets.
# Using the same simulation data from the question above, for what proportion of the 1,000 replicates 
# do we reject the null hypothesis at least once (more than 0 false positives)?

mean(plessthan>0) 

# Quiz 1 
# Suppose you plan to run an experiment screening a panel of 30,000 small molecules to determine which 
# ones increase expression of a fluorescent reporter gene. In untreated cells, the reporter gene expression 
# follows a normal distribution with a mean of 8 units and a standard deviation of 2 units. 
# There will be 100 untreated control wells, and each of the 30,000 molecules will be tested in 
# 10 technical replicates. You want to simulate the experiment to figure out how many hits would come 
# out of your screen if the null hypothesis is true for all 30,000 cases.

# Set the seed to 3, then generate the results of the 10 control wells:

set.seed(3)
ctrl = rnorm(100, 8, 2)

# This example code simulates 10 technical replicates for one compound for which the null hypothesis is true:

expt = rnorm(10, 8, 2)
t.test(ctrl, expt)$p.value

# Now set the seed to 4 and use replicate() to simulate 30,000 tests for which the null hypothesis is true. 
# The example code for one compound should go inside your replicate() call. Note that each test will compare 
# the same ctrl vector to a new simulated experimental vector.

# Make a histogram of the p-values from this experiment.
# Question 1 - Quiz 1: Which distribution do these p-values follow most closely?

set.seed(4)
B <- 30000
pval_distribution <- replicate(B, {
  expt = rnorm(10, 8, 2)
  t.test(ctrl, expt)$p.value})

hist((pval_distribution))

# Question 2 - Quiz 1: 

# What proportion of tests have a p-value below 0.05?

sum(pval_distribution < 0.05)/B

# Question 3 - Quiz 1:
# Since this simulation assumes that the null distribution is true for all compounds, 
# any results that have a p-value below a given cutoff will be false positives.

# How many compounds have a p-value below 0.05?

sum(pval_distribution < 0.05)


# Question 4 - Quiz 1: 

# If the p-value cutoff is lowered to 0.001, how many false positive results are there?

sum(pval_distribution < 0.001)

# Assume you are testing the effectiveness of 30 drugs on the white blood cell count of mice. 
# For each of the 30 drugs you run an experiment with 5 control mice and 5 treated mice. 
# Assume the null hypothesis that the drug has no effect is true for all 30 drugs and that white blood cell
# counts follow a normal distribution with mean 7.5 units and a standard deviation of 2.5 units.

# We will analyze the number of significant p-values expected by chance under the null distribution.


# Question 5 - Quiz 1: 

# Set the seed to 28, then run a Monte Carlo simulation for one of these studies by randomly generating 
# white blood cell counts for the cases and controls. 
# Use a t-test to compute the p-value for this simulated study.

# What is the p-value for the one simulated study?

set.seed(28)
cases <- rnorm(5,7.5,2.5)
control <- rnorm(5,7.5,2.5)
t.test(cases, control)$p.value

# Question 6 - Quiz 1:

# Now run a Monte Carlo simulation imitating the results for the experiment for all 30 drugs. 
# Set the seed to 51, set.seed(51), 
# then use your code from the previous question inside of replicate().

# How many of the 30 simulated p-values are below 0.05?

set.seed(51)

# Number of Drugs
B<- 30 
# Simulation of the drugs
Monte_carlos_H0_real <- replicate(B, {
  cases <- rnorm(5,7.5,2.5)
  control <- rnorm(5,7.5,2.5)
  t.test(cases, control)$p.value
})

sum(Monte_carlos_H0_real < 0.05)

# Question 7 - Quiz 1:

# Set the seed to 100, then repeat the simulated experiment 1000 times by using your code from the 
# previous question inside a second replicate() loop. 
# For each experiment, save the number of simulated p-values below 0.05.

# What is the average of the counts of p-values below 0.05?

set.seed(100)
Monte_carlos_null <- replicate(1000, { 
  pvals = replicate(B, {
  cases <- rnorm(5,7.5,2.5)
  control <- rnorm(5,7.5,2.5)
  t.test(cases, control)$p.value})
  sum(pvals < 0.05)})

# Notes: You have to define the variable when doing a replicate within a replicate first. See how
# you first define pvals and then were able to return the values that were less than 0.05. On the 
# previous question, you defined the function as Monte_Carlos_H0_real whereas in this case you defined 
# it as pvals.

# Question 8 - Quiz 1:

# Make a histogram of the p-value counts from Question 7.
# Which of the following is NOT true about the distribution of p-values?

hist(Monte_carlos_null)

# Question 9 - Quiz 1:

# What proportion of simulated experiments have more than 3 p-values below 0.05?

mean(Monte_carlos_null > 3)

# In the previous assessment we saw how the probability of incorrectly rejecting the null for at 
# least one of 20 experiments for which the null is true is well over 5%. 
# Now let's consider a case in which we run thousands of tests as 
# we would do in a high throughput experiment.

# We previously learned that under the null, the probability of a p-value < p is p. 
# If we run 8,793 independent tests, 
# what is the probability of incorrectly rejecting at least one of the null hypotheses?

n <- 8793
p_1 <- 0.95
p_2 <- 0.05

0.95*8793

# Let  𝑃1,…,𝑃8793  be the the p-values (that are random variables).

# Pr(at least one rejection)=1−Pr(no rejections)=1−∏8793𝑖=1Pr(𝑃𝑖>0.05)=1−0.958793≈1
1 - 0.5^8793
# Or if you want to use a simulation:
B<-1000
minpval <- replicate(B, min(runif(8793,0,1))< 0)
mean(minpval>=1)

# Suppose we need to run 8,793 statistical tests and we want to make the probability of a mistake very small,
# say 5%. Using the answer to exercise #2, how small do we have to change the cutoff, previously 0.05, 
# to lower our probability of at least one mistake to be 5%.

# Sidak Procedure: 1 - (1 - alpla)^(1/m)
alpha <- 1 - (1 - 0.05)^(1/8793)

##warning this can take several minutes
##and will only give an approximate answer

B=10000
cutoffs = 10^seq(-7,-4,0.1) ##we know it has to be small
prob = sapply(cutoffs,function(cutoff){
  minpval =replicate(B, min(runif(8793,0,1))<=cutoff)
  mean(minpval>=1)
})
cutoffs[which.min(abs(prob-0.05))]

#-----------------------------------------------------------------------------------------------

# On optimization of code via vectorization:

library(downloader) 
url <-"https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- "femaleControlsPopulation.csv"

if (!file.exists(filename)) download(url,destfile=filename)
set.seed(1)
(population = unlist( read.csv("femaleControlsPopulation.csv") ))

# To give an example of how we can simulate  𝑉  and  𝑆  we constructed a simulation with:

alpha <- 0.05
N <- 12
m <- 10000
p0 <- 0.90 ##10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
(nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1)))
delta <- 3

# We then ran a Monte Carlo simulation by repeating a procedure in which 10,000 tests were run one by 
# one using sapply().

B <- 10 ##number of simulations 
system.time(
  VandS <- replicate(B,{
    calls <- sapply(1:m, function(i){
      control <- sample(population,N)
      treatment <- sample(population,N)
      if(!nullHypothesis[i]) treatment <- treatment + delta
      t.test(treatment,control)$p.val < alpha
    })
    c(sum(nullHypothesis & calls),sum(!nullHypothesis & calls))
  })
)

# In each iteration we checked if that iteration was associated with the null or alternative hypothesis. 
# We did this with the line:

# if(!nullHypothesis[i]) treatment <- treatment + delta

# HOWEVER, in R, operations based on matrices are typically much faster than operations performed within 
# loops or sapply(). We can vectorize the code to make it go much faster. 
# This means that instead of using sapply() to run m tests, we will create a matrix with all data in one 
# call to sample.

# This code runs several times faster than the code above, which is necessary here due to the fact that 
# we will be generating several simulations. 
# Understanding this chunk of code and how it is equivalent to the code above using sapply() will 
# take a you long way in helping you code efficiently in R.

library(genefilter) ##rowttests is here
set.seed(1)

##Define groups to be used with rowttests

(g <- factor( c(rep(0,N),rep(1,N)) ))
B <- 10 ##number of simulations
system.time(  
  VandS <- replicate(B,{
    ##matrix with control data (rows are tests, columns are mice)
    (controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m))
    
    ##matrix with control data (rows are tests, columns are mice)
    (treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m))
    
    ##add effect to 10% of them
    (treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta)
  
    ##combine to form one matrix
    (dat <- cbind(controls,treatments))

    (calls <- rowttests(dat,g)$p.value < alpha)
    
    c(sum(nullHypothesis & calls),sum(!nullHypothesis & calls))
  })
)

#-----------------------------------------------------------------------------------------------

# On P-value Adjustments: Bonferroni & Sidaks Corrections

library(rafalib)
mypar(1,2)
m <- 10000
alphas <- seq(0,0.25,0.01)
bonferroni <- alphas/m
sidaks <- 1 - (1-alphas)^(1/m)

plot(alphas,bonferroni)
plot(alphas,sidaks)

# Rafa's Code
alphas <- seq(0,0.25,0.01)
par(mfrow=c(2,2))
for(m in c(2,10,100,1000)){
  plot(alphas,alphas/m - (1-(1-alphas)^(1/m)),type="l")
  abline(h=0,col=2,lty=2)
}

# As shown in both codes, Bonferroni's is more consertaive showing smaller p-values/cutoff
# values than Zidak's Procedure

# To simulate the p-value results of, say, 8,793 t-tests for which the null is true, 
# we don't actual have to generate the original data. As we learned in class, 
# we can generate p-values from a uniform distribution like this:

pvals <- runif(8793,0,1)

# Using what we have learned, set the cutoff using the Bonferroni correction that 
# guarantees an FWER lower than 0.05 and report back the FWER. 
# Set the seed at 1,set.seed(1), and run 10,000 simulations. 
# Report the Monte Carlo estimate of the FWER below.

# Running a Monte Carlo using a bonferroni correction
set.seed(2)
bonferroni <- 0.05/8793
B<- 10000
FWER_Bonferroni <- replicate(B, {
  pvals <- runif(8793,0,1)
  sum(pvals <= bonferroni)
})
mean(FWER_Bonferroni)

# Alternatively,
set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- alpha/m
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)

# Using the same seed repeat the above for Sidak's cutoff.
# Report the FWER below.
# Running a Monte Carlo using a sidak's correction
set.seed(2)
sidak <- 1 - (1 - 0.05)^(1/8793)
B<- 10000
FWER_Sidak <- replicate(B, {
  pvals <- runif(8793,0,1)
  sum(pvals <= bonferroni)
})
mean(FWER_Sidak)

##if pvals already defined no need to rerun this
set.seed(2)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
pvals
k <- (1-(1-alpha)^(1/m))
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)

# It must be noted that Dr. Irizarry code is much better. My answer showed Bonferroni rate bigger
# than Zidak's. This is not the case because Bonferrni should have a lower FWER. However,
# he accepted my answer because it is dependent on the seed. His isn't.

#-----------------------------------------------------------------------------------------------

# Using the Qvalue package and Introduction to Key BiocManager Packages:

library(devtools)
library(rafalib)
install_github("genomicsclass/GSE5859Subset")
BiocManager::install(c("genefilter", "qvalue"))

library(GSE5859Subset)
data(GSE5859Subset)

library(genefilter)
library(qvalue)

# Question 1
# Compute a p-value for each gene using the function rowttests() from the genefilter package in 
# Bioconductor.

g <- sampleInfo$group
first_rowttest <- rowttests(geneExpression, factor(g))
sum(first_rowttest$p.value<0.05)

# Question 2
# Now applying the bonferroni correction
# Apply the Bonferroni correction to the p-values obtained in question #1 to achieve a FWER of 0.05. 
# How many genes are called significant under this procedure?

sum(first_rowttest$p.value< 0.05/8793)

# Question 3
# Note that the FDR is a property of a list of features, not each specific feature. 
# The q-value relates FDR to an individual feature. To define the q-value we order features we tested 
# by p-value then compute the FDRs for a list with the most significant, 
# the two most significant, the three most significant, etc... 
# The FDR of the list with the, say, m most significant tests is defined as the q-value of the 
# m-th most significant feature. 
# In other words, the q-value of a feature, is the FDR of the biggest list that includes that gene.

# In R, we can compute the q-value using the p.adjust function with the FDR option. 
# Read the help file for p.adjust and then, for our gene expression dataset, 
# compute how many genes achieve an FDR < 0.05

p_adjust <- p.adjust(first_rowttest$p.value, method = "fdr")
sum(p.adjust(first_rowttest$p.value, method = "fdr") < 0.05)

# Question 4

# Now use the qvalue function, in the Bioconductor qvalue package, to estimate q-values using the 
# procedure described by Storey.

# Using this estimate how many genes have q-values below 0.05?

mypar(1,1)
qvalue_list <- qvalue(first_rowttest$p.value,fdr.level = 0.05)
sum(qvalue_list$qvalues < 0.05)
hist(qvalue_list$qvalues)

# Question 5
# Read the help file for qvalue and report the estimated proportion of genes for which the null hypothesis 
# is true  𝜋0=𝑚0/𝑚

qvalue_list$pi0*8793

# Question 6

plot(qvalue_list$qvalues,p_adjust)

# Question 7 

# Create a Monte Carlo Simulation in which you simulate measurements from 8,793 genes for 24 samples: 
# 12 cases and 12 controls.

n <- 24
m <- 8793
mat <- matrix(rnorm(n*m),m,n)

# Now for 500 genes, there is a difference of 2 between cases and controls:

delta <- 2
positives <- 500
mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta

# So the null hypothesis is true for 8793-500 genes. 
# Using the notation from the videos m=8793, m0=8293 and m1=500
# Set the seed at 1, set.seed(1), and run this experiment 1,000 times with a Monte Carlo simulation. 
# For each instance compute p-values using a t-test (using rowttests() in the genefilter package) and 
# create three lists of genes using:
  
# Bonferroni correction to achieve an FWER of 0.05,
# p.adjust() estimates of FDR to achieve an FDR of 0.05, and
# qvalue() estimates of FDR to to achieve an FDR of 0.05.

# For each of these three lists compute the number of false positives in the list and the number of 
# false negatives: genes not in the list that should have been because the null hypothesis is not true 
# (we added 2 to the controls to create the cases).

# What is the false positive rate (false positives divided by m0) if we use Bonferroni?

set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FP1
})
mean(result/(m-positives))

# From the same Monte Carlo simulation as in the question above, what is the false negative rate if we
# use Bonferroni?

set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FN1 <- sum(pvals[1:positives]>0.05/m)
  c(FP1,FN1)
})
mean(result[2,]/(positives))

# Adjusting the qvalue for the false positive rates/false negative rates

set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##p.adjust value
  pvals_adjusted = p.adjust(pvals, method ="fdr")
  FP1 <- sum(pvals_adjusted[-(1:positives)]<=0.05)  
  FN1 <- sum(pvals_adjusted[1:positives]>0.05)
  c(FP1,FN1)
})

mean(result[1,]/(m-positives))
mean(result[2,]/(positives))


# Adjusting the qvalue for the false positive rates/false negative rates using the qvalue formula
set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Q value formula
  pvals_adjusted = qvalue(pvals, fdr.level = 0.05)
  FP1 <- sum(pvals_adjusted$qvalues[-(1:positives)]<=0.05)  
  FN1 <- sum(pvals_adjusted$qvalues[1:positives]>0.05)
  c(FP1,FN1)
})

mean(result[1,]/(m-positives))
mean(result[2,]/(positives))

qvalue_list <- qvalue(first_rowttest$p.value,fdr.level = 0.05)
sum(qvalue_list$qvalues < 0.05)

# Question 1 was listed in the booklet

# Question 2

# A clinical trial of a diagnostic test is performed on 200 people. 
# The null hypothesis is that an individual does not have the disease. 
# 92 people with the disease are correctly labeled as having the disease. 
# 9 people with the disease are incorrectly labeled as healthy when the disease is actually present. 
# 16 healthy people are incorrectly labeled as having the disease. 
# 83 healthy people are correctly labeled as healthy.

# A. How many type I errors are there?

# 16

# B. How many type II errors are there?

# 9

# C. What percentage of healthy people are false positives?

# 16.16 %

# D. What percentage of people with the disease are false negatives?

# 8.91 %

# Question 3
# A certain RNA-seq experiment measures expression of  𝑚=6,319  features.

# Using the Bonferroni correction, what p-value cutoff  𝑘  would control the familywise error rate at
# 𝛼=0.05 ?

p_value_cutoff <- 0.05/6319

# Question 4

# Simulate the results of the RNA-seq experiment from Question 3 assuming the null distribution is true 
# for all features. Set the seed to 11. Use runif() to simulate  𝑚  p-values.

# How many p-values are below the cutoff  𝑘 ?

set.seed(11)
FWER_bonferroni <- runif(6319,0,1)
sum(FWER_bonferroni<= p_value_cutoff)

# Question 5

# Perform a Monte Carlo simulation of the familywise error rate for the RNA-seq experiment in question 3. 
# Set the seed to 12. Use runif() and replicate() to simulate 10,000 sets of  𝑚  p-values. For each set,
# determine how many p-values are below the cutoff  𝑘 . Under the assumption of the null distribution
# ,these are false positives.

# What proportion of simulated experiments have at least one false positive?

set.seed(12)
B <- 10000
FWER_bonferroni <- replicate(B, {
  pvals <- runif(6319,0,1)
  sum(pvals <= p_value_cutoff)
})
sum(FWER_bonferroni > 0)/10000

# Question 6 - 10 

# This is a dataset produced by Bottomly et al., performing RNA-sequencing on two strains of mouse 
# with many biological replicates.

# download Bottomly et al. data

if (!file.exists("bottomly_eset.RData")) download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData", 
                                                       "bottomly_eset.RData")
load("bottomly_eset.RData")

# also load these libraries, which we previously installed from Bioconductor
library(Biobase)
library(genefilter)
library(qvalue)

# These data are stored in an ExpressionSet object. We will learn how to work with these objects in future 
# courses, but for now we can manually extract the gene expression and strain information:

dat = exprs(bottomly.eset)    # gene expression matrix
strain = pData(bottomly.eset)$strain    # strain factor

# dat is a matrix with each row representing a gene, each column representing a sample, and the values 
# representing RNA-seq read counts for a given gene in a given sample. strain is a factor representing 
# the genetic strain of each sample column.

# Question 6

# Use the rowttests() function from the genefilter libaary to calculate p-values for every gene (row) in 
# dat based on strain

library(genefilter)
results <- rowttests(dat,strain)
pvals <- results$p.value
sum(pvals < 0.05, na.rm = pvals)

# Question 7

# Using the Bonferroni correction, what p-value cutoff would be required to ensure a FWER below 0.05?

(p_value_cutoff <- 0.05/nrow(dat))

# Question 8

# How many genes have a p-value below the cutoff determined by the Bonferroni correction?
  
sum(pvals < p_value_cutoff, na.rm = pvals)

# Question 9

# Use p.adjust() with the method="fdr" option to compute q-values in order to determine how many genes 
# are significant at an FDR cutoff of 0.05.

# How many genes have significant q-values at an FDR cutoff of 0.05 when using p.adjust()?

pvals_fdr <- p.adjust(pvals, method = "fdr")
sum(pvals_fdr < 0.05, na.rm = pvals)

# Question 10

# Now try computing q-values with an alternative method, using the qvalue() function from the qvalue package.
# How many genes have significant q-values at an FDR cutoff of 0.05 when using qvalue()?
# You may need to remove NAs from pvals before finding this value.

qvalue_pvals <- qvalue(pvals, fdr.level = 0.05)
qvalue_p_ad <- na.omit(qvalue_pvals$qvalues)
sum(qvalue_p_ad < 0.05)

#----------------------------------------------------------------------------------------------------------

# Binomial and Poisson Distribution:

# Suppose you have an urn with blue and red balls. If  𝑁  balls are selected at random with replacement
# (you put the ball back after you pick it) we can denote the outcomes as random variables  𝑋1,…,𝑋𝑁
# that are 1 or 0. If the proportion of red balls is  𝑝  then the distribution of each of these is:

# Pr(Xi = 1) = p 

# These are also called Bernoulli trials. Note that these random variables are independent because 
# we replace the balls. Flipping a coin is an example of this with p = 0.5.

# You can show that the mean and variance are  𝑝  and  𝑝(1−𝑝)  respectively. The binomial distribut
# ion gives us the distribution of the sum  𝑆𝑁  of these random variables. The probability that we see
# 𝑘  red balls is given by:

# Pr(Sn = k) = (N k) combinatorial * (p ^ k) * ((1 - p)^(N - k))

# In R the function dbinom() gives you this result. The function pbinom() gives us  Pr(𝑆𝑁≤𝑘) .

# This equation has many uses in the life sciences. We give some examples below.
# The probability of conceiving a girl is 0.49. 

# Question # 1

# What is the probability that a family with 4 children has 2 girls and 2 boys (you can assume no twins)?

(p_x <- choose(4,2) * 0.49^2 * 0.51^2)
(twogs_twobs <- dbinom(2, 4, 0.49, log = FALSE))

# Question # 2
# What is the probability that a family with 10 children has 4 girls and 6 boys (you can assume no twins)?

(fourgirl_sixboys <- dbinom(4,10, 0.49))

# Question # 3

# The genome has 3 billion bases. About 20% are C, 20% are G, 30% are T and 30% are A. 
# Suppose you take a random interval of 20 bases, what is the probability that the GC-content 
# (proportion of Gs or Cs) is strictly above 0.5 in this interval (you can assume independence)?

dbinom(1,20,0.4)
pbinom(10, 20, 0.4, lower.tail = FALSE, log.p = FALSE)

# Question # 4

# The following two questions are motivated by this event.

# The probability of winning the lottery is 1 in 175,223,510. If 189,000,000 randomly generated 
# (with replacement) tickets are sold, what is the probability that at least one winning tickets is sold?
# Give your answer as a proportion between 0 and 1, not a percentage.

pr_lottery <- 1/175223510
lambda <- pr_lottery*189000000
(1 - ppois(0,lambda,lower.tail = TRUE, log.p = FALSE))
1 - dbinom(0, 189000000, 1/175223510)

# Statistical Models Exercises #5
# Using the information from the previous question, what is the probability that two or more winning 
# tickets are sold?

1 - dbinom(0, 189000000, 1/175223510) -dbinom(1, 189000000, 1/175223510)

# Question # 6

# We can show that the binomial distribution is approximately normal with  𝑁  is large and  𝑝  is no
# t too close to 0 or 1. This means that

# (S(n) - E(S(n))/sqrt(var (s(n)))

# is approximately normal with mean 0 and SD 1. Using the results for sums of independent 
# random variables we learned in a previous course, we can show that  E(𝑆𝑁)=𝑁𝑝  and  Var(𝑆𝑁)=𝑁𝑝(1−𝑝) .

# The genome has 3 billion bases. About 20% are C, 20% are G, 30% are T and 30% are A. 
# Suppose you take a random interval of 20 bases, what is the exact probability that the 
# GC-content (proportion of Gs of Cs) is greater than 0.35 and smaller or equal to 0.45 in this interval?

(pb_between45 <- pbinom(9, 20, 0.4, lower.tail = TRUE, log.p = FALSE) - pbinom(7, 20, 0.4, lower.tail = TRUE, log.p = FALSE))

# For the question above, what is the normal approximation to the probability?

((pb_between45 * 20) - (20*0.4))/sqrt(20 * 0.4 * 0.6)
pnorm(-0.5527957,mean = 0,sd = 1,lower.tail = FALSE, log.p=FALSE)

b <- (9 - 20*.4)/sqrt(20*.4*.6)
a <- (7 - 20*.4)/sqrt(20*.4*.6)
pnorm(b)-pnorm(a)

# Statistical Models Exercises #8
# Repeat Statistical Models Exercises #3, but using an interval of 1000 bases.
# What is the difference (in absolute value) between the normal approximation and the exact probability 
# (using binomial) of the GC-content being greater than 0.35 and lesser or equal to 0.45?

(pb_between1000 <- pbinom(450, 1000, 0.4, lower.tail = TRUE, log.p = FALSE) - pbinom(350, 1000, 0.4, lower.tail = TRUE, log.p = FALSE))

b <- (450 - 1000*.4)/sqrt(1000*.4*.6)
a <- (350 - 1000*.4)/sqrt(1000*.4*.6)
normal_approximation <- pnorm(b)-pnorm(a)

abs(pb_between1000 - normal_approximation)

# The Cs in our genomes can be methylated or unmethylated. Suppose we have a large (millions) group of 
# cells in which a proportion  𝑝  of a C of interest are methylated. We break up the DNA of these cells
# and randomly select pieces and end up with  𝑁  pieces that contain the C we care about. This means 
# that the probability of seeing  𝑘  methylated Cs is binomial:

exact = dbinom(k,Ns,ps)

# We can approximate this with the normal distribution:

a <- (k+0.5 - Ns*ps)/sqrt(Ns*ps*(1-ps))
b <- (k-0.5 - Ns*ps)/sqrt(Ns*ps*(1-ps))
approx = pnorm(a) - pnorm(b)

# Let

Ns <- c(5,10,30,50, 100)
ps <- c(0.01,0.10,0.5,0.9,0.99)

# Question 9

# Compare the normal approximation and exact probability (from binomial) of the proportion of Cs being 
# 𝑘=1,…,𝑁−1 . Plot the exact versus approximate probability for each  𝑝  and  𝑁  combination

# Study the plots and tell us which of the following is NOT true.


mypar(2,2)
k <- c(1)
plot(exact, ps)
plot(exact,Ns)
plot(approx,ps)
plot(approx,Ns)

Ns <- c(5,10,30,100)
ps <- c(0.01,0.10,0.5,0.9,0.99)
library(rafalib)
mypar(4,5)
for(N in Ns){
  ks <- 1:(N-1)
  for(p in ps){
    exact = dbinom(ks,N,p)
    a = (ks+0.5 - N*p)/sqrt(N*p*(1-p))
    b = (ks-0.5 - N*p)/sqrt(N*p*(1-p))
    approx = pnorm(a) - pnorm(b)
    LIM <- range(c(approx,exact))
    plot(exact,approx,main=paste("N =",N," p = ",p),xlim=LIM,ylim=LIM,col=1,pch=16)
    abline(0,1)
  }
}

# Question 10

# We saw in the previous question that when  𝑝  is very small, the normal approximation breaks down.
# If  𝑁  is very large, then we can use the Poisson approximation.

# Earlier we computed the probability of 2 or more tickets winning the lottery when the odds of winning 
# were 1 in 175,223,510 and 189,000,000 tickets were sold. Using the binomial, we can run the code 
# below to compute the probability of exactly two people winning to be:

N <- 189000000
p <- 1/175223510
dbinom(2,N,p)

# If we were to use the normal approximation, we would overestimate this, as you can see by running this 
# code:

a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)

# To use the Poisson approximation here, use the rate  𝜆=𝑁𝑝  representing the number of tickets p
# er 189,000,000 that win the lottery. Run the following code and note how much better the approximation is:

dpois(2,N*p)

# In this case it is practically the same because  𝑁  is very very large and  𝑁𝑝  is not 0. These a
# re the assumptions needed for the Poisson to work.

# What is the Poisson approximation for the probability of two or more tickets winning?

1 - dpois(1,N*p) - dpois(0,N*p)

#----------------------------------------------------------------------------------------------------------

# Maximum Likelihood Exercises using the Human Cytomegalovirus Genome:

# In this assessment we are going to try to answer the question: 
# is there a section of the human cytomegalovirus genome in which the rate of palindromes is 
# higher than expected?
  
# Make sure you have the latest version of the dagdata library:

library(devtools)
install_github("genomicsclass/dagdata")

# and then load the palindrome data from the Human cytomegalovirus genome:
  
library(dagdata)
data(hcmv)  

# These are the locations of palindromes on the genome of this virus:
  
library(rafalib)
mypar()
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")

# These palindromes are quite rare,  𝑝  is very small. If we break the genome into bins of 4000 basepairs,
# then we have  𝑁𝑝  not so small and we might be able to use Poisson to model the number of palindromes 
# in each bin:

breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))

# So if our model is correct counts should follow a Poisson distribution. The distribution seems
# about right:

hist(counts)
  
# Let  𝑋1,…,𝑋𝑛  be the random variables representing counts then

# Pr(𝑋𝑖=𝑘)=𝜆𝑘/𝑘!exp(−𝜆) 
# So to fully describe this distribution we need  𝜆 . For this we will use MLE.

# To compute the Maximum Likelihood Estimate (MLE) we ask what is the probability of observing our data (which we denote with small caps) for a given  𝜆 :
  
# 𝐿(𝜆)=Pr(𝑋1=𝑥1 and 𝑋2=𝑥2 and …𝑋𝑛=𝑥𝑛;𝜆) 
# We assume that the  𝑋  are independent, thus the probabilities multiply:
  
# 𝐿(𝜆)=Pr(𝑋1=𝑥1)×Pr(𝑋2=𝑥2)×⋯×Pr(𝑋𝑛=𝑥𝑛) 
# Now we can write it in R. For example for  𝜆=4  we have:

probs <- dpois(counts,4)
likelihood <- prod(probs)
likelihood

# Run the code above to note that this is a tiny number. It is usually more convenient to compute 
# log-likelihoods

logprobs <- dpois(counts,4,log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood

# Now write a function that takes  𝜆  and the vector of counts as input, and returns the log-likelihood.
# Compute this log-likelihood for 

lambdas = seq(0,15,len=300)

# and make a plot.

loglikelihood_func <- function(L){
  logprobs <- dpois(counts,L,log=TRUE)
  (loglikelihood <- sum(logprobs))
}

values <- sapply(lambdas,loglikelihood_func)
(v <- plot(lambdas, values))

mle=optimize(loglikelihood_func,c(0,15),maximum=TRUE)
abline(v=mle$maximum)
mle$maximum

# Alternative code:

loglikelihood = function(lambda,x){
  sum(dpois(x,lambda,log=TRUE))
}
lambdas = seq(1,15,len=300)
l = sapply(lambdas,function(lambda) loglikelihood(lambda,counts))
plot(lambdas,l)
mle=lambdas[which.max(l)]
abline(v=mle)
print(mle)

# The point of collecting this dataset was to try to determine if there is a region of the genome that
# has higher palindrome rate than expected. We can create a plot and see the counts per location:

# Question 2: What is the center of the bin with the highest count?

breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab=)

binLocation[which.max(counts)]

# Question 3

# For the question above, what is the maximum count?

max(counts)

# Question 4

# Now that we have identified the location with the largest palindrome count, we want to know if by 
# chance we could see a value this big.

(lambda = mean(counts[ - which.max(counts) ]))

# If  𝑋  is a Poisson random variable with rate

1 - ppois(13,5)

# MLE Exercises #5

# From the question above, we obtain a p-value smaller than 0.001 for a count of 14. 
# Why is it problematic to report this p-value as strong evidence of a location that is different?

# We selected the highest region out of 57 and need to adjust for multiple testing.

# MLE Exercise # 6
# Use the Bonferroni correction to determine the p-value cut-off that guarantees a FWER of 0.05.

# What is this p-value cutoff?

0.05/57

# MLE Exercise # 7

# Create a qq-plot to see if our Poisson model is a good fit:

ps <- (seq(along=counts) - 0.5)/length(counts)
lambda <- mean( counts[ -which.max(counts)])
poisq <- qpois(ps,lambda)
qqplot(poisq,counts)
abline(0,1)

# Poisson is a very good approximation except for one point that we actually think is associated with a 
# region of interest.

#----------------------------------------------------------------------------------------------------------

# On Models for Variance:

# Install and load the following data library:
  
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)

# Now load this data and select the columns related to endometrium:
  
data("tissuesGeneExpression")
library(genefilter)
y = e[,which(tissue=="endometrium")]
dim(y)
# Question # 1

# Compute the across sample variance for the fifteen samples. 
# Then make a qq-plot to see if the variances follow a normal distribution.

# Which statement is true? (pick one)

mypar(3,5)
for (i in 1:15){ 
  qqnorm(y[,i])
  qqline(y[,i])     
}

# Rafa's Code

library(genefilter)
(s2 <- rowVars(y))
library(rafalib)
mypar(1,2)
qqnorm(s2)
qqline(s2)
##To see the square root transformation does not help much:
qqnorm(sqrt(s2))
qqline(sqrt(s2))

# Question 2
# Now fit an F-distribution with 14 degrees of freedom using the fitFDist() function in the limma package:
  
# What is estimated the estimated scale parameter?

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
?fitFDist()
estimates <- fitFDist(s2,14)

# Question 3
# Now create a qq-plot of the observed sample standard deviation versus the quantiles predicted by the 
# F-distribution (remember to take square root).

# Which of the following best describes the qq-plot?

theoretical<- sqrt(qf((seq(0,999)+0.5)/1000, 14, estimates$df2)*estimates$scale)
observed <- s2

mypar(1,2)
qqplot(theoretical,observed)
abline(0,1)

# Dr. Irizarry's Code:

ps <- (seq(along=s2)-0.5)/length(s2)
theoretical<- qf(ps,14,estimates$df2)*estimates$scale 
LIM <- sqrt( range(c(theoretical,s2)) )
mypar(1,2)
qqplot(sqrt( theoretical ), sqrt( s2 ),ylim=LIM,xlim=LIM)
abline(0,1)
##close up excluding the upper 5%
K <- sqrt( quantile(s2,0.95) )
qqplot( sqrt( theoretical ), sqrt( s2 ),ylim=c(0,K),xlim=c(0,K))
abline(0,1)

#----------------------------------------------------------------------------------------------------------

# Mammograms, RNA seq Data Statistics, & Monte Carlo Simulations in F Distributions

# Mammograms are important breast cancer screening tests that have contributed to an increase in early 
# diagnosis and decrease in breast cancer mortality. However, like most screening tests, they have a high 
# false positive rate - most women with an abnormal mammogram who are called back for additional testing 
# (like breast biopsies) do not have breast cancer. 
# The probability that a woman with a positive mammogram has cancer on follow-up testing is around 0.1.

# Suppose you are a pathologist evaluating 30 random breast biopsies from mammogram follow-up tests.

# A. What is the probability that none of the biopsies show cancer?
  
1 - pbinom(9, 30, 0.1)

# B. What is the probability that exactly 3 of the biopsies show cancer?

dbinom(3, 30, 0.1)

# C. What is the probability that at least 10 of the biopsies show cancer?

1 - pbinom(9, 30, 0.1)

# Question 2

# Suppose you are analyzing RNA-seq data and transcript X is expressed at a level such that it represents 
# 2 out of every 1,000,000 transcripts. This means the probability of observing transcript X in a random 
# read is 0.000002. Now suppose that you evaluate 3,000,000 reads in an experiment.

# A. What is the expected number of reads for transcript X?

(lambda <- (2/1000000)*3000000)

# B. Use the Poisson distribution to calculate the probability of observing exactly 1 read for transcript X.

ppois(1,6)

# C. What is the probability of observing more than 10 reads for transcript X?
  
1 - ppois(10,6)

# Question 3

# In the human genome, cytosines that are followed by guanines (CpGs) are methylated 80% of the time.
# A. Consider 30 CpG sites. Using the binomial distribution, what is the exact probability that between 70% and 90% of the CpGs are methylated?

.70*30
.9*30
pbinom(27,30,0.8) -pbinom(21,30,0.8)

# B. Using the normal distribution, what is the approximate probability that between 70% and 90% of 
# CpGs are methylated?

N <- 30
p <- 0.8
a <- (27 - N*p)/sqrt(N*p*(1-p))
b <- (21 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)

# C. What is the difference (in absolute value) between the normal approximation and the exact 
# probability (using binomial) of observing methylation between 70% and 90%?

abs(0.8290965 - 0.8271703)

# Question 4

# In a previous week, we performed 1000 simulations of a series of 30 mouse experiments under the 
# null distribution and, for each simulation, counted the number of p-values under 0.05 to generate a
# vector pval_counts:

set.seed(100)
pval_counts = replicate(1000,{
  pvals = replicate(30, {
    cases = rnorm(5,7.5,2.5)
    controls = rnorm(5,7.5,2.5)
    t.test(cases,controls)$p.value
  })
  sum(pvals < 0.05)
})

mean(pval_counts)

# This random sampling can be modeled as a Poisson process, and the Maximum Likelihood Estimate can 
# be used to determine the  𝜆  that best fits this process.

# This function takes a  𝜆  and a vector of counts as inputs and returns the log-likelihood for that  𝜆 :

loglikelihood = function(lambda,x){
  sum(dpois(x,lambda,log=TRUE))
}

# Compute this log-likelihood for:

# A. Which value of  𝜆  maximizes the log likelihood?

lambdas = seq(0,10,len=101)
l = sapply(lambdas,function(lambdas) loglikelihood(lambdas,pval_counts))
plot(lambdas,l)
mle=lambdas[which.max(l)]
abline(v=mle)
max(mle)

# B. Given that value of lambda, what is the probability of observing 3 or more p-values below 0.05.

1 - ppois(3,1.3)

# C. Compare the estimated value of  𝜆  from the simulated experiment to the theoretical expected value
# of  𝜆 . How many p-values are expected to be below 0.05 due to chance given  𝑁=30  tests with 
# a probability of success of  p = 0.05 ?

lambda <- 0.05*30

# Question 5

# You can generate a set of random variables from an F-distribution with the function rf(). 
# This line of code generates 100 random numbers from an F-distribution with parameters  


x = rf(100,df1=8,df2=16)

# Set the seed to 25, then generate an F-distributed list of random numbers x using the code above. 
# Use fitFDist() from the limma package to fit an F-distribution to x using df1 = 8.

set.seed(25)
library(limma)
(F <- fitFDist(x,8))
F$df2

# Question 6

# Set the seed to 28, then use replicate() to repeat the previous procedure 1000 times: 
# each time, generate 100 F-distributed random numbers with the code provided, then use fitFDist() 
# with a known value of df1=8and determine the estimated value of df2.

# A. What is the median value of df2 in this Monte Carlo simulation?

set.seed(28)
F_values_1 <- replicate(1000,{
  x = rf(100,df1=8,df2=16)
  values <- fitFDist(x,8)
  values$df2
})

median(F_values)

# B. What proportion of estimated df2 values are between 12 and 20 in this Monte Carlo simulation?


sum(F_values_1 > 12 & F_values_1 < 20)/1000

# Question 7 

# Set the seed to 28 again, then repeat the previous question except this time increase the number 
# of randomly generated values in rf() to 1000, representing a larger sample size. 
# Again, use fitFDist() with a known value of df1=8 and determine the estimated value of df2.

# A. What is the median value of df2 in the Monte Carlo simulation with a larger sample size?

set.seed(28)
F_values_2 <- replicate(1000,{
  x = rf(1000,df1=8,df2=16)
  values <- fitFDist(x,8)
  values$df2
})

median(F_values)

# B. What proportion of estimated df2 values are between 12 and 20 in the Monte Carlo simulation 
# with a larger sample size?

sum(F_values_2 > 12 & F_values_2 < 20)/1000

# Question 8

boxplot(F_values_1)
boxplot(F_values_2)

#---------------------------------------------------------------------------------------------------------

# On Hierarchal Models, Limma and Introduction tp Volcano Plots:

# Question 1

p <- 1/4000

(0.99*p)/(0.99*p+(0.01*(1-p)))

# Question 2

tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
##this shows us files
filenames <- unzip(tmpfile,list=TRUE)
players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

# Which of the following dplyr commands gives us the batting averages (AVG) for players with more 
# than 500 at bats (AB) in 2012:@

filter(players,yearID==2012) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)

# Question 3 

# Edit the command above to obtain all the batting averages from 2010, 2011, 2012 and
# removing rows with AB < 500.

# What is the average of these batting averages?

players_list <- filter(players,yearID>= 2010, yearID <= 2012) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)

mean_prior <- mean(players_list$AVG)

# What is the standard deviation of these batting averages?

sd_prior <- sd(players_list$AVG)

# Use exploratory data analysis to decide which of the following distributions approximates the 
# distribution of the average across players (hint: this is contained in the AVG component)?

library(rafalib)
mypar(1,2)
hist(players_list$AVG)
qqnorm(players_list$AVG)
qqline(players_list$AVG)

# It is April and after 20 at bats, Jose Iglesias is batting .450 (this is very good). We 
# can think of this as a binomial distribution with 20 trials with probability of success  𝑝 .
# Our sample estimate of  𝑝  is .450. What is our estimate of standard deviation?
# Hint: This AVG a sum of Bernoulli trials, that is binomial, divided by 20.

p <- 0.450
n <- 20
sd_sample <- sqrt(((p)*(1-p))/n)

# The sum (numerator of AVG) is binomial so it has SD sqrt(Np(1-p)) . 
# The SD of a random variable times a constant is the SD of the random variable times that constant. 
# For the AVG we divide by  𝑁  to get  𝑝(1−𝑝)/𝑁‾‾‾‾‾‾‾‾‾‾‾√ . This is

# The Binomial is approximated by normal when the sample size is large, so our sampling distribution 
# is approximately normal with mean  𝜃  = 0.45 and SD  𝜎=0.11 . Earlier we used a baseball databas
# e to determine that our prior distribution for  𝜃  is Normal with mean  𝜇=0.275  and SD  𝜏=0.027 
# We saw that this is the posterior mean prediction of the batting average.

# What is your estimate of Jose Iglesias' batting average going forward taking into account 
# his current batting average?


B <- sd_sample^2/(sd_prior^2 + sd_sample^2)
(E <- mean_prior + (1 - B)*(0.450 - mean_prior))

# Load the following data (you can install it from Bioconductor) and extract the data matrix 
# using exprs() (we will discuss this function in detail in a future course):

BiocManager::install("SpikeInSubset")
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)

# This dataset comes from an experiment in which RNA was obtained from the same background pool to 
# create six replicate samples. Then RNA from 16 genes were artificially added in different quantities 
# to each sample. These quantities (in picoMolars) and gene IDs are stored here:

pData(rma95)
y <- exprs(rma95)

# Note that these quantities were the same in the first three arrays and in the last three arrays. 
# So we define two groups like this:

g <- factor(rep(0:1,each=3))

# and create an index of which rows are associated with the artificially added genes:
  
spike <- rownames(y) %in% colnames(pData(rma95))

# Note that only these 16 genes are differentially expressed since these six samples differ only due 
# to random sampling (they all come from the same background pool of RNA).

# Perform a t-test on each gene using the rowttests() function in the genefilter package.

# What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?

library(genefilter)
rtt = rowttests(y,g)
index = rtt$p.value < 0.01 
print (mean( !spike[index] ))

## We can make a volcano plot to visualize this:
mask <- with(rtt, abs(dm) < .2 & p.value < .01)
cols <- ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))
with(rtt,plot(-dm, -log10(p.value), cex=.8, pch=16,
              xlim=c(-1,1), ylim=c(0,5),
              xlab="difference in means",
              col=cols))
abline(h=2,v=c(-.2,.2), lty=2)

# Now compute the within group sample standard deviation for each gene (you can use group 1). 
# Based on the p-value < 0.01 cut-off, split the genes into true positives, false positives, true negatives and false negatives. 
# Create a boxplot comparing the sample SDs for each group. 
# Which of the following best described the box-plot?
  
library(genefilter)
sds <- rowSds(y[,g==0])
index <- paste0( as.numeric(spike), as.numeric(rtt$p.value<0.01))
index <- factor(index,levels=c("11","01","00","10"),labels=c("TP","FP","TN","FN"))
boxplot(split(sds,index)) 

# Question 3

# In the previous two questions we observed results consistent with the fact that the random 
# variability associated with the sample standard deviation leads to t-statistics that are large by chance.

# Note that the sample standard deviation we use in the t-test is an estimate and that with just a 
# pair of triplicate samples, the variability associated with the denominator in the t-test can be large.

# The following three steps perform the basic limma analysis. 
# The eBayes step uses a hierarchical model that provides a new estimate of the gene specific standard error.

library(limma)
fit <- lmFit(y, design=model.matrix(~ g))
colnames(coef(fit))
fit <- eBayes(fit)

# Make a plot of the original new hierarchical models based estimate versus the sample based estimate.

sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)

LIM = range( c(posteriorSD,sampleSD))
plot(sampleSD, posteriorSD,ylim=LIM,xlim=LIM)
abline(0,1)
abline(v=sqrt(fit$s2.prior))

# Use these new estimates (computed in Question 4.6.3) of standard deviation in the denominator 
# of the t-test and compute p-values. You can do it like this:

library(limma)
fit = lmFit(y, design=model.matrix(~ g))
fit = eBayes(fit)
##second coefficient relates to diffences between group
pvals = fit$p.value[,2] 
index <- pvals < 0.01

# What proportion of genes with a p-value < 0.01 (no multiple comparison correction) 
# are not part of the artificially added (false positives)?
index = rtt$p.value < 0.01 
print (mean( !spike[index] ))
  
# Compare to the previous volcano plot and notice that we no longer have small p-values for genes with 
# small effect sizes.

#---------------------------------------------------------------------------------------------------------

# Explorative Data Analysis: Log Ratio Comparisons

# Download and install the Bioconductor package SpikeInSubset and then load the library and the mas133 data:
  
data(mas133)  
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

# What proportion of the points are inside the box?
  
(e_1 <- as.numeric(e[,1]))
(e_2 <- as.numeric(e[,2]))
(sum(e_1 <= 3000) + sum(e_2 <= 3000))/44600

#Now make the sample plot with log:

plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

# When you take the log, 95% of data is no longer in a tiny section of plot.

# The two samples we are plotting are replicates (they random samples from the same batch of RNA). 
# The correlation of the data was 0.997 in the original scale, 0.96 in the log-scale. 
# High correlations are sometimes confused for evidence of replication. But replication implies 
# we get very small difference between the observations, which is better measured with distance or 
# differences.

# What is the standard deviation of the log ratios for this comparison?
# Make an MA-plot:
  
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5) 
sd(e[,2] - e[,1])

# How many fold changes above 2 do we see? Note that these measures of log (base 2) of expression so a fold
# change of 2 translates into a difference, in absolute value, of 1.

sum(abs(e[,2] - e[,1]) >= 1)

# ------------------------------------------------------------------------------------------------

# Hierarchal Models

# The incidence of prostate cancer in men over the age of 50 is roughly 0.5%. 
# A prostate cancer screening test (PSA) exists, but it has recently fallen out of favor for several reasons. 
# The PSA test is positive for 51% of men with advanced prostate cancer and negative for 91% of men without
# prostate cancer. These probabilities are summarized below:

# Question # 1

# A. 

(p_neg <- (0.49*0.005)+( 0.91*0.995))
(p_pc_neg <- (0.49*0.005)/p_neg)
(p_neg_pc <- (p_pc_neg*p_neg)/0.005)

# B. 

(p_pos <- 0.51*0.005 + 0.09*0.995)
(p_nopc_pos <- (0.49*0.995)/p_pos)
(p_pos_nopc <- (p_pos*p_nopc_pos)/0.995)

# C.

(0.51*0.005)/p_pos

# The ChickWeight dataset included in base R contains weights of chicks on different diets over the first 
# 21 days of life.

# Question 2

# Suppose we want to evaluate the weight of a chick at day 21. Filter the ChickWeight data to 
# include only weights at day 21:
 
library(tidyverse)
day21 <- ChickWeight %>% filter(Time == 21)    # you don't need to load any packages to access ChickWeight

# A. What is the mean weight of chicks on diet 3 at day 21?
# Hint this will be the u

mu <- mean(day21$weight)

# B. What is the standard deviation of chick weight on diet 3 at day 21?
# Hint this will be your tau  

tau <- sd(day21$weight)

# Questio 3

# In general, it is fairly uncommon for a 21 day old chick to weigh over 300g. 
# However, different diets affect the chick weights. Suppose the chick is on diet 3.

day21_B <- ChickWeight %>% filter(Time == 21 & Diet == 3) 

# A. What is the mean weight of chicks on diet 3 at day 21?
  
Y <- mean(day21_B$weight)

# B. What is the standard deviation of chick weight on diet 3 at day 21?

s <- sd(day21_B$weight)

# C. Assume that chick weights on diet 3 follow a normal distribution. 
# What is the probability that a 21 day old chick on diet 3 weighs 300g or more?

1 - pnorm(300, 270.3,71.62254)

# Question 4

# Chicks on diet 3 have a higher probability of weighing over 300g than the general population of chicks 
# on all diets. However, note that we have less information about chicks on each individual diet than 
# we do about chicks on all diets - there are only 10 weights for chicks on diet 3. 
# This means it may be helpful to apply a hierarchical model to chick weights based on diet.

# A. Using a hierarchical model that combines the overall weight data with the diet 3 weight data, what is the expected weight of a chick on diet 3 at day 21?
  
B <- s^2/(s^2 + tau^2)
(E <- mu + (1 - B)*(Y - mu))

# B.  Using a hierarchical model that combines the overall weight data with the diet 3 weight data, 
# what is the standard error of chick weights on diet 3 at day 21?
  
(sterror <- sqrt(1 /(1/s^2 + 1/tau^2)))

# C. Given the expected value and standard error of this hierarchical model, and assuming a normal 
# distribution, what is the probability that a 21 day old chick on diet 3 weighs over 300g?

1 - pnorm(300,E,sterror)

# The probability has reduce when taking into account the prior statistics.

# Question 5

# Suppose that you use rowttests() from the genefilter library to compare gene expression across two 
# biological contexts and you assign the resulting output to results. 
# You want to create a volcano plot of the results.

# Which of these options gives the values for the x-axis of the volcano plot?

# results$dm

# Which of these options gives the values for the y-axis of the volcano plot?
  
# -log10(results$p.value)

# Question 6 - 8

# In previous exercises, we analyzed an RNA-seq experiment from Bottomly et al. comparing gene expression 
# across two different strains of mice:

if (!file.exists("bottomly_eset.RData")) download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData", 
                                                       "bottomly_eset.RData")
load("bottomly_eset.RData")

# also load these libraries, which we previously installed from Bioconductor
library(Biobase)
library(genefilter)
library(qvalue)

dat = exprs(bottomly.eset)    # gene expression matrix
strain = pData(bottomly.eset)$strain    # strain factor

results <- rowttests(dat,strain)
pvals <- results$p.value

# Set the seed to 1, then permute the strain information:
  
set.seed(1)
permut <- sample(strain)
results_1 <- rowttests(dat,permut)
pvals_1 <- results_1$p.value

# Question 6
# How many genes have a p-value below .05 in this simulated null distribution?
sum(pvals_1 < 0.05, na.rm = TRUE)

# Question 7

# Create a histogram of p-values for both the original results, pvals, and the permuted p-values.
# Which of the following is NOT true about the distribution of p values?
  
library(rafalib)
mypar(1,2)
hist(pvals)
hist(pvals_1)

# Because the permuted p-values do not follow a uniform distribution, this suggests unexpected 
# correlation between some samples.

# Question 8

# Samples 1 and 4 are both from the mouse strain C57BL/6J. If these biological replicates are 
# highly correlated to each other, then an MA-plot of these samples should be symmetrical around the 
# line x=y (equivalent to x - y = 0).

# Assign x and y as the log base 2 of samples 1 and 4 respectively. Note that adding 1 before taking the 
# log prevents problems related to zeros in the data:
  
x <- log2(dat[,1] + 1)
y <- log2(dat[,4] + 1)

plot(((x+y)/2),x-y)
abline(0,0)
