####################################################################################################################

# By: Santiago Taguado Menza
# High-Dimensional Data Analysis Section 3: Hierarchical Clustering, Heat Maps, & k-Nearest Neighbour Classification
# March 30th, 2021
# HarvardX

####################################################################################################################

# Hierarchical Clustering & Dendrograms

# Create a random matrix with no correlation in the following way:
set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n

# Run hierarchical clustering on this data with the hclust() function with default parameters to cluster 
# the columns. Create a dendrogram.

# From the dendrogram which pairs of samples are the furthest away from each other?

set.seed(1)
d <- dist(t(x))
hc <-hclust(d)
length(hc)
sd(cutree(hc,h = 143))/sqrt(24)
plot(hc)

# Samples 17 and 9 are far away from each other
# Set the seed at 1 with set.seed(1) and replicate the creation of this matrix 100 times:

set.seed(1)
xs <- replicate(100,{
  m = 10000
  n = 24
  x = matrix(rnorm(m*n),m,n)
  d1 <- dist(t(x))
  hc_d1 <- hclust(d1)
  return(sd(cutree(hc_d1, h = 143))/7)})

set.seed(1)
m = 10000
n = 24
nc = replicate(100,{
  x = matrix(rnorm(m*n),m,n)
  hc = hclust( dist( t(x)))
  length(unique(cutree(hc,h=143)))
})

plot(table(nc)) ## look at the distribution
popsd(nc)

####################################################################################################################

# K-means
# Run kmeans() with 5 centers for the blood RNA data:

library(GSE5859Subset)
data(GSE5859Subset)

# Set the seed to 10, set.seed(10), right before running kmeans() with 5 centers.
# Explore the relationship of clusters and information in sampleInfo. 
# Which of the following best describes what you find:

set.seed(10)
mypar(3,1)
km <- kmeans(t(geneExpression), centers=5)
plot(sampleInfo$ethnicity,km$cluster)
plot(sampleInfo$date,km$cluster)
plot(sampleInfo$group,km$cluster)

mds=cmdscale(dist(t(geneExpression)))
set.seed(10)
result=kmeans(t(geneExpression),5)
library(rafalib)
mypar(1,1)
plot(mds,bg=result$cl,pch=21)
table(sampleInfo$group,result$cluster)
table(sampleInfo$date,result$cluster)

##looks better if we re-order:
table(sampleInfo$date,result$cluster)[,c(4,1,5,3,2)]

# Date is driving the clusters

####################################################################################################################

# Heat Map Exercises

library(GSE5859Subset)
data(GSE5859Subset)

# Pick the 25 genes with the highest across sample variance. This function might help

install.packages("matrixStats")
library(matrixStats)
##we use mads due to a outlier sample

# While a heatmap function is included in R, we recommend the heatmap.2 function from the gplots package 
# on CRAN because it is a bit more customized. For example, it stretches to fill the window.

library(gplots)
rowMads(geneExpression)

# Use heatmap.2() to make a heatmap showing the sampleInfo$group with color, the date as labels, 
# the rows labelled with chromosome, and scaling the rows.
# What do we learn from this heatmap?

library(RColorBrewer)
library(genefilter)
library(rafalib)

# Identifying the top 25 most variable genes via rowVars & RowMads
rv <- rowVars(geneExpression)
(idx <- order(-rv)[1:25])
rmads <- rowMads(geneExpression)
(idx_1 <- order(-rmads)[1:25])

# Making the colors
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]

# Making Columns
labcol= gsub("2005-","",sampleInfo$date) 

# Making the heatmap
heatmap.2(geneExpression[idx_1,], 
          col = cols,
          trace = "none",
          scale ="row",
          labRow = geneAnnotation$CHR[idx_1],
          labCol = labcol,
          ColSideColors=gcol,
          key=TRUE)

# Second Type of Heat Map

##make colors
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]

##make lables: remove 2005 since it's common to all
labcol= gsub("2005-","",sampleInfo$date)  

##pick highly variable genes:
sds =rowMads(geneExpression)
(ind = order(sds,decreasing=TRUE)[1:25])

## make heatmap
heatmap.2(geneExpression[ind,],
          col=cols,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[ind],
          labCol=labcol,
          ColSideColors=gcol,
          key=FALSE)

# A group of chrY genes are higher in group 0 and appear to drive the clustering. Within those clusters 
# there appears to be clustering by month.

# Create a large data set of random data that is completely independent of sampleInfo$group like this:

set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$g )
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)

# Create two heatmaps with these data. Show the group g either with labels or colors.

# 1. Taking the 50 genes with smallest p-values obtained with rowttests

# 2. Taking the 50 genes with largest standard deviations.

(pvals <- rowttests(x,g))
(idx_2 <- order(pvals$p.value,decreasing = FALSE)[1:50]) 

sds_1 = rowMads(x)
(idx_3 = order(sds,decreasing=TRUE)[1:50])

# Heatmap # 1
heatmap.2(x[idx_2,],
          col=cols,
          trace="none",
          scale="row",
          labCol = g)

# Heatmap # 2
heatmap.2(x[idx_3,],
          col=cols,
          trace="none",
          scale="row",
          labCol = g,
          key= FALSE)

# The first heatmap shows a relationships between g and x, but with so many test some will appear significants.
# Selecting genes with t-test gives us a deceiving results.

# How to vectorize the two heatmaps above: 

ttest = rowttests(x,g)
sds = rowSds(x)
Indexes = list(t=order(ttest$p.value)[1:50], s=order(-sds)[1:50])
Indexes
for(ind in Indexes){
  heatmap.2(x[ind,],
            col=cols,
            trace="none",
            scale="row",
            labCol=g,
            key=FALSE)
}

####################################################################################################################

# Conditional Expectations

# Throughout this assessment it will be useful to remember that when our data are 0s and 1s, 
# probabilities and expectations are the same thing. We can do the math, but here is an example 
# in the form of R code:

n = 1000
y = rbinom(n,1,0.25)
##proportion of ones Pr(Y)
sum(y==1)/length(y)
##expectaion of Y
mean(y)

# Generate some random data to imitate heights for men (0) and women (1):

n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]

# Treating the data generated above as the population, if we know someone is 176 cm tall, 
# what it the probability that this person is a woman:  Pr(𝑌=1|𝑋=176)=E(𝑌|𝑋=176) ?

(sum(y[x== 176]==1))/ (sum(y[x==176]==1) + sum(y[x==176]==0))
mean(y[x==176])

# Conditional Expectation Exercises #2

# Now make a plot of  𝐸(𝑌|𝑋=𝑥)  for x=seq(160,178) using the data generated in Conditional Expectat
# ion Exercises #1.Suppose for each height  𝑥  you predict 1 (female) if  Pr(𝑌|𝑋=𝑥)>0.5  and 0 (male) otherwi
# se. What is the largest height for which you predict female ?

(x_vals = seq(160,178))

mean_func <- function(k){
  mean(y[x == k])
}

(y_vals <- sapply(x_vals,mean_func))
plot(x_vals,y_vals, xlab = "Height", ylab = "Probability of being a Women")
abline(h=0.5)
abline(v=168)

####################################################################################################################

# Bin Smoothing with Local Weighted Regression

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Use the data generated in a previous question about men's and women's heights:

RNGkind(sample.kind = "Rounding")

n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]

# Set the seed at 5, set.seed(5), and take a random sample of 250 individuals from the population like this:

set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]

fit <- loess(Y ~ X)
predict(fit,168)
predict(fit,168,se = TRUE)$se.fit

# The loess estimate above is a random variable thus we should compute its standard error. 
# Use Monte Carlo simulation to compute the standard error of your estimate of  𝑓(168) .

# Set the seed to 5, set.seed(5), and perform 1000 simulations of the computations performed in question 
# 2.7.1. Report the the SE of the loess based estimate.

library(rafalib)
set.seed(5)
vals <- replicate(1000, {
  N = 250
  ind = sample(length(y),N)
  Y = y[ind]
  X = x[ind]
  fit <- loess(Y ~ X)
  preds <- predict(fit,168)
  return(preds)})

popsd(vals)

####################################################################################################################

# K-nearest neighbors

# Changes in R since the creation of this material have altered the randomization code. 
# You will need to include the following line in your code before you call set.seed(N) in order 
# to obtain the correct answers:

RNGkind(sample.kind = "Rounding")

library(GSE5859Subset)
data(GSE5859Subset)

# And define the outcome and predictors. To make the problem more difficult, we will only consider 
# autosomal genes:

(y = factor(sampleInfo$group))
X = t(geneExpression)
out = which(geneAnnotation$CHR%in%c("chrX","chrY"))
X = X[,-out]

# Note, you will also need to load the following package:

library(caret)

# kNN and Cross Validation Exercises #1
# Set the seed to 1, set.seed(1), then use the createFolds() function in the caret package to 
# create 10 folds of y.

# What is the 2nd entry in the fold 3?

set.seed(1)
(idx <- createFolds(y, k = 10))
sapply(idx, function(i) table(y[i]))

# For the following questions we are going to use kNN. 
# We are going to consider a smaller set of predictors by filtering genes using t-tests. 
# Specifically, we will perform a t-test and select the  𝑚  genes with the smallest p-values.

# Let  𝑚=8  and  𝑘=5  and train kNN by leaving out the second fold, idx[[2]].

# How many mistakes do we make on the test set? 
# Remember it is indispensable that you perform the ttest on the training data.

library(genefilter)
library(class)
m = 8
k = 5

# Perform P Test
ind = idx[[2]]
(pvals <- rowttests(t(X[-ind,]),factor(y[-ind])))
(idx_2 <- order(pvals$p.value,decreasing = FALSE)[1:8])

# Run Predicted Algorithm
i = 2
k = 5
predict=knn(X[-ind,idx_2],X[ind,idx_2],y[-ind],k=k)
sum(predict != y[ind])

# Now run the code for kNN and Cross Validation Exercises #2 for all 10 folds and keep track of the errors. 
# What is our error rate (number of errors divided by number of predictions) ?

m=8
k=5

ks<- 1:10
loop_func<- function(k){
  ind = idx[[k]]
  (pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val)
  ind2 = order(pvals)[1:8]
  predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=5)
  sum(predict!=y[ind])
}
sum(sapply(ks,loop_func))/length(y)

# Now we are going to select the best values of  𝑘  and  𝑚 
# Use the expand.grid() function to try out the following values:

ms=2^c(1:11)
ks=seq(1,11,1)
(params = expand.grid(k=ks,m=ms))
params
# Now use sapply() or a for loop to obtain error rates for each of these pairs of parameters. 
# Which pair of parameters minimizes the error rate?

errors_1 = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})

params[which.min(errors_1),]

##make a plot and confirm its just one min:
errors = matrix(errors_1,11,11)
library(rafalib)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks),bty = "n")

# Repeat question kNN and Cross Validation Exercises #4 but now perform the t-test filtering before the cross validation.
# Note how this biases the entire result and gives us much lower estimated error rates.

# What is the minimum error rate?

pvals <- rowttests(t(X),factor(y))$p.val
errors_2 = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})

params[which.min(errors),]

(errors = matrix(errors_2,11,11))
library(rafalib)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topleft",as.character(ks),lty=seq_along(ks),col=seq_along(ks),bty = "n")

# Repeat the cross-validation we performed in question kNN and Cross Validation Exercises #4, 
# but now instead of defining y as sampleInfo$group, use:

y = factor(as.numeric(format( sampleInfo$date, "%m")=="06"))
k = 1:10  
ms=2^c(1:11)
ks=seq(1,11,1)
params = expand.grid(k=ks,m=ms)

errors_3 = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})

params[which.min(errors_3),]

min(errors_3)

(errors = matrix(errors_3,11,11))
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks),bty = "n")

####################################################################################################################

# Quiz 3

# The heights dataset from the dslabs package (available from CRAN) contains self-reported heights (in inches) 
# for male and female students from three Harvard Biostatistics classes:

library(dslabs)
data(heights)

# For simplicity, round heights to the nearest inch:
# Treat this as a population for all.

(heights$height <- round(heights$height))

# Calculate the conditional probability that a person 67 inches tall is female.

n = 812
m = 238

men = heights[which(heights$sex == "Male"),]
women = heights[which(heights$sex == "Female"),]
(y = c(rep(0,n),rep(1,m)))
x = round(c(men$height,women$height))

# Conditional Probability Formula
mean(y[x==67])

# Calculate the conditional probability that a person is female for the vector of heights hts = 60:80. 
# Make a plot of this conditional probability versus hts. 
# Suppose you predict female for any height for which the conditional probability of being female  
# E(𝑌=Female|𝑋=𝑥)  is > 0.5. What is the maximum height for which you predict a person is female?

hts = seq(60,80,1)
conditional_prob = function(hts){
  mean(y[x == hts])
}
probs = sapply(hts,conditional_prob)
plot(hts,probs)  
abline(h = 0.5, v = 64)  

# The leukemiasEset contains 60 sets of bone marrow gene expression data from patients with one 
# of the 4 main types of leukemia (ALL, AML, CLL, CML) as well as control patients without 
# leukemia (NoL).

# Install and load the leukemiasEset data from the leukemiasEset Bioconductor package:

BiocManager::install("leukemiasEset")    # install if needed
library(leukemiasEset)
data(leukemiasEset)

# These data are stored in a container called an ExpressionSet. In future courses, we will learn how 
# to work with ExpressionSets directly, but for now we can extract gene expression data as a 
# matrix dat (features are rows, columns are samples):
  
dat = exprs(leukemiasEset)

# We can also create a vector noting which type of leukemia is present in each sample:

leuk = leukemiasEset$LeukemiaType

# A. How many features are present in dat?

dim(dat)

# B. How many samples are present in dat?
 
dim(dat)

# C. How many samples are from patients with AML?
  
sum(leuk == "AML")

# Make an MDS plot of dat and color the points by leuk.
# Which of the following are TRUE?

mds=cmdscale(dist(t(dat)))
plot(mds,col = leuk)
legend("bottomright",levels(leuk),col=seq_along(leuk),pch = 1,bty = "n")

# Observations
# CLL samples tend to have higher values of mds[,2] than ALL
# The samples with the highest values of mds[,1] are all CML
# At a glance, CML samples are more similar to NoL samples than to other leukemias.

# Run hierarchical clustering on this data with the hclust() function with default parameters to 
# cluster the columns. Create a dendrogram and use the leukemia type leuk as labels.

# Suppose you want to cut the tree so that there are 5 clusters. 
# Which of these heights would be the best cutoff?

d = hclust(dist(t(dat)))
sd(cutree(d,h = 150))/sqrt(24)
plot(d)
leuk_type = cutree(d, h = 150)
plot(leuk_type)

# Using the cutoff height that generates 5 clusters in the previous problem, one cluster contains 
# exactly 12 samples that are all from the same leukemia type.
# Which two leukemia types have all samples of that type in a unique cluster?

mypar(1,1)  
plot(d, labels = leuk, lab.col= leuk, cex=0.5)
abline(h = 150)

set.seed(4)
hclusters <- cutree(d, h = 150)
table(true=leuk, cluster=hclusters)

# Using Kmeans

set.seed(4)
result=kmeans(t(dat),5)
table(result$cluster, leuk)

# Pick the 25 genes with the highest across sample variance using the rowMads() function from 
# matrixStats:
  
library(matrixStats)
sds =rowMads(dat)
ind = order(sds,decreasing=TRUE)[1:25]

# Use heatmap.2() from gplots to make a heatmap showing the leuk type with column colors as well as
# column labels, and scaling the rows. (In the future, we will learn how to convert gene IDs, 
# like "ENSG000…", into gene names.)

# Which of the following statements are TRUE about the heatmap?

library(RColorBrewer)
library(gplots)
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
gcol=brewer.pal(5,"Dark2")
gcol=gcol[as.numeric(leuk)]

## make heatmap
heatmap.2(dat[ind,],
          col=cols,
          trace="none",
          scale="row",
          labCol=leuk,
          ColSideColors=gcol)

# Observations:
# Over 20 of the genes with the highest across sample variance are upregulated in CML and NoL and 
# downregulated in other leukemias.

# The bottom 2 genes in the plot tend to be upregulated in ALL and CLL and downregulated in AML and CML.

# Based on these 25 genes, the type of leukemia with the closest expression pattern to normal (NoL) 
# bone marrow is CML.

# Suppose you want to design an algorithm that can predict whether a sample from the leukemia dataset 
# is normal ("NoL") versus any type of leukemia.

# Start by creating a vector leukTF that is TRUE when a sample is normal and FALSE when a 
# sample is leukemia:

leukTF = leuk == "NoL"

# Load the caret library and set the seed to 2. Use createFolds() on leuk2 to create 5 folds for 
# cross-validation. Save the indices for these folds as idx.

# Before running any machine learning algorithms on these folds, it is best to ensure that each fold 
# contains both normal and leukemia samples. Count the number of normal samples in each fold.
factor(leukTF)

library(caret)
set.seed(2)

(idx = createFolds(leukTF, k=5))

normal_counts = sapply(1:length(idx), function(x){
  fold_ind = idx[[x]]
  sum(leukTF[fold_ind]==TRUE)
})

sum(normal_counts > 0)
sum(normal_counts == 3)

# We are going to consider a smaller set of predictors by filtering genes using t-tests. 
# Specifically, we will perform a t-test and select the  𝑚  genes with the smallest p-values.

# Let  𝑚=3 . Leave out the first fold, idx[[1]], and perform rowttests() from the genefilter
# library on the remaining samples. Find the row numbers of the 3 genes with the lowest 
# p-values and save these as gene_ind.

# Which of these rows does not represent one of the three genes with the lowest p-values when 
# omitting the first fold, stored in gene_ind?

library(genefilter)
set.seed(1)
m = 3
k=5
(ind = idx[[1]])
X = t(dat)
pvals = rowttests(dat[,-ind],factor(leukTF[-ind]))$p.val
(gene_ind = order(pvals)[1:3])
predict=knn(X[-ind,gene_ind],X[ind,gene_ind],leukTF[-ind],k=5)
sum(predict!=leukTF[ind])

# Separate dat into a test set consisting of samples in the first fold and a training set consisting of
# samples in all other folds. Keep only genes from gene_ind in these sets. 
# (Your test set should be an 11x3 matrix and your training set should be a 49x3 matrix.)

# Train a kNN model and generate predictions for the test set using the knn function from the class 
# library and k=5.

# How many errors does this model make on the test set for the first fold?

library(caret)
library(genefilter)
library(class)

set.seed(4)
idx = createFolds(leukTF, k=5)
m = 3
k = 5

error_rates = function(i){
  fold_ind = idx[[i]]
  pvals = rowttests(dat[,-fold_ind],factor(leukTF[-fold_ind]))$p.val
  gene_ind = order(pvals)[1:m]
  train_set = t(dat[gene_ind, -fold_ind])
  test_set = t(dat[gene_ind, fold_ind])
  train_classes = leukTF[-fold_ind]
  pred = knn(train_set, test_set, train_classes, k)
  sum(pred!=leukTF[fold_ind])
}

# Repeat the steps from questions 8 and 9 above for each of the 5 folds.

# A. What is the total number of errors across all 5 folds?

i = 1:5
num_errors = sum(sapply(i,error_rates))

# B. What proportion of the 60 samples are classified incorrectly by this model?

num_errors/length(leukTF)

# C. Accuracy is defined as 1 minus the error rate. What is the accuracy of this kNN model?

1 - num_errors/length(leukTF)
