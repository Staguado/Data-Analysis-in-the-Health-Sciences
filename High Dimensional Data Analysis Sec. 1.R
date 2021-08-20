####################################################################################################################

# By: Santiago Taguado Menza
# High-Dimensional Data Analysis Section 1: On Distance & Differences
# March 21, 2021

####################################################################################################################
####################################################################################################################

# On Distance & Differences in R
# If you have not done so already, install the data package tissuesGeneExpression.

library(devtools)
install_github("genomicsclass/tissuesGeneExpression")

# The data represents RNA expression levels for seven tissues, each with several biological replicates.
# We call samples that we consider to be from the same population, such as liver tissue from different 
# individuals, biological replicates:

library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)

# How many biological replicates are there for hippocampus?

table(tissue)

# What is the distance between samples 3 and 45?

d <- dist(t(e))
class(d)
e_mat <- as.matrix(d)
e_mat[3,45]

# What is the distance between gene 210486_at and 200805_at?

x <- e["210486_at",]
y <- e["200805_at",]
sqrt( crossprod(x-y) )

# If I run the command (don't run it!):

d = as.matrix(dist(e))

# How many cells (number of rows times number of columns) would this matrix have?

22215*(22215/1000000)

# Compute the distance between all pairs of samples:
# Read the help file for dist().
# How many distances are stored in d? (Hint: What is the length of d)?

length(d)

# Why is the answer above not ncol(e)^2?
# Because R takes advantage of symmetry: only the lower triangular matrix is stored, thus there are only 
# ncol(e)*(ncol(e)-1)/2 values.

####################################################################################################################
####################################################################################################################

# Quiz Week One

library(GSE5859Subset)
data(GSE5859Subset)

# Question 1 
# A: How many samples are in the dataset?

dim(geneExpression)
nrow(geneExpression)

# B: How many features are in the dataset?

# 24

# Question 2

# Inspect the sampleInfo data frame.
# A: How many samples are from the ethnicity "ASN"?

(index = which(sampleInfo$ethnicity %in% c("ASN")))
lenght(index)

# B: Which sample is from the ethnicity "CEU"?

(index_1 = which(sampleInfo$ethnicity %in% c("CEU")))

# Question 3

# Inspect the sampleInfo data frame.
# A: What is the distance between samples 3 and 7?

d <- dist(t(geneExpression,),diag = TRUE)
class(d)
(e_mat <- as.matrix(d))
e_mat[3,7]

# B What is the distance between samples 4 and 14?

e_mat[4,14]

# This code finds the mean distance between the first sample (column 1) and all other samples:
# Add an extra sapply() loop to this code to check the mean distance between each sample (column) and 
# all other samples.

row = 1:8793
nrow(geneExpression)
mean_dist = sapply(column, function(row){
  x = 1:nrow(geneExpression)
  dists = sapply(x, function(x){
    test = geneExpression[,x]
    target = geneExpression[,column]
    sqrt(crossprod(target-test))
  })
  mean(dists)
})

# Which sample (column) has the largest mean distance from other samples?

which.max(mean_dist)

# Use dist() to calculate the distance between all pairs of samples.
# What is the maximum distance between any two samples?

max(d)

# A: What is the distance between features "1007_s_at" and "201371_s_at"?

e_mat[,"1007_s_at"]
x <- geneExpression["1007_s_at",]
y <- geneExpression["201371_s_at",]
sqrt( crossprod(x-y) )

# B: What is the distance between features "202138_x_at" and "202152_x_at"?

x <- geneExpression[800,]
y <- geneExpression[123,]
sqrt( crossprod(x-y) )
x = geneExpression["1007_s_at",]
y = geneExpression["1053_at",]

# Use dist() to calculate the distance between all pairs of features.
# What is the maximum distance between any two features?

d2 = dist(geneExpression)
max(d2)