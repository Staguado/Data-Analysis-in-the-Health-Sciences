####################################################################################################################

# By: Santiago Taguado Menza
# High-Dimensional Data Analysis Section 2: On Projections & Dimension Reduction
# March 21, 2021
# HarvardX

####################################################################################################################
####################################################################################################################

# On Projections & Dimensions

library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

# Suppose you want to make an MA plot of the first two samples y = geneExpression[,1:2]. 
# Which of the following projections of  𝑦  gives us new coordinates such that column 2 versus column 
# 1 is an MA plot?

y = geneExpression[,1:2]
diffs <- y[,1] - y[,2]
average <- (y[,1] + y[,2])/2
plot(average,diffs)

c_1 <- c(1,1)
c_2 <- c(1,-1)
(conv <- cbind(c_1,c_2))

y_2 <- y %*% conv
diffs_1 <- y_2[,1] - y_2[,2]
average_1 <- (y_2[,1] + y_2[,2])/2
plot(average_1,diffs_1)

# Say  𝑌  is  𝑀×𝑁 , in the SVD  𝑌=𝑈𝐷𝑉⊤  which of the following is NOT correct?
# D are the coordinates of the projection U^t &*& Y

####################################################################################################################

# On Singular Value Decomposition

library(tissuesGeneExpression)
data(tissuesGeneExpression)

# Important note: When using the SVD in practice it is important to note that the solution to SVD is not 
# unique. This is because  𝐔𝐃𝐕⊤=(−𝐔)𝐃(−𝐕)⊤ . In fact we can flip the sign of each column of  𝐔
# and as long as we also flip the respective column in  𝐕  the decompostion works. Here is R code
# demonstrating this:

s = svd(e)
(signflips = sample(c(-1,1),ncol(e),replace=TRUE))

# Now we switch the sign of each column and check that we get the same answer. We do this using the 
# function sweep(). If x is a matrix and a is a vector, then sweep(x,1,y,FUN="*") applies the function 
# FUN to each row i FUN(x[i,],a[i]), in this case x[i,]*a[i]. If instead of 1 we use 2, sweep() applies 
# this to columns. To learn about sweep(), read ?sweep.

newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

# Question 1
# Compute the SVD of e:
# Now compute the mean of each row:

m = rowMeans(e)

# What is the correlation between the first column of  𝐔  and m?

cor(s$u[,1],m)

# Question 2

# In the above question, we saw how the first column relates to the mean of the rows of e. 
# Note that if we change these means, the distances between columns do not change. Here is some R 
# code showing how changing the means does not change the distances:

(newmeans = rnorm(nrow(e))) ##random values we will add to create new means
newe = e+newmeans ##we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))

# So we might as well make the mean of each row 0 since it does not help us approximate the column distances.
# We will define y as the detrended e and recompute the SVD:

y = e - rowMeans(e)
s = svd(y) 

# We showed that  𝐔𝐃𝐕⊤  is equal to y up to numerical error:

resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

# The above can be made more efficient in two ways. First, using the crossprod() and second not creating 
# a diagonal matrix. Note that in R we can multiply a matrix x by vector a. The result is a matrix with row 
# i equal to x[i,]*a[i]. Here is an example to illustrate this.

x=matrix(rep(c(1,2),each=5),5,2)
x*c(1:5)

# Note that the above code is actually equivalent to:

sweep(x,1,1:5,"*")

# This means that we don't have to convert s$d into a matrix to obtain  𝐃𝐕⊤ .

# Question 2 
# Which of the following gives us the same as diag(s$d)%*%t(s$v)?

diag(s$d)  %*% t(s$v)
s$d * t(s$v)

# Question 3

# If we define vd = t(s$d * t(s$v)), then which of the following is not the same as  𝐔𝐃𝐕⊤  :
# Define VD

vd = t(s$d * t(s$v))
tcrossprod(s$u,vd)

# Correct Answer
s$u %*% s$d * t(s$v)
s$u %*% (s$d * t(s$v) )
s$u %*% t(vd)

# Order of operations in R go left to right and s$u %*% s$d is multiplying non-conformable matrices

# Question 4
# Let z = s$d * t(s$v). We showed a derivation demonstrating that because  𝐔  is orthogonal,
# the distance between e[,3] and e[,45] is the same as the distance between y[,3] and y[,45], 
# which is the same as z[,3] and z[,45]:

z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

# Note that the columns z have 189 entries, compared to 22,215 for e.
# What is the difference (in absolute value) between the actual distance "sqrt(crossprod(e[,3]-e[,45]))" 
# and the approximation using only two dimensions of z?

realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistance = sqrt(crossprod(z[1:2,3]-z[1:2,45]))
abs(realdistance - approxdistance)

# Question 5
# What is the minimum number of dimensions we need to use for the approximation in SVD Exercises #4 to 
# be within 10% or less?

# Crude Method: Trial and error until achieving less than 10%
realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistance = sqrt(crossprod(z[1:7,3]-z[1:7,45]))
(abs(realdistance - approxdistance))/realdistance

# Elegant Method: Creating a function 
ks = 1:189
realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistances = sapply(ks,function(k){
  sqrt(crossprod(z[1:k,3,drop = FALSE]-z[1:k,45,drop = FALSE] )) 
})
percentdiff = 100*abs(approxdistances - realdistance)/realdistance
plot(ks,percentdiff) ##take a look
min(ks[which(percentdiff < 10)])

# Question 6
# Compute distances between sample 3 and all other samples:

(distances = sqrt(apply(e[,-3]-e[,3],2,crossprod)))

# Recompute this distance using the 2 dimensional approximation.

approx = sqrt(apply(z[1:2,-3] - z[1:2,3],2,crossprod))

# What is the Spearman correlation between this approximate distance and the actual distance?

cor(distances, approx, method = "spearman")

####################################################################################################################

# On Multi Dimensional Scaling Plots
# For the following questions, use the data loaded with:

library(tissuesGeneExpression)
data(tissuesGeneExpression)  

# In these exercise we will demonstrate the relationship between the SVD and the output of cmdscale(), 
# the function in R that performs MDS.

# Using the z we computed in SVD Exercises #4:

y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)

# we can make an MDS plot:

library(rafalib)
ftissue = factor(tissue)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1,bty='n')

# Now run the function cmdscale() on the original data:

d = dist(t(e))
mds = cmdscale(d)

# What is the correlation between the first row of z and the first column in mds?

cor(z[1,],mds[,1])

# What is the correlation between the second row of z and the second column of mds?

cor(z[2,],mds[,2])  

# Note that the MDS plot is not the same:

library(rafalib)
ftissue = factor(tissue)
mypar(1,2)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1,bty='n')
plot(mds[,1],mds[,2],col=as.numeric(ftissue))

# Given the answer to MDS Exercises #1 and #2, what do we have to do to z[1,] and z[2,] 
# to get a practically identical plot?

library(rafalib)
ftissue = factor(tissue)
mypar(1,1)
new_z <- z[1,] * -1
new_z_1 <- z[2,] * -1
plot(new_z, new_z_1,col=as.numeric(ftissue))
legend("topright",inset = c(-0.2,0),levels(ftissue),col=seq_along(ftissue),pch = 1,xpd = TRUE)
plot(mds[,1],mds[,2],col=as.numeric(ftissue))

# Question 4

library(GSE5859Subset)
data(GSE5859Subset)
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)


month = as.numeric(factor(format( sampleInfo$date, "%m")))

cor_func <- function(x){
  cor(z[x,],month)
}

k <- c(1:24)
sapply(k,cor_func)

# Continue working with the z calculated from the GSE5859Subset data.
# What is this max correlation?

# 0.62368575

# Which dimension of z has the second highest correlation with the outcome sampleInfo$group?

which.max(cor(sampleInfo$g,t(z))[-1]) + 1

# Note these measurements were made during two months:

sampleInfo$date

# We can extract the month this way:

(month = format( sampleInfo$date, "%m"))
(month = (factor( month)))

# Which dimension of z has the highest correlation with the outcome month?

which.max(cor( as.numeric(month), t(z)))

# What is this correlation?

max(cor( as.numeric(month), t(z)))

# Note: this is an advanced question. Please feel free to research this question online.
# In MDS Exercises #7 we saw that that one of the dimensions was highly correlated to the sampleInfo$group. 
# Now take the 5th column of  𝐔  and stratify by the gene chromosome.
# Remove chrUn and make a boxplot of the values of  𝐔6  stratified by chromosome.

# Which chromosome looks different from the rest?

chromosomes <- factor(geneAnnotation$CHR)
u_6 <- s$u[,6]
one <- cbind(chromosomes,u_6)
results_2 <- split(u_6,geneAnnotation$CHR)
boxplot(results_2)


result = split(s$u[,6],geneAnnotation$CHR)
result = result[ which(names(result)!="chrUn") ]
boxplot(result,range=0)
boxplot(result,range=0,ylim=c(-0.025,0.025))
medians = sapply(result,median)
names(result)[ which.max(abs(medians)) ]

# chrY

####################################################################################################################

# Singular Value Decomposition & Distance Concept Review

library(GSE5859Subset)
data(GSE5859Subset)

# Question 1
# Which of the following are true about the singular value decomposition matrices  𝐘𝐕=𝐔𝐃 ?

# Y is the original matrix
# U and V are orthogonal matrices meaning that the inverse and the original are equal
# D is a diagonal with decreasing value
# UD are the new coordinates for the projection YV
# This equation can be written as Y = U^t*D*V

# Compute the SVD of geneExpression. Save it as the variable s.
# A: What is the first entry of s$d?
s = svd(geneExpression)

# The proportion of variability in the data explained by the xth column of  𝐔  is equal to s$d[x]^2 divided
# by the sum of all s$d values. What proportion of variability is explained by the first column of  𝐔 ?

s$d[1]^2/sum(s$d^2)

# C: Compute the mean of each row of geneExpression as a vector m. 
# What is the correlation between m and the first column of s$u?

nrow(geneExpression)
m  = rowMeans(geneExpression)
u_1 = s$u[,1]
cor(m,u_1)

# Which of the following are true about the singular value decomposition matrices  𝐘𝐕=𝐔𝐃 ?

# The row means are almost perfectly correlated with the first column of  𝐔  aside from a sign change.
# Most of the variability in the gene expression matrix is driven by average expression levels of each 
# gene rather than biological differences between samples.
# Removing the row means before computing the SVD would help reveal the underlying biological signal.

# Define y as geneExpression - rowMeans(geneExpression), then compute the SVD of y and save the result as s.

y = geneExpression - rowMeans(geneExpression)
s = svd(y)

# A: What is the first entry of s$d?
s$d

# What proportion of the variability is explained by the first column of U?

s$d[1]^2/sum(s$d^2)

# C: Calculate the proportion of variability explained by each column of  𝐔 .
# How many individual columns explain more than 5% of the variability?

l = 1:24
variability = sapply(l, function(l){
  s$d[l]^2/sum(s$d^2)})

variability
sum(variability > 0.05)

# D: What percent of variability is explained by the first 10 columns of

sum(variability[1:10])

# Confirm that  𝐘=𝐔𝐃𝐕⊤
# . First, multiply the matrices within s to regenerate the data matrix. 
# Then, subtract the regenerated matrix from the original matrix to find residuals. 
# Find and report the residual with the maximum absolute value to show the matrices 
# are identical to within numerical error.

# What is the value of the residual with the maximum absolute value?

actual =  s$u %*% (s$d * t(s$v))

residuals = abs(actual - y)
max(residuals)

# Let z = s$d * t(s$v). 
# Compare the distance between columns 1 and 2 in geneExpression (the original matrix), 
# y (the de-trended matrix with row means subtracted), and z.

# A: What is the distance between columns 1 and 2 in geneExpression?

d1 = dist(t(geneExpression))
d1 = as.matrix(d1)

sqrt(crossprod(geneExpression[,1] - geneExpression[,2]))

# B: What is the distance between columns 1 and 2 in y?

sqrt(crossprod(y[,1] - y[,2]))

# Now using dist
# C: What is the distance between columns 1 and 2 in z?

z = s$d * t(s$v)
sqrt(crossprod(z[,1] - z[,2]))

# Now using dist
# D: What is the distance between 1 and 2 in z using only the first 10 columns as an approximation?

sqrt(crossprod(z[1,1:10] - z[2,1:10]))

# Perform MDS on the original geneExpression data:

d = dist(t(geneExpression))
mds = cmdscale(d)

# Which of the following is true about the relationship of date to the MDS plot?
  
library(rafalib)
mypar(1,2)
fdate = factor(sampleInfo$date)
plot(mds[,1],fdate,,col=as.numeric(fdate))
legend("topleft",levels(fdate),col=seq_along(fdate),pch = 1,bty='n')
plot(mds[,2],fdate,,col=as.numeric(fdate))

# The first dimension of the MDS plot appears to correlate with date, with earlier dates to the left and later dates 
# to the right