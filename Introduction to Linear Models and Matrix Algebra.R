# By: Santiago Taguado Menza
# Introduction to Linear Models and Matrix Algebra
# Harvard Course

#-------------------------------------------------------------------------------------------------

# Packages Required
install.packages("UsingR")
install.packages("dplyr")
install.packages("contrast")
install.packages("downloader")

# Once packages are loaded, use the following to code to access Galton's father and son heights:
library(UsingR)
data("father.son",package="UsingR")

#-------------------------------------------------------------------------------------------------

# Data Exploration:

# What is the average height of the sons (don't round off)?

mean(father.son$fheight)
mean(father.son$sheight)

# One of the defining features of regression is that we stratify one variable based on others. 
# In Statistics, we use the verb "condition". 
# For example, the linear model for son and father heights answers the question how tall do I expect 
# a son to be if I condition on his father being x inches. 
# The regression line answers this question for any x.

# Using the father.son dataset described above, we want to know the expected height of sons if we condition 
# on the father being 71 inches. 

# Create a list of son heights for sons that have fathers with heights of 71 inches (round to the nearest inch).

# What is the mean of the son heights for fathers that have a height of 71 inches (don't round off your answer)?
# Hint: use the function round() on the fathers' heights)

father_son_refined <-  father.son %>% mutate((round(father.son$fheight))) %>% filter(round(father.son$fheight) == 71)
mean(father_son_refined$sheight)

# We say a statistical model is a linear model when we can write it as a linear combination of parameters 
# and known covariates plus random error terms. In the choices below,  𝑌  represents our observations,
# time  𝑡  is our only covariate, unknown parameters are represented with letters
# 𝑎,𝑏,𝑐,𝑑  and measurement error is represented by the letter e
# Note that if  𝑡  is known, then any transformation of𝑡is also known.

# So, for example, both  𝑌=𝑎+𝑏𝑡+𝑒  and  𝑌=𝑎+𝑏𝑓(𝑡)+𝑒  are linear models.
# Which of the following can't be written as a linear model?

# Answer => Y = a + b^t + e. This cannot be transformed to a linear model.

# Suppose you model the relationship between weight and height across individuals with a linear model. 
# You assume that the height of individuals for a fixed weight  𝑥  follows a linear model  𝑌=𝑎+𝑏𝑥+𝑒 .

# Which of the following best describes what  𝑒  represents?
# Between-individual variability: people of the same weight vary in their height.

#------------------------------------------------------------------------------------------------------------

# Matrix Notation & Extraction:

# Create the matrix from the vector 1:1000 like this:
  
X = matrix(1:1000,100,10)

# What is the entry in row 25, column 3 ?

X[25,3]

# Matrix Notation Exercises

# Using the function cbind(), create a 10 x 5 matrix with first column x=1:10. 
# Then use 2*x, 3*x, 4*x and 5*x respectively in columns 2 through 5.

# What is the sum of the elements of the 7th row?

x = 1:10
Z <- cbind(x,2*x,3*x,4*x,5*x)
sum(Z[7,])

# Which of the following creates a matrix with multiples of 3 in the third column?

matrix(1:60,20,3,byrow=TRUE)

# Whats the last element of the vector returned by seq(10,1,-2)?
  
seq(10,1,-2)
# 2

# Let Y be the following:

q = c(4,5)
w = c(9,0)
Y = cbind(q,w)
  
# Which of the following is NOT equivalent to X?

Y %*% matrix(1,ncol(Y) )
Y%*%diag(ncol(Y))

# Solve the following system of equations using R:
  
#3𝑎+4𝑏−5𝑐+𝑑=10 
#2𝑎+2𝑏+2𝑐−𝑑=5 
#𝑎−𝑏+5𝑐−5𝑑=7 
#5𝑎+𝑑=4

e = c(3,4,-5,1) 
r = c(2,2,2,-1)
t = c(1,-1,5,-5) 
y = c(5,0,0,1)
u = cbind(c(10,5,7,4))

U = rbind(e,r,t,y)
solve(U) %*% u

# Load the following two matrices into R:
  
a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)

# Note the dimensions of a and the dimensions of b.

# In the question below, we will use the matrix multiplication operator in R, %*%, 
# to multiply these two matrices.

# What is the value in the 3rd row and the 2nd column of the matrix product of a and b?

c = a%*%b
c[3,2]

# Multiply the 3rd row of a with the 2nd column of b, using the element-wise vector multiplication with *.
# What is the sum of the elements in the resulting vector?

a_1 <- a[3,]
b_1 <- b[,2]

sum(a_1*b_1)

#------------------------------------------------------------------------------------------------------------

# Linear Models: Galton's Father/Son Height & Gravity on The Tower of Pisa

# Fitting Galton's Data
data(father.son,package="UsingR")
x=father.son$fheight
y=father.son$sheight
X <- cbind(1,x)
betahat <- solve( t(X) %*% X ) %*% t(X) %*% y
# Alternatively,
betahat <- solve( crossprod(X) ) %*% crossprod( X, y )
newx <- seq(min(x),max(x),len=100)
X <- cbind(1,newx)
fitted <- X%*%betahat
plot(x,y,xlab="Father's height",ylab="Son's height")
lines(newx,fitted,col=2)

# The Tower of Pisa
set.seed(1)
g <- 9.8 #meters per second
n <- 25
tt <- seq(0,3.4,len=n) #time in secs, t is a base function
d <- 56.67  - 0.5*g*tt^2 + rnorm(n,sd=1)
X <- cbind(1,tt,tt^2)
y <- d
betahat <- solve(crossprod(X))%*%crossprod(X,y)
newtt <- seq(min(tt),max(tt),len=100)
X <- cbind(1,newtt,newtt^2)
fitted <- X%*%betahat
plot(tt,y,xlab="Time",ylab="Height")
lines(newtt,fitted,col=2)

# Data Anlysis of Treatment Groups
# Suppose we are analyzing a set of 4 samples. The first two samples are from a treatment group A 
# and the second two samples are from a treatment group B. This design can be represented with a model 
# matrix like so:
  
X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
rownames(X) <- c("a","a","b","b")

# Suppose that the fitted parameters for a linear model give us:
  
beta <- c(5, 2)

# Use the matrix multiplication operator, %*%, in R to answer the following questions:

samples_a <- X[1:2,]
samples_b <- X[3:4,]

samples_a %*% beta
samples_b %*% beta

#------------------------------------------------------------------------------------------------------------

# Monte Carlo on Gravity Model to Identify Estimates:

# We have shown how to find the least squares estimates with matrix algebra. 
# These estimates are random variables as they are linear combinations of the data. 
# For these estimates to be useful we also need to compute the standard errors.

# Here we review standard errors in the context of linear models.

# It is useful to think about where randomness comes from. 
# In our falling object example, randomness was introduced through measurement errors. 
# Every time we rerun the experiment a new set of measurement errors will be made which implies our 
# data will be random. This implies that our estimate of, for example, the gravitational constant will change. 
# The constant is fixed, but our estimates are not. 
# To see this we can run a Monte Carlo simulation. 
# Specifically we will generate the data repeatedly and compute the estimate for the quadratic term each time.

g = 9.8 ## meters per second
h0 = 56.67
v0 = 0
n = 25
tt = seq(0,3.4,len=n) ##time in secs, t is a base function
y = h0 + v0 *tt  - 0.5* g*tt^2 + rnorm(n,sd=1)

# Now we act as if we didn't know h0, v0 and -0.5*g and use regression to estimate these. 
# We can rewrite the model as y = b0 + b1 t + b2 t^2 + e and obtain the LSE we have used in this class. 
# Note that g = -2*b2.

# To obtain the LSE in R we could write:

X = cbind(1,tt,tt^2)
A = solve(crossprod(X))%*%t(X)

# Given how we have defined A, which of the following is the LSE of g, the acceleration due to gravity?
# Suggestion: try the code in R.

set.seed(1)
gs <- replicate(100000,
                     {y = h0 + v0 *tt  - 0.5* g*tt^2 + rnorm(n,sd=1) 
                       X = cbind(1,tt,tt^2)
                       A = solve(crossprod(X))%*%t(X)
                       2 * (A %*% y) [3]})
mean(gs)
sd(gs)

#---------------------------------------------------------------------------------------------------------

# Using the lm function in R, Effects of Random Sampling and Standard Errors of Beta:

library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)
N = 50
set.seed(1)
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef

# Note that the fitted values𝑌̂from a linear model can be obtained with:

fit = lm(y ~ x)
summary(fit)
fit$fitted.values

r <-  (sampledat$sheight - fit$fitted.values)^2 
sum(r)
ssr <- sum(r)/48
X = cbind(rep(1,N), x)
solve(t(X) %*% X)

# Now we are one step away from the standard error of beta-hat. 
# Take the diagonals from the  (𝑋𝑇𝑋)−1  matrix above, using the diag() funct
# Now multiply our estimate of  𝜎2  and the diagonals of this matrix.
# This is the estimated variance of beta-hat, so take the square root of this. 
# You should end up with two numbers, the standard error for the intercept and the standard error 
# for the slope.

# What is the standard error for the slope?

sqrt(diag(solve(t(X) %*% X))*ssr)

# Suppose now we are comparing two treatments B and C to a control group A, each with two samples. 
# This design can be represented with a model matrix like so:
  
X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
rownames(X) <- c("a","a","b","b","c","c")
# which results in a matrix that looks like

# a 1 0 0
# a 1 0 0
# b 1 1 0
# b 1 1 0
# c 1 0 1
# c 1 0 1

# Suppose that the fitted values for the linear model are given by:
beta <- c(10,3,-3)
fitted_values <- X  %*% beta 

# Question # 1 - What is the fitted value for the B samples?
# 13
X  %*% beta  

# Question # 2 - What is the fitted value for the C samples?
# 7

# Question # 3 
# In the father and son height examples we have randomness because we have a random sample of father and 
# son pairs. For the sake of illustration let's assume that this is the entire population:

library(UsingR)
X = father.son$fheight
Y = father.son$sheight
n = length(y)

# Now let's run a Monte Carlo simulation in which we take a sample of size 50 over and over again. 
# Here is how we obtain one sample:

N =  50
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat =  lm(y~x)$coef

# Use the function replicate() to take 10,000 samples.
# What is the standard error of the slope estimate? That is, calculate the standard deviation of the estimate
# from many random samples. Again, set the seed to 1 before using replicate().

set.seed(1)
galtons_replications <- replicate(10000, {
  index = sample(n,N)
  sampledat = father.son[index,]
  x = sampledat$fheight
  y = sampledat$sheight
  betahat =  lm(y~x)$coef
})

# Out of curiosity, here is the standard error for the intercept:
sd(galtons_replications[1,])
sd(galtons_replications[2,])

# Question # 4
# We are defining a new concept: covariance. The covariance of two lists of numbers  𝑋=𝑋1,...,𝑋𝑛  and  𝑌=𝑌1,..
# .,𝑌𝑛  is mean((Y - mean(Y))*(X-mean(X))).

# Which of the following is closest to the covariance between father heights and son heights?

(covariance <-((Y - mean(Y)) * (X - mean(X))))
mean(covariance)
cov(X,Y)
cov(Y,X)
# But what is covariance? Covaiance is measures the joint variability of two variables

# Question # 5 
# Look back at the Standard Errors exercises that used father and son heights. 
# In that problem set we quite sensibly created a model with the son's height as a function of the 
# father's height: fit = lm(y ~ x). 
# What if we had accidentally defined the reverse function: fit = lm(x ~ y), with the father's 
# height as a function of the son's?
# How would our model's slope change if we reversed the variables?

fit_1 = lm(y ~ x)
summary(fit_1)

fit_2 = lm(x ~ y)
summary(fit_2)

# Nothing would change. The estimate would lower because the mean average height would reduce.
# Or as Prof. Irizarry puts it, the only values that will be the same for both models will be items related 
# to the correlation coefficient: R2, t value, etc. 
# The slope, the intercept, and their standard errors will all be different, and they will not have 
# such a simple relationship to the old values.

#---------------------------------------------------------------------------------------------------------

# Experimental Design: Creating Model Matrices

# Suppose we have an experiment with the following design: 
# on three different days, we perform an experiment with two treated and two control samples. 
# We then measure some outcome  𝑦𝑖 , and we want to test the effect of condition
# n, while controlling for whatever differences might have occurred due to the the different day 
# (maybe the temperature in the lab affects the measuring device). 

# Assume that the true condition effect is the same for each day (no interaction between condition and day). 
# We then define factors in R for day and for condition.
# day: A B C
# condition: --------------
# treated | 2 2 2
# control | 2 2 2

# Given the factors we have defined above, and not defining any new ones, 
# which of the following R formula will produce a design matrix (model matrix) 
# that let's us analyze the effect of condition, controlling for the different days:

# Factors
# day: A B C
# condition: --------------
# treated | 2 2 2
# control | 2 2 2

# ~ day + condition
# Using the ~ and the names for the two variables we want in the model will produce a 
# design matrix controlling for all levels of day and all levels of condition, so ~ day + condition. 
# We do not use the levels A,B,C etc in the design formula.

# SE = sqrt(var(diff))
# var(diff) = (1/nx + 1/ny) (sum {(x_i - mu_x)^2} + sum{(y_i - mu_y)^2}) / (nx + ny - 2)

# We can make a design matrix X for a two group comparison either using model.matrix() or simply with:

ny = 7
nx = 5
X = cbind(rep(1,nx + ny),
          rep(c(0,1),
              c(nx, ny)))
solve(t(X)%*%X) 

# Alternatively, we can use model.matrix:
species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
X <- model.matrix(~ species + condition)

# Suppose we want to build a contrast of coefficients for the above experimental design.
# You can either figure this question out through logic, by looking at the design matrix, 
# or using the contrast() function from the contrast library. If you have not done so already, 
# you should download the contrast library. The contrast vector is returned as contrast()$X.

# What should the contrast vector be, for the contrast of (species=B and condition=control) vs (species=A and condition=treatment)? 
# Assume that the beta vector from the model fit by R is: Intercept, speciesB, condition treated.

y = rnorm(4)
fit = lm(y ~ species + condition)
contrast(fit, list(species="B",condition="control"), list(species="A",condition="treated"))$X

# notes: Create a random distribution and combine when binary numbers are present

#---------------------------------------------------------------------------------------------------------

# Spider Legs & Friction Study Analysis

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)
head(spider)

# Run the code from the Rmd script of the spider dataset. 
# Suppose we build a model using two variables: ~ type + leg.

# What is the t-value for the contrast of leg pair L4 vs leg pair L2?

U <- model.matrix(~ type + leg, data = spider)
fitTL <- lm(friction ~ type + leg, data=spider)
summary(fitTL)
(coefs <- coef(fitTL))

L4vsL2 <- contrast(fitTL,list(leg="L4",type="push"),list(leg="L2",type="push"))
L4vsL2

spider$log2friction <- log2(spider$friction)
boxplot(log2friction ~ type*leg, data=spider)

# What is the t-value for the interaction of type push and leg L4? 
# If this t-value is sufficiently large, we would reject the null hypothesis that the push vs pull effect 
# on log2(friction) is the same in L4 as in L1.

print(spider$log2friction)
U <- model.matrix(~ type + leg + type:leg, data = spider)
fitTL2 <- lm(spider$log2friction ~ type + leg + type:leg, data=spider)
summary(fitTL2)

# The t-value is sufficiently large and we would reject the null hypothesis.

# Anove Test
anova(fitTL2)

L2vsL1 <- contrast(fitTL2, list(leg = "L2", type = "push"), list(leg = "L1", type = "push"))
L2vsL1

#---------------------------------------------------------------------------------------------------------

# Collinearity Problems

Sex <- c(0,0,0,0,1,1,1,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")

# Which of the above design matrices does NOT have the problem of collinearity?

# You can check in R, the rank of the E matrix is equal to the number of columns, 
# so all of the columns are independent.

m = matrix(c(1,1,1,1,0,0,1,1,0,1,0,1,0,0,0,1),4,4)
qr(m)$rank
  
sex <- factor(rep(c("female","male"),each=4))
trt <- factor(c("A","A","B","B","C","C","D","D"))
X <- model.matrix( ~ sex + trt)
qr(X)$rank
Y <- 1:8
makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b
fitTheRest <- function(a,b) {
  Ystar <- makeYstar(a,b)
  Xrest <- X[,-c(2,5)]
  betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
  residuals <- Ystar - Xrest %*% betarest
  sum(residuals^2)
}
Xrest
fitTheRest(1,2)

expand.grid(1:3,1:3)
betas = expand.grid(-2:8,-2:8)
rss = apply(betas,1,function(x) fitTheRest(x[1],x[2]))
print(rss)

fitTheRest(8,-2)
fitTheRest(6,0)
fitTheRest(1,5)

# There is no minimum

library(rafalib)

## plot the pairs what are minimum

themin=min(rss)

rss
plot(betas[which(rss==2),]) 

#---------------------------------------------------------------------------------------------------------

# QR Decomposition Computation Exercises on Spider Study

fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)
Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)

QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)

solve(R) %*% (t(Q) %*% Y)
betahat

# Quiz 3

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))
C <- matrix(c(0,0,-1,0,1),1,5)
# Note the diagonal elements of the sigma function give the betahats. The covariances are the intersection
# between variables

# With the function below you can create the C function
L4vsL2 <- contrast(fitTL,list(leg = "L4", type = "push"), list(leg = "L2", type = "push"))
L4vsL2$X

# sqrt(Var(beta-hat_L4 - beta-hat_L2)) = sqrt(Var(beta-hat_L4) + Var(beta-hat_L2) 
# - 2 Cov(beta-hat_L4, beta-hat_L2))

L4vsL2$SE
# By rearranging the formula above we can actually find the covariance
-((C %*% Sigma %*% t(C)) - 0.0011819981 - 0.0020871318)/2 

#---------------------------------------------------------------------------------------------------------

# On the F-Test

N <- 40
p <- 4
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)
set.seed(1)

Fvalues = replicate(1000, {
Y <- rnorm(N,mean=42,7)
mu0 <- mean(Y)
initial.ss <- sum((Y - mu0)^2)
s <- split(Y, group)
after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
group.ss <- initial.ss - after.group.ss
group.ms <- group.ss / (p - 1)
after.group.ms <- after.group.ss / (N - p)
group.ms / after.group.ms
})
mean(Fvalues)
hist(Fvalues)

# On your own, you may wish to plot the distribution of the 1000 F-values:
  
# hist(Fs, col="grey", border="white", breaks=50, freq=FALSE)
# Overlay the theoretical F-distribution, with parameters df1=p - 1, df2=N - p.

# xs <- seq(from=0,to=6,length=100)
# lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")
# This is the distribution which is used to calculate the p-values for the ANOVA table produced by anova().

# Quiz 1 

# Question Number 3

8848 - 0.5*9.8*(100)

# Question Number 4

# On a matrix [row, columns]

# Question 5

2.5 : 6.5
seq(2.5,6.5)

# Question 6

c(seq(1,2),seq(3,4))
# This is evidently a vector

# Question 7

q = c(2,1)
w = c(6,3)
Y = cbind(q,w)
solve(Y)
# An error is displayed meaning it is not invertible

x + y + z = 6
2y + 5z = −4
2x + 5y − z = 27

q = c(1,0,2) 
w = c(1,2,5)
e = c(1,5,-1)
r = cbind(c(6,-4,27))
f = cbind(q,w,e)

solve(f)%*%r

# On your own, you may wish to plot the distribution of the 1000 F-values:

hist(Fvalues, col="grey", border="white", breaks=50, freq=FALSE)

# Overlay the theoretical F-distribution, with parameters df1=p - 1, df2=N - p.

xs <- seq(from=0,to=6,length=100)
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")

# This is the distribution which is used to calculate the p-values for the ANOVA table produced by anova(). 

#---------------------------------------------------------------------------------------------------------
