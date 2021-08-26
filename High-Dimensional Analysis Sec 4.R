####################################################################################################################

# By: Santiago Taguado Menza
# High-Dimensional Data Analysis Section 3: Confounding, Singular Value Decomposition, & Factor Analysis
# April 14th, 2021
# HarvardX

####################################################################################################################

# On Confounding Errors: Simpsons Paradox

library(dagdata)
data(admissions)

# UC Berkeley Admissions Data
print(admissions)
head(admissions)

# Men Data
(index = which(admissions$Gender==1))
(accepted= sum(admissions$Number[index] * admissions$Percent[index]/100))
(applied = sum(admissions$Number[index]))
(not_accepted = applied - accepted)

(totals = c(accepted,not_accepted))

(acceptance_rate_men = accepted/applied)

# Women Data
index_1 = which(admissions$Gender==0)
(accepted_1= sum(admissions$Number[index_1] * admissions$Percent[index_1]/100))/
applied_1 = sum(admissions$Number[index_1])
(not_accepted_1 = applied_1 - accepted_1)

(acceptance_rate_women = accepted_1/applied_1)

# Question 2
(totals_1 = c(accepted_1,not_accepted_1))
# Table:  
#       accepted   not accepted     
# men   1198.02      1492.98
# women 556.62       1278.38 
z = rbind(totals,totals_1)
chisq.test(z)

# Admissions 
# We can quantify how "hard" a major is using the percent of students that were accepted. 
# Compute the percent that were accepted (regardless of gender) to each major and call this vector H.

(H_0 = admissions$Number[1:12])
(H_1 = admissions$Percent[1:12]/100)
(H_2 = H_0 * H_1)

ex= function(x,y){(H_2[x] + H_2[y])/(H_0[x]+H_0[y])}
x = seq(1,6,1)
y = seq(7,12,1)
mapply(ex,x,y)

admissions_rate_per_major = c((H_2[1] + H_2[7])/(H_0[1]+H_0[7]),(H_2[2] + H_2[8])/(H_0[2]+H_0[8]),(H_2[3] + H_2[9])/(H_0[3]+H_0[9]),
         (H_2[4] + H_2[10])/(H_0[4]+H_0[10]),(H_2[5] + H_2[11])/(H_0[5]+H_0[11]),(H_2[6] + H_2[12])/(H_0[6]+H_0[12]))
major = c("A","B","C","D","E","F")

(H = cbind(major,admissions_rate_per_major))

# Alternatively:
major = admissions[1:6,1]
men = admissions[1:6,]
women =admissions[7:12,]
H = (men$Number*men$Percent/100 + women$Number*women$Percent/100) / (men$Number+women$Number)
major[which.min(H)]

# Confounding Exercise 5
# For men, what is the correlation between the number of applications across majors and H?

(number_across = as.numeric(admissions$Number[1:6])) 
cor(number_across,H)

# Confounding Exercise 6
# For women, what is the correlation between the number of applications across majors and H?

(number_across = as.numeric(admissions$Number[7:12]) ) 
cor(number_across,H)

####################################################################################################################

# Confounding in Genomics 

# Note that this is the original dataset from which we selected the subset used in GSE5859Subset.  
# You can obtain it from the genomicsclass GitHub repository:

library(Biobase)
library(GSE5859)
data(GSE5859)

library(devtools)
install_github("genomicsclass/GSE5859")

# We can extract the gene expression data and sample information table using the Bioconductor functions exprs() 
# and pData() like this:

geneExpression = exprs(e)
sampleInfo = pData(e)

# Familiarize yourself with the sampleInfo table. Note that some samples were processed at different times. 
# This is an extraneous variable and should not affect the values in geneExpression. 
# However, as we have seen in previous analyses, it does appear to have an effect, so we will explore this here.

# You can extract the year from each date like this:

year = format(sampleInfo$date,"%y")
tab = table(year,sampleInfo$ethnicity)
print(tab)

(year_sampleInfo = sampleInfo %>% mutate(Year = format(sampleInfo$date,"%y")))

# Note there are 5 unique years and 3 unique ethnicities for which we have data.

length( unique(year) )

(ethnicity = format(sampleInfo$ethnicity))
length(unique(ethnicity))

(q = (cbind(year, ethnicity)))

# Do not form a matrix if you want to look at the data

tab = table(year,sampleInfo$ethnicity)
print(tab)
x = rowSums( tab != 0)
sum( x >= 2)

# For how many of these years do we have more than one ethnicity represented?
# 2

# Repeat the above exercise but now instead of year consider the month as well. 
# Specifically, instead of the year variable defined above, use:

month.year = format(sampleInfo$date,"%m%y")
tab = table (month.year, sampleInfo$ethnicity)
print(tab)
x = rowSums(tab != 0)
mean(x >= 2)

# Perform a t-test (use rowttests() from the genefilter package) comparing CEU samples processed in 2002 
# to those processed in 2003. Then use the qvalue package to obtain q-values for each gene.

library(genefilter)
library(qvalue)
year = factor( format(sampleInfo$date,"%y"))
idx = which(year%in% c("02","03") & sampleInfo$ethnicity=="CEU")
year = droplevels(year[idx])
pvals = rowttests(geneExpression[,idx], year)$p.value
q_vals = qvalue(pvals)
sum(q_vals$qvalues < 0.05)

# Now perform a t-test (use rowttests()) comparing CEU samples processed in 2003 to CEU samples processed in 2004. 
# Then use the qvalue package to obtain q-values for each gene.
# How many genes have q-values < 0.05?

year = factor( format(sampleInfo$date,"%y"))
index = which(year%in% c("03","04") & sampleInfo$ethnicity=="CEU")              
year = droplevels(year[index])
year
pval = rowttests(geneExpression[ ,index], year)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)
               
#Now we are going to compare ethnicities as was done in the original publication in which these data were first presented.
# Use the rowttests() function to compare the ASN population to the CEU population. Once again, 
# use the qvalue() function to obtain q-values.
               
# How many genes have q-values < 0.05?
               
(ethnicity = factor( sampleInfo$ethnicity))
(index = which(sampleInfo$ethnicity %in% c("CEU","ASN")))
(ethnicity = droplevels(ethnicity[index]))
pval = rowttests(geneExpression[ ,index],ethnicity)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

# Note that over 80% of genes are called differentially expressed between ethnic groups. 
# However, due to the confounding with processing date, we need to confirm these differences are actually due to 
# ethnicity. This will not be easy due to the almost perfect confounding. 
# However, above we noted that two groups were represented in 2005. 
# Just like we stratified by majors to remove the "major effect" in our admissions example, 
# here we can stratify by year and perform a t-test comparing ASN and CEU, but only for samples processed in 2005.

# How many genes have q-values < 0.05?

year = factor( format(sampleInfo$date,"%y"))
(ethnicity = factor( sampleInfo$ethnicity))
(index = which(sampleInfo$ethnicity %in% c("CEU","ASN") & year %in% c("05")))
(ethnicity = droplevels(ethnicity[index]))
pval = rowttests(geneExpression[ ,index],ethnicity)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)
    
# To provide a more balanced comparison, we repeat the analysis but now by taking 3 random CEU samples from 2002. 
# Repeat the analysis above but comparing the ASN from 2005 to three random CEU samples from 2002. 
# Set the seed at 3, set.seed(3), before random sampling.
               
RNGkind(sample.kind = "Rounding") 
sampleInfo = pData(e)
library(qvalue)
library(genefilter)

year = factor( format(sampleInfo$date,"%y") )
(index_1 = which(sampleInfo$ethnicity == "ASN" & year == "05"))
set.seed(3)
(index_2 = sample(which(sampleInfo$ethnicity == "CEU" & year == "02"),3))
index = c(index_1,index_2)          
g = droplevels(sampleInfo$ethnicity[index])
pvals_rand = rowttests(geneExpression[,index],g)$p.value
qval = qvalue(pvals_rand)
sum(qval$qvalue < 0.05)

####################################################################################################################

# Modeling Batch Effects

# For the dataset we have been working with, models do not help due to the almost perfect confounding. 
# This is one reason we created the subset dataset:

library(GSE5859Subset)
data(GSE5859Subset)

library(class)
library(genefilter)
library(qvalue)
library(qvalue)
library(genefilter)

# Here we purposely confounded month and group (sex) but not completely:

sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)  

# Using the functions rowttests() and qvalue() compare the two groups, in this case males and females coded in sex. 
# Because this is a smaller dataset, which decreases our power, we will use a more lenient FDR cut-off of 10%.

# How many gene have q-values less than 0.1?

sex = factor( sex)
pval = rowttests(geneExpression, sex)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.1)

# Note that sampleInfo$group here represents males and females. 
# Thus we expect differences to be on chrY and, for genes that escape inactivation, chrX. 
# Note that we do not expect many autosomal genes to be different between males and females. 
# This gives us an opportunity to evaluate false and true positives with experimental data. 
# For example, we evaluate results using the proportion genes of the list that are on chrX or chrY.

# For the list of genes with q<0.1 calculated in Modeling Batch Effects Exercises #1, 
# what proportion of genes are on chrX or chrY?

index = which(qval$qvalue < 0.1)
index_2 = which(geneAnnotation$CHR[index] =="chrY")
index_3 = which(geneAnnotation$CHR[index]=="chrX")
index_4 = c(index_2, index_3)

length(index_4)/length(index)

# Now for the autosomal genes (not on chrX and chrY) for which q-value < 0.1 perform a t-test comparing samples 
# processed in June to those processed in October.
# What proportion of these have p-values < 0.05?

# First Step remove sex chromosomes from the list of genes.
# Then perform a rowttest with using this index and the month.

library(qvalue)
library(genefilter)
sex = factor( sex)
pval = rowttests(geneExpression, sex)$p.value
qval = qvalue(pval)
qvals = qval$qvalues
index = which(qvals<0.1 & !geneAnnotation$CHR%in%c("chrX","chrY"))
month = factor( format(sampleInfo$date,"%m"))
pval = rowttests(geneExpression[index,], month)$p.value
mean(pval<0.05)

# The above result shows that the great majority of the autosomal genes show differences due to processing data. 
# This provides further evidence that confounding is resulting in false positives. 
# So we are going to try to model the month effect to better estimate the sex effect. 
# We are going to use a linear model.

# Which of the following creates the appropriate design matrix?

X = model.matrix(~sex+month)

# Now use the X defined above to fit a regression model using lm for each gene. 
# Note that you can obtain p-values for estimated parameters using summary(). Here is an example:

X = model.matrix(~sex+month)
i = 234
y = geneExpression[i,]
fit = lm(y~X-1)
summary(fit)$coef
summary(fit)$coef[2,c(1,4)]

# How many of the q-values for the group comparison are <0.1 now?

res <- t( sapply(1:nrow(geneExpression),function(i){
  y <- geneExpression[i,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
} ) )

##turn into data.frame so we can use the same code for plots as above
res <- data.frame(res)
names(res) <- c("dm","p.value") 
qvals = qvalue(res$p.value)

sum(qvals$qvalues < 0.1) 
index = which(qvals$qvalues < 0.1)

# With this new list, what proportion of these are chrX and chrY?

(sum(geneAnnotation$CHR[index]=="chrY",na.rm=TRUE) +
sum(geneAnnotation$CHR[index]=="chrX",na.rm=TRUE) )/length(index)

# Now, from the linear model in Modeling Batch Effects Exercises #6, 
# extract the p-values related to the coefficient representing the October versus June differences using the same 
# linear model.

# How many of the q-values for the month comparison are < 0.1 now?

res <- t( sapply(1:nrow(geneExpression),function(i){
  y <- geneExpression[i,]
  fit <- lm(y~X-1)
  summary(fit)$coef[3,c(1,4)]
} ) )
summary(fit)$coef[2,c(1,4)]
res <- data.frame(res)
names(res) <- c("dm","p.value") 
qvals = qvalue(res$p.value)
sum(qvals$qvalues < 0.1) 
test = which(qvals$qvalues < 0.1 & !geneAnnotation$CHR %in% c("chrX","chrY"))

####################################################################################################################

# Factor Analysis

library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

y = geneExpression - rowMeans(geneExpression)

# Compute and plot an image of the correlation for each sample. 
# Make two image plots of these correlations. 
# In the first one, plot the correlation as image. 
# In the second, order the samples by date and then plot the an image of the correlation. 
# The only difference in these plots is the order in which the samples are plotted.
# ased on these plots, which of the following you would say is true:

# Simple Version
myplibrary(rafalib)
mypar(1,2)
image(cor(y))
data = order(sampleInfo$date)
(x = y[,data])
image(cor(x))

# Advanced Version
library(rafalib)
sex = sampleInfo$group
mypar(1,2)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
cors = cor(y)
image(1:ncol(y),1:ncol(y),cors,col=cols,zlim=c(-1,1),
      xaxt="n",xlab="",yaxt="n",ylab="")
axis(2,1:ncol(y),sex,las=2)
axis(1,1:ncol(y),sex,las=2)
o = order(sampleInfo$date)
image(1:ncol(y),1:ncol(y),cors[o,o],col=cols,zlim=c(-1,1),
      xaxt="n",xlab="",yaxt="n",ylab="")
label = gsub("2005-","",sampleInfo$date[o])
axis(2,1:ncol(y),label,las=2)
axis(1,1:ncol(y),label,las=2)

# Based on the correlation plots above, we could argue that there are at least two hidden factors. 
# Using PCA estimate these two factors. Specifically, apply the svd() to y and use the first two PCs as estimates.
# Which command gives us these estimates?

(pcs = svd(y)$v[,1:2])

# Plot each of the estimated factor ordered by date. Use color to denote month. 
# The first factor is clearly related to date.
# Which of the following appear to be most different according to this factor?

u = pcs[,1]
max(abs(u))
times <-sampleInfo$date
(o=order(times))
plot(u,times[o],pch=21,ylab="date")

u1 = pcs[,2]
max(abs(u1))
times <-sampleInfo$date
(o=order(times))
plot(u1,times[o],pch=21,ylab="date")

pcs = svd(y)$v[,1:2]
o = order(sampleInfo$date)
month = factor( format(sampleInfo$date,"%m"))
cols = as.numeric(month)[o]
mypar(2,1)
for(i in 1:2){
  plot(pcs[o,i],col="blue",xaxt="n",xlab="")
  label = gsub("2005-","",sampleInfo$date[o])
  axis(1,1:ncol(y),label,las=2)
}

# Use the svd() function to obtain the principal components (PCs) for our detrended gene expression data y.
# How many principal components (PCs) explain more than 10% each of the variability?

s = svd(y)
(varexplained = s$d^2/ sum(s$d^2))
mypar(1,1)
plot(varexplained)
sum(varexplained>0.10)

# Which PC most correlates (negative or positive correlation) with month?
# What is this correlation (in absolute value)?

s = svd(y)
month = factor( format(sampleInfo$date,"%m"))
cors = cor( as.numeric(month),s$v)
mypar(1,1)
plot(t(cors))
max(abs(cors))

# Which PC most correlates (negative or positive correlation) with sex?
# What is this correlation (in absolute value)?

s = svd(y)
sex = factor( sex)
cors = cor( as.numeric(sex),s$v)
plot(t(cors))
which.max(abs(cors))
max(abs(cors))

# Now instead of using month, which we have shown does not quite describe the batch, add the two estimated factors 
# in Factor Analysis Exercises #6 to the linear model we used in previous exercises:
# Apply this model to each gene, and compute q-values for the sex difference.
# How many q-values are <0.1 for the sex comparison?

# Replicate this model for every row in the GeneExpression File
X <- model.matrix(~sex+s$v[,1:2])
fit <- lm(y[1,]~X)
summary(fit)
summary(fit)$coef[2,4]

# Provided Code
X= model.matrix(~sex+s$v[,1:2])
pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X)
  summary(fit)$coef[2,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)

# What proportion of the genes are on chrX and chrY?  
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

####################################################################################################################

# Surrogate Variable Analysis Function

# In this section we will use the sva() function in the sva package and apply it to the following data:

library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

# In the previous section we estimated factors using PCA. But we noted that the first factor was correlated 
# with our outcome of interest:

s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])

# As in the previous questions we are interested in finding genes that are differentially expressed between the 
# two groups (males and females in this case). Here we learn to use SVA to estimate these effects while using 
# a factor analysis approach to account for batch effects.

# The svafit() function estimates factors, but downweighting the genes that appear to correlate with the outcome 
# of interest. It also tries to estimate the number of factors and returns the estimated factors like this:

sex = sampleInfo$group
mod = model.matrix(~sex)
(svafit = sva(geneExpression,mod))
head(svafit$sv) 

#Note that the resulting estimated factors are not that different from the PCs:

for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}

# Now fit a linear model to estimate the difference between males and females for each gene but that instead of 
# accounting for batch effects using month it includes the factors estimated by sva in the model. 
# Use the qvalue() function to estimate q-values.

# How many genes have q-value < 0.1?

library(qvalue)
library(sva)
X= model.matrix(~sex+svafit$sv)
pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)

# What proportion of the genes from SVA Exercises #1 are from chrY or chrX?

index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

# Ploting log10 p-values and Differences

res = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
})

qvals = qvalue(res[2,])$qvalue
pcutoff = max( res[2,qvals < .1] )
library(rafalib)
mypar2(1,1)

plot(res[1,],-log10(res[2,]),xlab="M",ylab="log10 p-value")

ind = which(geneAnnotation$CHR=="chrY")
points(res[1,ind],-log10(res[2,ind]),col=1,pch=16)

ind = which(geneAnnotation$CHR=="chrX")
points(res[1,ind],-log10(res[2,ind]),col=2,pch=16)

abline(h=-log10(pcutoff))
legend("bottomleft",c("chrX","chrY"),col=c(2,1),pch=16)

####################################################################################################################

# The bladderbatch dataset from Bioconductor is a collection of gene expression data on bladder cancers 
# from 5 different batches.

# Quiz # 4

# Run Code 
BiocManager::install("bladderbatch")

library(bladderbatch)
data(bladderdata)

# Get the expression data
edata = exprs(bladderEset)
# Get the pheno data
pheno = pData(bladderEset)

# Create a reduced dataset containing only batches 1-3. Save the subsetted expression data as expr and 
# save the subsetted sample data as pdata:

ind = which(pheno$batch %in% 1:3)
expr = edata[,ind]
pdata = data.frame(batch = factor(pheno$batch[ind]),
                   cancer = factor(pheno$cancer[ind]))

# Make a table of cancer status by batch.
# Which of the following are true?

table(pdata$cancer,pdata$batch)

# Compare gene expression in the normal samples from batches 2 and 3. Use this code to extract the 
# relevant subset of the data:

# Use rowttests() from the genefilter package to compare expression across the two batches and 
# extract p-values. Then use the qvalue() function from the qvalue package to obtain q-values 
# for each gene.

# What proportion of genes have an FDR less than 0.05 when comparing normal samples across batches?

(index = which(pdata$cancer == "Normal"))
expr_norm = edata[ ,index]
(batch_norm = factor(pdata$batch[index]))

library(genefilter)
test_1 = rowttests(expr_norm,batch_norm)
qvals = qvalue(test_1$p.value)
mean(qvals$qvalues < 0.05)

# Use rowttests() from the genefilter library to find which genes in expr appear to be differentially 
# expressed between cancer and normal samples. Do not include batch effects. 
# Then use the qvalue() function from the qvalue package to obtain q-values for each gene.

# What proportion of genes appear differentially expressed between cancer and normal samples at an 
# q-value cutoff of 0.05?

cancer = factor(pdata$cancer)
pvals = rowttests(expr,cancer)
qvals_var_2 = qvalue(pvals$p.value)
mean(qvals_var_2$qvalues < 0.05)

# 0.6458735. Notice the number of significant genes went down. There are probably some batch effects present. 
# Well, of course, because within the pdata only batches 1 through 3 are being utilized.

# The pdata sample information associated with this experiment includes a variable batch. 
# It is not immediately clear what these batches represent, whether they include all the major sources 
# of experimental variability, and whether they will be useful for improving interpreation of the data.

# Define a model matrix X that includes both cancer status and batch as variables.

# Which of these commands correctly defines X?

X = model.matrix(~pdata$cancer + pdata$batch)

# Now use the model matrix X defined above to fit a regression model using lm() for each gene. 
# Note that you can obtain p-values for estimated parameters using summary(). 
# Here is an example for the first gene:
  
i = 1
y = expr[i,]
fit = lm(y~X-1)
summary(fit)$coef[2,c(1,4)]

# Find the p-value (Pr(>|t|)) for the expression difference between cancer and normal samples for 
# each gene. You can do this by modifying the example code above and using sapply(). 
# Then use the qvalue() function from the qvalue package to obtain q-values for each gene.

# A. What proportion of genes appear to be differentially expressed between cancer and normal samples 
# at a q-value cutoff of 0.05 when including batch in the model matrix?

res <- t( sapply(1:nrow(expr),function(i){
  y <- expr[i,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
} )) 

res <- data.frame(res)
names(res) <- c("dm","p.value") 
qvals = qvalue(res$p.value)

mean(qvals$qvalues < 0.05) 

# B. What proportion of genes appear to be differentially expressed between batch 1 and batch 2? 

res_1 <- t( sapply(1:nrow(expr),function(i){
  y <- expr[i,]
  fit <- lm(y~X-1)
  summary(fit)$coef[3,c(1,4)]
} )) 

res_1 <- data.frame(res_1)
names(res_1) <- c("dm","p.value") 
qvals_1 = qvalue(res_1$p.value)

mean(qvals_1$qvalues < 0.05) 

# C. What proportion of genes appear to be differentially expressed between batch 1 and batch 3?

res_2 <- t( sapply(1:nrow(expr),function(i){
  y <- expr[i,]
  fit <- lm(y~X-1)
  summary(fit)$coef[4,c(1,4)]
} )) 

res_2 <- data.frame(res_2)
names(res_2) <- c("dm","p.value") 
qvals_2 = qvalue(res_2$p.value)

mean(qvals_2$qvalues < 0.05) 

# Subtract the average expression of each gene from expr and save these results as y:

y = expr - rowMeans(expr)

# Use the svd() function to obtain the principal components (PCs) for our detrended gene expression data y.

s = svd(y)

# How many principal components (PCs) explain more than 5% each of the variability?

(varexplained = s$d^2/ sum(s$d^2))
plot(varexplained)
sum(varexplained>0.05)

# Plot the first 2 principal components on the x and y axis respectively. 
# Try coloring the points by either cancer status or batch number.

pcs = s$v[,1:2]

(cols = as.numeric(pdata$batch))
plot(pcs[,1],pcs[,2],col = cols,pch = 16)
legend("bottomleft",c("Batch 1","Batch 2","Batch 3"),col= 1:3,pch = 16,bty = "n")

(cols_1 = as.numeric(pdata$cancer))
plot(pcs[,1],pcs[,2],col = cols_1, pch = 16)
legend("bottomleft",c("Cancer","Normal"),col= 1:2,pch=16,bty = "n")

# What is the absolute value of the correlation coefficient between the first principal component and 
# cancer status?

(cor(pcs[,1],pdata$cancer == "Cancer"))

# Load the sva library and use it to infer the surrogate variables in expr other than cancer status.
# Define mod as a model matrix including cancer status as a variable. 
# Do not include batch as a variable - we will infer the batch effects with this approach. 
# Then, use sva() to estimate the surrogate variables and store the output as sv.

# How many significant surrogate variables affect the data?

mod = model.matrix(~pdata$cancer)

library(sva)
(sv = sva(expr,mod = mod))

# Question 10
# Define mod0 as a null model matrix:
  
mod0 = model.matrix(~1, data=pdata)

# The f.pvalue() function from sva quickly calculates p-values for each gene (row) given a design 
# matrix mod with the variable of interest and a null matrix mod0 that contains all variables except
# the variable of interest:

(fpvals = f.pvalue(expr, mod, mod0)  )

# Note that the q-values from this function are the same as the results from using rowttests() 
# in question 3:

fqvals = qvalue(fpvals)$qvalue
mean(fqvals < 0.05)

# Now, alter the alternative and null model matrices to adjust for the surrogate variables:

modSv = cbind(mod,sv$sv)
mod0Sv = cbind(mod0,sv$sv)

fvals_adj = f.pvalue(expr,modSv,mod0Sv)  
fqvals = qvalue(fvals_adj)$qvalue
mean(fqvals < 0.05)
index_finals = which(fqvals < 0.05)
mean(expr[index_finals,])
