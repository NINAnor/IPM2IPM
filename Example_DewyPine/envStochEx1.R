
## IPM course Barcelona November 2017
### R code by Maria Paniw
### Exercise on environmental stochasticity

## This exercise is intended to show you how to construct IPMs
## using mixed effect models for vital rates 
## and how to incorporate random variation (random effect in vital rate models)
## into population porjections

### This exercise consists of 3 parts

### Part A: Vital rate functions 

### Part B: Construct IPM kernels and simulate log lambda.s using kernel selection

### Part C: Construct IPM kernels and simulate log lambda.s using parameter selection

###########################################################################
### PART A - VITAL RATE FUNCTIONS
##########################################################################

# set working directory (change to appropriate one)
setwd("/Users/maria/Dropbox/teaching/IPM2017")

data=read.csv("dataDroso.csv")

head(data)

data$year=as.factor(data$year)

# Now, the study species is adapted to fire regimes, but let's say we don't know this,
# and go to the field to collect data from habitats with TSF >3

data=droplevels(data[data$TSF==">three",])

# load lme4 package for mixed effect models
library(lme4)

#### Vital rate models

# Empty list to save model coefficients 
paramCont=list(NULL)

# Run and compare mixed model for survival and save parameters of best model (determined by LRT)
surv0=glmer(surv~1+(1|year),family=binomial,data=data)
surv1=glmer(surv~size+(1|year),family=binomial,data=data)

AIC(surv0,surv1)

paramCont[[1]]=as.matrix(coef(surv1)$year) # save coefficients 

# Run and compare mixed model for growth and save parameters of best model (determined by LRT)
gr0=lmer(sizeNext~1+(1|year),data=data[!is.na(data$size)&!is.na(data$sizeNext),])
gr1=lmer(sizeNext~size+(1|year),data=data[!is.na(data$size)&!is.na(data$sizeNext),])

anova(gr0,gr1,test="Chisq")

paramCont[[2]]=cbind(as.matrix(coef(gr1)$year),rep(sd(residuals(gr1)),4)) # the third column is for the standard deviation of growth 

# seedling size distribution is not a function of size; save parameters of NULL model
sds0=lmer(sds~1+(1|year),data=data)


paramCont[[3]]=cbind(as.matrix(coef(sds0)$year),rep(sd(residuals(sds0)),4))

# Run and compare mixed model for probability of flowering and save parameters of best model (determined by LRT)
fl0=glmer(fl~1+(1|year),family=binomial,data=data)
fl1=glmer(fl~size+(1|year),family=binomial,data=data)


anova(fl0,fl1, test="Chisq")

paramCont[[4]]=as.matrix(coef(fl1)$year)

# Run and compare mixed model for number of flowering stalks and save parameters of best model (determined by LRT)

fs0=glmer(fs~1+(1|year),family=poisson,data = data)
fs1=glmer(fs~size+(1|year),family=poisson,data = data)

anova(fs0,fs1)

paramCont[[5]]=as.matrix(coef(fs1)$year)

# Run and compare mixed model for number of flowers per stalk and save parameters of best model (determined by LRT)

fps0=glmer(fps~1+(1|year),family=poisson,data = data)
fps1=glmer(fps~size+(1|year),family=poisson,data = data)

anova(fps0,fps1)

paramCont[[6]]=as.matrix(coef(fps1)$year)

###########################################################################
### PART B - IPM and SIMULATIONS BASED ON KERNEL SELECTION
##########################################################################

# Construct an IPM kernel K using the parameters we obtained from the models

  
  # First define the elements that make up the IPM (vital rates):
  
  # SURVIVAL:
  
  S.fun <- function(z,year) {
    
    mu.surv=exp(paramCont[[1]][year,"(Intercept)"] + paramCont[[1]][1,"size"]*z)

    return(mu.surv/(1+mu.surv))
  }
  
  # GROWTH (we assume a contant variance)
  
  GR.fun <- function(z,zz,year){
    
    growth.mu=(paramCont[[2]][year,"(Intercept)"] + paramCont[[2]][1,"size"]*z)
 
    # Get residual variance 
    var.res= (paramCont[[2]][1,3])^2
    # Density distribution function of the normal distribution
    gr1 = sqrt(2*pi*var.res)
    gr2 = ((zz-growth.mu)^2)/(2*var.res)
    
    return(exp(-gr2)/gr1)
    
  }
  
  ## SEEDLING SIZES (same approach as in growth function)
  
  SDS.fun <- function(zz,year){
    
    sds.mu=(paramCont[[3]][year,"(Intercept)"])

    # Get residual variance 
    var.res= (paramCont[[3]][1,2])^2
    # Density distribution function of the normal distribution
    sds1 = sqrt(2*pi*var.res)
    sds2 = ((zz-sds.mu)^2)/(2*var.res)
    
    return(exp(-sds2)/sds1)
    
  }
  
  # PROBABILITY OF FLOWERING 
  
  FL.fun <- function(z,year) {
    
      mu.fl=exp(paramCont[[4]][year,"(Intercept)"] + paramCont[[4]][1,"size"]*z)

    return(mu.fl/(1+mu.fl))
  }
  
  # NUMBER OF FLOWERING STALKS 
  
  FS.fun <- function(z,year) {

    mu.fs=exp(paramCont[[5]][year,"(Intercept)"] + paramCont[[5]][1,"size"]*z)

  }
  
  # NUMBER OF FLOWERS PER STALK
  
  FPS.fun <- function(z,year) {
    
      mu.fps=exp(paramCont[[6]][year,"(Intercept)"] + paramCont[[6]][1,"size"]*z)

    return(mu.fps)
  }
  
  # Second, put together the kernels - four kernels for four years:
  
  # Define the lower and upper integration limit
  L=0.9*min(data$size[!is.na(data$size)]) # minimum size
  U=1.1*max(data$size[!is.na(data$size)]) # maximum size
  
  n=50 # bins
  
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width 
  
  # define years
  years=1:4
  
  # define K
  
  K=array(0,c(n,n,length(years)))
  
  for(i in 1:length(years)){
    
    
    S <- diag(S.fun(meshp,years[i])) # Survival  
    G <- h*t(outer(meshp,meshp,GR.fun,years[i])) # Growth
    
    #Recruits distribution
    R <- h*matrix(rep(SDS.fun(meshp,years[i]),n),n,n,byrow=F)
    
    #Probability of flowering
    Fec01 = (diag(FL.fun(meshp,years[i])))
    
    #Number of Flowering Stalks 
    Fec02 = (diag(FS.fun(meshp,years[i])))
    
    #Number of flowers per Stalk
    
    Fec03= (diag(FPS.fun(meshp,years[i])))
    
    #Number of seeds per flower that survive to become offspring 
    Fec04 = (diag(rep(9.8,n)))
    
    FecALL= Fec01*Fec02*Fec03*Fec04
    
    # Control for eviction:
    # this is equivalent to redistributing evictd sizes evenly among existing size classes 
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    R=R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
    Pkernel <- as.matrix(G%*%S)
    
    Fkernel <- as.matrix(R%*%FecALL)
    
    K[,,i] <-Pkernel+Fkernel
    
  }
  
## Simulate log lambda.s

# Length of simulations
trun=5000

# vector of environments (here random samples from 4 years)
states = sample(4, trun, replace=T)

# vector to hold your 1 step growth rates
growth=rep(NA,trun)

# Initial population vector:

vec1=rep(1,n) 
vec1 <- vec1/sum(vec1) # scale so that the pop size is not unbounded

# Calculate lambda.s for a total of trun simulations: 
for(i in 1:trun){
  # get the environmental state/year at each iteration
  s <- states[i] 
  mat <- K[,,s] # get the K associated with s
  vec1 <- mat%*%as.numeric(vec1) # multiply IPM by vector 
  # get log(lambda.s) for iteration i
  growth[i]  <- log(sum(vec1))
  vec1 <- vec1/sum(vec1) 
  
}
# the stochastic growth rate is the mean growth rate over trun iterations
log.lambda.s = sum(growth)/trun

plot(4500:trun,growth[4500:trun],type="l")

# QUICK EXERCISE:

# What is log.lambda.s if you discard the first 500 growth rates from simulations? Does it change much?
# Why is it advisable to discard the initial simulations? 

###########################################################################
### PART C - IPM and SIMULATIONS BASED ON PARAMETER SELECTION
##########################################################################

# The first step here is to create a variance-covariance matrix 
# of vital rate means - these are the ones we consider to vary between years

# 

# The variances can be obtained from the variance of the random effect for each model
variances=c(as.data.frame(VarCorr(surv1))[,4],as.data.frame(VarCorr(gr1))[1,4],
            as.data.frame(VarCorr(sds0))[1,4],as.data.frame(VarCorr(fl1))[1,4],as.data.frame(VarCorr(fs1))[,4],
            as.data.frame(VarCorr(fps1))[,4])

# The covariances are calculated from the year-specific intercept estimates for each vital rate

covMat=cov(cbind(paramCont[[1]][,1],paramCont[[2]][,1],paramCont[[3]][,1],
             paramCont[[4]][,1],paramCont[[5]][,1],paramCont[[6]][,1]))

diag(covMat)=variances

covMat

# As you can see, the number of flowering stalks and flowers per stalk does not vary among years
# so we can remove these from our matrix

covMat=covMat[1:4,1:4]

# We must also get the mean parametes for each intercept

mean.param=c(fixef(surv1)[1],fixef(gr1)[1],fixef(sds0)[1],fixef(fl1)[1])

# Function to construct an IPM kernel K using the parameters we obtained from the models

# To simulate parameters, we need to install MASS

library(MASS)

# This time, wrap the IPM construction into a function 
IPM.param=function(par.new,n,L,U){
  
  
  # First, define the elements that make up the IPM (vital rates):
  # This is almost the same as the functions for kernel selction
  # with the exception that the intercepts are taken from the new simulated paramCont
  # see below
 
   # SURVIVAL:
  
  S.fun <- function(z,par.new) {
    
    mu.surv=exp(par.new[1] + paramCont[[1]][1,"size"]*z)
    
    return(mu.surv/(1+mu.surv))
  }
  
  # GROWTH (we assume a contant variance)
  
  GR.fun <- function(z,zz,par.new){
    
    growth.mu=(par.new[2] + paramCont[[2]][1,"size"]*z)
    
    # Get residual variance 
    var.res= (paramCont[[2]][1,3])^2
    # Density distribution function of the normal distribution
    gr1 = sqrt(2*pi*var.res)
    gr2 = ((zz-growth.mu)^2)/(2*var.res)
    
    return(exp(-gr2)/gr1)
    
  }
  
  ## SEEDLING SIZES (same approach as in growth function)
  
  SDS.fun <- function(zz,par.new){
    
    sds.mu=(par.new[3])
    
    # Get residual variance 
    var.res= (paramCont[[3]][1,2])^2
    # Density distribution function of the normal distribution
    sds1 = sqrt(2*pi*var.res)
    sds2 = ((zz-sds.mu)^2)/(2*var.res)
    
    return(exp(-sds2)/sds1)
    
  }
  
  # PROBABILITY OF FLOWERING 
  
  FL.fun <- function(z,par.new) {
    
    mu.fl=exp(par.new[4] + paramCont[[4]][1,"size"]*z)
    
    return(mu.fl/(1+mu.fl))
  }
  
  # NUMBER OF FLOWERING STALKS (here, we don't simulate new paramCont)
  
  FS.fun <- function(z) {
    
    mu.fs=exp(paramCont[[5]][1,"(Intercept)"] + paramCont[[5]][1,"size"]*z)
    
  }
  
  # NUMBER OF FLOWERS PER STALK (here, we don't simulate new paramCont)
  
  FPS.fun <- function(z) {
    
    mu.fps=exp(paramCont[[6]][1,"(Intercept)"] + paramCont[[6]][1,"size"]*z)
    
    return(mu.fps)
  }
  
  # Second, put together the kernel for each time:

  # define K
  
  K=array(0,c(n,n))
  
    
    b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
    meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
    
    h=(U-L)/n # bin width 
    
    S <- diag(S.fun(meshp,par.new)) # Survival  
    G <- h*t(outer(meshp,meshp,GR.fun,par.new)) # Growth
    
    #Recruits distribution
    R <- h*matrix(rep(SDS.fun(meshp,par.new),n),n,n,byrow=F)
    
    # Control for eviction:
    # this is equivalent to redistributing evictd sizes evenly among existing size classes 
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    R=R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
    #Probability of flowering
    Fec01 = (diag(FL.fun(meshp,par.new)))
    
    #Number of Flowering Stalks 
    Fec02 = (diag(FS.fun(meshp)))
    
    #Number of flowers per Stalk
    
    Fec03= (diag(FPS.fun(meshp)))
    
    #Number of seeds per flower that survive to become offspring 
    Fec04 = (diag(rep(9.8,n)))
    
    FecALL= Fec01*Fec02*Fec03*Fec04
    
    
    Pkernel <- as.matrix(G%*%S)
    
    Fkernel <- as.matrix(R%*%FecALL)
    
    K <-Pkernel+Fkernel
    
  
  return(K)
}

# Define the lower and upper integration limit
L=0.9*min(data$size[!is.na(data$size)]) # minimum size
U=1.1*max(data$size[!is.na(data$size)]) # maximum size

n=50 # bins

## Simulate log lambda.s

# Length of simulations
trun=5000

# vector to hold your 1 step growth rates
growth=rep(NA,trun)

# Initial population vector:

vec1=rep(1,n) 
vec1 <- vec1/sum(vec1) # scale so that the pop size is not unbounded

# Calculate lambda.s for a total of trun simulations: 
for(i in 1:trun){
  # simulate paramters from a multivariate normal distribution
  par.new=mvrnorm(1, mean.param, covMat)
  
  mat <- IPM.param(par.new,n,L,U) # get the Kt
  vec1 <- mat%*%as.numeric(vec1) # multiply IPM by vector 
  # get log(lambda.s) for iteration i
  growth[i]  <- log(sum(vec1))
  vec1 <- vec1/sum(vec1) 
}
# the stochastic growth rate is the mean growth rate over trun iterations
log.lambda.sV2 = sum(growth)/trun


############# QUICK EXERCISE

# When taking random parameter samples from a Gaussian distribution,
# it is recommended  to bound the range of parameters generated, 
# which we can do in the following way:

# 1.Cholesky factorization of the variance-covariance matrix 

chol.covMat=chol(covMat)

# 2. during simulation: 
#    generating independent, truncated Gaussian deviates (0.001, 0.999 quantiles)
#    and transforming these using the Cholesky factorization

growth=rep(NA,trun)

# Initial population vector:

vec1=rep(1,n) 
vec1 <- vec1/sum(vec1) # scale so that the pop size is not unbounded

# Calculate lambda.s for a total of trun simulations: 
for(i in 1:trun){
  # generate deviates
  par.new=mean.param+matrix(qnorm(runif(4,0.001,0.999)),1,4) %*% chol.covMat
  
  mat <- IPM.param(par.new,n,L,U) # get the Kt
  vec1 <- mat%*%as.numeric(vec1) # multiply IPM by vector 
  # get log(lambda.s) for iteration i
  growth[i]  <- log(sum(vec1))
  vec1 <- vec1/sum(vec1) 
}
# the stochastic growth rate is the mean growth rate over trun iterations
log.lambda.sV3 = sum(growth)/trun

### Do log.lambda.s, log.lambda.sV2, log.lambda.sV3 differ?
### How can you decrease the differences? 
