## IPM course Barcelona November 2017
### R code by Maria Paniw
### Exercise on environmental stochasticity

## This exercise is intended to show you how to construct IPMs
## using mixed effect models for vital rates 
## and how to incorporate dependent ensvironmental states (TSF changes)
## and random environmental variation (random effect in vital rate models) 
## into stochastic population projections

## NOTE: The IPMs now include a discrete stage - the seed bank
##       for more detail, see Paniw et al. 2016, Oikos

### This exercise consists of 3 parts

### Part A: Vital rate functions 

### Part B: Construct IPM kernels 

### PART C: simulate log lambda.s using kernel selection

###########################################################################
### PART A - VITAL RATE FUNCTIONS
##########################################################################

# set working directory (change to appropriate one)
setwd("/Users/mariapaniw/Dropbox/teaching/IPM2017")

data=read.csv("dataDroso.csv")

head(data)

data$year=as.factor(data$year)

# Order the levels of TSF categories
data$TSF=factor(data$TSF,levels = c("one", "two", "three",">three"))

# load lme4 package for mixed effect models
library(lme4)

#### Vital rate models

# Empty list to save model coefficients for continuous transitions in
paramCont=list(NULL)

# Run and compare mixed model for survival and save parameters of best model (determined by LRT)
surv0=glmer(surv~1+(1|year),family=binomial,data=data)
surv1=glmer(surv~size+(1|year),family=binomial,data=data)
surv2=glmer(surv~size+TSF+(1|year),family=binomial,data=data)
surv3=glmer(surv~size+TSF+size*TSF+(1|year),family=binomial,data=data)

anova(surv0,surv1,surv2,surv3, test="Chisq")

paramCont[[1]]=as.matrix(coef(surv2)$year)

# Run and compare mixed model for growth and save parameters of best model (determined by LRT)
gr0=lmer(sizeNext~1+(1|year),data=data[!is.na(data$size)&!is.na(data$sizeNext),])
gr1=lmer(sizeNext~size+(1|year),data=data[!is.na(data$size)&!is.na(data$sizeNext),])
gr2=lmer(sizeNext~size+TSF+(1|year),data=data[!is.na(data$size)&!is.na(data$sizeNext),])
gr3=lmer(sizeNext~size+TSF+size*TSF+(1|year),data=data[!is.na(data$size)&!is.na(data$sizeNext),])

anova(gr0,gr1,gr2,gr3,test="Chisq")

paramCont[[2]]=cbind(as.matrix(coef(gr3)$year),rep(sd(residuals(gr3)),4))

# Run and compare mixed model for seedling size distribution and save parameters of best model (determined by LRT)
sds0=lmer(sds~1+(1|year),data=data)
sds1=lmer(sds~TSF+(1|year),data=data)

anova(sds0,sds1,test="Chisq")

paramCont[[3]]=cbind(as.matrix(coef(sds1)$year),rep(sd(residuals(sds1)),4))

# Run and compare mixed model for probability of flowering and save parameters of best model (determined by LRT)
fl0=glmer(fl~1+(1|year),family=binomial,data=data)
fl1=glmer(fl~size+(1|year),family=binomial,data=data)
fl2=glmer(fl~size+TSF+(1|year),family=binomial,data=data)
fl3=glmer(fl~size+TSF+size*TSF+(1|year),family=binomial,data=data)

anova(fl0,fl1,fl2,fl3, test="Chisq")

paramCont[[4]]=as.matrix(coef(fl3)$year)

# Run and compare mixed model for number of flowering stalks and save parameters of best model (determined by LRT)
fs0=glmer(fs~1+(1|year),family=poisson,data = data)
fs1=glmer(fs~size+(1|year),family=poisson,data = data)
fs2=glmer(fs~size+TSF+(1|year),family=poisson,data = data)
fs3=glmer(fs~size+TSF+size*TSF+(1|year),family=poisson,data = data)

anova(fs0,fs1,fs2,fs3, test="Chisq")

paramCont[[5]]=as.matrix(coef(fs2)$year)

# Run and compare mixed model for number of flowers per stalk and save parameters of best model (determined by LRT)

fps0=glmer(fps~1+(1|year),family=poisson,data = data)
fps1=glmer(fps~size+(1|year),family=poisson,data = data)
fps2=glmer(fps~size+TSF+(1|year),family=poisson,data = data)
fps3=glmer(fps~size+TSF+size*TSF+(1|year),family=poisson,data = data)

anova(fps0,fps1,fps2,fps3, test="Chisq")

paramCont[[6]]=as.matrix(coef(fps2)$year)

###########################################################################
### PART B - IPMs - for each TSF and year
##########################################################################

# Function to construct an IPM kernel K using the parameters we obtained from the models

### NOTICE that the IPMs now contain a simplified discrete stage: the seed bank

# The seed-bank parameters are

# goSB - proportion of seeds going into the seed bank vs. going to continuous stages (goCont; recruits)
# staySB - proportion of seeds surviving and staying in the seedbank
# outSB - proportion of seed surviving, leaving the seedbank and estabishing as recruits

########################### You can define the arguments for the functions later, here I do it in case you want to run function line by line 
# Define the lower and upper integration limit
L=0.9*min(data$size[!is.na(data$size)]) # minimum size
U=1.1*max(data$size[!is.na(data$size)]) # maximum size

n=50 # bins
parameters=paramCont
###########################

IPM.kernel=function(parameters,n,L,U){
  
  # First define the elements that make up the IPM (vital rates):
  # Because we are fitting ANCOVA models and using the default R effects parameterization,
  # extracting the parameters can be a bit tricky
  
  # SURVIVAL:
  
  S.fun <- function(z,tsf,year) {
    
    # For time since fire 1, vital rates are modeled as just the intercept + size
    if(tsf=="one"){mu.surv=exp(paramCont[[1]][year,"(Intercept)"] + paramCont[[1]][1,"size"]*z)}else{
      mu.surv=exp(paramCont[[1]][year,"(Intercept)"] + paramCont[[1]][1,"size"]*z + paramCont[[1]][1,tsf] )
    }
    
    
    return(mu.surv/(1+mu.surv))
  }
  
  # GROWTH (we assume a contant variance)
  
  GR.fun <- function(z,zz,tsf,year){
    
    if(tsf=="one"){growth.mu=(paramCont[[2]][year,"(Intercept)"] + paramCont[[2]][1,"size"]*z)}else{
      growth.mu=(paramCont[[2]][year,"(Intercept)"] + paramCont[[2]][1,"size"]*z + paramCont[[2]][1,tsf] + paramCont[[2]][1,paste("size:",tsf,sep="")]*z)
    }
    
    # Get residual variance 
    var.res= (paramCont[[2]][1,9])^2
    # Density distribution function of the normal distribution
    gr1 = sqrt(2*pi*var.res)
    gr2 = ((zz-growth.mu)^2)/(2*var.res)
    
    return(exp(-gr2)/gr1)
    
  }
  
  ## SEEDLING SIZES (same approach as in growth function)
  
  SDS.fun <- function(zz,tsf,year){
    
    if(tsf=="one"){sds.mu=(paramCont[[3]][year,"(Intercept)"])}else{
      sds.mu=(paramCont[[3]][year,"(Intercept)"] + paramCont[[3]][1,tsf])
    }
    
    # Get residual variance 
    var.res= (paramCont[[3]][1,5])^2
    # Density distribution function of the normal distribution
    sds1 = sqrt(2*pi*var.res)
    sds2 = ((zz-sds.mu)^2)/(2*var.res)
    
    return(exp(-sds2)/sds1)
    
  }
  
  # PROBABILITY OF FLOWERING 
  
  FL.fun <- function(z,tsf,year) {
    
    if(tsf=="TSFtwo"){mu.fl=exp(paramCont[[4]][year,"(Intercept)"] + paramCont[[4]][1,"size"]*z)}else{
      mu.fl=exp(paramCont[[4]][year,"(Intercept)"] + paramCont[[4]][1,"size"]*z + paramCont[[4]][1,tsf] +paramCont[[4]][1,paste("size:",tsf,sep="")]*z )
    }
    
    
    return(mu.fl/(1+mu.fl))
  }
  
  # NUMBER OF FLOWERING STALKS
  
  FS.fun <- function(z,tsf,year) {
    
    if(tsf=="TSFtwo"){mu.fs=exp(paramCont[[5]][year,"(Intercept)"] + paramCont[[5]][1,"size"]*z)}else{
      mu.fs=exp(paramCont[[5]][year,"(Intercept)"] + paramCont[[5]][1,"size"]*z + paramCont[[5]][1,tsf] )
    }
    
  }
  
  # NUMBER OF FLOWERS PER STALK
  
  FPS.fun <- function(z,tsf,year) {
    
    if(tsf=="TSFtwo"){mu.fps=exp(paramCont[[6]][year,"(Intercept)"] + paramCont[[6]][1,"size"]*z)}else{
      mu.fps=exp(paramCont[[6]][year,"(Intercept)"] + paramCont[[6]][1,"size"]*z + paramCont[[6]][1,tsf] )
    }
    return(mu.fps)
  }
  
  # Second, put together the kernels:
  
  # define TSF categories
  tsf.cat=c("zero","one","TSFtwo","TSFthree","TSF>three")
  
  # define years
  years=1:4
  
  # define the bins
  
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width 
  
  # define K
  
  K=array(0,c(n+1,n+1,length(tsf.cat),length(years)))
  
  for(t in 1:length(tsf.cat)){ # loop through TSF categories
    tsf=tsf.cat[t]
    
    for(i in 1:length(years)){ # loop through 4 years per TSF category
      
      ### FOR BURNING (all above-ground plants die, germination from seed bank occurs)
      
      if(tsf=="zero"){
        outSB <- 0.6 
        staySB <- 0.1 
        goCont=0.01 # so that we don?t devide by zero! (see below)
        goSB=0
        S <-matrix(0,n,n) # survival, growth (G), and fecundity (FecALL) are all zero
        G <- matrix(0,n,n)
        FecALL=matrix(0,n,n)
        R <- h* matrix(rep(SDS.fun(meshp,tsf="one",years[i]),n),n,n,byrow=F)# the relevant non-0 transition is the size distribution of seedlings
        
        # Control for eviction:
        # this is equivalent to redistributing evictd sizes evenly among existing size classes 
        R=R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        
      }else if(tsf=="one"){
        
        ### FOR TSF 1 (reproduction is still 0, but growth occurs)
        outSB <- 0.1 
        staySB <- 0.05 
        goCont=0.01 # so that we don?t devide by zero! (see below)
        goSB=0
        S <- diag(S.fun(meshp,tsf=tsf,years[i])) # Survival  
        G <- h*t(outer(meshp,meshp,GR.fun,tsf=tsf,years[i])) # Growth
        R <- h*matrix(rep(SDS.fun(meshp,tsf=tsf,years[i]),n),n,n,byrow=F)
        
        # Control for eviction:
        # this is equivalent to redistributing evictd sizes evenly among existing size classes 
        G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        R=R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        
        FecALL=matrix(0,n,n)
        
      
      }else{
        
        # for TSF >1
        outSB <- 0.01 
        staySB <- 0.6 
        goCont=0.05 
        goSB=0.45
        
        S <- diag(S.fun(meshp,tsf=tsf,years[i])) # Survival  
        G <- h*t(outer(meshp,meshp,GR.fun,tsf=tsf,years[i])) # Growth
        
        #Recruits distribution
        R <- h*matrix(rep(SDS.fun(meshp,tsf=tsf,years[i]),n),n,n,byrow=F)
        
        # Control for eviction:
        # this is equivalent to redistributing evictd sizes evenly among existing size classes 
        G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        R=R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        
        #Probability of flowering
        Fec01 = (diag(FL.fun(meshp,tsf=tsf,years[i])))
        
        #Number of Flowering Stalks 
        Fec02 = (diag(FS.fun(meshp,tsf=tsf,years[i])))
        
        #Number of flowers per Stalk
        
        Fec03= (diag(FPS.fun(meshp,tsf=tsf,years[i])))
        
        #Number of seeds per flower that survive to become offspring 
        Fec04 = (diag(rep(9,n)))
        
        FecALL= Fec01*Fec02*Fec03*Fec04*goCont
        
        
      }
      
      
      # construct the kernel per TSF category & year
      
      Pkernel.cont <- as.matrix(G%*%S) # continuous transitions
      Pkernel.discr = c(staySB,outSB*R[,2]/sum(R[,2]))
      
      Pkernel <- cbind(Pkernel.discr,rbind(rep(0,length(meshp)),Pkernel.cont))# discrete component
      
      Fkernel.cont <- as.matrix(R%*%FecALL)
      Fkernel.discr =rep(0,length(meshp)+1)
      
      Fkernel <- cbind(Fkernel.discr,rbind(diag(FecALL)*goSB/goCont,Fkernel.cont))
      
      K[,,t,i] <- Pkernel+Fkernel
      
    }
  }
  
  
  return(K)
}

K=IPM.kernel(paramCont,n,L,U)

## Simulate log lambda.s

# Length of projections
trun=5000

### Create an environmental transition matrix
# based on autocorrelation coefficient and frequency of burning

f=0.1
v1=-0.05

q=f*(1-v1)
p=v1+q

envS=matrix(c(p,1-p,rep(0,3),
           q,0,1-q,rep(0,2),
           q,0,0,1-q,0,
           q,0,0,0,1-q,
           q,0,0,0,1-q),nrow=5,ncol=5,byrow=F)


#### SIMULATIONS

# Unlike in Exercise 1, we now want to repeat each projection of 5000 years 
# 100 times. This is advisable if you want to calculate the variance of stochastic lambda 

nsim=100

# Create log.lambda.s vector

log.lambda.s=numeric(nsim)

# Loop through the 100 simulations and safe log.lambda.s for each

for(sim in 1:nsim){
  
  # vector of environments
  states = numeric(trun)
  
  # start with burning 
  
  env.at.t=1
  states[1]=1
  
  for(x in 2:trun)
    states[x] = env.at.t =sample(ncol(envS),1,pr=envS[,env.at.t])
  
  ### vector of random years
  
  years=sample(4, trun, replace=T)
  
  # vector to hold your 1 step growth rates
  growth=rep(NA,trun)
  
  # Initial population vector:
  
  vec1=c(1000,rep(1,n)) # start with 1000 seeds in the seed bank 
  vec1 <- vec1/sum(vec1) # scale so that the pop size is not unbounded
  
  # Calculate lambda.s for a total of trun simulations: 
  for(i in 1:trun){
    # get the environmental state/year at each iteration
    s <- states[i] 
    y <- years[i]
    mat <- K[,,s,y] # get the K associated with s
    vec1 <- mat%*%as.numeric(vec1) # multiply IPM by vector 
    # get log(lambda.s) for iteration i
    growth[i]  <- log(sum(vec1))
    vec1 <- vec1/sum(vec1) 
  }
  # the stochastic growth rate is the mean growth rate over trun iterations
  log.lambda.s[sim] = sum(growth)/trun
  
}

hist(log.lambda.s)

### Viability

# We may want to know if after a certain number of years (here 100 y), 
# the population is quasi-extinct (< 100 seeds in the seed bank in post fire >3)

nsim=100
trun=100 # now we are only projecting 100 years! (Transients)
# Create log.lambda.s vector

extinction=rep(0,nsim)

# Loop through the 100 simulations and safe log.lambda.s for each

for(sim in 1:nsim){
  
  # vector of environments
  states = numeric(trun)
  
  # start with burning 
  
  env.at.t=1
  states[1]=1
  
  for(x in 2:trun)
    states[x] = env.at.t =sample(ncol(envS),1,pr=envS[,env.at.t])
  
  ### vector of random years
  
  years=sample(4, trun, replace=T)
  
  # now in stead of growth rate, you want to hold the sum of seeds in the seed bank
  n.seeds=rep(NA,trun)
  
  # Initial population vector:
  
  vec1=c(1000,rep(1,n)) # start with 1000 seeds in the seed bank, no scaling 

  # Calculate lambda.s for a total of trun simulations: 
  for(i in 1:trun){
    # get the environmental state/year at each iteration
    s <- states[i] 
    y <- years[i]
    mat <- K[,,s,y] # get the K associated with s
    vec1 <- mat%*%as.numeric(vec1) # multiply IPM by vector 
    # get log(lambda.s) for iteration i
    n.seeds[i]  <- vec1[1]
     
  }
  # quasi extinction
  
  if(any(n.seeds[states==5]<100)) extinction[sim] <- 1
  
}

mean(extinction)
mean(extinction)*(1-mean(extinction)) # variability 
# QUICK EXERCISE:

# What happens to log.lambda.s if you decrease f to 0.01 and change v1 to -0.005
