# Following the results of stage 6(viral load heterogeneity) code 'Bayesian_IBM6.b.R'
# I now want to write up the entire model with comparable results from stage 1 to
# the final level of viral load. Since the data age age groups have changed over
# the growth of the model, I now want to rerun all the levels at once. Also the 
# definition of the stages will change a bit.
# 
# Stage 1: Basic model with eta and epsilon
# Stage 2: Heterogeneity in susceptibility by infection history
# Stage 3: Heterogeneity in susceptibility by age group
# Stage 4: Heterogeneity in exposure by age
# Stage 5: Household structure
#           a) Contact matrices by age
# Stage 6: Effects of HH size 
#           a) Density dependence Vs Frequency dependence transmission
#           b) Transmission in large Vs small HHs
# Stage 7: Heterogeneity in infectivity by viral load
#           a) Continuous
#           b) Categorical
# 
# 18th July 2017
################################################################################

#-------------------------------------------------------------------------------
# Stage 1
#-------------------------------------------------------------------------------
# Likelihood
IBM.LL.1<-function(theta){
    # 1) Relative susceptibility
    PI<-matrix(1,D,N,byrow=FALSE)
    PI[which(shedding>0)]=0
    # 1) Rate of exposure
      # HH
      eta<-exp(theta[1]);HH.risk<-t(Cc %*% t(shedding));in.house<-(status)*(eta*HH.risk) 
      # Community
      epsilon<-exp(theta[2]);Comm <- epsilon*matrix(comm.risk,D,N)
      # Total 
      Lambda<-in.house + Comm
    # 2) Probability of Infection
      Alpha<-(1-exp(-PI*Lambda))
    # 3) Probability of Starting Shedding
      PrShed<-matrix(0,D,N)
      for (t in 1:D){
        # for every day,t
        days<-t-seq(0,length(Dist.lat)-1)
        PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
      }
    # Likelihood for days of observed data, ie onset and at risk day(ie present and
    # not shedding)
      nPrShed<-1-PrShed;pple<-seq(N)
      logger<-function(p){
        return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
      }
      Like<-sapply(pple, logger)
    return(sum(Like))
}
# Prior
IBM.prior.1<-function(theta){
    Eta<-dunif(theta[1], min = -20, max = 0, log = TRUE)
    Epsilon<-dunif(theta[2], min = -20, max = 0, log = TRUE)
    ## the joint prior
    log.sum<-Eta+Epsilon
    ## default is to give the prior but otherwise give the log prior
    return(log.sum)
}
# Posterior
IBM.posterior.1<-function(theta){
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.1(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.1(theta)
  # return the logged posterior probability
  return(x+y)
}

#-------------------------------------------------------------------------------
# Stage 2
#-------------------------------------------------------------------------------
# Likelihood
IBM.LL.2<-function(theta){
    # 1) Relative susceptibility
    # Infection history component
    hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
    hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
    PI<-exp(hist.mat)
    PI[which(shedding>0)]=0
    # 2) Rate of exposure
    # HH
    eta<-exp(theta[2]);HH.risk<-t(Cc %*% t(shedding));in.house<-(status)*(eta*HH.risk) 
    # Community
    epsilon<-exp(theta[3]);Comm <- epsilon*matrix(comm.risk,D,N)
    # Total 
    Lambda<-in.house + Comm
    # 3) Probability of Infection
    Alpha<-(1-exp(-PI*Lambda))
    # 4) Probability of Starting Shedding
    PrShed<-matrix(0,D,N)
    for (t in 1:D){
      # for every day,t
      days<-t-seq(0,length(Dist.lat)-1)
      PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
    }
    # Likelihood for days of observed data, ie onset and at risk day(ie present and
    # not shedding)
    nPrShed<-1-PrShed;pple<-seq(N)
    logger<-function(p){
      return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
    }
    Like<-sapply(pple, logger)
    return(sum(Like))
}
# Prior
IBM.prior.2<-function(theta){
  inf.hist<-dunif(theta[1], min = -10, max = 10, log = TRUE)
  Eta<-dunif(theta[2], min = -20, max = 0, log = TRUE)
  Epsilon<-dunif(theta[3], min = -20, max = 0, log = TRUE)
  ## the joint prior
  log.sum<-inf.hist+Eta+Epsilon
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.2<-function(theta){
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.2(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.2(theta)
  # return the logged posterior probability
  return(x+y)
}

#-------------------------------------------------------------------------------
# Stage 3
#-------------------------------------------------------------------------------
# Likelihood
IBM.LL.3<-function(theta){
  # 1) Relative susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # 2) Rate of exposure
  # HH
  eta<-exp(theta[s.g+1]);HH.risk<-t(Cc %*% t(shedding));in.house<-(status)*(eta*HH.risk) 
  # Community
  epsilon<-exp(theta[s.g+2]);Comm <- epsilon*matrix(comm.risk,D,N)
  # Total 
  Lambda<-in.house + Comm
  # 3) Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # 4) Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior
IBM.prior.3<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[s.g+2], min = -20, max = 0, log = TRUE)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+Epsilon
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.3<-function(theta){
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.3(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.3(theta)
  # return the logged posterior probability
  return(x+y)
}
#-------------------------------------------------------------------------------
# Stage 4
#-------------------------------------------------------------------------------
# Likelihood
IBM.LL.4<-function(theta){
  # _________________1. Relative Susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # _________________2. Rate of Exposure
  # Within HH
  # General prob of in house transmission
  Eta<-exp(theta[s.g+1]) 
  HH.risk<-t(Cc %*% t(shedding))
  # Total within household risk
  in.house<-status*(Eta*HH.risk) 
  # Community
  x<-length(theta)-e.g
  Eps<-exp(theta[(x+1)])
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N);comm.time<-matrix(comm.risk,D,1)
  Comm <- Eps*(comm.time %*% comm.age)  # Time varying community exposure
  # Total exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior
IBM.prior.4<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  x<-length(theta)-e.g
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[(x+1)], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.4<-function(theta){
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.4(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.4(theta)
  # return the logged posterior probability
  return(x+y)
}
#-------------------------------------------------------------------------------
# Stage 5
#-------------------------------------------------------------------------------
# a) Contact matrices by age
#----------------------------
# Likelihood
IBM.LL.5<-function(theta){
  # _________________1. Relative Susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # _________________2. Rate of Exposure
  # Within HH
  # General prob of in house transmission
  Eta<-exp(theta[s.g+1]) 
  x<-length(theta)-e.g
  # Effect of different contact types
  con.par<-c(0,1,exp(theta[(s.g+2):(x)])) 
  Cij<-matrix(con.par[C],N,N,byrow=FALSE);diag(Cij)<-0 
  HH.risk<-t(Cij %*% t(shedding)) 
  # Total within household risk
  in.house<-status*(Eta*HH.risk) 
  # Community
  Eps<-exp(theta[(x+1)])
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N);comm.time<-matrix(comm.risk,D,1)
  Comm <- Eps*(comm.time %*% comm.age)  # Time varying community exposure
  # Total exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior
IBM.prior.5<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  x<-length(theta)-e.g
  ## contact type parameters
  con.par<-sapply(theta[(s.g+2):(x)],sus)
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[(x+1)], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+sum(con.par)+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.5<-function(theta){
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.5(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.5(theta)
  # return the logged posterior probability
  return(x+y)
}
#-------------------------------------------------------------------------------
# Stage 6
#-------------------------------------------------------------------------------
# a) FD vs DD
#----------------------------

# Likelihood
IBM.LL.6<-function(theta){
  # _________________1. Relative Susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # _________________2. Rate of Exposure
  # Within HH
  # General prob of in house transmission
  Eta<-exp(theta[s.g+1]) 
  HH.risk<-t(Cc %*% t(shedding))
  # FD vs DD
  x<-length(theta)-e.g
  omega<-round(exp(theta[x]),3)
  hh.eff<-matrix(0,D,N)
  hh.eff<-(Nh-1)^(-omega) 
  in.house<-status* Eta*hh.eff*HH.risk 
  # Community
  Eps<-exp(theta[(x+1)])
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N);comm.time<-matrix(comm.risk,D,1)
  Comm <- Eps*(comm.time %*% comm.age)  # Time varying community exposure
  # Total exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior
IBM.prior.6<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  x<-length(theta)-e.g
  omega<-sapply(theta[x],sus)
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[(x+1)], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+omega+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.6<-function(theta){
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.6(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.6(theta)
  # return the logged posterior probability
  return(x+y)
}

# b) Categorical HH size
#----------------------------

# Likelihood
IBM.LL.6.b<-function(theta){
  # _________________1. Relative Susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # _________________2. Rate of Exposure
  # Within HH
  # General prob of in house transmission
  Eta<-exp(theta[s.g+1]) 
  # Effect of HH size
  hh.par<-c(1,exp(theta[s.g+2]));hh.mat<-matrix(hh.par[H],D,N,byrow=TRUE)
  # Mapping viral load into infectivity
  x<-length(theta)-e.g
  HH.risk<-t(Cc %*% t(shedding))
  # Total within household risk
  in.house<-(status*hh.mat)*(Eta*HH.risk) 
  # Community
  Eps<-exp(theta[(x+1)])
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N);comm.time<-matrix(comm.risk,D,1)
  Comm <- Eps*(comm.time %*% comm.age)  # Time varying community exposure
  # Total exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior
IBM.prior.6.b<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  ## uniform prior on household size covariate
  hh.size<-sus(theta[s.g+2])
  x<-length(theta)-e.g
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[(x+1)], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+hh.size+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.6.b <- function(theta) {
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.6.b(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.6.b(theta)
  # return the logged posterior probability
  return(x+y)
}
#-------------------------------------------------------------------------------
# Stage 7
#-------------------------------------------------------------------------------
# a) Continuous Viral load
#----------------------------

# Likelihood function
IBM.LL.7.a<-function(theta){
  # _________________1. Relative Susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # _________________2. Rate of Exposure
  # Within HH
  # General prob of in house transmission
  Eta<-exp(theta[s.g+1]) 
  # Effect of HH size
  hh.par<-c(1,exp(theta[s.g+2]));hh.mat<-matrix(hh.par[H],D,N,byrow=TRUE)
  # Mapping viral load into infectivity
  x<-length(theta)-e.g
  if(max(shedding)>1){
    Vpar<-theta[(s.g+3):x]
    # 1st mapping
    V.par<-exp(c(theta[(s.g+3):x]))
    Infectivity<-V.map(shedding,V.par) 
    HH.risk<-t(Cc %*% t(Infectivity))
  }else{
    HH.risk<-t(Cc %*% t(shedding))
  }
  # Total within household risk
  in.house<-(status*hh.mat)*(Eta*HH.risk) 
  # Community
  Eps<-exp(theta[(x+1)])
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N);comm.time<-matrix(comm.risk,D,1)
  Comm <- Eps*(comm.time %*% comm.age)  # Time varying community exposure
  # Total exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior function
IBM.prior.7.a<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  ## uniform prior on household size covariate
  hh.size<-sus(theta[s.g+2])
  # Viral mapping function parameters
  x<-length(theta)-e.g
  if(max(shedding)>1){
    V.par<-theta[(s.g+3):x]
    # 1st mapping
    if(names(V.par)[1]=="beta"){
      v.pars<-c(dunif(V.par[1],min=-5,max=5,log=TRUE),
                dunif(V.par[2],min=-10,max=5,log=TRUE),
                dunif(V.par[3],min=-10,max=5,log=TRUE))
    }
  }else{
    v.pars=0
  }
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[(x+1)], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+hh.size+sum(v.pars)+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior function
IBM.posterior.7.a <- function(theta) {
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.7.a(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.7.a(theta)
  # return the logged posterior probability
  return(x+y)
}

# b) Categorical Viral load
#----------------------------

IBM.LL.7.b<-function(theta){
  # _________________1. Relative Susceptibility
  # Infection history component
  hist.par<-c(0,theta[1]) # i.e no modifying factor if no previous infection
  hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Age component
  age.par<-c(0,theta[2:s.g])
  age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  PI<-exp(hist.mat+age.mat)
  PI[which(shedding>0)]=0
  # _________________2. Rate of Exposure
  # Within HH
  # General prob of in house transmission
  Eta<-exp(theta[s.g+1]) 
  # Effect of HH size
  hh.par<-c(1,exp(theta[s.g+2]));hh.mat<-matrix(hh.par[H],D,N,byrow=TRUE)
  # Mapping viral load into infectivity
  x<-length(theta)-e.g
  if(max(shedding)>1){
    V.par<-c(0,1,exp(theta[x]))
    Infectivity<-matrix(V.par[Shedding.cat],D,N,byrow=FALSE)
    HH.risk<-t(Cc %*% t(Infectivity))
  }else{
    HH.risk<-t(Cc %*% t(shedding))
  }
  # Total within household risk
  in.house<-(status*hh.mat)*(Eta*HH.risk) 
  # Community
  Eps<-exp(theta[(x+1)])
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N);comm.time<-matrix(comm.risk,D,1)
  Comm <- Eps*(comm.time %*% comm.age)  # Time varying community exposure
  # Total exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding
  PrShed<-matrix(0,D,N)
  for (t in 1:D){
    # for every day,t
    days<-t-seq(0,length(Dist.lat)-1)
    PrShed[t,]<-Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],]
  }
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# Prior
IBM.prior.7.b<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[1])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[s.g+1], min = -20, max = 0, log = TRUE)
  ## uniform prior on household size covariate
  hh.size<-sus(theta[s.g+2])
  # Viral mapping function parameters
  x<-length(theta)-e.g
  if(max(shedding)>1){
    v.pars<-sus(theta[x])
  }else{
    v.pars=0
  }
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[(x+1)], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+hh.size+sum(v.pars)+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# Posterior
IBM.posterior.7.b <- function(theta) {
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.7.b(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.7.b(theta)
  # return the logged posterior probability
  return(x+y)
}


