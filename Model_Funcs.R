# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 13th September 2018
################################################################################
# This code contains functions for fitting the HH transmission model either for 
# a single pathogen (SECTION 1) or for two interacting pathogens (SECTION 2). 
################################################################################
# SECTION 1 (12 parameters)
# For a single pathogen the following parameters are estimated:
# a) Effect of previous infection on susceptibility
# b) Effect of age on susceptibility for age groups 1-5, 5-15 and >15 relative to <1
# c) Within household transmission coefficient
# d) Effect of household size for HHs with >=8 people relative to <8 
# e) Effect of the following combinations of viral load and symptoms:
#    High and Asymptomatic, low and symptomatic, high and symptomatic, relative
#    to low and asymptomatic
# f) Effect of on community exposure for age groups 1-5, >5 relative to <1
# g) Community transmission coefficient for each pathogen

# --- Likelihood function
IBM.LL<-function(theta){
  # _________________1. Relative Susceptibility ________________________________
  # Effect of previous infection on susceptibility
  hist.par<-c(0,theta[['Prev.inf']]);hist.mat<-matrix(hist.par[B],D,N,byrow=FALSE)
  # Effect of age on susceptibility
  age.par<-c(0,theta[2:s.g]);age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  # Combining age and infection history effect on susceptibility
  PI<-exp(hist.mat+age.mat)
  # If individuals are shedding then they are not susceptible
  PI[which(Shedding>0)]<-0
  # _________________2. Rate of Exposure _______________________________________
  # Within household rate of exposure
  # ---------------------------------
  # Within household transmission coefficient 
  Eta<-exp(theta[['eta']]) 
  # Effect of household size on within household transmission
  hh.par<-c(1,exp(theta[['hh.size']]));hh.mat<-matrix(hh.par[H],D,N,byrow=TRUE)
  # Effect of high viral load on infectivity
  V.par<-c(0,1,exp(theta[['HighAsym']]),exp(theta[['LowSym']]),exp(theta[['HighSym']]))
  Infectivity<-matrix(V.par[Load.n.sym],D,N,byrow=FALSE)
  HH.risk<-t(Cc %*% t(Infectivity))
  # Total within household risk, status identifies if an indiviual is present in the household
  in.house<-(status*hh.mat)*(Eta*HH.risk) 
  # Community rate of exposure
  # --------------------------
  # Community transmission coefficient
  Eps<-exp(theta[['epsilon']])
  # Effect of age on community level exposure
  x<-length(theta)-e.g
  Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N)
  # Time varying component of community exposure/background density funtion
  comm.time<-matrix(Comm.risk,D,1)
  # Combining the age effects and time varying function
  Comm <- Eps*(comm.time %*% comm.age)  
  # Total community exposure rate 
  Lambda<-in.house + Comm
  # _________________3. Probability of Infection _______________________________
  Alpha<-(1-exp(-PI*Lambda))
  # _________________4. Probability of Starting Shedding _______________________
  # Function that calculates the probability of onset given latency
  shedder<-function(t){
    days<-t-seq(0,length(Dist.lat)-1)
    return(Dist.lat[which(days>0)]  %*% Alpha[days[which(days>0)],])
  }
  PrShed<-t(sapply(seq(D),shedder))
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed<-1-PrShed;pple<-seq(N)
  logger<-function(p){
    return(sum(log(PrShed[(Onset[,p]==1),p])) + sum(log(nPrShed[(atRisk[,p]==1),p])))
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# --- Prior distribution
# Limits of the prior distribution for each parameters
prior.limits<-data.frame(lower=c(Prev.hom=-10,age.2=-10,
                                 age.3=-10,age.4=-10,eta=-20,hh.size=-10,
                                 HighAsym=-10,LowSym=-10,HighSym=-10,epsilon=-20,
                                 eps.age2=-10,eps.age3=-10),
                         upper=c(Prev.A=10,age.2=10,age.3=10,age.4=10,eta=0, 
                                 hh.size=10,HighAsym=10,LowSym=10,HighSym=10,
                                 epsilon=0,eps.age2=10,eps.age3=10))

IBM.prior<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <- sus(theta[['Prev.inf']])
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[2:s.g],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta<-dunif(theta[['eta']], min = -20, max = 0, log = TRUE)
  ## uniform prior on household size covariate
  hh.size<-sus(theta[['hh.size']])
  # Infectivity mapping function parameters
  x<-length(theta)-e.g
  Inf.pars<-c(sus(theta[['HighAsym']]),sus(theta[['LowSym']]),sus(theta[['HighSym']]))
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon<-dunif(theta[['epsilon']], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community exposure with age
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-inf.hist+sum(age.sus)+Eta+hh.size+sum(Inf.pars)+Epsilon+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# --- Posterior distribution
IBM.posterior <- function(theta) {
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL(theta)
  # return the logged posterior probability
  return(x+y)
}
################################################################################
# SECTION 2 (15 parameters)
# For interacting pathogens the following parameters are estimated:
# a) Effect of previous homologous infection on susceptibility
# b) Effect of previous heterologous infection on susceptibility
# c) Effect of age on susceptibility for age groups 1-5, 5-15 and >15 relative to <1
# d) Within household transmission coefficient for each pathogen
# e) Effect of household size for HHs with >=8 people relative to <8 
# f) Effect of the following combinations of viral load and symptoms:
#    High and Asymptomatic, low and symptomatic, high and symptomatic, relative
#    to low and asymptomatic
# g) Effect of on community exposure for age groups 1-5, >5 relative to <1
# h) Community transmission coefficient for each pathogen

# --- Likelihood function
IBM.LL.group<-function(theta){
  # _________________1. Relative Susceptibility_________________________________
  # Effect of previous homologous and heterologous infection on susceptibility
  hist.par.A<-c(0,theta[['Prev.hom']], theta[['Prev.het']],(theta[['Prev.hom']]+theta[['Prev.het']]))
  hist.mat.A<-matrix(hist.par.A[B],D,N,byrow=FALSE)
  hist.par.B<-c(0,theta[['Prev.het']], theta[['Prev.hom']],(theta[['Prev.hom']]+theta[['Prev.het']]))
  hist.mat.B<-matrix(hist.par.B[B],D,N,byrow=FALSE)
  # Effect of age on susceptibility
  age.par<-c(0,theta[3:(2+s.g-1)]);age.mat<-matrix(age.par[A],D,N,byrow=TRUE)
  # Combining age and infection history effect on susceptibility
  PI.A<-exp(hist.mat.A+age.mat);PI.A[which(Shedding.a.est.load>0)]=0
  PI.B<-exp(hist.mat.B+age.mat);PI.B[which(Shedding.b.est.load>0)]=0
  # _________________2. Rate of Exposure________________________________________
  # Within household rate of exposure
  # ---------------------------------
  # Within household transmission coefficient
  Eta.A<-exp(theta[['eta.A']]); Eta.B<-exp(theta[['eta.B']])
  # Effect of household size on within household transmission
  hh.par<-c(1,exp(theta[['hh.size']]));hh.mat<-matrix(hh.par[H],D,N,byrow=TRUE)
  # Effect of viral load and symptoms on infectivity
  V.par<-c(0,1,exp(theta[['HighAsym']]),exp(theta[['LowSym']]),exp(theta[['HighSym']]))
  Infectivity.A<-matrix(V.par[Load.n.sym.a],D,N,byrow=FALSE)
  HH.risk.A<-t(Cc %*% t(Infectivity.A))
  Infectivity.B<-matrix(V.par[Load.n.sym.b],D,N,byrow=FALSE)
  HH.risk.B<-t(Cc %*% t(Infectivity.B))
  # Total within household risk, status identifies if an indiviual is present in the household
  in.house.A<-(status*hh.mat)*(Eta.A*HH.risk.A) 
  in.house.B<-(status*hh.mat)*(Eta.B*HH.risk.B) 
  # Community rate of exposure
  # --------------------------
  # Community transmission coefficient
  Eps.A<-exp(theta[['epsilon.A']]);Eps.B<-exp(theta[['epsilon.B']])
  # Effect of age on community level exposure
  x<-length(theta)-e.g;Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
  comm.age<-matrix(Eps.age[E],1,N)
  # Time varying component of community exposure/background density funtion
  comm.time.A<-matrix(Comm.risk.a,D,1);comm.time.B<-matrix(Comm.risk.b,D,1)
  # Combining the age effects and time varying function
  Comm.A <- Eps.A*(comm.time.A %*% comm.age)  
  Comm.B <- Eps.B*(comm.time.B %*% comm.age)
  # Total exposure rate 
  Lambda.A<-in.house.A + Comm.A
  Lambda.B<-in.house.B + Comm.B
  # _________________3. Probability of Infection________________________________
  Alpha.A<-(1-exp(-PI.A*Lambda.A));Alpha.B<-(1-exp(-PI.B*Lambda.B))
  # _________________4. Probability of Starting Shedding _______________________
  # Function that calculates the probability of onset given latency
  shedder.a<-function(t){
    days<-t-seq(0,length(Dist.lat)-1)
    return(Dist.lat[which(days>0)]  %*% Alpha.A[days[which(days>0)],])
  }
  PrShed.A<-t(sapply(seq(D),shedder.a))
  shedder.b<-function(t){
    days<-t-seq(0,length(Dist.lat)-1)
    return(Dist.lat[which(days>0)]  %*% Alpha.B[days[which(days>0)],])
  }
  PrShed.B<-t(sapply(seq(D),shedder.b))
  
  # Likelihood for days of observed data, ie onset and at risk day(ie present and
  # not shedding)
  nPrShed.A<-1-PrShed.A;nPrShed.B<-1-PrShed.B;pple<-seq(N)
  logger<-function(p){
    log.A<-sum(log(PrShed.A[(Onset.a[,p]==1),p])) + sum(log(nPrShed.A[(atRisk.a[,p]==1),p]))
    log.B<-sum(log(PrShed.B[(Onset.b[,p]==1),p])) + sum(log(nPrShed.B[(atRisk.b[,p]==1),p]))
    return(log.A + log.B)
  }
  Like<-sapply(pple, logger)
  return(sum(Like))
}
# --- Prior distribution
# Limits of the prior distribution for each parameters
prior.limits.group<-data.frame(lower=c(Prev.hom=-10,Prev.het=-10,age.2=-10,
                                 age.3=-10,age.4=-10,eta.A=-20,eta.B=-20,hh.size=-10,
                                 HighAsym=-10,LowSym=-10,HighSym=-10,epsilon.A=-20,
                                 epsilon.B=-20,eps.age2=-10,eps.age3=-10),
                         upper=c(Prev.A=10,Prev.B=10,age.2=10,age.3=10,age.4=10,eta.A=0, eta.b=0,
                                 hh.size=10,HighAsym=10,LowSym=10,HighSym=10,
                                 epsilon.A=0,epsilon.B=0,eps.age2=10,eps.age3=10))

IBM.prior.group<-function(theta) {
  sus <- function(val) return(dunif(val, min = -10, max = 10, log = TRUE))
  ## uniform prior on factor reducing susceptibility by Inf history
  inf.hist <-sapply(theta[1:2],sus)
  ## uniform prior on factor reducing susceptibility by age(4 grps)
  age.sus <- sapply(theta[3:(2+s.g-1)],sus)
  ## uniform prior on x1 such that eta=exp(x1) U[-10,0]
  Eta.a<-dunif(theta[['eta.A']], min = -20, max = 0, log = TRUE)
  Eta.b<-dunif(theta[['eta.B']], min = -20, max = 0, log = TRUE)
  ## uniform prior on household size covariate
  hh.size<-sus(theta[['hh.size']])
  # Infectivity mapping function parameters
  Inf.pars<-c(sus(theta[['HighAsym']]),sus(theta[['LowSym']]),sus(theta[['HighSym']]))
  ## uniform prior on x2 such that epsilon=exp(x2): U[-10,0]; for 3 categories
  Epsilon.a<-dunif(theta[['epsilon.A']], min = -20, max = 0, log = TRUE)
  Epsilon.b<-dunif(theta[['epsilon.B']], min = -20, max = 0, log = TRUE)
  ## unifrom prior on community eaxposure with age
  x<-length(theta)-e.g
  epsilon.age<-sapply(theta[(x+2):length(theta)],sus)
  ## the joint prior
  log.sum<-sum(inf.hist)+sum(age.sus)+Eta.a+Eta.b+hh.size+sum(Inf.pars)+Epsilon.a+Epsilon.b+sum(epsilon.age)
  ## default is to give the prior but otherwise give the log prior
  return(log.sum)
}
# --- Posterior distribution
IBM.posterior.group <- function(theta) {
  # calculate the log-prior for parameter `theta`
  x<-IBM.prior.group(theta)
  # calculate the log-likelihood for parameter `theta` and
  y<-IBM.LL.group(theta)
  # return the logged posterior probability
  return(x+y)
}













