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
Pos.cut.off<-35             # Cut-off used to define positive sample
source("DataPrep2.R")       # Modifies mother data set to shedding matrices
# based on cut-off given
source("CommRisk.R")        # Function creates time varying comm risk data,also 
# depends on the cut-off for positivity
source("Sampler_fitR.R")      # Loads some functions from the fitR package
# Model Functions
source("/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/Bayesian_IBM.All.funcs.R")
# Diagnostics functions
source("/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/Bayesian_IBM.All.diag.R")

#--------------------------------Data variables---------------------------------
shedding<-Shedding
onset<-Onset
status<-Ct.Data$status
sampledates<-as.Date(rownames(rsva.shedding))
# pdf of duration of latency
Dist.lat<-c(0, 0, 4, 4, 3, 1)/12 
# Community risk
Risk<-CommRisk()
z<-which(Risk$Hosp==0);Risk$Hosp[z]<-0.0001
s<-which(Risk$HH==0);Risk$HH[s]<-0.0000001
# Used in FOI function to add up household shedders
houses<-unique(Demo.data$hhid)
Cc<-matrix(0,N,N)
for (i in 1:length(houses)){
  # Their housemates
  pp<-which(Demo.data$hhid==houses[i])
  Cc[pp,pp]<-1
}
diag(Cc)<-0

#-------------------------------------------------------------------------------
# Stage 1
#-------------------------------------------------------------------------------
Proposal.sd1<-c(eta=40,epsilon=40)
N.iter<-10000
prior.limits1<-data.frame(lower=c(eta=-20,epsilon=-20),upper=c(eta=0,epsilon=0))

#---------------------------
# a) Constant community risk
#---------------------------
comm.risk<-rep(1,D)
# Chain 1
Paras1.C.c1<-log(c(eta=0.01,epsilon=0.01))
T1.C.c1<-system.time(Res1.C.c1<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c1, proposal.sd=Proposal.sd1,
                                          n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                          adapt.size.cooling=0.999))
save(Res1.C.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.C.c1.RData")
# Chain 2
Paras1.C.c2<-log(c(eta=0.01,epsilon=0.001))
T1.C.c2<-system.time(Res1.C.c2<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c2, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.C.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.C.c2.RData")
# Chain 3
Paras1.C.c3<-log(c(eta=0.001,epsilon=0.01))
T1.C.c3<-system.time(Res1.C.c3<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c3, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.C.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.C.c3.RData")
# Diagnostics(b=3000,th=0)
Res1.C.diag<-IBM.all.diag(Res1.C.c1,Res1.C.c2,Res1.C.c3,prior.limits1)
Mean1.C<-apply(Res1.C.diag$Res[names(Paras1.C.c1)],2,mean)
Median1.C<-apply(Res1.C.diag$Res[names(Paras1.C.c1)],2,median)
Max.post1.C<-Res1.C.diag$Res[ min(which(Res1.C.diag$Res[,3]==max(Res1.C.diag$Res[,3]))) ,-3]
exp(cbind(Mean1.C,Median1.C,as.numeric(Max.post1.C)))
Dic1.C.diag<-DIC.calc(Res1.C.diag$Res,Paras1.C.c1,IBM.prior.1)

#-------------------------------------
# b) Hospital incidence community risk
#-------------------------------------
comm.risk<-Risk$Hosp
# Chain 1
T1.H.c1<-system.time(Res1.H.c1<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c1, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.H.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.H.c1.RData")

# Chain 2
T1.H.c2<-system.time(Res1.H.c2<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c2, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.H.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.H.c2.RData")

# Chain 3
T1.H.c3<-system.time(Res1.H.c3<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c3, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.H.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.H.c3.RData")
# Diagnostics(3000 and 0)
Res1.H.diag<-IBM.all.diag(Res1.H.c1,Res1.H.c2,Res1.H.c3,prior.limits1)
Mean1.H<-apply(Res1.H.diag$Res[names(Paras1.C.c1)],2,mean)
Median1.H<-apply(Res1.H.diag$Res[names(Paras1.C.c1)],2,median)
Max.post1.H<-Res1.H.diag$Res[ min(which(Res1.H.diag$Res[,3]==max(Res1.H.diag$Res[,3]))) ,-3]
exp(cbind(Mean1.H,Median1.H,as.numeric(Max.post1.H)))
Dic1.H.diag<-DIC.calc(Res1.H.diag$Res,Paras1.C.c1,IBM.prior.1)

#------------------------------
# c) Gaussian HH community risk
#------------------------------
comm.risk<-Risk$HH
# Chain 1
T1.Hh.c1<-system.time(Res1.Hh.c1<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c1, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.Hh.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.Hh.c1.RData")

# Chain 2
T1.Hh.c2<-system.time(Res1.Hh.c2<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c2, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.Hh.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.Hh.c2.RData")

# Chain 3
T1.Hh.c3<-system.time(Res1.Hh.c3<-mcmcMH(target=IBM.posterior.1, init.theta=Paras1.C.c3, proposal.sd=Proposal.sd1,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res1.Hh.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res1.Hh.c3.RData")

# Diagnostics(4000 and 0)
Res1.Hh.diag<-IBM.all.diag(Res1.Hh.c1,Res1.Hh.c2,Res1.Hh.c3,prior.limits1)
Mean1.Hh<-apply(Res1.Hh.diag$Res[names(Paras1.C.c1)],2,mean)
Median1.Hh<-apply(Res1.Hh.diag$Res[names(Paras1.C.c1)],2,median)
Max.post<-Res1.Hh.diag$Res[ min(which(Res1.Hh.diag$Res[,3]==max(Res1.Hh.diag$Res[,3]))) ,-3]
exp(cbind(Mean1.Hh,Median1.Hh,as.numeric(Max.post)))
Dic1.Hh.diag<-DIC.calc(Res1.Hh.diag$Res,Paras1.C.c1,IBM.prior.1)

#-------------------------------------------------------------------------------
# Stage 2
#-------------------------------------------------------------------------------
# Changing covariate(Bit), infection history with 2 groups, coding is 1=no,2=yes
B <- matrix(1,D,N);B[postInf==1] <- 2 
Proposal.sd2<-c(Inf.hist=40,eta=40,epsilon=40)
N.iter<-10000
prior.limits2<-data.frame(lower=c(Inf.hist=-10,eta=-20,epsilon=-20),
                          upper=c(Inf.hist=10,eta=0,epsilon=0))

#------------------------------
# c) Gaussian HH community risk
#------------------------------
comm.risk<-Risk$HH
# Chain 1
Paras2.Hh.c1<-c(Inf.hist=log(1),Median1.Hh)
T2.Hh.c1<-system.time(Res2.Hh.c1<-mcmcMH(target=IBM.posterior.2, init.theta=Paras2.Hh.c1, proposal.sd=Proposal.sd2,
                                       n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                       adapt.size.cooling=0.999))
save(Res2.Hh.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res2.Hh.c1.RData")
# Chain 2
Paras2.Hh.c2<-c(Inf.hist=log(1.5),Median1.Hh)
T2.Hh.c2<-system.time(Res2.Hh.c2<-mcmcMH(target=IBM.posterior.2, init.theta=Paras2.Hh.c2, proposal.sd=Proposal.sd2,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res2.Hh.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res2.Hh.c2.RData")
# Chain 3
Paras2.Hh.c3<-c(Inf.hist=log(0.5),Median1.Hh)
T2.Hh.c3<-system.time(Res2.Hh.c3<-mcmcMH(target=IBM.posterior.2, init.theta=Paras2.Hh.c3, proposal.sd=Proposal.sd2,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res2.Hh.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res2.Hh.c3.RData")

# Diagnostics(3000 and 0)
Res2.Hh.diag<-IBM.all.diag(Res2.Hh.c1,Res2.Hh.c2,Res2.Hh.c3,prior.limits2)
Mean2.Hh<-apply(Res2.Hh.diag$Res[names(Paras2.Hh.c1)],2,mean)
Median2.Hh<-apply(Res2.Hh.diag$Res[names(Paras2.Hh.c1)],2,median)
x=(length(Paras2.Hh.c1)+1)
Max.post2.Hh<-Res2.Hh.diag$Res[ min(which(Res2.Hh.diag$Res[,x]==max(Res2.Hh.diag$Res[,x]))) ,-x]
exp(cbind(Mean2.Hh,Median2.Hh,as.numeric(Max.post2.Hh)))
Dic2.Hh.diag<-DIC.calc(Res2.Hh.diag$Res,Paras2.Hh.c1,IBM.prior.2)

#-------------------------------------------------------------------------------
# Stage 3
#-------------------------------------------------------------------------------
# Age covariate for susceptibility
A<-rep(4,N);s.g<-4
A[Demo.data$ageyrs<1]<-1;A[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
A[(Demo.data$ageyrs>=5) & (Demo.data$ageyrs<15)]<-3

Proposal.sd3<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,epsilon=40)
N.iter<-10000
prior.limits3<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,epsilon=-20),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,epsilon=0))
#------------------------------
# c) Gaussian HH community risk
#------------------------------
comm.risk<-Risk$HH
# Chain 1
Paras3.Hh.c1<-c(Median2.Hh[1],log(c(age.2=1,age.3=1,age.4=1)),Median2.Hh[2],Median2.Hh[3])

T3.Hh.c1<-system.time(Res3.Hh.c1<-mcmcMH(target=IBM.posterior.3, init.theta=Paras3.Hh.c1, proposal.sd=Proposal.sd3,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res3.Hh.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res3.Hh.c1.RData")
# Chain 2
Paras3.Hh.c2<-c(Median2.Hh[1],log(c(age.2=0.5,age.3=0.5,age.4=0.5)),Median2.Hh[2],Median2.Hh[3])

T3.Hh.c2<-system.time(Res3.Hh.c2<-mcmcMH(target=IBM.posterior.3, init.theta=Paras3.Hh.c2, proposal.sd=Proposal.sd3,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res3.Hh.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res3.Hh.c2.RData")
# Chain 3
Paras3.Hh.c3<-c(Median2.Hh[1],log(c(age.2=1.5,age.3=1.5,age.4=1.5)),Median2.Hh[2],Median2.Hh[3])

T3.Hh.c3<-system.time(Res3.Hh.c3<-mcmcMH(target=IBM.posterior.3, init.theta=Paras3.Hh.c3, proposal.sd=Proposal.sd3,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res3.Hh.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res3.Hh.c3.RData")
# Diagnostics(5000 and 0)
Res3.Hh.diag<-IBM.all.diag(Res3.Hh.c1,Res3.Hh.c2,Res3.Hh.c3,prior.limits3)
Mean3.Hh<-apply(Res3.Hh.diag$Res[names(Paras3.Hh.c1)],2,mean)
Median3.Hh<-apply(Res3.Hh.diag$Res[names(Paras3.Hh.c1)],2,median)
x=(length(Paras3.Hh.c1)+1)
Max.post3.Hh<-Res3.Hh.diag$Res[ min(which(Res3.Hh.diag$Res[,x]==max(Res3.Hh.diag$Res[,x]))) ,-x]
exp(cbind(Mean3.Hh,Median3.Hh,as.numeric(Max.post3.Hh)))
Dic3.Hh.diag<-DIC.calc(Res3.Hh.diag$Res,Paras3.Hh.c1,IBM.prior.3)

#-------------------------------------------------------------------------------
# Stage 4
#-------------------------------------------------------------------------------
# Age Covariate for community risk/exposure
E<-rep(3,N);e.g<-3
E[Demo.data$ageyrs<1]<-1;E[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2

Proposal.sd4<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,epsilon=40,eps.age2=40,
                eps.age3=40)
N.iter<-20000
prior.limits4<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,epsilon=-20,eps.age2=-10,
                                  eps.age3=-10),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,epsilon=0,eps.age2=10,
                                  eps.age3=10))
#------------------------------
# c) Gaussian HH community risk
#------------------------------
comm.risk<-Risk$HH
# Chain 1
Paras4.Hh.c1<-c(Median3.Hh,log(c(eps.age2=1,eps.age3=1)))
T4.Hh.c1<-system.time(Res4.Hh.c1<-mcmcMH(target=IBM.posterior.4, init.theta=Paras4.Hh.c1, proposal.sd=Proposal.sd4,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res4.Hh.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res4.Hh.c1.RData")
# Chain 2
Paras4.Hh.c2<-c(Median3.Hh,log(c(eps.age2=0.5,eps.age3=0.5)))
T4.Hh.c2<-system.time(Res4.Hh.c2<-mcmcMH(target=IBM.posterior.4, init.theta=Paras4.Hh.c2, proposal.sd=Proposal.sd4,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res4.Hh.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res4.Hh.c2.RData")
# Chain 3
Paras4.Hh.c3<-c(Median3.Hh,log(c(eps.age2=1.5,eps.age3=1.5)))
T4.Hh.c3<-system.time(Res4.Hh.c3<-mcmcMH(target=IBM.posterior.4, init.theta=Paras4.Hh.c2, proposal.sd=Proposal.sd4,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res4.Hh.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res4.Hh.c3.RData")

# Diagnostics(1000 and 0)
Res4.Hh.diag<-IBM.all.diag(Res4.Hh.c1,Res4.Hh.c2,Res4.Hh.c3,prior.limits4)
Mean4.Hh<-apply(Res4.Hh.diag$Res[names(Paras4.Hh.c1)],2,mean)
Median4.Hh<-apply(Res4.Hh.diag$Res[names(Paras4.Hh.c1)],2,median)
x=(length(Paras4.Hh.c1)+1)
Max.post4.Hh<-Res4.Hh.diag$Res[ min(which(Res4.Hh.diag$Res[,x]==max(Res4.Hh.diag$Res[,x]))) ,-x]
exp(cbind(Mean4.Hh,Median4.Hh,as.numeric(Max.post4.Hh)))
Dic4.Hh.diag<-DIC.calc(Res4.Hh.diag$Res,Paras4.Hh.c1,IBM.prior.4)

#-------------------------------------------------------------------------------
# Stage 5
#-------------------------------------------------------------------------------
# a) Contact matrices by age(no convergence)
#----------------------------
# Age defined contact types covariate (groups:[0-5),>=5)
C.g<-rep(2,N);C.g[Demo.data$ageyrs<5]<-1;
Contacts1<-matrix(0,N,N)
for (i in 1:N){
  # Their housemates
  pp<-which(Demo.data$hhid==Demo.data$hhid[i])
  # HH members age group
  if(C.g[i]==1){
    Contacts1[i,pp]<-C.g[pp]
  } else{
    Contacts1[i,pp]<-C.g[pp]+2
  } 
}
C<-Contacts1+1

Proposal.sd5<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,
                C12=40,C21=40,C22=40,epsilon=40,eps.age2=40,eps.age3=40)
N.iter<-20000
prior.limits5<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,
                                  C12=-10,C21=-10,C22=-10,epsilon=-20,eps.age2=-10,
                                  eps.age3=-10),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,
                                  C12=10,C21=10,C22=10,epsilon=0,eps.age2=10,
                                  eps.age3=10))

# Chain 1
Paras5.Hh.c1<-c(Median4.Hh[1:5],log(c(C12=1,C21=1,C22=1)),Median4.Hh[6:8])
T5.Hh.c1<-system.time(Res5.Hh.c1<-mcmcMH(target=IBM.posterior.5, init.theta=Paras5.Hh.c1, proposal.sd=Proposal.sd5,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res5.Hh.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res5.Hh.c1.RData")
# Chain 2
Paras5.Hh.c2<-c(Median4.Hh[1:5],log(c(C12=0.5,C21=0.5,C22=0.5)),Median4.Hh[6:8])
T5.Hh.c2<-system.time(Res5.Hh.c2<-mcmcMH(target=IBM.posterior.5, init.theta=Paras5.Hh.c2, proposal.sd=Proposal.sd5,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res5.Hh.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res5.Hh.c2.RData")
# Chain 3
Paras5.Hh.c3<-c(Median4.Hh[1:5],log(c(C12=1.5,C21=1.5,C22=1.5)),Median4.Hh[6:8])
T5.Hh.c3<-system.time(Res5.Hh.c3<-mcmcMH(target=IBM.posterior.5, init.theta=Paras5.Hh.c2, proposal.sd=Proposal.sd5,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res5.Hh.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res5.Hh.c3.RData")

# Diagnostics(7500 and 3)
Res5.Hh.diag<-IBM.all.diag(Res5.Hh.c1,Res5.Hh.c2,Res5.Hh.c3,prior.limits5)

# trying different age groupings;<1,>=1
C.g<-rep(2,N);C.g[Demo.data$ageyrs<1]<-1;
Contacts1<-matrix(0,N,N)
for (i in 1:N){
  # Their housemates
  pp<-which(Demo.data$hhid==Demo.data$hhid[i])
  # HH members age group
  if(C.g[i]==1){
    Contacts1[i,pp]<-C.g[pp]
  } else{
    Contacts1[i,pp]<-C.g[pp]+2
  } 
}
C<-Contacts1+1
# Chain 1
T5.Hh.c1.2<-system.time(Res5.Hh.c1.2<-mcmcMH(target=IBM.posterior.5, init.theta=Paras5.Hh.c1, proposal.sd=Proposal.sd5,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res5.Hh.c1.2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res5.Hh.c1.2.RData")

Res5.Hh.diag<-IBM.all.diag(Res5.Hh.c1.2,Res5.Hh.c1.2,Res5.Hh.c1.2,prior.limits5)

#-------------------------------------------------------------------------------
# Stage 6
#-------------------------------------------------------------------------------
# a) DD vs FD (poor converegence with omega ->0)
#----------------------------
Nh<-matrix(Demo.data$hhsize,D,N, byrow = TRUE)
Proposal.sd6<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,
                omega=40,epsilon=40,eps.age2=40,eps.age3=40)
N.iter<-40000
prior.limits6<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,
                                  omega=-10,epsilon=-20,eps.age2=-10,eps.age3=-10),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,
                                  omega=10,epsilon=0,eps.age2=10,eps.age3=10))

# Chain 1
Paras6.Hh.c1<-c(Median4.Hh[1:5],log(c(omega=4e-4)),Median4.Hh[6:8])
T6.Hh.c1<-system.time(Res6.Hh.c1<-mcmcMH(target=IBM.posterior.6, init.theta=Paras6.Hh.c1, proposal.sd=Proposal.sd6,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res6.Hh.c1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c1.RData")
# Chain 2
Paras6.Hh.c2<-c(Median4.Hh[1:5],log(c(omega=1)),Median4.Hh[6:8])
T6.Hh.c2<-system.time(Res6.Hh.c2<-mcmcMH(target=IBM.posterior.6, init.theta=Paras6.Hh.c2, proposal.sd=Proposal.sd6,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res6.Hh.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c2.RData")
# Chain 3
Paras6.Hh.c3<-c(Median4.Hh[1:5],log(c(omega=2)),Median4.Hh[6:8])
T6.Hh.c3<-system.time(Res6.Hh.c3<-mcmcMH(target=IBM.posterior.6, init.theta=Paras6.Hh.c3, proposal.sd=Proposal.sd6,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res6.Hh.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c3.RData")
# Chain 4
Paras6.Hh.c4<-c(Median4.Hh[1:5],log(c(omega=1000)),Median4.Hh[6:8])
T6.Hh.c4<-system.time(Res6.Hh.c4<-mcmcMH(target=IBM.posterior.6, init.theta=Paras6.Hh.c4, proposal.sd=Proposal.sd6,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res6.Hh.c4,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c4.RData")

# Diagnostics()
Res6.Hh.diag<-IBM.all.diag(Res6.Hh.c1,Res6.Hh.c2,Res6.Hh.c3,prior.limits6)
# Res6.Hh.diag<-IBM.all.diag(Res6.Hh.c2,Res6.Hh.c3,Res6.Hh.c4,prior.limits6)
plot(density(exp(Res6.Hh.diag$Res[,6])),type='l')
Mean4.Hh<-apply(Res4.Hh.diag$Res[names(Paras4.Hh.c1)],2,mean)
Median4.Hh<-apply(Res4.Hh.diag$Res[names(Paras4.Hh.c1)],2,median)
x=(length(Paras4.Hh.c1)+1)
Max.post4.Hh<-Res4.Hh.diag$Res[ min(which(Res4.Hh.diag$Res[,x]==max(Res4.Hh.diag$Res[,x]))) ,-x]
exp(cbind(Mean4.Hh,Median4.Hh,as.numeric(Max.post4.Hh)))
Dic4.Hh.diag<-DIC.calc(Res4.Hh.diag$Res,Paras4.Hh.c1,IBM.prior.4)

#----------------------------
# b) HH size categorical (converged)
#----------------------------
# HH size covariate
H<-rep(2,N);H[Demo.data$hhsize<8]<-1

Proposal.sd6.b<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,
                  hh.size=40,epsilon=40,eps.age2=40,eps.age3=40)
N.iter<-20000
prior.limits6.b<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,
                                    hh.size=-10,epsilon=-20,eps.age2=-10,eps.age3=-10),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,
                                  hh.size=10,epsilon=0,eps.age2=10,eps.age3=10))

# Chain 1
Paras6.Hh.c1.b<-c(Median4.Hh[1:5],log(c(hh.size=1)),Median4.Hh[6:8])
T6.Hh.c1.b<-system.time(Res6.Hh.c1.b<-mcmcMH(target=IBM.posterior.6.b, init.theta=Paras6.Hh.c1.b, proposal.sd=Proposal.sd6.b,
                                         n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                         adapt.size.cooling=0.999))
save(Res6.Hh.c1.b,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c1.b.RData")

# Chain 2
Paras6.Hh.c2.b<-c(Median4.Hh[1:5],log(c(hh.size=0.5)),Median4.Hh[6:8])
T6.Hh.c2.b<-system.time(Res6.Hh.c2.b<-mcmcMH(target=IBM.posterior.6.b, init.theta=Paras6.Hh.c2.b, proposal.sd=Proposal.sd6.b,
                                             n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                             adapt.size.cooling=0.999))
save(Res6.Hh.c2.b,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c2.b.RData")
# Chain 3
Paras6.Hh.c3.b<-c(Median4.Hh[1:5],log(c(hh.size=1.5)),Median4.Hh[6:8])
T6.Hh.c3.b<-system.time(Res6.Hh.c3.b<-mcmcMH(target=IBM.posterior.6.b, init.theta=Paras6.Hh.c3.b, proposal.sd=Proposal.sd6.b,
                                             n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                             adapt.size.cooling=0.999))
save(Res6.Hh.c3.b,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/IBM_All_stages/BayIBM_Res6.Hh.c3.b.RData")
# Diagnostics()
Res6.Hh.diag.b<-IBM.all.diag(Res6.Hh.c1.b,Res6.Hh.c2.b,Res6.Hh.c3.b,prior.limits6.b)
Mean6.Hh.b<-apply(Res6.Hh.diag.b$Res[names(Paras6.Hh.c1.b)],2,mean)
Median6.Hh.b<-apply(Res6.Hh.diag.b$Res[names(Paras6.Hh.c1.b)],2,median)
x=(length(Paras6.Hh.c1.b)+1)
Max.post6.Hh.b<-Res6.Hh.diag.b$Res[ min(which(Res6.Hh.diag.b$Res[,x]==max(Res6.Hh.diag.b$Res[,x]))) ,-x]
exp(cbind(Mean6.Hh.b,Median6.Hh.b,as.numeric(Max.post6.Hh.b)))
Dic6.Hh.diag.b<-DIC.calc(Res6.Hh.diag.b$Res,Paras6.Hh.c1.b,IBM.prior.6.b)

#-------------------------------------------------------------------------------
# Stage 7
#-------------------------------------------------------------------------------
# a) continuous viral load (no convergence)
#----------------------------
# one
V.map<-function(v,para){
  beta=para[1];alpha=para[2];b=para[3]
  v[v!=0]<-b/(1+exp(-beta*(v[v!=0]-alpha)))
  return(v)
} 
shedding<-Shedding.est.load

Proposal.sd1<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,
                hh.size=40,beta=40,alpha=40,b=40,epsilon=40,eps.age2=40,
                eps.age3=40)
prior.limits7.a<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,
                                  hh.size=-10,beta=-5,alpha=-10,b=-10,epsilon=-20,eps.age2=-10,
                                  eps.age3=-10),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,
                                  hh.size=10,beta=5,alpha=5,b=5,epsilon=0,eps.age2=10,
                                  eps.age3=10))

# Chain 1
Paras1<-log(c(Inf.hist=0.601,age.2=1.092,age.3=0.370,age.4=0.196,eta=0.0237,
              hh.size=0.454,beta=100,alpha=V.load(Pos.cut.off),b=1,epsilon=0.00944,eps.age2=0.467,
              eps.age3=1.442))
T.6.1<-system.time(Res.6.a.1<-mcmcMH(target=IBM.posterior.7.a, init.theta=Paras1, proposal.sd=Proposal.sd1,
                                     n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                     adapt.size.cooling=0.999))
save(Res.6.a.1,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/BayIBM_Res.6.a.1.RData")
# Chain 2
Paras1.c2<-log(c(Inf.hist=0.601,age.2=1.092,age.3=0.370,age.4=0.196,eta=0.0237,
                 hh.size=0.454,beta=1,alpha=7,b=4,epsilon=0.00944,eps.age2=0.467,
                 eps.age3=1.442))
T.6.1.c2<-system.time(Res.6.a.1.c2<-mcmcMH(target=IBM.posterior.7.a, init.theta=Paras1.c2, proposal.sd=Proposal.sd1,
                                           n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                           adapt.size.cooling=0.999))
save(Res.6.a.1.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/BayIBM_Res.6.a.1.c2.RData")

# Chain 3 
Paras1.c3<-log(c(Inf.hist=0.601,age.2=1.092,age.3=0.370,age.4=0.196,eta=0.0237,
                 hh.size=0.454,beta=0.4,alpha=1e-4,b=4,epsilon=0.00944,eps.age2=0.467,
                 eps.age3=1.442))
T.6.1.c3<-system.time(Res.6.a.1.c3<-mcmcMH(target=IBM.posterior.7.a, init.theta=Paras1.c3, proposal.sd=Proposal.sd1,
                                           n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                           adapt.size.cooling=0.999))
save(Res.6.a.1.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/BayIBM_Res.6.a.1.c3.RData")

# Diagnostics()
Res7.Hh.diag.a<-IBM.all.diag(Res.6.a.1,Res.6.a.1.c2,Res.6.a.1.c3,prior.limits7.a)


# b) Categorical viral load 
#----------------------------
# Categorical variable for shedding
high.load<-6
Shedding.cat<-matrix(1,D,N)
Shedding.cat[which(shedding>0,TRUE)]<-2
Shedding.cat[which(shedding>high.load,TRUE)]<-3
N.iter<-500000
Proposal.sd1<-c(Inf.hist=40,age.2=40,age.3=40,age.4=40,eta=40,
                hh.size=40,Infec=40,epsilon=40,eps.age2=40,
                eps.age3=40)
prior.limits7.b<-data.frame(lower=c(Inf.hist=-10,age.2=-10,age.3=-10,age.4=-10,eta=-20,
                                  hh.size=-10,Infec=-10,epsilon=-20,eps.age2=-10,
                                  eps.age3=-10),
                          upper=c(Inf.hist=10,age.2=10,age.3=10,age.4=10,eta=0,
                                  hh.size=10,Infec=10,epsilon=0,eps.age2=10,
                                  eps.age3=10))
# Chain 1
Paras7.Hh.c1.b<-log(c(Inf.hist=0.601,age.2=1.092,age.3=0.370,age.4=0.196,eta=0.0237,
              hh.size=0.454,Infec=0.5,epsilon=0.00944,eps.age2=0.467,
              eps.age3=1.442))

T6.b.six<-system.time(Res.6.b.six<-mcmcMH(target=IBM.posterior.7.b, init.theta=Paras7.Hh.c1.b, proposal.sd=Proposal.sd1,
                                          n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                          adapt.size.cooling=0.999))
save(Res.6.b.six,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/BayIBM_Res.6.b.six 2.RData")

# Chain 2
Paras7.Hh.c2.b<-log(c(Inf.hist=0.601,age.2=1.092,age.3=0.370,age.4=0.196,eta=0.0237,
                 hh.size=0.454,Infec=1,epsilon=0.00944,eps.age2=0.467,
                 eps.age3=1.442))
T6.b.six.c2<-system.time(Res.6.b.six.c2<-mcmcMH(target=IBM.posterior.7.b, init.theta=Paras7.Hh.c2.b, proposal.sd=Proposal.sd1,
                                                n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                                adapt.size.cooling=0.999))
save(Res.6.b.six.c2,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/BayIBM_Res.6.b.six.c2 2.RData")

# Chain 3
Paras7.Hh.c3.b<-log(c(Inf.hist=0.601,age.2=1.092,age.3=0.370,age.4=0.196,eta=0.0237,
                 hh.size=0.454,Infec=1.5,epsilon=0.00944,eps.age2=0.467,
                 eps.age3=1.442))
T6.b.six.c3<-system.time(Res.6.b.six.c3<-mcmcMH(target=IBM.posterior.7.b, init.theta=Paras7.Hh.c3.b, proposal.sd=Proposal.sd1,
                                                n.iterations=N.iter,adapt.size.start = 1000,adapt.shape.start = 500,
                                                adapt.size.cooling=0.999))
save(Res.6.b.six.c3,file="/Users/ikombe/Desktop/IKombe/PhD/Codes/HH_model/BayIBM_Res.6.b.six.c3 2.RData")

# Diagnostics(30000 and 0)
Res7.Hh.diag.b<-IBM.all.diag(Res.6.b.six,Res.6.b.six.c2,Res.6.b.six.c3,prior.limits7.b)
Mean7.Hh.b<-apply(Res7.Hh.diag.b$Res[names(Paras7.Hh.c1.b)],2,mean)
Median7.Hh.b<-apply(Res7.Hh.diag.b$Res[names(Paras7.Hh.c1.b)],2,median)
x=(length(Paras7.Hh.c1.b)+1)
Max.post7.Hh.b<-Res7.Hh.diag.b$Res[ min(which(Res7.Hh.diag.b$Res[,x]==max(Res7.Hh.diag.b$Res[,x]))) ,-x]
exp(cbind(Mean7.Hh.b,Median7.Hh.b,as.numeric(Max.post7.Hh.b)))
Dic7.Hh.diag.b<-DIC.calc(Res7.Hh.diag.b$Res,Paras7.Hh.c1.b,IBM.prior.7.b)



