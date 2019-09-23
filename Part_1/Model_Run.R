# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 13th September 2018
################################################################################
# This is the main model fitting function. It calls the data, the model functions, 
# the model fitting functions and the diagnostic fuctions. The diagnosis part
# is interactive, trace plots are plotted and the user has to input the burn-in 
# and the thinning interval, there are no default values. The output of the 
# diagnosis is a caterpiller plot showing the distributions of the fitted parameters
# and a plot showing correlation patterns.
# 
# Parameter fitting is done in two parts:
# Part 1: Data distinguishes between RSV A and B and fits an interactive model
#         Given the posterior estimates of the parameters, distributions of the
#         per person household transmission rate and community transmission rate
#         are plotted. This is Figure 3. in the main article.
# Part 2: The data is fitted as a single pathogen, either only RSV A, only RSV B 
#         or a combination where the data is not separated into RSV A or B, just 
#         RSV.
# 
################################################################################

# --- Loading the data needed for fitting; shedding information and some demographics
source("Data_Prep.R")
# --- Loading the model functions: the likelihood, prior and posterior
source("Model_Funcs.R")
# --- Loading the fitting functions (from the fitR package) 
source("Sampler_fitR.R") 
# --- Loading the diagnostic functions; Takes in three MCMC chains
source("Model_Diag.R")
# --- Defining some variables needed for fiting that are not pathogen specific
# -pdf of duration of latency
Dist.lat<-c(0, 0, 4, 4, 3, 1)/12 
# -Age covariate for susceptibility
A<-rep(4,N);s.g<-4
A[Demo.data$ageyrs<1]<-1;A[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
A[(Demo.data$ageyrs>=5) & (Demo.data$ageyrs<15)]<-3
# -Age Covariate for community risk/exposure
E<-rep(3,N);e.g<-3
E[Demo.data$ageyrs<1]<-1;E[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
# -HH size covariate
H<-rep(2,N);H[Demo.data$hhsize<8]<-1
# -An NxN matrix the identifies housemates.Used in rate of exposure function 
# within the likelihood function to add up HH shedders
houses<-unique(Demo.data$hhid)
Cc<-matrix(0,N,N)
for (i in 1:length(houses)){
  # Their housemates
  pp<-which(Demo.data$hhid==houses[i])
  Cc[pp,pp]<-1
}
diag(Cc)<-0
################################### PART 1 #####################################

# --- Algorithm specifics
# The standard deviation for the multivariate proposal distribution
Proposal.sd1<-c(Prev.hom=40,Prev.het=40,age.2=40,age.3=40,age.4=40,eta.A=40,eta.B=40,
                hh.size=40,HighAsym=40,LowSym=40,HighSym=40,epsilon.A=40, epsilon.B=40,
                eps.age2=40,eps.age3=40)
# The number of iterations to run
N.iter<-250000
# Initial parameter values for the three MCMC chains to be run 
Paras1<-log(c(Prev.hom=0.658,Prev.het=0.658,age.2=1.113,age.3=0.364,age.4=0.195,
              eta.A=0.0175,eta.B=0.0175,hh.size=0.423,HighAsym=3,LowSym=0.5,HighSym=1.5,
              epsilon.A=0.00952,epsilon.B=0.00952,eps.age2=0.446,eps.age3=1.473))
Paras1.c2<-log(c(Prev.hom=0.158,Prev.het=0.658,age.2=1.113,age.3=0.364,age.4=0.195,
                 eta.A=0.035,eta.B=0.0175,hh.size=0.423,HighAsym=1,LowSym=1,HighSym=1,
                 epsilon.A=0.0476,epsilon.B=0.00952,eps.age2=0.446,eps.age3=1.473))
Paras1.c3<-log(c(Prev.hom=1,Prev.het=1,age.2=1,age.3=1,age.4=1,
                 eta.A=0.02,eta.B=0.02,hh.size=1,HighAsym=5,LowSym=2,HighSym=2,
                 epsilon.A=0.002,epsilon.B=0.002,eps.age2=1,eps.age3=1))

# --- Defining some variables needed for fiting that are pathogen specific
# -Extracting matrices from the list object output of 'Data_Prep.R'
# RSV A
Shedding.a<<-Shed.data$Shedding; Onset.a<<-Shed.data$Onset; Shedding.a.est.ct<<-Shed.data$Shedding.est.ct;
Shedding.a.est.load<<-Shed.data$Shedding.est.load;status<<-Shed.data$status;
atRisk.a<<-Shed.data$atRisk; postInf.a<<-Shed.data$postInf; Comm.risk.a<<-Shed.data$Comm.risk;
# RSV B
Shedding.b<<-Shed.data.2$Shedding; Onset.b<<-Shed.data.2$Onset; Shedding.b.est.ct<<-Shed.data.2$Shedding.est.ct;
Shedding.b.est.load<<-Shed.data.2$Shedding.est.load;
atRisk.b<<-Shed.data.2$atRisk; postInf.b<<-Shed.data.2$postInf; Comm.risk.b<<-Shed.data.2$Comm.risk;

# -Changing covariate for infection history with 3 groups
B <<- matrix(1,D,N);B[postInf.a==1]<-2;B[postInf.b==1]<-3
B[(postInf.a==1 & postInf.b==1)]<-4 

# -Cut-off for high viral load
high.load<<-6
# -Combined categories of viral load and symptoms
Load.n.sym.a<-matrix(1,D,N);Load.n.sym.b<-matrix(1,D,N)
# Reference group;low viral load and asymptomatic
Load.n.sym.a[which(Shedding.a.est.load>0 & ARI.estimated!=1,TRUE)]<-2
Load.n.sym.b[which(Shedding.b.est.load>0 & ARI.estimated.2!=1,TRUE)]<-2
# high viral load and asymptomatic
Load.n.sym.a[which(Shedding.a.est.load>high.load & ARI.estimated!=1,TRUE)]<-3
Load.n.sym.b[which(Shedding.b.est.load>high.load & ARI.estimated.2!=1,TRUE)]<-3
# low viral load and symptomatic
Load.n.sym.a[which(Shedding.a.est.load>0 & ARI.estimated==1,TRUE)]<-4
Load.n.sym.b[which(Shedding.b.est.load>0 & ARI.estimated.2==1,TRUE)]<-4
# high viral load and symptomatic
Load.n.sym.a[which(Shedding.a.est.load>high.load & ARI.estimated==1,TRUE)]<-5
Load.n.sym.b[which(Shedding.b.est.load>high.load & ARI.estimated.2==1,TRUE)]<-5

# --- Running 3 chains intiated at different points
# Chain 1
Res.group.c1<-mcmcMH(target=IBM.posterior.group,init.theta=Paras1, 
                         proposal.sd=Proposal.sd1,n.iterations=N.iter,
                         adapt.size.start = 1000,adapt.shape.start = 500,
                         adapt.size.cooling=0.999)
# Chain 2
Res.group.c2<-mcmcMH(target=IBM.posterior.group,init.theta=Paras1.c2, 
                         proposal.sd=Proposal.sd1,n.iterations=N.iter,
                         adapt.size.start = 1000,adapt.shape.start = 500,
                         adapt.size.cooling=0.999)
# Chain 3
Res.group.c3<-mcmcMH(target=IBM.posterior.group,init.theta=Paras1.c3, 
                         proposal.sd=Proposal.sd1,n.iterations=N.iter,
                         adapt.size.start = 1000,adapt.shape.start = 500,
                         adapt.size.cooling=0.999)

# --- Diagnostics
# Running the main diagnostic function in 'Model_Diag.R'
Res.group.diag<-IBM.all.diag(chain1=Res.group.c1,chain2=Res.group.c2,chain3=Res.group.c3,
                               prior.limits=prior.limits.group,desc='RSV A and B')
# Save the chain results so that they can be accessed by the validation code ***
# In this case we created a folder 'Para_distributions' that contains the results 
# of each chain run.
save(Res.group.diag,file='Para_distributions/Res.group.diag.RData')

# Calculating the median and 95% Credible intervals
Res.summary.RSV<-t(exp(apply(Res.group.diag$Res[names(Paras1)],2,Quntile.func)))

# --- Plotting the distributions of the per person household transmission rate 
#     and community transmission rate. This is Figure 3. in the main 
#     article.
Rates.plot(Par.dist=Res.group.diag$Res[,-16])
################################### PART 2 #####################################
# --- Algorithm specifics
# The standard deviation for the multivariate proposal distribution
Proposal.sd1<<-c(Prev.inf=40,age.2=40,age.3=40,age.4=40,eta=40,hh.size=40,
                 HighAsym=40,LowSym=40,HighSym=40,epsilon=40,eps.age2=40,eps.age3=40)
# The number of iterations to run
N.iter<-100000
# Initial parameter values for the three MCMC chains to be run 
Paras1<<-log(c(Prev.inf=0.5,age.2=0.8,age.3=0.4,age.4=0.2,eta=0.02,hh.size=0.5,
               HighAsym=3,LowSym=0.5,HighSym=1.5,epsilon=0.002,eps.age2=1,eps.age3=1.5))
Paras1.c2<-log(c(Prev.inf=1,age.2=1,age.3=1,age.4=1,eta=0.02,hh.size=1,HighAsym=1,
                 LowSym=1,HighSym=1,epsilon=0.02,eps.age2=1,eps.age3=1))
Paras1.c3<-log(c(Prev.inf=1.5,age.2=1.5,age.3=1.5,age.4=1.5,eta=0.01,hh.size=1.5,
                 HighAsym=1.5,LowSym=1.5,HighSym=1.5,epsilon=0.02,eps.age2=1.5,
                 eps.age3=1.5))
# --- RSV A fitting 
# --- Defining some variables needed for fiting that are pathogen specific
Shedding<<-Shed.data$Shedding; Onset<<-Shed.data$Onset; Shedding.est.ct<<-Shed.data$Shedding.est.ct;
Shedding.est.load<<-Shed.data$Shedding.est.load;status<<-Shed.data$status;
atRisk<<-Shed.data$atRisk; postInf<<-Shed.data$postInf; Comm.risk<<-Shed.data$Comm.risk;
# -Changing covariate tracking infection history
B <<- matrix(1,D,N);B[postInf==1]<-2
# -Cut-off for high viral load
high.load<<-6
# -Combined categories of viral load and symptoms
Load.n.sym<-matrix(1,D,N)
# Reference group;low viral load and asymptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load<=high.load) & ARI.estimated==0]<-2                      
# high viral load and asymptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load>high.load) & ARI.estimated==0]<-3                      
# low viral load and symptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load<=high.load) & ARI.estimated==1]<-4                      
# high viral load and symptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load>high.load) & ARI.estimated==1]<-5           

# --- Running 3 chains intiated at different points
# Chain 1
Res.RSV.a.c1<-mcmcMH(target=IBM.posterior,init.theta=Paras1, 
                     proposal.sd=Proposal.sd1,n.iterations=N.iter,
                     adapt.size.start = 1000,adapt.shape.start = 500,
                     adapt.size.cooling=0.999)
# Chain 2
Res.RSV.a.c2<-mcmcMH(target=IBM.posterior,init.theta=Paras1.c2, 
                     proposal.sd=Proposal.sd1,n.iterations=N.iter,
                     adapt.size.start = 1000,adapt.shape.start = 500,
                     adapt.size.cooling=0.999)
# Chain 3
Res.RSV.a.c3<-mcmcMH(target=IBM.posterior,init.theta=Paras1.c3, 
                     proposal.sd=Proposal.sd1,n.iterations=N.iter,
                     adapt.size.start = 1000,adapt.shape.start = 500,
                     adapt.size.cooling=0.999)

# --- RSV B fitting 
# --- Defining some variables needed for fiting that are pathogen specific
Shedding<<-Shed.data.2$Shedding; Onset<<-Shed.data.2$Onset; Shedding.est.ct<<-Shed.data.2$Shedding.est.ct;
Shedding.est.load<<-Shed.data.2$Shedding.est.load;status<<-Shed.data.2$status;
atRisk<<-Shed.data.2$atRisk; postInf<<-Shed.data.2$postInf; Comm.risk<<-Shed.data.2$Comm.risk;
# -Changing covariate tracking infection history
B <<- matrix(1,D,N);B[postInf==1]<-2
# -Cut-off for high viral load
high.load<<-6
# -Combined categories of viral load and symptoms
Load.n.sym<-matrix(1,D,N)
# Reference group;low viral load and asymptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load<=high.load) & ARI.estimated==0]<-2                      
# high viral load and asymptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load>high.load) & ARI.estimated==0]<-3                      
# low viral load and symptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load<=high.load) & ARI.estimated==1]<-4                      
# high viral load and symptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load>high.load) & ARI.estimated==1]<-5           

# --- Running 3 chains intiated at different points
# Chain 1
Res.RSV.b.c1<-mcmcMH(target=IBM.posterior,init.theta=Paras1, 
                     proposal.sd=Proposal.sd1,n.iterations=N.iter,
                     adapt.size.start = 1000,adapt.shape.start = 500,
                     adapt.size.cooling=0.999)
# Chain 2
Res.RSV.b.c2<-mcmcMH(target=IBM.posterior,init.theta=Paras1.c2, 
                     proposal.sd=Proposal.sd1,n.iterations=N.iter,
                     adapt.size.start = 1000,adapt.shape.start = 500,
                     adapt.size.cooling=0.999)
# Chain 3
Res.RSV.b.c3<-mcmcMH(target=IBM.posterior,init.theta=Paras1.c3, 
                     proposal.sd=Proposal.sd1,n.iterations=N.iter,
                     adapt.size.start = 1000,adapt.shape.start = 500,
                     adapt.size.cooling=0.999)

#----RSV fitting 
Shedding<<-Shed.data.3$Shedding; Onset<<-Shed.data.3$Onset; 
Shedding.est.load<<-Shed.data.3$Shedding.est.load;status<<-Shed.data.3$status;
atRisk<<-Shed.data.3$atRisk; B<<-Shed.data.3$B;B[B==0]<-1; Comm.risk<<-Shed.data.3$Comm.risk;
# -Cut-off for high viral load
high.load<<-6
# -Combined categories of viral load and symptoms
Load.n.sym<-matrix(1,D,N)
# Reference group;low viral load and asymptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load<=high.load) & ARI.estimated.3==0]<-2                      
# high viral load and asymptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load>high.load) & ARI.estimated.3==0]<-3                      
# low viral load and symptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load<=high.load) & ARI.estimated.3==1]<-4                      
# high viral load and symptomatic
Load.n.sym[(Shedding.est.load>0 & Shedding.est.load>high.load) & ARI.estimated.3==1]<-5           
# --- Running 3 chains intiated at different points
# Chain 1
Res.RSV.c1<-mcmcMH(target=IBM.posterior,init.theta=Paras1, 
                   proposal.sd=Proposal.sd1,n.iterations=N.iter,
                   adapt.size.start = 1000,adapt.shape.start = 500,
                   adapt.size.cooling=0.999)
# Chain 2
Res.RSV.c2<-mcmcMH(target=IBM.posterior,init.theta=Paras1.c2, 
                   proposal.sd=Proposal.sd1,n.iterations=N.iter,
                   adapt.size.start = 1000,adapt.shape.start = 500,
                   adapt.size.cooling=0.999)
# Chain 3
Res.RSV.c3<-mcmcMH(target=IBM.posterior,init.theta=Paras1.c3, 
                   proposal.sd=Proposal.sd1,n.iterations=N.iter,
                   adapt.size.start = 1000,adapt.shape.start = 500,
                   adapt.size.cooling=0.999)

# ---- Diagnostics
# RSV A
Res.diag.a<-IBM.all.diag(chain1=Res.RSV.a.c1,chain2=Res.RSV.a.c2,
                         chain3=Res.RSV.a.c3,prior.limits=prior.limits,desc='RSV A')
# Calculating the median and 95% Credible intervals
Res.summary.a<-t(exp(apply(Res.diag.a$Res[names(Paras1)],2,Quntile.func)))

# RSV B
Res.diag.b<-IBM.all.diag(chain1=Res.RSV.b.c1,chain2=Res.RSV.b.c2,
                         chain3=Res.RSV.b.c3,prior.limits=prior.limits,desc='RSV B')
# Calculating the median and 95% Credible intervals
Res.summary.b<-t(exp(apply(Res.diag.b$Res[names(Paras1)],2,Quntile.func)))

# RSV 
Res.diag.rsv<-IBM.all.diag(chain1=Res.RSV.c1,chain2=Res.RSV.c2,
                           chain3=Res.RSV.c3,prior.limits=prior.limits,desc='RSV')
# Calculating the median and 95% Credible intervals
Res.summary.rsv<-t(exp(apply(Res.diag.rsv$Res[names(Paras1)],2,Quntile.func)))





















