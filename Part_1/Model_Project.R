# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 13th September 2018
################################################################################
# This code tries to quantify the contributions of asymptomatic and symotomatic 
# individuals to transmission. This will be done in two parts as described:
# Part 1: Contribution of symptomatic individuals to transmission
#         Assume that symptomatic individuals are as infectious as asymptomatic
#         ones and simulate. The difference in the total numbers infected between
#         the data and the simulations should give the contribution of symptomatics
# Part 2: Contribution of asymptomatic individuals.
#         -If asymptomatic individuals were not important to transmission what
#          would the epidemic look like? Assume asymptomatics are not infectious
#          and simulate
#         -If asymptomatic individuals were not sampled, what inference would we
#         make on the parameters? Remove all asymptomatics are re-estimate the 
#         relevant parameters.
# At this point, all asymptomatic individuals ars treated as one group as opposed
# to high viral load asymptomatics and low viral load asymptomatics. This is due
# to the fact that the infered posterior distribution was very wide and hence
# not very informative
# 
################################################################################
# This function is for sampling p sets of parameters from a distribution Dist.
samp.par<<-function(Dist,p){
  Iter<-length(Dist[,1]); p.iter<-(sample(seq(Iter),size=p,replace=TRUE))
  return(Dist[p.iter,])
}
# This function extracts shedding profiles given shedding data and ARI information
Profile.shed<-function(shedding,shedding.est,ARI.est,Pathogen.name){
  prof.find<-function(shed.mat,shed.mat.load,ARI){
    shed.p<-which(colSums(shed.mat)>0)
    dur<-c();profile<-list();ari<-c();ari.profile<-list()
    for(i in 1:length(shed.p)){
      end<-which(diff(shed.mat[,shed.p[i]])==-1)
      strt<-which(diff(shed.mat[,shed.p[i]])==1) + 1
      strt<-strt[seq(length(end))]
      if(length(strt)==1){
        # Shedding profile
        profile[[length(dur)+1]]<-shed.mat.load[seq(strt,end),shed.p[i]]
        # ARI profile 
        ari.profile[[length(dur)+1]]<-ARI[seq(strt,end),shed.p[i]]
        ari<-c(ari,as.numeric(sum(ARI[seq(strt,end),shed.p[i]])>0))
      } 
      if(length(strt)>1){
        for(j in 1:length(strt)){
          # Shedding profile
          profile[[length(dur)+j]]<-shed.mat.load[seq(strt[j],end[j]),shed.p[i]]
          # ARI profile
          ari.profile[[length(dur)+j]]<-ARI[seq(strt[j],end[j]),shed.p[i]]
          ari<-c(ari,as.numeric(sum(ARI[seq(strt[j],end[j]),shed.p[i]])>0))
        }
      }
      dur<-c(dur,end-strt+1)
    }
    Res<-list(Profile=profile, Durations=dur, ARI=ari,ARI.profile=ari.profile)
    return(Res)
  }
  ############## All profiles
  Profiles.A<-prof.find(shedding,shedding.est,ARI.est)
  # visualizing
  par(mfrow=c(2,1))
  #symptomatic
  plot(unlist(Profiles.A$Profile[Profiles.A$ARI==1]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' symptomatic'))
  mapply(lines,Profiles.A$Profile[Profiles.A$ARI==1],col=seq_along(Profiles.A$Profile[Profiles.A$ARI==1]),lty=2)
  # asymptomatic
  plot(unlist(Profiles.A$Profile[Profiles.A$ARI==0]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' asymptomatic'))
  mapply(lines,Profiles.A$Profile[Profiles.A$ARI==0],col=seq_along(Profiles.A$Profile[Profiles.A$ARI==0]),lty=2)
  ############## Profiles by age group
  A<-rep(4,N);A[Demo.data$ageyrs<1]<-1
  A[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
  A[(Demo.data$ageyrs>=5) & (Demo.data$ageyrs<15)]<-3
  Profiles.A.1<-prof.find(shedding[,A==1],shedding.est[,A==1],ARI.est[,A==1])
  Profiles.A.2<-prof.find(shedding[,A==2],shedding.est[,A==2],ARI.est[,A==2])
  Profiles.A.3<-prof.find(shedding[,A==3],shedding.est[,A==3],ARI.est[,A==3])
  Profiles.A.4<-prof.find(shedding[,A==4],shedding.est[,A==4],ARI.est[,A==4])
  # visualizing
  par(mfrow=c(2,4))
  #symptomatic
  plot(unlist(Profiles.A.1$Profile[Profiles.A.1$ARI==1]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' symptomatic:<1s'))
  mapply(lines,Profiles.A.1$Profile[Profiles.A.1$ARI==1],col=seq_along(Profiles.A.1$Profile[Profiles.A.1$ARI==1]),lty=2)
  # symptomatic
  plot(unlist(Profiles.A.2$Profile[Profiles.A.2$ARI==1]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' symptomatic:1-5s'))
  mapply(lines,Profiles.A.2$Profile[Profiles.A.2$ARI==1],col=seq_along(Profiles.A.2$Profile[Profiles.A.2$ARI==1]),lty=2)
  # symptomatic
  plot(unlist(Profiles.A.3$Profile[Profiles.A.3$ARI==1]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' symptomatic:5-15s'))
  mapply(lines,Profiles.A.3$Profile[Profiles.A.3$ARI==1],col=seq_along(Profiles.A.3$Profile[Profiles.A.3$ARI==1]),lty=2)
  # symptomatic
  plot(unlist(Profiles.A.4$Profile[Profiles.A.4$ARI==1]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' symptomatic:>15s'))
  mapply(lines,Profiles.A.4$Profile[Profiles.A.4$ARI==1],col=seq_along(Profiles.A.4$Profile[Profiles.A.4$ARI==1]),lty=2)
  
  # asymptomatic
  plot(unlist(Profiles.A.1$Profile[Profiles.A.1$ARI==0]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' asymptomatic:<1s'))
  mapply(lines,Profiles.A.1$Profile[Profiles.A.1$ARI==0],col=seq_along(Profiles.A.1$Profile[Profiles.A.1$ARI==0]),lty=2)
  # asymptomatic
  plot(unlist(Profiles.A.2$Profile[Profiles.A.2$ARI==0]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' asymptomatic:1-5s'))
  mapply(lines,Profiles.A.2$Profile[Profiles.A.2$ARI==0],col=seq_along(Profiles.A.2$Profile[Profiles.A.2$ARI==0]),lty=2)
  # asymptomatic
  plot(unlist(Profiles.A.3$Profile[Profiles.A.3$ARI==0]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' asymptomatic:5-15s'))
  mapply(lines,Profiles.A.3$Profile[Profiles.A.3$ARI==0],col=seq_along(Profiles.A.3$Profile[Profiles.A.3$ARI==0]),lty=2)
  # asymptomatic
  plot(unlist(Profiles.A.4$Profile[Profiles.A.4$ARI==0]),type="n",xlim=c(1,50),ylim=c(2,8.5),ylab='V.load',xlab='Days',main=paste(Pathogen.name,' asymptomatic:>15s'))
  mapply(lines,Profiles.A.4$Profile[Profiles.A.4$ARI==0],col=seq_along(Profiles.A.4$Profile[Profiles.A.4$ARI==0]),lty=2)
  
  Res<-list(All=Profiles.A,age.1=Profiles.A.1,age.2=Profiles.A.2,age.3=Profiles.A.3,
            age.4=Profiles.A.4)
  return(Res)
}
# given a set of parameters(theta), and the number of simulations to run(runs),
# this function simulates epidemics
Set.sim<-function(theta,runs){
  # Given a set of parameters theta and participant demographics, this function 
  # simulates a single epidemic
  IBM.simulator<-function(){
    # This function calculates the probability of infection at a time t given parameter set theta
    Prob.Inf.func<-function(t){
      # _________________1. Relative Susceptibility_________________________________
      # Effect of previous homologous and heterologous infection on susceptibility
      hist.par.A<-c(0,theta[['Prev.hom']], theta[['Prev.het']],(theta[['Prev.hom']]+theta[['Prev.het']]))
      hist.mat.A<-matrix(hist.par.A[B[t,]],1,N,byrow=FALSE)
      hist.par.B<-c(0,theta[['Prev.het']], theta[['Prev.hom']],(theta[['Prev.hom']]+theta[['Prev.het']]))
      hist.mat.B<-matrix(hist.par.B[B[t,]],1,N,byrow=FALSE)
      # Effect of age on susceptibility
      age.par<-c(0,theta[3:(2+s.g-1)]);age.mat<-matrix(age.par[A],1,N,byrow=TRUE)
      # Combining age and infection history effect on susceptibility
      PI.A<-exp(hist.mat.A+age.mat);PI.A[S.a[t,]==0]<-0
      PI.B<-exp(hist.mat.B+age.mat);PI.B[S.b[t,]==0]<-0
      # _________________2. Rate of Exposure________________________________________
      # Within household transmission coefficients
      # ---------------------------------
      # Within household transmission coefficient
      Eta.A<-exp(theta[['eta.A']]); Eta.B<-exp(theta[['eta.B']])
      # Effect of household size on within household transmission
      hh.par<-c(1,exp(theta[['hh.size']]));hh.mat<-matrix(hh.par[H],1,N,byrow=TRUE)
      # Effect of viral load and symptoms on infectivity
      V.par<-c(0,1,exp(theta[['LowSym']]),exp(theta[['HighSym']]))
      Infectivity.A<-matrix(V.par[Load.n.sym.a[t,]],1,N,byrow=FALSE)*status[t,] # *status[t,] so that absent shedders are not included
      HH.risk.A<-t(Cc %*% t(Infectivity.A))
      Infectivity.B<-matrix(V.par[Load.n.sym.b[t,]],1,N,byrow=FALSE)*status[t,]
      HH.risk.B<-t(Cc %*% t(Infectivity.B))
      # Total within household risk, status identifies if an indiviual is present in the household
      in.house.A<-(status[t,]*hh.mat)*(Eta.A*HH.risk.A) 
      in.house.B<-(status[t,]*hh.mat)*(Eta.B*HH.risk.B) 
      # Community rate of exposure
      # --------------------------
      # Community transmission coefficient
      Eps.A<-exp(theta[['epsilon.A']]);Eps.B<-exp(theta[['epsilon.B']])
      # Effect of age on community level exposure
      x<-length(theta)-e.g;Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
      comm.age<-matrix(Eps.age[E],1,N)
      # Time varying component of community exposure/background density funtion
      comm.time.A<-matrix(Comm.risk.a[t],1,1);comm.time.B<-matrix(Comm.risk.b[t],1,1)
      # Combining the age effects and time varying function
      Comm.A <- Eps.A*(comm.time.A %*% comm.age)  
      Comm.B <- Eps.B*(comm.time.B %*% comm.age)
      # Total community exposure rate 
      Lambda.A<-in.house.A + Comm.A
      Lambda.B<-in.house.B + Comm.B
      # _________________3. Probability of Infection________________________________
      Alpha.A<-(1-exp(-PI.A*Lambda.A));Alpha.B<-(1-exp(-PI.B*Lambda.B))
      Res<-data.frame(Prob.a=Alpha.A[1,],Prob.b=Alpha.B[1,])
      return(Res)
    }
    # Given the number of transmission events, this function finds people involved 
    Person.finder<-function(No.events,cum.rate,pool.pple){
      # at some points number of proposed events exceeds the pool of pple available
      # to experience them, so need to change number of people to be 'found'
      if(No.events<length(pool.pple)){
        No.pple<-No.events
        pple.involved<-c();
        while(length(unique(pple.involved))<No.pple){
          Random.no<-runif(1,0,1);
          prsn<-which(((c(0,cum.rate[-length(pool.pple)]) < Random.no*cum.rate[length(pool.pple)])
                       & (Random.no*cum.rate[length(pool.pple)] <= cum.rate)) ==TRUE)
          pple.involved<-c(prsn,pple.involved)
        }
        return(unique(pple.involved))
        # 'unique' function is in case different random numbers support the same person
      }else{
        return(pool.pple)
      }
    }
    # This function randomly assigns shedding profiles, i.e duration, load and ARI status
    shed.profiler<-function(Profiles,age.grp){
      if(age.grp==1){
        take<-sample(seq(length(Profiles$age.1$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.1$Profile[[take]],            # viral load
                   ari=Profiles$age.1$ARI.profile[[take]])         # ARI
      }
      if(age.grp==2){
        take<-sample(seq(length(Profiles$age.2$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.2$Profile[[take]],            # viral load
                   ari=Profiles$age.2$ARI.profile[[take]])         # ARI
      }
      if(age.grp==3){
        take<-sample(seq(length(Profiles$age.3$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.3$Profile[[take]],            # viral load
                   ari=Profiles$age.3$ARI.profile[[take]])         # ARI
      }
      if(age.grp==4){
        take<-sample(seq(length(Profiles$age.4$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.4$Profile[[take]],            # viral load
                   ari=Profiles$age.4$ARI.profile[[take]])         # ARI
      }  
      return(prof)
    }
    # ---------- initiate system
    T.a<-rep(0,D);T.b<-rep(0,D)
    Nms<-list(rows=paste('t',seq(D),sep=""),cols=paste('p',seq(N),sep=""))
    # Binary indicator variable of susceptibility:1 => susceptible
    S.a<<-matrix(1,nrow = D,ncol = N,dimnames=Nms);S.b<<-S.a
    # Binary indicator variables of exposure state
    E.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms); E.b<-E.a; 
    # Continuous variables of viral load
    I.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);I.b<-I.a
    # Categorical variables of shedding status
    Sh.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);Sh.b<-Sh.a
    # Categorical variables of ARI status
    ARI.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);ARI.b<-I.a
    # Categorical variables of infectious group
    Load.n.sym.a<<-matrix(1,nrow = D,ncol = N,dimnames=Nms);Load.n.sym.b<<-Load.n.sym.a
    # Categorical variable tracking infection history
    postInf.a<<-matrix(0,nrow = D,ncol = N,dimnames=Nms);postInf.b<<-postInf.a
    B<<-matrix(1,nrow = D,ncol = N,dimnames=Nms)
    # integer variable for latency duration at exposure
    Inf.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);Inf.b<-Inf.a
    
    # Continuous variables for prob. infection (at time t=1, no one from HHs in shedding)
    alpha.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);alpha.b<-alpha.a;
    alpha.a[1,]<-Prob.Inf.func(1)[['Prob.a']];alpha.b[1,]<-Prob.Inf.func(1)[['Prob.b']]
    # ---------- Iteration per time step
    for(t in 2:D){
      # Probability of transmission at time t
      alpha.a[t,]<-Prob.Inf.func(t)[['Prob.a']];alpha.b[t,]<-Prob.Inf.func(t)[['Prob.b']]
      
      # 1. Number of TRANSMISSION events at time t
      T.a[t]<-rpois(n=1,lambda=sum(alpha.a[t,])) # RSV A
      T.b[t]<-rpois(n=1,lambda=sum(alpha.b[t,])) # RSV B
      # susceptibles at current time step (prior to any infection events)
      Sus.a<-which(S.a[t-1,]==1);Sus.b<-which(S.b[t-1,]==1)
      # if infection events did occur, who was involved
      if(T.a[t]>0){
        sum.alpha<-as.numeric(cumsum(alpha.a[t,Sus.a]))
        Trans.a<-Person.finder(T.a[t],sum.alpha,Sus.a)
        # 1. Assign latency durations
        y<-Sus.a[Trans.a];x<-A2[y]
        # Latency durations
        Inf.a[t,y]<-sample(Lat, size = length(Trans.a), replace = TRUE)
        # 2. Updating state variables over time
        for (i in 1:length(Trans.a)){
          # a) Exposure
          d.e<-seq(t,t+Inf.a[t,y[i]]-1);d.e<-d.e[d.e<=D];E.a[d.e,y[i]]<-1
          # b) Viral load and Infectiousness category
          # shedding profile                                                            
          prof<-shed.profiler(Profiles.rsv.a,age.grp=x[i]) # Random assignment of profiles
          d.sh<-seq((t+Inf.a[t,y[i]]),(t+Inf.a[t,y[i]]+length(prof$shed)-1)) # Time window of infection
          d.sh<-d.sh[d.sh<=D];# modifying window incase it goes beyond the study period
          prof$shed<-prof$shed[seq(length(d.sh))];prof$ari<-prof$ari[seq(length(d.sh))]
          # --viral load and ARI
          I.a[d.sh,y[i]]<-prof$shed;ARI.a[d.sh,y[i]]<-prof$ari   
          # --infectiousness category
          # Reference group; asymptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed<=high.load & prof$ari==0]<-2                      
          Load.n.sym.a[d.sh,y[i]][prof$shed>high.load & prof$ari==0]<-2                      
          # low viral load and symptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed<=high.load & prof$ari==1]<-3                      
          # high viral load and symptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed>high.load & prof$ari==1]<-4              
          # c) Susceptible 
          d.s<-seq(t,(t+Inf.a[t,y[i]]+length(prof$shed)-1));d.s<-d.s[d.s<=D];S.a[d.s,y[i]]<-0
        }
      }
      if(T.b[t]>0){
        sum.alpha<-as.numeric(cumsum(alpha.b[t,Sus.b]))
        Trans.b<-Person.finder(T.b[t],sum.alpha,Sus.b)
        # 1. Assign latency durations
        y<-Sus.b[Trans.b];x<-A2[y]
        # Latency durations
        Inf.b[t,y]<-sample(Lat, size = length(Trans.b), replace = TRUE)
        # 2. Updating state variables over time
        for (i in 1:length(Trans.b)){
          # a) Exposure
          d.e<-seq(t,t+Inf.b[t,y[i]]-1);d.e<-d.e[d.e<=D];E.b[d.e,y[i]]<-1
          # b) Viral load and Infectiousness category
          # shedding profile                                                            
          prof<-shed.profiler(Profiles.rsv.b,age.grp=x[i]) # Random assignment of profiles
          d.sh<-seq((t+Inf.b[t,y[i]]),(t+Inf.b[t,y[i]]+length(prof$shed)-1)) # Time window of infection
          d.sh<-d.sh[d.sh<=D];# modifying window incase it goes beyond the study period
          prof$shed<-prof$shed[seq(length(d.sh))];prof$ari<-prof$ari[seq(length(d.sh))]
          # --viral load and ARI
          I.b[d.sh,y[i]]<-prof$shed;ARI.b[d.sh,y[i]]<-prof$ari   
          # --infectiousness category
          # Reference group; asymptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed<=high.load & prof$ari==0]<-2                      
          Load.n.sym.b[d.sh,y[i]][prof$shed>high.load & prof$ari==0]<-2                      
          # low viral load and symptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed<=high.load & prof$ari==1]<-3                      
          # high viral load and symptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed>high.load & prof$ari==1]<-4             
          # c) Susceptible 
          d.s<-seq(t,(t+Inf.b[t,y[i]]+length(prof$shed)-1));d.s<-d.s[d.s<=D];S.b[d.s,y[i]]<-0
        }
      }
      # Updating infection history matrix 
      if(t>2){ # for computation purposes(either way B[=2] cannot be >1 since at t=1 everyone in susc)
        # previously shed A but have recovered
        b1<-which(colSums(I.a[seq(t-1),])>0 & (I.a[t,]==0))
        postInf.a[seq(t,D),b1]<-1
        # previously shed B but have recovered
        b2<-which(colSums(I.b[seq(t-1),])>0 & (I.b[t,]==0))
        postInf.b[seq(t,D),b2]<-1
        # The infection history matrix
        B[postInf.a==1]<-2;B[postInf.b==1]<-3
        B[(postInf.a==1 & postInf.b==1)]<-4 
      }
    }
    # Results of a single simulation
    Results<-list(RSV.A=I.a, RSV.B=I.b,Inf.hist.mat=B,Load.n.sym.a=Load.n.sym.a,
                  Load.n.sym.b=Load.n.sym.b)
    return(Results)
  }
  # -------------Set up
  # Set up the simulation environment;participant demographics and community epidemic curve
  # -Age covariate for susceptibility 
  A<-rep(4,N);s.g<-4
  A[Demo.data$ageyrs<1]<-1;A[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
  A[(Demo.data$ageyrs>=5) & (Demo.data$ageyrs<15)]<-3
  # -Age Covariate for community risk/exposure
  E<-rep(3,N);e.g<-3
  E[Demo.data$ageyrs<1]<-1;E[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
  # -HH size covariate
  H<-rep(2,N);H[Demo.data$hhsize<8]<-1
  # -An NxN matrix the identifies housemates.Used in FOI function to add up HH shedders
  houses<-unique(Demo.data$hhid)
  Cc<-matrix(0,N,N)
  for (i in 1:length(houses)){
    # Their housemates
    pp<-which(Demo.data$hhid==houses[i])
    Cc[pp,pp]<-1
  }
  diag(Cc)<-0
  # -Exposure duration(same for A and B)
  Lat<-c(rep(2,4),rep(3,4),rep(4,3),5)
  # -Age groups for simulation
  A2<-A
  # -Background community function
  Comm.risk.a<-Shed.data$Comm.risk;Comm.risk.b<-Shed.data.2$Comm.risk
  # -status variable
  status<<-Shed.data$status;
  # -cut-off for high viral load
  high.load<<-6
  # -------------Simulations
  Res.iter<-list()
  for(i in 1:runs){
    Res.iter<-c(Res.iter,IBM.simulator()) 
    message(i,' simulations done')
  } 
  return(Res.iter)
}
# Given the simulated data sets, this function calculates the following outcome
# measures; total number of people infected, total number of househols infected
# proportion of infected people with multiple onsets and timing of epidemic peak
Outcome.measures<-function(runs,Res.iter){
  # RSV A
  RSV.a<-array(0,c(D,N,runs));a.s<-c(seq(1,(runs*5)-4,5))
  # * status so that individuals who are shedding outside of the household are not
  # counted so that we compare more similar things (NB the data only saw shedding 
  # times if the person was present for sampling.Due to out imputaion method, there
  # are 3 onsets in RSV B that occurred while people were away, but that will not
  # give a very big difference in total numbers, but including absent simulated
  # shedders might )
  for(j in 1:length(a.s)) RSV.a[,,j]<-Res.iter[[a.s[j]]] * status 
  iter.a<-0*RSV.a;iter.a[which(RSV.a>0,TRUE)]<-1
  Shedding.iter.a<-apply(iter.a,3,rowSums) # epidemic curve
  # RSV B
  RSV.b<-array(0,c(D,N,runs));b.s<-c(seq(2,(runs*5)-3,5))
  for(j in 1:length(b.s)) RSV.b[,,j]<-Res.iter[[b.s[j]]] * status
  iter.b<-0*RSV.b;iter.b[which(RSV.b>0,TRUE)]<-1
  Shedding.iter.b<-apply(iter.b,3,rowSums) # epidemic curve
  # RSV
  iter.rsv<-iter.a+iter.b;iter.rsv[which(iter.rsv>0,TRUE)]<-1
  Shedding.iter.rsv<-apply(iter.rsv,3,rowSums) # epidemic curve
  # function counting total numbers for each iteration
  tot.inf<-function(it) {
    rsv.a.pple=which(colSums(iter.a[,,it])>0);
    rsv.b.pple=which(colSums(iter.b[,,it])>0)
    rsv.a.hh<-unique(Demo.data[['hhid']][rsv.a.pple])
    rsv.b.hh<-unique(Demo.data[['hhid']][rsv.b.pple])
    #Onsets
    strt.a<-which(diff(iter.a[,,it])==1,TRUE);strt.a[,1]<-strt.a[,1]+1
    Ons.a<-0*iter.a[,,it] ;Ons.a[strt.a]<-1
    prop.a.m<-length(which(colSums(Ons.a)>1))/length(rsv.a.pple)
    
    strt.b<-which(diff(iter.b[,,it])==1,TRUE);strt.b[,1]<-strt.b[,1]+1
    Ons.b<-0*iter.b[,,it] ;Ons.b[strt.b]<-1
    prop.b.m<-length(which(colSums(Ons.b)>1))/length(rsv.b.pple)
    
    strt.rsv<-which(diff(iter.rsv[,,it])==1,TRUE);strt.rsv[,1]<-strt.rsv[,1]+1
    Ons.rsv<-0*iter.rsv[,,it] ;Ons.rsv[strt.rsv]<-1
    prop.rsv.m<-length(which(colSums(Ons.rsv)>1))/length(union(rsv.a.pple,rsv.b.pple))
    
    Res<-c(a.pple=length(rsv.a.pple), # Toatl number of people with RSV A
           b.pple=length(rsv.b.pple), # Total number of people with RSV B
           rsv.pple=length(union(rsv.a.pple,rsv.b.pple)), # Total number of people with RSV
           a.hh=length(rsv.a.hh),  # Total number of households with RSV A
           b.hh=length(rsv.b.hh),  # Total number of households with RSV B
           rsv.hh=length(union(rsv.a.hh,rsv.b.hh)), # Total number of households with RSV
           prop.a.m=prop.a.m, # Proportion of people who had multiple onsets of RSV A
           prop.b.m=prop.b.m,  # Proportion of infected people who had multiple onsets of RSV B
           prop.rsv.m=prop.rsv.m) # Proportion of infected people who had multiple onsets of RSV
    
    return(Res)
  }
  total.infected<-t(sapply(seq(runs),tot.inf))
  # Funtion getting the peak epidemic timings
  peak.inf<-function(it) {
    p.a<-min(which(Shedding.iter.a[,it]==max(Shedding.iter.a[,it])))
    p.b<-min(which(Shedding.iter.b[,it]==max(Shedding.iter.b[,it])))
    p.rsv<-min(which(Shedding.iter.rsv[,it]==max(Shedding.iter.rsv[,it])))
    return(c(peak.a=p.a,peak.b=p.b,peak.rsv=p.rsv))
  }
  peak.tim<-t(sapply(seq(runs),peak.inf))
  Res<-as.data.frame(cbind(total.infected,peak.tim))
  return(Res)
}
# Counts total number infected by age group for <1, 1-5, 5-15, >15
Outcome.measures.2<-function(runs,Res.iter){
  # RSV A
  RSV.a<-array(0,c(D,N,runs));a.s<-c(seq(1,(runs*5)-4,5))
  for(j in 1:length(a.s)) RSV.a[,,j]<-Res.iter[[a.s[j]]] * status
  iter.a<-0*RSV.a;iter.a[which(RSV.a>0,TRUE)]<-1
  Shedding.iter.a<-apply(iter.a,3,rowSums) # epidemic curve
  # RSV B
  RSV.b<-array(0,c(D,N,runs));b.s<-c(seq(2,(runs*5)-3,5))
  for(j in 1:length(b.s)) RSV.b[,,j]<-Res.iter[[b.s[j]]] * status
  iter.b<-0*RSV.b;iter.b[which(RSV.b>0,TRUE)]<-1
  Shedding.iter.b<-apply(iter.b,3,rowSums) # epidemic curve
  # RSV
  iter.rsv<-iter.a+iter.b;iter.rsv[which(iter.rsv>0,TRUE)]<-1
  Shedding.iter.rsv<-apply(iter.rsv,3,rowSums) # epidemic curve
  # function counting total numbers for each iteration
  tot.inf<-function(it) {
    rsv.a.pple.1<-which(colSums(iter.a[,A==1,it])>0);rsv.b.pple.1<-which(colSums(iter.b[,A==1,it])>0)
    rsv.a.pple.2<-which(colSums(iter.a[,A==2,it])>0);rsv.b.pple.2<-which(colSums(iter.b[,A==2,it])>0)
    rsv.a.pple.3<-which(colSums(iter.a[,A==3,it])>0);rsv.b.pple.3<-which(colSums(iter.b[,A==3,it])>0)
    rsv.a.pple.4<-which(colSums(iter.a[,A==4,it])>0);rsv.b.pple.4<-which(colSums(iter.b[,A==4,it])>0)
    Res<-c(rsv.pple.1=length(union(rsv.a.pple.1,rsv.b.pple.1)),
           rsv.pple.2=length(union(rsv.a.pple.2,rsv.b.pple.2)),
           rsv.pple.3=length(union(rsv.a.pple.3,rsv.b.pple.3)),
           rsv.pple.4=length(union(rsv.a.pple.4,rsv.b.pple.4))) 
    return(Res)
  }
  
  total.infected<-t(sapply(seq(runs),tot.inf))
  
  Res<-as.data.frame(total.infected)
  return(Res)
}
# ---- Load the real data
source('Data_Prep.R')

###################################### PART 1: #################################
# ---- Given the real data, extract the shedding profiles
Profiles.rsv.a<-Profile.shed(Shed.data$Shedding,Shed.data$Shedding.est.load,ARI.estimated,'RSV A')
Profiles.rsv.b<-Profile.shed(Shed.data.2$Shedding,Shed.data.2$Shedding.est.load,ARI.estimated.2,'RSV B')
# -Age covariate for susceptibility
A<-rep(4,N);s.g<-4
A[Demo.data$ageyrs<1]<-1;A[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
A[(Demo.data$ageyrs>=5) & (Demo.data$ageyrs<15)]<-3

# ---- Load the posterior distribution 
load("Para_distributions/Res.group.diag.RData")
theta.names<-names(Res.group.diag$Res)[-c(9,16)]
Par.dist<-Res.group.diag$Res[,theta.names]
# ---- Parameters
# Sampling from the posterior
paras.set1<-(samp.par(Par.dist,1));paras.set1<-as.numeric(paras.set1);
names(paras.set1)<-names(Par.dist)
# Checking where the parameter selected falls in the distribution
par(mfrow=c(3,5))
for(i in 1:length(paras.set1)){
  plot(density(Par.dist[,i]),type='l',lwd=2,main=names(Par.dist)[i],
       xlab=names(Par.dist)[i])
  abline(v=paras.set1[i],col='green',lwd=2,lty=2)
}
# ---- Number of simulations to run
sims<-200

# ---- Simulations & outcome measures: Scenario 1, relative infectivity as estimated
# Original relative infectiousness
set1<-Set.sim(paras.set1,sims);outcome1<-Outcome.measures(sims,set1)
outcome1.age<-Outcome.measures.2(sims,set1)

# ---- Simulations & outcome measures: Scenario 2, Reducing infectivity of symptomatics
# Changing the relative infectiousness, such that symptomatic have the same baseline
# infectivity as asymptomatics
paras.set1.2<-paras.set1
paras.set1.2[['LowSym']]<-0;paras.set1.2[['HighSym']]<-0
set1.2<-Set.sim(paras.set1.2,sims);outcome1.2<-Outcome.measures(sims,set1.2)
outcome1.2.age<-Outcome.measures.2(sims,set1.2)

###################################### PART 2: #################################
# Assuming it is possible to get asymptomatic infections but that such individuals
# do not transmit
#-------------------------------------------------------------------------------
# Only difference from Set.sim is that the infectiousness of the reference group
# has been change to 0 rather than 1
Set.sim.2<-function(theta,runs){
  # Given a set of parameters theta and participant demographics, this function 
  # simulates a single epidemic
  IBM.simulator<-function(){
    # This function calculates the probability of infection at a time t given parameter set theta
    Prob.Inf.func<-function(t){
      # _________________1. Relative Susceptibility_________________________________
      # Effect of previous homologous and heterologous infection on susceptibility
      hist.par.A<-c(0,theta[['Prev.hom']], theta[['Prev.het']],(theta[['Prev.hom']]+theta[['Prev.het']]))
      hist.mat.A<-matrix(hist.par.A[B[t,]],1,N,byrow=FALSE)
      hist.par.B<-c(0,theta[['Prev.het']], theta[['Prev.hom']],(theta[['Prev.hom']]+theta[['Prev.het']]))
      hist.mat.B<-matrix(hist.par.B[B[t,]],1,N,byrow=FALSE)
      # Effect of age on susceptibility
      age.par<-c(0,theta[3:(2+s.g-1)]);age.mat<-matrix(age.par[A],1,N,byrow=TRUE)
      # Combining age and infection history effect on susceptibility
      PI.A<-exp(hist.mat.A+age.mat);PI.A[S.a[t,]==0]<-0
      PI.B<-exp(hist.mat.B+age.mat);PI.B[S.b[t,]==0]<-0
      # _________________2. Rate of Exposure________________________________________
      # Within household transmission coefficients
      # ---------------------------------
      # Within household transmission coefficient
      Eta.A<-exp(theta[['eta.A']]); Eta.B<-exp(theta[['eta.B']])
      # Effect of household size on within household transmission
      hh.par<-c(1,exp(theta[['hh.size']]));hh.mat<-matrix(hh.par[H],1,N,byrow=TRUE)
      # Effect of viral load and symptoms on infectivity
      V.par<-c(0,0,exp(theta[['LowSym']]),exp(theta[['HighSym']]))
      Infectivity.A<-matrix(V.par[Load.n.sym.a[t,]],1,N,byrow=FALSE)*status[t,] # *status[t,] so that absent shedders are not included
      HH.risk.A<-t(Cc %*% t(Infectivity.A))
      Infectivity.B<-matrix(V.par[Load.n.sym.b[t,]],1,N,byrow=FALSE)*status[t,]
      HH.risk.B<-t(Cc %*% t(Infectivity.B))
      # Total within household risk, status identifies if an indiviual is present in the household
      in.house.A<-(status[t,]*hh.mat)*(Eta.A*HH.risk.A) 
      in.house.B<-(status[t,]*hh.mat)*(Eta.B*HH.risk.B) 
      # Community rate of exposure
      # --------------------------
      # Community transmission coefficient
      Eps.A<-exp(theta[['epsilon.A']]);Eps.B<-exp(theta[['epsilon.B']])
      # Effect of age on community level exposure
      x<-length(theta)-e.g;Eps.age<-c(1,exp(theta[(x+2):length(theta)]))
      comm.age<-matrix(Eps.age[E],1,N)
      # Time varying component of community exposure/background density funtion
      comm.time.A<-matrix(Comm.risk.a[t],1,1);comm.time.B<-matrix(Comm.risk.b[t],1,1)
      # Combining the age effects and time varying function
      Comm.A <- Eps.A*(comm.time.A %*% comm.age)  
      Comm.B <- Eps.B*(comm.time.B %*% comm.age)
      # Total community exposure rate 
      Lambda.A<-in.house.A + Comm.A
      Lambda.B<-in.house.B + Comm.B
      # _________________3. Probability of Infection________________________________
      Alpha.A<-(1-exp(-PI.A*Lambda.A));Alpha.B<-(1-exp(-PI.B*Lambda.B))
      Res<-data.frame(Prob.a=Alpha.A[1,],Prob.b=Alpha.B[1,])
      return(Res)
    }
    # Given the number of transmission events, this function finds people involved 
    Person.finder<-function(No.events,cum.rate,pool.pple){
      # at some points number of proposed events exceeds the pool of pple available
      # to experience them, so need to change number of people to be 'found'
      if(No.events<length(pool.pple)){
        No.pple<-No.events
        pple.involved<-c();
        while(length(unique(pple.involved))<No.pple){
          Random.no<-runif(1,0,1);
          prsn<-which(((c(0,cum.rate[-length(pool.pple)]) < Random.no*cum.rate[length(pool.pple)])
                       & (Random.no*cum.rate[length(pool.pple)] <= cum.rate)) ==TRUE)
          pple.involved<-c(prsn,pple.involved)
        }
        return(unique(pple.involved))
        # 'unique' function is in case different random numbers support the same person
      }else{
        return(pool.pple)
      }
    }
    # This function randomly assigns shedding profiles, i.e duration, load and ARI status
    shed.profiler<-function(Profiles,age.grp){
      if(age.grp==1){
        take<-sample(seq(length(Profiles$age.1$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.1$Profile[[take]],            # viral load
                   ari=Profiles$age.1$ARI.profile[[take]])         # ARI
      }
      if(age.grp==2){
        take<-sample(seq(length(Profiles$age.2$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.2$Profile[[take]],            # viral load
                   ari=Profiles$age.2$ARI.profile[[take]])         # ARI
      }
      if(age.grp==3){
        take<-sample(seq(length(Profiles$age.3$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.3$Profile[[take]],            # viral load
                   ari=Profiles$age.3$ARI.profile[[take]])         # ARI
      }
      if(age.grp==4){
        take<-sample(seq(length(Profiles$age.4$Profile)),size = 1) # index
        prof<-list(shed=Profiles$age.4$Profile[[take]],            # viral load
                   ari=Profiles$age.4$ARI.profile[[take]])         # ARI
      }  
      return(prof)
    }
    # ---------- initiate system
    T.a<-rep(0,D);T.b<-rep(0,D)
    Nms<-list(rows=paste('t',seq(D),sep=""),cols=paste('p',seq(N),sep=""))
    # Binary indicator variable of susceptibility:1 => susceptible
    S.a<<-matrix(1,nrow = D,ncol = N,dimnames=Nms);S.b<<-S.a
    # Binary indicator variables of exposure state
    E.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms); E.b<-E.a; 
    # Continuous variables of viral load
    I.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);I.b<-I.a
    # Categorical variables of shedding status
    Sh.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);Sh.b<-Sh.a
    # Categorical variables of ARI status
    ARI.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);ARI.b<-I.a
    # Categorical variables of infectious group
    Load.n.sym.a<<-matrix(1,nrow = D,ncol = N,dimnames=Nms);Load.n.sym.b<<-Load.n.sym.a
    # Categorical variable tracking infection history
    postInf.a<<-matrix(0,nrow = D,ncol = N,dimnames=Nms);postInf.b<<-postInf.a
    B<<-matrix(1,nrow = D,ncol = N,dimnames=Nms)
    # integer variable for latency duration at exposure
    Inf.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);Inf.b<-Inf.a
    
    # Continuous variables for prob. infection (at time t=1, no one from HHs in shedding)
    alpha.a<-matrix(0,nrow = D,ncol = N,dimnames=Nms);alpha.b<-alpha.a;
    alpha.a[1,]<-Prob.Inf.func(1)[['Prob.a']];alpha.b[1,]<-Prob.Inf.func(1)[['Prob.b']]
    # ---------- Iteration per time step
    for(t in 2:D){
      # Probability of transmission at time t
      alpha.a[t,]<-Prob.Inf.func(t)[['Prob.a']];alpha.b[t,]<-Prob.Inf.func(t)[['Prob.b']]
      
      # 1. Number of TRANSMISSION events at time t
      T.a[t]<-rpois(n=1,lambda=sum(alpha.a[t,])) # RSV A
      T.b[t]<-rpois(n=1,lambda=sum(alpha.b[t,])) # RSV B
      # susceptibles at current time step (prior to any infection events)
      Sus.a<-which(S.a[t-1,]==1);Sus.b<-which(S.b[t-1,]==1)
      # if infection events did occur, who was involved
      if(T.a[t]>0){
        sum.alpha<-as.numeric(cumsum(alpha.a[t,Sus.a]))
        Trans.a<-Person.finder(T.a[t],sum.alpha,Sus.a)
        # 1. Assign latency durations
        y<-Sus.a[Trans.a];x<-A2[y]
        # Latency durations
        Inf.a[t,y]<-sample(Lat, size = length(Trans.a), replace = TRUE)
        # 2. Updating state variables over time
        for (i in 1:length(Trans.a)){
          # a) Exposure
          d.e<-seq(t,t+Inf.a[t,y[i]]-1);d.e<-d.e[d.e<=D];E.a[d.e,y[i]]<-1
          # b) Viral load and Infectiousness category
          # shedding profile                                                            
          prof<-shed.profiler(Profiles.rsv.a,age.grp=x[i]) # Random assignment of profiles
          d.sh<-seq((t+Inf.a[t,y[i]]),(t+Inf.a[t,y[i]]+length(prof$shed)-1)) # Time window of infection
          d.sh<-d.sh[d.sh<=D];# modifying window incase it goes beyond the study period
          prof$shed<-prof$shed[seq(length(d.sh))];prof$ari<-prof$ari[seq(length(d.sh))]
          # --viral load and ARI
          I.a[d.sh,y[i]]<-prof$shed;ARI.a[d.sh,y[i]]<-prof$ari   
          # --infectiousness category
          # Reference group; asymptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed<=high.load & prof$ari==0]<-2                      
          Load.n.sym.a[d.sh,y[i]][prof$shed>high.load & prof$ari==0]<-2                      
          # low viral load and symptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed<=high.load & prof$ari==1]<-3                      
          # high viral load and symptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed>high.load & prof$ari==1]<-4              
          # c) Susceptible 
          d.s<-seq(t,(t+Inf.a[t,y[i]]+length(prof$shed)-1));d.s<-d.s[d.s<=D];S.a[d.s,y[i]]<-0
        }
      }
      if(T.b[t]>0){
        sum.alpha<-as.numeric(cumsum(alpha.b[t,Sus.b]))
        Trans.b<-Person.finder(T.b[t],sum.alpha,Sus.b)
        # 1. Assign latency durations
        y<-Sus.b[Trans.b];x<-A2[y]
        # Latency durations
        Inf.b[t,y]<-sample(Lat, size = length(Trans.b), replace = TRUE)
        # 2. Updating state variables over time
        for (i in 1:length(Trans.b)){
          # a) Exposure
          d.e<-seq(t,t+Inf.b[t,y[i]]-1);d.e<-d.e[d.e<=D];E.b[d.e,y[i]]<-1
          # b) Viral load and Infectiousness category
          # shedding profile                                                            
          prof<-shed.profiler(Profiles.rsv.b,age.grp=x[i]) # Random assignment of profiles
          d.sh<-seq((t+Inf.b[t,y[i]]),(t+Inf.b[t,y[i]]+length(prof$shed)-1)) # Time window of infection
          d.sh<-d.sh[d.sh<=D];# modifying window incase it goes beyond the study period
          prof$shed<-prof$shed[seq(length(d.sh))];prof$ari<-prof$ari[seq(length(d.sh))]
          # --viral load and ARI
          I.b[d.sh,y[i]]<-prof$shed;ARI.b[d.sh,y[i]]<-prof$ari   
          # --infectiousness category
          # Reference group; asymptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed<=high.load & prof$ari==0]<-2                      
          Load.n.sym.b[d.sh,y[i]][prof$shed>high.load & prof$ari==0]<-2                      
          # low viral load and symptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed<=high.load & prof$ari==1]<-3                      
          # high viral load and symptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed>high.load & prof$ari==1]<-4             
          # c) Susceptible 
          d.s<-seq(t,(t+Inf.b[t,y[i]]+length(prof$shed)-1));d.s<-d.s[d.s<=D];S.b[d.s,y[i]]<-0
        }
      }
      # Updating infection history matrix 
      if(t>2){ # for computation purposes(either way B[=2] cannot be >1 since at t=1 everyone in susc)
        # previously shed A but have recovered
        b1<-which(colSums(I.a[seq(t-1),])>0 & (I.a[t,]==0))
        postInf.a[seq(t,D),b1]<-1
        # previously shed B but have recovered
        b2<-which(colSums(I.b[seq(t-1),])>0 & (I.b[t,]==0))
        postInf.b[seq(t,D),b2]<-1
        # The infection history matrix
        B[postInf.a==1]<-2;B[postInf.b==1]<-3
        B[(postInf.a==1 & postInf.b==1)]<-4 
      }
    }
    # Results of a single simulation
    Results<-list(RSV.A=I.a, RSV.B=I.b,Inf.hist.mat=B,Load.n.sym.a=Load.n.sym.a,
                  Load.n.sym.b=Load.n.sym.b)
    return(Results)
  }
  # -------------Set up
  # Set up the simulation environment;participant demographics and community epidemic curve
  # -Age covariate for susceptibility 
  A<-rep(4,N);s.g<-4
  A[Demo.data$ageyrs<1]<-1;A[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
  A[(Demo.data$ageyrs>=5) & (Demo.data$ageyrs<15)]<-3
  # -Age Covariate for community risk/exposure
  E<-rep(3,N);e.g<-3
  E[Demo.data$ageyrs<1]<-1;E[(Demo.data$ageyrs>=1) & (Demo.data$ageyrs<5)]<-2
  # -HH size covariate
  H<-rep(2,N);H[Demo.data$hhsize<8]<-1
  # -An NxN matrix the identifies housemates.Used in FOI function to add up HH shedders
  houses<-unique(Demo.data$hhid)
  Cc<-matrix(0,N,N)
  for (i in 1:length(houses)){
    # Their housemates
    pp<-which(Demo.data$hhid==houses[i])
    Cc[pp,pp]<-1
  }
  diag(Cc)<-0
  # -Exposure duration(same for A and B)
  Lat<-c(rep(2,4),rep(3,4),rep(4,3),5)
  # -Age groups for simulation
  A2<-A
  # -Background community function
  Comm.risk.a<-Shed.data$Comm.risk;Comm.risk.b<-Shed.data.2$Comm.risk
  # -status variable
  status<<-Shed.data$status;
  # -cut-off for high viral load
  high.load<<-6
  # -------------Simulations
  Res.iter<-list()
  for(i in 1:runs){
    Res.iter<-c(Res.iter,IBM.simulator()) 
    message(i,' simulations done')
  } 
  return(Res.iter)
}

# ---- Simulations & outcome measures: Scenario 3, asymptomatics who do not transmit

# Changing the relative infectiousness, such that symptomatic low shedders are
# the reference group
paras.set1.3<-paras.set1
paras.set1.3[['HighSym']]<-paras.set1.3[['HighSym']]-paras.set1.3[['LowSym']]
paras.set1.3[['LowSym']]<-0; # Symptomatic low load shedders are reference group

set1.3<-Set.sim.2(paras.set1.3,sims);outcome1.3<-Outcome.measures(sims,set1.3)
outcome1.3.age<-Outcome.measures.2(sims,set1.3)

# ---- Comparing 3 scenarios 
# Numbers infected per age group from the observed data
Res<-c()
Res[1]<-sum((colSums(Shed.data.3$Onset)[A==1])>0) 
Res[2]<-sum((colSums(Shed.data.3$Onset)[A==2])>0)
Res[3]<-sum((colSums(Shed.data.3$Onset)[A==3])>0)
Res[4]<-sum((colSums(Shed.data.3$Onset)[A==4])>0)
# Density plots by age (plots and saves as a pdf)
pdf('FigA.23.pdf',width = 10,height = 7)
par(mfrow=c(2,3))
# <1 year olds
plot(density(outcome1.age[,1]),type='l',xlim=c(10,60),lwd=2,main='<1 years',
     xlab='Total numbers infected in a single outbreak',
     ylim=c(0,max(density(outcome1.age[,1])$y,density(outcome1.2.age[,1])$y,
                  density(outcome1.3.age[,1])$y)))
lines(density(outcome1.2.age[,1]),col='red',lty=1,lwd=2)
lines(density(outcome1.3.age[,1]),col='blue',lty=1,lwd=2)
abline(v=Res[1],col='grey60',lty=2,lwd=2);
text(x=Res[1]-5,y=0.1,labels = c('Real data'),cex=0.7,col='grey60')
# 1-5 year olds
plot(density(outcome1.age[,2]),type='l',xlim=c(20,80),lwd=2,main='1-5 years',
     xlab='Total numbers infected in a single outbreak',
     ylim=c(0,max(density(outcome1.age[,2])$y,density(outcome1.2.age[,2])$y,
                  density(outcome1.3.age[,2])$y)))
lines(density(outcome1.2.age[,2]),col='red',lty=1,lwd=2)
lines(density(outcome1.3.age[,2]),col='blue',lty=1,lwd=2)
abline(v=Res[2],col='grey60',lty=2,lwd=2);
text(x=Res[2]-5,y=0.08,labels = c('Real data'),cex=0.7,col='grey60')
# 5-15 year olds
plot(density(outcome1.age[,3]),type='l',xlim=c(20,110),lwd=2,main='5-15 years',
     xlab='Total numbers infected in a single outbreak',
     ylim=c(0,max(density(outcome1.age[,3])$y,density(outcome1.2.age[,3])$y,
                  density(outcome1.3.age[,3])$y)))
lines(density(outcome1.2.age[,3]),col='red',lty=1,lwd=2)
lines(density(outcome1.3.age[,3]),col='blue',lty=1,lwd=2)
abline(v=Res[3],col='grey60',lty=2,lwd=2);
text(x=Res[3]-5,y=0.02,labels = c('Real data'),cex=0.7,col='grey60')
# >15 year olds
plot(density(outcome1.age[,4]),type='l',xlim=c(0,80),lwd=2,main='>15 years',
     xlab='Total numbers infected in a single outbreak',
     ylim=c(0,max(density(outcome1.age[,4])$y,density(outcome1.2.age[,4])$y,
                  density(outcome1.3.age[,4])$y)))
lines(density(outcome1.2.age[,4]),col='red',lty=1,lwd=2)
lines(density(outcome1.3.age[,4]),col='blue',lty=1,lwd=2)
abline(v=Res[4],col='grey60',lty=2,lwd=2);
text(x=Res[4]-5,y=0.02,labels = c('Real data'),cex=0.7,col='grey60')
# All the ages
plot(density(outcome1[,3]),type='l',xlim=c(70,270), main='All ages',lwd=2,
     xlab='Total numbers infected in a single outbreak',
     ylim=c(0,max(density(outcome1[,3])$y,density(outcome1.2[,3])$y,
                  density(outcome1.3[,3])$y)))
lines(density(outcome1.2[,3]),col='red',lty=1,lwd=2)
lines(density(outcome1.3[,3]),col='blue',lty=1,lwd=2)
abline(v=179,col='grey60',lty=2,lwd=2);
text(x=170,y=0.035,labels = c('Real data'),cex=0.7,col='grey60')

plot.new()
legend('topright',lty=c(1,1,1,2),col=c('black','red','blue','grey'),lwd=2,
       legend=c('Unaltered infectiousness','Reduced infectiousness of symptomatic',
                'Removed infectiousness of asymptomatic',
                'Observed data'),
       cex=1,bty='n')
dev.off()

#-------------------------------------------------------------------------------
# If only the symptomatic cases had been sampled, how different would the 
# parameter estimates be. We remove all the episodes that were not accompanied
# by an ARI. Symptomatic shedders can have days within the episode that were not
# accompanied by an ARI, these are still included in the model.
#-------------------------------------------------------------------------------
# Loading the model functions
source("Model_Funcs.R")
# Loading the fitting functions (from the fitR package) 
source("Sampler_fitR.R") 
# Loading the diagnostic functions
source("Model_Diag.R")
# Variables need for fiting that are not pathogen specific
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
# -An NxN matrix the identifies housemates.Used in FOI function to add up HH shedders
houses<-unique(Demo.data$hhid)
Cc<-matrix(0,N,N)
for (i in 1:length(houses)){
  # Their housemates
  pp<-which(Demo.data$hhid==houses[i])
  Cc[pp,pp]<-1
}
diag(Cc)<-0
# -Background community function
Comm.risk.a<<-Shed.data$Comm.risk;Comm.risk.b<<-Shed.data.2$Comm.risk
# -status variable
status<<-Shed.data$status

# ---- Algorithm specifics
# Standard deviation for the proposal distribution in the MCMC algorithm
Proposal.sd1<-c(Prev.hom=40,Prev.het=40,age.2=40,age.3=40,age.4=40,eta.A=40,eta.B=40,
                hh.size=40,HighAsym=40,LowSym=40,HighSym=40,epsilon.A=40, epsilon.B=40,
                eps.age2=40,eps.age3=40)
N.iter<-25000
# Initial parameter values for each chain
Paras1<-log(c(Prev.hom=0.658,Prev.het=0.658,age.2=1.113,age.3=0.364,age.4=0.195,
              eta.A=0.0175,eta.B=0.0175,hh.size=0.423,HighAsym=3,LowSym=0.5,HighSym=1.5,
              epsilon.A=0.00952,epsilon.B=0.00952,eps.age2=0.446,eps.age3=1.473))
Paras1.c2<-log(c(Prev.hom=0.158,Prev.het=0.658,age.2=1.113,age.3=0.364,age.4=0.195,
                 eta.A=0.035,eta.B=0.0175,hh.size=0.423,HighAsym=1,LowSym=1,HighSym=1,
                 epsilon.A=0.0476,epsilon.B=0.00952,eps.age2=0.446,eps.age3=1.473))

# ---- Infection data
# -RSV A
Shedding.a<<-Shed.data$Shedding;Onset.a<<-Shed.data$Onset
# 1) Remove all the infection data from people with no ARI
ari.a<-which(colSums(ARI.estimated)>0) # pple with ari RSV A
Shedding.a[,-ari.a]<-0;Onset.a[,-ari.a]<-0
# 2) For every ARI person with >1 onset, remove the non symptomatic episodes
x<-which(colSums(Onset.a[,ari.a])>1) # ARI RSV A pple with multiple onsets
for(i in 1:length(x)){
  # the episodes
  strt<-which(diff(Shedding.a[,ari.a[x[i]]])==1)+1
  end<-which(diff(Shedding.a[,ari.a[x[i]]])==-1)
  # Find if episode was symptomatic
  for(e in 1:length(strt)){
    ari.ep<-ARI.estimated[strt[e]:end[e],ari.a[x[i]]] 
    # if no ARI remove the data
    if (sum(ari.ep)==0){
      Shedding.a[strt[e]:end[e],ari.a[x[i]]] <-0
      Onset.a[strt[e],ari.a[x[i]]]<-0
    } 
  }
}
# Modify the rest of the data
Shedding.a.est.ct<-0*Onset.a;Shedding.a.est.ct[Shedding.a==1]<-Shed.data$Shedding.est.ct[Shedding.a==1]
Shedding.a.est.load<-0*Onset.a;Shedding.a.est.load[Shedding.a==1]<-Shed.data$Shedding.est.load[Shedding.a==1]
atRisk.a<-0*Onset.a;atRisk.a[Shedding.a==0 & status==1]<-1
postInf.a<-0*Shedding.a
for(p in 1:N){
  y<-which(diff(Shedding.a[,p])==-1)
  if(length(y)>0) postInf.a[((min(y)+1):D),p]<-1
}
rm(x)

# -RSV B
Shedding.b<<-Shed.data.2$Shedding; Onset.b<<-Shed.data.2$Onset
# 1) Remove all the infection data from people with no ARI
ari.b<-which(colSums(ARI.estimated.2)>0) # pple with ari RSV B
Shedding.b[,-ari.b]<-0;Onset.b[,-ari.b]<-0
# 2) For every ARI person with >1 onset, remove the non symptomatic episodes
y<-which(colSums(Onset.b[,ari.b])>1) # ARI RSV B pple with multiple onsets
for(i in 1:length(y)){
  # the episodes
  strt<-which(diff(Shedding.b[,ari.b[y[i]]])==1)+1
  end<-which(diff(Shedding.b[,ari.b[y[i]]])==-1)
  # Find if episode was symptomatic
  for(e in 1:length(strt)){
    ari.ep<-ARI.estimated.2[strt[e]:end[e],ari.b[y[i]]] 
    # if no ARI remove the data
    if (sum(ari.ep)==0){
      Shedding.b[strt[e]:end[e],ari.b[y[i]]] <-0
      Onset.b[strt[e],ari.b[y[i]]]<-0
    } 
  }
}
# Modify the rest of the data
Shedding.b.est.ct<-0*Onset.b;Shedding.b.est.ct[Shedding.b==1]<-Shed.data.2$Shedding.est.ct[Shedding.b==1]
Shedding.b.est.load<-0*Onset.b;Shedding.b.est.load[Shedding.b==1]<-Shed.data.2$Shedding.est.load[Shedding.b==1]
atRisk.b<-0*Onset.b;atRisk.b[Shedding.b==0 & status==1]<-1
postInf.b<-0*Shedding.b
for(p in 1:N){
  y<-which(diff(Shedding.b[,p])==-1)
  if(length(y)>0) postInf.b[((min(y)+1):D),p]<-1
}

# -Changing covariate(Bit), infection history with 2 groups, coding is 1=no,2=yes
B <<- matrix(1,D,N);B[postInf.a==1]<-2;B[postInf.b==1]<-3
B[(postInf.a==1 & postInf.b==1)]<-4 

# -Categorical variable for shedding
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

# Chain 1
T7.b.group<-system.time(Res.NoAsym.c1<-mcmcMH(target=IBM.posterior.group.sym,
                                                  init.theta=Paras1, proposal.sd=Proposal.sd1,n.iterations=N.iter,
                                                  adapt.size.start = 1000,adapt.shape.start = 500,adapt.size.cooling=0.999))

# Chain 2
T7.b.group.c2<-system.time(Res.NoAsym.c2<-mcmcMH(target=IBM.posterior.group.sym,
                                                     init.theta=Paras1.c2, proposal.sd=Proposal.sd1,n.iterations=N.iter,
                                                     adapt.size.start = 1000,adapt.shape.start = 500,adapt.size.cooling=0.999))

# ---- Diagnostics
prior.limits.sym<-data.frame(lower=c(Prev.hom=-10,Prev.het=-10,age.2=-10,
                                     age.3=-10,age.4=-10,eta.A=-20,eta.B=-20,hh.size=-10,
                                     HighAsym=-10,LowSym=-10,HighSym=-10,epsilon.A=-20,
                                     epsilon.B=-20,eps.age2=-10,eps.age3=-10),
                             upper=c(Prev.A=10,Prev.B=10,age.2=10,age.3=10,age.4=10,eta.A=0, eta.b=0,
                                     hh.size=10,HighAsym=10,LowSym=10,HighSym=10,
                                     epsilon.A=0,epsilon.B=0,eps.age2=10,eps.age3=10))
# RSV A and B (35k, 0, 35k, 0)
Res.diag.RSV.NoAsym<-IBM.all.diag(chain1=Res.NoAsym.c1,chain2=Res.NoAsym.c2,
                                      chain3=Res.NoAsym.c2,prior.limits=prior.limits.sym,
                                      desc='RSV A and B sym')

# Comparing to parameter distributions when using all the data
load("Para_distributions/Res.group.diag.RData")
theta.names<-names(Res.group.diag$Res)[-c(9,16)]
Par.dist1<-Res.group.diag$Res[,theta.names]

set.seed<-5
take<-sample(seq(dim(Par.dist1)[1]),500);
take2<-sample(seq(dim(Res.diag.RSV.sym.NoAsym$Res)[1]),500)
dens.plot<-function(j){
  x<-data.frame(All.data=Par.dist1[take,j],No.Asym=Res.diag.RSV.sym.NoAsym$Res[take2,j]);
  para.name<-names(Par.dist1)[j]
  data<- melt(x)
  ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
    labs(title=paste(para.name)) 
}
p1<-dens.plot(1);p2<-dens.plot(2);p3<-dens.plot(3);p4<-dens.plot(4);p5<-dens.plot(5)
p6<-dens.plot(6);p7<-dens.plot(7);p8<-dens.plot(8);p9<-dens.plot(9);p10<-dens.plot(10)
p11<-dens.plot(11);p12<-dens.plot(12);p13<-dens.plot(13);p14<-dens.plot(14)
p15<-dens.plot(15)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15)

Res.summary.All<-t(exp(apply(Par.dist1,2,Quntile.func)))

Res.summary.NoAsym<-t(exp(apply(Res.diag.RSV.sym.NoAsym$Res[names(Par.dist1)],2,Quntile.func)))

cbind(Res.summary.All,Res.summary.NoAsym)















