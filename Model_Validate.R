# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 13th September 2018
################################################################################
# Following the fitting procedure in 'Model_Run.R', these functions validate the 
# model in a two part process.
# 
# Part 1:Simulate epidemics using parameters sampled from the posterior and 
#        compare to the real data
# Part 2:Using simulated data sets, re-estimate parameters to check if parameters
#        used in simulation are contained in the re-estimated posterior
# 
################################################################################
# --- Funtions used in the validation 
# This function is for correlated sampling of p sets of parameters from a
# distribution Dist.
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
      V.par<-c(0,1,exp(theta[['HighAsym']]),exp(theta[['LowSym']]),exp(theta[['HighSym']]))
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
          # Reference group;low viral load and asymptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed<=high.load & prof$ari==0]<-2                      
          # high viral load and asymptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed>high.load & prof$ari==0]<-3                      
          # low viral load and symptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed<=high.load & prof$ari==1]<-4                      
          # high viral load and symptomatic
          Load.n.sym.a[d.sh,y[i]][prof$shed>high.load & prof$ari==1]<-5              
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
          # shedding profile                                                            <------------
          prof<-shed.profiler(Profiles.rsv.b,age.grp=x[i]) # Random assignment of profiles
          d.sh<-seq((t+Inf.b[t,y[i]]),(t+Inf.b[t,y[i]]+length(prof$shed)-1)) # Time window of infection
          d.sh<-d.sh[d.sh<=D];# modifying window incase it goes beyond the study period
          prof$shed<-prof$shed[seq(length(d.sh))];prof$ari<-prof$ari[seq(length(d.sh))]
          # --viral load and ARI
          I.b[d.sh,y[i]]<-prof$shed;ARI.b[d.sh,y[i]]<-prof$ari   
          # --infectiousness category
          # Reference group; asymptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed<=high.load & prof$ari==0]<-2                      
          # high viral load and asymptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed>high.load & prof$ari==0]<-3                      
          # low viral load and symptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed<=high.load & prof$ari==1]<-4                      
          # high viral load and symptomatic
          Load.n.sym.b[d.sh,y[i]][prof$shed>high.load & prof$ari==1]<-5              # <-------------------
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
# This function is for drawing box plots of the outcome measure
plot.box<-function(){
  # Total number of people infected
  boxplot(boxer(1),ylim=c(50,300),ylab='Total number infected',main='RSV A')
  abline(h=length(shed.a),col='red',lty=2)
  text(x=1,y =length(shed.a)+3,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(2),ylim=c(50,300),ylab='Total number infected',main='RSV B')
  abline(h=length(shed.b),col='red',lty=2)
  text(x=1,y =length(shed.b)+3,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(3),ylim=c(50,300),ylab='Total number infected',main='RSV')
  abline(h=length(shed.rsv),col='red',lty=2)
  text(x=1,y =length(shed.rsv)+3,labels = 'Data',col='red',cex=0.8)
  # Total number of households infected
  boxplot(boxer(4),ylim=c(20,47),ylab='Total number HH infected')
  abline(h=length(shed.a.hh),col='red',lty=2)
  text(x=1,y =length(shed.a.hh)+3,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(5),ylim=c(20,47),ylab='Total number HH infected')
  abline(h=length(shed.b.hh),col='red',lty=2)
  text(x=1,y =length(shed.b.hh)+3,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(6),ylim=c(20,47),ylab='Total number HH infected')
  abline(h=length(shed.rsv.hh),col='red',lty=2)
  text(x=1,y =length(shed.rsv.hh)+3,labels = 'Data',col='red',cex=0.8)
  # prop with >1 infections
  boxplot(boxer(7),ylim=c(0,0.6),ylab='Prop with >1 inf.')
  abline(h=prop.multi.a,col='red',lty=2)
  text(x=1,y =prop.multi.a+0.03,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(8),ylim=c(0,0.6),ylab='Prop with >1 inf.')
  abline(h=prop.multi.b,col='red',lty=2)
  text(x=1,y =prop.multi.b+0.03,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(9),ylim=c(0,0.6),ylab='Prop with >1 inf.')
  abline(h=prop.multi.rsv,col='red',lty=2)
  text(x=1,y =prop.multi.rsv+0.03,labels = 'Data',col='red',cex=0.8)
  
  # peak timings
  boxplot(boxer(10),ylim=c(40,180),ylab='Peak timing')
  abline(h=peak.a,col='red',lty=2)
  text(x=1,y =peak.a+3,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(11),ylim=c(40,180),ylab='Peak timing')
  abline(h=peak.b,col='red',lty=2)
  text(x=1,y =peak.b+3,labels = 'Data',col='red',cex=0.8)
  
  boxplot(boxer(12),ylim=c(40,180),ylab='Peak timing')
  abline(h=peak.rsv,col='red',lty=2)
  text(x=1,y =peak.rsv+3,labels = 'Data',col='red',cex=0.8)
  
}
# Given the simulation results, this function sorts them out and plots the 
# simulated curves against the real data 
curve.plot<-function(runs,Res.iter,description){
  # Real data
  real.a<-Shed.data$Shedding;real.b<-Shed.data.2$Shedding
  real.rsv<-real.a+real.b;real.rsv[real.rsv>0]<-1
  # RSV A
  RSV.a<-array(0,c(D,N,runs));a.s<-c(seq(1,(runs*5)-4,5))
  #  * status so that individuals who are shedding outside of the household are not
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
  
  # Everyone
  matplot(Shedding.iter.a,type='l',col='grey',lty=1,xlab='Days',main=paste('RSV A:',description),
          ylab='Numbers shedding',ylim=c(0,80));lines(rowSums(real.a),lwd=2)
  matplot(Shedding.iter.b,type='l',col='cyan',lty=1,xlab='Days',main=paste('RSV B:',description),
          ylab='Numbers shedding',ylim=c(0,80));lines(rowSums(real.b),lwd=2,col='blue')
  matplot(Shedding.iter.rsv,type='l',col='orange',lty=1,xlab='Days',main=paste('RSV:',description),
          ylab='Numbers shedding',ylim=c(0,80));lines(rowSums(real.rsv),lwd=2,col='green')
}

############################### PART 1: SIMULATION #############################
# We will used 5 different sets of parameters sampled dependently from the 
# posterior distribution and for each set we will run 200 simulations.

# Load the real data
source('Data_Prep.R')
# Given the real data, extract the shedding profiles
Profiles.rsv.a<-Profile.shed(Shed.data$Shedding,Shed.data$Shedding.est.load,ARI.estimated,'RSV A')
Profiles.rsv.b<-Profile.shed(Shed.data.2$Shedding,Shed.data.2$Shedding.est.load,ARI.estimated.2,'RSV B')
# Load the posterior distribution 
load("Para_distributions/Res.group.diag.RData")
theta.names<-names(Res.group.diag$Res)[-16]
Par.dist<-Res.group.diag$Res[,theta.names]
# Number of simulations to run per set
sims<-200
# Outcome measures from real data to compare to measures from simulated data
Shedding<-Shed.data$Shedding;Shedding[Shed.data.2$Shedding==1]<-1 
strt<-which(diff(as.matrix(Shedding))==1,TRUE);strt[,1]<-strt[,1]+1;
Onset<<-0*Shedding;Onset[strt]<-1
shed.a<-(which(colSums(Shed.data$Shedding)>0)) # pple who shed A (88)
shed.b<-(which(colSums(Shed.data.2$Shedding)>0)) # pple who shed B (113)
shed.rsv<-union(shed.a,shed.b)         # pple who shed any (from real data) (179)
shed.a.hh<-unique(Demo.data[['hhid']][shed.a])
shed.b.hh<-unique(Demo.data[['hhid']][shed.b])
shed.rsv.hh<-union(shed.a.hh,shed.b.hh)
peak.a<-min(which(rowSums(Shed.data$Shedding)==max(rowSums(Shed.data$Shedding))))
peak.b<-min(which(rowSums(Shed.data.2$Shedding)==max(rowSums(Shed.data.2$Shedding))))
peak.rsv<-min(which(rowSums(Shedding)==max(rowSums(Shedding))))
prop.multi.a<-length(which(colSums(Shed.data$Onset)>1))/length(shed.a)
prop.multi.b<-length(which(colSums(Shed.data.2$Onset)>1))/length(shed.b)
prop.multi.rsv<-length(which(colSums(Onset)>1))/length(shed.rsv)
# ---------------------------- Simulations -------------------------------------
# --- set1
# Parameters
paras.set1<-(samp.par(Par.dist,1));paras.set1<-as.numeric(paras.set1);
names(paras.set1)<-names(Par.dist)
# Simulations & outcome measures
set1<-Set.sim(paras.set1,sims);outcome1<-Outcome.measures(sims,set1)
# --- set2
# Parameters
paras.set2<-(samp.par(Par.dist,1));paras.set2<-as.numeric(paras.set2);
names(paras.set2)<-names(Par.dist)
# Simulations & outcome measures
set2<-Set.sim(paras.set2,sims);outcome2<-Outcome.measures(sims,set2)
# --- set3
# Parameters
paras.set3<-(samp.par(Par.dist,1));paras.set3<-as.numeric(paras.set3);
names(paras.set3)<-names(Par.dist)
# Simulations & outcome measures
set3<-Set.sim(paras.set3,sims);outcome3<-Outcome.measures(sims,set3)
# --- set4
# Parameters
paras.set4<-(samp.par(Par.dist,1));paras.set4<-as.numeric(paras.set4);
names(paras.set4)<-names(Par.dist)
# Simulations & outcome measures
set4<-Set.sim(paras.set4,sims);outcome4<-Outcome.measures(sims,set4)
# --- set5
# Parameters
paras.set5<-(samp.par(Par.dist,1));paras.set5<-as.numeric(paras.set5);
names(paras.set5)<-names(Par.dist)
# Simulations & outcome measures
set5<-Set.sim(paras.set5,sims);outcome5<-Outcome.measures(sims,set5)

Sim.res<-list(par1=paras.set1,runs1=set1, par2=paras.set2,runs2=set2,
              par3=paras.set3,runs3=set3, par4=paras.set4,runs4=set4,
              par5=paras.set5,runs5=set5)
# Saving simulation results
save(Sim.res, file= "Para_distributions/Sim.res.RData")

# ------------- Visualizing simulated epidemics and outcome measures -----------
# Box plot. Figure A.13 in the supplementary
boxer<-function(j)return(list(set1=outcome1[,j],set2=outcome2[,j],
                              set3=outcome3[,j],set4=outcome4[,j],
                              set5=outcome5[,j]))
par(mfrow=c(4,3),mar=c(4,4.5,2,4),cex.axis=0.7,cex.lab=1.1,cex.main=1.5,bg='white')
plot.box()
# Epidemic curve. Figure 4 in the main article
par(mfrow=c(3,3),mar=c(4,4.5,2,4),cex.axis=0.7,cex.lab=1.1,cex.main=1.5,bg='white')
curve.plot(sims,set1,'200 simsumations:set1');curve.plot(sims,set2,'200 simsumations:set2')
curve.plot(sims,set3,'200 simsumations:set3');curve.plot(sims,set4,'200 simsumations:set4')
curve.plot(sims,set5,'200 simsumations:set5')
# Parameters used in estimation relative to the posterior distribution. Figure A.12
par(mfrow=c(3,3),cex.main=1.1)
for(i in 1:15){
  hist(Par.dist[,i],col='grey',xlab='log value',main=names(Par.dist)[i],border = 'grey')
  abline(v=paras.set1[i],col='red',lty=2,lwd=2)
  abline(v=paras.set2[i],col='magenta',lty=2,lwd=2)
  abline(v=paras.set3[i],col='black',lty=2,lwd=2)
  abline(v=paras.set4[i],col='yellow',lty=2,lwd=2)
  abline(v=paras.set5[i],col='cyan',lty=2,lwd=2)
}
legend('topleft',lwd=c(7,rep(2,5)),col=c('grey','red','magenta','black','yellow','cyan'),
       lty=c(1,rep(2,5)),legend = c('Posterior','set1','set2','set3','set4','set5'),cex=0.6,
       bty='n')
# Single epidemics from each run
A.1<-Sim.res$runs1[[1 + 5*1]];A.1[A.1>0]<-1;B.1<-Sim.res$runs1[[2 + 5*1]];B.1[B.1>0]<-1
A.2<-Sim.res$runs2[[1 + 5*5]];A.2[A.2>0]<-1;B.2<-Sim.res$runs2[[2 + 5*5]];B.2[B.2>0]<-1
A.3<-Sim.res$runs3[[1 + 5*10]];A.3[A.3>0]<-1;B.3<-Sim.res$runs3[[2 + 5*10]];B.3[B.3>0]<-1
A.4<-Sim.res$runs4[[1 + 5*20]];A.4[A.4>0]<-1;B.4<-Sim.res$runs4[[2 + 5*20]];B.4[B.4>0]<-1
A.5<-Sim.res$runs5[[1 + 5*35]];A.5[A.5>0]<-1;B.5<-Sim.res$runs5[[2 + 5*35]];B.5[B.5>0]<-1
par(mfrow=c(2,3))
plot(rowSums(Shed.data$Shedding),type='l',col='black',lwd=2,main='Real data',ylim=c(0,35),
     xlab='Days',ylab='Numbers shedding')
lines(rowSums(Shed.data.2$Shedding),col='blue',lwd=2)
legend('topleft',lwd=2,col=c('black','blue'),legend=c('RSV A','RSV B'),bty='n',cex=0.8)
plot(rowSums(A.1),type='l',col='grey',lwd=2,main='set1',ylim=c(0,40),
     xlab='Days',ylab='Numbers shedding')
lines(rowSums(B.1),col='cyan',lwd=2)
legend('topleft',lwd=2,col=c('grey','cyan'),legend=c('RSV A','RSV B'),bty='n',cex=0.8)
plot(rowSums(A.2),type='l',col='grey',lwd=2,main='set2',ylim=c(0,35),
     xlab='Days',ylab='Numbers shedding')
lines(rowSums(B.2),col='cyan',lwd=2)
legend('topleft',lwd=2,col=c('grey','cyan'),legend=c('RSV A','RSV B'),bty='n',cex=0.8)
plot(rowSums(A.3),type='l',col='grey',lwd=2,main='set3',ylim=c(0,35),
     xlab='Days',ylab='Numbers shedding')
lines(rowSums(B.3),col='cyan',lwd=2)
legend('topleft',lwd=2,col=c('grey','cyan'),legend=c('RSV A','RSV B'),bty='n',cex=0.8)
plot(rowSums(A.4),type='l',col='grey',lwd=2,main='set4',ylim=c(0,35),
     xlab='Days',ylab='Numbers shedding')
lines(rowSums(B.4),col='cyan',lwd=2)
legend('topleft',lwd=2,col=c('grey','cyan'),legend=c('RSV A','RSV B'),bty='n',cex=0.8)
plot(rowSums(A.5),type='l',col='grey',lwd=2,main='set5',ylim=c(0,35),
     xlab='Days',ylab='Numbers shedding')
lines(rowSums(B.5),col='cyan',lwd=2)
legend('topleft',lwd=2,col=c('grey','cyan'),legend=c('RSV A','RSV B'),bty='n',cex=0.8)

############################### PART 2: RE-ESTIMATING ##########################
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
# Initial parameter values for each chain
Paras1<-log(c(Prev.hom=0.658,Prev.het=0.658,age.2=1.113,age.3=0.364,age.4=0.195,
              eta.A=0.0175,eta.B=0.0175,hh.size=0.423,HighAsym=3,LowSym=0.5,HighSym=1.5,
              epsilon.A=0.00952,epsilon.B=0.00952,eps.age2=0.446,eps.age3=1.473))
Paras1.c2<-log(c(Prev.hom=0.158,Prev.het=0.658,age.2=1.113,age.3=0.364,age.4=0.195,
                 eta.A=0.035,eta.B=0.0175,hh.size=0.423,HighAsym=1,LowSym=1,HighSym=1,
                 epsilon.A=0.0476,epsilon.B=0.00952,eps.age2=0.446,eps.age3=1.473))

# Length of each of the 2 MCMC chains
N.iter<-50000
# ------------------- Re-estimation using data generated by set1 ---------------
# Re-estimating using data generated by set1
Epi1<-Sim.res$runs1[seq(5)]
# Pairwise pathogen data needed for fitting
# Viral load and infection history matrices
Shedding.a.est.load<<-Epi1$RSV.A; Shedding.b.est.load<<-Epi1$RSV.B; B<<-Epi1$Inf.hist.mat
# Infectiousness category matrices
Load.n.sym.a<<-Epi1$Load.n.sym.a;Load.n.sym.b<<-Epi1$Load.n.sym.b
# binary shedding status
Shedding.a<<-0*B;Shedding.a[Shedding.a.est.load>0]<-1
Shedding.b<<-0*B;Shedding.b[Shedding.b.est.load>0]<-1
# Onset matrices
strt.a<-which(diff(Shedding.a)==1,TRUE);strt.a[,1]<-strt.a[,1]+1;
Onset.a<<-0*B;Onset.a[strt.a]<-1
strt.b<-which(diff(Shedding.b)==1,TRUE);strt.b[,1]<-strt.b[,1]+1;
Onset.b<<-0*B;Onset.b[strt.b]<-1
# At risk matrices
atRisk.a<<-0*B;atRisk.a[status==1 & Shedding.a==0]<-1
atRisk.b<<-0*B;atRisk.b[status==1 & Shedding.b==0]<-1
# Comparing the outbreak to be used for re-estimation to the real outbreak. Figure A.14
par(mfrow=c(1,2))
plot(rowSums(Shedding.a),type='l',col='black',lwd=2,main='RSV A',xlab='Days',
     ylab='Numbers shedding')
lines(rowSums(Shed.data$Shedding),col='red',lwd=2)
legend('topleft',col=c('black','red'),lty=c(1,1),lwd=c(2,2),legend=c('Epi1','Real'))
plot(rowSums(Shedding.b),type='l',col='black',lwd=2,main='RSV B',xlab='Days',
     ylab='Numbers shedding')
lines(rowSums(Shed.data.2$Shedding),col='red',lwd=2)
legend('topleft',col=c('black','red'),lty=c(1,1),lwd=c(2,2),legend=c('Epi1','Real'))
# ---- Running multiple chains
Time.group.c1<-system.time(Res.epi1.c1<-mcmcMH(target=IBM.posterior.group,
                                               init.theta=Paras1, proposal.sd=Proposal.sd1,n.iterations=N.iter,
                                               adapt.size.start = 1000,adapt.shape.start = 500,adapt.size.cooling=0.999))
Time.group.c2<-system.time(Res.epi1.c2<-mcmcMH(target=IBM.posterior.group,
                                               init.theta=Paras1.c2, proposal.sd=Proposal.sd1,n.iterations=N.iter,
                                               adapt.size.start = 1000,adapt.shape.start = 500,adapt.size.cooling=0.999))
# ---- Diagnostics
# Limits for the prior distributions
prior.limits.sym<-data.frame(lower=c(Prev.hom=-10,Prev.het=-10,age.2=-10,
                                     age.3=-10,age.4=-10,eta.A=-20,eta.B=-20,hh.size=-10,
                                     HighAsym=-10,LowSym=-10,HighSym=-10,epsilon.A=-20,
                                     epsilon.B=-20,eps.age2=-10,eps.age3=-10),
                             upper=c(Prev.A=10,Prev.B=10,age.2=10,age.3=10,age.4=10,eta.A=0, eta.b=0,
                                     hh.size=10,HighAsym=10,LowSym=10,HighSym=10,
                                     epsilon.A=0,epsilon.B=0,eps.age2=10,eps.age3=10))
# As only 2 chains were run, but the diagnostic function takes in 3, feed chain 3
# the same data as chain 1 but burn at the last iteration.
Res.epi1.diag<-IBM.all.diag(chain1=Res.epi1.c1,chain2=Res.epi1.c2,
                            chain3=Res.epi1.c1,prior.limits=prior.limits.sym,desc='RSV A and B sym')
Res.summary.epi1<-t(exp(apply(Res.epi1.diag$Res[names(Paras1)],2,Quntile.func)))

# ---- Visualisation of the re-estimation results
# Comparing parameter densities. Figure A.15
comp.dens.plot<-function(original,re.estimate){
  set.seed<-5
  take<-sample(seq(dim(original)[1]),500);take2<-sample(seq(dim(re.estimate$Res)[1]),500)
  dens.plot<-function(j){
    x<-data.frame(Original=original[take,j],Re.est=re.estimate$Res[take2,j]);para.name<-names(original)[j]
    data<- melt(x)
    ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
      labs(title=paste(para.name,':ESS of Re-est=',round(re.estimate$ESS[j],0)))+
      geom_vline(data=data,aes(xintercept = theta[j]),linetype = "dashed",col='red')    
  }
  p1<-dens.plot(1);p2<-dens.plot(2);p3<-dens.plot(3);p4<-dens.plot(4);p5<-dens.plot(5)
  p6<-dens.plot(6);p7<-dens.plot(7);p8<-dens.plot(8);p9<-dens.plot(9);p10<-dens.plot(10)
  p11<-dens.plot(11);p12<-dens.plot(12);p13<-dens.plot(13)
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13)
  
}
theta<-Sim.res$par1
comp.dens.plot(original=Par.dist,re.estimate = Res.epi1.diag)
# Infection in the largest household
Shedding<-Shedding.a;Shedding[Shedding.b==1]<-1;colnames(Shedding.a)<-Demo.data$sid
Dot.plot<-function(pple,pnt.size,lab.size){
  lng<-colSums(Shedding[,pple])
  Tm<-c();Pr<-c();Gr<-c()
  for(p in 1:length(pple)){
    # Time points with shedding
    t=which(Shedding[,pple[p]]==1);Tm<-c(Tm,t)
    Pr<-c(Pr,rep(p,length(t)))
    # A
    t.a<-which(Shedding.a[,pple[p]]==1)
    g<-rep('B',length(t));g[which(t %in% t.a)]<-'A'
    # B
    t.b<-which(Shedding.b[,pple[p]]==1)
    g[which(t %in% intersect(t.b,t.a))]<-'AB'
    Gr<-c(Gr,g)
  }
  
  Inf.plot<-data.frame(time=Tm,Person=Pr,RSV=Gr)
  p<-ggplot(Inf.plot, aes(time,Person, color=factor(RSV)))+geom_point(size=pnt.size)
  p + scale_y_continuous(breaks=seq(1,length(pple)),labels=colnames(Shedding.a)[pple]) +
    labs(title='Shedding patterns for Infected individuals') +
    theme(axis.text.y = element_text(size=lab.size,colour ='cyan'))+
    scale_colour_manual(values=c("#999999", "#56B4E9", "#CC79A7"))
  
}
hh.5<-which(Demo.data$hhid==5)
Dot.plot(hh.5,2,10)

