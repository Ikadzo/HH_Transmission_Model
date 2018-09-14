# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 13th September 2018
################################################################################
# This function modifies the raw data into matrices that are then used to estimate
# parameters of interest. It is called by the main script 'Model_Run.R'.
# There are two types of raw data that are modified in this script:
# 1) Household (HH) incidence data. This has 11 variables as described in the variable 
#    dictionary. It is modified to give the following data:
#    a) Binary matrix of estimated shedding durations : "Shedding"
#    b) Binary matrix of onset dates based on estimated shedding durations: "Onset"
#    c) Numeric matrix of imputed Ct values: "Shedding.est.ct"
#    d) Numeric matrix of imputed viral loads: "Shedding.est.load"
#    e) Binary matrix of precence or absence in the household: "status"
#    f) Binary matrix indicating if an individual is at risk of within HH exposure:
#       "atRisk"
#    g) Binary matrix indicating if an individual has recovered from an infection 
#       within the epidemic: "postInf"
#    h) Numeric vector of time varying community density curves based on primary 
#       household incidence: "Comm.risk"
#    g) Binary matrix of imputed ARI status within a shedding episode: "ARI.estimated"
#    
#    These matrices are created for RSV A and B and stored as objects in a list 
#    for each RSV group
#    
# 2) Hospital incidence. This has 3 variables also described in the variable 
#    dictionary.It is modified to give density curves for RSV A and B that can 
#    be used as an alternative to those generated from primary household incidence. 
#    These curves are contained in a dataframe named "Risk.hosp".
################################################################################
################################################################################
library('ggplot2')
# ------- These are the functions that are used to modify the raw data ---------

# This function extracts matrices of ct values on given dates 
Mat.gen<-function(path.data){
  Nms<-list(row=seq(D),cols=Demo.data$sid)
  load.shedding<-data.frame(matrix(-3,D,N),row.names = Sample.dates);colnames(load.shedding)<-Demo.data$sid
  for( i in 1:N){
    # days of  samples
    y<-which(Data$sid %in% Demo.data$sid[i])# ct data for person i
    yy<-which(Sample.dates %in% Data$sampledate[y]) # dates of i's ct data on new matrix
    if (length(y)!=length(yy)){
      y.2<-y[!duplicated(Data$sampledate[y])] # removing repeated rows of no ct
      load.shedding[yy,i]<-path.data[y.2]# fill given rows with Ct value
    }else{
      load.shedding[yy,i]<-path.data[y]# fill given rows with Ct value
    }
  }
  load.shedding[which(is.na(load.shedding)==TRUE,TRUE)]<--3
  return(load.shedding)
}
# This function generates the data needed for fitting based on matrices of Ct values
# In addition, it also plots the background community function against the data 
# that waas used to estimate it (primary household incidence) and the epidemic 
# curve for the given data based on imputed shedding durations.
Data.gen<-function(virus.shedding,rsv.group){
  # Funtion for counting episodes and estimating shedding durations
  Episode.finder<-function(person.data){
    shedding.p<-rep(0,length(person.data));onset.p<-rep(0,length(person.data))
    censored.p<-rep(0,length(person.data))
    s<-which(person.data>=0) # samples
    # Are there samples
    if (length(s)>0){
      p<-which((person.data>0) & (person.data<=Pos.cut.off)) # positives
      n<-setdiff(s,p) # negatives
      # Are there positive samples
      if (length(p)>0){
        d<-diff(p);e<-which(d>14)
        E<-c(1,e+1) #index of first sample for every episode 
        L<-c((E[-1]-1),length(p)) #index of last sample for every episode 
        for (i in 1:length(E)){
          # First +ve sample and last negative before 1st +ve
          f.p<-p[E[i]]
          yy<-which(n<f.p);if(length(yy)>0)l.n<-n[max(yy)] else l.n<- 0
          # last +ve sample and 1st -ve after episode  
          l.p<-p[L[i]]
          xx<-which(n>l.p);if(length(xx)>0)f.n<-n[min(xx)] else f.n<-length(person.data)+1
          # Start of episode
          if(((f.p-l.n)>7)) {
            start<-f.p-half.mean.SI;censored.p[floor(start)]<-1
          }else{
            start<-l.n+0.5*(f.p-l.n)
          }
          # End of episode
          if(((f.n-l.p)>7)){
            end<-l.p+half.mean.SI;censored.p[floor(end)]<-2
          }else{
            end<-l.p+0.5*(f.n-l.p)
          }
          #shedding.p[floor(start):floor(end)]<-1;onset.p[floor(start)]<-1
          dur<-ceiling(end-start)
          # incase of short interval between samples i.e. a day 
          if((floor(start)+dur-1)<l.p) lastday<-l.p else lastday<-(floor(start)+dur-1)
          shedding.p[floor(start):lastday]<-1;onset.p[floor(start)]<-1
          
        }
      }
    }
    Res<-data.frame(shedding=shedding.p,onset=onset.p,censored=censored.p)
  }
  # Ct values Interpolating function
  Ct.Interpolate<-function(i){
    a.est<-0*Onset.a[,i]
    # Did they have an onset?
    O.a<-which(Onset.a[,i]>0)
    if(length(O.a)>0){
      epis<-length(O.a) # No. of episodes
      for(e in 1:epis){
        X<-which(Shedding.a[,i]==0);end<-X[min(which(X>O.a[e]))]-1 # end of episode
        epi<-c(O.a[e]:end) # shedding days
        # First and last days of epi,Ct=boderline
        a.est[c(O.a[e],end)]<-Pos.cut.off
        # days with +ve samples in episode + 1st and last days
        load<-which((virus.shedding[epi,i]>0) & (virus.shedding[epi,i]<=Pos.cut.off))
        a.est[epi[load]]<-virus.shedding[epi[load],i]
        load<-c(1,load,length(epi))
        # Interpolating(assuming 1st and last days of epi have data=boderline +ve)
        gap<-diff(load)
        for(g in 1:length(gap)){ # for each gap
          if(gap[g]>1){
            for(gg in 1:(gap[g]-1)) # for each day in the gap
              a.est[(epi[load[g]]+gg)]=a.est[epi[load[g]]] +
                gg*((a.est[epi[load[g+1]]]-a.est[epi[load[g]]])/gap[g])
          }
        }
      }
    }
    return(a.est)
  }
  # Function for converting Ct values to viral densities/load
  V.load<-function(Ct){
    (Ct-42.9)/-3.308  
  }
  # Function for interpolatinf the status matrix. It need info presence during NFS collection
  Stat.fill<-function(p){
    Status<-rep(-1,D)
    # Data from the person
    p.d<-which(Data$sid==Demo.data$sid[p])
    # Dates of data 
    days<-as.numeric(Data$sampledate[p.d]-min(Data$sampledate)+1)
    stat<-Data$awaynfs[p.d]
    # filling in the data gaps
    for(g in 1:(length(days)-1)){
      f.p<-days[g];l.p<-days[g+1]
      if((l.p-f.p)>1){ # if there is a discontinuity fill in gaps
        if(stat[g]==stat[(g+1)]){ 
          Status[(f.p:l.p)]<-stat[g]
        }else{
          mid.point<-f.p+round((l.p-f.p)/2,0)
          Status[(f.p:mid.point)]<-stat[g]
          Status[((mid.point+1):l.p)]<-stat[(g+1)]
        }
      }else{ # if there is no discontinuity
        Status[c(f.p,l.p)]<-stat[c(g,(g+1))]
      }
    }
    return(Status)
  }
  # This function generates the background community density curves from primary 
  # household incidence for a given pathogen
  CommRisk.func<<-function(virus.shedding,rsv.group){
    # Scaling between 0 and 1
    scaler<-function(x) return(((x-min(x))/diff(range(x))))
    # infection data dates
    sampledates<-range(as.Date(rownames(virus.shedding)))
    #_______________________________________________________________________________
    #_________________________1. Gaussian fitting___________________________________
    # a. infection data hhid
    houses<-unique(Demo.data$hhid);weeks<-seq(sampledates[1],sampledates[2],7)
    days<-seq(sampledates[1],sampledates[2],1)
    HH.on.a=matrix(0,D,length(houses));colnames(HH.on.a)<-houses
    # b. For every household, get the date of primary onset
    for (j in 1:length(houses)){
      pple<-which(Demo.data$hhid==houses[j])
      on.days.a<-which(rowSums(Onset.a[,pple])>0)
      if(length(on.days.a)>0) HH.on.a[min(on.days.a),j]<-1 
    } 
    HH.inc<-data.frame(inc.a=rowSums(HH.on.a))
    # c. Average weekly incidence
    HH.inc.wk<-data.frame(inc.a=rep(0,length(weeks)))
    rownames(HH.inc.wk)<-weeks
    for(i in 1:length(weeks)){
      if(i==length(weeks)){
        # Weekly numbers
        y=which(days>=weeks[i])
        HH.inc.wk[i,]<-sum(HH.inc[y,])/(as.numeric(diff(c(weeks[i],sampledates[length(sampledates)])))+1)
      }else{
        # Weekly averages
        y=which((days>=weeks[i]) & (days<weeks[(i+1)]))
        HH.inc.wk[i,]<-sum(HH.inc[y,])/as.numeric(diff(weeks[c(i,i+1)])) 
      }
    }  
    
    # d. Fitting the curves
    # Number of weeks
    model_time<-seq(length(weeks))
    # Incidence function
    incidence_fn <- function ( hazard_fn, coeffs ) {
      #returns the incidence by integrating the hazard and calculating the cumulative
      #distribution function
      integ_haz <- c() #the integrated hazard from 0 to t
      cdf <- c(0) #the cumulative distribution function
      incid <- 0 #the model incidence in each week
      for (t in 1:length(model_time)){
        xx <- integrate(hazard_fn,lower=0,upper=model_time[t],coeffs)
        integ_haz[t] <- xx$value
        cdf[t] <- length(houses) * (1 - exp ( -integ_haz[t]))
        if (t==1){
          incid[t] <- max(1e-200,cdf[t])    
        } else {
          incid[t] <- max(1e-200,(cdf[t]-cdf[t-1])) # added the max bit so incid !<0
        }
      }
      return (incid)
    }
    # -log likelihood function(using total RSV data)
    calc_neg_ll <- function (hazard_coeffs, hazard_fn, rsv.data ) {
      # Model
      model_incid<-incidence_fn( hazard_fn, hazard_coeffs )
      # Log likelihood
      Res<- sum((rsv.data*log(model_incid))-model_incid)
      return(-Res)
    }
    # scaled Gaussin hazard function
    hazard_gauss <- function ( tt, cffs ) {
      a1 <- cffs[1] ; b1 <- cffs[2] ; c1 <- cffs[3]
      haz<-a1*exp(-((tt-b1)/c1)^2)
      return (haz)
    }
    # Fitting 
    guess<-c(0.5,12,2)
    model.gauss.a<- optim(guess,calc_neg_ll, hazard_fn = hazard_gauss, rsv.data=HH.inc.wk$inc.a )
    model_incid.a<-incidence_fn( hazard_gauss, model.gauss.a$par)
    # e. Stretching the model time curves
    HH.inc.dy.mod<-data.frame(inc.a=rep(0,length(days)))
    for(i in 1:length(weeks)){
      if(i==length(weeks)){
        # Daily numbers
        HH.inc.dy.mod$inc.a[which((days>=weeks[i])&(days<=sampledates[2]))]<-model_incid.a[i]
      }else{
        # Daily numbers
        HH.inc.dy.mod$inc.a[which((days>=weeks[i])&(days<weeks[i+1]))]<-model_incid.a[i]
      }
    }
    #_______________________________________________________________________________
    # Visualizing output of process
    #_______________________________________________________________________________
    if(rsv.group=='A') plot.col<-'grey';if(rsv.group=='B') plot.col<-'magenta'
    plot(weeks,HH.inc.wk$inc.a,type='b',pch=20,col='black',xlab='Time in weeks',
         ylab='Weekly average number of primary outbeaks',
         main=paste('RSV',rsv.group,': Incidence of primary outbreaks at household level'),lwd=2)
    lines(weeks,model_incid.a,lty=1,col=plot.col,lwd=3)
    legend('topright',col=c('black',plot.col),lty=c(NA,1),pch=c(20,NA),bty='n',
           legend=c('Data','Fitted curve'),cex=1.2,lwd=c(NA,3))
    # Data for model fitting
    Risk.virus<-scaler(HH.inc.dy.mod$inc.a)
    z<-which(Risk.virus==0);Risk.virus[z]<-0.0000001
    
    return(Risk.virus)
  }
  
  D<-dim(virus.shedding)[1];N<-dim(virus.shedding)[2]
  #--------------------1. Extracting shedding profiles ---------------------------
  # a) Mean sampling interval for samples within an episode
  # -------------------------------------------------------
  virus<-0*virus.shedding
  a.pos<-which(((virus.shedding>0) & (virus.shedding<=Pos.cut.off)),TRUE)
  virus[a.pos]<-1;all.means<-matrix(-1,N,5)
  for(j in 1:N){
    p<-which(virus[,j]==1) # +ve samples
    if (length(p)>0){
      d<-diff(p);e<-which(d>14) # Identifying episodes as samples that are >14 days apart
      E<-c(1,e+1) #index of first sample for every episode 
      L<-c((E[-1]-1),length(p)) #index of last sample for every episode 
      for (i in 1:length(E)){
        days<-(p[E[i]]:p[L[i]])# days in episode
        s<-which(virus.shedding[days,j]>=0) # days with data in the episode
        all.means[j,i]<-mean(diff(days[s]))
      }
    }
  }
  rm(virus)
  # Total Number of episodes
  test<-as.numeric(all.means);test2<-test[-which(test==-1)]
  # Remove single sample episodes from calculation of overall mean
  x<-which(is.na(test2)==TRUE);half.mean.SI<-0.5*mean(test2[-x])
  
  message('Using a cut-off of ',Pos.cut.off,',',length(test2),' episodes were 
          observed, of which ',length(x),' were excluded from calculation of the mean sampling
          interval on account of having a single sample')
  # b) Episodes
  # -----------
  Virus.A<-apply(virus.shedding,2,Episode.finder)
  Virus.A<-as.data.frame(Virus.A)
  # 0/1 matrices of shedding and onset days
  Shedding.a<-Virus.A[,seq(1,length(Virus.A),3)];colnames(Shedding.a)<-colnames(virus.shedding)
  Onset.a<-Virus.A[,seq(2,length(Virus.A),3)];colnames(Onset.a)<-colnames(virus.shedding)
  #----------------2. Linear Interplolation(of Ct values)-----------------------
  Shedding.a.est.ct<-apply(matrix(seq(N)),1,Ct.Interpolate)
  #----------------3. converting interpolated Cts to V.loads---------------------
  Shedding.a.est.load<-0*Shedding.a.est.ct
  x<-which(Shedding.a.est.ct>0,TRUE)
  Shedding.a.est.load[x]<-V.load(Shedding.a.est.ct[x]);rm(x)
  #----------------4. Additional Data needing for model fitting-------------------
  # 2. Status matrix
  # Interpolating status like shedding duration calculation
  status<-apply(matrix(seq(N)),1,Stat.fill)
  # if shedding then they were present
  status[Shedding.a==1]<-1 
  # Filling in the rest as away
  x<-which(status==-1,TRUE); status[x]<-0
  
  # 3. At risk
  atRisk.a<-0*Shedding.a;atRisk.a[((status==1) & (Shedding.a==0))]<-TRUE
  
  # 4. PostInf    
  postInf.a<-0*Shedding.a
  for(p in 1:N){
    y<-which(diff(Shedding.a[,p])==-1)
    if(length(y)>0) postInf.a[((min(y)+1):D),p]<-1
  }
  # 5. Community densty function 
  par(mfrow=c(1,1),oma = c(0, 0, 1, 0),mar=c(5,5,2,4)+.1,cex=1.2,cex.lab=1.2,
      cex.main=1.1,bg='white')
  plot(rowSums(Shedding.a),type='l',xlab='Days in the data',main=Pathogen.name,
       ylab='Total number of people shedding',lwd=2)
  
  Comm.risk<-CommRisk.func(virus.shedding,rsv.group)
  
  Results<-list(Shedding=Shedding.a,Onset=Onset.a,Shedding.est.ct=Shedding.a.est.ct,
                Shedding.est.load=Shedding.a.est.load,status=status,atRisk=atRisk.a,
                postInf=postInf.a,Comm.risk=Comm.risk)
  return(Results)
}
# This function estimates duration of ARI episodes within shedding episodes
ARI.est<-function(Shedding,Onset,ARI){
  # For every person, this fuction finds the shedding episode and estimates ARI
  # episode duration within it.
  epi.finder<-function(j){
    ari.p<-rep(0,D)
    Strt<-which(Onset[,j]==1) # start days of each episode
    Fin<-which(diff(Shedding[,j])==-1) # end days of each episode
    # If there are shedding episodes
    if(length(Strt)>0){
      for (i in 1:length(Strt)){
        x<-Strt[i]:Fin[i] # Days of shedding in the episode
        # per episode
        person.data<-ARI[x,j]
        symptom.p<-rep(0,length(person.data));onset.p<-rep(0,length(person.data))
        censored.p<-rep(0,length(person.data))
        s<-which(person.data>=0) # samples
        if (length(s)>0){
          p<-which(person.data==1) # positives
          n<-setdiff(s,p) # negatives
          # Are there positive samples
          if (length(p)>0){
            d<-diff(p);e<-which(d>14)
            E<-1 #index of first sample for ARI episode
            L<-length(p) #index of last sample for ARI episode
            # First +ve sample and last negative before 1st +ve
            f.p<-p[E]
            yy<-which(n<f.p);if(length(yy)>0)l.n<-n[max(yy)] else l.n<- 0
            # last +ve sample and 1st -ve after episode  
            l.p<-p[L]
            xx<-which(n>l.p);if(length(xx)>0)f.n<-n[min(xx)] else f.n<-length(person.data)+1
            # Start of episode
            if(((f.p-l.n)>7)) {
              start<-f.p-half.mean.SI;censored.p[floor(start)]<-1
            }else{
              if(f.p==1)start<-1 else start<-l.n+0.5*(f.p-l.n)
            }
            # End of episode
            if(((f.n-l.p)>7)){
              end<-l.p+half.mean.SI;censored.p[floor(end)]<-2
            }else{
              end<-l.p+0.5*(f.n-l.p)
            }
            #shedding.p[floor(start):floor(end)]<-1;onset.p[floor(start)]<-1
            dur<-ceiling(end-start)
            # incase of short interval between samples i.e. a day 
            if((floor(start)+dur-1)<l.p) lastday<-l.p else lastday<-(floor(start)+dur-1)
            symptom.p[floor(start):lastday]<-1;onset.p[floor(start)]<-1
            ari.p[x]<-symptom.p
          }
        }
      }
      
    }
    return(ari.p)
  }
  # -------------------------------------------------------
  # a) Mean sampling interval for samples within an episode
  # -------------------------------------------------------
  all.means<-matrix(-1,N,10)
  for(j in 1:N){
    p<-which(ARI[,j]==1) # +ve samples
    if (length(p)>0){
      d<-diff(p);e<-which(d>14) # Identifying episodes as samples that are >14 days apart
      E<-c(1,e+1) #index of first sample for every episode 
      L<-c((E[-1]-1),length(p)) #index of last sample for every episode 
      for (i in 1:length(E)){
        days<-(p[E[i]]:p[L[i]])# days in episode
        s<-which(ARI[days,j]>=0) # days with data in the episode
        all.means[j,i]<-mean(diff(days[s]))
      }
    }
  }
  # Total Number of episodes
  test<-as.numeric(all.means);test2<-test[-which(test==-1)]
  # Remove single sample episodes from calculation of overall mean
  x<-which(is.na(test2)==TRUE);half.mean.SI<-0.5*mean(test2[-x])
  # b) Episodes
  # -----------
  ARI.est<-sapply(seq(N),epi.finder)
  return(ARI.est)
}
# This function generates data for fitting RSV as a single pathogen, i.e. not 
# making a distinction between RSV A and B
RSV.gen<-function(){
  CommRisk.func<<-function(virus.shedding){
    # Scaling between 0 and 1
    scaler<-function(x) return(((x-min(x))/diff(range(x))))
    # infection data dates
    sampledates<-range(as.Date(rownames(virus.shedding)))
    #_______________________________________________________________________________
    #_________________________1. Gaussian fitting___________________________________
    # a. infection data hhid
    houses<-unique(Demo.data$hhid);weeks<-seq(sampledates[1],sampledates[2],7)
    days<-seq(sampledates[1],sampledates[2],1)
    HH.on.a=matrix(0,D,length(houses));colnames(HH.on.a)<-houses
    # b. For every household, get the date of primary onset
    for (j in 1:length(houses)){
      pple<-which(Demo.data$hhid==houses[j])
      on.days.a<-which(rowSums(Onset.a[,pple])>0)
      if(length(on.days.a)>0) HH.on.a[min(on.days.a),j]<-1 
    } 
    HH.inc<-data.frame(inc.a=rowSums(HH.on.a))
    # c. Average weekly incidence
    HH.inc.wk<-data.frame(inc.a=rep(0,length(weeks)))
    rownames(HH.inc.wk)<-weeks
    for(i in 1:length(weeks)){
      if(i==length(weeks)){
        # Weekly numbers
        y=which(days>=weeks[i])
        HH.inc.wk[i,]<-sum(HH.inc[y,])/(as.numeric(diff(c(weeks[i],sampledates[length(sampledates)])))+1)
      }else{
        # Weekly averages
        y=which((days>=weeks[i]) & (days<weeks[(i+1)]))
        HH.inc.wk[i,]<-sum(HH.inc[y,])/as.numeric(diff(weeks[c(i,i+1)])) 
      }
    }  
    
    # d. Fitting the curves
    # Number of weeks
    model_time<-seq(length(weeks))
    # Incidence function
    incidence_fn <- function ( hazard_fn, coeffs ) {
      #returns the incidence by integrating the hazard and calculating the cumulative
      #distribution function
      integ_haz <- c() #the integrated hazard from 0 to t
      cdf <- c(0) #the cumulative distribution function
      incid <- 0 #the model incidence in each week
      for (t in 1:length(model_time)){
        xx <- integrate(hazard_fn,lower=0,upper=model_time[t],coeffs)
        integ_haz[t] <- xx$value
        cdf[t] <- length(houses) * (1 - exp ( -integ_haz[t]))
        if (t==1){
          incid[t] <- max(1e-200,cdf[t])    
        } else {
          incid[t] <- max(1e-200,(cdf[t]-cdf[t-1])) # added the max bit so incid !<0
        }
      }
      return (incid)
    }
    # -log likelihood function(using total RSV data)
    calc_neg_ll <- function (hazard_coeffs, hazard_fn, rsv.data ) {
      # Model
      model_incid<-incidence_fn( hazard_fn, hazard_coeffs )
      # Log likelihood
      Res<- sum((rsv.data*log(model_incid))-model_incid)
      return(-Res)
    }
    # scaled Gaussin hazard function
    hazard_gauss <- function ( tt, cffs ) {
      a1 <- cffs[1] ; b1 <- cffs[2] ; c1 <- cffs[3]
      haz<-a1*exp(-((tt-b1)/c1)^2)
      return (haz)
    }
    # Fitting 
    guess<-c(0.5,12,2)
    model.gauss.a<- optim(guess,calc_neg_ll, hazard_fn = hazard_gauss, rsv.data=HH.inc.wk$inc.a )
    model_incid.a<-incidence_fn( hazard_gauss, model.gauss.a$par)
    # e. Stretching the model time curves
    HH.inc.dy.mod<-data.frame(inc.a=rep(0,length(days)))
    for(i in 1:length(weeks)){
      if(i==length(weeks)){
        # Daily numbers
        HH.inc.dy.mod$inc.a[which((days>=weeks[i])&(days<=sampledates[2]))]<-model_incid.a[i]
      }else{
        # Daily numbers
        HH.inc.dy.mod$inc.a[which((days>=weeks[i])&(days<weeks[i+1]))]<-model_incid.a[i]
      }
    }
    #_______________________________________________________________________________
    # Visualizing output of process
    #_______________________________________________________________________________
    par(mfrow=c(1,1))
    plot(weeks,HH.inc.wk$inc.a,type='b',pch=20,col='black',xlab='Time in weeks',
         ylab='Weekly average number of primary outbeaks',lwd=2,
         main='RSV Incidence of primary outbreaks at household level')
    lines(weeks,model_incid.a,lty=1,col='red',lwd=2)
    legend('topright',col=c('black','red'),lty=c(NA,1),pch=c(20,NA),legend=c('Data','Fitted curve'),
           cex=1,bty='n')

    # Data for model fitting
    Risk.virus<-scaler(HH.inc.dy.mod$inc.a)
    z<-which(Risk.virus==0);Risk.virus[z]<-0.0000001
    
    return(Risk.virus)
  }
  # Viral load matrix
  Shedding.est.load<-Shed.data$Shedding.est.load+Shed.data.2$Shedding.est.load
  # Binary shedding matrix
  Shedding<-0*Shedding.est.load;Shedding[Shedding.est.load>0]<-1
  # Infection history matrix
  B<-0*Shedding
  for(t in 3:D){
    b1<-which(colSums(Shedding[seq(t-1),])>0 & (Shedding[t,]==0));B[seq(t,D),b1]<-2
  }
  # Onset
  strt<-which(diff(Shedding)==1,TRUE);strt[,1]<-strt[,1]+1;
  Onset<-0*Shedding;Onset[strt]<-1
  # Presence in HH
  status<-Shed.data$status
  # At risk
  atRisk<-0*Shedding;atRisk[Shedding==0 & status==1]<-1
  # Community risk
  Onset.a<-Onset;Comm.risk<-CommRisk.func(Mat.gen(Data$rsva))
  
  Res=list(Shedding=Shedding,Onset=Onset,Shedding.est.load=Shedding.est.load,
           status=status,atRisk=atRisk,B=B,Comm.risk=Comm.risk)
  return(Res)
}
# This function generates the hospital incidence based background community curves
Comm.risk.hosp<-function(){
  # Scaling between 0 and 1
  scaler<-function(x) return(((x-min(x))/diff(range(x))))
  #_______________________________________________________________________________
  #_________________________1. Hospital incidence_________________________________
  sampledates<-range(Sample.dates)
  # Dates
  Hosp.data$date_admission<-as.Date(as.character(Hosp.data$date_admission),"%d/%m/%y")
  x<-which((Hosp.data$date_admission>=sampledates[1]) & (Hosp.data$date_admission<=sampledates[2]))
  # Data to keep
  #Hosp.inc<-Hosp.data[x,c(3,5,6,7)];rm(Hosp.data)
  Hosp.inc<-Hosp.data[x,];rm(Hosp.data)
  # Modifying factor variables
  Hosp.inc$rsva_result<-as.numeric(Hosp.inc$rsva_result)
  Hosp.inc$rsva_result[Hosp.inc$rsva_result!=2]<-0
  Hosp.inc$rsva_result[Hosp.inc$rsva_result==2]<-1
  Hosp.inc$rsvb_result<-as.numeric(Hosp.inc$rsvb_result)
  Hosp.inc$rsvb_result[Hosp.inc$rsvb_result!=2]<-0
  Hosp.inc$rsvb_result[Hosp.inc$rsvb_result==2]<-1
  # Getting weekly averages and filling in subsequent daily nunbers
  weeks<-seq(sampledates[1],sampledates[2],7)
  Hosp.inc.wk<-data.frame(inc.a=rep(0,length(weeks)),inc.b=rep(0,length(weeks)))
  rownames(Hosp.inc.wk)<-weeks
  
  days<-seq(sampledates[1],sampledates[2],1)
  Hosp.inc.dy<-data.frame(inc.a=rep(0,length(days)),inc.b=rep(0,length(days)))
  rownames(Hosp.inc.dy)<-days
  for(i in 1:length(weeks)){
    if(i==length(weeks)){
      # Weekly numbers
      y=which(Hosp.inc$date_admission>=weeks[i])
      Hosp.inc.wk[i,]<-colSums(Hosp.inc[y,c("rsva_result","rsvb_result")])/(as.numeric(diff(c(weeks[i],sampledates[length(sampledates)])))+1)
      # Daily numbers
      Hosp.inc.dy[which((days>=weeks[i])&(days<=sampledates[2])),]<-Hosp.inc.wk[i,]
    }else{
      # Weekly averages
      y=which((Hosp.inc$date_admission>=weeks[i]) & (Hosp.inc$date_admission<weeks[(i+1)]))
      Hosp.inc.wk[i,]<-colSums(Hosp.inc[y,c("rsva_result","rsvb_result")])/as.numeric(diff(weeks[c(i,i+1)])) 
      # Daily numbers
      Hosp.inc.dy[which((days>=weeks[i])&(days<weeks[i+1])),]<-Hosp.inc.wk[i,]
    }
  }  
  Hosp.inc.dy$inc<-Hosp.inc.dy$inc.a+Hosp.inc.dy$inc.b
  Hosp.inc.wk$inc<-Hosp.inc.wk$inc.a+Hosp.inc.wk$inc.b
  # Total cases
  sum(Hosp.inc$rsva_result)+sum(Hosp.inc$rsvb_result)   
  Risk.rsv<-data.frame(Hosp.a=scaler(Hosp.inc.dy$inc.a),Hosp.b=scaler(Hosp.inc.dy$inc.b),Hosp=scaler(Hosp.inc.dy$inc))
  
  z<-which(Risk.rsv$Hosp.a==0);Risk.rsv$Hosp.a[z]<-0.0000001
  z<-which(Risk.rsv$Hosp.b==0);Risk.rsv$Hosp.b[z]<-0.0000001
  z<-which(Risk.rsv$Hosp==0);Risk.rsv$Hosp[z]<-0.0000001
  return(Risk.rsv)
  
}
# This plots the shedding episodes
Dot.plot<-function(Shedding.a, Shedding.b,pple,pnt.size,lab.size){
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
  lab.col<-colnames(Shedding.a)[pple]
  Lcol<-round((as.numeric(lab.col)/100),0) # different colours
  
  Inf.plot<-data.frame(Days=Tm,Person=Pr,RSV=Gr)
  p<-ggplot(Inf.plot, aes(Days,Person, color=factor(RSV)))+geom_point(size=pnt.size)
  p + scale_y_continuous(breaks=seq(1,length(pple)),labels=colnames(Shedding.a)[pple]) +
    labs(title='') +
    theme(axis.text.y = element_text(size=lab.size,colour = Lcol))+
    scale_colour_manual(values=c("#999999", "#56B4E9", "#CC79A7")) + 
    theme( panel.background = element_rect(fill = "white",colour = "white"),
           panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                           colour = "white"),
           panel.grid.major.y = element_line(size = 0.1, linetype = 'solid',
                                             colour = "white"),
           panel.border = element_blank()) + 
    geom_hline(yintercept=hh.bound, linetype="solid", color = "grey", size=0.3)
  
}
# This plots the shedding and ARI episodes within the shedding episodes
Dot.plot2<-function(Shedding,Symptom,pple,pnt.size,lab.size){
  lng<-colSums(Shedding[,pple])
  Tm<-c();Pr<-c();Gr<-c()
  for(p in 1:length(pple)){
    # Time points with shedding
    t=which(Shedding[,pple[p]]==1);Tm<-c(Tm,t)
    Pr<-c(Pr,rep(p,length(t)))
    # Virus
    t.a<-which(Shedding[,pple[p]]==1)
    g<-rep('Virus',length(t));
    # ARI 
    t.b<-which(Symptom[,pple[p]]==1)
    g[which(t %in% intersect(t.b,t.a))]<-'ARI'
    Gr<-c(Gr,g)
  }
  lab.col<-colnames(Shedding)[pple]
  Lcol<-round((as.numeric(lab.col)/100),0) # different colours
  
  Inf.plot<-data.frame(time=Tm,Person=Pr,RSV=Gr)
  p<-ggplot(Inf.plot, aes(time,Person, color=factor(RSV)))+geom_point(size=pnt.size)
  p + scale_y_continuous(breaks=seq(1,length(pple)),labels=colnames(Shedding)[pple]) +
    labs(title=paste('Shedding and ARI patterns for ',Pathogen.name)) +
    theme(axis.text.y = element_text(size=lab.size,colour = Lcol))
  
}
################################################################################
library('lubridate')
# --- Define the Ct value cut-off for a positive sample
Pos.cut.off<-35
# --- Loading raw household data
Data<-read.csv('Household_incidence_IBM.csv')
# --- Changing date format into a form the R can run operations on
Data$sampledate<- as.Date(as.character(Data$sampledate),"%d-%B-%y")
# --- Creating a vector of sample dates based on the first and last date in the data
Sample.dates<-seq(min(Data$sampledate),max(Data$sampledate),1) 
# --- Re-coding some data variables from factor to numeric for ease of manupulation
# Away during sampling
Data$awaynfs<-as.numeric(Data$awaynfs);Data$awaynfs[Data$awaynfs==2]<-1  # not away
Data$awaynfs[Data$awaynfs==3]<-0 # away
# Cough
Data$cough<-as.numeric(Data$cough);Data$cough[Data$cough==1]<--1  
Data$cough[Data$cough==2]<-0;Data$cough[Data$cough==3]<-1 
# Nasal discharge or blocked nose
Data$runnynose<-as.numeric(Data$runnynose);Data$runnynose[Data$runnynose==1]<--1  
Data$runnynose[Data$runnynose==2]<-0;Data$runnynose[Data$runnynose==3]<-1 
# Difficulty in breathing
Data$dib<-as.numeric(Data$dib);Data$dib[Data$dib==1]<--1  
Data$dib[Data$dib==2]<-0;Data$dib[Data$dib==3]<-1 
# --- Extracting participant demographics
Demo.data<-Data[!duplicated(Data$sid),c("sid","hhid","hhsize","ageyrs") ]
# --- Extracting the number of individuals and the length of thr study period in days
N<-length(Demo.data[,1]);D<-as.numeric(max(Data$sampledate)-min(Data$sampledate)+1)     
# --- Generating ARI matrices
# DxN matrices of any of the three sypmtoms
Cough<-Mat.gen(Data$cough);Nasal<-Mat.gen(Data$runnynose);DiB<-Mat.gen(Data$dib)
# Matrix of ARI
ARI<-Cough;ARI[which(Nasal==1,TRUE)]<-1;ARI[which(DiB==1,TRUE)]<-1 # The +ves
ARI[which(Nasal==0 & ARI<0,TRUE)]<-0;ARI[which(DiB==0 & ARI<0,TRUE)]<-0 # more -ves

# --- Filling in the dates where no data was collected (Imputing) 
# RSV A
# converting data to matrix form
virus.shedding<-Mat.gen(Data$rsva);Pathogen.name<<-'RSV A'
# imputing shedding durations & generating all the other matrices needed for fitting
Shed.data<-Data.gen(virus.shedding,'A')
# Imputing ARI episodes within shedding episodes
ARI.estimated<-ARI.est(Shed.data$Shedding,Shed.data$Onset,ARI)

# RSV B
# converting data to matrix form
virus.shedding.2<-Mat.gen(Data$rsvb);Pathogen.name<<-'RSV B'
# imputing shedding durations & generating all the other matrices needed for fitting
Shed.data.2<-Data.gen(virus.shedding.2,'B')
# Imputing ARI episodes within shedding episodes
ARI.estimated.2<-ARI.est(Shed.data.2$Shedding,Shed.data.2$Onset,ARI)

# Plotting the estimated background community density curves
layout(matrix(c(0,0,1,1),nrow=2,ncol=2,byrow=T))
plot(Shed.data$Comm.risk,type='l',col='grey',lwd=2,xlab='Time in days',
     ylab='Relative rate',main='Estimated background community rate functions')
lines(Shed.data.2$Comm.risk,col='magenta',lwd=2)
legend('topright',col=c('grey','magenta'),lty=1,legend=c('RSV A','RSV B'),bty='n',
       cex=1)

# RSV as a single pathogen
Shed.data.3<-RSV.gen()
ARI.estimated.3<-ARI.estimated;ARI.estimated.3[ARI.estimated.2==1]<-1

# --- Hospital incidence based community curve
# Loading the raw data
Hosp.data<-read.csv('Hospital_incidence_IBM.csv')
# Coverting it into a form that can be used for fitting
Risk.hosp<-Comm.risk.hosp()

# --- Visualizing shedding and ARI patterns 
# RSV A shedding an ARI
shed.a<-which(colSums(Shed.data$Onset)>0);Pathogen.name<<-'RSV A'
Dot.plot2(Shedding=Shed.data$Shedding,Symptom=ARI.estimated,pple=shed.a,pnt.size=1,lab.size=5)
# RSV B shedding and ARI
shed.b<-which(colSums(Shed.data.2$Onset)>0);Pathogen.name<<-'RSV B'
Dot.plot2(Shedding=Shed.data.2$Shedding,Symptom=ARI.estimated.2,pple=shed.b,pnt.size=1,lab.size=5)
# RSV A and B shedding
Shedding<-Shed.data$Shedding;Shedding[Shed.data.2$Shedding==1]<-1
shed.rsv<-which(colSums(Shedding)>0)
house.rsv<-(Demo.data$hhid[shed.rsv]);hh.bound<-which(diff(house.rsv)>0)+1
Dot.plot(Shed.data$Shedding,Shed.data.2$Shedding,shed.rsv,1,4)
# --- Removing functions and data not needed for fitting
rm(list=c('Data','Data.gen','Mat.gen','Pos.cut.off','Pathogen.name','Hosp.data',
          'virus.shedding','virus.shedding.2','Cough','Nasal','DiB','Sample.dates',
          'ARI.est','CommRisk.func','RSV.gen','shed.a','shed.b','shed.rsv','house.rsv',
          'Dot.plot','Dot.plot2','Shedding','hh.bound','Comm.risk.hosp','ARI'))

