# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 26th September 2019
################################################################################
# This code contains function needed to modify the data and is called by 
# 'Data_Mod.R'
################################################################################
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
Data.gen<-function(virus.shedding){
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
  # Function for interpolating the status matrix. It need info presence during NFS collection
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
  # Dimensions of teh data
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
  
  Results<-list(Shedding=Shedding.a,Onset=Onset.a,Shedding.est.ct=Shedding.a.est.ct,
                Shedding.est.load=Shedding.a.est.load,status=status,atRisk=atRisk.a,
                postInf=postInf.a)
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
  
  Res=list(Shedding=Shedding,Onset=Onset,Shedding.est.load=Shedding.est.load,
           status=status,atRisk=atRisk,B=B)
  return(Res)
}

# Function 1: Function for imputing clusters for shedding episodes with some 
#             sequence info and cases with no sequence info but are part of a 
#             household outbreak with some sequences
# ------------------------------------------------------------------------------
# Needs:
#  - Cluster id for the sequences (Seq.cluster.id)
#  - The position of the sequences in the shedding matrix (cl.pos)
#  - Matrix identifying start time, index and end time of cases(shed.times) 
#  - Matrix identifying each household outbreak
Cluster.fill<-function(Seq.cluster.id,cl.pos,shed.times){
  # Identifying cluster group for positions(people and dates) with sequences
  Seq.clusters<-matrix(0,length(Sample.dates),length(Sample.sid))
  Seq.clusters[cl.pos]<-Seq.cluster.id
  # Predefine matrix of imputed clusters
  Seq.clusters.imp<-0*Seq.clusters
  # Keep track of the episodes that have been assigned
  seq.epis<-c()
  
  # 1) # Fill in the gaps for the episodes with some sequence data
  for(e in 1:length(shed.times$start.T)){
    # Cluster info for the shedding duration
    Shed.cl<-Seq.clusters[shed.times$start.T[e]:shed.times$stop.T[e],shed.times$Person[e]]
    if(any(Shed.cl>0)){
      # Keep track of the episode
      seq.epis<-c(seq.epis,e)
      # Cluster ids within that episode
      epi<-Seq.clusters[shed.times$start.T[e]:shed.times$stop.T[e],shed.times$Person[e]]
      Cl.id<-unique(epi[epi>0])
      # If only one cluster, assign that id to the whole episode
      if(length(Cl.id)==1) Seq.clusters.imp[shed.times$start.T[e]:shed.times$stop.T[e],shed.times$Person[e]]<-Cl.id
      # If more that one cluster, distribute the ids across the episode
      if(length(Cl.id)>1){
        # Time indices for the episode
        epi.ind<-shed.times$start.T[e]:shed.times$stop.T[e]
        # Relative to the episode, start and end times for each cluster
        strt<-c(1,(which(epi %in% Cl.id)+1));strt<-strt[-length(strt)]
        end<-which(epi %in% Cl.id);end<-c(end[-length(end)],length(epi))
        # Assign cluster id
        for(cl in 1:length(strt)){
          Seq.clusters.imp[epi.ind[strt[cl]:end[cl]],shed.times$Person[e]] <- epi[epi>0][cl]
        } 
      }
    } 
  }
  
  # Episodes still unclustered
  comp.missing<-setdiff(seq(1:length(shed.times$start.T)),seq.epis)
  # 2) Impute cluster information for episodes with no sequences but part of a HH
  #    outbreak with a unique cluster id. Assignment is done in order of time
  #    i.e. the earliest case is assigned first. The gaps for cases with some 
  #    genetic info had to be filled first, hence these remaining cases are in a
  #    separate loop. 
  for(i in comp.missing){
    # Case details
    x<-shed.times[i,]
    # Find their HH and any cluster info
    hh.case<-Demo.data$hhid[x$Person]
    # If the HH has cluster info
    if (any(rowSums(Seq.clusters.imp[,Demo.data$hhid==hh.case])>0)){
      # For the given case, define the time window within which to look for 
      # cluster info, i.e. the 5 days pre-onset to 5 days post end of shedding
      win<-(x$start.T-5):(x$stop.T+5)
      # Days in the window with cluster ids
      days<-which(rowSums(Seq.clusters.imp[win,Demo.data$hhid==hh.case])>0)
      # if there is any info on cluster id proceed
      if (length(days)>0){
        # Earliest day in that window with a cluster id
        day<-win[min(days)]
        # cluster option(s)
        cluster<-unique(Seq.clusters.imp[day,Demo.data$hhid==hh.case]);
        cluster<-cluster[cluster>0]
        # If there is only one cluster option, assign the case the cluster id for 
        # the chosen day
        if(length(cluster)==1) {
          Seq.clusters.imp[x$start.T:x$stop.T,x$Person]<-cluster
          # Keep track of assigned case
          seq.epis<-c(seq.epis,i)
        }
      }
    }
  }
  # Return imputed shedding clusters
  return(Seq.clusters.imp)
}

# ------------------------------------------------------------------------------
# Function 2: Function for obtaining the pair-wise genetic distance between 
#             episodes with sequence data.
# ------------------------------------------------------------------------------
# Needs:
#  - Matrix of genetic distances between sequences (gen.dist)
#  - The position of the sequences in the shedding matrix (cl.pos)
#  - Matrix identifying start time, index and end time of cases(shed.times)
Dist.gen.cases<-function(gen.dist,cl.pos,shed.times){
  # Start by identifying episodes with genetic info
  shed.times$seq<-0
  for (e in 1:length(shed.times$start.T)){
    # Shedding window
    dur.e<-shed.times$start.T[e]:shed.times$stop.T[e]
    # any sequences in the window?
    x<-which((cl.pos[,1] %in% dur.e) & cl.pos[,2]==(shed.times$Person[e]))
    if(length(x)>0) shed.times$seq[e]<-1
  }
  sum(shed.times$seq)
  # For every episode get the pairwise distance
  # - pre-allocate matrix
  Dist.gen.epi<-matrix(NA,length(shed.times$seq),length(shed.times$seq))
  for(i in 1:(length(shed.times$seq)-1)){
    for(j in (i+1):length(shed.times$seq)){
      # If both episodes have some sequence data
      if(shed.times$seq[i]==1 & shed.times$seq[j]==1){
        # - Find the first sequence from i
        # Shedding window
        dur.i<-shed.times$start.T[i]:shed.times$stop.T[i]
        # any sequences in the window?
        x<-which((cl.pos[,1] %in% dur.i) & cl.pos[,2]==(shed.times$Person[i]))
        # pick the first one in order of time
        x<-x[which(cl.pos[x,1]==min(cl.pos[x,1]))]
        
        # - Find the closest temporal sequence from j
        # Shedding window
        dur.j<-shed.times$start.T[j]:shed.times$stop.T[j]
        # any sequences in the window?
        y<-which((cl.pos[,1] %in% dur.j) & cl.pos[,2]==(shed.times$Person[j]))
        # time between j's sequences and i's
        t.btwn<-abs(cl.pos[y,1]-cl.pos[x,1])
        # pick the one with the shortest time(incase of equal forward and backward
        # temporal distance, prioritise backward i.e. j's seq prior to i's)
        y<-y[which(t.btwn==min(t.btwn))[1]]
        
        # - Find the genetic distance between the relevant sequences
        Dist.gen.epi[i,j]<-gen.dist[x,y]
        
      }else{
        # Assign an indicator value out of the range of possible genetic distances
        # so that the model will know to randomly assign a distance from the 
        # relevant cluster
        Dist.gen.epi[i,j]<--1
      }
    }
  }
  diag(Dist.gen.epi)<-0
  Dist.gen.epi<-forceSymmetric(Dist.gen.epi)
  
  return(Dist.gen.epi)
}

# Need another function for randomly assigning the clusters for HH outbreaks
# in the likelihood function for cases left unclustered after running the
# Cluster.fill function. This is used to initialize the cluster pattern prior to
# running the MCMC where 'Cluster.update' is used in the likelihood function
Cluster.fill.random<-function(Shedding.cluster,Shedding,Onset){
  # Available cluster ids
  cl.id<-unique(Shedding.cluster[Shedding.cluster>0])
  # All the episodes
  all.epis<-which(Onset>0,T)
  # The clustered episodes
  cl.epis<-which(Shedding.cluster>0 & Onset==1,T)
  # The difference
  x<-which(Onset==1 & Shedding.cluster==0,T)
  # Identify the household outbreaks in the unclustered cases
  # Unclustered shedding patterns
  Shedding.un<-Shedding;Shedding.un[Shedding.cluster>0]<-0
  # HHs
  hhs<-unique(Demo.data$hhid[x[,2]])
  # For each HH, cluster the outbreaks
  for(h in hhs){
    # Shedding times
    times<-which(rowSums(Shedding.un[,Demo.data$hhid==h])>0)
    # Time when each HH outbreak started and ended
    start<-times[c(1,(which(diff(times)>5)+1))]
    end<-times[c(which(diff(times)>5),length(times))]
    # Pick at random the cluster to assign to each outbreak
    cl.assign<-sample(x = cl.id,size = length(start),replace = T)
    for(e in 1:length(end)){
      Shedding.cluster[start[e]:end[e],Demo.data$hhid==h][Shedding.un[start[e]:end[e],Demo.data$hhid==h]==1]<-cl.assign[e]
    }
  }
  
  return(Shedding.cluster)
}
# This function gets the indices for time and person for one case from each
# HH outbreak. Doing this so that I can check how the ids are being infered
# in the MCMC chain
one.case.per.hh.outbreak<-function(x,Shedding,Shedding.cluster.certain){
  # matrix for keeping track of 1 case from every outbreak
  one.case<-data.frame(times=c(),pple=c())
  # their HHs
  hhs<-Demo.data$hhid[x[,2]]
  hh.shed<-unique(hhs)
  # for each HH, find the number of outbreaks
  for(hh in hh.shed){
    # Shedding times
    times<-which(rowSums((Shedding==1 & Shedding.cluster.certain==0)[,Demo.data$hhid==hh])>0)
    # if differences between shedding times >5 those are diff outbreaks
    x<-which(diff(times)>5)
    # end times of each outbreak
    times<-times[c(x,length(times))]
    # find at least one person shedding at these times
    y<-which((Shedding==1 & Shedding.cluster.certain==0)[times,Demo.data$hhid==hh]==1,T)
    if(length(y)==1) {
      pple<-which(Demo.data$hhid==hh)[y]
    }else{
      pple<-which(Demo.data$hhid==hh)[y[,2]]
      times<-times[y[,1]]
    }
    # store the info
    # if one person has > epi in separate outbreaks
    one.case<-rbind(one.case,cbind(times,pple)) 
  }
  return(one.case)
  
}
# Find the initial number of HH outbreaks in each cluster
outbreak.counter<-function(Shedding.cluster){
  No.cls<-max(Shedding.cluster)
  outbreak.count<-c(cluster=rep(0,No.cls))
  # for each cluster, count the HH outbreaks
  for(i in 1:No.cls){
    # Pple and shedding times
    pple.time<-as.data.frame(which(Shedding.cluster==i,T));colnames(pple.time)<-c('time','pple')
    # pple shedding
    pple.shed<-unique(pple.time$pple)
    # their HHs
    pple.time$hhs<-Demo.data$hhid[pple.time$pple]
    hh.shed<-unique(pple.time$hhs)
    # for each HH, find the number of outbreaks
    for(hh in hh.shed){
      # Shedding times
      times<-which(rowSums((Shedding.cluster==i)[,Demo.data$hhid==hh])>0)
      # if differences between shedding times >5 those are diff outbreaks
      x<-which(diff(times)>5)
      # Add on to count
      outbreak.count[i]<-outbreak.count[i] + length(x) + 1
    }
  }
  return(outbreak.count)
}
