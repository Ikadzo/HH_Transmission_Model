# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 25th September 2019
################################################################################
# This code modifies the original data into structures that can be used in the 
# model(s). The modified structured are saved as CSV files that can then be 
# used in the julia environment. 
# There are three main part of data modification:
# Part 1: Extracting the shedding patterns and participant demographics
# Part 2: Extracting the spatial distance matrix
# Part 3: Extracting the patterns of the genetic clusters
################################################################################
#     Part 1:Extracting the shedding patterns and participant demographics
################################################################################
# --- Loading the functions required
source('Data_Mod_Funcs.R')
# --- Loading complete household data
Data<-read.csv('Household_incidence_IBM.csv')
# --- Loading the required libraries 
library('lubridate')

# --- Define the Ct value cut-off for a positive sample
Pos.cut.off<-35
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
Shed.data<-Data.gen(virus.shedding)
# Imputing ARI episodes within shedding episodes
ARI.estimated<-ARI.est(Shed.data$Shedding,Shed.data$Onset,ARI)

# RSV B
# converting data to matrix form
virus.shedding.2<-Mat.gen(Data$rsvb);Pathogen.name<<-'RSV B'
# imputing shedding durations & generating all the other matrices needed for fitting
Shed.data.2<-Data.gen(virus.shedding.2)
# Imputing ARI episodes within shedding episodes
ARI.estimated.2<-ARI.est(Shed.data.2$Shedding,Shed.data.2$Onset,ARI)

# RSV as a single pathogen
Shed.data.3<-RSV.gen()
ARI.estimated.3<-ARI.estimated;ARI.estimated.3[ARI.estimated.2==1]<-1

# --- Dimensions of the data, time and population size
Time<-dim(Shed.data$Shedding)[1]
Population<-dim(Shed.data$Shedding)[2]
Sample.sid<-Demo.data$sid

# --- Defining some variables needed for fiting that are group specific
# -status variable
status<-Shed.data$status;
# -Extracting matrices from the list object output of 'Data_Prep.R'
# RSV A
Shedding.a<<-as.matrix(Shed.data$Shedding);colnames(Shedding.a)<-Demo.data$sid
Onset.a<<-Shed.data$Onset;Shedding.a.est.load<<-Shed.data$Shedding.est.load;
postInf.a<<-Shed.data$postInf; 
# RSV B
Shedding.b<<-as.matrix(Shed.data.2$Shedding);colnames(Shedding.b)<-Demo.data$sid
Onset.b<<-Shed.data.2$Onset; Shedding.b.est.load<<-Shed.data.2$Shedding.est.load;
postInf.b<<-Shed.data.2$postInf; 
# -Changing covariate for infection history with 3 groups
B <<- matrix(1,Time,Population);B[postInf.a==1]<-2;B[postInf.b==1]<-3
B[(postInf.a==1 & postInf.b==1)]<-4 
# -Combined categories of viral load and symptoms
# cut-off for high viral load
high.load<-6
Load.n.sym.a<-matrix(1,Time,Population);Load.n.sym.b<-matrix(1,Time,Population)
# Reference group;shedding and asymptomatic
Load.n.sym.a[which(Shedding.a.est.load>0 & ARI.estimated!=1,TRUE)]<-2
Load.n.sym.b[which(Shedding.b.est.load>0 & ARI.estimated.2!=1,TRUE)]<-2
# low viral load and symptomatic
Load.n.sym.a[which(Shedding.a.est.load>0 & ARI.estimated==1,TRUE)]<-3
Load.n.sym.b[which(Shedding.b.est.load>0 & ARI.estimated.2==1,TRUE)]<-3
# high viral load and symptomatic
Load.n.sym.a[which(Shedding.a.est.load>high.load & ARI.estimated==1,TRUE)]<-4
Load.n.sym.b[which(Shedding.b.est.load>high.load & ARI.estimated.2==1,TRUE)]<-4
# -All the people with an onset
shed.pple.a<-which(colSums(Onset.a)>0)
shed.pple.b<-which(colSums(Onset.b)>0)
# -All the onsets, store the start, end date and the person with an onset
# RSV A
shed.times.a<-which(Onset.a==1,T);colnames(shed.times.a)<-c('start.T','Person')
shed.times.a<-as.data.frame(shed.times.a);shed.times.a$stop.T<-0
for(i in 1:length(shed.pple.a)){
  # stop time
  stop<-which(diff(Shedding.a[,shed.pple.a[i]])==-1)
  # update matrix
  shed.times.a$stop.T[which(shed.times.a$Person %in% shed.pple.a[i])]<-stop
}
# RSV B
shed.times.b<-which(Onset.b==1,T);colnames(shed.times.b)<-c('start.T','Person')
shed.times.b<-as.data.frame(shed.times.b);shed.times.b$stop.T<-0
for(i in 1:length(shed.pple.b)){
  # stop time
  stop<-which(diff(Shedding.b[,shed.pple.b[i]])==-1)
  # update matrix
  shed.times.b$stop.T[which(shed.times.b$Person %in% shed.pple.b[i])]<-stop
}

# --- Shedding data for RSV without group identification
B.rsv<-Shed.data.3$B;Shedding.rsv<-Shed.data.3$Shedding
Load.n.sym.rsv<-matrix(1,Time,Population)
# Reference group;shedding and asymptomatic
Load.n.sym.rsv[which(Shed.data.3$Shedding.est.load>0 & ARI.estimated.3!=1,TRUE)]<-2
# low viral load and symptomatic
Load.n.sym.rsv[which(Shed.data.3$Shedding.est.load>0 & ARI.estimated.3==1,TRUE)]<-3
# high viral load and symptomatic
Load.n.sym.rsv[which(Shed.data.3$Shedding.est.load>high.load & ARI.estimated.3==1,TRUE)]<-4


# --- Removing functions and data not needed for fitting
rm(list=c('Data','Data.gen','Mat.gen','Pos.cut.off','Pathogen.name',
          'virus.shedding','virus.shedding.2','Cough','Nasal','DiB',
          'ARI.est','RSV.gen','ARI','D','N','stop','i'))

################################################################################
#                 Part 2:Extracting the spatial distance matrix
################################################################################
library('Matrix');library('geosphere')
# --- Load the incomplete location data
Data.location<-read.csv("HH_members_location.csv")
# Ordering by SID
Data.location<-Data.location[order(Data.location$sid),]
# Matching the people with location info to the overall data
x<-which( Demo.data$sid %in% Data.location$sid)
Data.location$hhid<-Demo.data$hhid[x];rownames(Data.location)<-NULL
# Distance matrix for the people with location info
test<-sapply((1:dim(Data.location)[1]),
             function(x) distHaversine(Data.location[x,c(3,4)],Data.location[,c(3,4)]))
# Putting it in a larger matrix
Distance.matrix<-matrix(0,Population,Population);Distance.matrix[x,x]<-test
rownames(Distance.matrix)<-Demo.data$sid;colnames(Distance.matrix)<-Demo.data$sid
# Fill in for the people with no location info such that for pairs where one or
# both is missing location then the distance between them is:
# - average within HH pair-wise dist if i and j are in the same HH
# - distance between the infants in the respective HH if i and j are in diff HHs
missing.loc<-setdiff(1:Population,x)
for(i in missing.loc){
  hh.i<-Demo.data$hhid[i];
  for(j in setdiff(1:Population,i)){
    hh.j<-Demo.data$hhid[j]
    # If either i or j do not have location info and are housemates
    if(hh.i==hh.j){
      # Find their housemates with location info
      xx<-which(Data.location$hhid==hh.i)
      # Assign the average pair-wise distance in the HH as the missing distance
      Distance.matrix[i,j]<-sum(test[xx,xx][upper.tri(test[xx,xx])])/sum(upper.tri(test[xx,xx]))
      Distance.matrix[j,i]<-Distance.matrix[i,j]
    }else{
      # Find the infant in the respective household(ordered by age for just the 
      # infant is the first person in the household)
      yy.i<-which(Data.location$hhid==hh.i);yy.i<-yy.i[1]
      yy.j<-which(Data.location$hhid==hh.j);yy.j<-yy.j[1]
      # Distance between i and j is the distance between the infants
      Distance.matrix[i,j]<-test[yy.i,yy.j];Distance.matrix[j,i]<-Distance.matrix[i,j]
    }
  }
}
Distance.matrix<-as.matrix(forceSymmetric(Distance.matrix,'L'))
Distance.matrix<-Distance.matrix/1000 # converting to km
# --- Removing functions and data not needed for fitting
rm(list=c('Data.location','x','test','missing.loc','i','j','hh.i','hh.j','yy.i',
          'yy.j','xx'))

################################################################################
#            Part 3:Extracting the patterns of the genetic clusters
################################################################################

# --- Loading the genetic data 
# Sequences in matrix format 
load('rsv.sequences.final.RData')
# Cluster ids
rsv.cluster<-read.csv("rsv.clusters.csv")

seq.id=rsv.cluster$tree_id
# --- Matching the data
# - RSV A
seq.id.a=rownames(rsv.sequences$rsva)
cl.ind.a=match(seq.id.a,seq.id)
# - RSV B
seq.id.b=rownames(rsv.sequences$rsvb)
cl.ind.b=match(seq.id.b,seq.id)

# --- Extracting cluster id for the sequences 
# - RSV A
# Extracting cluster id
temp.a=rsv.cluster$clade[cl.ind.a]
cl.name.a=as.character(unique(temp.a)[order(unique(temp.a))])
# replacing cluster labels with numbers
Seq.cluster.id.a=rep(0,length(temp.a));
for(i in 1:length(cl.name.a)){
  Seq.cluster.id.a[temp.a==cl.name.a[i]]=i
}
# - RSV B
# Extracting cluster id
temp.b=rsv.cluster$clade[cl.ind.b]
cl.name.b=as.character(unique(temp.b)[order(unique(temp.b))])
# replacing cluster labels with numbers
Seq.cluster.id.b=rep(0,length(temp.b));
for(i in 1:length(cl.name.b)){
  Seq.cluster.id.b[temp.b==cl.name.b[i]]=i
}

# --- Extracting SID and dates
# - SID
seqs.sid.a<-as.numeric(strtrim(labels(rsv.sequences$rsva),4)) # RSV A
seqs.sid.b<-as.numeric(strtrim(labels(rsv.sequences$rsvb),4)) # RSV B
# -Dates
strReverse <- function(x)sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
Rev.labs<-strReverse(labels(rsv.sequences$rsva))
seqs.dates.a<-as.Date(strReverse(strtrim(Rev.labs,11)),"%d-%B-%Y") 
Rev.labs<-strReverse(labels(rsv.sequences$rsvb))
seqs.dates.b<-as.Date(strReverse(strtrim(Rev.labs,11)),"%d-%B-%Y") 

# --- Distance matrix from raw genetic distances. Deletion of missing sites is
#  pair-wise
Nucleotide.D.a <- as.data.frame(as.matrix(dist.dna(rsv.sequences$rsva, model = "N",pairwise.deletion=T)))
Nucleotide.D.b <- as.data.frame(as.matrix(dist.dna(rsv.sequences$rsvb, model = "N",pairwise.deletion=T)))

# --- People and dates with sequence information relative to the entire data set
cl.pos.a<-cbind(day=match(seqs.dates.a, Sample.dates),sid=match(seqs.sid.a, Sample.sid))
cl.pos.b<-cbind(match(seqs.dates.b, Sample.dates),match(seqs.sid.b, Sample.sid))

# --- Imputing clusters for episodes and household outbreaks with some genetic info
Shedding.cluster.a<-Cluster.fill(Seq.cluster.id.a,cl.pos.a,shed.times.a)
Shedding.cluster.b<-Cluster.fill(Seq.cluster.id.b,cl.pos.b,shed.times.b)

# --- Calculating the pair-wise genetic distances between cases
dist.gen.cases.a<-Dist.gen.cases(Nucleotide.D.a,cl.pos.a,shed.times.a)
dist.gen.cases.b<-Dist.gen.cases(Nucleotide.D.b,cl.pos.b,shed.times.b)

# --- For every cluster, generate vectors of unique pair-wise distance
dist.cl.a<-list()
for(cl in 1:max(Seq.cluster.id.a)){
  dist.cl<-Nucleotide.D.a[Seq.cluster.id.a==cl,Seq.cluster.id.a==cl]
  if(length(dist.cl)>1) dist.cl.a[[cl]]<-unique(dist.cl[upper.tri(dist.cl)])
  if(length(dist.cl)==1) dist.cl.a[[cl]]<- dist.cl
}
dist.cl.b<-list()
for(cl in 1:max(Seq.cluster.id.b)){
  dist.cl<-Nucleotide.D.b[Seq.cluster.id.b==cl,Seq.cluster.id.b==cl]
  if(length(dist.cl)>1) dist.cl.b[[cl]]<-unique(dist.cl[upper.tri(dist.cl)])
  if(length(dist.cl)==1) dist.cl.b[[cl]]<- dist.cl
}

# --- One case per HH outbreak
# unclustered shedding episodes
x.a<-which(Onset.a==1 & Shedding.cluster.a==0,T)
x.b<-which(Onset.b==1 & Shedding.cluster.b==0,T)
# One case from each outbreak
one.case.a<-one.case.per.hh.outbreak(x.a,Shedding.a,Shedding.cluster.a)
one.case.b<-one.case.per.hh.outbreak(x.b,Shedding.b,Shedding.cluster.b)

# --- Initial random cluster assignment for missing clusters
Shedding.cluster.A<-Cluster.fill.random(Shedding.cluster.a,Shedding.a,Onset.a)
Shedding.cluster.B<-Cluster.fill.random(Shedding.cluster.b,Shedding.b,Onset.b)

# --- Initial number of HH outbreaks in each cluster
init.rsva.outbreaks<-outbreak.counter(Shedding.cluster.A)
init.rsvb.outbreaks<-outbreak.counter(Shedding.cluster.B)

rm(list=c('temp.a','temp.b','cl.name.a','cl.name.b','cl','x.a','x.b'))

################################################################################
#                        Saving the data as CSV files
################################################################################
# 0/1 matrices of shedding episodes
write.table(Shedding.a, file = "Data_from_R/Shedding.a.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
write.table(Shedding.b, file = "Data_from_R/Shedding.b.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
write.table(Shedding.rsv, file = "Data_from_R/Shedding.rsv.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
# Categorical variable for infection history (3 groups)
write.table(B, file = "Data_from_R/B.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
write.table(B.rsv, file = "Data_from_R/B.rsv.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
# Combined categories of viral load and symptoms
write.table(Load.n.sym.a, file = "Data_from_R/Load.n.sym.a.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
write.table(Load.n.sym.b, file = "Data_from_R/Load.n.sym.b.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
write.table(Load.n.sym.rsv, file = "Data_from_R/Load.n.sym.rsv.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
# 0/1 matrix of presence or absence in the household
write.table(status, file = "Data_from_R/status.csv",row.names=FALSE, na="",col.names=Demo.data$sid, sep=",")
# Demographic information
write.table(Demo.data, file = "Data_from_R/Demo.data.csv",row.names=FALSE, na="", sep=",")
# Spatial distance matrix
write.table(Distance.matrix, file = "Data_from_R/Distance.matrix.csv",row.names=FALSE,col.names=Demo.data$sid, na="", sep=",")
# For every onset, the start and end day and the person
write.table(shed.times.a, file = "Data_from_R/shed.times.a.csv",row.names=FALSE, na="", sep=",")
write.table(shed.times.b, file = "Data_from_R/shed.times.b.csv",row.names=FALSE, na="", sep=",")
# The onset time and person index of one case from every HH outbreak
write.table(one.case.a, file = "Data_from_R/one.case.a.csv",row.names=FALSE, na="", sep=",")
write.table(one.case.b, file = "Data_from_R/one.case.b.csv",row.names=FALSE, na="", sep=",")
# People and dates with sequence information relative to the entire data set
write.table(cl.pos.a, file = "Data_from_R/cl.pos.a.csv",row.names=FALSE, na="", sep=",")
write.table(cl.pos.b, file = "Data_from_R/cl.pos.b.csv",row.names=FALSE, na="", sep=",")
# Pair-wise genetic distances between cases with seq data
write.table(as.matrix(dist.gen.cases.a), file = "Data_from_R/dist.gen.cases.a.csv",row.names=FALSE, na="", sep=",")
write.table(as.matrix(dist.gen.cases.b), file = "Data_from_R/dist.gen.cases.b.csv",row.names=FALSE, na="", sep=",")
# Genetic distance matrices
write.table(Nucleotide.D.a, file = "Data_from_R/Nucleotide.D.a.csv",row.names=FALSE, na="", sep=",")
write.table(Nucleotide.D.b, file = "Data_from_R/Nucleotide.D.b.csv",row.names=FALSE, na="", sep=",")
# Clustering of sequences using a cut-off of 10 nucleotides
write.table(Seq.cluster.id.a, file = "Data_from_R/Seq.cluster.id.a.csv",row.names=FALSE, na="", sep=",")
write.table(Seq.cluster.id.b, file = "Data_from_R/Seq.cluster.id.b.csv",row.names=FALSE, na="", sep=",")
# Imputed clusters for episodes and household outbreaks with some genetic info
write.table(Shedding.cluster.a, file = "Data_from_R/Shedding.cluster.fixedA.csv",row.names=FALSE, na="", sep=",")
write.table(Shedding.cluster.b, file = "Data_from_R/Shedding.cluster.fixedB.csv",row.names=FALSE, na="", sep=",")
# Initial number of HH outbreaks in each cluster
write.table(init.rsva.outbreaks, file = "Data_from_R/init.rsva.outbreaks.csv",row.names=FALSE, na="", sep=",")
write.table(init.rsvb.outbreaks, file = "Data_from_R/init.rsvb.outbreaks.csv",row.names=FALSE, na="", sep=",")
# Initial random cluster assignment for missing clusters
write.table(Shedding.cluster.A, file = "Data_from_R/Shedding.cluster.A.csv",row.names=FALSE, na="", sep=",")
write.table(Shedding.cluster.B, file = "Data_from_R/Shedding.cluster.B.csv",row.names=FALSE, na="", sep=",")


# ---------- FUNCTIONS  to be replicated
nLL=rep(0,10)
system.time(for(i in 1:10)nLL[i]<-IBM.cluster.LL(theta,Shedding.cluster.A,Shedding.cluster.B))

Shedding.cluster.A=read.csv("Julia_Code/Data_from_R/Shedding.cluster.A.csv")
Shedding.cluster.B=read.csv("Julia_Code/Data_from_R/Shedding.cluster.B.csv")
status2=read.csv("Julia_Code/Data_from_R/status.csv")
#


