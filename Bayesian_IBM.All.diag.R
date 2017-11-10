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
# 19th July 2017
# This is the diagnostic function code. It assumes three chains were run
################################################################################
library('coda');library('psych');library('ggplot2');library('gridExtra')
library('MCMCvis');library('mcmcse')
trace.plotter<-function(res,res2,res3,prior.limits,cols){
  for(i in 1:length(res[1,])){
    # Traceplots
    plot(as.numeric(res[,i]),type='l',xlab='Iteration',ylab=names(res)[i],
         xlim=c(0,max(length(res[,i]),length(res2[,i]),length(res3[,i]))),
         ylim=as.numeric(c(prior.limits[i,])),main=paste('Trace plot:HH'))
    lines(as.numeric(res2[,i]),col=cols[1],lty=2);lines(as.numeric(res3[,i]),col=cols[2],lty=2)
    legend('bottomright',lty=c(1,1,1),col=c('black',cols),cex=0.5,
           legend=c('Chain1','Chain2','Chain3'))
  }
}
BnT<-function(RES,b,t){
  func<-function(res){
    # Burn-in
    Res<-res[-seq(b)]
    # Thining
    if(t>0){
      th<-seq(t,length(Res),t);Res<-Res[th];return(Res) 
    }else{
      return(Res)  
    }
  }
  return(apply(RES,2,func)) 
}
scat.h<-function(Res,X,Y,col){
  scatter.hist(Res[,X],Res[,Y],pch=20,col=col,ellipse=FALSE, smooth=FALSE,xlab=names(Res)[X],
               ylab=names(Res)[Y],title='',x.breaks=round(length(Res[,X])/500),
               y.breaks=round(length(Res[,Y])/500))
} 
mean.fn<-function(X) return(round(mean(X),3))
median.fn<-function(X)return(round(median(X),3))
Dens.plot<-function(Res,pars){
  Means<-apply(Res[,pars],2,mean.fn)
  para.name<-c(names(Res)[pars])
  dat <- data.frame(Estimates = as.vector(as.matrix(Res[,pars])),Paras = rep(c(para.name), each =  length(Res[,1])) )
  
  p<-ggplot(dat, aes(x = Estimates, fill = Paras)) + geom_density(alpha = 0.3) +  
    labs(title=paste('Parameters',paste(para.name,sep="", collapse=","),
                     'means=',paste(Means,sep="", collapse=","),'respectively')) +
    geom_vline(data=dat,aes(xintercept = rep(c(Means),each=length(Res[,1])),col=Paras))
}
Dens.plot3<-function(Res,pars){
  Means<-apply(Res[,pars],2,mean.fn)
  para.name<-c(par.names[pars])
  dat <- data.frame(Estimates = as.vector(as.matrix(Res[,pars])),Paras = rep(c(para.name), each =  length(Res[,1])) )
  
  p<-ggplot(dat, aes(x = Estimates, fill = Paras)) + geom_density(alpha = 0.3) +  
    labs(title=paste('Parameters',paste(para.name,sep="", collapse=","),
                     'means=',paste(Means,sep="", collapse=","),'respectively')) +
    geom_vline(data=dat,aes(xintercept = rep(c(Means),each=length(Res[,1])),col=Paras))
}
IBM.all.diag<-function(chain1,chain2,chain3,prior.limits){
  par.names<-names(chain1$trace[-length(chain1$trace[1,])])
  res1<-data.frame(apply(chain1$trace,2,unlist))
  res2<-data.frame(apply(chain2$trace,2,unlist))
  res3<-data.frame(apply(chain3$trace,2,unlist))
  # Trace plots
  par(mfrow=c(3,3))
  trace.plotter(res1[par.names],res2[par.names],res3[par.names],prior.limits,c('blue','cyan'))
  # Burn-in 
  ReBurn<-"Y"
  while(ReBurn=="Y"){
    burn.at<- readline("What is the burn-in?")
    burn.at<- as.numeric(unlist(strsplit(burn.at, ",")))
    thin.every<- readline("What is the thinning interval?")
    thin.every<- as.numeric(unlist(strsplit(thin.every, ",")))
    
    Res1.burn<-BnT(res1,burn.at,thin.every);Res2.burn<-BnT(res2,burn.at,thin.every)
    Res3.burn<-BnT(res3,burn.at,thin.every)
    # chain combination
    Res.all<-data.frame(rbind(Res1.burn,Res2.burn,Res3.burn),row.names=NULL);
    autocorr.plot(Res.all[,par.names])  
    ESS<-effectiveSize(Res.all[par.names]) # univariate ESS (factors in autocorrelation)
    mESS<-multiESS(as.matrix(Res.all[par.names])) # multivariate ESS(factors in correlation)
    print('The ESS per parameter=');print(ESS)
    print('The multivariate ESS =');print(mESS)
    
    ReBurn<- readline("Change burn-in and thinning?Y or N")
  }
  # Acceptance rate
 acceptance<-((chain1$acceptance.rate*(dim(chain1$trace)[1])) + (chain2$acceptance.rate*(dim(chain2$trace)[1]))
  + (chain3$acceptance.rate*(dim(chain3$trace)[1])))/
    ((dim(chain1$trace)[1]) + (dim(chain2$trace)[1]) + (dim(chain3$trace)[1]))
  # Density and scatter plots
  Res.all.plot<-exp(Res.all[par.names])
  #-------All the density plots
  par(mfrow=c(3,4))
  for(j in 1:length(par.names)){
    Mn<-round(exp(mean.fn(Res.all[,j])),3);Md<-round(exp(median.fn(Res.all[,j])),3)
    plot(density(Res.all.plot[,j]),xlab=names(Res.all.plot)[j],main=names(Res.all.plot)[j],lty=1)
    legend('topright',cex=0.7,
           legend=c(paste('Mean=',Mn),paste('Median='),Md),bty='n')
  }
  #--------Caterpillar plot
  RSV.plot<-as.matrix(Res.all[par.names]);
  #RSV.plot<-as.matrix(Res.all.plot)
  Labels<-paste(par.names,'(ESS =',round(ESS,0),')')
  Min<-min(apply(RSV.plot,2,min));Max<-max(apply(RSV.plot,2,max))
  par(mfrow=c(1,1))
  MCMCplot(RSV.plot,labels = Labels,labels_sz = 1,med_sz = 2,thick_sz = 7,
           thin_sz = 3,ax_sz = 4,xlim=c(Min,Max),xlab = "Parameter Estimate(log scale)")
  #--------Correlation plots
  #x<-which(((cor(Res.all[par.names])>0.5) & (cor(Res.all[par.names])<1)),arr.ind=TRUE)
  #y<-which((cor(Res.all[par.names])< -0.5),arr.ind=TRUE)
  #Scat.plot(3,4,c('red'))
  
  
  return(list(Res=Res.all,ESS=ESS,mESS=mESS,acceptance.rate=acceptance))
}
IBM.all.diag2<-function(chain1,chain2,chain3,prior.limits){
  par.names<-names(chain1$trace[-length(chain1$trace[1,])])
  res1<-data.frame(apply(chain1$trace,2,unlist))
  res2<-data.frame(apply(chain2$trace,2,unlist))
  res3<-data.frame(apply(chain3$trace,2,unlist))
  # Trace plots
  par(mfrow=c(3,3))
  trace.plotter(res1[par.names],res2[par.names],res3[par.names],prior.limits,c('blue','cyan'))
  # Burn-in 
  ReBurn<-"Y"
  while(ReBurn=="Y"){
    # chain 1
    burn.at<- readline("What is the chain 1 burn-in?");burn.at<- as.numeric(unlist(strsplit(burn.at, ",")))
    thin.every<- readline("What is the thinning interval?");thin.every<- as.numeric(unlist(strsplit(thin.every, ",")))
    # chain 2
    burn.at2<- readline("What is the chain2 burn-in?");burn.at2<- as.numeric(unlist(strsplit(burn.at2, ",")))
    thin.every2<- readline("What is the thinning interval?");thin.every2<- as.numeric(unlist(strsplit(thin.every2, ",")))
    # chain 3
    burn.at3<- readline("What is the chain3 burn-in?");burn.at3<- as.numeric(unlist(strsplit(burn.at3, ",")))
    thin.every3<- readline("What is the thinning interval?");thin.every3<- as.numeric(unlist(strsplit(thin.every3, ",")))
    
    Res1.burn<-BnT(res1,burn.at,thin.every);Res2.burn<-BnT(res2,burn.at2,thin.every2)
    Res3.burn<-BnT(res3,burn.at3,thin.every3)
    # chain combination
    Res.all<-data.frame(rbind(Res1.burn,Res2.burn,Res3.burn),row.names=NULL);
    autocorr.plot(Res.all[,par.names])  
    ESS<-effectiveSize(Res.all[par.names]) # univariate ESS (factors in autocorrelation)
    mESS<-multiESS(as.matrix(Res.all[par.names])) # multivariate ESS(factors in correlation)
    print('The ESS per parameter=');print(ESS)
    print('The multivariate ESS =');print(mESS)
    
    ReBurn<- readline("Change burn-in and thinning?Y or N")
  }
  # Acceptance rate
  acceptance<-((chain1$acceptance.rate*(dim(chain1$trace)[1])) + (chain2$acceptance.rate*(dim(chain2$trace)[1]))
               + (chain3$acceptance.rate*(dim(chain3$trace)[1])))/
    ((dim(chain1$trace)[1]) + (dim(chain2$trace)[1]) + (dim(chain3$trace)[1]))
  # Density and scatter plots
  Res.all.plot<-exp(Res.all[par.names])
  #-------All the density plots
  par(mfrow=c(3,4))
  for(j in 1:length(par.names)){
    Mn<-round(exp(mean.fn(Res.all[,j])),3);Md<-round(exp(median.fn(Res.all[,j])),3)
    plot(density(Res.all.plot[,j]),xlab=names(Res.all.plot)[j],main=names(Res.all.plot)[j],lty=1)
    legend('topright',cex=0.7,
           legend=c(paste('Mean=',Mn),paste('Median='),Md),bty='n')
  }
  #--------Caterpillar plot
  RSV.plot<-as.matrix(Res.all[par.names]);
  #RSV.plot<-as.matrix(Res.all.plot)
  Labels<-paste(par.names,'(ESS =',round(ESS,0),')')
  Min<-min(apply(RSV.plot,2,min));Max<-max(apply(RSV.plot,2,max))
  par(mfrow=c(1,1))
  MCMCplot(RSV.plot,labels = Labels,labels_sz = 1,med_sz = 2,thick_sz = 7,
           thin_sz = 3,ax_sz = 4,xlim=c(Min,Max),xlab = "Parameter Estimate(log scale)")
  #--------Correlation plots
  #x<-which(((cor(Res.all[par.names])>0.5) & (cor(Res.all[par.names])<1)),arr.ind=TRUE)
  #y<-which((cor(Res.all[par.names])< -0.5),arr.ind=TRUE)
  #Scat.plot(3,4,c('red'))
  
  
  return(list(Res=Res.all,ESS=ESS,mESS=mESS,acceptance.rate=acceptance))
}

DIC.calc<-function(All.res,paras,prior.func){
  Res<-All.res[names(paras)]
  posterior<-All.res["log.density"]
  prior<-apply(Res,1,prior.func)
  logL<-posterior-prior
  Deviance<--2*logL
  Pv<-0.5*var(Deviance)
  DIC<-Pv + apply(Deviance,2,mean)
  return(list(Pv=as.numeric(Pv),DIC=as.numeric(DIC)))
}

Quntile.func<-function(XX){
  quantile(XX,c(0.025,0.5,0.975))
}



