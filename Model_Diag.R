# This function is called by 'Model_Run.R' to run diagnostics on the MCMC runs
# Convergence is assessed visually and hence the user has to input the burn-in 
# and thinning. This works for 3 MCMC chains. 
# The main diagnostic function returns the estimated posterior distribution, the 
# univariate effective sample size for each parameter, the multivariate effective
# sample size and the acceptance rate 
# ################################################################################
# Required libraries
library('coda');library('psych');library('ggplot2');library('gridExtra')
library('MCMCvis');library('mcmcse');library('mcmcplots');library('GGally')
library('reshape2')
# Function for plotting trace plots
trace.plotter<-function(res1,res2,res3,prior.limits,cols,desc){
  for(i in 1:length(res1[1,])){
    # Traceplots
    plot(as.numeric(res1[,i]),type='l',xlab='Iteration',ylab=names(res1)[i],
         xlim=c(0,max(length(res1[,i]),length(res2[,i]))),
         ylim=as.numeric(c(prior.limits[i,])),main=desc)
    lines(as.numeric(res2[,i]),col=cols[1],lty=2);
    lines(as.numeric(res3[,i]),col=cols[2],lty=2);
    legend('topright',lty=c(1,1),col=c('black',cols),cex=0.6,bty='n',
           legend=c('Chain1','Chain2','Chain3'))
  }
}
# Function for burning and thinning the MCMC results
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
# Overall diagnostic function
IBM.all.diag<-function(chain1,chain2,chain3,prior.limits,desc){
  par.names<-names(chain1$trace[-length(chain1$trace[1,])])
  res1<-data.frame(apply(chain1$trace,2,unlist));names(res1)<-c(par.names,"log.density")
  res2<-data.frame(apply(chain2$trace,2,unlist));names(res2)<-c(par.names,"log.density")
  res3<-data.frame(apply(chain3$trace,2,unlist));names(res3)<-c(par.names,"log.density")
  
  # Trace plots
  par(mfrow=c(2,2))
  trace.plotter(res1[par.names],res2[par.names],res3[par.names],prior.limits,c('red','cyan'),desc)
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
  acceptance<-((chain1$acceptance.rate*(dim(chain1$trace)[1])) + 
                 (chain2$acceptance.rate*(dim(chain2$trace)[1])) + 
                 (chain3$acceptance.rate*(dim(chain3$trace)[1]))) /
    ((dim(chain1$trace)[1]) + (dim(chain2$trace)[1]) + (dim(chain3$trace)[1]))
  # Displaying results
  #--------Caterpillar plot
  RSV.plot<-as.matrix(Res.all[par.names]);
  Labels<-paste(par.names,'(ESS =',round(ESS,0),')')
  Min<-min(apply(RSV.plot,2,min));Max<-max(apply(RSV.plot,2,max))
  par(mfrow=c(1,1))
  MCMCplot(RSV.plot,labels = Labels,labels_sz = 1,med_sz = 2,thick_sz = 7,
           thin_sz = 3,ax_sz = 4,xlim=c(Min,Max),xlab = "Parameter Estimate(log scale)")
  #--------Correlation plots
  iter<-dim(Res.all)[1];take<-sample(seq(iter),500)
  pars<-dim(Res.all)[2]-1
  pairs.panels(Res.all[take,seq(pars)], 
               method = "pearson", # correlation method
               hist.col = "#00AFBB",
               density = TRUE,  # show density plots
               ellipses = TRUE # show correlation ellipses
  )
  
  return(list(Res=Res.all,ESS=ESS,mESS=mESS,acceptance.rate=acceptance))
}
# Function for calculating the median and 95% credible interval
Quntile.func<-function(XX){
  quantile(XX,c(0.025,0.5,0.975))
}
# This function is passed the posterior distributions of the individual parameters
# and plots the distributions of the per person household transmission rate and
# community transmission rate. This is Figure 3. in the main article.
Rates.plot<-function(Par.dist){
  # Taking a sample from the posterior rather than every value
  set.seed(5);take<-sample(seq(dim(Par.dist)[1]),1000)
  
  # RSV A: Calculating the per person pair-wise within household risk of exposure
  # for different combinations of categories
  # small HH and asymptomatic infectious contact
  small.asym<-exp(Par.dist[['eta.A']][take]) 
  # small HH and symptomatic low viral load infectious contact
  small.sym.low<-exp(Par.dist[['eta.A']][take]) * exp(Par.dist[['LowSym']][take]) 
  # small HH and symptomatic high viral load infectious contact
  small.sym.high<-exp(Par.dist[['eta.A']][take]) * exp(Par.dist[['HighSym']][take]) 
  # large HH and asymptomatic infectious contact
  large.asym<-exp(Par.dist[['eta.A']][take]) 
  # large HH and symptomatic low viral load infectious contact
  large.sym.low<-exp(Par.dist[['eta.A']][take]) * exp(Par.dist[['LowSym']][take]) *  exp(Par.dist[['hh.size']][take])
  # large HH and symptomatic high viral load infectious contact
  large.sym.high<-exp(Par.dist[['eta.A']][take]) * exp(Par.dist[['HighSym']][take]) * exp(Par.dist[['hh.size']][take])
  # Putting the above calculations together for plotting purposes
  A.data<-c(small.asym, small.sym.low, small.sym.high,
            large.asym, large.sym.low, large.sym.high)
  HH.size<-rep(c('Small','Large'),each=3000)
  Infectiousness=rep(c('Asymptomatic','Symptomatic & low load','Symptomatic & high load'),each=1000)
  data.a<-data.frame(HH.size,Infectiousness,A.data)
  # The plot for RSV A within HH rates
  p1<-ggplot(data.a, aes(x=HH.size, y=A.data, fill=Infectiousness)) +
    geom_boxplot(outlier.shape = NA) + theme(legend.position="none") + 
    labs(y = "Pair-wise rate of within household exposure /susceptible /day", x= "Household size",title='RSV A') + 
    scale_y_continuous(limits = quantile(data.a$A.data, c(0.025, 0.975))) +
    theme(axis.title.y=element_text(size=rel(1), angle=90))
  
  # RSV B: Calculating the per person pair-wise within household risk of exposure
  # for different combinations of categories
  # small HH and asymptomatic infectious contact
  small.asym<-exp(Par.dist[['eta.B']][take]) 
  # small HH and symptomatic low viral load infectious contact
  small.sym.low<-exp(Par.dist[['eta.B']][take]) * exp(Par.dist[['LowSym']][take]) 
  # small HH and symptomatic high viral load infectious contact
  small.sym.high<-exp(Par.dist[['eta.B']][take]) * exp(Par.dist[['HighSym']][take]) 
  # large HH and asymptomatic infectious contact
  large.asym<-exp(Par.dist[['eta.B']][take]) 
  # large HH and symptomatic low viral load infectious contact
  large.sym.low<-exp(Par.dist[['eta.B']][take]) * exp(Par.dist[['LowSym']][take]) *  exp(Par.dist[['hh.size']][take])
  # large HH and symptomatic high viral load infectious contact
  large.sym.high<-exp(Par.dist[['eta.B']][take]) * exp(Par.dist[['HighSym']][take]) * exp(Par.dist[['hh.size']][take])
  # Putting the above calculations together for plotting purposes
  B.data<-c(small.asym, small.sym.low, small.sym.high,
            large.asym, large.sym.low, large.sym.high)
  data.b<-data.frame(HH.size,Infectiousness,B.data)
  # The plot for RSV B within HH rates
  p2<-ggplot(data.b, aes(x=HH.size, y=B.data, fill=Infectiousness)) +
    geom_boxplot(outlier.shape = NA) + theme(legend.position="none") + 
    labs(y = "", x= "Household size",title='RSV B') +
    scale_y_continuous(limits = quantile(data.a$A.data, c(0.025, 0.975))) 
  
  # Calculating the per person community risk of exposure for different
  # categories at different time points
  Comm.age1.a<-matrix(0,1000,D);Comm.age2.a<-matrix(0,1000,D);Comm.age3.a<-matrix(0,1000,D)
  Comm.age1.b<-matrix(0,1000,D);Comm.age2.b<-matrix(0,1000,D);Comm.age3.b<-matrix(0,1000,D)
  Comm.risk.a<-Shed.data$Comm.risk;Comm.risk.b<-Shed.data.2$Comm.risk
  for(i in 1:D){
    # RSV A
    # First age category
    Comm.age1.a[,i]<-Comm.risk.a[i]*exp(Par.dist[['epsilon.A']][take])
    # Second age category
    Comm.age2.a[,i]<-Comm.risk.a[i]*exp(Par.dist[['epsilon.A']][take])*exp(Par.dist[['eps.age2']][take])
    # Third age category
    Comm.age3.a[,i]<-Comm.risk.a[i]*exp(Par.dist[['epsilon.A']][take])*exp(Par.dist[['eps.age3']][take])
    # RSV B
    # First age category
    Comm.age1.b[,i]<-Comm.risk.b[i]*exp(Par.dist[['epsilon.B']][take])
    # Second age category
    Comm.age2.b[,i]<-Comm.risk.b[i]*exp(Par.dist[['epsilon.B']][take])*exp(Par.dist[['eps.age2']][take])
    # Third age category
    Comm.age3.b[,i]<-Comm.risk.b[i]*exp(Par.dist[['epsilon.B']][take])*exp(Par.dist[['eps.age3']][take])
  }
  # From the calculated communty risk, finding the lower are upper bound
  x<-which(Comm.age1.a==min(Comm.age1.a),TRUE)[1,1];y<-which(Comm.age1.a==max(Comm.age1.a),TRUE)[1,1]
  min.comm.age1.a<-Comm.age1.a[x,];max.comm.age1.a<-Comm.age1.a[y,]
  x<-which(Comm.age2.a==min(Comm.age2.a),TRUE)[1,1];y<-which(Comm.age2.a==max(Comm.age2.a),TRUE)[1,1]
  min.comm.age2.a<-Comm.age2.a[x,];max.comm.age2.a<-Comm.age2.a[y,]
  x<-which(Comm.age3.a==min(Comm.age3.a),TRUE)[1,1];y<-which(Comm.age3.a==max(Comm.age3.a),TRUE)[1,1]
  min.comm.age3.a<-Comm.age3.a[x,];max.comm.age3.a<-Comm.age3.a[y,]
  
  
  x<-which(Comm.age1.b==min(Comm.age1.b),TRUE)[1,1];y<-which(Comm.age1.b==max(Comm.age1.b),TRUE)[1,1]
  min.comm.age1.b<-Comm.age1.b[x,];max.comm.age1.b<-Comm.age1.b[y,]
  x<-which(Comm.age2.b==min(Comm.age2.b),TRUE)[1,1];y<-which(Comm.age2.b==max(Comm.age2.b),TRUE)[1,1]
  min.comm.age2.b<-Comm.age2.b[x,];max.comm.age2.b<-Comm.age2.b[y,]
  x<-which(Comm.age3.b==min(Comm.age3.b),TRUE)[1,1];y<-which(Comm.age3.b==max(Comm.age3.b),TRUE)[1,1]
  min.comm.age3.b<-Comm.age3.b[x,];max.comm.age3.b<-Comm.age3.b[y,]
  # Putting the above values together for plotting purposes
  Days <- c(1:D)
  Dat <- data.frame(Days, min.comm.age1.a, max.comm.age1.a, min.comm.age2.a,
                    max.comm.age2.a, min.comm.age3.a, max.comm.age3.a,
                    min.comm.age1.b, max.comm.age1.b, min.comm.age2.b,
                    max.comm.age2.b, min.comm.age3.b, max.comm.age3.b)
  # The plot for RSV A community rates
  p1.c <- ggplot(Dat, aes(x= Days)) + 
    geom_ribbon(aes(ymin= min.comm.age3.a, ymax= max.comm.age3.a,alpha=0.1), fill="black") +
    geom_ribbon(aes(ymin= min.comm.age1.a, ymax= max.comm.age1.a,alpha=0.1), fill="cyan") +
    geom_ribbon(aes(ymin= min.comm.age2.a, ymax= max.comm.age2.a,alpha=0.1), fill="red") +
    labs(y = "Community rate of exposure /susceptible /day",title='RSV A') + 
    theme(legend.position="none") +
    scale_y_continuous(limits = c(0,0.035)) 
  # The plot for RSV B community rates
  p2.c <- ggplot(Dat, aes(x= Days)) + 
    geom_ribbon(aes(ymin= min.comm.age3.b, ymax= max.comm.age3.b,alpha=0.1), fill="black") +
    geom_ribbon(aes(ymin= min.comm.age2.b, ymax= max.comm.age2.b,alpha=0.1), fill="red") +
    geom_ribbon(aes(ymin= min.comm.age1.b, ymax= max.comm.age1.b,alpha=0.1), fill="cyan") +
    labs(y = "",title='RSV B') + 
    theme(legend.position="none")+
    scale_y_continuous(limits = c(0,0.035)) 
  grid.arrange(p1,p2,p1.c,p2.c,nrow=2)
  
}