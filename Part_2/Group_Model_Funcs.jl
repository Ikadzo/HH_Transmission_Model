# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 24th September 2019
################################################################################
# This code contains functions needed for running the model that uses data
# identified at the pathogen group level for the purpose of making comparisons
# to the model in 'Cluster_Model_Run.jl'. Unlike the model in the publication
# 'Model based estimates of transmission of respiratory syncytial virus within
# households' this version modifies susceptibility based on current shedding
# status and uses an emperical based background community function.
# It is called by the function 'Group_Model_Run.jl'. There are 5 functions in
# this file:
#     1) CommRisk_func
#     2) IBM_LL
#     3) IBM_prior
#     4) IBM_posterior
#     5) MH_MCMC
################################################################################
################################################################################
# This function generates the empirical based background community function. It
# takes in 5 arguments:
# delta = the parameter for the base background rate prior to any observed onset
# beta = the rate of decay parameter for background community function
# Seq_clusters = an array of the range of cluster ids
# Shedding_cluster = A matrix of all shedding episodes identified by cluster type/id
# Shedding_group = 0/1 matrix of shedding episodes
function CommRisk_func(delta,beta,Seq_clusters,Shedding_cluster,Shedding_group)
  # This function works out fc(t) for one cluster at a time
  function comm_func2(cluster,Shedding_cluster)
    # - Find all cases for relevant cluster
    shedding_c = zeros(Time,Population)
    shedding_c[Shedding_cluster.==cluster].=1;
    pple_c = findall(sum(shedding_c,dims=1).>0);
    # - Period from onset onwards for the relevant cluster cases
    # Pre-define matrix
    Tau=zeros(Time,Population);Tau=Tau .-1;
    for i in 1:length(pple_c)
      # Start of period post onset
      start_i = findall(in(1),diff(shedding_c[:,pple_c[i][2]])) .+1
      # End of period post onset
      if length(start_i)==1
        end_i=Time
      else
        end_i=push!((start_i[2:end].-1),Time)
      end
      # for post onset period, fill in the time since infection
      for j in 1:length(start_i)
        span_a=start_i[j]:end_i[j]
        Tau[span_a,pple_c[i][2]]=span_a.-start_i[j]
      end
    end
    # - Prior to onset, then the function is 0
    xx=exp.(-beta.*Tau);xx[Tau.<0].=0
    # - Time curve for the cluster
    yy=delta.+sum(xx,dims=2)
    return(yy)
  end
  # This function is for weighing the cluster curves per time step
  # (cluster curves have to be weighed so they add up to the group curve)
  function cluster_weigh(t,cluster_funcs,group_func)
    (cluster_funcs[t,:]./sum(cluster_funcs[t,:]))*group_func[t]
  end
  # - Establish unweighted background functions by cluster
  Cluster_func=zeros(Time,length(Seq_clusters))
  for cluster in Seq_clusters
    Cluster_func[:,cluster] = comm_func2(cluster,Shedding_cluster)
  end
  # - Establish background functions by group
  Group_func=comm_func2(1,Shedding_group)
  # - Weigh the cluster curves if clusters != groups
  # RSV A
  if length(Seq_clusters)>1
    Cluster_func_w=zeros(Time,length(Seq_clusters))
    for t in 1:Time
      Cluster_func_w[t,:] = cluster_weigh(t,Cluster_func,Group_func)
    end
  else
    Cluster_func_w=Group_func
  end
  return(Cluster_func_w)
end
# ------------------------------------------------------------------------------
# The likelihood
# ------------------------------------------------------------------------------
# This function takes in theta, a named array of 18 parameters values
function IBM_LL(theta)
  # ___1. Generating background community rates ____
  # RSV A

  Comm_risk_a=CommRisk_func(exp(theta["delta"]),exp(theta["beta"]),1,
                             Shedding_cluster_A,Shedding_a)
  # RSV B
  Comm_risk_b=CommRisk_func(exp(theta["delta"]),exp(theta["beta"]),1,
                             Shedding_cluster_B,Shedding_b)
  # _______________2. Bits that do not change at cluster level _________________
  # Effect of previous homologous and heterologous(group level) infection on susceptibility
  histParA=[0,theta["PrevHom"], theta["PrevHet"],(theta["PrevHom"]+theta["PrevHet"])];
  histMatA = reshape(histParA[B[:]],Time, Population);
  histParB=[0,theta["PrevHet"], theta["PrevHom"],(theta["PrevHom"]+theta["PrevHet"])];
  histMatB = reshape(histParB[B[:]],Time, Population);

  # Effect of current heterologous infection on susceptibility
  currPar=[0,theta["CurrHet"]];
  currMatA = reshape(currPar[Shedding_b[:].+1],Time,Population);
  currMatB = reshape(currPar[Shedding_a[:].+1],Time,Population);

  # Effect of age on susceptibility.
  agePar=[0,theta["age2"],theta["age3"],theta["age4"]];
  ageMat= transpose(reshape(repeat(agePar[A],Time),Population,Time))

  # Combining age and infection history effect on susceptibility
  PI_A = exp.(histMatA .+ currMatA + ageMat);PI_A[Shedding_a.>0].=0;
  PI_B = exp.(histMatB + currMatB + ageMat);PI_B[Shedding_b.>0].=0;

  # Within household transmission coefficient (group level)
  EtaA=exp(theta["etaA"]); EtaB=exp(theta["etaB"])
  # Effect of household size on within household transmission
  hhPar = [1,exp(theta["hhSize"])];
  hhMat = transpose(reshape(repeat(hhPar[H],Time),Population,Time));

  # Effect of viral load and symptoms on infectivity
  V_par = [0,1,exp(theta["LowSym"]),exp(theta["HighSym"])];

  # Community transmission coefficient
  EpsA = exp(theta["epsilonA"]);EpsB = exp(theta["epsilonB"]);

  # Effect of age on community level exposure
  EpsAge = [1,exp(theta["epsAge2"]),exp(theta["epsAge3"])];
  commAge = transpose(reshape(repeat(EpsAge[E],Time),Population,Time));

  # Transmission from close neighbours
  distPar = exp(theta["distRate"]);
  exp_neighbours=exp.(-distPar*Distance_matrix); # exponential distance function

  # _________________3. Given the group level bits, this function ______________
  # Gets the rate or exposure
  function Lambda_c(cluster,group)
    # Predefining some variables
    Load_n_sym = ones(Int64,Time,Population);Shedding= zeros(Int64,Time,Population);

    # Identify the RSV group and select the appropriate variables and data
    # --------------------------
    if group=='A'
      PI = copy(PI_A);Eta=EtaA;Eps=EpsA;
      Load_n_sym[Shedding_cluster_A.==cluster].= Load_n_sym_a[Shedding_cluster_A.==cluster];
      Onset = copy(Onset_a);Onset[Shedding_cluster_A.!=cluster].=0;
      Comm_risk=Comm_risk_a[:,cluster];

      Shedding_cluster=Shedding_cluster_A;
      Shedding[Shedding_cluster_A.==cluster].=1;
    else
      PI = copy(PI_B);Eta=EtaB;Eps=EpsB;
      Load_n_sym[Shedding_cluster_B.==cluster].= Load_n_sym_b[Shedding_cluster_B.==cluster];
      Onset = copy(Onset_b);Onset[Shedding_cluster_B.!=cluster].=0;
      Comm_risk=Comm_risk_b[:,cluster];

      Shedding_cluster=Shedding_cluster_B;
      Shedding[Shedding_cluster_B.==cluster].=1;
    end
    # Within HH rate of exposure
    # --------------------------
    Infectivity= reshape(V_par[Load_n_sym[:]],Time,Population) .* status;
    HHrisk = Infectivity * Cc;

    # Total within household risk, status identifies if an indiviual is present in the household
    in_house= (status.*hhMat).*(Eta*HHrisk);
    # Community rate of exposure
    # --------------------------
    # Time varying component of community exposure/background density funtion
    comm_time=reshape(repeat(Comm_risk,Population),Time,Population);
    # Transmission from close neighbours
    Neighbour_risk=(Infectivity * (Dd.*exp_neighbours))
    # Recalculate
    Neighbour_risk=Neighbour_risk .* status; # RSV *status[t,] so that absent susceptibles are not included
    # Total community risk
    Comm = (Eps*commAge).*(Neighbour_risk .+ comm_time);
    # Total exposure rate (unmodified)
    Lambda=PI.*(in_house .+ Comm);
    return(Lambda)
  end

  # Cluster specific RoE
  Lambda_a = Lambda_c(1,'A');
  Lambda_b= Lambda_c(1,'B');

  # Total RoE
  sum_Lambda_a = Lambda_a;
  sum_Lambda_b = Lambda_b;
  # Gets the likelihood per cluster
  function Like_c(cluster,Lambda,sum_Lambda,group)
    # Identify the RSV group and select the appropriate variables and data
    # --------------------------
    if group=='A'
      Onset = copy(Onset_a);Onset[Shedding_cluster_A.!=cluster].=0;
      atRisk = copy(status);atRisk[Shedding_cluster_A.==cluster].=0;
    else
      Onset = copy(Onset_b);Onset[Shedding_cluster_B.!=cluster].=0;
      atRisk = copy(status);atRisk[Shedding_cluster_B.==cluster].=0;
    end
    # _________________3. Probability of Infection________________________________
    Alpha = (1 .- exp.(-sum_Lambda)) .* (Lambda[:,:,cluster]./sum_Lambda)
    # Because some bits of sum.Lambda will be 0 because of PI being 0, we need
    # do ensure we do not end up with undefined values in Alpha. If PI is 0, so
    # should Alpha
    Alpha[isnan.(Alpha)].=0;
    # _________________4. Probability of Starting Shedding _______________________
    # Function that calculates the probability of onset given latency
    function shedder(t)
      days = t.-Lat_dist.days;
      return(Lat_dist.prob[days.>0]' * Alpha[days[days.>0],:])
    end
    PrShed = zeros(Time, Population);
    for t in 1:Time
      PrShed[t,:]= shedder(t);
    end
    # Likelihood for days of observed data, ie onset and at risk day(ie present and
    # not shedding)
    nPrShed = 1 .- PrShed;
    function logger(p)
      logP= sum(log.(PrShed[(Onset[:,p].==1),p])) + sum(log.(nPrShed[(atRisk[:,p].==1),p]));
      return(logP)
    end
    Like = zeros(Population)
    for p in 1:Population
      Like[p]=logger(p);
    end
    return(sum(Like))
  end

  # Likelihood calculated for all the clusters
  LikeA= Like_c(1,Lambda_a,sum_Lambda_a,'A');
  LikeB= Like_c(1,Lambda_b,sum_Lambda_b,'B');
  # The total likelihood, summed across clusters then across RSV groups
  return(sum(LikeA) + sum(LikeB))
end

# ------------------------------------------------------------------------------
#  The Prior function
# ------------------------------------------------------------------------------
# This function takes in theta, a named array of 18 parameters values
function IBM_prior(theta)
  dnorm(x, μ=0., σ=1.) = logpdf(Normal(μ, σ), x)

  # Factor modifying susceptibility to homologous re-infection
  PreHom = dnorm(theta["PrevHom"], 0, 5)
  # Factor modifying susceptibility to heterologous re-infection
  PreHet = dnorm(theta["PrevHet"], 0, 5)
  # Factor modifying susceptibility to heterologous super-infection
  CurrHet = dnorm(theta["CurrHet"], 0, 5)
  # Factor modifying susceptibility by age group for ages 1-5 years
  age2 = dnorm(theta["age2"], 0, 5)
  # Factor modifying susceptibility by age group for ages 5-15 years
  age3 = dnorm(theta["age3"], 0, 5)
  # Factor modifying susceptibility by age group for ages >15 years
  age4 = dnorm(theta["age4"], 0, 5)
  # Within HH transmission coefficient for RSV A
  etaA =  dnorm(theta["etaA"], -5, 2)
  # Within HH transmission coefficient for RSV B
  etaB =  dnorm(theta["etaB"], -5, 2)
  # Factor modifying within HH transmission by HH size
  hhSize = dnorm(theta["hhSize"], 0, 5)
  # Factor modifying infectiousness by v,load and symptoms for low load symptomatis
  LowSym = dnorm(theta["LowSym"], 0, 5)
  # Factor modifying infectiousness by v,load and symptoms for high load symptomatis
  HighSym = dnorm(theta["HighSym"], 0, 5)
  # Parameter for spatial density kernel: rate of exponential decrease
  distRate = dnorm(theta["distRate"], 0, 5)
  # Community rate transmission coefficient for RSV A
  epsilonA = dnorm(theta["epsilonA"], -5, 2)
  # Community rate transmission coefficient for RSV B
  epsilonB = dnorm(theta["epsilonB"], -5, 2)
  # Factor modifying community transmission by age for ages 1-5
  epsAge2 = dnorm(theta["epsAge2"], 0, 5)
  # Factor modifying community transmission by age for ages >5
  epsAge3 = dnorm(theta["epsAge3"], 0, 5)
  # base rate for the community rate curves: fc(t)=delta + sum(exp(-beta*(t-tau_i)))
  delta = dnorm(theta["delta"], 0, 5)
  # rate of exponential decay for the community rate curves
  beta = dnorm(theta["beta"], -5, 2)
  # the joint prior
  logSum = PreHom + PreHet + CurrHet + age2 + age3 + age4 + etaA + etaB + hhSize +
  + LowSym + HighSym + distRate +  epsilonA + epsilonB + epsAge2 +
  epsAge3 + delta + beta

  return(logSum)
end
# ------------------------------------------------------------------------------
#  The Posterior function
# ------------------------------------------------------------------------------
# This function takes in theta, a named array of 18 parameters values
function IBM_posterior(theta)
  # calculate the log-prior for parameter `theta`
  x = IBM_prior(theta)
  # calculate the log-likelihood for parameter `theta` and
  y = IBM_LL(theta)
  # return the logged posterior probability
  return(x+y)
end
# ------------------------------------------------------------------------------
#  The adaptive MH-MCMC function
# ------------------------------------------------------------------------------
# This function takes in 7 arguments, the last 4 are optional
# target_dist = The posterior function
# init_theta = a named array of 19 parameters values to initiate the MCMC algorithm
# nIter = the number of interations to be run
# proposal_sd = an array of values giving the standard deviation for the joint
#               proposal distribution for the parameters (optional,default is NaN)
# adapt_at = an integer value giving the iteration where the algorithm should
#            start adapting (optional,default is NaN)
# long_runs = Either a true or false value indicating if temporary results of the
#             MCMC algorithm should be saved every (nIter/5) iterations
#            (optional, default is false)
# chain = Where multiple long chains are being run (i.e. temporary results for
#         each are being saved in the same environment) this is an integer value
#         indicating the chain number (optional,default is NaN)
function MH_MCMC(target_dist,init_theta,nIter,proposal_sd=NaN,adapt_at=NaN,long_runs = false,chain=NaN)

  # --- Number of times to print output
  print_every = Tuple(floor(nIter/10):floor(nIter/10):nIter)
  # --- Number of times to save temp files of results (incase of a random crash)
  if long_runs
    save_every=Tuple(floor(nIter/4):floor(nIter/4):nIter)
  end
  # --- Pre-allocation of matrices and random numbers
  # Parameters( keeps track of accepted parameters)
  Params=[transpose(init_theta);zeros(nIter,length(init_theta))]
  # random numbers
  R_acc = log.(rand(Uniform(),(nIter+1)))# for acceptance probabilities when parameters
  if !isnan(adapt_at)
    R_cov=rand(Uniform(),(nIter+1)) # for adapting covariance matrix
  end
  # Log posterior of values in the chain
  target_trace = repeat([target_dist(init_theta)],(nIter+1))
  # initialize values keeping track of number of accepted parameter proposals & acceptance rate
  accepted = 0;acceptance_rate = repeat([0.0],(nIter+1))
  # --- SD and Covariance matrix of proposal distribution( prior to any adaptation)
  if any(isnan.(proposal_sd))
     proposal_sd = collect(init_theta)./10
  end
  covmat_proposal=(0.01/length(init_theta)) * Diagonal(proposal_sd.^2);

  # Per iteration
  for i in 2:(nIter+1)
    # Current parameter values
    theta_current=Params[i-1,:]
    # -------------------------------------------------------------------------
    # 1) UPDATE TRANSMISSION PARAMETERS (simultaneously for all parameters)
    # Adapt the proposal distribution
    if (!isnan(adapt_at) & (i >= adapt_at))
      if R_cov[i]<0.05
        covmat_proposal=(0.01/length(init_theta)) * Diagonal(proposal_sd.^2);
      else
        covmat_proposal= ((2.38^2)/length(init_theta)) * cov(Params[1:i,:])
        covmat_proposal=Matrix(reshape(covmat_proposal,length(init_theta),length(init_theta)))

      end
    end
    # propose a new set of parameters from the proposal distribution
    theta_proposed = NamedArray(reshape(rand(MvNormal(theta_current,covmat_proposal),1),length(init_theta)),(theta_names,))
    # compute the acceptance probability = log((Posterior*/Posterior)(Proposal/Proposal*))
    # but for a symmetric proposal, = log(Posterior*/Posterior)
    target_proposed=target_dist(theta_proposed)
    target_current=target_trace[i-1]
    # accept or reject
    if R_acc[i]< (target_proposed - target_current)
      # update parameter matrix,value of target and acceptance vector of parameters
      Params[i,:]=theta_proposed;
      target_trace[i:(nIter+1)].=target_proposed;
      accepted=accepted + 1;acceptance_rate[i]=accepted/i;
    else
      Params[i,:]=theta_current;acceptance_rate[i]=accepted/i
    end

    # -------------------------------------------------------------------------
    # Print output
    if any(print_every.==i)
      print("      Iteration: ",i,", Acceptance rate = ",acceptance_rate[i], " log target density = ",target_trace[i],"     ")
    end
    # Save intermediate results
    if long_runs
      if any(save_every.==i)
        # Save intermediate results
        writedlm( string("TempRes/ParameterTrace",chain,".csv"),  Params[1:i,:], ',')
        writedlm( string("TempRes/TargetTrace",chain,".csv"),  target_trace[1:i], ',')
        writedlm( string("TempRes/ParaAcceptance",chain,".csv"),  acceptance_rate[1:i], ',')

      end
    end
  end # end of loop going from 2:(nIter+1)
  # Create an object of all the results
  Results = Dict([("Parameters",Params), ("para.acceptance.rate", acceptance_rate),
  ("target.trace",target_trace)]);
  return(Results)
end
