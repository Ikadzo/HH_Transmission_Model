# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 20th September 2019
################################################################################
# This code contains functions needed for running the model that uses genetic
# cluster. It is called by the function 'Cluster_Model_Run.jl'.
# There are 7 functions in this file:
#     1) gen_weight
#     2) CommRisk_func
#     3) ClusterUpdate
#     4) IBM_cluster_LL
#     5) IBM_cluster_prior
#     6) IBM_cluster_posterior
#     7) MH_MCMC_cluster
################################################################################
################################################################################
# This function calculates the genetic weight(a value btwn 0 and 1) associated
# with case.i and case.j at time t. For case pairs that are missing sequences a
# value is assigned by randomly sampling for the pool of genetic distances in
# that cluster. The value therefore changes everytime the function is called.
# This function takes in 5 more arguments in addition to case indices case.i,
# case.j and time t:
# shed_times = a matrix containing the start and end of each shedding episode
#              and the index of the case
# rate_decay = the rate of decay of the exponential function used to model the
#              probability of transmission with increasing genetic distance
# dist_gen_cases = A matrix of the shortest genetic distance between each
#                  shedding episode, if sequence data is available
# dist_cl = distribution of the genetic distances for each genetic cluster
# Shedding_cluster = A matrix of shedding episodes identified by cluster type/id
function gen_weight(case_j,case_i,t,shed_times,rate_decay,dist_gen_cases,dist_cl,
                     Shedding_cluster)
  # - if case.i had a shedding episode
  # Given time t, what episodes are we interested in comparing?
  # Episode(s) from i. look for onsets that are within a 5 day window of t
  o=findall((shed_times[:,2].==case_i) .& (((shed_times[:,1].-t).<=5) .& ((shed_times[:,1].-t).>=0)));
  # Find the relavant episode from j
  w= findall((shed_times[:,2].==case_j) .& ((shed_times[:,1].<=t) .& (t.<=shed_times[:,3])));
  # gen weight between the relevant episodes
  d_gen=dist_gen_cases[o,w];
  # if not available, sample from cluster distances
  if d_gen[1] == -1
    # find cluster
    cl=Shedding_cluster[t,case_j];
    # sample from distances and calculate weight
    d_gen=rand(dist_cl[cl])
  end
  # Calculate the weight
  w_gen = exp(-rate_decay*d_gen)
  # Return the genetic weight
  return(w_gen)
end
# ------------------------------------------------------------------------------
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
# This function updates/changes the cluster ID of the unclustered household
# outbreaks for one household outbreak at a time.
# It takes in 4 arguments:
# Shedding_cluster = A matrix of all shedding episodes identified by cluster type/id
# Shedding_cluster_imputed = A matrix of shedding episodes identified by cluster
#                            type/id for episodes and household outbreaks with some genetic info
# Onset = 0/1 matrix of onset days
# Shedding = 0/1 matrix of shedding episodes
function ClusterUpdate(Shedding_cluster, Shedding_cluster_imputed, Onset,Shedding)
  # Available cluster ids
  cl_id=unique(Shedding_cluster[Shedding_cluster.>0]);
  # All the episodes
  all_epis=findall(Onset.>0);
  # The clustered episodes
  cl_epis=findall((Onset.== 1) .& (Shedding_cluster .> 0))
  # The difference
  x=setdiff(all_epis,cl_epis)
  # If there are unclustered episodes
  if length(x)>0
    # Make a copy otherwise it changes the original matrix being passed even when
    # outside the function environment
    Shedding_cluster_imputed2=copy(Shedding_cluster_imputed)
    # Unclustered shedding patterns
    Shedding_un=copy(Shedding);Shedding_un[Shedding_cluster.>0].= 0;
    # Randomly select a case( and hence a household outbreak) to update
    case= rand(x);
    # HHs
    hh=Demo_data.hhid[case[2]];ppleHH=findall(Demo_data.hhid.==hh)
    # Find all the cases in the same outbreak
    # Shedding times
    times=findall(Tuple(sum(Shedding[:,ppleHH],dims=2) .>0))
    # Time when each HH outbreak started and ended
    startT = times[prepend!((findall(Tuple(diff(times).>5)).+1),1)]
    endT = append!(times[findall(Tuple(diff(times).>5))],times[end])
    # find the appropriate outbreak ( for some reason hh_o does not get assigned
    # if I do not make the process happen within a function)
    function outfind(startT,endT,case)
      hh_o=0
      for o in 1:length(startT)
        if any((startT[o]:endT[o]).== case[1])
          #print(hh_o)
          hh_o = o
        end
      end
      return(hh_o)
    end
    hh_o=outfind(startT,endT,case)
    # Pick at random the cluster to assign to each outbreak that is different from
    # the present cluster
    old_cl=Shedding_cluster_imputed[case[1],case[2]]
    cl_assign= rand(cl_id[cl_id.!=old_cl])
    # Assign the cluster
    timeWin=collect(startT[hh_o]:endT[hh_o])
    Shedding_cluster_imputed2[timeWin,ppleHH]=replace(Shedding_cluster_imputed2[timeWin,ppleHH],old_cl=>cl_assign)

    return(Shedding_cluster_imputed2)
  else
  return(Shedding_cluster)
  end
end
# ------------------------------------------------------------------------------
# The likelihood
# ------------------------------------------------------------------------------
# This function takes in 3 arguments:
# theta = a named array of 19 parameters values
# Shedding_cluster_A = A matrix of all RSV A shedding episodes identified by cluster type/id
# Shedding_cluster_B = A matrix of all RSV B shedding episodes identified by cluster type/id
function IBM_cluster_LL(theta,Shedding_cluster_A,Shedding_cluster_B)
  # ___1. Assigning genetic clusters & generating background community rates ____
  # RSV A
  Seq_clusters_a=Tuple(1:maximum(Shedding_cluster_A))

  Comm_risk_a=CommRisk_func(exp(theta["delta"]),exp(theta["beta"]),Seq_clusters_a,
                             Shedding_cluster_A,Shedding_a)
  # RSV B
  Seq_clusters_b=Tuple(1:maximum(Shedding_cluster_B))
  Comm_risk_b=CommRisk_func(exp(theta["delta"]),exp(theta["beta"]),Seq_clusters_b,
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

      dist_gen_cases=dist_gen_cases_a;
      dist_cl=dist_cl_a;
      Seq_cluster_id = Shedding_cluster_A[cl_pos_a];
      Shedding_cluster=Shedding_cluster_A;
      shed_times=shed_times_a;
      Shedding[Shedding_cluster_A.==cluster].=1;
    else
      PI = copy(PI_B);Eta=EtaB;Eps=EpsB;
      Load_n_sym[Shedding_cluster_B.==cluster].= Load_n_sym_b[Shedding_cluster_B.==cluster];
      Onset = copy(Onset_b);Onset[Shedding_cluster_B.!=cluster].=0;
      Comm_risk=Comm_risk_b[:,cluster];

      dist_gen_cases=dist_gen_cases_b;
      dist_cl=dist_cl_b;
      Seq_cluster_id = Shedding_cluster_B[cl_pos_b];
      Shedding_cluster=Shedding_cluster_B;
      shed_times=shed_times_b;
      Shedding[Shedding_cluster_B.==cluster].=1;
    end
    # Within HH rate of exposure
    # --------------------------
    Infectivity= reshape(V_par[Load_n_sym[:]],Time,Population) .* status;
    HHrisk = Infectivity * Cc;

    # People and onsets in the cluster for which to recalcluate the HH RoE
    pple_in_cluster = findall(Onset.>0);

    # Modifying the HH rate of exposure to take into account within cluster variation
    # This function recalculates HH.risk for the 5 days leading up to an onset in i,
    # if there are infectious housemates
    for i in 1:length(pple_in_cluster)
      # Find onset time and person
      epis=pple_in_cluster[i][1];p=pple_in_cluster[i][2]
      # time period to recalculate risk
      tWind= collect((epis-5):epis);
      # If in this time window there are no infectious housemates, skip this case
      if sum(Shedding[tWind,Cc[p,:].==1])>0
        tWind = tWind[(sum(Shedding[tWind,Cc[p,:].==1],dims=2).>0)[:,1]]
        # For every day of the window that has an infectious housemate, recalculate the risk
        for t in tWind
          shedders = findall((Shedding[t,:].==1) .& (Cc[p,:].==1));
          if length(shedders)>0
            risk=zeros(length(shedders))
            for j in 1:length(shedders)
              risk[j] = (gen_weight(shedders[j],p,t,shed_times,exp(theta["genRate"]),
                                        dist_gen_cases, dist_cl,
                                        Shedding_cluster)[1]) * Infectivity[t,shedders[j]]
            end
            HHrisk[t,p]= sum(risk)
          end
        end
      end
    end
    # Total within household risk, status identifies if an indiviual is present in the household
    in_house= (status.*hhMat).*(Eta*HHrisk);
    # Community rate of exposure
    # --------------------------
    # Time varying component of community exposure/background density funtion
    comm_time=reshape(repeat(Comm_risk,Population),Time,Population);
    # Transmission from close neighbours
    Neighbour_risk=(Infectivity * (Dd.*exp_neighbours))
    # Modifying the comm rate of exposure to take into account within cluster variation
    # This function recalculates Neighbour.risk for the 5 days leading up to an onset in i,
    # if there are infectious neighbors close by
    for i in 1:length(pple_in_cluster)
      # Find onset time and person
      epis=pple_in_cluster[i][1];p=pple_in_cluster[i][2]
      # time period to recalculate risk
      tWind= collect((epis-5):epis);
      # If in this time window there are no infectious neighbour, skip this case
      if sum(Shedding[tWind,Dd[p,:].==1])>0
        tWind = tWind[(sum(Shedding[tWind,Dd[p,:].==1],dims=2).>0)[:,1]]
        # For every day of the window that has an infectious housemate, recalculate the risk
        for t in tWind
          shedders = findall((Shedding[t,:].==1) .& (Dd[p,:].==1));
          if length(shedders)>0
            risk=zeros(length(shedders))
            for j in 1:length(shedders)
              risk[j] = (gen_weight(shedders[j],p,t,shed_times,exp(theta["genRate"]),
                                        dist_gen_cases, dist_cl,
                                        Shedding_cluster)[1] * exp_neighbours[p,shedders[j]]) *
                                        Infectivity[t,shedders[j]] *
                                        exp_neighbours[t,shedders[j]]
            end
            Neighbour_risk[t,p]= sum(risk)
          end
        end
      end
    end
    # Recalculate
    Neighbour_risk=Neighbour_risk .* status; # RSV *status[t,] so that absent susceptibles are not included
    # Total community risk
    Comm = (Eps*commAge).*(Neighbour_risk .+ comm_time);
    # Total exposure rate (unmodified)
    Lambda=PI.*(in_house .+ Comm);
    return(Lambda)
  end

  # Cluster specific RoE
  Lambda_a = zeros(Time,Population,maximum(Shedding_cluster_A));
  for i in 1:maximum(Shedding_cluster_A)
    Lambda_a[:,:,i] = Lambda_c(i,'A');
  end
  Lambda_b= zeros(Time,Population,maximum(Shedding_cluster_B));
  for i in 1:maximum(Shedding_cluster_B)
    Lambda_b[:,:,i] = Lambda_c(i,'B');
  end

  # Total RoE
  sum_Lambda_a = sum(Lambda_a,dims=3);
  sum_Lambda_b = sum(Lambda_b,dims=3);
  # Gets the likelihood per cluster
  function Like_c(cluster,Lambda,sum_Lambda,group)
    # Identify the RSV group and select the appropriate variables and data
    # --------------------------
    if group=='A'
      Onset = copy(Onset_a);Onset[Shedding_cluster_A.!=cluster].=0;
      # If pple are shedding any cluster, then they are not at risk since the
      # clusters are competing for hosts
      atRisk = copy(status);atRisk[Shedding_cluster_A.>0].=0;
    else
      Onset = copy(Onset_b);Onset[Shedding_cluster_B.!=cluster].=0;
      # If pple are shedding any cluster, then they are not at risk since the
      # clusters are competing for hosts
      atRisk = copy(status);atRisk[Shedding_cluster_B.>0].=0;
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
  LikeA= zeros(maximum(Seq_clusters_a))
  for c in 1:maximum(Seq_clusters_a)
    LikeA[c]=Like_c(c,Lambda_a,sum_Lambda_a,'A');
  end

  LikeB= zeros(maximum(Seq_clusters_b))
  for c in 1:maximum(Seq_clusters_b)
    LikeB[c]=Like_c(c,Lambda_b,sum_Lambda_b,'B');
  end

  # The total likelihood, summed across clusters then across RSV groups
  return(sum(LikeA) + sum(LikeB))
end

# ------------------------------------------------------------------------------
#  The Prior function
# ------------------------------------------------------------------------------
# This function takes in theta, a named array of 19 parameters values
function IBM_cluster_prior(theta)
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
  # Parameter for genetic distance kernel: rate of exponential decrease
  genRate = dnorm(theta["genRate"], 0, 5)
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
  + LowSym + HighSym + distRate + genRate + epsilonA + epsilonB + epsAge2 +
  epsAge3 + delta + beta

  return(logSum)
end
# ------------------------------------------------------------------------------
#  The Posterior function
# ------------------------------------------------------------------------------
# This function takes in the same 3 arguments as the likelihood function
function IBM_cluster_posterior(theta,Shedding_cluster_A,Shedding_cluster_B)
  # calculate the log-prior for parameter `theta`
  x = IBM_cluster_prior(theta)
  # calculate the log-likelihood for parameter `theta` and
  y = IBM_cluster_LL(theta,Shedding_cluster_A,Shedding_cluster_B)
  # return the logged posterior probability
  return(x+y)
end
# ------------------------------------------------------------------------------
#  The adaptive Reversible Jump MH-MCMC function
# ------------------------------------------------------------------------------
# This function takes in 11 arguments, the last 4 are optional
# Shedding_cluster_A = A matrix of all RSV A shedding episodes identified by cluster type/id
# Shedding_cluster_B = A matrix of all RSV B shedding episodes identified by cluster type/id
# target_dist = The posterior function
# init_theta = a named array of 19 parameters values to initiate the MCMC algorithm
# init_rsva_outbreaks = Initial number of RSV A HH outbreaks in each cluster
# init_rsvb_outbreaks = Initial number of RSV B HH outbreaks in each cluster
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
function MH_MCMC_cluster(Shedding_cluster_A,Shedding_cluster_B,target_dist,init_theta,init_rsva_outbreaks,init_rsvb_outbreaks,
                         nIter,proposal_sd=NaN,adapt_at=NaN,long_runs = false,chain=NaN)

  # --- Number of times to print output
  print_every = Tuple(floor(nIter/10):floor(nIter/10):nIter)
  # --- Number of times to save temp files of results (incase of a random crash)
  if long_runs
    save_every=Tuple(floor(nIter/5):floor(nIter/5):nIter)
  end
  # --- Pre-allocation of matrices and random numbers
  # Parameters( keeps track of accepted parameters)
  Params=[transpose(init_theta);zeros(nIter,length(init_theta))]
  # random numbers
  R_acc = log.(rand(Uniform(),(nIter+1)))# for acceptance probabilities when parameters
  R_accA = log.(rand(Uniform(),(nIter+1)))# for acceptance probabilities when updating rsv a data
  R_accB = log.(rand(Uniform(),(nIter+1)))# for acceptance probabilities when updating rsv b data
  if !isnan(adapt_at)
    R_cov=rand(Uniform(),(nIter+1)) # for adapting covariance matrix
  end
  # Log posterior of values in the chain
  target_trace = repeat([target_dist(init_theta,Shedding_cluster_A,Shedding_cluster_B)],(nIter+1))
  # initialize values keeping track of number of accepted parameter proposals & acceptance rate
  accepted = 0;acceptance_rate = repeat([0.0],(nIter+1))
  # Keeping track of accepted chages to the augmented data
  acceptedA = 0;acceptance_rateA = repeat([0.0],(nIter+1))
  acceptedB = 0;acceptance_rateB = repeat([0.0],(nIter+1))
  # Keep track of the number of HH outbreaks in each cluster
  RSVA_outbreaks=transpose(reshape(repeat(init_rsva_outbreaks,(nIter+1)),length(init_rsva_outbreaks),(nIter+1)));
  RSVB_outbreaks=transpose(reshape(repeat(init_rsvb_outbreaks,(nIter+1)),length(init_rsvb_outbreaks),(nIter+1)));
  # Keep track of the cluster ids for the uniformed outbreaks
  Uninformed_outbreak_a= transpose(reshape(repeat(Shedding_cluster_A[one_case_a],(nIter+1)),length(one_case_a),(nIter+1)));
  Uninformed_outbreak_b= transpose(reshape(repeat(Shedding_cluster_B[one_case_b],(nIter+1)),length(one_case_b),(nIter+1)));
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
        covmat_proposal=Matrix(reshape(covmat_proposal,19,19))

      end
    end
    # propose a new set of parameters from the proposal distribution
    theta_proposed = NamedArray(reshape(rand(MvNormal(theta_current,covmat_proposal),1),19),(theta_names,))
    # compute the acceptance probability = log((Posterior*/Posterior)(Proposal/Proposal*))
    # but for a symmetric proposal, = log(Posterior*/Posterior)
    target_proposed=target_dist(theta_proposed,Shedding_cluster_A,Shedding_cluster_B)
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
    # 2) UPDATE CLUSTER ID (for a single outbreak at a time, sequentially for RSV A and B)
    # ----- RSV A
    # Proposed changes
    Shedding_cluster_A_prop = ClusterUpdate(Shedding_cluster_fixedA,Shedding_cluster_A,
                                            Onset_a,Shedding_a)
    # Difference between the proposed and current
    # Counts the frequency of each cluster
    function cl_count(Shedding_cluster,cl)
      return(sum(Shedding_cluster[:].==cl))
    end
    diff_cl=repeat([0],maximum(Shedding_cluster_A))
    for k in 1:maximum(Shedding_cluster_A)
      diff_cl[k]=cl_count(Shedding_cluster_A_prop,k)-cl_count(Shedding_cluster_A,k)
    end

    # which clusters has lost an outbreak
    Ck = findmin(diff_cl)[2];
    # which cluster has gained an outbreaks
    Cm = findmax(diff_cl)[2];
    # Proposal ratio for the change in cluster confirguration
    prob_aug = RSVA_outbreaks[i-1,Ck]/(RSVA_outbreaks[i-1,Cm]+1)
    # Acceptance probability
    target_proposed=target_dist(Params[i,:],Shedding_cluster_A_prop,Shedding_cluster_B)
    target_current=target_trace[i]
    r = (target_proposed - target_current) .+ log(prob_aug);
    # accept or reject
    if R_accA[i]<r
      # Update matrix keeping track of number on outbreak in each cluster,
      # value of target, acceptance vector and the clustering pattern
      RSVA_outbreaks[i,[Ck,Cm]] = RSVA_outbreaks[i-1,[Ck,Cm]].+[-1,1];
      target_trace[i:(nIter+1)] .= target_proposed;
      acceptedA = acceptedA + 1;acceptance_rateA[i]=acceptedA/i;
      Shedding_cluster_A = Shedding_cluster_A_prop;
    else
      RSVA_outbreaks[i,[Ck,Cm]] = RSVA_outbreaks[i-1,[Ck,Cm]]
      acceptance_rateA[i] = acceptedA/i
    end
    Uninformed_outbreak_a[i,:] = Shedding_cluster_A[one_case_a]
    # ----- RSV B
    # Proposed changes
    Shedding_cluster_B_prop = ClusterUpdate(Shedding_cluster_fixedB,Shedding_cluster_B,
                                            Onset_b,Shedding_b)
    # Difference between the proposed and current
    xx=countmap(Shedding_cluster_B_prop[:])
    yy=countmap(Shedding_cluster_B[:])
    diff_cl=repeat([0],maximum(Shedding_cluster_B))
    for k in 1:maximum(Shedding_cluster_B)
      diff_cl[k]=xx[k]-yy[k]
    end
    # which clusters has lost an outbreak
    Ck = findmin(diff_cl)[2];
    # which cluster has gained an outbreaks
    Cm = findmax(diff_cl)[2];
    # Proposal ratio for the change in cluster confirguration
    prob_aug = RSVB_outbreaks[i-1,Ck]/(RSVB_outbreaks[i-1,Cm]+1)
    # Acceptance probability
    target_proposed=target_dist(Params[i,:],Shedding_cluster_A,Shedding_cluster_B_prop)
    target_current=target_trace[i]
    r = (target_proposed - target_current) .+ log(prob_aug);
    # accept or reject
    if R_accB[i]<r
      # Update matrix keeping track of number on outbreak in each cluster,
      # value of target, acceptance vector and the clustering pattern
      RSVB_outbreaks[i,[Ck,Cm]] = RSVB_outbreaks[i-1,[Ck,Cm]].+[-1,1];
      target_trace[i:(nIter+1)] .= target_proposed;
      acceptedB = acceptedB + 1;acceptance_rateB[i]=acceptedB/i;
      Shedding_cluster_B = Shedding_cluster_B_prop;
    else
      RSVB_outbreaks[i,[Ck,Cm]] = RSVB_outbreaks[i-1,[Ck,Cm]]
      acceptance_rateB[i] = acceptedB/i
    end
    Uninformed_outbreak_b[i,:] = Shedding_cluster_B[one_case_b]

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
        writedlm( string("TempRes/Uninformed_outbreak_a",chain,".csv"),  Uninformed_outbreak_a[1:i,:], ',')
        writedlm( string("TempRes/Uninformed_outbreak_b",chain,".csv"),  Uninformed_outbreak_b[1:i,:], ',')
        writedlm( string("TempRes/AugAccA",chain,".csv"),  acceptance_rateA[1:i], ',')
        writedlm( string("TempRes/AugAccB",chain,".csv"),  acceptance_rateB[1:i], ',')
      end
    end
  end # end of loop going from 2:(nIter+1)
  # Create an object of all the results
  Results = Dict([("Parameters",Params), ("para.acceptance.rate", acceptance_rate),
  ("target.trace",target_trace),("rsva.outbreaks",RSVA_outbreaks),("aug.a.acc",acceptance_rateA),
  ("rsvb.outbreaks",RSVB_outbreaks),("aug.b.acc",acceptance_rateB),("uninf.a",Uninformed_outbreak_a),
  ("uninf.b",Uninformed_outbreak_b)]);
  return(Results)
end
