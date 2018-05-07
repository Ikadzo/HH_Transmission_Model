# fitR will not load on the latest version or R
# so I got these from github...need to code my own sampler!!!
# 
# 24th May 2017
library('tmvtnorm')
updateCovmat <- function(covmat,theta.mean,theta,i) {
  
  if(is.null(names(theta))){
    stop("Argument ",sQuote("theta")," must be named.",.call=FALSE)
  }
  if(is.null(names(theta.mean))){
    stop("Argument ",sQuote("theta.mean")," must be named.",.call=FALSE)
  }
  if(is.null(rownames(covmat))){
    stop("Argument ",sQuote("covmat")," must have named rows.",.call=FALSE)
  }
  if(is.null(colnames(covmat))){
    stop("Argument ",sQuote("covmat")," must have named columns.",.call=FALSE)
  }
  
  covmat <- covmat[names(theta),names(theta)]
  theta.mean <- theta.mean[names(theta)]
  
  residual <- as.vector(theta-theta.mean)
  covmat <- (covmat*(i-1)+(i-1)/i*residual%*%t(residual))/i
  theta.mean <- theta.mean + residual/i
  
  return(list(covmat=covmat,theta.mean=theta.mean))
}
printNamedVector <- function(x, fmt="%.2f", sep=" | ") {
  
  paste(paste(names(x),sprintf(fmt,x),sep=" = "),collapse=sep)
  
}
mcmcMH <- function(target, init.theta, proposal.sd = NULL,
                   n.iterations, covmat = NULL,
                   limits=list(lower = NULL, upper = NULL),
                   adapt.size.start = NULL, adapt.size.cooling = 0.99,
                   adapt.shape.start = NULL,
                   print.info.every = n.iterations/100,
                   verbose = FALSE, max.scaling.sd = 50,
                   acceptance.rate.weight = NULL,
                   acceptance.window = NULL) {
  
  # initialise theta
  theta.current <- init.theta
  theta.propose <- init.theta
  
  # extract theta of gaussian proposal
  covmat.proposal <- covmat
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  
  # reorder vector and matrix by names, set to default if necessary
  theta.names <- names(init.theta)
  if (!is.null(proposal.sd) && is.null(names(proposal.sd))) {
    names(proposal.sd) <- theta.names
  }
  
  if (is.null(covmat.proposal)) {
    if (is.null(proposal.sd)) {
      proposal.sd <- init.theta/10
    }
    covmat.proposal <-
      matrix(diag(proposal.sd[theta.names]^2, nrow = length(theta.names)),
             nrow = length(theta.names),
             dimnames = list(theta.names, theta.names))
  } else {
    covmat.proposal <- covmat.proposal[theta.names,theta.names]
  }
  
  if (is.null(lower.proposal)) {
    lower.proposal <- init.theta
    lower.proposal[] <- -Inf
  } else {
    lower.proposal <- lower.proposal[theta.names]
  }
  
  if (is.null(upper.proposal)) {
    upper.proposal <- init.theta
    upper.proposal[] <- Inf
  } else {
    upper.proposal <- upper.proposal[theta.names]
  }
  
  # covmat init
  covmat.proposal.init <- covmat.proposal
  
  adapting.size <- FALSE # will be set to TRUE once we start
  # adapting the size
  
  adapting.shape <- FALSE # will be set to TRUE once we start
  # adapting the shape
  
  # find estimated theta
  theta.estimated.names <- names(which(diag(covmat.proposal) > 0))
  
  # evaluate target at theta init
  target.theta.current <- target(theta.current)
  
  # if return value is a vector, set log.density and trace
  if (class(target.theta.current) == "numeric") {
    target.theta.current <-
      list(log.density = target.theta.current,
           trace = theta.current)
  }
  
  if (!is.null(print.info.every)) {
    message("Init: ", printNamedVector(theta.current[theta.estimated.names]),
            ", target: ", target.theta.current[["log.density"]])
  }
  
  trace <- data.frame(t(c(target.theta.current[["trace"]],
                          target.theta.current["log.density"])))
  
  # acceptance rate
  acceptance.rate <- 0
  if (!is.null(acceptance.window)) {
    acceptance.series <- c()
  }
  
  # scaling factor for covmat size
  scaling.sd  <- 1
  
  # scaling multiplier
  scaling.multiplier <- 1
  
  # empirical covariance matrix (0 everywhere initially)
  covmat.empirical <- covmat.proposal
  covmat.empirical[,] <- 0
  
  # empirical mean vector
  theta.mean <- theta.current
  
  # if print.info.every is null never print info
  if (is.null(print.info.every)) {
    print.info.every <- n.iterations + 1
  }
  
  start_iteration_time <- Sys.time()
  
  for (i.iteration in seq_len(n.iterations)) {
    
    # adaptive step
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
        (is.null(adapt.shape.start) || acceptance.rate*i.iteration < adapt.shape.start)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      # adapt size of covmat until we get enough accepted jumps
      scaling.multiplier <- exp(adapt.size.cooling^(i.iteration-adapt.size.start) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd,max.scaling.sd))
      # only scale if it doesn't reduce the covariance matrix to 0
      covmat.proposal.new <- scaling.sd^2*covmat.proposal.init
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] <
                .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      }
      
    } else if (!is.null(adapt.shape.start) &&
               acceptance.rate*i.iteration >= adapt.shape.start) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        # flush.console()
        adapting.shape <- TRUE
      }
      # adapt shape of covmat using optimal scaling factor for multivariate target distributions
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    }
    
    # print info
    if (i.iteration %% ceiling(print.info.every) == 0) {
      ## end_iteration_time <- Sys.time()
      state.mcmc <- target.theta.current$trace
      ## suppressMessages(time.estimation <- round(as.period((end_iteration_time-start_iteration_time)*10000/round(print.info.every))))
      ## message("Iteration: ",i.iteration,"/",n.iterations,", ETA: ",time.estimation,", acceptance rate: ",sprintf("%.3f",acceptance.rate),appendLF=FALSE)
      message("Iteration: ",i.iteration,"/", n.iterations,
              ", acceptance rate: ",
              sprintf("%.3f",acceptance.rate), appendLF=FALSE)
      if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd),
                ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier),
                appendLF=FALSE)
      }
      message(", state: ",printNamedVector(state.mcmc))
      message(", logdensity: ", target.theta.current$log.density)
      ## start_iteration_time <- end_iteration_time
    }
    
    # propose another parameter set
    if (any(diag(covmat.proposal)[theta.estimated.names] <
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names,theta.estimated.names])
      stop("non-positive definite covmat",call.=FALSE)
    }
    if (length(theta.estimated.names) > 0) {
      theta.propose[theta.estimated.names] <-
        as.vector(rtmvnorm(1,
                           mean =
                             theta.current[theta.estimated.names],
                           sigma =
                             covmat.proposal[theta.estimated.names,theta.estimated.names],
                           lower =
                             lower.proposal[theta.estimated.names],
                           upper = upper.proposal[theta.estimated.names]))
    }
    
    # evaluate posterior of proposed parameter
    target.theta.propose <- target(theta.propose)
    # if return value is a vector, set log.density and trace
    if (class(target.theta.propose) == "numeric") {
      target.theta.propose <-
        list(log.density = target.theta.propose,
             trace = theta.propose)
    }
    
    if (!is.finite(target.theta.propose$log.density)) {
      # if posterior is 0 then do not compute anything else and don't accept
      log.acceptance <- -Inf
      
    }else{
      
      # compute Metropolis-Hastings ratio (acceptance probability)
      log.acceptance <- target.theta.propose$log.density -
        target.theta.current$log.density
      log.acceptance <- log.acceptance +
        dtmvnorm(x = theta.current[theta.estimated.names],
                 mean =
                   theta.propose[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      log.acceptance <- log.acceptance -
        dtmvnorm(x = theta.propose[theta.estimated.names],
                 mean = theta.current[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      
    }
    
    if (verbose) {
      message("Propose: ", theta.propose[theta.estimated.names],
              ", target: ", target.theta.propose[["log.density"]],
              ", acc prob: ", exp(log.acceptance), ", ",
              appendLF = FALSE)
    }
    
    if (is.accepted <- (log(runif (1)) < log.acceptance)) {
      # accept proposed parameter set
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      if (verbose) {
        message("accepted")
      }
    } else if (verbose) {
      message("rejected")
    }
    trace <- rbind(trace,c(target.theta.current[["trace"]],
                           target.theta.current["log.density"]))
    
    # update acceptance rate
    if (i.iteration == 1) {
      acceptance.rate <- is.accepted
    } else {
      if (is.null(acceptance.rate.weight)) {
        if (is.null(acceptance.window)) {
          acceptance.rate <- acceptance.rate +
            (is.accepted - acceptance.rate) / i.iteration
        } else {
          acceptance.series <- c(is.accepted, acceptance.series)
          if (length(acceptance.series) > acceptance.window) {
            acceptance.series <- acceptance.series[-length(acceptance.series)]
          }
          acceptance.rate <- mean(acceptance.series)
        }
      } else {
        acceptance.rate <- acceptance.rate * (1 - acceptance.rate.weight) +
          is.accepted * acceptance.rate.weight
      }
    }
    
    # update empirical covariance matrix
    tmp <- updateCovmat(covmat.empirical, theta.mean,
                        theta.current, i.iteration)
    covmat.empirical <- tmp$covmat
    theta.mean <- tmp$theta.mean
    
  }
  
  return(list(trace = trace,
              acceptance.rate = acceptance.rate,
              covmat.empirical = covmat.empirical))
}
