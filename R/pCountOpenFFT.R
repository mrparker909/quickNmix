`%dopar%` <- foreach::`%dopar%`

#' @title Fit FFT N-mixture Model
#' @description Fit an open population N-mixture model using the FFT method of computing the Transition Probability matrix. The four parameters are mean initial site abundance lambda, mean recruitments gamma, survival probability omega, and probability of detection pdet. Parameters can be made to vary over sites and over times by including parameter covariates. Note that this function is essentially a wrapper for optim acting on the nll_FFT function.
#' @param nit Matrix of counts data. Rows represent sites, columns represent sampling occasions. Note that if the data is a vector, then it will be converted to a matrix with a single row.
#' @param K Upper bound on summations in the likelihood function. K should be chosen large enough that the negative log likelihood function is stable (unchanging as K increases). If K=NULL, K=5*max(nit) will be used as default. Default: NULL
#' @param starts Either NULL for default starting values, or a vector of parameter values: `c(log(lambda), log(gamma), logit(omega), logit(pdet))`. Note that the parameter vector will need to be longer by one for each parameter coefficient if covariate values are supplied. The order of coefficients is: `c(lambda, l_s_c, gamma, g_s_c, g_t_c, omega, o_s_c, o_t_c, pdet, p_s_c, p_t_c)`
#' @param l_s_c List of lambda site covariates, Default: NULL
#' @param g_s_c List of gamma site covariates, Default: NULL
#' @param g_t_c List of gamma time covariates, Default: NULL
#' @param o_s_c List of omega site covariates, Default: NULL
#' @param o_t_c List of omega time covariates, Default: NULL
#' @param p_s_c List of pdet site covariates, Default: NULL
#' @param p_t_c List of pdet time covariates, Default: NULL
#' @param VERBOSE If TRUE, will print additional information during model fitting, Default: FALSE
#' @param outfile Location of csv file to write/append parameter values, can be used to checkpoint long running model fits. Default: NULL (no csv file created).
#' @param method Optimization method, passed to optim function, options include: "BFGS", "Nelder-Mead", "CG". Default: "BFGS"
#' @param ... Additional arguments passed to the optimization function optim. For example: `control = list(trace=1, REPORT=1, reltol=1e-10)`
#' @return Returns the fitted model object.
#' @examples 
#' if (interactive()) {
#' # No Covariates
#' nit = matrix(c(1,1,0,2,3), nrow=1) # observations for 1 site, 5 sampling occassions
#' model1 = pCountOpenFFT(nit, K=10)  # fit the model with population upper bound K=10
#' 
#' # Site Covariates
#' o_s_c = list(cov1=c(0,0,1)) # omega site covariates, cov1 is categorical
#' nit = matrix(c(1,1,0,2,3, 
#'                1,0,1,3,2, 
#'                4,1,3,2,0), nrow=3, byrow=T) # 3 sites, 5 sampling occassions
#' model2 = pCountOpenFFT(nit, K=20, o_s_c=o_s_c) # fit the model with population upper bound K=20
#' 
#' # Time Covariates
#' g_t_c = list(temp=c(0.5,0.3,0.6,0.7,NA)) # transition covariates: only first T-1=4 values used 
#' model3 = pCountOpenFFT(nit, K=10, g_t_c=g_t_c)  # fit the model with population upper bound K=10
#' }
#' @importFrom stats optim plogis
#' @rdname pCountOpenFFT
#' @export 
pCountOpenFFT <- function(nit, K=NULL, starts=NULL, 
                          l_s_c=NULL, g_s_c=NULL, g_t_c=NULL, 
                          o_s_c=NULL, o_t_c=NULL, p_s_c=NULL, p_t_c=NULL, 
                          VERBOSE=FALSE, outfile=NULL, method = "BFGS", ...) {
  if(is.null(nrow(nit))) {
    warning("WARNING: converting vector nit to a matrix with one row.")
    nit = matrix(nit, nrow=1)
  }
  
  lamb_starts <- c(1)
  gamm_starts <- c(1)
  omeg_starts <- c(0)
  pdet_starts <- c(0)
  
  lamb_names  <- c("B_l_0")
  gamm_names  <- c("B_g_0")
  omeg_names  <- c("B_o_0")
  pdet_names  <- c("B_p_0")
  
  if(!is.null(l_s_c)) {
    if(!is.list(l_s_c)) {stop("invalid lambda site covariates (l_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")}
    if(any( !unlist(lapply(X = l_s_c, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid lambda site covariates (l_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")
    }
    # update default starting values
    lamb_starts <- rep(1, times=length(l_s_c)+1)
    numCov      <- length(lamb_starts)-1
    lamb_names  <- c(lamb_names, paste0("B_l_s_",1:numCov))
  }
  
  if(!is.null(g_s_c)) {
    if(!is.list(g_s_c)) {stop("invalid gamma site covariates (g_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")}
    if(any( !unlist(lapply(X = g_s_c, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid gamma site covariates (g_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")
    }
    # update default starting values
    gamm_starts <- rep(1, times=length(g_s_c)+1)
    numCov      <- length(gamm_starts)-1
    gamm_names  <- c(gamm_names, paste0("B_g_s_",1:numCov))
  }
  
  if(!is.null(g_t_c)) {
    if(!is.list(g_t_c)) {stop("invalid gamma time covariates (g_t_c) - must be either NULL or a list of vectors, where each vector has length T = ncol(nit), the number of sampling occasions.")}
    if(any( !unlist(lapply(X = g_t_c, FUN = function(X) {is.vector(X) && length(X)==(ncol(nit))})) )) {
      stop("invalid gamma time covariates (g_t_c) - must be either NULL or a list of vectors, where each vector has length T = ncol(nit), the number of sampling occasions.")
    }
    # update default starting values
    gamm_starts <- c(gamm_starts, rep(1, times=length(g_t_c)))
    numCov      <- length(gamm_starts)-length(gamm_names)
    gamm_names  <- c(gamm_names, paste0("B_g_t_",1:numCov))
  }
  
  if(!is.null(o_s_c)) {
    if(!is.list(o_s_c)) {stop("invalid omega site covariates (o_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")}
    if(any( !unlist(lapply(X = o_s_c, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid omega site covariates (o_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")
    }
    # update default starting values
    omeg_starts <- rep(0, times=length(o_s_c)+1)
    omeg_starts[1] <- 0
    numCov      <- length(omeg_starts)-1
    omeg_names  <- c(omeg_names, paste0("B_o_s_",1:numCov))
  }
  
  if(!is.null(o_t_c)) {
    if(!is.list(o_t_c)) {stop("invalid omega time covariates (o_t_c) - must be either NULL or a list of vectors, where each vector has length T = ncol(nit), the number of sampling occasions.")}
    if(any( !unlist(lapply(X = o_t_c, FUN = function(X) {is.vector(X) && length(X)==(ncol(nit))})) )) {
      stop("invalid omega time covariates (o_t_c) - must be either NULL or a list of vectors, where each vector has length T = ncol(nit), the number of sampling occasions.")
    }
    # update default starting values
    omeg_starts <- c(omeg_starts, rep(0, times=length(o_t_c)))
    numCov      <- length(omeg_starts)-length(omeg_names)
    omeg_names  <- c(omeg_names, paste0("B_o_t_",1:numCov))
  }
  
  if(!is.null(p_s_c)) {
    if(!is.list(p_s_c)) {stop("invalid pdet site covariates (p_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")}
    if(any( !unlist(lapply(X = p_s_c, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid pdet site covariates (p_s_c) - must be either NULL or a list of vectors, where each vector has length R = nrow(nit), the number of sites.")
    }
    # update default starting values
    pdet_starts <- rep(0, times=length(p_s_c)+1)
    pdet_starts[1] <- 0
    numCov      <- length(pdet_starts)-1
    pdet_names  <- c(pdet_names, paste0("B_p_s_",1:numCov))
  }
  
  if(!is.null(p_t_c)) {
    if(!is.list(p_t_c)) {stop("invalid pdet time covariates (p_t_c) - must be either NULL or a list of vectors, where each vector has length T = ncol(nit), the number of sampling occasions.")}
    if(any( !unlist(lapply(X = p_t_c, FUN = function(X) {is.vector(X) && length(X)==(ncol(nit))})) )) {
      stop("invalid pdet time covariates (p_t_c) - must be either NULL or a list of vectors, where each vector has length T = ncol(nit), the number of sampling occasions.")
    }
    # update default starting values
    pdet_starts <- c(pdet_starts, rep(0, times=length(p_t_c)))
    numCov      <- length(pdet_starts)-length(pdet_names)
    pdet_names  <- c(pdet_names, paste0("B_p_t_",1:numCov))
  }
  
  if(is.null(starts)) {
    starts <- c(lamb_starts, gamm_starts, omeg_starts, pdet_starts)
  }
  
  NAMES <- c(lamb_names, gamm_names, omeg_names, pdet_names)
  
  if(length(NAMES) != length(starts)) {
    stop(paste0("ERROR: starts must be length: ", length(NAMES), ", not: ", length(starts)))
  }
  
  if(is.null(K)) {
    K = 5*max(nit)
  }
  
  opt = optim(par = starts,
              fn  = nll_FFT, 
              nit = nit,
              K   = K, 
              l_s_c = l_s_c, 
              g_s_c = g_s_c, 
              g_t_c = g_t_c, 
              o_s_c = o_s_c, 
              o_t_c = o_t_c, 
              p_s_c = p_s_c, 
              p_t_c = p_t_c, 
              VERBOSE = VERBOSE, 
              outfile = outfile,
              method  = method,
              hessian = F,
              ...     = ...)
  
  names(opt$par) <- NAMES
  model_data = list(K = K, 
                    nit = nit,  
                    l_s_c = l_s_c, 
                    g_s_c = g_s_c, 
                    g_t_c = g_t_c, 
                    o_s_c = o_s_c, 
                    o_t_c = o_t_c, 
                    p_s_c = p_s_c, 
                    p_t_c = p_t_c)
  
  AIC = 2*opt$value+2*length(opt$par)
  
  param_length = 1
  
  par = opt$par
  # extract coefficients from par, setup covariate matrices (eg:lamb), and coefficient vectors (eg:B_l)
  lamb <- matrix(rep(1,times=nrow(nit)),ncol=1)
  B_l <- par[param_length]
  if(!is.null(l_s_c)) { # Design matrix: rows are sites, cols are covariate values
    lamb <- cbind(lamb, do.call(cbind, l_s_c)) 
    B_l  <- sapply(X = param_length - 1 + 1:(length(l_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_l)
  
  gamm <- matrix(rep(1,times=nrow(nit)),ncol=1)
  B_g <- par[param_length]
  if(!is.null(g_s_c)) { # Design matrix: rows are sites, cols are covariate values
    gamm <- cbind(gamm, do.call(cbind, g_s_c)) 
    B_g  <- sapply(X = param_length - 1 + 1:(length(g_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_g)
  
  gamt <- matrix(rep(0,times=ncol(nit)),ncol=1)
  B_gt <- NULL
  if(!is.null(g_t_c)) { # Design matrix: rows are times, cols are covariate values
    gamt <- do.call(cbind, g_t_c)
    B_gt <- sapply(X = param_length + 1:(length(g_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_gt)
  
  omeg <- matrix(rep(1,times=nrow(nit)),ncol=1)
  B_o <- par[param_length]
  if(!is.null(o_s_c)) { # Design matrix: rows are sites, cols are covariate values
    omeg <- cbind(omeg, do.call(cbind, o_s_c)) 
    B_o  <- sapply(X = param_length - 1 + 1:(length(o_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_o)
  
  omet <- matrix(rep(0,times=ncol(nit)),ncol=1)
  B_ot <- NULL
  if(!is.null(o_t_c)) { # Design matrix: rows are times, cols are covariate values
    omet <- do.call(cbind, o_t_c) 
    B_ot <- sapply(X = param_length + 1:(length(o_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_ot)
  
  pdet <- matrix(rep(1,times=nrow(nit)),ncol=1)
  B_p <- par[param_length]
  if(!is.null(p_s_c)) { # Design matrix: rows are sites, cols are covariate values
    pdet <- cbind(pdet, do.call(cbind, p_s_c)) 
    B_p  <- sapply(X = param_length - 1 + 1:(length(p_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_p)
  
  pdett <- matrix(rep(0,times=ncol(nit)),ncol=1)
  B_pt <- NULL
  if(!is.null(p_t_c)) { # Design matrix: rows are times, cols are covariate values
    pdett <- do.call(cbind, p_t_c) 
    B_pt <- sapply(X = param_length + 1:(length(p_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  
  lambdaM = matrix(0, nrow=nrow(nit), ncol=1)
  for(row in 1:nrow(nit)) {
    lambdaM[row,1] = exp(sum(B_l * lamb[row,], na.rm = T))
  }
  
  gammaM = matrix(0, nrow=nrow(nit), ncol=ncol(nit)-1)
  for(row in 1:nrow(nit)) {
    for(col in 1:(ncol(nit)-1)) {
      gammaM[row,col] = exp(sum(B_g * gamm[row,], na.rm = T)+sum(B_gt * gamt[col,], na.rm = T))
    }
  }
  
  omegaM = matrix(0, nrow=nrow(nit), ncol=ncol(nit)-1)
  for(row in 1:nrow(nit)) {
    for(col in 1:(ncol(nit)-1)) {
      omegaM[row,col] = plogis(sum(B_o * omeg[row,], na.rm = T)+sum(B_ot * omet[col,], na.rm = T))
    }
  }
  
  pdetM = matrix(0, nrow=nrow(nit), ncol=ncol(nit))
  for(row in 1:nrow(nit)) {
    for(col in 1:ncol(nit)) {
      pdetM[row,col] = plogis(sum(B_p * pdet[row,], na.rm = T)+sum(B_pt * pdett[col,], na.rm = T))
    }
  }
  
  estimate_matrices = list(lambda = lambdaM,
                           gamma  = gammaM,
                           omega  = omegaM,
                           pdet   = pdetM)
  
  model_results = list(NLL=opt$value,
                       AIC=AIC,
                       estimate_matrices = estimate_matrices)
  
  model_object = list(optim_results = opt, 
                      model_results = model_results,
                      model_data = model_data)
  
  return(model_object)
}


###########################################
# Likelihood Function:
###########################################
nll_FFT <- function(par, nit, K, l_s_c=NULL, g_s_c=NULL, g_t_c=NULL, o_s_c=NULL, o_t_c=NULL, p_s_c=NULL, p_t_c=NULL, VERBOSE=FALSE, outfile=NULL) {
  if(any(nit > K)) {
    stop("ERROR: K must be larger than any observed data in nit.")
  }
  
  if(is.null(nrow(nit))) {
    warning("WARNING: converting vector nit to a matrix with one row.")
    nit = matrix(nit, nrow=1)
  }
  
  TT <- ncol(nit)
  R <- nrow(nit)
  if(is.null(dim(K)) || !all.equal(dim(K),dim(nit))) {
    K = matrix(K[1], nrow=R, ncol=TT)
  }
  
  param_length = 1
  
  # extract coefficients from par, setup covariate matrices (eg:lamb), and coefficient vectors (eg:B_l)
  lamb <- matrix(rep(1,times=R),ncol=1)
  B_l <- par[param_length]
  if(!is.null(l_s_c)) { # Design matrix: rows are sites, cols are covariate values
    lamb <- cbind(lamb, do.call(cbind, l_s_c)) 
    B_l  <- sapply(X = param_length - 1 + 1:(length(l_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_l)
  
  gamm <- matrix(rep(1,times=R),ncol=1)
  B_g <- par[param_length]
  if(!is.null(g_s_c)) { # Design matrix: rows are sites, cols are covariate values
    gamm <- cbind(gamm, do.call(cbind, g_s_c)) 
    B_g  <- sapply(X = param_length - 1 + 1:(length(g_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_g)
  
  gamt <- matrix(rep(0,times=TT),ncol=1)
  B_gt <- NULL
  if(!is.null(g_t_c)) { # Design matrix: rows are times, cols are covariate values
    gamt <- do.call(cbind, g_t_c)
    B_gt <- sapply(X = param_length + 1:(length(g_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_gt)
  
  omeg <- matrix(rep(1,times=R),ncol=1)
  B_o <- par[param_length]
  if(!is.null(o_s_c)) { # Design matrix: rows are sites, cols are covariate values
    omeg <- cbind(omeg, do.call(cbind, o_s_c)) 
    B_o  <- sapply(X = param_length - 1 + 1:(length(o_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_o)
  
  omet <- matrix(rep(0,times=TT),ncol=1)
  B_ot <- NULL
  if(!is.null(o_t_c)) { # Design matrix: rows are times, cols are covariate values
    omet <- do.call(cbind, o_t_c) 
    B_ot <- sapply(X = param_length + 1:(length(o_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_ot)
  
  pdet <- matrix(rep(1,times=R),ncol=1)
  B_p <- par[param_length]
  if(!is.null(p_s_c)) { # Design matrix: rows are sites, cols are covariate values
    pdet <- cbind(pdet, do.call(cbind, p_s_c)) 
    B_p  <- sapply(X = param_length - 1 + 1:(length(p_s_c)+1), FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_p)
  
  pdett <- matrix(rep(0,times=TT),ncol=1)
  B_pt <- NULL
  if(!is.null(p_t_c)) { # Design matrix: rows are times, cols are covariate values
    pdett <- do.call(cbind, p_t_c) 
    B_pt <- sapply(X = param_length + 1:(length(p_t_c)) - 1, FUN = function(X,par) { par[X] }, par=par)
  }
  param_length = param_length + length(B_pt)
  
  Y <- nit
  g1_t_star <- list()
  g1_t      <- list()
  g1        <- list()
  g2        <- list()
  g_star    <- list()
  g3        <- matrix(-Inf, nrow = K[1]+1, ncol = K[1]+1)
  
  # only need to calculate tp_MAT one time if no covariates for omega and gamma
  if(length(B_g) == 1 && length(B_o) == 1 &&
     length(B_gt) == 0 && length(B_ot) == 0) {
    t_omeg = plogis(sum(B_o*omeg[1,1], na.rm = T))
    t_gamm = exp(sum(B_g*gamm[1,1], na.rm = T))
    g3 <- tp_MAT(M = g3, omeg = t_omeg, gamm= t_gamm)
  }
  
  for(i in 1:R) {
    g1_t_star[[i]] <- numeric(K[1]+1)
    g1_t[[i]]      <- numeric(K[1]+1)
    g_star[[i]]    <- rep(1, K[1]+1)
  }
  
  # apply over sites (1 to R), TODO: this is a prime candidate for parallel since each site i is independent
  #ll_i  <- vapply(X = 1:R, FUN = function(i, K, T, Y, lamb, B_l, pdet, B_p, pt, B_pt, g3, g1_t_star, g1_t,g1,g2, g_star) {
  ll_i = numeric(R)
  for(i in 1:R) {
    ig1_t_star <- g1_t_star[[i]]
    ig1_t      <- g1_t[[i]]
    ig_star    <- g_star[[i]]
    iK         <- K[i]

    # only need to calculate tp_MAT once for each site if no time covariates for omega and gamma
    if((length(B_g) != 1 || length(B_o) != 1) &&
       (length(B_gt) == 0 && length(B_ot) == 0)) {
      t_omeg = plogis(sum(B_o*omeg[i,], na.rm = T))
      t_gamm = exp(sum(B_g*gamm[i,], na.rm = T))
      g3 <- tp_MAT(M = g3, omeg = t_omeg, gamm= t_gamm)
    }
    
    # loop backwards over times t, stopping at t==2
    for(t in TT:2) {
      # calculate tp_MAT for each site and time
      if((length(B_gt) != 0 || length(B_ot) != 0)) {
        t_omeg = plogis(sum(B_o*omeg[i,], na.rm = T) + sum(B_ot*omet[t,], na.rm = T))
        t_gamm = exp(sum(B_g*gamm[i,], na.rm = T) + sum(B_gt*gamt[t,], na.rm = T))
        g3 <- tp_MAT(M = g3, omeg = t_omeg, gamm= t_gamm)
      }
      
      t_pdet = plogis(sum(B_p * pdet[i,], na.rm = T) + sum(B_pt * pdett[t,], na.rm = T))
      ig1_t <- dbinom(x = Y[i,t], size = (0:iK), prob = t_pdet, log=F)
      ig1_t_star <- ig1_t * ig_star
      
      # update g_star
      ig_star = g3 %*% ig1_t_star
    }
    
    # size takes possible values of N (0 to K) at time t==1
    g1 <- dbinom(x = Y[i,1], size = (0:iK), prob = plogis(sum(B_p * pdet[i,], na.rm = T) + sum(B_pt * pdett[1,], na.rm = T)), log=F)
    g2 <- dpois(x = 0:iK, lambda = exp(sum(B_l * lamb[i,], na.rm = T)), log=F)
    
    # apply recursive definition of likelihood
    ll_i[i] = log(sum(g1 * g2 * ig_star))
  }
  
  ll <- sum(ll_i)
  ###########################
  
  if(is.nan(ll)) { ll <- -Inf }
  return(-1*ll)
}


# compute convolution in log space using FFT
#' @importFrom stats convolve
conv_log_FFT <- function(log_x,log_y) {
  m = max(log_x)
  n = max(log_y)
  x = exp(log_x-m)
  y = exp(log_y-n)
  
  conv = convolve(x, rev(y), type = "open")
  conv[which(conv < 0)] = 0
  lconv = log(conv) + m + n
}

tp_MAT <- function(M, omeg, gamm) {
  K <- nrow(M)
  
  M <- foreach(row = 1:K, .combine = "rbind", .export = c("conv_log_FFT")) %dopar% {
    Mrow = M[row,]
    if(row-1==0) { # only recruitments
      Mrow = dpois(0:(K-1),gamm,log=F)
    } else { # recruitments + survivals
      # fast convolve
      Mrow = exp(conv_log_FFT(dbinom(x = 0:(K-1),size = row-1, prob = omeg, log=TRUE), dpois(x = 0:(K-1),lambda = gamm, log=TRUE)))[1:K]
    }
  }
  return(M)
}
