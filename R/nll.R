`%dopar%` <- foreach::`%dopar%`
`%do%`    <- foreach::`%do%`

#' @title Negative Log Likelihood
#' @description Computes the negative log likelihood function for the open population asymptotic N-mixtures model.
#' @param par Vector of parameter values: c(log(lambda), log(gamma), logit(omega), logit(pdet)). Note that the parameter vector will need to be longer if covariate values are supplied.
#' @param nit Matrix of counts data. Rows represent sites, columns represent sampling occasions. Note that if the data is a vector, then it will be converted to a matrix with a single row.
#' @param K Upper bound on summations in the likelihood function. K should be chosen large enough that the negative log likelihood function is stable (unchanging as K increases).
#' @param l_s_c List of lambda site covariates, Default: NULL
#' @param g_s_c List of gamma site covariates, Default: NULL
#' @param g_t_c List of gamma time covariates, Default: NULL
#' @param o_s_c List of omega site covariates, Default: NULL
#' @param o_t_c List of omega time covariates, Default: NULL
#' @param p_s_c List of pdet site covariates, Default: NULL
#' @param p_t_c List of pdet time covariates, Default: NULL
#' @param SMALL_a_CORRECTION If TRUE will apply the small a correction when calculating the transition probability matrix, Default: FALSE
#' @param VERBOSE If TRUE, will print additional information, Default: FALSE
#' @param outfile Location of csv file to write/append parameter values, Default: NULL
#' @return Returns the negative log likelihood function evaluated at par.
#' @details DETAILS
#' @importFrom stats plogis dbinom dnorm dpois
#' @importFrom utils read.csv write.table
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname nll
#' @export 
nll <- function(par, nit, K, l_s_c=NULL, g_s_c=NULL, g_t_c=NULL, o_s_c=NULL, o_t_c=NULL, p_s_c=NULL, p_t_c=NULL, SMALL_a_CORRECTION=FALSE, VERBOSE=FALSE, outfile=NULL) {
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
  
  if(length(B_g) == 1 && length(B_o) == 1 &&
     length(B_gt) == 0 && length(B_ot) == 0) {
    t_omeg = plogis(sum(B_o*omeg[1,1]))
    t_gamm = exp(sum(B_g*gamm[1,1]))
    g3 <- log_tp_MAT_lse(M = g3, omeg = t_omeg, gamm= t_gamm, corrections = SMALL_a_CORRECTION)
  }
  
  for(i in 1:R) {
    g1_t_star[[i]] <- numeric(K[1]+1)
    g1_t[[i]]      <- numeric(K[1]+1)
    g_star[[i]]    <- rep(0, K[1]+1)
  }
  
  ll_i = numeric(R)
  for(i in 1:R) {
    ig1_t_star <- g1_t_star[[i]]
    ig1_t      <- g1_t[[i]]
    ig_star    <- g_star[[i]]
    iK         <- K[i]
    
    if((length(B_g) != 1 || length(B_o) != 1) &&
       (length(B_gt) == 0 && length(B_ot) == 0)) {
      t_omeg = plogis(sum(B_o*omeg[i,]))
      t_gamm = exp(sum(B_g*gamm[i,]))
      g3 <- log_tp_MAT_lse(M = g3, omeg = t_omeg, gamm= t_gamm, corrections = SMALL_a_CORRECTION)
    }
    
    # loop backwards over times t, stopping at t==2
    for(t in TT:2) {
      if((length(B_gt) != 0 || length(B_ot) != 0)) {
        t_omeg = plogis(sum(B_o*omeg[i,]) + sum(B_ot*omet[t,]))
        t_gamm = exp(sum(B_g*gamm[i,]) + sum(B_gt*gamt[t,]))
        g3 <- log_tp_MAT_lse(M = g3, omeg = t_omeg, gamm= t_gamm, corrections = SMALL_a_CORRECTION)
      }
      
      t_pdet = plogis(sum(B_p * pdet[i,]) + sum(B_pt * pdett[t,]))
      
      # size takes possible value of N (0 to K) at time t (for t > 1)
      ig1_t <- dbinom(x = Y[i,t], size = (0:iK), prob = t_pdet, log=T)
      ig1_t_star <- ig1_t + ig_star
      
      # update g_star
      ig_star = Ax_log(g3,ig1_t_star)
    }
    
    # size takes possible values of N (0 to K) at time t==1
    g1 <- dbinom(x = Y[i,1], size = (0:iK), prob = plogis(sum(B_p * pdet[i,]) + sum(B_pt * pdett[1,])), log=T)
    
    # probability for initial population sizes given lambda
    g2 <- dpois(x = 0:iK, lambda = exp(sum(B_l * lamb[i,])), log=T)
    g2[2:Y[i,1]-1] = -Inf # initial population cannot be less than observed initial population
    
    # apply recursive definition of likelihood
    ll_i[i] = logSumExp(g1 + g2 + ig_star)
  }
  
  ll <- sum(ll_i)
  if(is.nan(ll)) { ll <- -Inf }
  
  if(VERBOSE) {
    print(paste0("nll=",-1*ll))
    print("par=")
    print(par)
  }
  
  if(!is.null(outfile)) {
    parnames = c(paste0("par", 1:length(par)), "nll")
    parstate = c(par, -1*ll)
    names(parstate) = parnames
    datcsv = data.frame(t(parstate))
    if(file.exists(outfile)) {
      checknames = colnames(read.csv(outfile))
      if(!all(names(parstate)==checknames)) { warning(paste("oufile column names do not match new column names: ", outfile)) }
      write.table(datcsv, file = outfile, append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
    } else {
      write.table(datcsv, file = outfile, col.names = TRUE, row.names = FALSE, sep = ",")
    }
  }
  
  return(-1*ll)
}

# numerically stable transition probability matrix
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param M PARAM_DESCRIPTION
#' @param omeg PARAM_DESCRIPTION
#' @param gamm PARAM_DESCRIPTION
#' @param corrections PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[foreach]{foreach}}
#' @rdname log_tp_MAT_lse
#' @importFrom foreach foreach
log_tp_MAT_lse <- function(M, omeg, gamm, corrections) {
  K <- nrow(M)
  
  M <- foreach::foreach(row = 1:K, .combine = "rbind",.export = c("Pab", "Pab_gamma", "Pab_omega", "Pab_asymptotic", "logSumExp", "logSubtractExp")) %dopar% {
    Mrow = M[row,]
    b=1
    while(b < K+1) {
      Mrow[b] = Pab_asymptotic(row-1,b-1,omeg,gamm,corrections)
      b = b + 1
    }
    return(Mrow)
  }
  
  return(M)
}

# asymptotic transition probability
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @param omega PARAM_DESCRIPTION
#' @param gamma PARAM_DESCRIPTION
#' @param corrections PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Pab_asymptotic
#' @importFrom stats pnorm
Pab_asymptotic = function(a,b,omega,gamma, corrections = FALSE) {
  if(corrections) {
    # small a correction
    if(a < 0.5*(b-0.25*gamma) & a < 0.5*(1.75*gamma-b)) {
      return(Pab(a,b,omega,gamma))
    }
  }
  
  k = min(a,b)
  
  s1 = a*omega*(1-omega)
  s2 = gamma
  
  # no survival term!
  if(s1==0) { return(Pab_gamma(a,b,gamma))}
  
  # no recruitment term!
  if(s2==0) { return(Pab_omega(a,b,omega))}
  
  mu1 = a*omega
  mu2 = gamma
  
  mstar = (mu1*s2 + (b-mu2)*s1)/(s1+s2)
  sstar = s1*s2/(s1+s2) 
  
  L = 0.5*log(sstar/(2*pi*s1*s2)) - 0.5*(mu1+mu2-b)^2/(s1+s2)
  
  L + logSubtractExp(c(pnorm(q = k, mean = mstar-0.5, sd = sqrt(sstar), log.p=T),
                       pnorm(q = 0, mean = mstar, sd = sqrt(sstar), log.p=T)))
}

# numerically stable log space addition
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname logSumExp
logSumExp <- function(x) {
  if(all(is.infinite(x))) { return(x[1]) }
  x = x[which(is.finite(x))]
  ans = x[1]
  for(i in seq_along(x)[-1]) {
    ma = max(ans,x[i])
    mi = min(ans,x[i])
    ans = ma + log1p(exp(mi-ma))
  }
  return(ans)
}

# numerically stable log space subtraction
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname logSubtractExp
logSubtractExp <- function(x) {
  ma = max(x)
  mi = min(x)
  ma + log1p(-exp(mi-ma))
}

# Multiply a matrix times a vector y=Ax in log space
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param logA PARAM_DESCRIPTION
#' @param logx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Ax_log
Ax_log <- function(logA,logx) {
  xl = length(logx)
  if(is.null(xl)) { xl = nrow(logx) }
  Ac = ncol(logA)
  if(Ac != xl) {stop("columns of A must match length of x")}
  
  y = matrix(-Inf, nrow=nrow(logA), ncol=1)
  for(r in 1:nrow(logA)) {
    temp = logA[r,] + logx
    y[r,1] = logSumExp(temp)
  }
  return(y)
}

# no survival term, omega=0 or a=0
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @param gamma PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Pab_gamma
#' @importFrom stats pnorm
Pab_gamma = function(a, b, gamma) {
  if(gamma <= 0) { return(-Inf) }
  if(b==0) { return(dnorm(x = b, mean = gamma-0.5, sd = sqrt(gamma), log = T)) } # no summation
  logSubtractExp(c(pnorm(q = b, mean = gamma-0.5, sd = sqrt(gamma), log.p=T), 
                   pnorm(q = 0, mean = gamma-0.5, sd = sqrt(gamma), log.p=T)))
}

# no recruitment term, gamma=0 implies b <= a
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @param omega PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Pab_omega
#' @importFrom stats pnorm
Pab_omega = function(a, b, omega) {
  if(omega*a*(1-omega) <= 0) { ifelse(b==0, return(0), return(-Inf)) }
  if(b > a) { return(-Inf) }
  if(b==0) { return(dnorm(x = 0, mean = omega*a-0.5, sd = sqrt(omega*a*(1-omega)), log = T)) }
  logSubtractExp(c(pnorm(q = b, mean = omega*a-0.5, sd = sqrt(omega*a*(1-omega)), log.p=T),
                   pnorm(q = 0, mean = omega*a-0.5, sd = sqrt(omega*a*(1-omega)), log.p=T)))
}

# traditional Pab in log space:
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @param omega PARAM_DESCRIPTION
#' @param gamma PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Pab
Pab <- function(a,b,omega,gamma) {
  s = log(0)
  for(c in 0:min(a,b)) {
    s = logSumExp(c(s, dbinom(x = c,size = a, prob = omega,log=T)+dpois(x = b-c,lambda = gamma, log=T)))
  }
  
  return(s)
}
