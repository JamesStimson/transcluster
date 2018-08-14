##################################
#
# James Stimson 05/03/18
#
# uses integrate function from stats
#
##################################

ProbHeight = function(h, N, lambda, dd, dterm_list){
  psum <- 0.0

  for (i in seq(0, N))
  {
    hterm = lambda^{i+1} * h^i * exp(-lambda*h)/factorial(i)
    dterm = dterm_list[[i+1]]
    psum <- psum + hterm*dterm
  }
  return(psum)
}

# Gamma version of ProbHeight() -TO DO: Wire this up-
ProbHeightGamma = function(h, N, a, b, dd, dterm_list){
  psum <- 0.0

  for (i in seq(0, N))
  {
    hterm = b^{a*(i+1)} * h^{a*(i+1)-1} * exp(-b*h)/gamma(a*(i+1))
    dterm = dterm_list[[i+1]]
    psum <- psum + hterm*dterm
  }
  return(psum)
}

ProbTrans = function(t, k, beta){
  return({beta*t}^k * exp(-beta*t)/factorial(k))
}

Integrand = function(h, delta_time, N, k, lambda, beta, dterm_list){
  return(ProbHeight(h, N, lambda, delta_time, dterm_list) * ProbTrans(h+delta_time, k, beta))
}

ProbKTransmissions = function(N, k, t1, t2, lambda, beta){
  upper_limit = Inf
  delta_time = abs(t1-t2)
  dsum <- 0.0
  dterm_list = list()
  for (i in seq(0, N)){
    dterm = {lambda*delta_time}^{N-i} * exp(-lambda*delta_time)/factorial(N-i)
    dsum = dsum + dterm
    dterm_list[[i+1]] = dterm
  }
  result = tryCatch(
    stats::integrate(Integrand, lower=0, upper=upper_limit, delta_time,N,k,lambda,beta,dterm_list),
    error=function(error_message) {
      message("Integration error in ProbKTransmissions().")
      print(paste0(t1,':',t2))
      message(error_message)
      return(NA)
    }
  )
  return(result$val/dsum)
}

