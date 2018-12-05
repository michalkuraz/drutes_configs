my_optim <- function(ini_vals, alpha, beta, gamma, maxit, abstol){
  source("drut_eval.R")
lower = c(2e-3 , 1.1, 0.3, 0.8)
upper = c(3e-3, 2.3, 0.6, 1.3)
ndeps = c(1e-4, 1e-2, 0.05, 1e-1)
  
  out <- optim( 
    par = ini_vals,
    fn = eval_fun_mead,
    lower = lower,
    upper = upper,
    control = list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      maxit = maxit,
      abstol = abstol,
      ndeps = ndeps
    )
  )
  print(out)
  if(out$convergence != 0){
    print("something is fucked up!!!!!!")
  }
  write(c(out$par, out$value, out$counts[1], out$convergence), "results/results.out", )
}
