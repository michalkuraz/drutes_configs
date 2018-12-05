my_optim <- function(ini_vals, alpha, beta, gamma, maxit, abstol){
  source("drut_eval.R")
  out <- optim( 
    par = ini_vals,
    fn = eval_fun_mead,
    control = list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      maxit = maxit,
      abstol = abstol
    )
  )
  print(out)
  if(out$convergence != 0){
    print("something is fucked up!!!!!!")
  }
  write(c(out$par, out$value, out$counts[1], out$convergence), "results/results.out", )
}