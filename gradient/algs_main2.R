# Optimization set-up to be used in combination with DRUtES

source('sourcealgs.R') # creates connection to algorithms
### algorithm options
## single objective
# 1 - pso - linear time adaptive particle swarm optimiyation with some strategies to leave local optima
# 2 - pso_bn - same as pso with a bad neighbourhood approach, where bad neighbourhoods are assigned and population drawn towards best
# 3 - pso_sce - same as pso with shuffling complexes mechanism
# 4 - pso_sce_bn - pso with bn and sce
# 5 - tlbo - learning experience teaching learning based algorithm 
# 6 - tlbo_bn - tlbo with bn
# 7 - tlbo_sce - tlbo with sce
# 8 - tlbo_sce_bn - tlbo with sce and bn
# 9 - Nelder mead - uses optim()

## biobjective (two objectives)
# 21 - mo_PSO - particle swarm optimization with dynamic neighborhood approach. Only two objecttives can be solved at this time.

## choose algorithm
alg <- 9

### set up

## Nelder-MEad simplex
# Only set-up for Nelder-Mead, optim() function. Initial values, for nelder mead, an initial values are needed, not used for other algorithms
ini_vals <- c(a = 3.802e-3 , n = 1.279, ths = 0.594, Ks = 1.165)
# alpha is the reflection coefficient
alpha <- 1
# beta is the contraction
beta <- 0.5
# gamma is the expansion coefficient
gamma <- 2
# maximum iteration
maxit <- 1000
#
abstol <- 1e-5

# other
# ranges, need to be in the same order as in drut_opti.sh script, equal min and max values will be ignored
mins <- c(amin = 0.01 , nmin = 2, thsmin = 0.4, Ksmin = 1)#
maxs <- c(amax = 0.01 , nmax = 2, thsmax = 0.4, Ksmax = 1)# 
logscale <- T # sample uniform between min and max, but use exponent of sampled value as input



# optimization
pop <- 4 # population 
complexes <- 2 # only important for algorithms with sce, will be ignored otherwise

# reinitilization
reini_prop <- 0.95 # reinitilization probability for population to not be reinitialized. Set to value higher than 1 if population should not be reinitialized. 
red_fac <- 0.99 # reduction fac of reiniprop. After each algorithm the reinitilization probability will be reduced to make reinitilization more liklely. Set to 1 or higher if not wanted.

# termination criteria
gen <- 1 # maximum number of generations, maximum number of function calls = pop*gen
optimum <- 0 # ideal value of function value of the optimum, e.g. using RMSE is 0, does not need to be realistic 
conv <- 1e-6 # convergence criteria, if the global best has not changed more than this.
conv_gen <- 100 # if global best has not changed more than conv for conv_gen consecutive generations, the optimization terminates

# parallel executions
para <- 1

## Minimization or Maximization problem?
# all problems will be evaluated with minimzation, but for maximization will be converted to minimzation problem
minimize <- TRUE # TRUE= Minimization, FALSE=Maximization

## output options
# "all" - entire population is written into seperate files every generation, after optimization finishes, all global best will be printed into one file
# "gbest" - only global best of each generation will be printed into file after optimization finishes, for mo_PSO only final non-dominated solutions will be printed
output <- "all" 


## restart optimization from old file
restart <- F
filename <- "solution_files/result_out.csv"

##
# function call in sourcealgs.R
callopti(  alg = alg, pop = pop, complexes = complexes,
           mins = mins, maxs = maxs,
           gen = gen,
           optimum = optimum,
           conv = conv, conv_gen = conv_gen,
           para = para,
           minimize = minimize,
           reini_prop = reini_prop,
           red_fac = red_fac,
           output = output,
           restart = restart, 
           filename = filename,
           logscale = logscale,
           ini_vals = ini_vals,
           alpha = alpha,
           beta = beta,
           gamma = gamma,
           maxit = maxit,
           abstol = abstol)
