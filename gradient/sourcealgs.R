source('mo_PSO.R')
source('pso.R')
source('tlbo.R')
source('my_optim.R')

callopti=function(alg, pop, complexes,
                  mins, maxs,
                  gen, optimum=0, conv , conv_gen,
                  para, minimize=T,
                  reini_prop, red_fac,
                  output,
                  restart=F, filename, logscale,
                  ini_vals, alpha, beta, gamma,
                  maxit, abstol){
  
  # check if alg option exist
  alg_opts=c(1:9,21)
  if(!any(alg==alg_opts)){
    stop('This algorithm option is invalid. Choose 1-8 for single objective optimization and 21 for biobjective optimization')
  }
  
  # set-up parallel computation
  system(paste('./opti_setup.sh', para))
  
  dim=length(mins)
  printall <- switch(output,
         "all"=TRUE,
         "gbest"=FALSE
        )
  if(!dir.exists('results')){
    system('mkdir results')
  }
 

  # select algorithm
  switch(as.character(alg),
         "1"=return(PSO_all(pop,complexes = 1,dim,
                            mins,maxs,gen,printall=printall,
                            maxeval=para,start_shuffle_prob=reini_prop,
                            red_fac=red_fac,optimum,conv,conv_gen,
                            minimize,bn=F,sce=F,restart=restart,
                            filename=filename,logscale=logscale)),
         "2"=return(PSO_all(pop,complexes = 1,dim,
                            mins,maxs,gen,printall=printall,
                            maxeval=para,start_shuffle_prob=reini_prop,
                            red_fac=red_fac,optimum,conv,conv_gen,
                            minimize,bn=T,sce=F,restart=restart,
                            filename=filename,logscale=logscale)),
         "3"=return(PSO_all(pop,complexes,dim,
                            mins,maxs,gen,printall=printall,
                            maxeval=para,start_shuffle_prob=reini_prop,
                            red_fac=red_fac,optimum,conv,conv_gen,
                            minimize,bn=F,sce=T,restart=restart,
                            filename=filename,logscale=logscale)),
         "4"=return(PSO_all(pop,complexes,dim,
                            mins,maxs,gen,printall=printall,
                            maxeval=para,start_shuffle_prob=reini_prop,
                            red_fac=red_fac,optimum,conv,conv_gen,
                            minimize,bn=T,sce=T,restart=restart,
                            filename=filename,logscale=logscale)),
         "5"=return(TLBO_all(pop,1,dim,
                             mins,maxs,gen,printall=printall,
                             maxeval=para,start_shuffle_prob=reini_prop,
                             red_fac=red_fac,optimum,conv,conv_gen,
                             minimize,bn=F,sce=F,restart=restart,
                             filename=filename,logscale=logscale)),
         "6"=return(TLBO_all(pop,1,dim,
                             mins,maxs,gen,printall=printall,
                             maxeval=para,start_shuffle_prob=reini_prop,
                             red_fac=red_fac,optimum,conv,conv_gen,
                             minimize,bn=T,sce=F,restart=restart,
                             filename=filename,logscale=logscale)),
         "7"=return(TLBO_all(pop,complexes,dim,
                             mins,maxs,gen,printall=printall,
                             maxeval=para,start_shuffle_prob=reini_prop,
                             red_fac=red_fac,optimum,conv,conv_gen,
                             minimize,bn=F,sce=T,restart=restart,
                             filename=filename,logscale=logscale)),
         "8"=return(TLBO_all(pop,complexes,dim,
                             mins,maxs,gen,printall=printall,
                             maxeval=para,start_shuffle_prob=reini_prop,
                             red_fac=red_fac,optimum,conv,conv_gen,
                             minimize,bn=T,sce=T,restart=restart,
                             filename=filename,logscale=logscale)),
         "9"=return(my_optim(ini_vals = ini_vals, alpha = alpha, beta = beta, gamma = gamma,
                             maxit = maxit, abstol = abstol)),
         "21"=return(mo_PSO(pop,complexes=1,dim,
                            mins,maxs,gen,printall=printall,
                            maxeval=para,start_shuffle_prob=reini_prop,
                            red_fac=red_fac,minimize,restart=restart,
                            filename=filename,logscale=logscale))
      )
}
