
mo_PSO=function(pop,complexes,dim,xmin,xmax,gen,printall=T,maxeval,start_shuffle_prob=0.995,red_fac=0.99,minimize,restart=F,filename=NULL,logscale=F){
  source("drut_eval.R")
  if(minimize){
    facmin=1
  }else{
    facmin=-1
  }
  
  fit_func<-function(ln_id){
    result=facmin*eval_fun(ln_id,obj=2)
    return(result) # this should be maxeval rows if combined opti is used
  }
  
  reshuffle_prob=start_shuffle_prob
  if(length(xmax)!=dim){
    stop('xmax: wrong number of boundaries defined. Dimension and length(xmax) differ')
  }
  if(length(xmin)!=dim){
    stop('xmin: wrong number of boundaries defined. Dimension and length(xmin) differ')
  }
  ###################
  #### Initialisation
  ###################
  
  vmax=0.3*(xmax-xmin)
  swarm=matrix(t(runif(pop*dim,xmin,xmax)),nrow=pop,byrow=T) 
  vel=matrix(rep(0,dim*pop),nrow=pop,byrow=T) 
  xmax_mat=matrix(rep(xmax,pop),nrow=pop,byrow=T)
  xmin_mat=matrix(rep(xmin,pop),nrow=pop,byrow=T)
  
  if(restart){
    pbest_loc=read.csv(filename,header=T)
    pbest_loc=as.matrix(pbest_loc[,2:(dim+1)])
    swarm=pbest_loc
  }
  
  # parallel evaluation
  index=sort(rep(1:ceiling((pop/maxeval)),maxeval))
  index=index[1:pop]
  for(i in 1:ceiling(pop/maxeval)){
    ln_id=(length(index[index==i]))
    pars_in=cbind(rep('p',ln_id),matrix(swarm[index==i,],nrow=ln_id))
    if(logscale){
            pars_in[,2:(dim+1)]=exp(as.numeric(pars_in[,2:(dim+1)]))
    }
    write(t(pars_in),'pars.in',append = F,ncol=dim+1)
    if(i>1){
      result=rbind(result,fit_func(ln_id))
    }else{
      result=fit_func(ln_id)
    }
  }
  pbest=result
  pbest_loc=swarm
  ranks_ob1=rank(result[,1],ties.method = "random")
  ranks_ob2=rank(result[,2],ties.method = "random")
  write.csv(cbind(pbest_loc,pbest),paste('results/population_gen_init.csv',sep=""))
  
  ## Neighbourhoods
  # 1. each particle only has two neighbours to exchange information, which are neighhbours based on objective function 1
  # 2. The particle is going to move towards the best particle in that small neighbourhood
  # based on objective function 2
  # 3. later pbest is only updated when non-dominated solution presetn, e.g. solution has improved for both objective functions
  # 4. can I do this vectorized as well?
  nhood=function(i){
    ind=ranks[i]
    if(ind==1){
      n1=which(ranks==(ind+1))
      n2=which(ranks==(ind+1))
    }else{
      if(ind==pop){
        n1=which(ranks==(ind-1))
        n2=which(ranks==(ind-1))
      }else{
        n1=which(ranks==(ind-1))
        n2=which(ranks==(ind+1))
      }
    }
    nhood=c(i,n1,n2)
    nmin=which.min(ranks2[nhood])
    nbests=nhood[nmin]
    return(nbests)
  }
  
  mode1=T
  switch=0
  if(mode1){
    ranks=ranks_ob1
    ranks2=ranks_ob2
  }else{
    ranks=ranks_ob2
    ranks2=ranks_ob1
  }
  pops=matrix(1:pop,ncol=1)
  nbests=apply(pops,1,nhood)
  nbests_loc=swarm[nbests,]

  ##############
  # Update
  ##############
  k=1 # generation index
   # to be returned in the end
  while(k<=gen){
    print(paste("generation",k))
    mode1=!mode1
    wmax=(0.9-0.2)*(gen-k)/gen+0.2 # based on Suganthan, Roshida and yoshida et al. #0.9
    #wmin=0.2
    vmax=(0.3*(xmax-xmin)-(xmax-xmin)/20)*(gen-k)/gen+(xmax-xmin)/20#changed from 0.5 to 0.3 and from min 10 to 20 (run 18). changed 0.1 to 0.5  
    #particle specific vmax, low when gbest (improve exploitation) and high when ranked badly (go on exploration)
    vmax_part=matrix(rep(vmax,pop),nrow=pop,byrow=T)/(pop-ranks+1)
    w=wmax
    c2=0.5+(2.5-0.5)*k/gen ## increasing attraction to global best
    c1=0.1+(1.5-0.1)*(gen-k)/gen # decreasing attraction to personal best
    vel=w*vel+c1*runif(pop)*(pbest_loc-swarm)+c2*runif(pop)*(nbests_loc-swarm)
    
    too_fast=which(vel>vmax_part)
    vel[too_fast]=vmax_part[too_fast]
    too_fast=which(vel<(-vmax_part))
    vel[too_fast]=-vmax_part[too_fast]
    
    swarm=swarm+vel
    # reflection back into the space when 'hitting' the boundary
    bound_max=which(swarm > xmax_mat)
    bound_min=which(swarm < xmin_mat)
    lnmin=length(bound_min[bound_min])
    lnmax=length(bound_min[bound_min])
    while(lnmin>0 | lnmax>0){
      swarm[bound_max] = xmax_mat[bound_max]-(swarm[bound_max]-xmax_mat[bound_max])
      swarm[bound_min] = xmin_mat[bound_min]-(swarm[bound_min]-xmin_mat[bound_min])
      bound_max=which(swarm > xmax_mat)
      bound_min=which(swarm < xmin_mat)
      lnmin=length(bound_min[bound_min])
      lnmax=length(bound_min[bound_min])
    }
    
    
    # parallel evaluation
    for(i in 1:ceiling(pop/maxeval)){
      ln_id=(length(index[index==i]))
      pars_in=cbind(rep('p',ln_id),matrix(swarm[index==i,],nrow=ln_id))
      if(logscale){
              pars_in[,2:(dim+1)]=exp(as.numeric(pars_in[,2:(dim+1)]))
      }
      write(t(pars_in),'pars.in',append = F,ncol=dim+1)
      if(i>1){
        result=rbind(result,fit_func(ln_id))
      }else{
        result=fit_func(ln_id)
      }
    }
    ranks_ob1=rank(result[,1],ties.method = "random")
    ranks_ob2=rank(result[,2],ties.method = "random")
    
    if(mode1){
      ranks=ranks_ob1
      ranks2=ranks_ob2
    }else{
      ranks=ranks_ob2
      ranks2=ranks_ob1
    }
    
    reshuffle_prob=reshuffle_prob*start_shuffle_prob
    a=runif(1)
    # reinitilize worst half if prob reached
    if(a>reshuffle_prob){
      worst=length(ranks[which(ranks>floor(pop/2))])
      swarm[which(ranks>floor(pop/2)),]=matrix(t(runif(worst*dim,xmin,xmax)),nrow=worst,byrow=T) 
      reshuffle_prob=start_shuffle_prob
      for(i in 1:ceiling(pop/maxeval)){
        ln_id=(length(index[index==i]))
        pars_in=cbind(rep('p',ln_id),matrix(swarm[index==i,],nrow=ln_id))
        if(logscale){
                pars_in[,2:(dim+1)]=exp(as.numeric(pars_in[,2:(dim+1)]))
        }
        write(t(pars_in),'pars.in',append = F,ncol=dim+1)
        if(i>1){
          result=rbind(result,fit_func(ln_id))
        }else{
          result=fit_func(ln_id)
        }
      }
      ranks_ob1=rank(result[,1],ties.method = "random")
      ranks_ob2=rank(result[,2],ties.method = "random")
      
      if(mode1){
        ranks=ranks_ob1
        ranks2=ranks_ob2
      }else{
        ranks=ranks_ob2
        ranks2=ranks_ob1
      }
    }
    #new pbest, only if both are better
    newpbest=apply(result<=pbest,1,all)
    pbest[newpbest,]=result[newpbest,]
    pbest_loc[newpbest,]=swarm[newpbest,]
    nbests=apply(pops,1,nhood)
    nbests_loc=swarm[nbests,]
    if(printall){
      write.csv(cbind(pbest_loc,pbest),paste('results/population_gen',k,'.csv',sep=""))
    }

    k=k+1  
  }
  # pareto optimal front
  testdom=function(testobj){
    test=matrix(rep(testobj,pop),ncol=2,byrow=T)
    bool=!any(apply(test>pbest,1,all))
  }
  nond=apply(pbest,1,testdom)
  pareto=cbind(pbest_loc[nond,],pbest[nond,])
  results=pareto
  write.csv(pareto,'results/non_dom_solutions.csv')
  return(results)
}
