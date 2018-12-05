
eval_fun <- function(ln_id,obj=1){
  system2("./drut_opti.sh",wait=T)
  result=matrix(ncol=obj,nrow=ln_id)
  for(i in 1:ln_id){ 
    test=read.table(paste(i,"/out/objfnc.val",sep=''), quote="\"",comment.char="#", sep="",skip=1)
    result[i,1]=test$V1[1]
    if(obj>1){
      #result[i,2]=test$V1[2]
    }
  }
  return(result)
}

eval_fun_mead <- function(pars_in){
  write(c("p",pars_in),'pars.in', ncol = length(pars_in)+1)
  system2("./drut_opti.sh",wait=T)
  obj.val <- read.table(paste("1/out/objfnc.val",sep=''), quote="\"",comment.char="#", sep="",skip=1)
  result <- obj.val$V1[1]
  return(result)
}
