library('DEoptim')

setwd("/home/miguel/drutes/R")


drutes <- function(pp){
  system(paste('./startopt.sh', pp[1], " ",  pp[2]))
  model = read.table("../out/obspt_ADER-1.out")
  experiment = read.table("data.in")
  
  val = 0.0
  pos = 1
  for (i in 1:dim(model)[1]) {
    for (j in pos:(dim(experiment)[1]-1)) {
      if ( model[i,1] >= experiment[j,1] & model[i,1] < experiment[j+1,1]) {
	pos_mem = j
	break
      }
    }
    pos = pos_mem
    slope_exp = (experiment[pos+1,2] - experiment[pos,2])/(experiment[pos+1,1] - experiment[pos,1])
    val4exp = slope_exp * (model[i,1] - experiment[pos,1]) + experiment[pos,2]
    val = val + abs(model[i,2]-val4exp)
  }
 return(val)
}



upper <- c(10.0)
lower <- c(0.0001)



outDEoptim<-DEoptim(drutes, lower, upper, DEoptim.control(NP = 10, itermax = 100,  trace =T))
summary(outDEoptim)
## plot the output

plot(outDEoptim)
