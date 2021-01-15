# library('GA')
# library('DEoptim')
# library('foreach')
# library('doMC')

ReferenceInfl <- function(t) {

  A = 0.113173248618806
  k = 0.00000515892375196475
  s = 0.000855461248783132

  -(s/(2*sqrt(t*3600))*exp(-A*sqrt(t*3600))+k)*360000 
  
}

ReferenceCumInfl <- function(t) {

  A = 0.113173248618806
  k = 0.00000515892375196475
  s = 0.000855461248783132
  
  (s/A*(1-exp(-A*sqrt(t*3600)))+k*t*3600)*100 
  
}


DRU = function(){
  
  data = read.table('out/obspt_RE_matrix-1.out')
  errorCumi = mean(abs(data[,5] - ReferenceCumInfl(data[,1])))
  
  dt=1e-2
  
  dimen = dim(data)[1]
  
  t2 = data[dimen,1]
  
  t1 = data[dimen-1,1]
  
  drefdt = (ReferenceCumInfl(t2) - ReferenceCumInfl(t1))/(t2-t1)
  
  ddatadt = (data[dimen,5] - data[dimen-1,5])/(t2-t1)
  
  error_der = abs(ddatadt-drefdt)
  
  for (i in 1:(dimen-1)) {
    if (data[i,1] < 0.075 & data[i+1,1] >= 0.075) {
      error3 = abs(data[i+1,5] - ReferenceCumInfl(data[i+1,1]))
    }
  }
      
  toterr = errorCumi + error_der + error3
  write(toterr, file="objfnc.val")
 # write(error_der, file="objfnc.val", append=TRUE)
 # write(error3, file="objfnc.val", append=TRUE)


}


DRU()




