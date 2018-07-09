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

  write(errorCumi, file="objfnc.val")

}


DRU()




