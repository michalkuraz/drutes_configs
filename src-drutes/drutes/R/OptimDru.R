library('GA')
library('DEoptim')
library('foreach')
library('doMC')

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


DRU = function(pars){
  
  samp = c(0:9, letters, LETTERS)
  dir = paste(sample(samp, 16), collapse="")
  parVec = paste(as.vector(pars), collapse = " ")
  system(paste('./bash.sh', dir, parVec))
  Sys.sleep(0.01)data
  data = read.table(paste0(dir, '/out/obspt_RE_matrix-1.out'))
  system(paste('rm -rf', dir))  
#   errorInfl = mean(abs(data[,4] - ReferenceInfl(data[,1])))
  errorCumi = mean(abs(data[,5] - ReferenceCumInfl(data[,1])))

#   MAE = errorInfl + 200 * errorCumi
  MAE = errorCumi

  MAE

}


## NASTAVENI UMISTENI
setwd('/home/vojta/OptimalizaceV1/opt1')

nPar = 5
ranges = matrix(c(0.005,  0.05,
                  1.050,  2.10,
                  0.500,  0.90,
                  0.050, 36.00,
                  0.000,  0.10), 
                  nPar, 2, byrow = T)
# rownames(ranges) = c('a', 'n', 'Ths', 'Ks', 'Ss')

## Paralelizace

nProc = 12

registerDoMC(nProc)

ptm = proc.time()

DE = DEoptim(DRU, ranges[ ,1], ranges[ ,2], DEoptim.control(itermax = 150, parallelType = 2))

message(paste('elapsed time:',(proc.time()-ptm)['elapsed']/3600), 'h')

saveRDS(DE, 'DEoptimalizace.rds')
# DE = readRDS('DEoptimalizace.rds')

txtout = list(criterion = DE[[1]]$bestval, parameters = DE[[1]]$bestmem)
unlisted = unlist(txtout)
write('Criterion:', "DEtxtOut", ncolumns=1000)
write(unlisted[1], "DEtxtOut", append=TRUE, ncolumns=1000)
write('Parameters:', "DEtxtOut", append=TRUE, ncolumns=1000)
write(unlisted[2:length(unlisted)], "DEtxtOut", append=TRUE, ncolumns=1000)


bestPar = DE[[1]][['bestmem']]
parVec = paste(as.vector(bestPar), collapse = " ")
if(file.exists('Recomputation')) system('rm -rf Recomputation') 
system(paste('./bash.sh', 'Recomputation', parVec))

data = read.table('Recomputation/out/obspt_RE_matrix-1.out')

pdf('./Pictures/OptV1Infl.pdf')
plot(data[,1], ReferenceInfl(data[,1]), type = 'p', main = 'V1 opt. Infl')
lines(data[,1], data[,4], col = 'red')
legend('bottomright', c('measured', 'computed'), pch = c(1,-1), lty = c(-1,1),col = c('black','red'))
dev.off()

pdf('./Pictures/OptV1CumInfl.pdf')
plot(data[,1], ReferenceCumInfl(data[,1]), type = 'p', main = 'V1 opt. CumInfl')
lines(data[,1], data[,5], col = 'red')
legend('bottomright', c('measured', 'computed'), pch = c(1,-1), lty = c(-1,1),col = c('black','red'))
dev.off()





