a = 0.113173248618806
k = 0.00000515892375196475
s = 0.000855461248783132

cumi <- function(t) {
  return(s/a*(1-exp(-a*sqrt(t)))+k*t)
  }

cumii <- function(t) {
  return( ((-2*s*(-a^(-2) - sqrt(t)/a))/exp(1)^(a*sqrt(t)) + s*t + (a*k*t^2)/2)/a )
  }

data = read.table("out/obspt_RE_matrix-1.out")

val1 = data[1,4]
val2 = data[1,5]

dim = length(data[,1])

for (i in 2:dim) {
  val1 = val1 + data[i,4]*(data[i,1]-data[i-1,1])
  val2 = val2 + data[i,5]*(data[i,1]-data[i-1,1])
}

val1 = -val1

error1 = abs(val1 - cumi(data[dim,1]))

error2 = abs(val2 - cumii(data[dim,1]))