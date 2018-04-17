 A = 0.113173248618806
 k = 0.00000515892375196475
 s = 0.000855461248783132

 f(x) =  (s/A*(1-exp(-A*sqrt(x*3600)))+k*x*3600)*100
 
 set xlabel "time [hrs]"
 
 set ylabel "cumulative infiltration [cm]" 
 
 set xrange [0:0.8]
 
 set grid

 set terminal postscript colour "Helvetica" 20 lw 2

 set key bottom right
 
 set output "plot.eps"
 
 plot f(x) lt 1 title "experimental data"  , "out/obspt_RE_matrix-1.out" u 1:5 w l lt 3 title "RE model data" 



