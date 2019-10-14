
 
 set xlabel "time [s]"
 
 set ylabel "infiltration rate [m/s]" 
 
 set grid

 set terminal postscript colour "Helvetica" 20 lw 2

 set key bottom right
 
 set output "plot.eps"
 
 plot "obspt_runoff-1.out" u 1:3 w l lc rgb "blue" lw 2 title "model data" , "drutes.conf/kinwave/inputs.dat" u 1:2 title "experimental data"

set terminal png

set output "plot.png"

replot

