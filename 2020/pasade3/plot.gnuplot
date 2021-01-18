
 
 set xlabel "time [s]"
 
 set ylabel "surface runoff [m^3.s^{-1}.m^{-1}]" 
 
 set grid

 set terminal postscript enhanced colour "Helvetica" 20 lw 2

 set key top right
 
 set output "plot.eps"
 
 plot "out/obspt_runoff-1.out" u 1:3 w l lc rgb "blue" lw 2 title "model data" , "drutes.conf/inverse_modeling/data-water3" u 1:2 title "experimental data"

set terminal png

set output "plot.png"

replot

