
 
 set xlabel "time [s]"
 
 set ylabel "conc. flux [kg.m^{3}.s^{-1}]" 
 
 set grid

 set terminal postscript enhanced colour "Helvetica" 20 lw 2

 set key bottom right
 
 set output "plot.eps"
 
 plot "out/obspt_debris_flow-1.out" u 1:4 w l lc rgb "blue" lw 2 title "model data" ,"drutes.conf/inverse_modeling/data2model" u 1:2 title "experimental data"

set terminal png

set output "plot.png"

replot

