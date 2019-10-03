
 
 set xlabel "time [min]"
 
 set ylabel "pressure head [cm]" 
 
 set xrange [0:300]
 
 set grid

 set terminal postscript colour "Helvetica" 20 lw 2

 set key bottom right
 
 set output "plot.eps"
 
 plot "out/obspt_RE_matrix-1.out" u 1:2 w l lc rgb "blue" lw 2 title "point 1 - model" , "drutes.conf/inverse_modeling/experiment.dat" u 1:2 title "point 1 - experiment",  "out/obspt_RE_matrix-2.out" u 1:2 w l lc rgb "red" lw 2 title "point 2 - model" , "drutes.conf/inverse_modeling/experiment.dat" u 1:3 title "point 2 - experiment", "out/obspt_RE_matrix-3.out" u 1:2 w l lc rgb "violet" lw 2 title "point 3 - model" , "drutes.conf/inverse_modeling/experiment.dat" u 1:4 title "point 3 - experiment"

set terminal png

set output "plot.png"

replot

