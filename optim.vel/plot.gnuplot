set terminal png
set output "plot.png"
plot "drutes.conf/ADE/SOIL/A/control-a-obs-bottom.dat" u 1:2 , "out/obspt_ADER_in_liquid-1.out" u 1:2

