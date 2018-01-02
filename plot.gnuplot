set terminal png
set output "1.png"
plot "out/obspt_ADER_in_liquid-1.out" u 1:2 w l , "drutes.conf/inverse_modeling/biochar.dat" u 1:2
set output "2.png" 
plot "out/obspt_ADER_in_liquid-2.out" u 1:2 w l , "drutes.conf/inverse_modeling/biochar.dat" u 1:3
set output "3.png" 
plot "out/obspt_ADER_in_liquid-3.out" u 1:2 w l , "drutes.conf/inverse_modeling/biochar.dat" u 1:4 
