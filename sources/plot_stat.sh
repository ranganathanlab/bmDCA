#!/bin/bash
gnuplot << EOF

plot "errorlist.txt" u 3:4 w p pt 0 title "1 point marginals"
replot "my_corr.dat" u 1:2 w p pt 0 title "correlations"
set grid
set xlabel "MSA"
set ylabel "MC "
set term jpeg
set out "cfr_stat.jpeg"
rep

EOF
