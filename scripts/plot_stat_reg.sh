#! /bin/bash
set -eu

gnuplot << EOF

plot "errorlist.txt" u 3:8 w p pt 0 title "1 point marginals"
replot "my_corr.dat" u 1:3 w p pt 0 title "correlations"
set grid
set xlabel "MSA"
set ylabel "MCMC"
set term jpeg
set out "cfr_stat_reg.jpeg"
rep

EOF
