#!/bin/bash
gnuplot << EOF

plot "overlap.txt" u 1:2 w lp title "overlap"
replot "overlap.txt" u 1:2:3 w ye notitle
replot "overlap_inf.txt" u 1:2 w l title "independent init"
set grid
set xlabel "MC steps"
set ylabel "1-Hamming "
set term jpeg
set logscale x
set out "overlap.jpeg"
rep

EOF
