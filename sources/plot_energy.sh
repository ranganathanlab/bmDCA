#!/bin/bash
gnuplot << EOF

plot "energy.dat" u 1:2 w lp title "Average energy"
replot "energy.dat" u 1:2:3 w ye notitle
set xlabel "MC steps"
set grid
set ylabel "Energy"
set term jpeg
set out "energy.jpeg"
rep

EOF
