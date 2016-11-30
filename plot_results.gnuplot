#! /usr/bin/gnuplot -p

set style data lines
set terminal png size 800,600

set log y
set output 'jfnk_residual.png'
plot 'jfnk.txt' using 1:2 title 'residual'
unset log y

set output 'jfnk_alpha.png'
plot 'jfnk.txt' using 1:3 title 'alpha'

set log y
set output 'gmres_convergence.png'
plot 'gmres.txt' using 2:3:1 palette
unset log y
