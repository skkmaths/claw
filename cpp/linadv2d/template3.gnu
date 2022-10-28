# colored plot

set term postscript eps enhanced color
set output 'colorplot.eps'



set pm3d map

set xlabel 'x-values'
set ylabel 'y-values'


splot "solution.dat" u 1:2:3 


