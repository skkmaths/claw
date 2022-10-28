#contour plot

set term postscript eps enhanced color
set output 'template2.eps'



set contour
unset surface
set view map
set hidden3d
set cntrparam level 40
set xlabel 'x-values'
set ylabel 'y-values'
set zlabel 'f(x,y)'
#set cntrparam levels incr -8,1,10



splot "solution.dat" u 1:2:3 t 'f-values' w l  lw 2


