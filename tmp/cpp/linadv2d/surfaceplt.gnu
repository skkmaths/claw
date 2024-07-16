#surface plot

set term postscript eps enhanced color
set output '3dplot.eps'



set surface 
set hidden3d
set xlabel 'x-values'
set ylabel 'y-values'
set zlabel 'f(x,y)'



splot "solution-100.dat" u 1:2:3 t'f-values' w l


