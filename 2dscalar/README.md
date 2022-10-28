## This code is written for u_t+ a(x,y)u_x + b(x,y)u_y = 0, 
minmod+ssprk22 and Lx Wendroff schemes are implemented for constant and variable advection speeds

Run as follows for minmod+ssprk22 scheme
```
python clw2d.py -scheme fv -Tf 1.0 -pde linear -ic sin2pi -ncellx 100 -ncelly 100
```

Run as follows for LxW scheme
```
python clw2d.py -scheme lw -Tf 8.0 -pde varadv -ic expo -ncellx 100 -ncelly 100 
```
## To test the order of accuracy, run the script file
```
sh runconvergence "20 40 80 160 320 640"
python plotrate.py
```
Order of accuracy test works for sin2pi with linear pde