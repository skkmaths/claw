# This code is written for u_t+ a(x,y)u_x + b(x,y)u_y = 0, and u_t + (u^2/2)_x + (u^2/2)_y = 0

# 2D Linear advection equation and Burgers equation

Following schemes are implemented

1. Lax Wendroff scheme for Linear advection 
2. Second order MUSCL scheme with, ssprk22 time discretization for Linear advection and Burgers equation

Convergence test for order of accuracy works for 

1. LW scheme ---> Linear advection ( including variable advection)

2. MUSCL scheme---> Linear advection ( including variable advection) and Burgers equation

Default code runs the first order FV scheme

```
python clw2d.py -Tf 1.0 -cfl 0.9 -save_freq 1 -plot_freq 0
```
For others

```
python clw2d.py -Tf 1.0  -pde linear -scheme rk2 -limit mmod
```

```
python clw2d.py -scheme lw -Tf 1.0 -pde linear -ic sin2pi 
```

## To test the order of accuracy, run the script file
```
sh runconvergence "20 40 80 160 320 640"
python plotrate.py
```
