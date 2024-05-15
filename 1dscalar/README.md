## This code is written for u_t+ f(u)_x = 0, to generate First order FVM and Second order minimod+ssprk22 scheme 
To run first order scheme enter the following

```
python clw.py -pde linear -ic smooth -nc 50 -time_scheme euler -limit no
```
To run the second order scheme, MUSCL with ssprk22 scheme enter the following
```
python clw.py -pde linear -ic smooth -nc 50 -time_scheme ssprk22 -limit mmod -Tf 4 -cfl 0.5
```
## To test the order of accuracy, run the script file
```
sh runconvergence "20 40 80 160 320 640"
python plotrate.py
```
Order of accuracy test works for both advection and burger equation with periodic bc
For discontinuous flux use the -pde dflux
