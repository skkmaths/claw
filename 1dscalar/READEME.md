To run first order scheme enter the following

python clw.py -pde linear -ic smooth -nc 50 -time_scheme euler -theta 0.0

To run the second order scheme, MUSCL with ssprk22 scheme enter the following

$python clw.py -pde linear -ic smooth -nc 50 -time_scheme ssprk22 -theta 0.5 -Tf 4 -cfl 0.5