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

Solutions are saved to .plt files: use Visit to plot

Instruction to install visit in centos-lynux machine

1. Download the tar.gz file from https://wci.llnl.gov/simulation/computer-codes/visit/executables

2. Downlad install script https://github.com/visit-dav/visit/releases/download/v3.1.4/visit-install3_1_4

3. Rename the file to visit-install

4. Chmod u+x visit-install3_1_4

5. ./visit-install3_1_4 3.1.4 linux-x86_64-centos8 /usr/local/

6  You can install to /usr/local/visit as well

7. The directory should contain the .tar.gz file 
