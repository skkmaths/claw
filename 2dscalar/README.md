# 2D Linear advection equation and Burgers equation
Following schemes are implemented
1. Lax Wendroff scheme for Linear advection 
2. Second order MUSCL scheme with, ssprk22 time discretization for Linear advection and Burgers equation
3. MUSCL Hancock Scheme 

Kappa Reconstruction can be implemented and one can choose between minmod limiter and mc limiter


Convergence test for order of accuracy works for 
1. LW scheme ---> Linear advection ( including variable advection)
2. MUSCL scheme---> Linear advection ( including variable advection) and Burgers equation
3. MUSCL Scheme with Kappa Reconstruction for linear advection

Default code runs the first order FV scheme

```
python clw2d.py -Tf 1.0 -cfl 0.9 -save_freq 1 -plot_freq 0
```
For others

```
python clw2d.py -Tf 1.0  -pde linear -scheme rk2 
```

```
python clw2d.py -scheme lw -Tf 1.0 -pde linear -ic sin2pi 
```

## To test the order of accuracy, run the script file
```
sh runconvergence "20 40 80 160 320 640"
python plotrate.py is included in runconvergence.sh file and will automatically execute to print the table.

For the linear advection: set cfl = 1 (for fo scheme) and cfl =0.5 (for rk2 scheme) while computing EOC and for cfl = 0.4 (for muscl Hancock)

Note: Not much improvement in rk2 muscl and muscl hancock while using minmod limiter, issue of low EOC with respect to minmod limiter and minor improvement in result with use of mc limiter. Still working on that front. 

Note: Interesting Result are seen with use of kappa reconstruction with kappa = 1/3 (default setting) The EOC is reached upto 2 and beyond. (Not reaching EOC of three is most probably due to numerical solution is strictly not being a cell averged solution)

Solutions are saved to .plt files: use Visit to plot

Instructions to install visit in centos-lynux machine

1. Download the tar.gz file from https://wci.llnl.gov/simulation/computer-codes/visit/executables

2. Download install script https://github.com/visit-dav/visit/releases/download/v3.1.4/visit-install3_1_4

3. Rename the file to visit-install

4. Chmod u+x visit-install3_1_4

5. ./visit-install3_1_4 3.1.4 linux-x86_64-centos8 /usr/local/

6.  You can install to /usr/local/visit as well

7. The directory should contain the .tar.gz file 
