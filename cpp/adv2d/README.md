## TwoDProblem Solver
This project implements a 2 dimensional finite volume scheme for linear advection equation, with constant coefficients. This also uses OpenMP for parallelization. But, seems to be not efficient. The program takes several command-line arguments to configure the problem's parameters and performs computations accordingly.
## Features
- Solves two-dimensional problems with configurable grid size and time step.
- Uses OpenMP for parallel computation to enhance performance.
## Requirements
- C++ compiler (e.g., g++)
- OpenMP support
- This code is for the discretization of u_t+u_x+u_y = 0
## Compile 
- $make clean
- $make 
- $./run
- By default it uses number of cells, nx, ny, and final time Tf and cfl.
- You can also run as, eg: $./run -nx 50 -ny 50 -Tf 1.0 -cfl 0.4 -save_freq 10 -scheme so

## E.O.C. computation 
- E.O.C. works only for smooth periodic test case
- $sh runcongs.sh "20 40 80 160"
- While computing E.O.C. set -save_freq 0 to save computational time
- This will produce and error.txt file and can print the order of convergence as follows
- $python printrate -f erro.txt
- Both first and second -order schemes gives the expected E.O.C
## Visualization
-You can use visit for visualization
## To do
- Avoid print data in each iteration while computing E.O.C

