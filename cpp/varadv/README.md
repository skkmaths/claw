## TwoD Problem Solver
This project implements two dimensional finite volume scheme for linear advection equation, with varying coefficients. The program takes several command-line arguments to configure the problem's parameters and performs computations accordingly.
## Reference
- The example considered here is taken from page 460, solid body rotation,cha: Mulitidimensional scalar equations, of the book;  Randall J. Leveque, Finite volume methods for Hyperobolic problems, Cambridge text in applied Mathematics.
## Requirements
- C++ compiler (e.g., g++)
- This code is for the discretization of $u_t + div(a(x,y)u) = 0, a = (a_1,a_2)$
## Compile 
- $make clean
- $make 
- $./run
- By default it uses number of cells, nx, ny, and final time Tf and cfl.
- You can also run as, eg: $./run -nx 50 -ny 50 -Tf 1.0 -cfl 0.4 -save_freq 10 -scheme so
- Set the  domain by assigning boundary nodes in .run function
- Set the type of boundary conditions  .run
- Exact solution is written for a=(-y,x)
- You can adjust the advection_speed function to get linear advection equations
## E.O.C. computation 
- For E.O.C. computation set bc as periodic test case
- Final time = 2*pi for ic = "expo"
- Final time = 3  for ic = "sin"
- Check if we are using the full domain [-1,1] X [-1,1]
- $sh runcongs.sh "20 40 80 160"
- While computing E.O.C. set -save_freq 0 to save computational time
- This will produce and error.txt file and can print the order of convergence as follows
- $python printrate -f erro.txt
- Both first and second -order schemes gives the expected E.O.C, with the "expo" test case.
## Visualization
-You can use visit for visualization
## To do
- Avoid print data in each iteration while computing E.O.C

