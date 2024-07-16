# TwoDProblem Solver

This project implements a two-dimensional problem solver using OpenMP for parallelization. The program takes several command-line arguments to configure the problem's parameters and performs computations accordingly.

## Features

- Solves two-dimensional problems with configurable grid size and time step.
- Uses OpenMP for parallel computation to enhance performance.

## Requirements

- C++ compiler (e.g., g++)
- OpenMP support

# Run instructions

- This code is for the discretization of u_t+u_x+u_y = 0
# Compile 

- $make clean
- $make 
- $./twodproblem
- By defaulty it uses number of cells, nx, ny, and final time Tf and cfl.
- For your choice you can use $./twodproblem -nx 50 -ny 50 -Tf 1.0 -cfl 0.4

# E.O.C computation 
-$sh runcongs.sh "20 40 80 160"
-This will produce and error.txt file and can print the order of convergence as follows
-$python printrate -f erro.txt
