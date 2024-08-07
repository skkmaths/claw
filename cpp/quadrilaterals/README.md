## TwoD FVM code for linear advection equation on unstructured quadrilateral grids

Grid files are written in grid.h, it creates an object of Nodes, Faces and Cells. Nodes and Cells are assigned id

ranging from 0 to size of Nodes and Cells respectively. For each Face, we can access the Left and Right cells. 

If a face is a part of the boundary, then right cell is always marked as -1 as there is no right cell. 

##Things to do

1. Implementation of boundary condition: Periodic, Dirichlet, Inflow, Outflow etc.
2. Fixing CFL and use of various numfluxes.
3. Need to assign the ID s for Cells and update the vtk file.
4. Initialize the solution by finding the average through quadrature

## How to run

$ make clean

$ gmsh -2 meshfile.geo -o mesh.msh

$ make

$ ./run

Solutions are saved in ./sol directory. You can visualize using ViSiT or Paraview