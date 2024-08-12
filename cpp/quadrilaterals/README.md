## TwoD FVM code for linear advection equation on unstructured quadrilateral grids

Grid files are written in grid.h, it creates an object of Nodes, Faces and Cells. Nodes and Cells are assigned id

ranging from 0 to size of Nodes/Cells.  For each Face, we can access the Left and Right cells. 

If a face is a part of the boundary, then right cell is always marked as -1, as there is no right cell. 

## Things to do

1. Implementation of boundary condition: Periodic, Dirichlet, Inflow, Outflow etc.
2. Need to assign the ID s for Cells and update the vtk file.
3. Initialize the solution by finding the average through quadrature

## How to run
Use one of this in the make file: 

INCLUDES = -I/usr/local/include 

LDFLAGS = -L/usr/local/lib  -lgmsh

or 

INCLUDES = -I/opt/homebrew/include 

LDFLAGS = -L/opt/homebrew/lib  -lgmsh

Make sure that boundaries defined in the .geo file and grid.h files coincide

$ make clean

$ gmsh -2 meshfile.geo -o mesh.msh 

$ make

$ ./run

You can also use the following to specify the mesh size

$gmsh -2 strumesh.geo -setnumber lc 0.01 -o mesh.msh

To verify the EOC do the following

$ sh runcongs.sh "0.1 0.05 0.025 0.0125 0.00625"

$python errorplot.py

While doing the EOC, make sure that, ic = "expo", domain [-1,1]X[-1,1], set the boundary flux computation to zero, set time = 2*$\pi$

Solutions are saved in ./sol directory. You can visualize using ViSiT or Paraview