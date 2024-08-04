// Define the characteristic length of the mesh
lc = 0.01;

// Define the points of the geometry
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Define the lines of the geometry
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the curve loop and plane surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define the transfinite lines and surface for structured quadrilateral meshing
Transfinite Line {1, 2, 3, 4} = 1 / lc + 1;
Transfinite Surface {1};

// Optionally, you can specify the recombination algorithm to generate quadrilateral elements
Recombine Surface {1};

// Optionally, you can specify the mesh algorithm and other parameters
// Setting the mesh algorithm to 8 (Frontal-Delaunay for quadrangles)
Mesh.Algorithm = 8; // 8 corresponds to Frontal-Delaunay for quadrangles

// Generate the mesh
Mesh 2;
