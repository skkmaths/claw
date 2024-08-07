// Define the characteristic length of the mesh
lc = 0.01;

// Define the points of the geometry
Point(1) = {-1, -1, 0, lc};
Point(2) = {1, -1, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {-1, 1, 0, lc};

// Define the lines of the geometry
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the curve loop and plane surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Generate the triangular mesh using the Delaunay algorithm
Mesh.Algorithm = 5; // 5 corresponds to the Delaunay algorithm

// Recombine the triangular elements into quadrilateral elements
Recombine Surface {1};

// Generate the mesh
Mesh 2;
