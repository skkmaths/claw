// Define the unit square
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Mesh settings
Mesh.Algorithm = 8; // Delquad for quadrilateral meshing
Mesh.RecombineAll = 1;
Mesh.SubdivisionAlgorithm = 1;
Mesh.ElementOrder = 1;

// Define the mesh size
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.1;

// Generate the mesh
Mesh 2;
