\// This script generates an unstructured quadrilateral mesh for a simple square domain

// Specify the file format version
Mesh.MshFileVersion = 4.1;

// Define the points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

// Define the lines (edges) connecting the points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the line loop and the surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Mesh the surface using quadrangles
Mesh.Algorithm = 8; // Use the "Frontal-Delaunay" for quads
Mesh.RecombineAll = 1; // Recombine triangles into quads
Mesh.SubdivisionAlgorithm = 1; // Apply quad-dominant algorithm

// Set the element size
Characteristic Length {1, 2, 3, 4} = 0.1;

// Generate the mesh
Mesh 2;

// Save the mesh to a .gmh file
Save StrFormat('unstructured_quad_mesh.gmh');
