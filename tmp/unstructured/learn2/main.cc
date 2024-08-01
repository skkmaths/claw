#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// Node structure to represent a mesh node
struct Node {
    double x, y, z;  // Coordinates of the node
};

// Triangle structure to represent a mesh element
struct Triangle {
    std::vector<std::size_t> nodeIndices; // Indices of the nodes that form the triangle
};

// Face structure to represent a mesh face
struct Face {
    std::vector<std::size_t> nodeIndices; // Indices of the nodes that form the face
    bool isBoundary;                      // Flag to indicate if it is a boundary face
    std::size_t leftTriangleIndex;        // Index of the left triangle
    std::size_t rightTriangleIndex;       // Index of the right triangle
};

// Mesh class to store nodes, triangles, and faces
class Mesh {
public:
    std::vector<Node> nodes;          // List of nodes in the mesh
    std::vector<Triangle> triangles;  // List of triangles in the mesh
    std::vector<Face> faces;          // List of faces in the mesh
    std::vector<double> solution;     // Solution vector to store initial condition values

    // Function to read mesh data from Gmsh
    void readFromGmsh() {
        std::vector<std::size_t> nodeTags;
        std::vector<double> coord, parametricCoord;
        gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

        nodes.resize(nodeTags.size());
        for (std::size_t i = 0; i < nodeTags.size(); ++i) {
            nodes[i].x = coord[3 * i];
            nodes[i].y = coord[3 * i + 1];
            nodes[i].z = coord[3 * i + 2];
        }

        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

        for (std::size_t i = 0; i < elementTypes.size(); ++i) {
            if (elementTypes[i] == 2) { // Type 2 corresponds to 2D triangular elements
                for (std::size_t j = 0; j < elementTags[i].size(); ++j) {
                    Triangle tri;
                    tri.nodeIndices = {elementNodeTags[i][3 * j] - 1, elementNodeTags[i][3 * j + 1] - 1, elementNodeTags[i][3 * j + 2] - 1};
                    triangles.push_back(tri);
                }
            } else if (elementTypes[i] == 1) { // Type 1 corresponds to 1D line elements (faces)
                for (std::size_t j = 0; j < elementTags[i].size(); ++j) {
                    Face face;
                    face.nodeIndices = {elementNodeTags[i][2 * j] - 1, elementNodeTags[i][2 * j + 1] - 1};
                    face.isBoundary = true; // Initially set as boundary
                    faces.push_back(face);
                }
            }
        }

        // Associate faces with triangles and determine if they are boundary faces
        for (auto &face : faces) {
            face.leftTriangleIndex = face.rightTriangleIndex = -1; // Initialize as invalid index
            for (std::size_t t = 0; t < triangles.size(); ++t) {
                const auto &tri = triangles[t];
                if (std::find(tri.nodeIndices.begin(), tri.nodeIndices.end(), face.nodeIndices[0]) != tri.nodeIndices.end() &&
                    std::find(tri.nodeIndices.begin(), tri.nodeIndices.end(), face.nodeIndices[1]) != tri.nodeIndices.end()) {
                    if (face.leftTriangleIndex == -1) {
                        face.leftTriangleIndex = t; // First triangle is the left triangle
                    } else {
                        face.rightTriangleIndex = t; // Second triangle is the right triangle
                        face.isBoundary = false; // If there are two triangles, it's not a boundary face
                    }
                }
            }
        }

        // Resize solution vector to match the number of triangles
        solution.resize(triangles.size());
    }

    // Function to write initial condition to a VTK file
    void writeInitialConditionToVTK(const std::string &filename) const {
        std::ofstream file(filename);
        file << "# vtk DataFile Version 3.0\n";
        file << "Advection Initial Condition\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Write points
        file << "POINTS " << nodes.size() << " float\n";
        for (const auto &node : nodes) {
            file << node.x << " " << node.y << " " << node.z << "\n";
        }

        // Write cells
        file << "CELLS " << triangles.size() << " " << 4 * triangles.size() << "\n";
        for (const auto &tri : triangles) {
            file << "3 " << tri.nodeIndices[0] << " " << tri.nodeIndices[1] << " " << tri.nodeIndices[2] << "\n";
        }

        // Write cell types
        file << "CELL_TYPES " << triangles.size() << "\n";
        for (std::size_t i = 0; i < triangles.size(); ++i) {
            file << "5\n"; // VTK_TRIANGLE
        }

        // Write cell data
        file << "CELL_DATA " << triangles.size() << "\n";
        file << "SCALARS scalar_field float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (const auto &value : solution) {
            file << value << "\n";
        }

        file.close();
    }

    // Function to apply initial condition at the centroids of the triangles
    void applyInitialCondition(double (*initialCondition)(double, double)) {
        for (std::size_t t = 0; t < triangles.size(); ++t) {
            const auto &tri = triangles[t];
            std::size_t i1 = tri.nodeIndices[0];
            std::size_t i2 = tri.nodeIndices[1];
            std::size_t i3 = tri.nodeIndices[2];

            double x1 = nodes[i1].x, y1 = nodes[i1].y;
            double x2 = nodes[i2].x, y2 = nodes[i2].y;
            double x3 = nodes[i3].x, y3 = nodes[i3].y;

            double centroidX = (x1 + x2 + x3) / 3.0;
            double centroidY = (y1 + y2 + y3) / 3.0;

            solution[t] = initialCondition(centroidX, centroidY);
        }
    }

    // Function to print the centroid of each triangle
    void printCentroidsOfAllTriangles() const {
        for (std::size_t t = 0; t < triangles.size(); ++t) {
            const Triangle &tri = triangles[t];
            std::size_t i1 = tri.nodeIndices[0];
            std::size_t i2 = tri.nodeIndices[1];
            std::size_t i3 = tri.nodeIndices[2];

            double x1 = nodes[i1].x, y1 = nodes[i1].y;
            double x2 = nodes[i2].x, y2 = nodes[i2].y;
            double x3 = nodes[i3].x, y3 = nodes[i3].y;

            double centroidX = (x1 + x2 + x3) / 3.0;
            double centroidY = (y1 + y2 + y3) / 3.0;

            std::cout << "Centroid of triangle " << t << ": (" << centroidX << ", " << centroidY << ")" << std::endl;
        }
    }

    // Function to print the midpoint of each face
    void printMidpointsOfAllFaces() const {
        for (std::size_t f = 0; f < faces.size(); ++f) {
            const Face &face = faces[f];
            std::size_t i1 = face.nodeIndices[0];
            std::size_t i2 = face.nodeIndices[1];

            double x1 = nodes[i1].x, y1 = nodes[i1].y;
            double x2 = nodes[i2].x, y2 = nodes[i2].y;

            double midpointX = (x1 + x2) / 2.0;
            double midpointY = (y1 + y2) / 2.0;

            std::cout << "Midpoint of face " << f << ": (" << midpointX << ", " << midpointY << ")" << std::endl;
        }
    }

    // Function to print the unit vector perpendicular to the face pointing to the left triangle
    void printPerpendicularUnitVectors() const {
        for (std::size_t f = 0; f < faces.size(); ++f) {
            const Face &face = faces[f];
            std::size_t i1 = face.nodeIndices[0];
            std::size_t i2 = face.nodeIndices[1];

            double x1 = nodes[i1].x, y1 = nodes[i1].y;
            double x2 = nodes[i2].x, y2 = nodes[i2].y;

            // Midpoint of the face
            double midpointX = (x1 + x2) / 2.0;
            double midpointY = (y1 + y2) / 2.0;

            // Vector from node 1 to node 2
            double dx = x2 - x1;
            double dy = y2 - y1;

            // Perpendicular vector
            double perpX = -dy;
            double perpY = dx;

            // Normalize the perpendicular vector
            double length = std::sqrt(perpX * perpX + perpY * perpY);
            perpX /= length;
            perpY /= length;

            // Determine the centroid of the left triangle
            if (face.leftTriangleIndex != -1) {
                const auto &leftTri = triangles[face.leftTriangleIndex];
                std::size_t li1 = leftTri.nodeIndices[0];
                std::size_t li2 = leftTri.nodeIndices[1];
                std::size_t li3 = leftTri.nodeIndices[2];

                double lx1 = nodes[li1].x, ly1 = nodes[li1].y;
                double lx2 = nodes[li2].x, ly2 = nodes[li2].y;
                double lx3 = nodes[li3].x, ly3 = nodes[li3].y;

                double centroidX = (lx1 + lx2 + lx3) / 3.0;
                double centroidY = (ly1 + ly2 + ly3) / 3.0;

                // Check the direction of the perpendicular vector
                double dotProduct = (centroidX - midpointX) * perpX + (centroidY - midpointY) * perpY;
                if (dotProduct < 0) {
                    // Reverse the direction if it's not pointing to the left triangle
                    perpX = -perpX;
                    perpY = -perpY;
                }

                std::cout << "Perpendicular unit vector of face " << f << " pointing to left triangle: (" << perpX << ", " << perpY << ")" << std::endl;
            }
        }
    }
};

// Initial condition function
double initialCondition(double x, double y) {
    return exp(-10 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
}

int main(int argc, char **argv) {
    try {
        gmsh::initialize();
        gmsh::model::add("unit_square");

        // Define the geometry: a unit square
        double lc = 1; // Decreased characteristic length for more triangles
        gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
        gmsh::model::geo::addPoint(1, 0, 0, lc, 2);
        gmsh::model::geo::addPoint(1, 1, 0, lc, 3);
        gmsh::model::geo::addPoint(0, 1, 0, lc, 4);

        gmsh::model::geo::addLine(1, 2, 1);
        gmsh::model::geo::addLine(2, 3, 2);
        gmsh::model::geo::addLine(3, 4, 3);
        gmsh::model::geo::addLine(4, 1, 4);

        gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
        gmsh::model::geo::addPlaneSurface({1}, 1);

        gmsh::model::geo::synchronize();
        gmsh::model::mesh::generate(2);

        Mesh mesh;
        mesh.readFromGmsh();
        mesh.applyInitialCondition(initialCondition);
        mesh.writeInitialConditionToVTK("initial_condition.vtk");
        mesh.printCentroidsOfAllTriangles();
        mesh.printMidpointsOfAllFaces();
        mesh.printPerpendicularUnitVectors();

        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        gmsh::finalize();
        return -1;
    }

    return 0;
}
