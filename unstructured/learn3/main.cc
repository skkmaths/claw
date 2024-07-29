#include <gmsh.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <functional>

// Node structure to represent a mesh node
struct Node {
    double x, y, z;  // Coordinates of the node
};

// Triangle structure to represent a mesh element
struct Triangle {
    std::vector<std::size_t> nodeIndices; // Indices of the nodes that form the triangle
    double area; // Area of the triangle
    Node centroid; // Centroid of the triangle

    // Constructor to compute area and centroid
    Triangle() : area(0.0), centroid{0.0, 0.0, 0.0} {}

    void computeAreaAndCentroid(const std::vector<Node>& nodes) {
        const auto& p0 = nodes[nodeIndices[0]];
        const auto& p1 = nodes[nodeIndices[1]];
        const auto& p2 = nodes[nodeIndices[2]];

        // Calculate area using the cross product of vectors
        double x0 = p1.x - p0.x;
        double y0 = p1.y - p0.y;
        double x1 = p2.x - p0.x;
        double y1 = p2.y - p0.y;
        area = std::abs(x0 * y1 - y0 * x1) / 2.0;

        // Calculate centroid
        centroid.x = (p0.x + p1.x + p2.x) / 3.0;
        centroid.y = (p0.y + p1.y + p2.y) / 3.0;
        centroid.z = (p0.z + p1.z + p2.z) / 3.0;
    }
};

// Face structure to represent a mesh face
struct Face {
    std::vector<std::size_t> nodeIndices; // Indices of the nodes that form the face
    std::size_t leftTriangle;             // Index of the left triangle
    std::size_t rightTriangle;            // Index of the right triangle (set to -1 if boundary)
    bool isBoundary;                      // Flag to indicate if the face is a boundary face
    Node normalLeft;                     // Normal vector for left triangle
    Node normalRight;                    // Normal vector for right triangle
    Node midpoint;

    Face(const std::vector<std::size_t>& nodes, std::size_t left, std::size_t right)
        : nodeIndices(nodes), leftTriangle(left), rightTriangle(right), isBoundary(right == -1) {}
};

// Custom hash function for std::pair
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ hash2;
    }
};

// Mesh class to store nodes, triangles, and faces
class Mesh {
public:
    std::vector<Node> nodes;              // List of nodes in the mesh
    std::vector<Triangle> triangles;      // List of triangles in the mesh
    std::vector<Face> faces;              // List of faces in the mesh
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t, pair_hash> faceMap; // Mapping from node pairs to face indices

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
                    tri.computeAreaAndCentroid(nodes);
                    triangles.push_back(tri);
                }
            }
        }

        createFaces();
    }

    // Function to create faces from triangles
    void createFaces() {
        for (std::size_t i = 0; i < triangles.size(); ++i) {
           const auto& tri = triangles[i];
             for (std::size_t j = 0; j < 3; ++j) {
                std::pair<std::size_t, std::size_t> faceKey = std::minmax(tri.nodeIndices[j], tri.nodeIndices[(j + 1) % 3]);
                if (faceMap.find(faceKey) == faceMap.end()) {
                // New face
                std::size_t node1 = faceKey.first;
                std::size_t node2 = faceKey.second;
                Node midpoint;
                midpoint.x = (nodes[node1].x + nodes[node2].x) / 2.0;
                midpoint.y = (nodes[node1].y + nodes[node2].y) / 2.0;
                midpoint.z = (nodes[node1].z + nodes[node2].z) / 2.0;
                faces.emplace_back(std::vector<std::size_t>{node1, node2}, i, -1);
                faces.back().midpoint = midpoint;
                faceMap[faceKey] = faces.size() - 1;
            } else {
                // Existing face, update right triangle
                auto& face = faces[faceMap[faceKey]];
                face.rightTriangle = i;
                face.isBoundary = false;
            }
        }
    }

    computeFaceNormals();
}

    // Function to compute normal vectors for faces
    void computeFaceNormals() {
        for (auto& face : faces) {
            const auto& leftTri = triangles[face.leftTriangle];
            Node leftNormal = computeNormalToFace(face.nodeIndices, leftTri.centroid);
            face.normalLeft = leftNormal;
            if (!face.isBoundary) {
                const auto& rightTri = triangles[face.rightTriangle];
                Node rightNormal = computeNormalToFace(face.nodeIndices, rightTri.centroid);
                face.normalRight = rightNormal;
            }
        }
    }

    // Function to compute normal vector to a face from the centroid of a triangle
   Node computeNormalToFace(const std::vector<std::size_t>& faceNodes, const Node& centroid) const {
        const auto& p0 = nodes[faceNodes[0]];
        const auto& p1 = nodes[faceNodes[1]];

        // Vector from centroid to a point on the face
        Node foot ;
        Node vector;
        double xbar = p1.x-p0.x;
        double ybar = p1.y-p0.y;
        double lam = ( xbar*(centroid.x-p1.x) + ybar*(centroid.y-p1.y) )/ (xbar*xbar +ybar*ybar);
        foot.x = p1.x +(p1.x-p0.x)*lam;
        foot.y = p1.y + (p1.y-p0.y)*lam;
        vector.x = centroid.x -foot.x;
        vector.y = centroid.y - foot.y;
        vector.z = centroid.z - foot.y;
        // Calculate the normal (perpendicular to the face) in 2D
        double length = std::sqrt(vector.x * vector.x + vector.y * vector.y);
        if (length != 0) {
            vector.x /= length;
            vector.y /= length;
        }
        /*
        std::cout<<"points"<< "( "<< p0.x <<","<<p0.y<<")"<< " and"<< "( "<< p1.x <<","<<p1.y<<")" <<std::endl;
        std::cout<<"foot" <<  "( "<<foot.x <<","<<foot.y<<")"<<std::endl;
        std::cout<<"vector"<< "( "<< vector.x <<","<<vector.y<<")"<<std::endl;
        */

        return vector;
    }

    // Function to print centroids of all triangles
    void printCentroidsOfTriangles() const {
        for (const auto& tri : triangles) {
            std::cout << "Triangle centroid: (" << tri.centroid.x << ", " << tri.centroid.y << ", " << tri.centroid.z << ")\n";
        }
    }
    // Function to print face information
    // Function to print face information
    void printFaceInfo() const {
        for (const auto& face : faces) {
        std::cout << "Face " << (&face - &faces[0]) << ":\n"; // Calculate the index of the face
        std::cout << "  Boundary: " << (face.isBoundary ? "Yes" : "No") << "\n";
        std::cout << "  Left Triangle Index: " << face.leftTriangle << "\n";
        std::cout << "  Right Triangle Index: " << (face.isBoundary ? "N/A" : std::to_string(face.rightTriangle)) << "\n";
        std::cout << "  Normal Left: (" << face.normalLeft.x << ", " << face.normalLeft.y << ")\n";
        if (!face.isBoundary) {
            std::cout << "  Normal Right: (" << face.normalRight.x << ", " << face.normalRight.y << ")\n";
        }
        std::cout << "  Midpoint: (" << face.midpoint.x << ", " << face.midpoint.y << ", " << face.midpoint.z << ")\n";
    }
}

};

 // Function to write initial condition to a VTK file
    void savesol(const std::string &filename, const Mesh& mesh, std::vector<double> solution) {
        std::ofstream file(filename);
        file << "# vtk DataFile Version 3.0\n";
        file << "Advection Initial Condition\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Write points
        file << "POINTS " << mesh.nodes.size() << " float\n";
        for (const auto &node : mesh.nodes) {
            file << node.x << " " << node.y << " " << node.z << "\n";
        }

        // Write cells
        file << "CELLS " << mesh.triangles.size() << " " << 4 * mesh.triangles.size() << "\n";
        for (const auto &tri : mesh.triangles) {
            file << "3 " << tri.nodeIndices[0] << " " << tri.nodeIndices[1] << " " << tri.nodeIndices[2] << "\n";
        }

        // Write cell types
        file << "CELL_TYPES " << mesh.triangles.size() << "\n";
        for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
            file << "5\n"; // VTK_TRIANGLE
        }

        // Write cell data
        file << "CELL_DATA " << mesh.triangles.size() << "\n";
        file << "SCALARS scalar_field float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (const auto &value : solution) {
            file << value << "\n";
        }

        file.close();
    }

// Initial condition function
double initialCondition(double x, double y) {
    return exp(-30 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
}

// Function to initialize the solution vector
void initialize(const Mesh& mesh, std::vector<double>& solution) {
    solution.resize(mesh.triangles.size());
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        solution[i] = initialCondition(mesh.triangles[i].centroid.x, mesh.triangles[i].centroid.y);
    }
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
        mesh.printCentroidsOfTriangles();
        mesh.printFaceInfo();
        std::vector<double> solution;
        initialize(mesh, solution);
        // Save initial conditions to VTK files
        savesol("solution.vtk", mesh, solution);

        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        gmsh::finalize();
        return -1;
    }

    return 0;
}
