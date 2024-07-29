#include <gmsh.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>


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
    std::size_t leftTriangle;             // Index of the left triangle
    std::size_t rightTriangle;            // Index of the right triangle (set to -1 if boundary)
    bool isBoundary;                      // Flag to indicate if the face is a boundary face

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
                    faces.emplace_back(std::vector<std::size_t>{faceKey.first, faceKey.second}, i, -1);
                    faceMap[faceKey] = faces.size() - 1;
                } else {
                    // Existing face, update right triangle
                    auto& face = faces[faceMap[faceKey]];
                    face.rightTriangle = i;
                    face.isBoundary = false;
                }
            }
        }
    }

    // Function to print centroids of all triangles
    void printCentroidsOfTriangles() const {
        for (const auto& tri : triangles) {
            double centroidX = 0, centroidY = 0, centroidZ = 0;
            for (const auto& nodeIndex : tri.nodeIndices) {
                centroidX += nodes[nodeIndex].x;
                centroidY += nodes[nodeIndex].y;
                centroidZ += nodes[nodeIndex].z;
            }
            centroidX /= 3;
            centroidY /= 3;
            centroidZ /= 3;
            std::cout << "Triangle centroid: (" << centroidX << ", " << centroidY << ", " << centroidZ << ")\n";
        }
    }

    // Function to print midpoints of all faces and the centroids of the left and right triangles
    void printFaceMidpointsAndTriangleCentroids() const {
        for (const auto& face : faces) {
            double midpointX = (nodes[face.nodeIndices[0]].x + nodes[face.nodeIndices[1]].x) / 2.0;
            double midpointY = (nodes[face.nodeIndices[0]].y + nodes[face.nodeIndices[1]].y) / 2.0;
            double midpointZ = (nodes[face.nodeIndices[0]].z + nodes[face.nodeIndices[1]].z) / 2.0;
            std::cout << "Face midpoint: (" << midpointX << ", " << midpointY << ", " << midpointZ << ")\n";

            const auto& leftTri = triangles[face.leftTriangle];
            double leftCentroidX = 0, leftCentroidY = 0, leftCentroidZ = 0;
            for (const auto& nodeIndex : leftTri.nodeIndices) {
                leftCentroidX += nodes[nodeIndex].x;
                leftCentroidY += nodes[nodeIndex].y;
                leftCentroidZ += nodes[nodeIndex].z;
            }
            leftCentroidX /= 3;
            leftCentroidY /= 3;
            leftCentroidZ /= 3;
            std::cout << "  Left triangle centroid: (" << leftCentroidX << ", " << leftCentroidY << ", " << leftCentroidZ << ")\n";

            if (!face.isBoundary) {
                const auto& rightTri = triangles[face.rightTriangle];
                double rightCentroidX = 0, rightCentroidY = 0, rightCentroidZ = 0;
                for (const auto& nodeIndex : rightTri.nodeIndices) {
                    rightCentroidX += nodes[nodeIndex].x;
                    rightCentroidY += nodes[nodeIndex].y;
                    rightCentroidZ += nodes[nodeIndex].z;
                }
                rightCentroidX /= 3;
                rightCentroidY /= 3;
                rightCentroidZ /= 3;
                std::cout << "  Right triangle centroid: (" << rightCentroidX << ", " << rightCentroidY << ", " << rightCentroidZ << ")\n";
            }
        }
    }
};


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
        mesh.printFaceMidpointsAndTriangleCentroids();
         // Create solution vector and save initial condition at centroids

        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        gmsh::finalize();
        return -1;
    }

    return 0;
}
