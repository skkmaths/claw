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
    int id;
    // Operator for comparing Node pointers (for unordered_map key)
    bool operator==(const Node& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// Triangle structure to represent a mesh element
struct Triangle {
    int id;
    std::vector<Node*> nodes; // Pointers to the nodes that form the triangle
    double area; // Area of the triangle
    Node centroid; // Centroid of the triangle

    // Constructor to compute area and centroid
    Triangle() : area(0.0), centroid{0.0, 0.0, 0.0} {}

    void computeAreaAndCentroid() {
        const auto& p0 = *nodes[0];
        const auto& p1 = *nodes[1];
        const auto& p2 = *nodes[2];

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
    std::vector<Node*> nodes; // Pointers to the nodes that form the face
    Triangle* leftTriangle;   // Pointer to the left triangle
    Triangle* rightTriangle;  // Pointer to the right triangle (nullptr if boundary)
    bool isBoundary;          // Flag to indicate if the face is a boundary face
    Node normalLeft;          // Normal vector for left triangle
    Node normalRight;         // Normal vector for right triangle
    Node midpoint;

    Face(const std::vector<Node*>& nodes, Triangle* left, Triangle* right)
        : nodes(nodes), leftTriangle(left), rightTriangle(right), isBoundary(right == nullptr) {}
};

// Custom hash function for std::pair of indices
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ (hash2 << 1);
    }
};

// Mesh class to store nodes, triangles, and faces
class Mesh {
public:
    std::vector<Node> nodes;              // List of nodes in the mesh
    std::vector<Triangle> triangles;      // List of triangles in the mesh
    std::vector<Face> faces;              // List of faces in the mesh
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t, pair_hash> faceMap; // Mapping from node indices to face indices

    // Function to read mesh data from Gmsh
    void readFromGmsh() {
        std::vector<std::size_t> nodeTags;
        std::vector<double> coord, parametricCoord;
        gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

        nodes.resize(nodeTags.size());
        for (std::size_t i = 0; i < nodeTags.size(); ++i) {
            nodes[i].id = static_cast<int>(i);
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
                    tri.id = static_cast<int>(triangles.size()); // Assign ID based on current size of triangles
                    tri.nodes = {&nodes[elementNodeTags[i][3 * j] - 1], &nodes[elementNodeTags[i][3 * j + 1] - 1], &nodes[elementNodeTags[i][3 * j + 2] - 1]};
                    tri.computeAreaAndCentroid();
                    triangles.push_back(tri);
                }
            }
        }

        createFaces();
    }

    // Function to create faces from triangles
    void createFaces() {
        for (std::size_t i = 0; i < triangles.size(); ++i) {
            auto& tri = triangles[i];
            for (std::size_t j = 0; j < 3; ++j) {
                std::size_t nodeIndex1 = std::find(nodes.begin(), nodes.end(), *tri.nodes[j]) - nodes.begin();
                std::size_t nodeIndex2 = std::find(nodes.begin(), nodes.end(), *tri.nodes[(j + 1) % 3]) - nodes.begin();
                std::pair<std::size_t, std::size_t> faceKey = std::minmax(nodeIndex1, nodeIndex2);

                if (faceMap.find(faceKey) == faceMap.end()) {
                    // New face
                    Node* node1 = tri.nodes[j];
                    Node* node2 = tri.nodes[(j + 1) % 3];
                    Node midpoint;
                    midpoint.x = (node1->x + node2->x) / 2.0;
                    midpoint.y = (node1->y + node2->y) / 2.0;
                    midpoint.z = (node1->z + node2->z) / 2.0;
                    
                    faces.emplace_back(std::vector<Node*>{node1, node2}, &tri, nullptr);
                    faces.back().midpoint = midpoint;
                    faceMap[faceKey] = faces.size() - 1;
                } else {
                    // Existing face, update the right triangle
                    std::size_t faceId = faceMap[faceKey];
                    faces[faceId].rightTriangle = &tri;
                    faces[faceId].isBoundary = false;
                }
            }
        }

        computeFaceNormals();
    }

    // Function to compute normal vectors for faces
    void computeFaceNormals() {
        for (auto& face : faces) {
            Node leftNormal = computeNormalToFace(face.nodes, face.leftTriangle->centroid);
            face.normalLeft = leftNormal;
            if (!face.isBoundary) {
                Node rightNormal = computeNormalToFace(face.nodes, face.rightTriangle->centroid);
                face.normalRight = rightNormal;
            }
        }
    }

    // Function to compute normal vector to a face from the centroid of a triangle
    Node computeNormalToFace(const std::vector<Node*>& faceNodes, const Node& centroid) const {
        const auto& p0 = *faceNodes[0];
        const auto& p1 = *faceNodes[1];

        // Vector from centroid to a point on the face
        Node foot;
        Node vector;
        double xbar = p1.x - p0.x;
        double ybar = p1.y - p0.y;
        double lam = (xbar * (centroid.x - p1.x) + ybar * (centroid.y - p1.y)) / (xbar * xbar + ybar * ybar);
        foot.x = p1.x + (p1.x - p0.x) * lam;
        foot.y = p1.y + (p1.y - p0.y) * lam;
        vector.x = centroid.x - foot.x;
        vector.y = centroid.y - foot.y;
        vector.z = centroid.z - foot.z;

        // Calculate the normal (perpendicular to the face) in 2D
        double length = std::sqrt(vector.x * vector.x + vector.y * vector.y);
        if (length != 0) {
            vector.x /= length;
            vector.y /= length;
        }

        return vector;
    }
};


// Function to write initial condition to a VTK file
    void savesol(const std::string &filename, const Mesh& mesh, std::vector<double>& solution) {
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
            file << "3 " << tri.nodes[0]->id << " " << tri.nodes[1]->id << " " << tri.nodes[2]->id << "\n";
        }

        // Write cell types
        file << "CELL_TYPES " << mesh.triangles.size() << "\n";
        for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
            file << "5\n"; // VTK_TRIANGLE
        }

        // Write cell data
        file << "CELL_DATA " << mesh.triangles.size() << "\n";
        file << "SCALARS sol float 1\n";
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
// intialize the solution
void initialize(const std::vector<Triangle> triangles, std::vector<double>& solution) {
    solution.resize(triangles.size());
    for (const auto &tri : triangles)
    {
    int i = tri.id;
    solution[i] = initialCondition(tri.centroid.x,tri.centroid.y);
    }
    
}

void compute_residue(const std::vector<double> &sol, std::vector<double> &res, const Mesh &mesh)
{
std::fill(res.begin(), res.end(), 0.0);
}

// Main function
int main() {
    gmsh::initialize();
    gmsh::model::add("t1");

    // Define geometry
    double lc = 0.01;
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

    // Generate mesh
    gmsh::model::mesh::generate(2);

    Mesh mesh;
    mesh.readFromGmsh();
    double dt=0.12;
    double time= 0.0;
    double Tf = 1.0;
    std::vector<double> solution, res;
    res.resize(mesh.triangles.size(),0.0); // intialize to zero
    initialize(mesh.triangles, solution); // initialize solution vector and res
    while (time < Tf)
    {
       if (time + dt >Tf){ dt = Tf-time;}
       compute_residue(solution, res, mesh);
       // update solution
       for (auto& tri : mesh.triangles)
          {
          unsigned int i = tri.id;
          solution[i] = solution[i]-dt*res[i];
          }
        time +=dt;
        std::cout<<" time = "<< time<<std::endl;

    }
    savesol("solution.vtk", mesh, solution);


    gmsh::finalize();
    return 0;
}
