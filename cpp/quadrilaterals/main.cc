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

    // Function to read mesh data from GMSH
    void readFromGmsh(const std::string &filename) {
        try {
            gmsh::open(filename);

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
        } catch (const std::exception &e) {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
            gmsh::finalize();
            throw;
        }
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
        foot.x = p1.x + lam * xbar;
        foot.y = p1.y + lam * ybar;
        foot.z = centroid.z;
        vector.x = foot.x - centroid.x;
        vector.y = foot.y - centroid.y;
        vector.z = foot.z - centroid.z;

        // Calculate the normal (perpendicular to the face) in 2D
        double length = std::sqrt(vector.x * vector.x + vector.y * vector.y);
        if (length != 0) {
            vector.x /= length;
            vector.y /= length;
        }

        return vector;
    }

    void printTriangles() const {
        for (const auto& tri : triangles) {
            std::cout << "Triangle ID: " << tri.id << "\n";
            std::cout << "Centroid: (" << tri.centroid.x << ", " << tri.centroid.y << ", " << tri.centroid.z << ")\n";
            std::cout << "Vertices:\n";
            for (const auto& node : tri.nodes) {
                std::cout << "  (" << node->x << ", " << node->y << ", " << node->z << ")\n";
            }
            std::cout << "\n";
        }
    }

    // Function to print details of all faces
    void printFaces() const {
        for (std::size_t i = 0; i < faces.size(); ++i) {
            const auto& face = faces[i];
            std::cout << "Face ID: " << i << "\n";
            std::cout << "Nodes: (" << face.nodes[0]->x << ", " << face.nodes[0]->y << ", " << face.nodes[0]->z << ") and ("
                      << face.nodes[1]->x << ", " << face.nodes[1]->y << ", " << face.nodes[1]->z << ")\n";
            std::cout << "Left Triangle ID: " << (face.leftTriangle ? face.leftTriangle->id : -1) << "\n";
            std::cout << "Right Triangle ID: " << (face.rightTriangle ? face.rightTriangle->id : -1) << "\n";
            std::cout << "Left Normal: (" << face.normalLeft.x << ", " << face.normalLeft.y << ", " << face.normalLeft.z << ")\n";
            if (!face.isBoundary) {
                std::cout << "Right Normal: (" << face.normalRight.x << ", " << face.normalRight.y << ", " << face.normalRight.z << ")\n";
            }
            std::cout << "Boundary: " << (face.isBoundary ? "Yes" : "No") << "\n";
            std::cout << "\n";
        }
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
    double r = std::sqrt(pow(x - 0.2, 2) + pow(y - 0.2, 2));
    if (r < 0.1) return 1.0;
    else return 0.0;
}

// Initialize the solution
void initialize(const std::vector<Triangle> triangles, std::vector<double>& solution) {
    solution.resize(triangles.size());
    for (const auto &tri : triangles) {
        int i = tri.id;
        solution[i] = initialCondition(tri.centroid.x, tri.centroid.y);
    }
}

void compute_residue(const std::vector<double> &sol, std::vector<double> &res, const Mesh &mesh, const double& dt) {
    std::fill(res.begin(), res.end(), 0.0);
    for (const auto& face : mesh.faces) {
        if (!face.isBoundary) { // Check if the face is not a boundary
            Node n = face.normalLeft; // Normal vector pointing towards right triangle
            Triangle* L = face.leftTriangle;
            Triangle* R = face.rightTriangle;
            double length_face = std::sqrt(std::pow(face.nodes[0]->x - face.nodes[1]->x, 2) + std::pow(face.nodes[0]->y - face.nodes[1]->y, 2));
            double theta = std::acos(n.x); // Angle between right normal and positive x axis

            if (n.y > 0.0) theta = std::acos(n.x);
            else theta = 2.0 * M_PI - std::acos(n.x);
            double speed_xi = std::cos(theta) + std::sin(theta); // Speed in the xi direction of the transformed PDE
            double splus = std::max(speed_xi, 0.0);
            double sminus = std::min(speed_xi, 0.0);
            double flux;
            if (n.x + n.y > 0.0) flux = sol[L->id] * (n.x + n.y);
            else flux = sol[R->id] * (n.x + n.y);
            res[L->id] += length_face * flux / L->area;
            res[R->id] -= length_face * flux / R->area;
        }
    }
}

double solmax(const std::vector<double>& vec) {
    auto maxIt = std::max_element(vec.begin(), vec.end());
    if (maxIt != vec.end()) {
        return *maxIt;
    } else {
        throw std::runtime_error("Vector is empty");
    }
}

double solmin(const std::vector<double>& vec) {
    auto minIt = std::min_element(vec.begin(), vec.end());
    if (minIt != vec.end()) {
        return *minIt;
    } else {
        throw std::runtime_error("Vector is empty");
    }
}

// Main function
int main() {
    try {
        gmsh::initialize();
        
        Mesh mesh;
        mesh.readFromGmsh("mesh.msh");
        mesh.printFaces();
        
        double dt = 0.01;
        double time = 0.0;
        double Tf = 0.3;
        double h = 1e20;
        for (auto &tri : mesh.triangles)h = std::min( h, tri.area);   
    
        dt = 0.5* h;
        unsigned int save_freq = 1;
        unsigned int iter = 0;
        std::vector<double> solution, res;
        res.resize(mesh.triangles.size(),0.0);
        initialize(mesh.triangles, solution); // initialize solution vector
        while (time < Tf)
        {
        if (time + dt >Tf){ dt = Tf-time;}
        compute_residue(solution, res, mesh, dt);
        // update solution
        for (auto& tri : mesh.triangles)
          {
          unsigned int i = tri.id;
          //std::cout<<"id"<< i <<std::endl;
          solution[i] = solution[i] - dt*res[i];
          }
        time +=dt;
        iter +=1;
        if (save_freq > 0){       
            //if (iter % save_freq == 0){savesol(time, sol);}
            std::cout << std::left;
            std::cout << "iter = " << std::setw(8) << iter 
            << "time = " << std::setw(10) << time 
            << "Max = " << std::setw(15) << solmax(solution)
            << "Min = " << std::setw(15) << solmin(solution) << std::endl;
        }
        

    }
        savesol("solution.vtk", mesh, solution);
        
        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }

    return 0;
}
