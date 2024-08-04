#include <gmsh.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <sys/stat.h>

int fileid = 0;

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
    double perimeter;
    // Constructor to compute area and centroid
    Triangle() : area(0.0), centroid{0.0, 0.0, 0.0} {}

    void computeAreaAndCentroid() {
        const auto& p0 = *nodes[0];
        const auto& p1 = *nodes[1];
        const auto& p2 = *nodes[2];

        // Calculate area using the cross product of vectors
        area = std::abs( p0.x*( p1.y-p2.y) + p1.x*(p2.y - p0.y) + p2.x*(p0.y-p1.y) ) / 2.0;

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
    double length;            // Length of the face

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
            perimeterTriangle();
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
        computeFaceLength();
    }
    void computeFaceLength()
    {  for (auto& face : faces)
       {
          Node* p0 = face.nodes[0];
          Node* p1 = face.nodes[1];
          face.length = std::sqrt( std::pow (p0->x-p1->x,2) + std::pow (p0->y - p1->y,2) );
       }
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
    void perimeterTriangle()
    {
        for( auto &tri : triangles)
        {
            Node* p1 = tri.nodes[0];
            Node* p2 = tri.nodes[1];
            Node* p3 = tri.nodes[2];

            double p = std::sqrt( std::pow( p1->x -p2->x,2) +std::pow(p1->y-p2->y,2) )+
                       std::sqrt( std::pow( p1->x -p3->x,2) +std::pow(p1->y-p3->y,2) )+
                        std::sqrt( std::pow( p2->x -p3->x,2) +std::pow(p2->y-p3->y,2) );
            tri.perimeter = p;
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
            std::cout << "Triangle perimeter = " << tri.perimeter << "\n";
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
            std::cout << "Face length: " << face.length << "\n";
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
double length_face(Face& face)
{
    Node* p0 = face.nodes[0];
    Node* p1 = face.nodes[1];
    return std::sqrt( std::pow (p0->x-p1->x,2) + std::pow (p0->y - p1->y,2) );
}

double minfacelength(const Mesh& mesh)
{
    double h = 1e20;
    for (auto &face : mesh.faces)
    {
        h = std::min( h, face.length);
    }
    return h;
}
// For file name
void createDirectory(const std::string& dirname) 
{
    struct stat info;
    if(stat(dirname.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        if(mkdir(dirname.c_str(), 0777) != 0) {
            std::cerr << "Error creating directory " << dirname << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Directory " << dirname << " is created" << std::endl;
    }
}
std::string getFilename(const std::string& basename, int id) {
    std::ostringstream oss;
    oss << basename << "_" << std::setfill('0') << std::setw(4) << id << ".vtk";
    return oss.str();
}
// Function to write initial condition to a VTK file
void savesol(const Mesh& mesh, std::vector<double>& solution, const double& t) {
    std::string dirname = "sol";
    createDirectory(dirname);
    if(fileid == 0) {
        std::cout << "The directory \"sol\" is going to be formatted!" << std::endl;
        std::string pattern = "./sol/*";
        // Remove existing files
        system(("rm -f " + pattern).c_str());
    }
    std::string filename = getFilename("sol/sol", fileid);
    std::ofstream file(filename);
    file << "# vtk DataFile Version 3.0\n";
    file << "Advection Solution at time " << t << "\n";
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

    // Write field data for time
    file << "FIELD FieldData 1\n";
    file << "TIME 1 1 float\n";
    file << t << "\n";

    // Write cell data
    file << "CELL_DATA " << mesh.triangles.size() << "\n";
    file << "SCALARS sol float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto &value : solution) {
        file << value << "\n";
    }

    file.close();
    fileid++;
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

// exact the solution
void exact(const std::vector<Triangle> triangles, std::vector<double>& ue, const double& t) {
    ue.resize(triangles.size());
    for (const auto &tri : triangles) {
        int i = tri.id;
        ue[i] = initialCondition(tri.centroid.x-t, tri.centroid.y-t);
    }
}
void compute_residue(const std::vector<double> &sol, std::vector<double> &res, const Mesh &mesh, const double& dt) {
    std::fill(res.begin(), res.end(), 0.0);
    for (const auto& face : mesh.faces) {
        if (!face.isBoundary) { // Check if the face is not a boundary
            Node n = face.normalLeft; // Normal vector pointing towards right triangle
            Triangle* L = face.leftTriangle;
            Triangle* R = face.rightTriangle;
            //double lengthface = std::sqrt(std::pow(face.nodes[0]->x - face.nodes[1]->x, 2) + std::pow(face.nodes[0]->y - face.nodes[1]->y, 2));
            double theta = std::acos(n.x); // Angle between right normal and positive x axis
            if (n.y > 0.0) theta = std::acos(n.x);
            else theta = 2.0 * M_PI - std::acos(n.x);
            //double speed_xi = std::cos(theta) + std::sin(theta); // Speed in the xi direction of the transformed PDE
            //double splus = std::max(speed_xi, 0.0);
            //double sminus = std::min(speed_xi, 0.0);
            double flux;
            //if (n.x + n.y > 0.0) flux = sol[L->id] * (n.x + n.y);
            //else flux = sol[R->id] * (n.x + n.y);
            //double lam = std::max(L->perimeter/L->area, R->perimeter/R->area)*dt;
            flux = 0.5*( (n.x + n.y)*(sol[L->id]+ sol[R->id])-0.001*(sol[R->id] -sol[L->id])/dt );
            res[L->id] += face.length * flux / L->area;
            res[R->id] -= face.length * flux / R->area;
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
        //mesh.printFaces();
        //mesh.printTriangles();
        double area = 0.0;
        for(auto &tri : mesh.triangles)
        area += tri.area;
        std::cout<<"area="<<area<<std::endl;
        double dt = 0.001;
        double time = 0.0;
        double Tf = 0.5;
        double h = 0.0;
        double alpha = 0.0;
        //h = minfacelength(mesh);
        for (auto &tri : mesh.triangles)  alpha= std::max( h, tri.perimeter/tri.area);   
       
        dt = 1/ alpha;
        unsigned int save_freq = 1;
        unsigned int iter = 0;
        std::vector<double> solution, res, ue;
        res.resize(mesh.triangles.size(),0.0);
        initialize(mesh.triangles, solution); // initialize solution vector
        //exact(mesh.triangles, ue, time);
        savesol(mesh, solution, time);
        while (time < Tf)
        {
            if (time + dt >Tf){ dt = Tf-time;}
            compute_residue(solution, res, mesh, dt);
            // update solution
            for (auto& tri : mesh.triangles)
            {
                unsigned int i = tri.id;
                solution[i] = solution[i] - dt * res[i];
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
            if( save_freq >0)
               if ( iter % save_freq == 0) 
               { 
                   savesol(mesh, solution, time);
               }

        }
        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }

    return 0;
}
