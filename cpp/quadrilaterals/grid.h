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
struct Cell {
    int id;
    std::vector<Node*> nodes; // Pointers to the nodes that form the triangle
    double area; // Area of the triangle
    Node centroid; // Centroid of the triangle
    double perimeter;
    // Constructor to compute area and centroid
    Cell() : area(0.0), centroid{0.0, 0.0, 0.0} {}

    void computeAreaAndCentroid() {
        const auto& p0 = *nodes[0];
        const auto& p1 = *nodes[1];
        const auto& p2 = *nodes[2];
        const auto& p3 = *nodes[3];

        // Calculate area using the cross product of vectors
        area = std::abs( (p0.x*p1.y  + p1.x * p2.y + p2.x * p3.y + p3.x * p0.y ) -  
               (p0.y * p1.x + p1.y * p2.x + p2.y * p3.x + p3.y * p0.x)) / 2.0;


        // Calculate centroid
        centroid.x = (p0.x + p1.x + p2.x + p3.x) / 4.0;
        centroid.y = (p0.y + p1.y + p2.y + p3.y) / 4.0;
        centroid.z = (p0.z + p1.z + p2.z + p3.z) / 4.0;
    }
};

// Face structure to represent a mesh face
struct Face {
    std::vector<Node*> nodes; // Pointers to the nodes that form the face
    Cell* leftCell;   // Pointer to the left triangle
    Cell* rightCell;  // Pointer to the right triangle (nullptr if boundary)
    bool isBoundary;          // Flag to indicate if the face is a boundary face
    Node normalLeft;          // Normal vector for left triangle
    Node normalRight;         // Normal vector for right triangle
    Node midpoint;
    double length;            // Length of the face

    Face(const std::vector<Node*>& nodes, Cell* left, Cell* right)
        : nodes(nodes), leftCell(left), rightCell(right), isBoundary(right == nullptr) {}
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
    std::vector<Cell> cells;      // List of triangles in the mesh
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
                if (elementTypes[i] == 3) { // Type 2 corresponds to 2D quadrilateral elements
                    for (std::size_t j = 0; j < elementTags[i].size(); ++j) {
                        Cell cell;
                        cell.id = static_cast<int>(cells.size()); // Assign ID based on current size of triangles
                        cell.nodes = {&nodes[elementNodeTags[i][4 * j] - 1],
                        &nodes[elementNodeTags[i][4 * j + 1] - 1],
                        &nodes[elementNodeTags[i][4 * j + 2] - 1],
                        &nodes[elementNodeTags[i][4 * j + 3] - 1]};
                        cell.computeAreaAndCentroid();
                        cells.push_back(cell);
                    }
                }
            }

            createFaces();
            perimeterCells();
        } catch (const std::exception &e) {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
            gmsh::finalize();
            throw;
        }
    }
   
    // Function to create faces from triangles
    void createFaces() {
    for (std::size_t i = 0; i < cells.size(); ++i) {
        auto& cell = cells[i];
        for (std::size_t j = 0; j < 4; ++j) {
            std::size_t nodeIndex1 = std::find(nodes.begin(), nodes.end(), *cell.nodes[j]) - nodes.begin();
            std::size_t nodeIndex2 = std::find(nodes.begin(), nodes.end(), *cell.nodes[(j + 1) % 4]) - nodes.begin();
            std::pair<std::size_t, std::size_t> faceKey = std::minmax(nodeIndex1, nodeIndex2);

            if (faceMap.find(faceKey) == faceMap.end()) {
                // New face
                Node* node1 = cell.nodes[j];
                Node* node2 = cell.nodes[(j + 1) % 4];
                Node midpoint;
                midpoint.x = (node1->x + node2->x) / 2.0;
                midpoint.y = (node1->y + node2->y) / 2.0;
                midpoint.z = (node1->z + node2->z) / 2.0;
                
                faces.emplace_back(std::vector<Node*>{node1, node2}, &cell, nullptr);
                faces.back().midpoint = midpoint;
                faceMap[faceKey] = faces.size() - 1;
            } else {
                // Existing face, update the right quadrilateral
                std::size_t faceId = faceMap[faceKey];
                faces[faceId].rightCell = &cell;
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
        void perimeterCells()
    {
        for( auto &cell : cells)
        {
            Node* p1 = cell.nodes[0];
            Node* p2 = cell.nodes[1];
            Node* p3 = cell.nodes[2];

            double p = std::sqrt( std::pow( p1->x -p2->x,2) +std::pow(p1->y-p2->y,2) )+
                       std::sqrt( std::pow( p1->x -p3->x,2) +std::pow(p1->y-p3->y,2) )+
                        std::sqrt( std::pow( p2->x -p3->x,2) +std::pow(p2->y-p3->y,2) );
            cell.perimeter = p;
        }
    }
    // Function to compute normal vectors for faces
    void computeFaceNormals() {
        for (auto& face : faces) {
            Node leftNormal = computeNormalToFace(face.nodes, face.leftCell->centroid);
            face.normalLeft = leftNormal;
            if (!face.isBoundary) {
                Node rightNormal = computeNormalToFace(face.nodes, face.rightCell->centroid);
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

    void printCells() const {
        for (const auto& cell : cells) {
            std::cout << "Cell ID: " << cell.id << "\n";
            std::cout << "Cell perimeter = " << cell.perimeter << "\n";
            std::cout << "Cell area = " << cell.area << "\n";
            std::cout << "Centroid: (" << cell.centroid.x << ", " << cell.centroid.y << ", " << cell.centroid.z << ")\n";
            std::cout << "Vertices:\n";
            for (const auto& node : cell.nodes) {
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
            std::cout << "Left Cell ID: " << (face.leftCell ? face.leftCell->id : -1) << "\n";
            std::cout << "Right Triangle ID: " << (face.rightCell ? face.rightCell->id : -1) << "\n";
            std::cout << "Left Normal: (" << face.normalLeft.x << ", " << face.normalLeft.y << ", " << face.normalLeft.z << ")\n";
            if (!face.isBoundary) {
                std::cout << "Right Normal: (" << face.normalRight.x << ", " << face.normalRight.y << ", " << face.normalRight.z << ")\n";
            }
            std::cout << "Boundary: " << (face.isBoundary ? "Yes" : "No") << "\n";
            std::cout << "\n";
        }
    }
    
};