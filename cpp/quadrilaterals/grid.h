#include <cassert>

// Boundary struct to store connected boundary segments and their IDs
struct Boundary {
    std::string id;
    std::vector<std::tuple<double, double>> nodes;  // Vector of tuples to store (x, y) coordinates
};

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
    Cell* leftLeftCell;          // Pointer to the left-left cell
    Cell* rightRightCell;        // Pointer to the right-right cell
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

// Function to check if a nodes lies on a line defined by two points
bool isNodesOnLine(const Node& node1, const Node& node2, const std::tuple<double, double>& point1, const std::tuple<double, double>& point2) {
    double x1 = std::get<0>(point1), y1 = std::get<1>(point1);
    double x2 = std::get<0>(point2), y2 = std::get<1>(point2);
    
    // Check if the first node lies on the line
    double lhs1 = (y2 - y1) * (node1.x - x1);
    double rhs1 = (node1.y - y1) * (x2 - x1);
    bool isNode1OnLine = std::abs(lhs1 - rhs1) < 1e-9;
    
    // Check if the second node lies on the line
    double lhs2 = (y2 - y1) * (node2.x - x1);
    double rhs2 = (node2.y - y1) * (x2 - x1);
    bool isNode2OnLine = std::abs(lhs2 - rhs2) < 1e-9;
    
    // Return true only if both nodes lie on the line
    return isNode1OnLine && isNode2OnLine;
}

// Mesh class to store nodes, triangles, and faces
class Mesh {
public:
    std::vector<Node> nodes;              // List of nodes in the mesh
    std::vector<Cell> cells;      // List of triangles in the mesh
    std::vector<Face> faces;              // List of faces in the mesh
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t, pair_hash> faceMap; // Mapping from node indices to face indices
    std::vector<Node*> nodeMap; // Mapping from node ID to corresponding Node pointer
    std::vector<Cell*> cellMap; // Mapping from cell ID to corresponding Cell
    std::vector<Boundary> boundaries;  // To store the boundaries
    std::string get_boundary_id(const Face& face);
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
                //Ensure that the assigned ID is within bounds before adding to the vector
                assert(nodes[i].id < nodes.size() && "Nodes ID is out of bounds.");
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
                        //Ensure that the assigned ID is within bounds before adding to the vector
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
            generateNodeMap();
            generateCellid();
            defineboundaries();
        } catch (const std::exception &e) {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
            gmsh::finalize();
            throw;
        }
    }
    void generateCellid(){
        static int i=0;
        for(auto &cell : cells)
        {
            cell.id = static_cast<int>(i);
            assert(cell.id < cells.size() && "Cell ID is out of bounds.");
            i++;
        }
    }
    // Function to generate the node ID to Node pointer map
    void generateNodeMap() {
        nodeMap.resize(nodes.size());
        for (auto& node : nodes) {
            nodeMap[node.id] = &node;
        }
    }
    // Function to generate the cell ID to Cell pointer map
    void generateCellMap() {
        cellMap.resize(cells.size());
        for (auto& cell : cells) {
            cellMap[cell.id] = &cell;
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
    find_nbr_cells();
}
// Function to find the cell adjacent to the given cell across the face opposite to the given face
    Cell* findAdjacentCell(Cell* cell, const Face* givenFace) {
        if (!cell || !givenFace) {
            std::cerr << "Error: Null cell or face pointer." << std::endl;
            return nullptr;
        }

        for (std::size_t i = 0; i < 4; ++i) {
            Node* node1 = cell->nodes[i];
            Node* node2 = cell->nodes[(i + 1) % 4];
        
           if ((node1 == givenFace->nodes[0] && node2 == givenFace->nodes[1]) ||
                (node1 == givenFace->nodes[1] && node2 == givenFace->nodes[0])) {
                Node* oppositeNode1 = cell->nodes[(i + 2) % 4];
                Node* oppositeNode2 = cell->nodes[(i + 3) % 4];
            
                for (const auto& face : faces) {
                    if (areNodesShared(face.nodes, {oppositeNode1, oppositeNode2})) {
                        return (face.leftCell == cell) ? face.rightCell : face.leftCell;
                    }
                }
            }
        }

        return nullptr;
    }

// Function to find neighboring cells for each face
    void find_nbr_cells() {
        for (auto& face : faces) {
            if (face.leftCell) {
                face.leftLeftCell = findAdjacentCell(face.leftCell, &face);
            } else {
                face.leftLeftCell = nullptr;
            }

            if (face.rightCell) {
                face.rightRightCell = findAdjacentCell(face.rightCell, &face);
            } else {
                face.rightRightCell = nullptr;
            }
        }
    }

// Helper function to check if two sets of nodes are shared (unordered comparison)
bool areNodesShared(const std::vector<Node*>& nodes1, const std::vector<Node*>& nodes2) {
    return (std::find(nodes1.begin(), nodes1.end(), nodes2[0]) != nodes1.end() &&
            std::find(nodes1.begin(), nodes1.end(), nodes2[1]) != nodes1.end());
}

    // To compute face length 
    void computeFaceLength()
    {  for (auto& face : faces)
       {
          Node* p0 = face.nodes[0];
          Node* p1 = face.nodes[1];
          face.length = std::sqrt( std::pow (p0->x-p1->x,2) + std::pow (p0->y - p1->y,2) );
       }
    }
    // To compute perimeter of each cell
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
        Node vnormal;
        double dx = p1.x - p0.x;
        double dy = p1.y - p0.y;
        double length = std::sqrt(dx * dx + dy * dy);
        vnormal.x = dy / length;
        vnormal.y = -dx / length;
        if(  (p0.x-centroid.x ) * vnormal.x + (p0.y - centroid.y) * vnormal.y <0 ){
        vnormal.x *=-1;
        vnormal.y *=-1;
        }
        return vnormal;
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
        std::cout << "Right Cell ID: " << (face.rightCell ? face.rightCell->id : -1) << "\n";
        std::cout << "LeftLeft Cell ID: " << (face.leftLeftCell ? face.leftLeftCell->id : -1) << "\n";
        std::cout << "RightRight Cell ID: " << (face.rightRightCell ? face.rightRightCell->id : -1) << "\n";
        std::cout << "Left Normal: (" << face.normalLeft.x << ", " << face.normalLeft.y << ", " << face.normalLeft.z << ")\n";
        if (!face.isBoundary) {
            std::cout << "Right Normal: (" << face.normalRight.x << ", " << face.normalRight.y << ", " << face.normalRight.z << ")\n";
        }
        std::cout << "Boundary: " << (face.isBoundary ? "Yes" : "No") << "\n";
        std::cout << "\n";
    }
}

    void defineboundaries()
    {
        // Create boundaries with  names
        Boundary boundary;

        boundary.id = "left";
        boundary.nodes.push_back(std::make_tuple(0.0, 0.0)); // Add a node with coordinates (1.0, 2.0)
        boundary.nodes.push_back(std::make_tuple(0.0, 1.0)); // Add another node with coordinates (3.0, 4.0)
        boundaries.push_back(boundary);

        boundary = Boundary(); // Clear the boundary object
        boundary.id = "right";
        boundary.nodes.push_back(std::make_tuple(1.0, 0.0)); // Add a node with coordinates (1.0, 2.0)
        boundary.nodes.push_back(std::make_tuple(1.0, 1.0)); // Add another node with coordinates (3.0, 4.0)
        boundaries.push_back(boundary);

        boundary = Boundary(); // Clear the boundary object
        boundary.id = "bottom";
        boundary.nodes.push_back(std::make_tuple(0.0, 0.0)); // Add a node with coordinates (1.0, 2.0)
        boundary.nodes.push_back(std::make_tuple(1.0, 0.0)); // Add another node with coordinates (3.0, 4.0)
        boundaries.push_back(boundary);

        boundary = Boundary(); // Clear the boundary object
        boundary.id = "top";
        boundary.nodes.push_back(std::make_tuple(1.0, 1.0)); // Add a node with coordinates (1.0, 2.0)
        boundary.nodes.push_back(std::make_tuple(0.0, 1.0)); // Add another node with coordinates (3.0, 4.0)
        boundaries.push_back(boundary);
    }
};
// function to get the boundary id of a given face
std::string Mesh::get_boundary_id(const Face& face)
{   
    Node* p0 = face.nodes[0];
    Node* p1 = face.nodes[1];
    for(const auto &bdry : boundaries)
        if( isNodesOnLine(*p0, *p1, bdry.nodes[0], bdry.nodes[1]) )
        return bdry.id; 
    return "not_on_boundary"; // Return this if no boundary matches
};