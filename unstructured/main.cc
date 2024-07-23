#include <iostream>
#include <vector>
#include <fstream>
#include <gmsh.h>

// Define constants
const double dt = 0.01; // Time step
const int num_steps = 100; // Number of time steps
const double vx = 1.0; // Constant velocity in x
const double vy = 0.0; // Constant velocity in y

struct Cell {
    std::vector<int> nodes; // Node indices
    double area; // Cell area
    double u; // Scalar field
};

struct Node {
    double x, y; // Coordinates
};

// Read the mesh
void readMesh(const std::string &filename, std::vector<Node> &nodes, std::vector<Cell> &cells) {
    gmsh::initialize();
    gmsh::open(filename);
    
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    nodes.resize(nodeTags.size());
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        nodes[i] = {nodeCoords[3 * i], nodeCoords[3 * i + 1]};
    }

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t > > elementTags, elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    for (std::size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == 3) { // Quadrilateral elements
            for (std::size_t j = 0; j < elementTags[i].size(); ++j) {
                Cell cell;
                cell.nodes = {elementNodeTags[i][4 * j], elementNodeTags[i][4 * j + 1], elementNodeTags[i][4 * j + 2], elementNodeTags[i][4 * j + 3]};
                cells.push_back(cell);
            }
        }
    }

    gmsh::finalize();
}

// Calculate cell area (for quadrilateral cells)
double calculateCellArea(const Cell &cell, const std::vector<Node> &nodes) {
    const Node &n0 = nodes[cell.nodes[0] - 1];
    const Node &n1 = nodes[cell.nodes[1] - 1];
    const Node &n2 = nodes[cell.nodes[2] - 1];
    const Node &n3 = nodes[cell.nodes[3] - 1];
    
    return 0.5 * std::abs(n0.x * n1.y + n1.x * n2.y + n2.x * n3.y + n3.x * n0.y -
                          (n1.x * n0.y + n2.x * n1.y + n3.x * n2.y + n0.x * n3.y));
}

// Initialize the scalar field
void initializeField(std::vector<Cell> &cells) {
    for (auto &cell : cells) {
        cell.u = 1.0; // Initial condition
    }
}

// Update the scalar field using finite volume scheme
void updateField(std::vector<Cell> &cells, const std::vector<Node> &nodes) {
    for (auto &cell : cells) {
        // Simple upwind scheme (for demonstration purposes)
        double u_new = cell.u - dt * (vx * (cell.u - cell.u) / cell.area); // Placeholder for actual update
        cell.u = u_new;
    }
}

int main() {
    std::vector<Node> nodes;
    std::vector<Cell> cells;

    readMesh("unitsquare.msh", nodes, cells);

    for (auto &cell : cells) {
        cell.area = calculateCellArea(cell, nodes);
    }

    initializeField(cells);

    for (int step = 0; step < num_steps; ++step) {
        updateField(cells, nodes);
    }

    // Output results
    std::ofstream output("results.txt");
    for (const auto &cell : cells) {
        output << cell.u << std::endl;
    }
    output.close();

    return 0;
}
