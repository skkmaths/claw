#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Structure for storing node information
struct Node {
    double x, y;
};

// Structure for storing element information
struct Element {
    std::vector<int> nodeIds;
};

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::open("mesh.geo");

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);
    gmsh::write("example_mesh.msh");

    // Read nodes
    std::vector<size_t> nodeTags;
    std::vector<double> nodeCoords;
    std::vector<double> nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

    std::vector<Node> nodes;
    for (size_t i = 0; i < nodeTags.size(); ++i) {
        nodes.push_back(Node{nodeCoords[3 * i], nodeCoords[3 * i + 1]});
    }

    // Read elements
    std::vector<int> elementTypes;
    std::vector<std::vector<size_t>> elementTags, elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    std::vector<Element> elements;
    for (size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == 2) { // 2D triangles
            for (size_t j = 0; j < elementTags[i].size(); ++j) {
                Element elem;
                for (size_t k = 0; k < 3; ++k) {
                    elem.nodeIds.push_back(elementNodeTags[i][j * 3 + k] - 1); // Convert 1-based to 0-based index
                }
                elements.push_back(elem);
            }
        }
    }

    gmsh::finalize();

    // Write VTK file
    std::ofstream vtkFile("example_mesh.vtk");
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Gmsh output with constant initial condition\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points
    vtkFile << "POINTS " << nodes.size() << " float\n";
    for (const auto &node : nodes) {
        vtkFile << node.x << " " << node.y << " 0.0\n";
    }
    
    // Write cells
    vtkFile << "CELLS " << elements.size() << " " << elements.size() * 4 << "\n";
    for (const auto &elem : elements) {
        vtkFile << "3 " << elem.nodeIds[0] << " " << elem.nodeIds[1] << " " << elem.nodeIds[2] << "\n";
    }

    // Write cell types
    vtkFile << "CELL_TYPES " << elements.size() << "\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        vtkFile << "5\n"; // VTK_TRIANGLE
    }

    // Write cell data with constant initial condition
    vtkFile << "CELL_DATA " << elements.size() << "\n";
    vtkFile << "SCALARS initial_condition float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        vtkFile << "1.0\n"; // Constant value of 1
    }

    vtkFile.close();
    
    std::cout << "VTK file with constant initial condition saved as example_mesh.vtk" << std::endl;

    return 0;
}
