#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>

// Node structure
struct Node {
    std::size_t tag;
    double x, y, z;
    double value;
};

// Triangle structure
struct Triangle {
    std::size_t tag;
    std::vector<std::size_t> nodeTags;
};

// Mesh class
class Mesh {
public:
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;

    void readFromGmsh() {
        std::vector<std::size_t> nodeTags;
        std::vector<double> coord, parametricCoord;
        gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

        nodes.resize(nodeTags.size());
        for (std::size_t i = 0; i < nodeTags.size(); ++i) {
            nodes[i].tag = nodeTags[i];
            nodes[i].x = coord[3 * i];
            nodes[i].y = coord[3 * i + 1];
            nodes[i].z = coord[3 * i + 2];
            nodes[i].value = 0.0; // Initialize with zero value
        }

        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

        for (std::size_t i = 0; i < elementTypes.size(); ++i) {
            if (elementTypes[i] == 2) { // Type 2 corresponds to 2D triangular elements
                for (std::size_t j = 0; j < elementTags[i].size(); ++j) {
                    Triangle tri;
                    tri.tag = elementTags[i][j];
                    tri.nodeTags = {elementNodeTags[i][3 * j], elementNodeTags[i][3 * j + 1], elementNodeTags[i][3 * j + 2]};
                    triangles.push_back(tri);
                }
            }
        }
    }

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
            file << "3 " << tri.nodeTags[0] - 1 << " " << tri.nodeTags[1] - 1 << " " << tri.nodeTags[2] - 1 << "\n";
        }

        // Write cell types
        file << "CELL_TYPES " << triangles.size() << "\n";
        for (std::size_t i = 0; i < triangles.size(); ++i) {
            file << "5\n"; // VTK_TRIANGLE
        }

        // Write point data
        file << "POINT_DATA " << nodes.size() << "\n";
        file << "SCALARS scalar_field float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (const auto &node : nodes) {
            file << node.value << "\n";
        }

        file.close();
    }

    void applyInitialCondition(double (*initialCondition)(double, double)) {
        for (auto &node : nodes) {
            node.value = initialCondition(node.x, node.y);
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
        double lc = 0.2; // Decreased characteristic length for more triangles
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

        // Refine the mesh by setting the characteristic length
        std::vector<std::pair<int, int>> entities;
        gmsh::model::getEntities(entities);
        for (auto &entity : entities) {
            gmsh::model::mesh::setSize({entity}, lc);
        }

        // Generate 2D triangular mesh
        gmsh::model::mesh::generate(2);

        // Read mesh into Mesh object
        Mesh mesh;
        mesh.readFromGmsh();

        // Apply initial condition
        mesh.applyInitialCondition(initialCondition);

        // Write initial condition to VTK file
        mesh.writeInitialConditionToVTK("initial_condition.vtk");

        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }
    return 0;
}
