#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <gmsh.h>

const double dt = 0.01; // Time step
const int num_steps = 100; // Number of time steps
const double vx = 1.0; // Constant velocity in x
const double vy = 0.0; // Constant velocity in y

struct Node {
    double x, y; // Coordinates
};

struct Face {
    int n1, n2; // Node indices
    double nx, ny; // Normal vector
    double length; // Length of the face
    int leftCell, rightCell; // Indices of the cells sharing this face
};

struct Cell {
    std::vector<int> nodes; // Node indices
    std::vector<int> faces; // Face indices
    double area; // Cell area
    double u; // Scalar field
    double cx, cy; // Cell center coordinates
};

void readMesh(const std::string &filename, std::vector<Node> &nodes, std::vector<Cell> &cells, std::vector<Face> &faces) {
    gmsh::initialize();
    gmsh::open(filename);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    nodes.resize(nodeTags.size());
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        nodes[i].x = nodeCoords[3 * i];
        nodes[i].y = nodeCoords[3 * i + 1];
    }

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags, elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    for (std::size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == 3) { // Quadrilateral elements
            for (std::size_t j = 0; j < elementTags[i].size(); ++j) {
                Cell cell;
                cell.nodes.push_back(static_cast<int>(elementNodeTags[i][4 * j]));
                cell.nodes.push_back(static_cast<int>(elementNodeTags[i][4 * j + 1]));
                cell.nodes.push_back(static_cast<int>(elementNodeTags[i][4 * j + 2]));
                cell.nodes.push_back(static_cast<int>(elementNodeTags[i][4 * j + 3]));
                cells.push_back(cell);
            }
        }
    }

    // Compute face connectivity
    for (auto &cell : cells) {
        for (size_t i = 0; i < 4; ++i) {
            int n1 = cell.nodes[i];
            int n2 = cell.nodes[(i + 1) % 4];
            bool faceFound = false;
            for (auto &face : faces) {
                if ((face.n1 == n1 && face.n2 == n2) || (face.n1 == n2 && face.n2 == n1)) {
                    face.rightCell = static_cast<int>(&cell - &cells[0]);
                    cell.faces.push_back(static_cast<int>(&face - &faces[0]));
                    faceFound = true;
                    break;
                }
            }
            if (!faceFound) {
                Face face = {n1, n2, 0, 0, 0, static_cast<int>(&cell - &cells[0]), -1};
                faces.push_back(face);
                cell.faces.push_back(static_cast<int>(faces.size() - 1));
            }
        }
    }

    // Compute face normals and lengths
    for (auto &face : faces) {
        const Node &n1 = nodes[face.n1 - 1];
        const Node &n2 = nodes[face.n2 - 1];
        face.length = std::sqrt(std::pow(n2.x - n1.x, 2) + std::pow(n2.y - n1.y, 2));
        face.nx = (n2.y - n1.y) / face.length;
        face.ny = -(n2.x - n1.x) / face.length;
    }

    gmsh::finalize();
}

double calculateCellArea(const Cell &cell, const std::vector<Node> &nodes) {
    const Node &n0 = nodes[cell.nodes[0] - 1];
    const Node &n1 = nodes[cell.nodes[1] - 1];
    const Node &n2 = nodes[cell.nodes[2] - 1];
    const Node &n3 = nodes[cell.nodes[3] - 1];

    return 0.5 * std::abs(n0.x * n1.y + n1.x * n2.y + n2.x * n3.y + n3.x * n0.y -
                          (n1.x * n0.y + n2.x * n1.y + n3.x * n2.y + n0.x * n3.y));
}

void calculateCellCenters(std::vector<Cell> &cells, const std::vector<Node> &nodes) {
    for (auto &cell : cells) {
        const Node &n0 = nodes[cell.nodes[0] - 1];
        const Node &n1 = nodes[cell.nodes[1] - 1];
        const Node &n2 = nodes[cell.nodes[2] - 1];
        const Node &n3 = nodes[cell.nodes[3] - 1];
        cell.cx = 0.25 * (n0.x + n1.x + n2.x + n3.x);
        cell.cy = 0.25 * (n0.y + n1.y + n2.y + n3.y);
    }
}

void initializeField(std::vector<Cell> &cells) {
    for (auto &cell : cells) {
        cell.u = 1.0; // Initial condition
    }
}

void computeFluxes(const std::vector<Cell> &cells, const std::vector<Face> &faces, std::vector<double> &fluxes) {
    fluxes.resize(faces.size());
    for (size_t i = 0; i < faces.size(); ++i) {
        const Face &face = faces[i];
        double un_left = cells[face.leftCell].u;
        double un_right = (face.rightCell != -1) ? cells[face.rightCell].u : un_left; // For boundary faces

        double flux = 0.0;
        double vn = vx * face.nx + vy * face.ny;

        if (vn >= 0) {
            flux = vn * un_left;
        } else {
            flux = vn * un_right;
        }

        fluxes[i] = flux * face.length;
    }
}

void updateField(std::vector<Cell> &cells, const std::vector<Face> &faces, const std::vector<double> &fluxes) {
    for (auto &cell : cells) {
        double net_flux = 0.0;
        for (size_t j = 0; j < cell.faces.size(); ++j) {
            int face_index = cell.faces[j];
            const Face &face = faces[face_index];
            if (face.leftCell == &cell - &cells[0]) {
                net_flux -= fluxes[face_index];
            } else {
                net_flux += fluxes[face_index];
            }
        }
        cell.u -= dt * net_flux / cell.area;
    }
}

void writeToTecplot(const std::string &filename, const std::vector<Cell> &cells) {
    std::ofstream file(filename);
    file << "VARIABLES = \"X\", \"Y\", \"U\"\n";
    file << "ZONE N=" << cells.size() << ", E=" << cells.size() << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n";

    for (const auto &cell : cells) {
        file << cell.cx << " " << cell.cy << " " << cell.u << "\n";
    }

    for (const auto &cell : cells) {
        file << cell.nodes[0] << " " << cell.nodes[1] << " " << cell.nodes[2] << " " << cell.nodes[3] << "\n";
    }

    file.close();
}

int main() {
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;

    readMesh("unit_square.msh", nodes, cells, faces);

    for (auto &cell : cells) {
        cell.area = calculateCellArea(cell, nodes);
    }

    calculateCellCenters(cells, nodes);
    initializeField(cells);

    for (int step = 0; step < num_steps; ++step) {
        std::vector<double> fluxes;
        computeFluxes(cells, faces, fluxes);
        updateField(cells, faces, fluxes);
    }

    writeToTecplot("results.plt", cells);

    return 0;
}

