#include <gmsh.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <sys/stat.h>
#include"grid.h"

int fileid = 0;
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
    file << "CELLS " << mesh.cells.size() << " " << 5 * mesh.cells.size() << "\n";
    for (const auto &cell : mesh.cells) {
        file << "4 " << cell.nodes[0]->id << " " << cell.nodes[1]->id << " " << cell.nodes[2]->id << " " << cell.nodes[3]->id << "\n";
    }

    // Write cell types
    file << "CELL_TYPES " << mesh.cells.size() << "\n";
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        file << "9\n"; // VTK_QUAD
    }

    // Write field data for time
    file << "FIELD FieldData 1\n";
    file << "TIME 1 1 float\n";
    file << t << "\n";

    // Write cell data
    file << "CELL_DATA " << mesh.cells.size() << "\n";
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
void initialize(const std::vector<Cell> cells, std::vector<double>& solution) {
    solution.resize(cells.size());
    for (const auto &cell : cells) {
        int i = cell.id;
        solution[i] = initialCondition(cell.centroid.x, cell.centroid.y);
    }
}

// exact the solution
void exact(const std::vector<Cell> cells, std::vector<double>& ue, const double& t) {
    ue.resize(cells.size());
    for (const auto &cell : cells) {
        int i = cell.id;
        ue[i] = initialCondition(cell.centroid.x-t, cell.centroid.y-t);
    }
}
void compute_residue(const std::vector<double> &sol, std::vector<double> &res, const Mesh &mesh, const double& dt) {
    std::fill(res.begin(), res.end(), 0.0);
    for (const auto& face : mesh.faces) {
        if (!face.isBoundary) { // Check if the face is not a boundary
            Node n = face.normalLeft; // Normal vector pointing towards right triangle
            Cell* L = face.leftCell;
            Cell* R = face.rightCell;
            //double lengthface = std::sqrt(std::pow(face.nodes[0]->x - face.nodes[1]->x, 2) + std::pow(face.nodes[0]->y - face.nodes[1]->y, 2));
            double theta = std::acos(n.x); // Angle between right normal and positive x axis
            if (n.y > 0.0) theta = std::acos(n.x);
            else theta = 2.0 * M_PI - std::acos(n.x);
            //double speed_xi = std::cos(theta) + std::sin(theta); // Speed in the xi direction of the transformed PDE
            //double splus = std::max(speed_xi, 0.0);
            //double sminus = std::min(speed_xi, 0.0);
            double flux;
            if (n.x + n.y > 0.0) flux = sol[L->id] * (n.x + n.y);
            else flux = sol[R->id] * (n.x + n.y);
            //double lam = std::max(L->perimeter/L->area, R->perimeter/R->area)*dt;
            //flux = 0.5*( (n.x + n.y)*(sol[L->id]+ sol[R->id])-(sol[R->id] -sol[L->id])/lam);
            res[L->id] += face.length * flux / L->area;
            res[R->id] -= face.length * flux / R->area;
        }
    }
}
double solmax(const std::vector<double>& vec) {
    double smax = -1e20;
    for(auto &value : vec)
    {
        smax = std::max(smax, value );
    }
    return smax;
}
double solmin(const std::vector<double>& vec) {
     double smin = 1e20;
    for(auto &value : vec)
    {
        smin = std::min(smin, value );
    }
    return smin;
}
// Main function
int main() {
    try {
        gmsh::initialize();
        Mesh mesh;
        mesh.readFromGmsh("mesh.msh");
        //mesh.printFaces();
        //mesh.printCells();
        double dt ;
        double time = 0.0;
        double Tf = 0.5;
        double cfl = 0.5;
        double speed = -1e-20;
        double dx = 0.01;
        double dy = 0.01;
        for(auto &cell : mesh.cells)
        speed = std::max( speed, 1.0/dx + 1.0/dy + 1e-14 );
        dt = cfl/speed;       
        unsigned int save_freq = 1;
        unsigned int iter = 0;
        std::vector<double> solution, res, ue;
        res.resize(mesh.cells.size(),0.0);
        initialize(mesh.cells, solution); // initialize solution vector
        //exact(mesh.triangles, ue, time);
        savesol(mesh, solution, time);
        while (time < Tf)
        {
            if (time + dt >Tf){ dt = Tf-time;}
            compute_residue(solution, res, mesh, dt);
            // update solution
            for (auto& cell : mesh.cells)
            {
                unsigned int i = cell.id;
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
        //savesol(mesh, solution, time);
        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }

    return 0;
}
