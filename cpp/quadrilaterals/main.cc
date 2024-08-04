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

// Main function
int main() {
    try {
        gmsh::initialize();
        Mesh mesh;
        mesh.readFromGmsh("mesh.msh");
        mesh.printFaces();
        mesh.printCells();
        std::vector<double> solution(mesh.cells.size(),0.0);
        savesol(mesh, solution, 0.0);
        /*
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
       
        dt = 0.8/ alpha;
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
       
        } */
        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }

    return 0;
}
