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
#include"vis.h"

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
            //double theta = std::acos(n.x); // Angle between right normal and positive x axis
            //if (n.y > 0.0) theta = std::acos(n.x);
            //else theta = 2.0 * M_PI - std::acos(n.x);
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
        double cfl = 0.9;
        double time = 0.0;
        double Tf = 0.5;
        double speed = 0.0;
        //h = minfacelength(mesh);
        for (auto &tri : mesh.triangles)  speed= std::max(speed, tri.perimeter/tri.area);   
        dt = cfl/ speed;
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
