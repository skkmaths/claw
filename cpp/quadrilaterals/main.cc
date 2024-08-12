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
// Define the flux type
std::string flux_type = "lf";
// Advection velocity
Node velocity(const double& x, const  double& y)
{
    Node v;
    v.x = -y;
    v.y = x;
    return v;
}
double initialCondition(const double& x,const  double& y){
    std::string ic = "smooth";
    if (ic == "nonsmooth" ) // solid body rotation non smooth
    {    double r = sqrt( pow(x+0.45,2)+pow(y,2));
         if ( x> 0.1 & x<0.6 & y>-0.25 & y< 0.25  ) 
         {
            return 1.0;
         }
         else if ( r < 0.35)
         {
            return 1-r/0.35;
         }
    else return 0.0;
    }
    // Smooth initial data
    else if ( ic == "smooth") // solid body rotation
    {   double r = sqrt( pow(x-0.3,2)+pow(y-0.3,2));
        if ( r < 0.2){ return 1-r/0.2; }
        else return 0.0;
    }
    else if( ic == "expo")
    {
        return exp(-100.0*( pow(x+0.3,2) + pow(y+0.3,2) ));
    }
    else 
    {    std::cout<<" Unknown ic"<<std::endl;
         abort();
    }
}
// Initialize the solution
void initialize(const std::vector<Cell> cells, std::vector<double>& solution) {
    solution.resize(cells.size());
    for (const auto &cell : cells) {
        int i = cell.id;
        solution[i] = initialCondition(cell.centroid.x, cell.centroid.y);
    }
}
// Exact solution for c1 = y, c2 = -x, of ut+c1ux+c2uy = 0
double exactvaradv(const Node& p, const double& t) {
   return initialCondition( std::cos(t) * p.x + std::sin(t) * p.y, -std::sin(t)* p.x + std::cos(t) *p.y);
}
// exact the solution
void exact(const std::vector<Cell> cells, std::vector<double>& ue, const double& t) {
    ue.resize(cells.size());
    for (const auto &cell : cells) {
        int i = cell.id;
        ue[i] = initialCondition(cell.centroid.x-t, cell.centroid.y-t);
    }
}
// Compute residue
void compute_residue(const std::vector<double> &sol, std::vector<double> &res,  Mesh &mesh, const double& dt, const double& time) {
    std::fill(res.begin(), res.end(), 0.0);
    for (const auto& face : mesh.faces) {
        
        double flux;
        Node n = face.normalLeft; // Outward unit normal to leftcell
        Cell* L = face.leftCell;
        Cell* R = face.rightCell;
        Node vel = velocity(face.midpoint.x, face.midpoint.y); // Advection velocity at face mid point
        double velnormal = vel.x * n.x + vel.y * n.y;
        double speed = std::abs(velnormal);
        if (!face.isBoundary) { // Check if the face is not a boundary
            if ( flux_type == "lf") 
            flux = 0.5 * ( velnormal * (sol[L->id] + sol[R->id] ) - (sol[R->id] - sol[L->id])*speed); // Lax-Friedrich Flux
            else if( flux_type == "upwind")
            flux = (velnormal>0)? velnormal * sol[L->id] : velnormal * sol[R->id];
            else {std::cout<<"Unknown flux type"<<std::endl; abort();}
            res[L->id] += face.length * flux / L->area;
            res[R->id] -= face.length * flux / R->area;
        }
        else if(face.isBoundary)
        {   std::string facebdryid = mesh.get_boundary_id(face);
            if(facebdryid == "left")
            {
                flux =  velnormal * sol[L->id]; 
                res[L->id] += face.length * flux / L->area;
            }
            if(facebdryid == "bottom")
            {
                flux =  velnormal * exactvaradv(face.midpoint, time); 
                res[L->id] += face.length * flux / L->area;
            }
        }
        else 
        {
            std::cout<<"Face type is not known !"<<std::endl;
             abort();
        }
    }
}
// To find max of solution vector
double solmax(const std::vector<double>& vec) {
    double smax = -1e20;
    for(auto &value : vec)
    {
        smax = std::max(smax, value );
    }
    return smax;
}
// To find min of solution vector
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
        std::cout<<"Reading mesh....."<<std::endl;
        mesh.readFromGmsh("mesh.msh");
        std::cout<<"Reading mesh completed"<<std::endl;
        double dt ;
        double time = 0.0;
        double Tf = 2.0*M_PI; // final time
        double cfl = 0.9;
        double speed = -1e-20;
        for(auto &cell : mesh.cells)
        {   Node vel = velocity(cell.centroid.x, cell.centroid.y);
            speed= std::max(speed, std::sqrt(vel.x*vel.x + vel.y*vel.y )* cell.perimeter/cell.area);   
        }
        dt = cfl/speed;       
        unsigned int save_freq = 10;
        unsigned int iter = 0;
        std::vector<double> solution, res, ue;
        res.resize(mesh.cells.size(),0.0);
        initialize(mesh.cells, solution); // initialize solution vector
        //exact(mesh.triangles, ue, time);
        savesol(mesh, solution, time);
        while (time < Tf)
        {
            if (time + dt >Tf){ dt = Tf-time;}
            compute_residue(solution, res, mesh, dt, time);
            // update solution
            for (auto& cell : mesh.cells)
            {
                unsigned int i = cell.id;
                solution[i] = solution[i] - dt * res[i];
            }
            time +=dt;
            iter +=1;
            if (save_freq > 0){       
                if (iter % save_freq == 0){savesol(mesh, solution, time);}
                std::cout << std::left;
                std::cout << "iter = " << std::setw(8) << iter 
                << "time = " << std::setw(10) << time 
                << "Max = " << std::setw(15) << solmax(solution)
                << "Min = " << std::setw(15) << solmin(solution) << std::endl;
            }
        } 
        std::cout<<"Total number of cells = "<<mesh.cells.size()<<std::endl;
        savesol(mesh, solution, time);
        gmsh::finalize();
    } catch (const std::exception &e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }
    return 0;
}
