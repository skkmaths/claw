#include "twodproblem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <iomanip> // for std::setprecision
#include <chrono>
#include <algorithm> // for std::min
#include <cmath>
using namespace std;
// You should set the CFL according to your problem.
#define SIGN(a) (((a)<0) ? -1:1)
//------------------------------------------------------------------------------
// Minmod limiter
// minmod ( 2(u0-ul), (ur-ul)/2, 2(ur-u0) )
//------------------------------------------------------------------------------
double minmod (const double& ul, const double& u0, const double& ur)
{
   double result;
   double beta=1.0;

   double db = beta*(u0 - ul);         // backward difference
   double df = beta*(ur - u0);         // forward difference
   double dc = 0.5 * (ur - ul); // central difference

   if ( (db*df > 0.0) & (dc*df >0.0) )
   {
      result = std::min(min(fabs(db), fabs(dc)), fabs(df) );
      result *= SIGN(db);
   }
   else
      result = 0.0;
   return result;
}
double TwoDProblem::reconstruct(const double& sol_ll,const double& sol_l,const double& sol_r)
{   
    if( scheme == "so"){
       double theta = 0.5;
       return  sol_l + 0.5 * 2.0 * theta * minmod(sol_ll, sol_l, sol_r);
    }
    else return sol_l;
}
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
    oss << basename << "_" << std::setfill('0') << std::setw(4) << id << ".plt";
    return oss.str();
}
vector<double> advection_velocity(const double& x, const double& y)
{
    vector<double> v(2);
    v[0]  = -y;
    v[1]  = x;
    return v;
}

//-------------------------------------------------------------------------------
// flux functions
//---------------------------------------------------------------------------------
double TwoDProblem::xflux(const double& x, const double& y, const double& u)
{
	return  advection_velocity(x,y)[0]*u;
}
double TwoDProblem::yflux(const double& x, const double& y, const double& u)
{
	return advection_velocity(x,y)[1]*u;
}
double  TwoDProblem::initial_data( const double& x, const double& y)
{
    // Discontinuous  initial data
    if (ic == "nonsmooth" ) // solid body rotation
    {double r = sqrt( pow(x+0.45,2)+pow(y,2));
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
    {double r = sqrt( pow(x+0.45,2)+pow(y,2));
     if ( r < 0.35){ return 1-r/0.35; }
     else return 0.0;
    }
    else if (ic == "sin"){
       return sin(2.0 * M_PI * x)* sin(2.0 * M_PI * y);
    }
    else if( ic == "expo")
    {
        return  exp(-200.0*( pow(x-0.3,2) + pow(y,2) ));
    }
    else 
    {cout<<" Unknown ic"<<endl;
     abort();
    }

}

// Exact solution for c1 = -y, c2 = x, of ut+c1ux+c2uy = 0
double TwoDProblem::exact(const double& x, const double& y, const double& t) {
   return initial_data( std::cos(t) * x + std::sin(t) * y, -std::sin(t)* x + std::cos(t) * y);
}
// Numerical  flux in the x direction   across vertical wall
double TwoDProblem::xnumflux(const double& x, const double& y, const double& ul,
       const double& ur)
{
    vector<double> v(2);
    v = advection_velocity(x,y);
    double vplus = std::max(0.0,v[0]);
    double vminus = std::min(0.0,v[0]);
    return vplus * ul + vminus * ur;
    /*
    vector<double> v(2);
    v = advection_velocity(x,y);
    return  0.5*( xflux(x,y,ul) + xflux(x,y,ur) - 0.5 *fabs(v[0])*(ur-ul) );
    */
}
// numerical flux in the y direction  across  horizontal wall
double TwoDProblem::ynumflux(const double& x, const double& y, const double& ul,
       const double& ur)
{
    vector<double> v(2);
    v = advection_velocity(x,y);
    double vplus = std::max(0.0, v[1]);
    double vminus = std::min(0.0,v[1]);
    return vplus * ul + vminus * ur;
    /*
    vector<double> v(2);
    v = advection_velocity(x,y);
    return 0.5*( yflux(x,y,ul) + yflux(x,y,ur) - 0.5* fabs(v[1])*(ur-ul));
    */
}

//------------------------------------------------------------------------------
// Create cartesian grid
//------------------------------------------------------------------------------
void TwoDProblem::make_grid ()
{
   grid.nx = nx;
   grid.ny = ny;
   grid.dx = (grid.xmax-grid.xmin)/grid.nx;
   grid.dy = (grid.ymax-grid.ymin)/grid.ny;
   grid.allocate();
   cout<<"Making grid for 2Dclw problem ..." << endl;
   // write cell centers starting from left bottom corner
   for (unsigned int i = 0; i < grid.nx ; ++i) 
    {
        for (unsigned int j = 0; j < grid.ny ; ++j) 
        {
            grid.xc(i,j) = grid.xmin + (i + 0.5 ) * grid.dx;
            grid.yc(i,j) = grid.ymin + (j + 0.5) * grid.dy ;
        }
    }
   // write cell vertices
   for (unsigned int i = 0; i < grid.nx +1 ; ++i) 
    {
        for (unsigned int j = 0; j < grid.ny+1 ; ++j) 
        {
            grid.x(i,j) = grid.xmin + i * grid.dx;
            grid.y(i,j) = grid.ymin + j* grid.dy ;
        }
    }
}
//------------------------------------------------------------------------------
// allocate memory and set initial condition
//------------------------------------------------------------------------------
void TwoDProblem::initialize ()
{ 
   sol.allocate(grid.nx+4, grid.ny+4); // with two ghost cells each side
   // initialize only real cells
   for (unsigned int i = 2; i < grid.nx + 2; ++i) 
      {
        for (unsigned int j = 2; j < grid.ny + 2; ++j) 
          {
            sol(i,j) = initial_data( grid.xc(i-2,j-2), grid.yc(i-2,j-2));
          }
      }        
}
//------------------------------------------------------------------------------
// residual
//------------------------------------------------------------------------------
void TwoDProblem::compute_residual(Matrix& res, const double& time)
{ 
    // set residual to zero; this is a potential point
    res = 0.0;
    lam_x = dt / grid.dx;
    lam_y = dt / grid.dy;
    // Loop over interior vertical faces 
    //#pragma omp parallel for collapse(2)
    for (unsigned int i = 2; i < grid.nx + 1; ++i){// face between (i,j) and (i+1,j)
       for (unsigned int j = 2; j < grid.ny + 2; ++j){
          // get the cordinate of center of cell face
          double x = grid.x(i-1,j-2);
          double y = grid.y(i-1,j-2)+grid.dy/2.0;
          double sl = reconstruct(sol(i-1,j),sol(i,j), sol(i+1,j));
          double sr = reconstruct(sol(i+2,j), sol(i+1,j), sol(i,j));
          double Fn = xnumflux(x, y, sl, sr);
          res(i,j) +=  lam_x * Fn;
          res(i+1,j)-=  lam_x * Fn;
        }
    }
      // Loop over horizontal interior  faces 
    // face between (i,j) and (i,j+1)
    //#pragma omp parallel for collapse(2)
    for (unsigned int j = 2; j < grid.ny + 1; ++j){
       for (unsigned int i = 2; i < grid.nx + 2; ++i){
          double x = grid.x(i-2,j-1)+grid.dx/2.0;
          double y = grid.y(i-2,j-1);
          double sl = reconstruct(sol(i,j-1), sol(i,j), sol(i,j+1));
          double sr = reconstruct(sol(i,j+2), sol(i,j+1), sol(i,j));
          double Gn = ynumflux(x, y, sl,sr);
          res(i,j) += lam_y * Gn;
          res(i,j+1) -= lam_y * Gn;
        }
    }
    // Boundary faces
    // left vertical
    for (unsigned int j = 2; j < grid.ny + 2; ++j){
        int i = 1; // update res only in real cells
        // get coordinates of face center
        double x = grid.x(i-1,j-2);
        double y = grid.y(i-1,j-2)+grid.dy/2.0;
        if(bc.left == "outflow"){
            double sl = sol(i+1,j); 
            double sr = sol(i+1,j);
            double Fn = xnumflux(x, y, sl, sr);
            res(i+1,j)-=  lam_x * Fn;
        }
        else if ( bc.left == "dc"){
            double Fn = xflux( x, y, exact (x, y, time )); // flux at exact solution
            res(i+1,j) -=lam_x * Fn;
        }
        else if (bc.left == "periodic"){
            double sl = reconstruct(sol(i-1,j),sol(i,j), sol(i+1,j));
            double sr = reconstruct(sol(i+2,j), sol(i+1,j), sol(i,j));
            double Fn = xnumflux(x, y, sl, sr);
            res(i+1,j) -=  lam_x * Fn;
        }
        else {std::cout<<"Unkown bc at left face: "<<std::endl; abort();}   
    }
    // right vertical face
    for (unsigned int j = 2; j < grid.ny + 2; ++j){
        int i = grid.nx+1;
        double x = grid.x(i-1,j-2);
        double y = grid.y(i-1,j-2)+grid.dy/2.0;
        if(bc.right == "outflow"){
            double sl = sol(i,j); 
            double sr = sol(i,j); 
            double Fn = xnumflux(x, y, sl, sr);
            res(i,j)+=  lam_x * Fn;
        }
        else if ( bc.right == "dc"){
            double Fn = xflux(x, y, exact (x, y , time ));
            res(i,j) +=lam_x * Fn;
        }
        else if (bc.right == "periodic"){
            double sl = reconstruct(sol(i-1,j),sol(i,j), sol(i+1,j));
            double sr = reconstruct(sol(i+2,j), sol(i+1,j), sol(i,j));
            double Fn = xnumflux(x, y, sl, sr);
            res(i,j) +=  lam_x * Fn;
        }
        else {std::cout<<"Unkown bc at left face: "<<std::endl; abort();}
    }
    // Bottom horizontal faces
    for (unsigned int i = 2; i < grid.nx + 2; ++i){
        int j = 1;
        double x = grid.x(i-2,j-1)+grid.dx/2.0;
        double y = grid.y(i-2,j-1);
        if ( bc.bottom == "outflow"){
            double sl = sol(i,j+1);  
            double sr = sol(i,j+1); 
            double Gn = ynumflux(x, y, sl,sr);
            res(i,j+1) -= lam_y * Gn;
        }
        else if( bc.bottom == "dc"){
            double Gn = yflux(x, y, exact(x, y, time ) );
            res(i,j+1) -= lam_y*Gn;
        }
        else if ( bc.bottom == "periodic")
        {
            double sl = reconstruct(sol(i,j-1), sol(i,j), sol(i,j+1));
            double sr = reconstruct(sol(i,j+2), sol(i,j+1), sol(i,j));
            double Gn = ynumflux(x, y, sl,sr);
            res(i,j+1) -= lam_y * Gn;
        }
        else {
            std::cout<<"Unkown bc at left face: "<<std::endl; abort();
        }
    }
    // top horizontal  face
    for (unsigned int i = 2; i < grid.nx + 2; ++i){
        int j = grid.ny+1;
        double x = grid.x(i-2,j-1)+grid.dx/2.0;
        double y = grid.y(i-2,j-1);
        if ( bc.top == "outflow"){
            double sl = sol(i,j); 
            double sr = sol(i,j); 
            double Gn = ynumflux(x, y, sl,sr);
            res(i,j) += lam_y * Gn;
        }
        else if( bc.top == "dc"){
            double Gn = yflux(x, y, exact(x, y, time ) );
            res(i,j) += lam_y * Gn;
        }
        else if ( bc.top == "periodic")
        {
            double sl = reconstruct(sol(i,j-1), sol(i,j), sol(i,j+1));
            double sr = reconstruct(sol(i,j+2), sol(i,j+1), sol(i,j));
            double Gn = ynumflux(x, y, sl,sr);
            res(i,j) += lam_y * Gn;
        }
        else {
            std::cout<<"Unkown bc at left face: "<<std::endl; abort();
        }
    }
}
//------------------------------------------------------------------------------
// Update solution in ghost cells
// Change this according to your test case
//------------------------------------------------------------------------------
void TwoDProblem::updateGhostCells ()
{ 
    // left vertical 
    for (unsigned int j = 0; j <= grid.ny + 3; ++j) {
        sol(0,j) = sol(grid.nx, j );
        sol(1,j) = sol(grid.nx+1,j) ; 
    }
    // right ghost cell
    for (unsigned int j = 0; j <= grid.ny + 3; ++j) 
    {   sol(grid.nx+3,j) = sol(3,j);
        sol(grid.nx+2,j) = sol(2,j);
    }
    // bottom ghost cell
    for (unsigned int i = 0; i <= grid.nx + 3; ++i) {
        sol(i,0) = sol( i, grid.ny);
        sol(i,1) = sol(i, grid.ny+1);
    }
    // top ghost cell
    for (unsigned int i = 0; i <= grid.nx + 3; ++i) {
        sol(i, grid.ny+2) = sol( i, 2);
        sol(i, grid.ny+3) = sol(i,3);
    }
}
//------------------------------------------------------------------
// Save solution to file in the .dat format for gnuplot purpose.
//------------------------------------------------------------------
void TwoDProblem::savesol(double t, Matrix& sol) 
{   
    std::string dirname = "sol";
    createDirectory(dirname);
    if(fileid == 0) {
        std::cout << "The directory \"sol\" is going to be formatted!" << std::endl;
        /* use this to prevent erasing the existing directory
        std::string response;
        std::cout << "Do You Want To Continue? [y/n] ";
        std::cin >> response;
        if(response != "y") {
            std::cerr << "Execution is terminated" << std::endl;
            exit(EXIT_FAILURE);
        }
        */
        std::string pattern = "./sol/*";
        // Remove existing files
        system(("rm -f " + pattern).c_str());
    }
    std::string filename = getFilename("sol/sol", fileid);
    std::ofstream file(filename);
    file << "TITLE = \"Linear advection equation\"" << std::endl;
    file << "VARIABLES = \"x\", \"y\", \"sol\"" << std::endl;
    file << "ZONE STRANDID=1, SOLUTIONTIME=" << t << ", I=" << grid.nx << ", J=" << grid.ny << ", DATAPACKING=POINT" << std::endl;
    for(unsigned int j = 0; j < grid.ny ; ++j) {
        for(unsigned int i = 0; i < grid.nx ; ++i) {
            file << std::setprecision(8) << std::fixed << grid.xc(i,j) << ", "
             << grid.yc(i,j) << ", " 
             << sol(i+2,j+2) << std::endl;
        }
    }
    file.close();
    fileid++;
}
void TwoDProblem::compute_dt()
{
 double speed = -1.0e20;
 vector<double> v(2);
 for (unsigned int i = 0; i < grid.nx; ++i)
      {
        for (unsigned int j = 0; j < grid.ny; ++j)
        {   
          v = advection_velocity(grid.xc(i,j), grid.yc(i,j));
          speed = max( speed, fabs(v[0])/grid.dx + fabs(v[1])/grid.dy + 1e-14 );
        }
      }
 dt = cfl/speed;
}
std::vector<double> TwoDProblem::findMinMax() 
{
    std::vector<double> maxmin(2);
    double sol_max = -1.0e20;
    double sol_min = 1.0e20;
    //#pragma omp parallel for collapse(2)
    for (unsigned int i = 2; i < grid.nx+2; ++i)
      {
        for (unsigned int j = 2; j < grid.ny+2; ++j)
        {
         sol_max = max(sol(i,j), sol_max);
         sol_min = min( sol(i,j), sol_min);
        }
      }
    maxmin[0] = sol_max;
    maxmin[1] = sol_min;
    return maxmin;
}
//-----------------------------------------------------------------------------
// Compute error in the case of smooth test case
// Final time solution and initial condition
// are assumed to be the same here
//-----------------------------------------------------------------------------
void TwoDProblem::compute_error(double& l1error)
{    l1error = 0.0;
    for (unsigned int i = 2; i < grid.nx+2; ++i)
      {
        for (unsigned int j = 2; j < grid.ny+2; ++j)
        {
            l1error += abs( sol(i,j)- initial_data( grid.xc(i-2,j-2), grid.yc(i-2,j-2) ) );
        } 
      } 
    l1error*=(grid.dx*grid.dy);
}
// ssprk2 time stepping
void TwoDProblem::apply_ssprk2(const double& time){
     sol_old = sol;
     updateGhostCells();
     compute_residual(res,time);
     sol = sol- res;
     updateGhostCells();
     compute_residual(res, time);
     sol = sol - res;
     sol = (sol_old + sol) * 0.5;  
}
// Forward euler time stepping
void TwoDProblem::apply_euler(const double& time){
     updateGhostCells();
     compute_residual(res, time);
     sol = sol- res;
    }
//------------------------------------------------------------------------------
// perform time stepping
//------------------------------------------------------------------------------
void TwoDProblem::solve(){
    double time  = 0.0; // initial time
    unsigned iter = 0;
    fileid = 0;
    compute_dt(); // dt independent of solution, so need to compute only once
    res.allocate(grid.nx+4, grid.ny+4);
    sol_old.allocate(grid.nx+4, grid.ny+4);
	savesol(0.0, sol); // save initial condition
    while (time < Tf)
    {
     if (time + dt >Tf){ dt = Tf-time;}
     //Update time solution
     if (scheme =="so") apply_ssprk2(time); // update solution by rk time stepping 
     else apply_euler(time); // update solution by fo scheme
     time +=dt;
     iter +=1;
     if (save_freq > 0){       
        if (iter % save_freq == 0){savesol(time, sol);}
        std::cout << std::left;
        std::cout << "iter = " << std::setw(8) << iter 
              << "time = " << std::setw(10) << time 
              << "Max = " << std::setw(15) << findMinMax()[0] 
              << "Min = " << std::setw(15) << findMinMax()[1] << std::endl;
     }
    }
    // save final time solution
    savesol(time,sol);
}
//------------------------------------------------------------------------------
// solve the whole problem
//------------------------------------------------------------------------------
void TwoDProblem::run ()
{ 
  // set the domain vertices
  grid.xmax = 1.0;
  grid.xmin = -1.0;
  grid.ymax = 1.0;
  grid.ymin = -1.0;
  // set boundary conditions: choose from "periodic, dc, outflow"
  bc.left = "periodic";
  bc.right = "periodic";
  bc.bottom = "periodic";
  bc.top = "periodic" ;
  // choose initial condition from "expo, sin, smooth, nonsmooth"
  ic =  "expo"; 
  make_grid();
  initialize();
  auto start_wall = std::chrono::system_clock::now();
  solve();
  auto end_wall = std::chrono::system_clock::now();
  std::chrono::duration<double> duration_wall = end_wall - start_wall;
  double l1error;
  compute_error(l1error);
  cout<<"#cells, h, l1error, WCT"<<endl;
  cout<< grid.nx * grid.ny <<" "<< grid.dx <<" "<< l1error<<" "<< duration_wall.count()<< endl;
  // WCT- Wall clock time
}
