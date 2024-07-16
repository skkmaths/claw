#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>  // Include for sin and cos functions
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <iomanip> // for std::setprecision


using namespace std;

int nx  = 100;  // number of cells in the x direction
int ny  = 100; // number of cells in the y direction
double xmax = 1.0; // xlimit
double xmin = 0.0;
double ymax = 1.0; // ylimit
double ymin = 0.0;
double dx =  (xmax -xmin)/nx;
double dy =  ( ymax - ymin)/ny;
int fileid = 0;

double initial_condition(const double& x, const double& y)
 {
    return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y);
 }

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

void savesol(double t, std::vector<std::vector<double>>& var_u) 
{
    std::string dirname = "sol";
    createDirectory(dirname);

    if(fileid == 0) {
        std::cout << "The directory \"sol\" is going to be formatted!" << std::endl;
        std::string response;
        std::cout << "Do You Want To Continue? [y/n] ";
        std::cin >> response;
        if(response != "y") {
            std::cerr << "Execution is terminated" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string pattern = "./sol/*";
        // Remove existing files
        system(("rm -f " + pattern).c_str());
    }

    std::string filename = getFilename("sol/sol", fileid);
    std::ofstream file(filename);

    file << "TITLE = \"Linear advection equation\"" << std::endl;
    file << "VARIABLES = \"x\", \"y\", \"sol\"" << std::endl;
    file << "ZONE STRANDID=1, SOLUTIONTIME=" << t << ", I=" << nx << ", J=" << ny << ", DATAPACKING=POINT" << std::endl;

    for(int j = 2; j < ny + 2; ++j) {
        for(int i = 2; i < nx + 2; ++i) {
            double x = xmin + (i - 2) * dx + 0.5 * dx;
            double y = ymin + (j - 2) * dy + 0.5 * dy;
            file << std::setprecision(8) << std::fixed << x << ", " << y << ", " << var_u[i][j] << std::endl;
        }
    }

    file.close();
    fileid++;
}
double xflux(double x, double y, double v)
{
    return v;
}
double xnumflux(double x, double y, double Fl, double Fr, double vl, double vr)
{
    return 0.5*(Fl + Fr) - 0.5*(vr - vl);
}
double yflux(double x, double y, double v)
{
    return v;
}
double ynumflux(double x, double y, double Gl, double Gr, double vl, double vr)
{
    return 0.5*(Gl + Gr) - 0.5*(vr - vl);
}
// compute residual
std::vector<std::vector<double>> compute_residual(double t, double lam_x, double lam_y, 
                                                   const std::vector<std::vector<double>>& v, 
                                                   std::vector<std::vector<double>>& vres) {
    for (int i = 0; i < nx + 4; ++i) {
        std::fill(vres[i].begin(), vres[i].end(), 0.0);
    }

    // Compute the inter-cell fluxes
    // Loop over interior vertical faces
    #pragma omp parallel for
    for (int i = 1; i < nx + 2; ++i) {  // face between (i,j) and (i+1,j)
        double xf = xmin + (i - 1) * dx;  // x location of this face
        #pragma omp parallel for
        for (int j = 2; j < ny + 2; ++j) {
            double y = ymin + (j - 2) * dy + 0.5 * dy; // center of vertical face
            double vl =  v[i][j];// reconstruct(v[i - 1][j], v[i][j], v[i + 1][j]);
            double vr =  v[i + 1][j]; //reconstruct(v[i + 2][j], v[i + 1][j], v[i][j]);
            double Fl = xflux(xf, y, vl);
            double Fr = xflux(xf, y, vr);
            double Fn = xnumflux(xf, y, Fl, Fr, vl, vr);
            vres[i][j] += lam_x * Fn;
            vres[i + 1][j] -= lam_x * Fn;
        }
    }

    // Loop over interior horizontal faces
    #pragma omp parallel for
    for (int j = 1; j < ny + 2; ++j) {
        double yf = ymin + (j - 1) * dy;
        #pragma omp parallel for
        for (int i = 2; i < nx + 2; ++i) {
            double x = xmin + (i - 2) * dx + 0.5 * dx;
            double vl = v[i][j]; //reconstruct(v[i][j - 1], v[i][j], v[i][j + 1]);
            double vr = v[i][j + 1]; // reconstruct(v[i][j + 2], v[i][j + 1], v[i][j]);
            double Gl = yflux(x, yf, vl);
            double Gr = yflux(x, yf, vr);
            double Gn = ynumflux(x, yf, Gl, Gr, vl, vr);
            vres[i][j] += lam_y * Gn;
            vres[i][j + 1] -= lam_y * Gn;
        }
    }

    return vres;
}
// Periodic boundary condition

void update_ghost_cells(std::vector<std::vector<double>>& u, int nx, int ny) {
    // left ghost cell
    for (int j = 0; j <= ny + 3; ++j) {
        u[0][j] = u[nx][j];
        u[1][j] = u[nx + 1][j];
    }

    // right ghost cell
    for (int j = 0; j <= ny + 3; ++j) {
        u[nx + 3][j] = u[3][j];
        u[nx + 2][j] = u[2][j];
    }

    // bottom ghost cell
    for (int i = 0; i <= nx + 3; ++i) {
        u[i][0] = u[i][ny];
        u[i][1] = u[i][ny + 1];
    }

    // top ghost cell
    for (int i = 0; i <= nx + 3; ++i) {
        u[i][ny + 2] = u[i][2];
        u[i][ny + 3] = u[i][3];
    }
}

// Main program

int main() {

double Tf = 20.0;
int save_freq = 5;
// allocate solution variables
vector<vector<double>> u(ny+4, vector<double>(nx+4));
vector<vector<double>> ures(ny+4, vector<double>(nx+4));


// Set the initial condition by interpolation
   //#pragma omp parallel for
   for (int i = 0; i < nx + 4; ++i) 
      {
        for (int j = 0; j < ny + 4; ++j) 
          {
            double x = xmin + (i - 2) * dx + 0.5 * dx;
            double y = ymin + (j - 2) * dy + 0.5 * dy;
            double val = initial_condition(x, y);
            u[i][j] = val;
           }
       }
   
   double t = 0.0; // initial time
   int iter = 0 ;
   savesol(t,u);
 
   while (t < Tf)
   {
    double dt = 0.01; //
    if (t+dt > Tf)
    {
        dt = Tf-t;
    }
    double lam_x = dt / dx ;
    double lam_y = dt / dy ;
    update_ghost_cells(u, nx, ny);
    ures = compute_residual( t, lam_x, lam_y, u, ures);
    #pragma omp parallel for
    for (size_t i = 0; i < nx+4; ++i) 
    {
        #pragma omp parallel for
        for (size_t j = 0; j < ny+4; ++j)
        {
            u[i][j] -= ures[i][j];
        }
    }
    //u  = u - ures;
    t += dt ;
    iter +=1;
    if ( iter % save_freq == 0)
    {
        savesol(t, u);
    }
    //std::cout<<"time ="<<t <<endl;

   }



    return 0;
}