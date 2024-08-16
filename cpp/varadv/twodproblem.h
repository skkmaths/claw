#include <string>
#include "matrix.h"
#include "grid.h"
#define SZERO 0.001
#include <cmath>  // Include for sin and cos functions
#include <chrono>
using namespace std;

class BoundaryConditions {
public:
    // Constructor
    BoundaryConditions(const std::string& left, const std::string& right, 
                       const std::string& bottom, const std::string& top) 
        : left(left), right(right), bottom(bottom), top(top) {}
   BoundaryConditions() 
        : left(""), right(""), bottom(""), top("") {}  // Initializes all strings to empty
std::string left;
std::string right;
std::string bottom;
std::string top;
};

// Class for reservoir problem
class TwoDProblem
{
   public:
      TwoDProblem (int nx_val, int ny_val, double Tfinal, double cfl, unsigned int save_freq, string scheme):
       nx(nx_val), ny(ny_val), Tf(Tfinal), cfl(cfl), save_freq(save_freq), scheme(scheme) {}
      ~TwoDProblem () {};
      void run ();

   private:
      int nx;
      int ny;
      Grid    grid;
      Matrix  sol;
      Matrix  sol_old;
      Matrix  res;
      int fileid ;
      double Tf; 
      double dt;
      double cfl;
      unsigned int save_freq;
      string scheme;
      string ic ;
      double lam_x;
      double lam_y;
      BoundaryConditions bc;
      void make_grid ();
      double initial_data( const double& x, const double& y);
      double xflux(const double& x, const double& y,const double& u);
      double yflux (const double& x, const double& y,const double& u);
      double xnumflux(const double& x, const double& y,const double& state_left,
    	       const double& state_right);
      double ynumflux(const double& x, const double& y, const double& state_left,
    	       const double& state_right);
      void initialize();
      void compute_residual(Matrix& res, const double& time);
      void updateGhostCells ();
      void solve();
      void savesol(double t, Matrix& sol);
      vector<double> findMinMax();
      void compute_error(double& l1error);
      void apply_ssprk2(const double& time);
      void apply_euler(const double& time);
      void compute_dt();
      double reconstruct(const double& sol_ll,const double& sol_l,const double& sol_r);
      double exact(const double& x, const double& y, const double& t);
};