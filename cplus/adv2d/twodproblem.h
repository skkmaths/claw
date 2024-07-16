#include <string>
#include "matrix.h"
#include "grid.h"
#define SZERO 0.001
#include <cmath>  // Include for sin and cos functions
#include <chrono>
using namespace std;
// Class for reservoir problem
class TwoDProblem
{
   public:
      TwoDProblem (int nx_val, int ny_val, double Tfinal, double cfl) : nx(nx_val), ny(ny_val), Tf(Tfinal), cfl(cfl) {}
      ~TwoDProblem () {};
      void run ();

   private:
      int nx;
      int ny;
      unsigned int save_freq;
      Grid    grid;
      Matrix  sol;
      Matrix  res;
      int fileid ;
      double Tf; 
      double dt;
      double cfl;
      double lam_x;
      double lam_y;
      void make_grid ();
      double initial_data( const double& x, const double& y);
      double xflux(const double& u);
      double yflux (const double& u);
      double xnumflux(const double& state_left,
    	       const double& state_right);
      double ynumflux(const double& state_left,
    	       const double& state_right);
      void initialize();
      void compute_residual(Matrix& res);
      void updateGhostCells ();
      void solve();
      void savesol(double t, Matrix& sol);
      vector<double> findMinMax();
      void compute_error(double& l1error);
};