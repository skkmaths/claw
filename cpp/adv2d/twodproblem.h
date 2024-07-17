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
      void apply_ssprk2();
      void apply_euler();
      double reconstruct(const double& sol_ll,const double& sol_l,const double& sol_r);
};