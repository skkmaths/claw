#include <string>
#include "matrix.h"
#include "grid.h"
#define SZERO 0.001
#include <cmath>  // Include for sin and cos functions

using namespace std;
// Class for reservoir problem
class TwoDProblem
{
   public:
      TwoDProblem () {};
      ~TwoDProblem () {};
      void run ();

   private:
      unsigned int save_freq;
      Grid    grid;
      Matrix  sol;

      void make_grid ();

      double initial_data( const double& x, const double& y);
      double xflux(const double& u);
      double yflux (const double& u);

      double xnumflux(const double& state_left,
    	       const double& state_right);

      double ynumflux(const double& state_left,
    	       const double& state_right);
      void initialize();
      void compute_residual(Matrix& ures);
      void updateGhostCells (const double& time);
      void solve();
      void savesol(double t, Matrix& sol);
};