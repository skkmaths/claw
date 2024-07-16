#ifndef __RESERVOIR_H__
#define __RESERVOIR_H__

#include <string>
#include "matrix.h"
#include "grid.h"
#include "quadrature.h"
#define SZERO 0.001

using namespace std;
// Class for reservoir problem
class TwoDProblem
{
   public:
      TwoDProblem () {};
      ~TwoDProblem () {};
      void run (unsigned int& j, double& L1error);

   private:
      unsigned int max_iter;
      unsigned int save_freq;
      unsigned int nrk;
      unsigned int order,file;
      std::string  flux_type;
      std::string  test_case;
      double  ark[3], brk[3];
      double  cfl, final_time, dt, time, current_time;
      double  min_velocity;
      double  max_velocity;
      double  c1inlet;
      double  c2inlet;
      Grid    grid;
      Matrix  solution;
      Matrix  c1,c2;
      Matrix  pressure;
      Matrix  permeability;
      double  abs_max_sol;
      

      void read_input (unsigned int& n);
      void make_grid ();
      void initialize ();
      void residual (const unsigned int&, Matrix&);
      void solve (double& l1error);
      void output (const unsigned int, const double  l1error ) const;
      void output_dat (const unsigned int, const double  l1error ) const;

      double initial_data( const double& x, const double& y);
      double exact_sol(const double& x, const double& y,const double& t);
      double reconstruct
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const;

      double reconstruct1
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const;

      double reconstruct2
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const;

     double num_flux_x(const double& state_left,
    	       const double& state_right);

     double num_flux_y(const double& state_left,
    	       const double& state_right);


      void updateGhostCells (const double& time);
      void findMinMax () const;
      void  compute_maximum();

      void compute_error( double& l1error);

};



double minmod (const double& ul, const double& u0, const double& ur);

std::vector<double> dflu_flux
       (
       const double& velocity,
       const std::vector<double>& state_left,
       const std::vector<double>& state_right,
       const double& g
       );
       
std::vector<double> um_flux
       (
       const double& velocity,
       const std::vector<double>& state_left,
       const std::vector<double>& state_right,
       const double& g
       );
       
double argmin_flux(const double& c1,const double& c2,
                   const double& permeability,
                   const double& velocity);
#endif
