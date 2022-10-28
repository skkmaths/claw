#include "twodproblem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
// You should set the CFL according to your problem.


#define SIGN(a) (((a)<0) ? -1:1)
extern int n_interior_min;
const double beta = 2.0; // factor in minmod limiter
using namespace std;
unsigned int  n_qpoints=3;
unsigned  int trape_points=15;

double vel_x=-2.0;
double vel_y=-3.0;
//------------------------------------------------------------------------------
// Minmod limiter
// minmod ( 2(u0-ul), (ur-ul)/2, 2(ur-u0) )
//------------------------------------------------------------------------------
double minmod (const double& ul, const double& u0, const double& ur)
{
   double result;
   double theta=1.0;


   /*
   double db = theta*(u0 - ul);         // backward difference
   double df = theta*(ur - u0);         // forward difference
   double dc = 0.5 * (ur - ul); // central difference

   if (db*dc > 0.0 && dc*df > 0.0)
   {
      result = min( min(fabs(beta*db), fabs(dc)), fabs(beta*df) );
      result *= SIGN(db);
   }
   else
      result = 0.0;
   */

   double db = theta*(u0 - ul);         // backward difference
   double df = theta*(ur - u0);         // forward difference

   if (db*df > 0.0 )
   {
      result = min(fabs(db), fabs(df) );
      result *= SIGN(db);
   }
   else
      result = 0.0;
   return result;
}

//------------------------------------------------------------------------------
// Find left state at interface between (il,jl) and (ir,jr)
// (ill,jll) is to the left of (il,jl)
//------------------------------------------------------------------------------
double TwoDProblem::reconstruct
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const
{
   if (order==1)
      return reconstruct1(ill, jll, il, jl, ir, jr);
   else
      return reconstruct2(ill, jll, il, jl, ir, jr);
}

//------------------------------------------------------------------------------
// First order reconstruction
//------------------------------------------------------------------------------
double TwoDProblem::reconstruct1
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const
{
   double state;

   state= solution(il, jl);
   return state;
}
//------------------------------------------------------------------------------
// Second order reconstruction
//------------------------------------------------------------------------------
double TwoDProblem::reconstruct2
       (
       const unsigned int& ill,
       const unsigned int& jll,
       const unsigned int& il,
       const unsigned int& jl,
       const unsigned int& ir,
       const unsigned int& jr
       ) const
{
   double state;
   double ds;
   // solution
   ds = minmod (solution(ill,jll),
                solution(il,jl),
                solution(ir,jr));
   state= solution (il,jl) + 0.5*ds;
   return state;
}
//-------------------------------------------------------------------------------
// flux functions
//---------------------------------------------------------------------------------
double f(const double& u)
{
	return  u;
}
double g(const double& u)
{
	return u;
}

// source function  h  of  u_t+f(u)_x+g(u)_y=h(x,y,t);

double source_function(const double& x, const double& y, const double& t)
{


//	return cos(2.0*M_PI*(x+t))*2.0*M_PI*2.0
//			+cos(2.0*M_PI*(y+t))*2.0*M_PI*2.0;

	double  c=2.0*M_PI;


	return  c*(cos(c*(x+t))+cos(c*(y+t)))
			+c*cos(c*(x+t))+c*cos(c*(y+t)) ;


}
// Exact solution
double TwoDProblem::exact_sol(const double& x, const double& y,const double& t)
{
	double  value;
	//value=sin(2.0*M_PI*(x-t))+cos(2.0*M_PI*(y-t));
	value=sin(2.0*M_PI*(x+t))+sin(2.0*M_PI*(y+t));
	return value;
}
double  TwoDProblem::initial_data( const double& x, const double& y)
{
	double value=0.0;

    value=exact_sol(x,y,0.0);


	return value;
}

// Numerical  flux in the x direction   across vertical wall
double TwoDProblem::num_flux_x( const double& left_state,
       const double& right_state)
{
 double flux;

//LF  flux
// double lambda=max(abs(left_state),abs(right_state));
 //flux=0.5*(f(left_state)+f(right_state)-lambda*(right_state-left_state));
 return left_state;
}


// numerical flux in the y direction  across  horizontal wall
double TwoDProblem::num_flux_y( const double& left_state,
       const double& right_state)
{
//	 double flux;
	//LF  flux
//	 double lambda=max(abs(left_state),abs(right_state));

//	 flux=0.5*(g(left_state)+g(right_state)-lambda*(right_state-left_state));
return left_state;

}

//------------------------------------------------------------------------------
// Read some input from file
//------------------------------------------------------------------------------
void TwoDProblem::read_input (unsigned int& step)
{
   cout << "Reading input from file data.in ..." << endl;

   ifstream inp;
   string input;

  inp.open ("data.in");
/*
   if(step==0)
   inp.open ("data0.in");
   else if(step==1)
   inp.open("data1.in");
   else if(step==2)
   inp.open("data2.in");
   else if(step==3)
   inp.open("data3.in");
   else
   {
    std::cout<<"Number cycles exceeds the limit check the read_input  function"<<endl;
    abort();
   }
*/


  inp >> input >> test_case;
  assert (input=="testcase");
  assert (test_case=="smooth" || test_case=="nonsmooth" );


  
   inp >> input >> flux_type;      
   assert (input=="flux");
   assert (flux_type=="dflu" || flux_type=="um" );

   inp >> input >> order;      
   assert(input == "order");
   assert(order == 1 || order == 2);

   inp >> input >> max_iter;   
   assert(input == "max_iter");
   assert(max_iter > 0);

   inp >> input >> save_freq;   
   assert(input == "save_freq");
   assert(save_freq > 0);

   inp >> input >> cfl;        
   assert(input == "cfl");
   assert(cfl > 0.0 && cfl <= 1.0);

   inp >> input >> final_time;
   assert(input == "final_time");


   inp >> input >> grid.xmin >> grid.xmax;
   assert(input == "xrange");

   inp >> input >> grid.ymin >> grid.ymax;
   assert(input == "yrange");

   inp >> grid.nx >> grid.ny;
   inp >> grid.n_boundary;
   // This part of the code is for smooth test case
    // it reduces the mesh size by two, four,etch according to the variable step.

     if(test_case=="smooth")
     {
     grid.nx=(grid.nx-1)*pow(2,step)+1;
     grid.ny=(grid.ny-1)*pow(2,step)+1;
     }


   // allocate memory for grid
   grid.allocate ();

   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      inp >> grid.ibeg[n] >> grid.iend[n]
          >> grid.jbeg[n] >> grid.jend[n]
          >> grid.boundary_condition[n];
      assert (grid.ibeg[n] >= 1 && grid.ibeg[n] <= grid.nx);
      assert (grid.iend[n] >= 1 && grid.iend[n] <= grid.nx);
      assert (grid.jbeg[n] >= 1 && grid.jbeg[n] <= grid.ny);
      assert (grid.jend[n] >= 1 && grid.jend[n] <= grid.ny);
      if(grid.ibeg[n] == grid.iend[n] &&
         grid.jbeg[n] == grid.jend[n])
      {
         cout << "Boundary " << n 
              << " is not a surface !!!" << endl;
         abort ();
      }
   }


   // it  redefines the boundary according to step size
   if(test_case=="smooth")
   {
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {

       grid.ibeg[n]=(grid.ibeg[n]-1)*pow(2,step)+1;
       grid.iend[n]=(grid.iend[n]-1)*pow(2,step)+1;
       grid.jbeg[n]=(grid.jbeg[n]-1)*pow(2,step)+1;
       grid.jend[n]=(grid.jend[n]-1)*pow(2,step)+1;
   }
   }



   if(test_case!="smooth"&&step>1)
   {
   cout<<"For non smooth test case order cannot be computed, set n_cycles to 1"<<endl;
   abort();

   }

/*
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if( grid.ibeg[n]!= grid.iend[n])
    	  grid.iend[n]=(grid.nx-1)*(pow(2,step))+1;
       if(grid.jbeg[n]!= grid.jend[n])
       grid.jend[n]=(grid.ny-1)*(pow(2,step))+1;

   }

   grid.nx=(grid.nx-1)*(pow(2,step))+1;
   grid.ny=(grid.ny-1)*(pow(2,step))+1;
*/
   inp.close ();

   cout << "Scheme order           = " << order << endl;
   cout << "Max no. of time steps  = " << max_iter << endl;
   cout << "Solution save freq.    = " << save_freq << endl;
   cout << "CFL number             = " << cfl << endl;
   cout << "nx x ny                = " << grid.nx << " x " << grid.ny << endl;
   cout << "Number of boundaries   = " << grid.n_boundary << endl;
   cout << "Number of cells        = " << grid.n_cells << endl;
   cout << "Number of actual cells = " << (grid.nx-1)*(grid.ny-1) << endl;

}

//------------------------------------------------------------------------------
// Create cartesian grid
//------------------------------------------------------------------------------
void TwoDProblem::make_grid ()
{
   cout << "Making grid for 2Dclw problem ..." << endl;

   // set location of boundary
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      if(grid.ibeg[n] == grid.iend[n])
      {
         if(grid.ibeg[n] == 0)
            grid.b_type[n] = imin;
         else
            grid.b_type[n] = imax;
      }

      if(grid.jbeg[n] == grid.jend[n])
      {
         if(grid.jbeg[n] == 0)
            grid.b_type[n] = jmin;
         else
            grid.b_type[n] = jmax;
      }
   }

   grid.dx   = (grid.xmax - grid.xmin)/(grid.nx - 1);
   grid.dy   = (grid.ymax - grid.ymin)/(grid.ny - 1);

   // This is implicitly assumed in the flux computations
   assert (grid.dx == grid.dy);

   // grid vertex coordinates
   for(unsigned int i=0; i<=grid.nx+1; ++i)
      for(unsigned int j=0; j<=grid.ny+1; ++j)
      {
         grid.x (i,j) = grid.xmin + (i-1) * grid.dx;
         grid.y (i,j) = grid.ymin + (j-1) * grid.dy;
      }

   // cell center coordinates
   for(unsigned int i=0; i<=grid.nx; ++i)
      for(unsigned int j=0; j<=grid.ny; ++j)
      {
         grid.xc (i,j) = 0.25 * ( grid.x(i,j)     + grid.x(i+1,j) + 
                                  grid.x(i+1,j+1) + grid.x(i,j+1) );
         grid.yc (i,j) = 0.25 * ( grid.y(i,j)     + grid.y(i+1,j) + 
                                  grid.y(i+1,j+1) + grid.y(i,j+1) );
      }

   // Ghost cell centers   are  found  correctly in the following
   for(unsigned int i=1; i<=grid.nx-1; ++i)
   {

	   grid.xc(i,0)=grid.xc(i,1);
	   grid.yc(i,0)=grid.yc(i,1)-grid.dy;

	   grid.xc(i,grid.nx)=grid.xc(i,1);
	   grid.yc(i,grid.nx)=grid.yc(i,grid.nx-1)+grid.dy;
   }
   for(unsigned int j=1; j<=grid.nx-1; ++j)
   {

	   grid.xc(grid.nx,j)=grid.xc(grid.dx-1,j)+grid.dx;
	   grid.yc(grid.nx,j)=grid.yc(1,j);

	   grid.xc(0,j)=grid.xc(1,j)-grid.dx;
	   grid.yc(0,j)=grid.yc(1,j);

   }
}

//------------------------------------------------------------------------------
// allocate memory and set initial condition
// Change the following function according to your test case
//------------------------------------------------------------------------------
void TwoDProblem::initialize ()
{
   solution.allocate    (grid.nx+1, grid.ny+1);
   // initialize only real cells, not for ghost cells
   for(unsigned int i=1; i<=grid.nx-1; ++i)
      for(unsigned int j=1; j<=grid.ny-1; ++j)
      {
         //-----------------
    	//-----------------------------------
    	  //for non smooth test case in [0 1]X[0 1]
        if(test_case=="smooth")
        	 solution (i,j) = exact_sol(grid.xc(i,j),grid.yc(i,j),0.0);
        else
        {
    	  double dist = (grid.xc(i,j)-0.5) * (grid.xc(i,j)-0.5) +
    	                (grid.yc(i,j)-0.5) * (grid.yc(i,j)-0.5);
    	   if(sqrt(dist)<= 0.2)
    	   solution(i,j)=2.0;
    	   else
    	   solution(i,j)=0.0;
         //----------------------------------------------------
        }
       // for smooth and exact solution known test case follow the below one



      }

   updateGhostCells (0.0);
   output (0,0.0);

   // RK scheme parameters
   if( order==1 )
      nrk = 1;
   else
      nrk = 3;
   ark[0] = 0.0; ark[1] = 3.0/4.0; ark[2] = 1.0/3.0;
   for (unsigned int i=0; i<nrk; ++i) brk[i] = 1.0 - ark[i];
}
//------------------------------------------------------------------------------
// residual
//------------------------------------------------------------------------------
void TwoDProblem::residual (const unsigned int& rk,
                                 Matrix& s_residual)
{ 

   Matrix residual_dummy (grid.nx+1, grid.ny+1);
   unsigned int i, j;
   double velocity;
   double state_left, state_right, flux;
   s_residual = 0.0;
   min_velocity = 1.0e20;
   max_velocity = 0.0;
   // interior vertical faces

   for(unsigned int i=2; i<=grid.nx-1; ++i)

      for(unsigned int j=1; j<=grid.ny-1; ++j)
      {
         state_left  = reconstruct (i-2, j, i-1, j, i, j);

         state_right = reconstruct (i+1, j, i, j, i-1, j);

          flux=num_flux_x(state_left,state_right);

         s_residual (i-1,j) +=flux*grid.dy;
         s_residual (i,  j) -=flux*grid.dy;
      }

   // interior horizontal faces
   for(unsigned int j=2; j<=grid.ny-1; ++j)
      for(unsigned int i=1; i<=grid.nx-1; ++i)
      {
         state_left  = reconstruct (i, j-2, i, j-1, i, j);
         state_right = reconstruct (i, j+1, i, j, i, j-1);

        flux=num_flux_y(state_left,state_right);
         s_residual (i,j)   -= flux * grid.dx;
         s_residual (i,j-1) += flux * grid.dx;

        }
   // Boundary  periodic is set
      // inlet/outlet boundaries
   for(unsigned int n=0; n<grid.n_boundary; ++n)
   {
      int bc = grid.boundary_condition[n];

      if (grid.ibeg[n] == grid.iend[n])//&& bc != SOLID)
      {
         i = grid.ibeg[n];
         for(j=grid.jbeg[n]; j<grid.jend[n]; ++j)
         {

            if (i == 1) //  Left vertical side
            {
               state_left  = reconstruct (grid.nx-2, j, i-1, j, i, j);
               state_right = reconstruct (i+1, j, i, j, i-1, j);

            flux=num_flux_x(state_left,state_right);

               s_residual(i,j) -= flux * grid.dy;
              }
            else // right vertical side
            {
               state_left  = reconstruct (i-2, j, i-1, j, i, j);
               state_right = reconstruct (2, j, i, j, i-1, j);
               double xl=grid.xc(i-1,j);
               double xr=grid.xc(i,j);
               double y=grid.yc(i-1,j);
             flux=num_flux_x(state_left,state_right);
               s_residual(i-1,j) += flux * grid.dy;
             }
         }
      }
    //  if (grid.jbeg[n] == grid.jend[n] && bc != SOLID)
      if (grid.jbeg[n] == grid.jend[n])
      {
         j = grid.jbeg[n];
         for(i=grid.ibeg[n]; i<grid.iend[n]; ++i)
         {

            if(j == 1) // inlet-horizontal side
            {
               state_left  = reconstruct (i, grid.ny-2, i, j-1, i, j);
               state_right = reconstruct (i, j+1, i, j, i, j-1);

           flux=num_flux_y(state_left,state_right);
             s_residual(i,j)  -= flux * grid.dx;
            }
            else // outlet-horizontal side
            {
               state_left  = reconstruct (i, j-2, i, j-1, i, j);
               state_right = reconstruct (i, 2, i, j, i, j-1);


            flux=num_flux_y(state_left,state_right);
                s_residual(i,j-1)  += flux * grid.dx;
             }
         }
      }
   }

   if(rk == 0)
   {

	   dt = cfl * max (grid.dx, grid.dy)/2.0; //  f'(u)+g'(u)
      if(time<final_time && time>final_time-dt)
      dt=final_time-time;
    }

   double lambda = dt / (grid.dx * grid.dy);
   s_residual *= lambda;

     // Change the following according to your test case
    // Warning !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // For usual splitting scheme, Contribution of source terms through mid point rule
   // should be commented the following appropriately

   for( unsigned int j=1;j<grid.nx; j++)
   {
	   for( unsigned int i=1;i<grid.nx;i++)
		{
		   s_residual(i,j)-=dt*source_function(grid.xc(i,j),grid.yc(i,j),current_time);

		}
   }



}
//------------------------------------------------------------------------------
// Update solution in ghost cells
// Change this according to your test case
//------------------------------------------------------------------------------
void TwoDProblem::updateGhostCells (const double& Time)
{
   unsigned int i, j;



   if(test_case=="smooth")
   {
	   // top/bottom ghost cells
	    for (i=1; i<grid.nx; ++i)
	    {
	       j = 0;
	    solution    (i,j)= solution(i,grid.ny-1);
	       j = grid.ny;
	      solution    (i,j) = solution    (i,1);
	    }

	    // left/right ghost cells
	    for (j=1; j<grid.ny; ++j)
	    {
	       i = 0;
	     solution    (i,j) = solution    (grid.nx-1,j);
	       i = grid.nx;
	      solution    (i,j) = solution    (1,j);

	    }


   }
   else
   {
   // top/bottom ghost cells
    for (i=1; i<grid.nx; ++i)
    {
       j = 0;
   //  solution    (i,j)= solution(i,j+1);
    solution    (i,j) = exact_sol(grid.xc(i,j),0,Time);
       j = grid.ny;
    //  solution    (i,j) = solution    (i,j-1);
       solution    (i,j) = exact_sol(grid.xc(i,j),1,Time);
    }

    // left/right ghost cells
    for (j=1; j<grid.ny; ++j)
    {
       i = 0;
   //  solution    (i,j) = solution    (i+1,j);
    solution    (i,j) = exact_sol(0,grid.yc(i,j),Time);
       i = grid.nx;
  //    solution    (i,j) = solution    (i-1,j);
   solution    (i,j) = exact_sol(1,grid.yc(i,j),Time);

    }
   }
}
//------------------------------------------------------------------------------
// Find min and max values of solution
//------------------------------------------------------------------------------
void TwoDProblem::findMinMax () const
{
   double s_min = 1.0e20;
   double s_max =-1.0e20;

   for (unsigned int i=1; i<grid.nx; ++i)
      for (unsigned int j=1; j<grid.ny; ++j)
      {
         s_min = min (s_min, solution(i,j));
         s_max = max (s_max, solution(i,j));
      }

   cout << "solution    = " << s_min << " " << s_max << endl;
   cout << "dt            = " << dt << endl;
}
//------------------------------------------------------------------------------
// perform time stepping
//------------------------------------------------------------------------------
void TwoDProblem::solve (double& l1error)
{
	 current_time=0.0;
   unsigned int iter = 0;
   time = 0.0;
   dt=0.0;
   Matrix s_residual (grid.nx+1, grid.ny+1);
   Matrix s_old      (grid.nx+1, grid.ny+1);
   output_dat(0,0.0);
   vector<double>  rkt(3);

   while (iter < max_iter)
   { 

	   rkt[0]=0.0;
	   rkt[1]=dt;
	   rkt[2]=0.5*dt;
      s_old = solution;
      // Runge-Kutta stages
      for (unsigned int irk=0; irk<nrk; ++irk)
      {
          current_time=time+rkt[irk];
         n_interior_min = 0; // reset counter
         compute_maximum();
         // computing  residual
         residual (irk, s_residual);
               // update solution
         solution  = s_old * ark[irk] + (solution - s_residual) * brk[irk];
         // update solution in ghost cells
        updateGhostCells (time); //required in the case of more than 1 order
      }
      // find solution range: to check for stability
      findMinMax ();
      
      time += dt;
      ++iter;
  //    updateGhostCells (time);
     // save solution to file
     // if (iter % save_freq == 0 || iter == max_iter)
    //     output (iter,time);

      cout << "Time= " << time << " iter= " << iter << endl;
 //     cout << "No. of interior min flux = " << n_interior_min << endl;
      cout << endl;
      
     if(time==final_time)
      {// output(iter,time);
       output_dat(iter,time);
      compute_error(l1error);

      std::cout<<"L1error="<<l1error<<std::endl;
 	  break;
      }

   }
}

// To compute the maximum of absolute value of the average of solution

void TwoDProblem::compute_maximum ()
{

	abs_max_sol=0.0;

	   for(unsigned int i=1; i<=grid.nx-1; ++i)
	   {
	      for(unsigned int j=1; j<=grid.ny-1; ++j)
	      {
	    	  if( abs(solution(i,j))>abs_max_sol)
	    		  abs_max_sol=abs(solution(i,j));
          }
	   }

}

// computes the error in the solution against the exact solution at current time
void TwoDProblem::compute_error (double& L1_error)

{

L1_error=0.0;

	   for(unsigned int i=1; i<=grid.nx-1; ++i)
	   {
	      for(unsigned int j=1; j<=grid.ny-1; ++j)
	      {
	    	  L1_error+=abs(solution(i,j)-exact_sol(grid.xc(i,j),grid.yc(i,j),time));
	      }
	   }

	   L1_error*=grid.dx*grid.dy;
}


//------------------------------------------------------------------
// Savbe solution to file in the .dat format for gnuplot purpose.
//------------------------------------------------------------------
void TwoDProblem::output_dat (const unsigned int iter, const double time) const
{

   unsigned int i, j;
   ofstream dat;
   ostringstream filename;
   filename << "solution-" << iter << ".dat";

   dat.open (filename.str().c_str());


   for(unsigned int j=1; j<=grid.ny-1; ++j)
   {
      for(unsigned int i=1; i<=grid.nx-1; ++i)
      {
    	  dat<<grid.xc(i,j)<<" " << grid.yc(i,j)<<" "<<solution(i,j)<<std::endl;
       }
    	  dat<<" "<<" " <<" "<<" "<< " "<<std::endl;
    }

   dat.close();
}

//------------------------------------------------------------------------------
// save solution to file
// only interior cells are written, ghost cells are not written
//------------------------------------------------------------------------------
void TwoDProblem::output (const unsigned int iter, const double time) const
{

   unsigned int i, j;
   ofstream vtk;
   ostringstream filename;
   filename << "solution-" << iter << ".vtk";

   vtk.open (filename.str().c_str());

   vtk << "# vtk DataFile Version 2.0" << endl;
   vtk << " Two Dimensional  Problem: iter = " << iter << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET STRUCTURED_GRID" << endl;
   vtk << "FIELD FieldData 1"<<endl;
   vtk<< "TIME"<<" "<<"1"<<" "<<"1"<<" "<<"double"<<endl;
   vtk<<time<<endl;
   
   vtk << "DIMENSIONS " << grid.nx << " " << grid.ny << " 1" << endl;
   // write coordinates
   vtk << "POINTS " << grid.nx * grid.ny << " float" << endl;
   for(j=1; j<=grid.ny; ++j)
      for(i=1; i<=grid.nx; ++i)
      {
         vtk << grid.x (i,j) << "  "
             << grid.y (i,j) << "  "
             << 0.0 << endl;
      }

   vtk << "CELL_DATA " << (grid.nx-1)*(grid.ny-1) << endl;

   // write solution
   vtk << "SCALARS solution float 1" << endl;
   vtk << "LOOKUP_TABLE default" << endl;
   for(j=1; j<grid.ny; ++j)
      for(i=1; i<grid.nx; ++i)
         vtk << fixed << solution (i,j) << endl;

      vtk.close ();
}
//------------------------------------------------------------------------------
// solve the whole problem
//------------------------------------------------------------------------------
void TwoDProblem::run ( unsigned int& j, double& L1error)
{
    read_input (j);
    make_grid ();
    initialize ();
    solve (L1error);
  }
