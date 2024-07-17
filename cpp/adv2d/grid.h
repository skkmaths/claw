#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include "matrix.h"

enum BType { imin, imax, jmin, jmax };

// class to hold the grid object
class Grid
{
   public:

      Grid () { nx = ny = 0;};
      ~Grid () {};

      unsigned int nx, ny;
      Matrix x, y;   // grid vertices
      Matrix xc, yc; // cell centers
      double xmin, xmax, ymin, ymax;
      double dx, dy;
      void allocate ();

};

#endif
