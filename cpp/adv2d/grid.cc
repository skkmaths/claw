#include <iostream>
#include <cassert>
#include <vector>
#include "grid.h"

using namespace std;

void Grid::allocate ()
{
   assert (nx > 1);
   assert (ny > 1);
   // only real cell vertices
   x.allocate (nx+1,ny+1);
   y.allocate (nx+1,ny+1);
   
   // only real cell centers
   xc.allocate (nx,ny);
   yc.allocate (nx,ny);

}

