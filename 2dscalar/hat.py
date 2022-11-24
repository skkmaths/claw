import numpy as np
from constants import *

xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0

#lvbc, rvbc, bhbc, thbc = periodic, periodic, periodic, periodic
#lvbc, rvbc, bhbc, thbc = absorbing, absorbing, absorbing, absorbing

def initial_condition(x,y):
   if (x >0.2 and x <0.5):
       if (y >-0.15 and y < 0.15):   
           return 1.0
       else:
        return 0.0
   else:
        return 0.0     