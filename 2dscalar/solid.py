import numpy as np
from constants import *

xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0

lvbc, rvbc, bhbc, thbc = periodic, periodic, periodic, periodic

def initial_condition(x,y):
    r = np.sqrt((x+0.45)**2 + y**2)
    if (((x > 0.1) and (x < 0.6)) and ((y > -0.25) and (y < 0.25))):
        return 1.0
    if r < 0.35:
        return  1.0- (r/0.35)
    else:
        return 0.0     