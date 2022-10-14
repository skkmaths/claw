import numpy as np
from constants import *

xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0

lvbc, rvbc, bhbc, thbc = periodic, periodic, periodic, periodic

def initial_condition(x,y):
    return 1.0 + np.exp(-100.0*((x-0.5)**2 + y**2))
