import numpy as np
from constants import *

xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0

lvbc, rvbc, bhbc, thbc = periodic, periodic, periodic, periodic

def initial_condition(x,y):
    return np.sin(2*np.pi*x) * np.sin(2*np.pi*y)
