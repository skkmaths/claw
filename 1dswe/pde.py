from os import minor
import numpy as np
from scipy import optimize

# set domain
xmin, xmax = -3.0, 3.0
epsilon = 0.01
g = 9.8

def set_epsilon(new_epsilon):
    global epsilon
    epsilon = new_epsilon  # Update global epsilon from main.py

def P(v, h):
    return epsilon * v * h * g

def flux_f2(v, w, h):
    return (w**2.0/v) + P(v, h) 
    
def B(x):
    return 1.0

def numflux(xf, ul, ur, h):
    density_l = ul[0]/ul[2]
    density_r = ur[0]/ur[2]

    velocity_l = ul[1]/ul[0]
    velocity_r = ur[1]/ur[0]

    height_l, height_r = ul[2], ur[2]

    flux = np.zeros_like(ul)

    k1 = flux_f2( ul[0], max(ul[1], 0.0), ul[2])
    k2 = flux_f2( ur[0], min(ur[1], 0.0), ur[2])
    flux[0]  = max(velocity_l,0.0)* density_l * height_l  if k1 >= k2 else  min(velocity_r,0.0)* density_r * height_r 
    flux[1]  = max( k1, k2 )
    flux[2]  = max(velocity_l,0) * height_l  if k1 >= k2 else  min(velocity_r,0)*  height_r 
    
    return flux


