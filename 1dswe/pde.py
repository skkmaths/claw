from os import minor
import numpy as np
from scipy import optimize

# set domain
xmin, xmax = -1.0, 1.0

epsilon = 0.01
g = 9.8

def P(v, h):
    return epsilon * v * h * g

def flux_f2(v, w, h):
    return (w**2.0/v) + P(v, h) 

def numflux(xf, ul, ur, h):
    density_l = ul[0]/ul[2]
    density_r = ur[0]/ur[2]

    velocity_l = ul[1]/ul[0]
    velocity_r = ur[1]/ur[0]

    height_l, height_r = ul[2], ur[2]

    flux = np.like(ul)
    k1 = flux_f2( ul[0], max(ul[1], 0.0), ul[2])
    k2 = flux_f2( ur[0], min(ur[1], 0.0), ur[2])
    flux[0]  = max(velocity_l,0)* density_l * height_l  if k1 >= k2 else  min(velocity_r,0)* density_r * height_r 
    flux[1]  = max( k1, k2 )
    flux[3]  = max(velocity_l,0) * height_l  if k1 >= k2 else  min(velocity_r,0)*  height_r 
    
    return flux
def uexact(x, t, u0):

    ue = np.zeros(np.size(x))
    def imp_eqn(u):
        return u - u0(xx-t*u)
    for i,xx in enumerate(x):
        seed_value = u0(xx)
        ue[i] = optimize.fsolve(imp_eqn, seed_value)
    return ue

