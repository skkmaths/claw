from os import minor
import numpy as np
from scipy import optimize

def xflux(x, y, u):
    return 0.5*(u**2)

def yflux(x, y, u):
    return 0.5*(u**2)

# Upwind flux
def xnumflux(x, y, Fl, Fr, ul, ur):
    return max(xflux(x, y, max(ul,0.0)),xflux(x, y, min(ur,0.0)) )
    
# Upwind flux
def ynumflux(x, y, Fl, Fr, ul, ur):
    return max(yflux(x, y, max(ul,0.0)),yflux(x, y, min(ur,0.0)) )

# Max speed on whole grid: 
# x, y = cell centers
# u    = cell average
def max_speed(u):
    return (np.abs(u).max())


# works only for smooth solution
def uexact(x, y, t, u0):
    def imp_eqn(u):
        return u - u0(xx-t*u, yy-t*u)
    ue = np.zeros((np.size(x), np.size(y)))
    for i,xx in enumerate(x):
        for j, yy in enumerate(y):
            seed_value = u0(xx, yy)
            ue[i,j] = optimize.fsolve(imp_eqn, seed_value)
    return ue

# thing to do
# set the proper upwind flux


