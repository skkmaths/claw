import numpy as np

# (f, g) = (a*u, a*u) with a = 1
def xflux(x, y, u):
    return u

def yflux(x, y, u):
    return u

def advection_velocity(x,y):
    return (1,1)

# Upwind flux
def xnumflux(x, y, Fl, Fr, ul, ur):
    return 0.5*(Fl + Fr) - 0.5*(ur - ul)

def ynumflux(x, y, Fl, Fr, ul, ur):
    return 0.5*(Fl + Fr) - 0.5*(ur - ul)

# Return (1,1) at all points (x,y)
def local_speed(x, y, u):
    return (np.ones(x.shape), np.ones(x.shape))

def max_speed(x, y, u):
    return (1.0, 1.0)
