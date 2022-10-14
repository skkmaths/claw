import numpy as np

# (f, g) = (-y*u, x*u)
def xflux(x, y, u):
    return -y*u

def yflux(x, y, u):
    return x*u

def advection_velocity(x,y):
    return (-y,x)

# Upwind flux
def xnumflux(x, y, Fl, Fr, ul, ur):
    sx, sy = advection_velocity(x,y)
    return 0.5*(Fl + Fr) - 0.5*np.abs(sx)*(ur - ul)

# Upwind flux
def ynumflux(x, y, Fl, Fr, ul, ur):
    sx, sy = advection_velocity(x,y)
    return 0.5*(Fl + Fr) - 0.5*np.abs(sy)*(ur - ul)

# Return advection velocity
def local_speed(x, y, u):
    return advection_velocity(x,y)

# Max speed on whole grid: 
# x, y = cell centers
# u    = cell average
def max_speed(x, y, u):
    sx, sy = advection_velocity(x,y)
    return (np.abs(sx).max(), np.abs(sy).max())
