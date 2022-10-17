import numpy as np

# f = u
def flux(x,u):
    return u

# rusanov flux
def rusanov(x, ul, ur, fl, fr):
    # lam = max |f'(u)| for u between [ual, uar]
    return 0.5*(fl + fr) - 0.5*(ur - ul)

def godunov(x, ul, ur, fl, fr):
    return ul
# Max speed based on cell average values
def max_speed(u):
    return 1.0

# works for any initial data    
def uexact(x, t, u0):
    return  u0(x-t)

numfluxes = ['rusanov','godunov', 'nt']