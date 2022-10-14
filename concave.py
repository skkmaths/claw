from os import minor
import numpy as np
from scipy import optimize

# f = u(1-u)
def flux(x,u):
    return u*(1-u)

# rusanov flux
def rusanov(x, ul, ur, fl, fr): # need to fix this
    # lam = max |f'(u)| for u between [ual, uar]
    a = np.abs([ul, ur]).max()
    return 0.5*(fl + fr) - 0.5*a*(ur - ul)

def upwind(x, ul, ur, fl , fr ):
    return ul*(1.0-ur)

def godunov(x, ul, ur, fl, fr):
    return min(flux(x,min(ul,0.5)),flux(x,max(ur,0.5)) )

# Max speed based on cell average values
def max_speed(u):
    return np.abs(u).max()
# numflux for Nessyahu Tadmore scheme
def nt(x,ul, ur, fl, fr, dul, dur, lam):
    return 0.5*(fl+fr)-(0.5/lam)*(ur-ul) + (0.25/lam)*(dul+dur)


def uexact(x, t, u0): 
    ue = np.zeros(np.size(x))
    return ue

numfluxes = ['rusanov','nt', 'godunov', 'upwind']