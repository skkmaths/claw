from os import minor
import numpy as np
from scipy import optimize

# f = u^2/2
def flux(x,u):
    return 0.5*u**2

# rusanov flux
def rusanov(x, ul, ur, fl, fr):
    # lam = max |f'(u)| for u between [ual, uar]
    a = np.abs([ul, ur]).max()
    return 0.5*(fl + fr) - 0.5*a*(ur - ul)

def godunov(x, ul, ur, fl, fr):
    return max(flux(x,max(ul,0.0),0.0),flux(x,min(ur,0.0),0.0) )

# Max speed based on cell average values
def max_speed(u):
    return np.abs(u).max()
# numflux for Nessyahu Tadmore scheme
def nt(x,ul, ur, fl, fr, dul, dur, lam):
    return 0.5*(fl+fr)-(0.5/lam)*(ur-ul) + (0.25/lam)*(dul+dur)


def uexact(x, t, u0):
    ue = np.zeros(np.size(x))
    def imp_eqn(u):
        return u - u0(xx-t*u)
    for i,xx in enumerate(x):
        seed_value = u0(xx)
        ue[i] = optimize.fsolve(imp_eqn, seed_value)
    return ue

numfluxes = ['rusanov','nt', 'godunov']