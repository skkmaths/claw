from os import minor
import numpy as np
from scipy import optimize

# set domain
xmin, xmax = -1.0, 1.0

# f = u^2/2
def flux(x,u):
    return 0.5*u**2
def dxf(x,u):
    return u

# rusanov flux
def rusanov(x, ul, ur, fl, fr,lamda, h, sl, sr):
    # lam = max |f'(u)| for u between [ual, uar]
    a = np.abs([ul, ur]).max()
    return 0.5*(fl + fr) - 0.5*a*(ur - ul)

def lxf(x, ul, ur, fl, fr,lamda, h, sl, sr ):
    return 0.5*(fl + fr) - 0.5*(ur - ul)/lamda

def godunov(x, ul, ur, fl, fr, lamda,h, sl, sr):
    return max(flux(x,max(ul,0.0)),flux(x,min(ur,0.0)) )

# Max speed based on cell average values
def max_speed(u):
    return np.abs(u).max()
# numflux for Nessyahu Tadmore scheme
def nt(x, ul, ur, fl, fr, lamda,h, sl, sr ):
    ui   = ul-sl/2.0
    uip1 = ur + sr/2.0
    Fl = flux(x, ui-0.5 * lamda *dxf(x, ui)* sl ) + 0.5 * sl/lamda
    Fr = flux(x, uip1-0.5 * lamda * dxf(x, uip1)* sr) + 0.5 * sr/lamda
    return  0.5*(Fl+Fr)-0.5*( uip1-ui )/lamda

# works only for smooth solution
def uexact(x, t, u0):
    ue = np.zeros(np.size(x))
    def imp_eqn(u):
        return u - u0(xx-t*u)
    for i,xx in enumerate(x):
        seed_value = u0(xx)
        ue[i] = optimize.fsolve(imp_eqn, seed_value)
    return ue

numfluxes = ['rusanov','nt', 'lxf','godunov']