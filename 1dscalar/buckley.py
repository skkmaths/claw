from os import minor
import numpy as np
from scipy import optimize

# set domain
xmin, xmax = -1.0, 1.0

# Here M = Max( max |f'|,max|g'|)
M = 2.4

def flux(x,u):
    return 4.0*u**2/( 4.0*u**2 + (1.0-u)**2)

def dxf(x,u):
    return 8.0*u*(1.0-u)/(5.0*u**2 - 2.0*u +1)**2

# rusanov flux
def rusanov(x, ul, ur, fl, fr,lamda, h, sl, sr):
    # lam = max |f'(u)| for u between [ual, uar]
    #a =1.452608221822056
    a = np.abs([dxf(0,ul), dxf(0,ur)]).max()
    return 0.5*(fl + fr) - 0.5*a*(ur - ul)

def lxf(x, ul, ur, fl, fr,lamda, h, sl, sr ):
    return 0.5*(fl + fr) - 0.5*(ur - ul)/lamda

def godunov(x, ul, ur, fl, fr, lamda,h, sl, sr):
    return flux(0,ul)

# Max speed based on cell average values
def max_speed(u):
    return np.abs(u).max()
# numflux for Nessyahu Tadmore scheme
def nt(x, ul, ur, fl, fr, lamda,h, sl, sr ):
    #ui   = ul-sl/2.0
    #uip1 = ur + sr/2.0
    Fl = flux(x, ul-lamda *dxf(x, ul)* sl/4.0 ) + sl/(lamda*4.0)
    Fr = flux(x, ur-lamda * dxf(x, ur)* sr/4.0) + sr/(lamda * 4.0)
    return  0.5*(Fl+Fr)-0.5*( ur-ul )/lamda
# works only for smooth solution

# exact solution of rare1 initial condition


def uexact(x,t,u0):
    if x < 0.5*t:
        return 1.0
    else:
        return 0.0
uexact = np.vectorize(uexact)

'''
def uexact(x,t,u0):
    if x < -t:
        return -1.0
    elif x > t:
        return 1.0
    else:
        return x/t
uexact = np.vectorize(uexact)
def uexact(x, t, u0):

    ue = np.zeros(np.size(x))
    def imp_eqn(u):
        return u - u0(xx-t*u)
    for i,xx in enumerate(x):
        seed_value = u0(xx)
        ue[i] = optimize.fsolve(imp_eqn, seed_value)
    return ue
'''
numfluxes = ['rusanov','nt', 'lxf', 'godunov', 'muscl']