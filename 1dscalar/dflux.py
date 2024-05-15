from os import minor
import numpy as np
from scipy import optimize

# discontinuous flux conservation law
# u_t + F(x, u)_x = 0

def g(u):
    return 2.0*u*(1.0-u)/(1.0+u)
def f(u):
    return  2.0*u*(1.0-u)/(2.0-u)
def H(x):
    if x <= 0:
        return 0.0
    else:
        return 1.0

def flux(x,u):
    return  (1.0 - H(x)) * g(u) + H(x) * f(u) 
# Here M = Max( max |f'|,max|g'|)
M = 2.0

theta_g = -1 + np.sqrt(2)
theta_f = 2.0 - np.sqrt(2)

def dflu(x, ul, ur, fl, fr, lamda, h):
    if x < 0:
        return min(flux(x,min(ul,theta_g)),flux(x,max(ur,theta_g)) )
    elif x > 0.0:
        return min(flux(x,min(ul,theta_f)),flux(x,max(ur,theta_f)) )
    else:
        return min(g(min(ul,theta_g)),f(max(ur,theta_f)) )
# Lax Friedrich of Karlsen 2004
def lxf(x, ul, ur, fl, fr, lamda, h):
    # center of left and right cells
    xcl = x -0.5*h
    xcr = x+0.5*h
    return 0.5*( flux(xcl, ul) + flux(xcr,ur) - (ur -ul)/lamda  )

# numflux for Nessyahu Tadmore scheme


# works only for smooth solution
def uexact(x, t, u0):
    ue = np.zeros(np.size(x))
    def imp_eqn(u):
        return u - u0(xx-t*u)
    for i,xx in enumerate(x):
        seed_value = u0(xx)
        ue[i] = optimize.fsolve(imp_eqn, seed_value)
    return ue

numfluxes = ['dflu', 'lxf']