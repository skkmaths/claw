from os import minor
import numpy as np
from scipy import optimize
import sympy as sp

# discontinuous flux conservation law
# u_t + F(x, u)_x = 0
# set domain
xmin, xmax = -4, 4.0
u = sp.symbols('u')
g =  2.0*u*(1.0-u)/(1.0+u)
f =  2.0*u*(1.0-u)/(2.0-u)

dxg = sp.diff(g,u)
dxf = sp.diff(f,u)

g = sp.lambdify(u,g)
f = sp.lambdify(u,f)

dxg = sp.lambdify(u,dxg)
dxf = sp.lambdify(u,dxf)

def H(x):
    if x <= 0:
        return 0.0
    else:
        return 1.0

def flux(x,u):
    return  (1.0 - H(x)) * g(u) + H(x) * f(u) 
def dxflux(x,u):
    return  (1.0 - H(x)) * dxg(u) + H(x) * dxf(u) 
# Here M = Max( max |f'|,max|g'|)
M = 2.0

theta_g = -1 + np.sqrt(2)
theta_f = 2.0 - np.sqrt(2)

def dflu(x, ul, ur, fl, fr, lamda, h, sl, sr):
    if x < 0:
        return min(flux(x,min(ul,theta_g)),flux(x,max(ur,theta_g)) )
    elif x > 0.0:
        return min(flux(x,min(ul,theta_f)),flux(x,max(ur,theta_f)) )
    else:
        return min(g(min(ul,theta_g)),f(max(ur,theta_f)) )
# Lax Friedrich of Karlsen 2004
def lxf(x, ul, ur, fl, fr, lamda, h, sl, sr):
    # center of left and right cells
    xcl = x -0.5*h
    xcr = x+0.5*h
    return 0.5*( flux(xcl, ul) + flux(xcr,ur) - (ur -ul)/lamda  )

# numflux for Nessyahu Tadmore scheme
def nt(x, ul, ur, fl, fr, lamda, h, sl, sr ):
    ui   = ul-sl/2.0
    uip1 = ur + sr/2.0
    xcl = x -0.5*h
    xcr = x + 0.5*h
    Fl = flux(xcl, ui-0.5 * lamda * dxflux(xcl, ui)* sl ) + 0.5 * sl/lamda
    Fr = flux(xcr, uip1-0.5 * lamda * dxflux(xcr, uip1)* sr) + 0.5 * sr/lamda
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

numfluxes = ['dflu', 'nt','lxf']