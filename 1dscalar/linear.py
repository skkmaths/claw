import numpy as np
#set domain
xmin, xmax = -1.0, 1.0
# f = u
def flux(x,u):
    return u
def dxf(x,u):
    return 1.0

# rusanov flux
def rusanov(x, ul, ur, fl, fr,lamda,h, sl, sr):
    # lam = max |f'(u)| for u between [ual, uar]
    return 0.5*(fl + fr) - 0.5*(ur - ul)
def lxf(x, ul, ur, fl, fr,lamda,h, sl, sr ):
    return 0.5*(fl + fr) - 0.5*(ur - ul)/lamda
def godunov(x, ul, ur, fl, fr, lamda, h,  sl, sr):
    return ul
# Max speed based on cell average values
def max_speed(u):
    return 1.0

# works for any initial data    
def uexact(x, t, u0):
    return  u0(x-t)

#NT numerical flux 
def nt(x, ul, ur, fl, fr, lamda,h, sl, sr ):
    ui   = ul-sl/2.0
    uip1 = ur + sr/2.0
    Fl = flux(x, ui-0.5 * lamda *dxf(x, ui)* sl ) + 0.5 * sl/lamda
    Fr = flux(x, uip1-0.5 * lamda * dxf(x, uip1)* sr) + 0.5 * sr/lamda
    return  0.5*(Fl+Fr)-0.5*( uip1-ui )/lamda

numfluxes = ['rusanov','godunov', 'nt', 'lxf']