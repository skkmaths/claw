from os import minor
import numpy as np
from scipy import optimize
import sympy as sp
from scipy.optimize import minimize_scalar
#-----------------set dflux --------------------------------------------------------------
# discontinuous flux conservation law
# u_t + F(x, u)_x = 0
# set domain
xmin, xmax = -4, 4.0
u = sp.symbols('u')
#g =  2.0*u*(1.0-u)/(1.0+u)
#f =  2.0*u*(1.0-u)/(2.0-u)
#theta_g = -1 + np.sqrt(2)
#theta_f = 2.0 - np.sqrt(2)
g = 20.0 * u * u * (1.0-u)**2 / ( u**2 + 2.0 * ( 1.0-u)**2)
f = 50.0 * u * u * (1.0-u)**2 / ( 10.0 * u**2 + ( 1.0-u)**2)
theta_f = 0.557506665975558
theta_g = 0.317014
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

#-----------------set dflux --------------------------------------------------------------
# Here M = Max( max |f'|,max|g'|)
#------------------------compute M ------------------------------------------------------
# maximize |g'(u)|
res_g = minimize_scalar(
    lambda x: -abs(dxg(x)),
    bounds=(0.0, 1.0),
    method='bounded'
)
# maximize |f'(u)|
res_f = minimize_scalar(
    lambda x: -abs(dxf(x)),
    bounds=(0.0, 1.0),
    method='bounded'
)
Mg = abs(dxg(res_g.x))
Mf = abs(dxf(res_f.x))
M = max(Mg, Mf)
#------------------------compute M ------------------------------------------------------
# Variables for computing the diffusion coefficients
y, z = sp.symbols('y z', real=True)
A, B, kappa = sp.symbols('A B kappa', positive=True)

# Left branch
aL = sp.Piecewise(
    (kappa*z/A, z <= A),
    (kappa + (z-A)/(1-A)*(1-kappa), True)
)

# Right branch
aR = sp.Piecewise(
    (kappa*z/B, z <= B),
    (kappa + (z-B)/(1-B)*(1-kappa), True)
)

# Complete function
a = sp.Piecewise(
    (aL, y < 0.0),
    (aR, True)
)
# Define the vaues of A and B and Kappa here
a_num = a.subs({
    A: 0.63839972,
    B: 0.317014,
    #A: 0.317014,
    #B: 0.472372,
    kappa: 0.5
})
# define the function as a variable in y and z
diffusion = sp.lambdify((y, z), a_num, "numpy")



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

def llfd(x, ul, ur, fl, fr, lamda, h, sl, sr):
    # center of left and right cells
    xcl = x -0.5*h
    xcr = x+0.5*h
    # find local speed
    speed = max ( abs( dxflux(xcl, ul) ), abs( dxflux(xcr,ur) ) )
    if x < 0.0 or x >0.0:
        return 0.5*( flux(x, ul) + flux(x, ur) - speed * ( ur  - ul )  )
    else:
        return 0.5*( flux(xcl, ul) + flux(xcr,ur) - speed * ( diffusion(xcr,ur)  - diffusion(xcl,ul) )  )

def llf(x, ul, ur, fl, fr, lamda, h, sl, sr):
    # center of left and right cells
    xcl = x -0.5*h
    xcr = x+0.5*h
    # find local speed
    speed = max ( abs( dxflux(xcl, ul) ), abs( dxflux(xcr,ur) ) )
    return 0.5*( flux(xcl, ul) + flux(xcr,ur) - speed * ( diffusion(xcr,ur)  - diffusion(xcl,ul) )  )

# numflux for Nessyahu Tadmore scheme
def nt(x, ul, ur, fl, fr, lamda, h, sl, sr ):
    xcl = x -0.5*h
    xcr = x + 0.5*h
    Fl = flux(xcl, ul-lamda *dxflux(xcl, ul)* sl/4.0 ) + sl/(lamda*4.0)
    Fr = flux(xcr, ur-lamda * dxflux(xcr, ur)* sr/4.0) + sr/(lamda * 4.0)
    return  0.5*(Fl+Fr)-0.5*( ur-ul )/lamda

# works only for smooth solution
def uexact(x, t, u0):
    ue = np.zeros(np.size(x))
    def imp_eqn(u):
        return u - u0(xx-t*u)
    for i,xx in enumerate(x):
        seed_value = u0(xx)
        ue[i] = optimize.fsolve(imp_eqn, seed_value)
    return ue

numfluxes = ['dflu', 'nt','lxf', 'llf', 'llfd']