"""
Solve u_t + f(u)_x = 0  for f(u) 
Finite volume scheme
"""
from xml import dom
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-nc', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-ic',
                    choices=('smooth','shock','dflu1','dflu2','rare1','composite','hat','buckley1','rare','expo','slope'),
                    help='Initial condition', default='smooth')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

parser.add_argument('-pde', choices=('linear','dflux','varadv','burger','buckley','concave',
                                     'burger_adv',
                                     'oreqn1','oreqn2','oreqn3','oreqn4'),
                    help='PDE', default='linear')
parser.add_argument('-numflux', help='Numerical Flux',choices=('rusanov','dflu','godunov','upwind','nt', 'lxf'), default='rusanov')
parser.add_argument('-compute_error', choices=('no','yes'),
                    help='Compute error norm', default='no')
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
parser.add_argument('-time_scheme', default = 'euler',
                    help = 'Chosen by degree if unspecified',
                    choices = ('euler','ssprk22'))
parser.add_argument('-bc', default = 'periodic',
                    help = 'Chose the boundary condition',
                    choices = ('periodic','dirichlet'))
parser.add_argument('-limit', choices=('no', 'mmod'), help='Apply limiter',
                    default='no')
parser.add_argument('-alpha', type=float, help='Slope parameter', default=0.5)
args = parser.parse_args()
# Select PDE

if args.pde == 'burger':
    from burger import *
elif args.pde == 'linear':
    from linear import *
elif args.pde == 'concave':
    from concave import *
elif args.pde == 'dflux':
    from dflux import *
elif args.pde == 'buckley':
    from buckley import *

# Select Numerical Flux
if args.numflux not in numfluxes:
    print("Incorrect numerical flux chosen, choices are ",numfluxes)
    exit()
if args.numflux == 'rusanov':
    numflux = rusanov
elif args.numflux == 'godunov':
    numflux = godunov
elif args.numflux == 'upwind':
    numflux = upwind
elif args.numflux == 'nt':
    numflux = nt
elif args.numflux == 'lxf':
    numflux = lxf
elif args.numflux == 'dflu':
    numflux = dflu
elif args.numflux == 'buckley1':
    numflux = buckley1
# constants
Tf    = args.Tf
cfl   = args.cfl
uinit = args.ic
nc    = args.nc
alpha = args.alpha  # parameter in the slopes
beta = 1.0 # parameter in minmod
time_scheme = args.time_scheme

if args.ic == 'smooth':
    uinit = smooth
elif args.ic == 'shock':
    uinit = shock
elif args.ic == 'hat':
    uinit = hat
elif args.ic == 'dflu1':
    uinit = dflu1
elif args.ic == 'dflu2':
    uinit = dflu2
elif args.ic == 'rare1':
    uinit = rare1
elif args.ic == 'composite':
    uinit = composite
elif args.ic == 'buckley1':
    uinit = buckley1
x   = np.zeros(nc)
h = (xmax - xmin)/nc
Mdx2 = args.tvbM*h**2.0
  
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h

u = np.zeros(nc+4)    # with 4 ghost cells, 2 each sides 
u[2:nc+2] = uinit(x)   # initialize solution variable
res = np.zeros(nc+4)
s_u = np.zeros(nc+4)
'''
def minmod(a,b,c,Mdx2):
    if np.abs(a) < Mdx2:
        return a
    sa = np.sign(a)
    sb = np.sign(b)
    sc = np.sign(c)
    if sa==sb and sb==sc:
        return sa * np.abs([a,b,c]).min()
    else:
        return 0.0
'''
def minmod(a,b,c):
    sa = np.sign(a)
    sb = np.sign(b)
    sc = np.sign(c)
    if sa==sb and sb==sc:
        return sa * np.abs([a,b,c]).min()
    else:
        return 0.0
t = 0.0 # time
# plot initial condition
if args.plot_freq >0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1,line2 = ax.plot(x, u[2:nc+2], 'ro',x, u[2:nc+2], 'b')
    #line1, = ax.plot(x, u, 'o')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.title('nc='+str(nc)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
    plt.legend(('Numerical','Exact'))
    #plt.ylim(0.3,0.7)
    plt.grid(True); plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

#Error computation
def compute_error(u1,t):
    error_norm1 = 0.0; error_norm2 = 0.0
    ue = uexact(x, t , uinit)
    dom_len = xmax - xmin
    error_norm1 = h*np.sum(np.abs(u1-ue))
    error_norm2 = np.sqrt(h*np.sum((u1-ue)**2))
    return error_norm1, error_norm2

def reconstruct(ujm1, uj, ujp1):
    if args.limit == 'no':
        return uj
    elif args.limit == 'mmod':
        dj = 2.0* alpha * minmod( beta*(uj-ujm1), \
                                0.5*(ujp1-ujm1), \
                                beta*(ujp1-uj))
        ul = uj + 0.5 * dj
        return ul
    else:
        print('limit type not defined')
        exit()
def update_ghost(u1):
    if args.bc == 'periodic':
        # left ghost cell
        u1[0] = u1[nc]
        u1[1] = u1[nc+1]
        u1[nc+3] = u1[3]
        u1[nc+2] = u1[2]
    elif args.bc == 'dirichlet':
        # left ghost cell
        u1[0] = uinit(x[0])
        u1[1] = uinit(x[0])
        u1[nc+3] = uinit(x[nc-1])
        u1[nc+2] = uinit(x[nc-1])
    else:
        print('unknown boundary condition')
        exit()
    
# First order Euler forward step
def apply_euler(t,lam, u_old, u, ures ):
    #first stage
    ts  = t
    update_ghost(u)
    ures = compute_residual(ts, lam, u, ures)
    u = u - lam * ures
    return u

def apply_ssprk22(t,lam, u_old, u, ures ):
    #first stage
    ts  = t
    update_ghost(u)
    ures = compute_residual(ts, lam, u, ures)
    u = u - lam * ures
    #second stage
    ts = t + dt
    update_ghost(u)
    ures = compute_residual(ts, lam, u, ures)
    u = 0.5 * u_old + 0.5 *(u - lam * ures)
    return u

def compute_residual(ts, lam, u, res):
    res[:] = 0.0    
    for i in range(1,nc+2): # face between i and i+1
        xf = xmin+(i-1)*h # location of the face
        ul, ur  = reconstruct(u[i-1], u[i], u[i+1]), reconstruct(u[i+2], u[i+1], u[i])
        #compute slopes for NT scheme written in non-stagered form
        if args.numflux == 'nt':
            uim1 = u[i-1]
            ui   = u[i]
            uip1 = u[i+1]
            uip2 = u[i+2]
            if ( i ==1):
                uim2 = u[-5]
            else:
                uim2 = u[i-2]
            if i== nc+1:
                uip3 = u[4]
            else: 
                uip3 = u[i+3]

            sl = 2.0* alpha * minmod( beta*(uip2-ui), \
                                0.5*(uip2-uim2), \
                                beta*(ui-uim2))
            sr = 2.0* alpha * minmod( beta*(uip3-uip1), \
                                0.5*(uip3-uim1), \
                                beta*(uip1-uim1))
        else: # dummy variables
            sl = 2.0 * (ul-u[i])
            sr = 2.0 * (u[i+1]-ur)
        fl, fr = flux(xf,ul), flux(xf,ur)
        if args.numflux == 'nt':
            ul = u[i]
            ur = u[i+1]
        fn = numflux(xf, ul, ur, fl, fr, lam, h, sl, sr)
        res[i] += fn
        res[i+1] -= fn
    return res
time_schemes = {'euler': apply_euler, 'ssprk22' : apply_ssprk22 }
t, it = 0.0, 0
while t < Tf:
    if args.pde == 'linear':
        dt= cfl * h
    elif args.pde == 'burger':
        dt= cfl * h /max_speed(u)
    elif args.pde == 'dflux':
        dt= cfl * h /M
    elif args.pde == 'buckley':
        dt= cfl * h /M
    else:
        print('dt is not set')
        exit()
    lam = dt/h
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/h
    u_old = u    
    u = time_schemes[time_scheme](t, lam, u_old, u, res)  # update solution
    t += dt; it += 1 # update time step
    if args.plot_freq >0:
        ue = uexact(x, t , uinit)
        line1.set_ydata(u[2:nc+2])
        line2.set_ydata(ue)
        plt.title('nc='+str(nc)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
        plt.draw(); plt.pause(0.1)

# save exact solution at final time
xe = np.linspace(xmin, xmax, 10000)
ue = uexact(xe,0, uinit)
fname = 'exact.txt'
np.savetxt(fname, np.column_stack([xe, ue]))
print('Saved file ', fname)

# save final time solution to a file
fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u[2:nc+2]]))
print('Saved file ', fname)

if args.compute_error == 'yes':
    er1, er2 = compute_error(u[2:nc+2],t)
    print('h, L1 error norm, L2 error norm = ')
    print(h, er1, er2)
if args.plot_freq >0:
    plt.show()



