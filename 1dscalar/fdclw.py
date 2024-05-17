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
                    choices=('smooth','shock','dflu1','dflu2','rare1','hat','rare','expo','slope'),
                    help='Initial condition', default='smooth')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

parser.add_argument('-pde', choices=('linear','dflux','varadv','burger','bucklev','concave',
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
                    choices = ('periodic','dc'))
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

x   = np.zeros(nc)
h = (xmax - xmin)/nc
Mdx2 = args.tvbM*h**2.0
  
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h

u = np.zeros(nc+4)    # with 4 ghost cells, 2 each sides 
u[2:nc+2] = uinit(x)   # initialize solution variable
res = np.zeros(nc+4)
slope = np.zeros(nc+4)
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
    elif args.bc == 'dc':
        # left ghost cell
        u1[0] = uinit(x[0])
        u1[1] = uinit(x[0])
        u1[nc+3] = uinit(x[nc-1])
        u1[nc+2] = uinit(x[nc-1])
    else:
        print('unknown boundary condition')
        exit()
    return u1
def compute_slopes(slope1, u1):
    for i in range(2,nc+2):
        slope1[i] = 2.0* alpha * minmod( beta*(u1[i]-u1[i-1]), \
                                0.5*(u1[i+1]-u1[i-1]), \
                                beta*(u1[i+1]-u1[i]))
    return slope1

def update_lf(lam, u ):
    if args.pde == 'burger':
        f = 0.5*u*u
    elif args.pde == 'linear':
        f = u
    unew = np.empty_like(u)
    for i in range(2, len(u)-2):
        unew[i] = 0.5*(u[i-1]+u[i+1]) - 0.5*lam*(f[i+1] - f[i-1])
    return unew
def update_nt(lam, slope, u):
    slope = compute_slopes(slope,u)
    unew = np.empty_like(u)
    for i in range(2, len(u)-2):
        unew[i] = (u[i+1]+u[i-1])/2.0 +  (slope[i-1]-slope[i+1])/4.0 \
                -( flux(x[1], u[i+1] -dxf(x[1], u[i+1])*slope[i+1]*lam/2.0) - \
                       flux(x[1], u[i-1] -dxf(x[1], u[i-1])*slope[i-1]*lam/2.0  ))*lam/2.0
    return unew
    
t, it = 0.0, 0
while t < Tf:
    if args.pde == 'linear':
        dt= cfl * h
    elif args.pde == 'burger':
        dt= cfl * h /max_speed(u)
    elif args.pde == 'dflux':
        dt= cfl * h /M
    else:
        print('dt is not set')
        exit()
    lam = dt/h
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/h
    u = update_ghost(u)
    if args.numflux == 'lxf':
        u = update_lf(lam, u)
    elif args.numflux == 'nt':
        u = update_nt(lam, slope, u)
    else:
        print('unknown flux')
        exit()
    t += dt; it += 1 # update time step
    if args.plot_freq >0:
        ue = uexact(x, t , uinit)
        line1.set_ydata(u[2:nc+2])
        line2.set_ydata(ue)
        plt.title('nc='+str(nc)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
        plt.draw(); plt.pause(0.1)
ue = uexact(x,t, uinit)
fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u[2:nc+2], ue]))
print('Saved file ', fname)

if args.compute_error == 'yes':
    er1, er2 = compute_error(u[2:nc+2],t)
    print('h, L1 error norm, L2 error norm = ')
    print(h, er1, er2)
if args.plot_freq >0:
    plt.show()



