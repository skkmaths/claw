"""
Solve shallow water equations
Finite volume scheme
"""
from xml import dom
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *
from pde import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-nc', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-ic',
                    choices=('smooth','shock','dflu1','dflu2','rare1','composite','hat','buckley1','rare','expo','slope'),
                    help='Initial condition', default='smooth')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

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


# Select Numerical Flux
if args.numflux not in numfluxes:
    print("Incorrect numerical flux chosen, choices are ",numfluxes)
    exit()
if args.numflux == 'godunov':
    numflux = godunov

# constants
Tf    = args.Tf
cfl   = args.cfl
nc    = args.nc
alpha = args.alpha  # parameter in the slopes
beta = 1.0 # parameter in minmod
time_scheme = args.time_scheme

x   = np.zeros(nc)
h = (xmax - xmin)/nc
Mdx2 = args.tvbM*h**2.0
  
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h

u = np.zeros((3,nc+4))    # with 4 ghost cells, 2 each sides 
# Conserved varialbe
# u[0,:] = rho * height , u[1, :] = rho * height * velocity, u[2,:] = height
# initialize solution variable
u[0,2:nc+2] = density(x) * height (x)   
u[1,2:nc+2] = density(x) * height (x) * velocity(x)
u[2,2:nc+2] = height(x)
res = np.zeros((3,nc+4))

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
    line1,line2 = ax.plot(x, u[0,2:nc+2], 'ro',x, u[0,2:nc+2], 'b')
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
        u1[:,0] = u1[:,nc]
        u1[:,1] = u1[:,nc+1]
        u1[:,nc+3] = u1[:,3]
        u1[:,nc+2] = u1[:,2]
    elif args.bc == 'dirichlet':
        # left ghost cell
        u1[:,0] = uinit(x[0])
        u1[:,1] = uinit(x[0])
        u1[:,nc+3] = uinit(x[nc-1])
        u1[:,nc+2] = uinit(x[nc-1])
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
    res[:,:] = 0.0    
    for i in range(1,nc+2): # face between i and i+1
        xf = xmin+(i-1)*h # location of the face
        ul, ur  = u[:,i], u[:,i+1] 
        fn = numflux(xf, ul, ur, h)
        res[:,i] += fn
        res[:,i+1] -= fn
    return res
time_schemes = {'euler': apply_euler, 'ssprk22' : apply_ssprk22 }
t, it = 0.0, 0
while t < Tf:
    
    dt= cfl * h 
    lam = dt/h
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/h
    u_old = u.copy()    
    u = time_schemes[time_scheme](t, lam, u_old, u, res)  # update solution
    t += dt; it += 1 # update time step
    if args.plot_freq >0:
        ue = uexact(x, t , uinit)
        line1.set_ydata(u[0,2:nc+2])
        line2.set_ydata(ue)
        plt.title('nc='+str(nc)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
        plt.draw(); plt.pause(0.1)

# save final time solution to a file
fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u[0,2:nc+2]], u[1,2:nc+2]], u[2,2:nc+2]]))
print('Saved file ', fname)

if args.compute_error == 'yes':
    er1, er2 = compute_error(u[0,2:nc+2],t)
    print('h, L1 error norm, L2 error norm = ')
    print(h, er1, er2)
if args.plot_freq >0:
    plt.show()



