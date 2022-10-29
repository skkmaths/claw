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
parser.add_argument('-scheme',
                    choices=('C','LF','GLF','LLF','LW','ROE','EROE','GOD'),
                    help='Scheme', default='LF')
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-ic',
                    choices=('smooth','shock','rare1','hat','rare','expo','slope'),
                    help='Initial condition', default='smooth')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

parser.add_argument('-pde', choices=('linear','varadv','burger','bucklev','concave',
                                     'burger_adv',
                                     'oreqn1','oreqn2','oreqn3','oreqn4'),
                    help='PDE', default='linear')
parser.add_argument('-numflux', help='Numerical Flux',choices=('rusanov','godunov','upwind'), default='rusanov')
parser.add_argument('-compute_error', choices=('no','yes'),
                    help='Compute error norm', default='no')
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
parser.add_argument('-time_scheme', default = 'euler',
                    help = 'Chosen by degree if unspecified',
                    choices = ('euler','ssprk22'))
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

# constants
Tf    = args.Tf
cfl   = args.cfl
uinit = args.ic
nc    = args.nc
alpha = args.alpha  # parameter in the slopes
beta = 2.0 # parameter in minmod
time_scheme = args.time_scheme

if args.ic == 'smooth':
    uinit = smooth
elif args.ic == 'shock':
    uinit = shock
elif args.ic == 'hat':
    uinit = hat

xmin, xmax = 0.0, 1.0
x   = np.zeros(nc)
h = (xmax - xmin)/nc
Mdx2 = args.tvbM*h**2.0
  
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h

u = np.zeros(nc+4)    # with 4 ghost cells, 2 each sides 
u[2:nc+2] = uinit(x)   # initialize solution variable
res = np.zeros(nc+4)
s_u = np.zeros(nc+4)
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
    
# plot initial condition
if args.plot_freq >0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1,line2 = ax.plot(x, u[2:nc+2], 'ro',x, u[2:nc+2], 'b')
    #line1, = ax.plot(x, u, 'o')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.title('nc='+str(nc)+', CFL='+str(cfl))
    plt.legend(('Numerical','Exact'))
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

def compute_slopes():
    s_u[:] = 0.0
    for i in range(1,nc+3):
        vl, vr = u[i-1], u[i+1]
        dvl = u[i] - vl
        dvr = vr - u[i]
        dvc = vr - vl
        s_u[i] =2.0* alpha * minmod(beta*dvl, 0.5*dvc, beta*dvr, Mdx2)

def reconstruct(ujm1, uj, ujp1):
    if args.limit == 'no':
        return uj
    elif args.limit == 'mmod':
        dj = 2.0* alpha * minmod( beta*(uj-ujm1), \
                                0.5*(ujp1-ujm1), \
                                beta*(ujp1-uj), Mdx2 )
        ul = uj + 0.5 * dj
        return ul
    else:
        print('limit type not defined')
        exit()
def update_ghost(u1):
    # left ghost cell
    u1[0] = u1[nc]
    u1[1] = u1[nc+1]

    u1[nc+3] = u1[3]
    u1[nc+2] = u1[2]
    
# First order Euler forward step
def apply_euler(t,lam, u_old, u, ures ):
    #first stage
    ts  = t
    compute_slopes(u)
    ures = compute_residual(ts, lam, u, ures)
    u = u - lam * ures
    update_ghost(u)
    return u

def apply_ssprk22(t,lam, u_old, u, ures ):
    #first stage
    ts  = t
    compute_slopes()
    ures = compute_residual(ts, lam, u, ures)
    u = u - lam * ures
    update_ghost(u)

    #second stage
    ts = t + dt
    compute_slopes()
    ures = compute_residual(ts, lam, u, ures)
    u = 0.5 * u_old + 0.5 *(u - lam * ures)
    update_ghost(u)
    return u

def compute_residual(ts, lam, u, res):
    res[:] = 0.0    
    for i in range(1,nc+2): # face between i and i+1
        xf = xmin+(i-1)*h # location of the face
        #ul, ur  = reconstruct(u[i-1], u[i], u[i+1]), reconstruct(u[i+2], u[i+1], u[i])
        ul, ur = u[i] + 0.5* s_u[i], u[i+1] - 0.5* s_u[i+1]
        fl, fr = flux(xf,ul), flux(xf,ur)
        fn = numflux(xf, ul, ur, fl, fr)
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
        plt.draw(); plt.pause(0.1)

fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u[2:nc+2]]))
print('Saved file ', fname)

if args.compute_error == 'yes':
    er1, er2 = compute_error(u[2:nc+2],t)
    print('h, L1 error norm, L2 error norm = ')
    print(h, er1, er2)
if args.plot_freq >0:
    plt.show()



