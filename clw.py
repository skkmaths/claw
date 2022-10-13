"""
Solve u_t + f(u)_x = 0  for f(u) = u^2/2
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
parser.add_argument('-ic',
                    choices=('smooth','shock','rare1','rare','expo','slope'),
                    help='Initial condition', default='smooth')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

parser.add_argument('-pde', choices=('linear','varadv','burger','bucklev','concave',
                                     'burger_adv',
                                     'oreqn1','oreqn2','oreqn3','oreqn4'),
                    help='PDE', default='burger')
parser.add_argument('-numflux', help='Numerical Flux',choices=('rusanov','godunov','upwind'), default='rusanov')
parser.add_argument('-compute_error', choices=('no','yes'),
                    help='Compute error norm', default='no')
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
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

if args.ic == 'smooth':
    uinit = smooth
elif args.ic == 'shock':
    uinit = shock

xmin, xmax = -1.0, 1.0
x   = np.zeros(nc)
h = (xmax - xmin)/nc
  
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h
    
u = uinit(x)   # solution variable
res = np.zeros(nc)
ue = uexact(x, 0.0 , uinit)

# plot initial condition
if args.plot_freq >0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #line1,line2 = ax.plot(x, u, 'o',x, ue, '*')
    line1, = ax.plot(x, u, 'o')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.title('nc='+str(nc)+', CFL='+str(cfl))
    plt.grid(True); plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

#Error computation
def compute_error(u1):
    error_norm1 = 0.0; error_norm2 = 0.0
    ue = uexact(x, t , uinit)
    dom_len = xmax - xmin
    error_norm1 = h*np.sum(np.abs(u1-ue))
    error_norm2 = np.sqrt(h*np.sum((u1-ue)**2))
    return error_norm1, error_norm2

t, it = 0.0, 0
while t < Tf:
    res *=0
    dt= cfl * h #/max_speed(u)
    lam = dt/h
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/h
    # interior face
    for i in range(1,nc):
        xf = xmin+i*h # location of the face
        ul, ur = u[i-1], u[i]
        fl, fr = flux(xf,ul), flux(xf,ur)
        fn = numflux(xf, ul, ur, fl, fr)
        res[i-1] += fn
        res[i]   -= fn

    #left boundary face 
    #ul, ur = u[-1], u[0]  # periodic bc
    ul, ur = u[0], u[0]   # Neumann bc
    fl, fr = flux(xmin, ul), flux(xmin, ur)
    fn = numflux(xmin, ul, ur, fl, fr)
    res[0] -= fn

    # right right boundary face
    #ul, ur = u[-1], u[0]  # periodic bc
    ul, ur = u[-1], u[-1] # Neumann bc
    fl, fr = flux(xmax,ul), flux(xmax, ur)
    fn = numflux(xmax, ul, ur, fl, fr)
    res[-1] += fn

    u = u - lam*res  # update solution

    t += dt; it += 1 # update time step
    if args.plot_freq >0:
        ue = uexact(x, t , uinit)
        line1.set_ydata(u)
        #line2.set_ydata(ue)
        plt.draw(); plt.pause(0.1)

fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u]))
print('Saved file ', fname)

if args.compute_error == 'yes':
    er1, er2 = compute_error(u)
    print('h, L1 error norm, L2 error norm = ')
    print(h, er1, er2)
if args.plot_freq >0:
    plt.show()



