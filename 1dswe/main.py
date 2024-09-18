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
parser.add_argument('-cfl', type=float, help='CFL number', default=0.5)
parser.add_argument('-eps', type=float, help='Epsilon', default=0.01)

parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
parser.add_argument('-time_scheme', default = 'euler',
                    help = 'Chosen by degree if unspecified',
                    choices = ('euler','ssprk22'))
parser.add_argument('-bc', default = 'dirichlet',
                    help = 'Chose the boundary condition',
                    choices = ('periodic','dirichlet'))
args = parser.parse_args()


# constants
Tf    = args.Tf
cfl   = args.cfl
nc    = args.nc
time_scheme = args.time_scheme
set_epsilon(args.eps)
x   = np.zeros(nc)
h = (xmax - xmin)/nc

# create real cells, 
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h

u = np.zeros((3,nc+4))    # with total 4 ghost cells, 2 each sides 
# Conserved varialbe
# u[0,:] = rho * height , u[1, :] = rho * height * velocity, u[2,:] = height
# initialize solution variable
u[0,2:nc+2] = density(x) * height (x)   
u[1,2:nc+2] = density(x) * height (x) * velocity(x)
u[2,2:nc+2] = height(x)
res = np.zeros((3,nc+4))

t = 0.0 # time
# plot initial condition
if args.plot_freq > 0:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # Create 1 row, 3 columns of subplots
    
    # Plot on the first subplot (axes[0])
    line1, = axes[0].plot(x, u[0, 2:nc+2]/u[2, 2:nc+2], '-')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('density')
    #axes[0].set_title('Plot 1: nc='+str(nc)+', CFL='+str(cfl)+', time='+str(np.round(t, 3)))
    axes[0].legend(('Numerical'))
    axes[0].grid(True)
    
    # Plot on the second subplot (axes[1])
    line2, = axes[1].plot(x, u[1, 2:nc+2]/u[0, 2:nc+2], '-')  # Example plot for second data set
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('velocity')
    #axes[1].set_title('Plot 2: nc='+str(nc)+', CFL='+str(cfl)+', time='+str(np.round(t, 3)))
    axes[1].legend(('Numerical'))
    axes[1].grid(True)

    # Plot on the third subplot (axes[2])
    line3, = axes[2].plot(x, u[2, 2:nc+2], '-')  # Example plot for third data set
    axes[2].set_xlabel('x')
    axes[2].set_ylabel('height')
    #axes[2].set_title('Plot 3: nc='+str(nc)+', CFL='+str(cfl)+', time='+str(np.round(t, 3) ) )
    axes[2].legend('Numerical')
    axes[2].grid(True)
    axes[2].set_ylim([np.min(u[2,2:nc+2])-0.1, np.max(u[2,2:nc+2])+0.1])
    fig.suptitle('nc='+str(nc)+', CFL='+str(cfl)+', time='+str(np.round(t, 3)))
    plt.tight_layout()  # Adjust subplots to fit in the figure area
    plt.pause(0.1)
    wait = input("Press enter to continue ")

def update_ghost(u1):
    if args.bc == 'periodic':
        # left ghost cell
        u1[:,0] = u1[:,nc]
        u1[:,1] = u1[:,nc+1]
        u1[:,nc+3] = u1[:,3]
        u1[:,nc+2] = u1[:,2]
    elif args.bc == 'dirichlet':
        # left ghost cell
        u1[:,0] = u1[:,2]
        u1[:,1] = u1[:,2]
        u1[:,nc+3] = u1[:,nc+1]
        u1[:,nc+2] = u1[:,nc+1]
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

def compute_residual(ts, lam, u, res):
    res[:,:] = 0.0    
    # loop over the faces and add to res
    for i in range(1,nc+2): # face between i and i+1
        xf = xmin+(i-1)*h # location of the face
        ul, ur  = u[:,i], u[:,i+1] 
        fn = numflux(xf, ul, ur, h)
        res[:,i] += fn
        res[:,i+1] -= fn
    # add source terms to res
    for i in range (2, nc+2):
        k1 = flux_f2( u[0,i], max(u[1,i], 0.0), u[2,i])
        k2 = flux_f2( u[0,i+1], max(u[1,i+1], 0.0), u[2,i+1])
        s = g*u[0,i]* (B(x[i-2]) - B(x[i-3]) )/h if k1 >= k2 else g*u[0,i+1]* (B(x[i-1]) - B(x[i-2]) )/h
        res[1,i] += h * s
    return res
time_schemes = {'euler': apply_euler}
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
        line1.set_ydata(u[0,2:nc+2]/u[2,2:nc+2])
        line2.set_ydata(u[1,2:nc+2]/u[0,2:nc+2])
        line3.set_ydata(u[2,2:nc+2])
        fig.suptitle('nc='+str(nc)+', CFL='+str(cfl)+', time='+str(np.round(t, 3)))
        axes[2].set_ylim([np.min(u[2,2:nc+2]-0.1), np.max(u[2,2:nc+2])+0.1])
        plt.draw(); plt.pause(0.1)

# save final time solution to a file
fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u[0,2:nc+2]/u[2,2:nc+2], u[1,2:nc+2]/u[0,2:nc+2], u[2,2:nc+2]]))
print('Saved file ', fname)

if args.plot_freq >0:
    plt.show()



