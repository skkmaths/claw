"""
number of cells =nx X ny
Solve scalar conservation law with periodic bc
To get help, type
    python lwfr.py -h
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib import cm
from matplotlib.ticker import LinearLocator

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pde', choices=('linear', 'varadv', 'burger', 'bucklev'),
                    help='PDE', default='linear')
parser.add_argument('-scheme', choices=('fr', 'dfr'), help='FR or DFR',
                    default='fr')
parser.add_argument('-corr', choices=('radau', 'g2'), help='Correction function',
                    default='radau')
parser.add_argument('-points', choices=('gl', 'gll'), help='Solution points',
                    default='gl')
parser.add_argument('-ncellx', type=int, help='Number of x cells', default=50)
parser.add_argument('-ncelly', type=int, help='Number of y cells', default=50)
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-cfl', type=float, help='CFL number', default=1.0)
parser.add_argument('-diss', type=int, choices=(1, 2), help='Dissipation type',
                    default=1)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
parser.add_argument('-ic', choices=('sin2pi', 'expo'),
                    help='Initial condition', default='sin2pi')
parser.add_argument('-limit', choices=('no', 'tvb', 'blend'), help='Apply limiter',
                    default='no')
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-compute_error', choices=('no', 'yes'),
                    help='Compute error norm', default='no')
args = parser.parse_args()

# Select PDE
if args.pde == 'linear':
    from linadv import *
elif args.pde == 'varadv':
    from varadv import *
else:
    print('PDE not implemented')
    exit()

# Select initial condition
if args.ic == 'sin2pi':
    from sin2pi import *
elif args.ic == 'expo':
    from expo import *
else:
    print('Unknown initial condition')
    exit()

# Select cfl
cfl = args.cfl


nx = args.ncellx       # number of cells in the x-direction
ny = args.ncelly       # number of cells in the y-direction

dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny

# Allocate solution variables
v = np.zeros((nx+2, ny+2))  # solution at n+1
vres = np.zeros((nx+2, ny+2))  # residual

# To store the cell averages only in real cells.

# Set initial condition by interpolation
for i in range(nx+2):
    for j in range(ny+2):
        x = xmin + (i-1)*dx+0.5 * dx     # transform gauss points to real cell
        y = ymin + (j-1)*dy + 0.5 * dy 
        val = initial_condition(x, y)
        v[i, j] = val

# it stores the coordinates of real cell centre
xgrid = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nx)
ygrid = np.linspace(ymin+0.5*dy, ymax-0.5*dy, ny)
ygrid, xgrid = np.meshgrid(ygrid, xgrid)

# it stores the coordinates of vertices of the real cells

Xgrid = np.linspace(xmin, xmax, nx+1)
Ygrid = np.linspace(ymin, ymax, ny+1)
Ygrid, Xgrid = np.meshgrid(Ygrid, Xgrid)

def minmod(a,b,c):
    sa = np.sign(a)
    sb = np.sign(b)
    sc = np.sign(c)
    if sa==sb and sb==sc:
        return sa * np.abs([a,b,c]).min()
    else:
        return 0.0

# Initialize plot
def init_plot(ax1, ax2, u0):
    '''
    sp = ax1.plot_surface(xgrid, ygrid, u0, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax1.set_title('Initial condition')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_ylabel('v')
    '''
    cp = ax2.contour(xgrid, ygrid, u0, levels=16)
    ax2.set_title('Initial condition')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(cp)

    plt.draw()
    plt.pause(0.1)

# Update plot
def update_plot(fig, t, u1):
    plt.clf()
    '''
    ax1 = fig.add_subplot(121,projection='3d')
    sp = ax1.plot_surface(xgrid, ygrid, u1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax1.set_title(str(nx)+'X'+str(ny)+' cells, CFL = '+str(round(cfl, 3)) +
              ', Diss = '+str(args.diss)+', t = '+str(round(t, 3)))
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_ylabel('v')
    '''

    #ax2 = fig.add_subplot(122)
    ax2 = fig.add_subplot(111)
    cp = ax2.contour(xgrid, ygrid, u1, levels=16)
    ax2.set_title(str(nx)+'X'+str(ny)+' cells, CFL = '+str(round(cfl, 3)) +
              ', Diss = '+str(args.diss)+', t = '+str(round(t, 3)))
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(cp)

    plt.draw()
    plt.pause(0.1)
    plt.clf()



# Fill ghost cells using periodicity
def update_ghost():
    # left ghost cell
    v[0,:] = v[nx,:]
    # right ghost cell
    v[nx+1,:] = v[1,:]
    # bottom ghost cell
    v[:,0]= v[:,ny]
    # top ghost cell
    v[:,ny+1] = v[:,1]
fig = plt.figure()
#ax1 = fig.add_subplot(121,projection='3d')
#ax2 = fig.add_subplot(122)
ax2 = fig.add_subplot(111)
#init_plot(ax1, ax2, v[1:nx+1,1:ny+1])
init_plot(ax2, ax2, v[1:nx+1,1:ny+1])
wait = input("Press enter to continue ")
t =0.1

# Find dt once since cfl does not depend on u or time
sx, sy = local_speed(xgrid,ygrid,v[1:nx+1,1:ny+1])
# |sigma_x| + |sigma_y| = cfl
dt = cfl/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()
#dt = 0.3 * dx
# function to compute the residual v^n+1_i = v^n_i -res_i
def compute_residual(t, dt,lam_x, lam_y, v, vres):
    vres[:,:] = 0.0
    update_ghost()  # Fill the ghost cell with values.
    # compute the inter-cell fluxes
    # loop over interior  vertical faces
    for i in range(0, nx+1):
        xf = (xmin + i*dx)  # location of this face
        for j in range(1, ny+1):
            y = ymin + (j-1)*dy+ 0.5*dy
            vl, vr = v[i, j], v[i+1, j]
            Fn = xnumflux(xf, y, vl, vr, vl , vr)
            vres[i, j] += lamx*Fn
            vres[i+1, j] -= lamx*Fn
    # loop over interior horizontal faces
    for j in range(0, ny+1):
        yf = (ymin + j*dy)
        for i in range(1, nx+1):
            x = xmin + (i-1)*dx + 0.5 * dx
            vl, vr = v[i, j], v[i, j+1]
            Gn = ynumflux(x, yf, vl, vr, vl, vr)
            vres[i, j] += lamy*Gn
            vres[i, j+1] -= lamy*Gn
    return vres

it, t = 0, 0.0
Tf = args.Tf
while t < Tf:
    if t+dt > Tf:
        dt = Tf - t
    lamx, lamy = dt/dx,  dt/dy
    # Loop over real cells (no ghost cell) and compute cell integral
    v_old = v
    vres = compute_residual(t,dt, lamx, lamy, v, vres )
    v = v - vres
    t, it = t+dt, it+1
    print('it,t,min,max =', it, t, v[1:nx+1,1:ny+1].min(), v[1:nx+1,1:ny+1].max())
    if it % args.plot_freq == 0 or np.abs(Tf-t) < 1.0e-13:
        update_plot(fig, t, v[1:nx+1,1:ny+1])

plt.show()
# Things to do
# --------------------------------
# Save cell average solution to file
# Dirichlet boundary condition
