"""
number of cells =nx X ny
Solve scalar conservation law with periodic bc
To get help, type
    python lwfr.py -h
"""
from this import d
import numpy as np
import matplotlib.pyplot as plt
import argparse
from basis import *
from constants import *
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
parser.add_argument('-compute_error', choices=('no', 'yes'),
                    help='Compute error norm', default='no')
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
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
theta = 0.5  # constant in the minmods


nx = args.ncellx       # number of cells in the x-direction
ny = args.ncelly       # number of cells in the y-direction

dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny
Mdx2 = args.tvbM*dx**2

# Allocate solution variables
v = np.zeros((nx+2, ny+2))  # solution at n+1
vres = np.zeros((nx+2, ny+2))  # residual
s_x = np.zeros((nx+2, ny+2)) # slopes in x directin in each cell
s_y = np.zeros((nx+2, ny+2)) # slopes in y direction in each cell

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


def minmod(a,b,c,Mdx2):
    val = 0.0
    if np.abs(a) < Mdx2:
            val = a
    else:
        sa = np.sign(a)
        sb = np.sign(b)
        sc = np.sign(c)
        if sa==sb and sb==sc:
            val = sa * np.abs([a, b, c]).min()
        else:
            val = 0.0
    return val

def compute_slopes():
    s_x[:,:] = 0.0
    s_y[:,:] = 0.0
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            vl, vr = v[i,j-1], v[i,j+1]
            dvl = v[i,j]-vl
            dvr = vr-v[i,j]
            dvc = vr-vl
            s_y[i,j] = 2.0*theta*minmod(dvl, 0.5*dvc, dvr,Mdx2)
    
    for j in range (1,ny+1):
        for i in range(1,nx+1):
            vl, vr = v[i-1,j], v[i+1,j]
            dvl = v[i,j] - v[i-1,j]
            dvr = vr - v[i,j]
            dvc = vr-vl
            s_x[i,j] = 2.0*theta*minmod(dvl, 0.5*dvc, dvr,Mdx2)
    
    # slopes in ghost cell fro periodic case
    # x direction
    # left ghost cell
    s_x[0,:] = s_x[nx,:]
    # right ghost cell
    s_x[nx+1,:] = s_x[1,:]
    # bottom ghost cell
    s_x[:,0]= s_x[:,ny]
    # top ghost cell
    s_x[:,ny+1] = s_x[:,1]
   
    # y direction
    # left ghost cell
    s_y[0,:] = s_y[nx,:]
    # right ghost cell
    s_y[nx+1,:] = s_y[1,:]
    # bottom ghost cell
    s_y[:,0]= s_y[:,ny]
    # top ghost cell
    s_y[:,ny+1] = s_y[:,1]

# Initialize plot
def init_plot(ax1, ax2, u0):
    sp = ax1.plot_surface(xgrid, ygrid, u0, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax1.set_title('Initial condition')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_ylabel('v')

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
    ax1 = fig.add_subplot(121,projection='3d')
    sp = ax1.plot_surface(xgrid, ygrid, u1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax1.set_title(str(nx)+'X'+str(ny)+' cells, CFL = '+str(round(cfl, 3)) +
              ', Diss = '+str(args.diss)+', t = '+str(round(t, 3)))
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_ylabel('v')

    ax2 = fig.add_subplot(122)
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
ax1 = fig.add_subplot(121,projection='3d')
ax2 = fig.add_subplot(122)
init_plot(ax1, ax2, v[1:nx+1,1:ny+1])
wait = input("Press enter to continue ")
t =0.1

# Find dt once since cfl does not depend on u or time
sx, sy = local_speed(xgrid,ygrid,v[1:nx+1,1:ny+1])
# |sigma_x| + |sigma_y| = cfl
dt = cfl/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()
#dt = 0.3 * dx
# function to compute the residual v^n+1_i = v^n_i -res_i
def compute_residual(v, vres):
    update_ghost()  # Fill the ghost cell with values.
    compute_slopes()
    # compute the inter-cell fluxes
    # loop over interior  vertical faces
    for i in range(1, nx+2):
        xf = (xmin + (i-1)*dx)  # location of this face
        for j in range(1, ny+1):
            y = ymin + (j-1)*dy+ 0.5*dy
            vl, vr = v[i-1, j] + 0.5 * s_x[i-1,j], v[i, j] -0.5 * s_x[i,j]
            Fn = xnumflux(xf, y, vl, vr, vl, vr)
            vres[i-1, j] += lamx*Fn
            vres[i, j] -= lamx*Fn
    # loop over interior horizontal faces
    for j in range(1, ny+2):
        yf = (ymin + (j-1)*dy)
        for i in range(1, nx+1):
            x = xmin + (i-1)*dx + 0.5 * dx
            vl, vr = v[i, j-1] + 0.5*s_y[i,j-1], v[i, j] -0.5 * s_y[i,j]
            Gn = ynumflux(x, yf, vl, vr, vl, vr)
            vres[i, j-1] += lamy*Gn
            vres[i, j] -= lamy*Gn
    return vres

def apply_ssprk22(t, dt, v_old, v, vres ):
    #first stage
    ts  = t
    vres = compute_residual(v, vres)
    v = v - vres
    #second stage
    ts = t + dt
    vres = compute_residual(v, vres)
    v = 0.5 * v_old + 0.5 *(v - vres)
    return  v

it, t = 0, 0.0
Tf = args.Tf
while t < Tf:
    if t+dt > Tf:
        dt = Tf - t
    lamx, lamy = dt/dx,  dt/dy
    # Loop over real cells (no ghost cell) and compute cell integral
    vres[:,:] = 0.0
    v_old = v
    v = apply_ssprk22(t, dt, v_old, v, vres )
    t, it = t+dt, it+1
    print('it,t,min,max =', it, t, v[1:nx+1,1:ny+1].min(), v[1:nx+1,1:ny+1].max())
    if it % args.plot_freq == 0 or np.abs(Tf-t) < 1.0e-13:
        update_plot(fig, t, v[1:nx+1,1:ny+1])
    
plt.show()
# Things to do
# --------------------------------
# Save cell average solution to file
# Dirichlet boundary condition
