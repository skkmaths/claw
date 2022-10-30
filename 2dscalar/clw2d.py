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
parser.add_argument('-scheme', choices=('lw','rk2','fo' ), help='lw',
                    default='fo')
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
parser.add_argument('-limit', choices=('no', 'mmod'), help='Apply limiter',
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
beta = 2.0 # parameter in minmod
nx = args.ncellx       # number of cells in the x-direction
ny = args.ncelly       # number of cells in the y-direction

dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny

# Allocate solution variables
v = np.zeros((nx+4, ny+4))  # 2 ghost cells each side
vres = np.zeros((nx+4, ny+4))  # 2 ghost cells each sideresidual

# To store the cell averages only in real cells.

# Set initial condition by interpolation
for i in range(nx+4):
    for j in range(ny+4):
        x = xmin + (i-2)*dx+0.5 * dx     # transform gauss points to real cell
        y = ymin + (j-2)*dy + 0.5 * dy 
        val = initial_condition(x, y)
        v[i, j] = val
# copy the initial condition
# index 2 to nx+1 and 2 to ny+1
v0 = v[2:nx+2, 2:ny+2].copy()
# it stores the coordinates of real cell centre
xgrid = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nx)
ygrid = np.linspace(ymin+0.5*dy, ymax-0.5*dy, ny)
x1, y1 = xgrid,  ygrid # to plot 1d graph
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

def reconstruct(conjm1, conj, conjp1):
    if args.limit == 'no':
        return conj
    elif args.limit == 'mmod':
        conl = conj + 0.5 * minmod( beta*(conj-conjm1), \
                                0.5*(conjp1-conjm1), \
                                beta*(conjp1-conj) )
        return conl
    else:
        print('limit type not define')
        exit()
# Initialize plot
def lwflux(u,v,dx,dy,dt,ql,qr,qdl,qdr,qul,qur):
    qt = - u * (qr - ql)/dx - v * (qul - qdl + qur - qdr)/(4.0*dy)
    flux = u * (ql + qr)/2.0  + 0.5 * dt * u * qt
    return flux

def init_plot(ax1, ax2,ax3, u0):
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

    line1,line2 = ax3.plot(y1,u0[ny-1,0:ny],'ro', y1, u0[ny-1,0:ny],'b')
    plt.grid(True);
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
    ax2 = fig.add_subplot(121)
    cp = ax2.contour(xgrid, ygrid, u1, levels=16)
    ax2.set_title(str(nx)+'X'+str(ny)+' cells, CFL = '+str(round(cfl, 3)) +
              ', Diss = '+str(args.diss)+', t = '+str(round(t, 3)))
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(cp)
    ax3 = fig.add_subplot(122)
    line1,line2 = ax3.plot(y1,u1[ny-1,0:ny],'ro', y1, v0[ny-1,0:ny],'b')
    plt.grid(True);
    plt.draw()
    plt.pause(0.1)
    plt.clf()


# Fill ghost cells using periodicity
def update_ghost(v1):
    # left ghost cell
    v1[0,:] = v1[nx,:]
    v1[1,:] = v1[nx+1,:]
    # right ghost cell
    v1[nx+3,:] = v1[3,:]
    v1[nx+2,:] = v1[2,:]
    # bottom ghost cell
    v1[:,0]= v1[:,ny]
    v1[:,1]= v1[:,ny+1]
    # top ghost cell
    v1[:,ny+2] = v1[:,2]
    v1[:,ny+3] = v1[:,3]

if args.plot_freq > 0:
    fig = plt.figure()
    #ax1 = fig.add_subplot(121,projection='3d')
    #ax2 = fig.add_subplot(122)
    ax2 = fig.add_subplot(121)
    ax3 = fig.add_subplot(122)
    #init_plot(ax1, ax2, v[1:nx+1,1:ny+1])
    init_plot(ax2, ax2, ax3, v[2:nx+2,2:ny+2])
    wait = input("Press enter to continue ")



#Update solution using RK time scheme
def apply_ssprk22 ( t, dt, lam_x, lam_y, v_old, v, vres):
    #first stage
    ts  = t
    vres = compute_residual(ts, lam_x, lam_y, v, vres)
    v = v - vres
    update_ghost(v)
    
    #second stage
    ts = t + dt
    vres = compute_residual(ts, lam_x, lam_y, v, vres)
    v = 0.5 * v_old + 0.5 * (v - vres)
    update_ghost(v)
    return v
# Residual for Lax-Wendroff scheme in fv conservative form
def compute_residual_lw(t, dt, lam_x, lam_y, v, vres):
    vres[:,:] = 0.0
    # compute the inter-cell fluxes
    # loop over interior  vertical faces
    for i in range(1, nx+2):  # face between (i,j) and (i+1,j)
        xf = (xmin + (i-1)*dx)  # x location of this face
        for j in range(2, ny+2):
            y = ymin + (j-2)*dy+ 0.5*dy # cetre of vertical face
            speed = advection_velocity(xf,y)
            Fn = lwflux(speed[0], speed[1], dx, dy, dt, \
                        v[i,j],  v[i+1,j], \
                        v[i,j-1], v[i+1,j-1], \
                        v[i,j+1],v[i+1,j+1])
            vres[i, j] += lamx*Fn
            vres[i+1, j] -= lamx*Fn
    # loop over interior horizontal faces
    for j in range(1, ny+2):
        yf = (ymin + (j-1)*dy)
        for i in range(2, nx+2):
            x = xmin + (i-2)*dx + 0.5 * dx
            speed = advection_velocity(x,yf)
            Gn = lwflux(speed[1], -speed[0], dy, dx, dt, \
                        v[i,j],  v[i,j+1], \
                        v[i+1,j],v[i+1,j+1], \
                        v[i-1,j],v[i-1,j+1])
            vres[i, j] += lamy*Gn
            vres[i,j+1] -= lamy*Gn
    return vres
# Residual for Lax-Wendroff scheme in fd form
def compute_residual_lxw(t, dt, lam_x, lam_y, v, res):
    for i in range(2, nx+2): # 2 to nx+1
        for j in range(2, ny+2):
            vres[i, j] = - 0.5 * lam_x * ( v[i+1,j]- v[i-1,j])- 0.5 * lam_y * (v[i,j+1]-v[i, j-1]) \
                          + 0.5 * lam_x**2 * (v[i-1,j] - 2.0*v[i,j] + v[i+1,j]) \
                          + 0.25 * lam_x * lam_y * (v[i+1,j+1] - v[i+1,j-1] - v[i-1,j+1] + v[i-1,j-1] ) \
                          + 0.5 * lam_y**2 * ( v[i,j-1] - 2.0*v[i,j] + v[i,j+1])
    return -vres
# Compute residual of fv  scheme
def compute_residual(t,lam_x, lam_y, v, vres):
    vres[:,:] = 0.0
    # compute the inter-cell fluxes
    # loop over interior  vertical faces
    for i in range(1, nx+2):  # face between (i,j) and (i+1,j)
        xf = (xmin + (i-1)*dx)  # x location of this face
        for j in range(2, ny+2):
            y = ymin + (j-2)*dy+ 0.5*dy # cetre of vertical face
            vl = reconstruct(v[i-1, j], v[i, j], v[i+1, j])
            vr = reconstruct(v[i+2, j], v[i+1, j], v[i, j])
            Fl, Fr = xflux(xf, y, vl), xflux(xf, y, vr)
            Fn = xnumflux(xf, y, Fl, Fr, vl, vr)
            vres[i, j] += lamx*Fn
            vres[i+1, j] -= lamx*Fn
    # loop over interior horizontal faces
    for j in range(1, ny+2):
        yf = (ymin + (j-1)*dy)
        for i in range(2, nx+2):
            x = xmin + (i-2)*dx + 0.5 * dx
            vl = reconstruct(v[i,j-1], v[i,j], v[i,j+1])
            vr = reconstruct(v[i,j+2], v[i,j+1],v[i,j])
            Gl, Gr = yflux(x,yf, vl), yflux(x,yf, vr)
            Gn = ynumflux(x, yf, Gl, Gr, vl, vr)
            vres[i, j] += lamy*Gn
            vres[i,j+1] -= lamy*Gn
    return vres

# Find dt once since cfl does not depend on u or time
sx, sy = local_speed(xgrid, ygrid, v[2:nx+2,2:ny+2])
# |sigma_x| + |sigma_y| = cfl
#dt = cfl/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()
if ( args.scheme == 'lw'):
    dt = 0.72/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()
elif (args.scheme == 'fo' or args.scheme == 'rk2'):
    dt = cfl/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()
dt = dx
it, t = 0, 0.0
Tf = args.Tf
while t < Tf:
    if t+dt > Tf:
        dt = Tf - t
    lamx, lamy = dt/dx,  dt/dy
    # Loop over real cells (no ghost cell) and compute cell integral
    v_old = v
    if args.scheme == 'lw':
        vres = compute_residual_lw(t, dt, lamx, lamy, v, vres)
        v = v - vres
        update_ghost(v)
    elif args.scheme == 'rk2':
        v = apply_ssprk22 ( t, dt, lamx, lamy, v_old, v, vres)
    elif args.scheme == 'fo':
        vres = compute_residual(t,lamx, lamy, v, vres)
        v = v - vres
        update_ghost(v)
    

    t, it = t+dt, it+1
    if args.plot_freq > 0:
        print('it,t,min,max =', it, t, v[2:nx+2,2:ny+2].min(), v[2:nx+2,2:ny+2].max())
        if it% args.plot_freq == 0:
            update_plot(fig, t, v[2:nx+2,2:ny+2])


# Compute error norm: initial condition is exact solution
if args.compute_error == 'yes':
    l1_err = np.sum(np.abs(v[2:nx+2,2:ny+2]-v0)) / (nx*ny)
    l2_err = np.sqrt(np.sum((v[2:nx+2,2:ny+2]-v0)**2) / (nx*ny))
    li_err = np.abs(v[2:nx+2,2:ny+2]-v0).max()
    print('dx,dy,l1,l2,linf error =')# %10.4e %10.4e %10.4e %10.4e %10.4e' % 
    print(dx,dy,l1_err,l2_err,li_err)
if args.plot_freq > 0: 
    plt.show()

# Things to do
# --------------------------------
# Save cell average solution to file
# Dirichlet boundary condition
