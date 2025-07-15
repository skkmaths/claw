"""
number of cells =nx X ny
Solve scalar conservation law with periodic bc
To get help, type
    python lwfr.py -h
"""
import os, glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys
# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pde', choices=('linear', 'varadv', 'burger', 'bucklev'),
                    help='PDE', default='linear')
parser.add_argument('-scheme', choices=('lw','rk2','fo' ), help='lw',
                    default='fo')
parser.add_argument('-corr', choices=('radau', 'g2'), help='Correction function',
                    default='radau')
parser.add_argument('-ncellx', type=int, help='Number of x cells', default=50)
parser.add_argument('-ncelly', type=int, help='Number of y cells', default=50)
parser.add_argument('-cfl', type=float, help='CFL number', default=1.0)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
parser.add_argument('-ic', choices=('sin2pi', 'expo','hat', 'solid'),
                    help='Initial condition', default='sin2pi')
parser.add_argument('-limit', choices=('no', 'mmod', 'mc'), help='Apply limiter',
                    default='no')
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-compute_error', choices=('no', 'yes'),
                    help='Compute error norm', default='no')
parser.add_argument('-save_freq', type=int, help='Frequency to save solution',
                    default=0)
args = parser.parse_args()

# Select PDE
if args.pde == 'linear':
    from linadv import *
elif args.pde == 'varadv':
    from varadv import *
elif args.pde == 'burger':
    from burger import *
else:
    print('PDE not implemented')
    exit()

# Select initial condition
if args.ic == 'sin2pi':
    from sin2pi import *
elif args.ic == 'expo':
    from expo import *
elif args.ic == 'hat':
    from hat import *
elif args.ic == 'solid':
    from solid import *
else:
    print('Unknown initial condition')
    exit()

# Select cfl
cfl = args.cfl
beta = 1.0 # parameter in minmod

nx = args.ncellx       # number of cells in the x-direction
ny = args.ncelly       # number of cells in the y-direction
global fileid
fileid = 0
dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny
# Allocate solution variables
v = np.zeros((nx+4, ny+4))  # 2 ghost cells each side
vres = np.zeros((nx+4, ny+4))  # 2 ghost cells each sideresidual
# To store the cell averages only in real cells.
# Set initial condition by interpolation
for i in range(nx+4):
    for j in range(ny+4):
        x = xmin + (i-2)*dx+0.5 * dx     
        y = ymin + (j-2)*dy + 0.5 * dy 
        val = initial_condition(x, y)
        v[i, j] = val
# copy the initial condition
# index 2 to nx+1 and 2 to ny+1
v0 = v[2:nx+2, 2:ny+2].copy()
# it stores the coordinates of real cell centre
xgrid1 = np.linspace(xmin+0.5*dx, xmax-0.5*dx, nx)
ygrid1 = np.linspace(ymin+0.5*dy, ymax-0.5*dy, ny)
ygrid, xgrid = np.meshgrid(ygrid1, xgrid1)

# it stores the coordinates of vertices of the real cells
Xgrid = np.linspace(xmin, xmax, nx+1)
Ygrid = np.linspace(ymin, ymax, ny+1)
Ygrid, Xgrid = np.meshgrid(Ygrid, Xgrid)
#------------To save solution--------------------------------------------
def getfilename(file, fileid):
    if fileid <10:
        file = file +"00"+str(fileid)+".plt"
    elif fileid <99:
        file = file +"0"+str(fileid)+".plt"
    else:
        file =file+str(fileid)+".plt"
    return file
# save solution to a file
def savesol(t, var_u):
    global fileid
    if not os.path.isdir("sol"): # creat a dir if not
       os.makedirs("sol")
       print('Directory "sol" is created')
    if fileid == 0: # remove the content of the folder
        print('The directory "sol" is going to be formated!')
        if input('Do You Want To Continue? [y/n] ') != 'y':
            sys.exit('Execution is treminated')
        fs = glob.glob('./sol/*')
        for f in fs:
           os.remove(f)
    filename = "sol"
    filename = getfilename(filename, fileid)
    file = open("./sol/"+filename,"a")
    file.write('TITLE = "Linear advectino equation" \n')
    file.write('VARIABLES = "x", "y", "sol" \n')
    file.write("ZONE STRANDID=1, SOLUTIONTIME= "+ str(t)+ ", I= "+str(nx)+", J ="+str(ny)+", DATAPACKING=POINT \n")
    for j in range(2, ny+2):
        for i in range(2, nx+2):
            x = xmin + (i-2)*dx + 0.5*dx
            y = ymin + (j-2)*dy + 0.5*dy
            file.write( str(x) + ", " + str(y) +"," + str(var_u[i,j])+"\n")
    file.close()
    fileid = fileid + 1
#-----------------------------------------------------------------------------
def minmod(a,b,c):
    sa = np.sign(a)
    sb = np.sign(b)
    sc = np.sign(c)
    if sa==sb and sb==sc:
        return sa * np.abs([a,b,c]).min()
    else:
        return 0.0
# mc limiter
def mc(a, b, c):
    if a * b <= 0.0:
        return 0.0
    else:
        min_val = min(2.0 * abs(a), 2.0 * abs(b), abs(c))
        return min_val if c >= 0 else -min_val

def reconstruct(conjm1, conj, conjp1):
    if args.limit == 'no':
        return conj
    elif args.limit == 'mmod':
        conl = conj + 0.5 * minmod( beta*(conj-conjm1), \
                                0.5*(conjp1-conjm1), \
                                beta*(conjp1-conj) )
        return conl
    elif args.limit == 'mc':
        conl = conj + 0.5 * mc( beta*(conj-conjm1), \
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

def init_plot(ax1, ax2, ax3, u0):
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
    
    line1,line2 = ax3.plot(xgrid1, np.diag(u0),'ro-', xgrid1, np.diag(u0),'b')
    plt.legend(('exact','approx'))
    ax3.set_xlabel('x')
    ax3.set_ylabel('v')
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
              ', t = '+str(round(t, 3)))
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(cp)
    ax3 = fig.add_subplot(122)
    line1,line2 = ax3.plot(xgrid1, np.diag(v0),'ro-', xgrid1, np.diag(u1),'b')
    plt.legend(('exact','approx'))
    ax3.set_xlabel('x')
    ax3.set_ylabel('v')
    ax3.set_title('Solution along the line x=y')
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
    fig = plt.figure(figsize = (16, 7))
    #ax1 = fig.add_subplot(121,projection='3d')
    #ax2 = fig.add_subplot(122)
    ax2 = fig.add_subplot(121)
    ax3 = fig.add_subplot(122)
    #init_plot(ax1, ax2, v[1:nx+1,1:ny+1])
    init_plot(ax2, ax2, ax3, v0)
    wait = input("Press enter to continue ")

#Update solution using RK time scheme
def apply_ssprk22 ( t, dt, lam_x, lam_y, v_old, v, vres):
    #first stage
    ts  = t
    update_ghost(v)
    vres = compute_residual(ts, lam_x, lam_y, v, vres)
    v = v - vres

    #second stage
    ts = t + dt
    update_ghost(v)
    vres = compute_residual(ts, lam_x, lam_y, v, vres)
    v = 0.5 * v_old + 0.5 * (v - vres)
    
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

if (args.pde == 'linear' or args.pde == 'varadv'):
    # Find dt once since cfl does not depend on u or time
    sx, sy = local_speed(xgrid, ygrid, v[2:nx+2,2:ny+2])
    #  |sigma_x| + |sigma_y| = cfl
    if ( args.scheme == 'lw'):
        dt = 0.72/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()
    elif (args.scheme == 'fo' or args.scheme == 'rk2'):
        dt = cfl/(np.abs(sx)/dx + np.abs(sy)/dy + 1.0e-14).max()

it, t = 0, 0.0
Tf = args.Tf
#save initial data
if args.save_freq > 0:
    savesol(t, v)
while t < Tf:
    if (args.pde == 'burger'):
        dt = cfl/(max_speed(v)/dx +max_speed(v)/dy)
    if t+dt > Tf:
        dt = Tf - t
    lamx, lamy = dt/dx,  dt/dy
    # Loop over real cells (no ghost cell) and compute cell integral
    v_old = v.copy()
    if args.scheme == 'lw':
        update_ghost(v)
        vres = compute_residual_lw(t, dt, lamx, lamy, v, vres)
        v = v - vres
    elif args.scheme == 'rk2':
        v = apply_ssprk22 ( t, dt, lamx, lamy, v_old, v, vres)
    elif args.scheme == 'fo':
        update_ghost(v)
        vres = compute_residual(t,lamx, lamy, v, vres)
        v = v - vres
        
    t, it = t+dt, it+1
    if args.save_freq > 0:
        if it % args.save_freq == 0:
            savesol(t, v)
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
# print final data
if args.save_freq > 0:
    savesol(t,v)
    print('it,t,min,max =', it, t, v[2:nx+2,2:ny+2].min(), v[2:nx+2,2:ny+2].max())
    print('solution saved to .plt files')


if args.plot_freq > 0: 
    plt.show()

# Things to do
# --------------------------------
# Save cell average solution to file
# Dirichlet boundary condition
