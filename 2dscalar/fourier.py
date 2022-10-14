"""
2d fourier stability analysis
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from basis import *
from scipy.linalg import eigvals

# 1d CFL for radau case
cfl2_radau= [1.0, 0.333, 0.170, 0.103, 0.069] # for diss=2

# 1d CFL for g2 case
cfl2_g2 = [1.0, 1.0, 0.333, 0.170, 0.103] # for diss=2

ncfl  = 100   # No of cfl in each direction
nwave = 100   # No of wave numbers in each direction
TOL   = 1.0e-4 # tolerance to check eigenvalue

def R1(a):
    nr, nc = a.shape[0], a.shape[1]
    if nr != nc:
        print('The argument of R1 is not a square matrix')
        exit()
    return np.kron(np.eye(nr), a)

def R2(a):
    nr, nc = a.shape[0], a.shape[1]
    if nr != nc:
        print('The argument of R2 is not a square matrix')
        exit()
    return np.kron(a.T, np.eye(nr))

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-degree', type=int, choices=(1,2,3,4),
                    help='Polynomial degree', default=1)
parser.add_argument('-points', choices=('gl', 'gll'), help='Solution points',
                    default='gl')
parser.add_argument('-corr', choices=('radau', 'g2'), help='Correction function',
                    default='radau')
args = parser.parse_args()

if args.corr == 'radau':
    from gradau import *
    print('Correction function: radau')
elif args.corr == 'g2':
    from g2 import *
    print('Correction function: g2')
else:
    print('Corrector function not implemented.')
    exit()

k  = args.degree      # polynomial degree
print('Degree: ',k)

nd = k + 1 # dofs per cell

if args.corr == 'radau':
    cfl_max = 1.1 * cfl2_radau[k]
else:
    cfl_max = 1.1 * cfl2_g2[k]

bl, br = np.zeros((nd, 1)), np.zeros((nd, 1))
Vl, Vr = np.zeros((nd, 1)), np.zeros((nd, 1))

xg,wg,Vl[:,0],Vr[:,0],bl[:,0],br[:,0],Dm,D1 = fr_operators('fr',args.points,k,dgl,dgr)
I = np.eye(nd*nd)

# Absolute maximum value of the eigen values of the amplification matrix
# Uses global vars: k, bl, br, Vl, Vr, Dm, D1, I
def max_eig(sigma_x, sigma_y, k1dx, k2dy):
    H1 = - sigma_x * R1(Dm) - sigma_y * R2(Dm.T)
    if k == 1:
        T = I + 0.5 * H1
    elif k == 2:
        T = I + 0.5 * H1 + (1.0/6.0) * (H1 @ H1)
    elif k == 3:
        T = I + 0.5 * H1 + (1.0/6.0) * (H1 @ H1) + (1.0/24.0) * (H1 @ H1 @ H1)
    elif k == 4:
        T = I + 0.5 * H1 + (1.0/6.0) * (H1 @ H1) + (1.0/24.0) * (H1 @ H1 @ H1) \
              + (1.0/120.0) * (H1 @ H1 @ H1 @ H1)
    else:
        print('Not implemented for degree =', k)
        exit(0)
    Al = - sigma_x * R1(np.outer(bl, Vr.T)) @ T
    Ab = - sigma_y * R2(np.outer(Vr, bl.T)) @ T
    Ae = I - sigma_x * R1(D1) @ T - sigma_y * R2(D1.T) @ T \
           - sigma_x * R1(np.outer(br, Vr.T)) @ T \
           - sigma_y * R2(np.outer(Vr, br.T)) @ T
    A = Al * np.exp(-1j*k1dx) + Ae + Ab * np.exp(-1j*k2dy)
    return np.abs(eigvals(A)).max()

def is_stable(sigma_x, sigma_y):
    K = np.linspace(0,2.0*np.pi,nwave)
    for k1dx in K:
        for k2dy in K:
            if max_eig(sigma_x,sigma_y,k1dx,k2dy) - 1.0 > TOL: # unstable
                return False
    return True

sigma_x_range = np.linspace(0.0, cfl_max, ncfl)
sigma_y_range = np.linspace(0.0, cfl_max, ncfl)
X, Y = [], []
A = [] # All sigma for which (sigma,sigma) is a stable pair.
for sigma_x in sigma_x_range:
    for sigma_y in sigma_y_range:
        status = is_stable(sigma_x,sigma_y)
        if status == True: # stable
            X.append(sigma_x); Y.append(sigma_y)
            if np.abs(sigma_x - sigma_y) < 1e-4:
                A.append(sigma_x)
        else: # unstable
            # No need to check bigger values of sigma_y
            # break out of the sigma_y loop
            break

fname = 'fourier_k'+str(k)+'_'+args.corr+'.txt'
np.savetxt(fname, np.column_stack([X,Y]))
print('Saved ',fname)

cfl = round(2.0 * np.max(A), 3) # 3 decimal places
print('Highest sigma for which (sigma,sigma) is stable pair is ', np.max(A))
plt.figure(figsize=(5,5))
plt.scatter(X,Y,c='y',label='Stable Region')
#plt.plot(A,A,label="$\sigma_x + \sigma_y$ = CFL")
plt.plot([cfl,0],[0,cfl],label="$\sigma_x + \sigma_y$ = CFL")
plt.title('LWFR with '+str(args.corr)+', Degree = '+str(k)+', CFL ='+str(cfl))
plt.xlabel('$\sigma_x$'); plt.ylabel('$\sigma_y$')
plt.axis('equal'); plt.legend(); plt.grid(True)
fname = 'fourier_k'+str(k)+'_'+args.corr+'.pdf'
plt.savefig(fname); print('Saved ',fname)
