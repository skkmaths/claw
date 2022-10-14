#from numpy import sqrt
import numpy as np
from lgl import *

# Legendre polynomials on [-1,+1]
def Legendre(n, x):
    if n==0:
        value = 1.0
    elif n==1:
        value = x
    else:
        value = (2.0*n-1.0)/n * x * Legendre(n-1,x) - (n-1.0)/n * Legendre(n-2,x)
    return value

# Derivative of Legendre
def dLegendre(n, x):
    if n==0:
        value = 0.0
    elif n==1:
        value = 1.0
    else:
        value = n * Legendre(n-1,x) + x * dLegendre(n-1,x)
    return value

# Normalized Legendre polynomials
def nLegendre(n,x):
    return np.sqrt(2*n+1) * Legendre(n,x)

# xp : set of grid points
# Returns i'th Lagrange polynomial value at x
def Lagrange(i, xp, x):
    value = 1.0
    n = len(xp)
    for j in range(n):
        if j != i:
            value = value * (x - xp[j]) / (xp[i] - xp[j])
    return value

# Vandermonde matrix for Lagrange polynomials
# xp: grid points
# x : evaluation points
def Vandermonde_lag(xp, x):
    n = len(xp)
    m = len(x)
    V = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            V[i,j] = Lagrange(j, xp, x[i])
    return V

# Vandermonde matrix for Legendre polynomials
# k : degree
# x : evaluation points, in [0,1]
def Vandermonde_leg(k, x):
    n = k + 1
    m = len(x)
    V = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            V[i, j] = nLegendre(j, 2.0*x[i]-1.0)
    return V

# Projection matrix: nodal --> modal
# xg must be in [0,1]
def nodal2modal(xg):
    n = len(xg)
    k = n - 1

    nq = k + 1
    x, w = np.polynomial.legendre.leggauss(nq)
    x = 0.5*(x + 1.0)
    w = 0.5*w
    
    Vleg = Vandermonde_leg(k, x)
    Vlag = Vandermonde_lag(xg, x)

    # Mass matrix
    M = np.zeros(n)
    for i in range(n):
        M[i] = np.sum(Vleg[:,i] * Vleg[:,i] * w)
    err = np.abs(M - np.ones(n)).max()
    if err > 1e-10:
        print('Legendre mass matrix =',M)
        exit()

    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = np.sum(Vleg[:,i] * Vlag[:,j] * w)

    for i in range(n):
        A[i,:] = A[i,:] / M[i]
    
    return A

# Barycentric weights in lagrange interpolation
def barycentric_weights(x):
    n = len(x)
    w = np.ones(n)

    for j in range(1,n):
        for k in range(0,j):
            w[k] = w[k] * (x[k] - x[j])
            w[j] = w[j] * (x[j] - x[k])

    return 1.0/w

# Differentiation matrix
# D[i,j] = l_j'(x_i)
def diff_mat(x):
    w = barycentric_weights(x)
    n = len(x)
    D = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            if j != i:
                D[i,j] = (w[j]/w[i]) * 1.0/(x[i] - x[j])
                D[i,i] = D[i,i] - D[i,j]
    return D

# scheme = fr, dfr
# points = gl, gll
# degree = 1,2,3,4,5
def fr_operators(scheme,points,degree,dgl,dgr):
    k = degree
    nd = k + 1
    if points == 'gl':
        # k+1 point gauss rule, integrates exactly upto degree 2*k+1
        xg, wg = np.polynomial.legendre.leggauss(nd)
    elif points == 'gll':
        if scheme == 'dfr':
            print('Cannot use GLL with DFR')
            exit()
        xg, wg = lglnodes(nd-1)
    else:
        print('Unknown points = ',points)
        exit()

    # Transform to [0,1]
    xg = 0.5*(xg + 1.0)
    wg = 0.5*wg

    # Required to evaluate solution at face
    Vl, Vr = np.zeros(nd), np.zeros(nd)
    for i in range(nd):
        Vl[i] = Lagrange(i, xg, 0.0)
        Vr[i] = Lagrange(i, xg, 1.0)

    # Differentiation matrix
    Dm = diff_mat(xg)

    if scheme == 'dfr':
        # Extended points for flux in DFR
        xe = np.zeros(k+3)
        xe[0] = 0.0
        xe[1:-1] = xg
        xe[-1] = 1.0
        Dt = diff_mat(xe)[1:-1, :]  # derivative not needed at x=0 and x=1
        bl = Dt[:, 0]   # first col
        br = Dt[:, -1]  # last co
        D1 = Dt[:, 1:-1]
    elif scheme == 'fr':
        bl, br = np.zeros(nd), np.zeros(nd)
        for i in range(nd):
            bl[i] = 2.0 * dgl(k, 2*xg[i]-1)
            br[i] = 2.0 * dgr(k, 2*xg[i]-1)
        D1 = Dm - np.outer(bl, Vl) - np.outer(br, Vr)
    else:
        print('Unknown scheme = ',scheme)
        exit()

    return xg,wg,Vl,Vr,bl,br,Dm,D1
