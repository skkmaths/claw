from basis import Legendre, dLegendre

# FR correction functions
# x is in [-1,+1]
def gl(k, x):
    return 0.5 * (-1)**k * (Legendre(k,x) - Legendre(k+1,x))

def gr(k, x):
    return 0.5 * (Legendre(k,x) + Legendre(k+1,x))

# Derivatives of FR correction functions
# x is in [-1,+1]
def dgl(k, x):
    return 0.5 * (-1)**k * (dLegendre(k,x) - dLegendre(k+1,x))

def dgr(k, x):
    return 0.5 * (dLegendre(k,x) + dLegendre(k+1,x))
