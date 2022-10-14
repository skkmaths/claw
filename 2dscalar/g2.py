from basis import Legendre, dLegendre

# FR correction functions for g2
# x is in [-1,+1]

# Left corrector function of degree k+1
def gl(k, x):
    value = 0.5 * (-1)**k * (Legendre(k,x) - ((k+1.0)*Legendre(k-1,x) + 
                                              k*Legendre(k+1,x))/(2.0*k+1.0))
    return value


# Right corrector function of degree k+1 defined by symmetry around origin
def gr(k, x):
    return gr(k,-x)

# Derivatives of FR correction functions
# x is in [-1,+1]

# Derivative of left corrector function of degree k+1
def dgl(k, x):
    value = 0.5 * (-1)**k * (1.0 - x) * dLegendre(k,x)
    return value

def dgr(k, x):
    return -dgl(k,-x)
