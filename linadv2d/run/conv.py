import os
import numpy as np
from matplotlib import pyplot as plt

# path to executable
exe = '../src/advect'
# check that executable is present
if os.path.isfile(exe)==False:
    print ("Could not find ", exe)
    exit()

N = [50, 100, 200, 400]
itmax = 50000
itsave = 100
fluxtype = "lw"
limtype  = "minmod"
method   = "method1"
cfl = 0.5
testcase = 1
dt = -1.0

os.system("rm -f err.dat")
for n in N:
    print ("n =",n)
    # create input file
    f = open('inp.dat','w')
    f.write(str(n)+'\n')
    f.write(str(n)+'\n')
    f.write(str(itmax)+'\n')
    f.write(str(itsave)+'\n')
    f.write(fluxtype+'\n')
    f.write(limtype+'\n')
    f.write(method+'\n')
    f.write(str(cfl)+'\n')
    f.write(str(dt)+'\n')
    f.write(str(testcase)+'\n')
    f.close()

    # run solver
    os.system(exe+" > out")
    os.system("tail -n1 out")
    os.system("tail -n1 out >> err.dat")

# read error from file
def f(x):
    return x**2

data = np.loadtxt('err.dat')
plt.loglog(data[:,0],data[:,4],'o-',  linewidth=2)
plt.loglog(data[:,0],data[:,5],'*--',linewidth=2)
plt.loglog(data[:,0],data[:,6],'s-.',linewidth=2)
plt.loglog(data[:,0],f(data[:,0]),'s-.',linewidth=2)

plt.xlabel('N')
plt.ylabel('Error')
plt.legend(('L1','L2','Linf'))
plt.savefig('error.pdf')
plt.show()
