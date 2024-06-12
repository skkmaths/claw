from xml import dom
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from functools import partial
from scipy import optimize

def smooth(x):
    return np.sin(np.pi*x)
def dsmooth(x): # derivative of initial smooth function
    return np.cos(np.pi*x)*np.pi

def impf(u):
    return u- smooth(xx - tt*u)
def dimpf(u): # derivative of implicity function
    return 1.0 + tt* dsmooth(xx-tt*u)

x  = np.linspace(-1,1,200)
u  = smooth(x)

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, u)
plt.grid(True); plt.draw(); plt.pause(0.1)
\
Tf = 1.0
deltat = 0.01
time = 0 
uold = u
it = 0
while time < Tf:
    uold = u
    time += deltat
    tt = time
    for i in range(0, len(x)):
        guess = smooth(x[i])
        xx = x[i]
        u0 = u1 = guess #uold[i]  # or = guess
        for j in range ( 200):
            u0 = u1
            u1 = u0 - impf(u0)/dimpf(u0)
        u[i] = u1
    it += 1
    print('Iteration = ', it)
    line1.set_ydata(u)
    plt.draw(); plt.pause(0.1)

plt.show()

   