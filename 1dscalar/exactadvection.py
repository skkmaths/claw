from xml import dom
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from functools import partial
from scipy import optimize

def smooth(x):
    return np.sin(np.pi*x)

def impf(u):
    return u- smooth(xx - tt*u)

x  = np.linspace(-1,1,200)
u  = smooth(x)

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, u)
plt.grid(True); plt.draw(); plt.pause(1)

Tf = 2.0;
deltat = 0.01
t = 0 

while t < Tf:
    t
    for i in range(0, len(x)):
        u[i]  = smooth(x[i]-t)
    line1.set_ydata(u)
    plt.draw(); plt.pause(0.1)
    t += deltat

plt.show()

   