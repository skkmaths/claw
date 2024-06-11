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
plt.grid(True); plt.draw(); plt.pause(0.1)

Tf = 0.4
deltat = 0.01
time = 0 

while time < Tf:
    time += deltat
    tt = time
    for i in range(0, len(x)):
        guess = smooth(x[i])
        xx = x[i]
        u[i]  = fsolve(impf, guess)
    line1.set_ydata(u)
    plt.draw(); plt.pause(0.1)

plt.show()

   