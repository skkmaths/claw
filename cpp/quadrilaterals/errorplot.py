import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 10
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.autoscale(enable=True, axis='x', tight=True)
#plt.rcParams['text.usetex'] = True    # This will use Latex fonts (slow)
plt.axis('tight')
plt.autoscale(enable=True, axis='x', tight=True)


fo = np.loadtxt("error.txt")

plt.figure()
plt.plot(fo[:,0],fo[:,1],'--',label='FO', c='r', marker = 'o')

# Add a reference line with slope -1, adjusted to be close to the data
x_vals = np.array([min(fo[:,0]), max(fo[:,0])])
y_vals = fo[0, 1] * (x_vals / fo[0, 0])**(-1)  # y = k * x^(-1) where k is chosen to match the first data point
plt.plot(x_vals, y_vals, label='Slope -1', color='b', linestyle='-.')

plt.xscale("log")
plt.yscale("log")
#plt.xlim(0,1)
plt.xlabel('Number of cells')
plt.ylabel('$||u_h-u_s||_1$')
plt.legend(fontsize=12); plt.grid(True, linestyle = '--', linewidth = 0.5)
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('eoc.pdf')