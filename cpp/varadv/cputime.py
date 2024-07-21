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


fo = np.loadtxt("error_fo.txt")
so = np.loadtxt("error_so.txt")

plt.figure()
plt.plot(fo[:,3],fo[:,2],'--',label='FO', c='r', marker = 'o')
plt.plot(so[:,3],so[:,2],':',label='SO', c='k',marker = '+') # marker = 'o', fillstyle='none', markersize = 5)

plt.xscale("log")
plt.yscale("log")
#plt.xlim(0,1)
plt.xlabel('CPU time [s]')
plt.ylabel('$||u_h-u_s||_1$')
plt.legend(fontsize=12); plt.grid(True, linestyle = '--', linewidth = 0.5)
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('cputime.pdf')

