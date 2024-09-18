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


u1 = np.loadtxt("sol1.txt") # eps  0.01
u2 = np.loadtxt("sol2.txt") # eps  0.001
u3 = np.loadtxt("sol3.txt") # eps  0.0001
u4 = np.loadtxt("solep0.txt") # eps  0.0


plt.figure()
plt.plot(u1[:,0],u1[:,1],'-',label=r'$\epsilon =10^{-2}$',color='black')
plt.plot(u2[:,0],u2[:,1],'-',label=r'$\epsilon =10^{-3}$',color='blue')
plt.plot(u3[:,0],u3[:,1],'-',label=r'$\epsilon =10^{-4}$',color='green')
plt.plot(u4[:,0],u4[:,1],'--',label=r'$\epsilon= 0$',color='red')
plt.xlim()
plt.xlabel('x')
plt.ylabel('density')
plt.legend(fontsize=12); plt.grid(True, linestyle = '-', linewidth = 0.5)
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('density.pdf')

plt.figure()
plt.plot(u1[:,0],u1[:,2],'-',label=r'$\epsilon= 10^{-2}$',color='black')
plt.plot(u2[:,0],u2[:,2],'-',label=r'$\epsilon= 10^{-3}$',color='blue')
plt.plot(u3[:,0],u3[:,2],'-',label=r'$\epsilon= 10^{-4}$',color='green')
plt.plot(u4[:,0],u4[:,2],'--',label=r'$\epsilon= 0 $',color='red')
plt.xlim()
plt.xlabel('x')
plt.ylabel('velocity')
plt.legend(fontsize=12); plt.grid(True, linestyle = '-', linewidth = 0.5)
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('velocity.pdf')

plt.figure()
plt.plot(u1[:,0],u1[:,3],'-',label=r'$\epsilon= 10^{-2}$',color='black')
plt.plot(u2[:,0],u2[:,3],'-',label=r'$\epsilon= 10^{-3}$',color='blue')
plt.plot(u3[:,0],u3[:,3],'-',label=r'$\epsilon= 10^{-4}$',color='green')
plt.plot(u4[:,0],u4[:,3],'-',label=r'$\epsilon= 0 $',color='red')
plt.xlim()
plt.xlabel('x')
plt.ylabel('height')
plt.legend(loc='upper left',fontsize=12); plt.grid(True, linestyle = '-', linewidth = 0.5)
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('height.pdf')
