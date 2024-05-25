from mpl import *
import numpy as np
import matplotlib.pyplot as plt


rusanov = np.loadtxt("rusanov.txt")
lxf = np.loadtxt("lxf.txt")
dflu = np.loadtxt("dflu.txt")
nt   = np.loadtxt("nt.txt")
muscl  = np.loadtxt("muscl.txt")
godunov = np.loadtxt("godunov.txt")
ex = np.loadtxt('exact.txt')
bex = np.loadtxt('buckleyexact.txt')

#--------------------------density--------------------

fig = plt.figure()
#plt.plot(dflu[:,0],dflu[:,2],label='Exact', c='k')
#plt.plot(lxf[:,0],lxf[:,2],label='Exact', c='k')
#plt.plot(lxf[:,0],lxf[:,1],'ro',fillstyle='none',label='LF')
#plt.plot(rusanov[:,0],rusanov[:,1],'r+',fillstyle='none',label='LLF')
#plt.plot(godunov[:,0],godunov[:,1],'b*',fillstyle='none',label='GD')
#plt.plot(nt[:,0], nt[:,1],'c+',fillstyle='none',label='NT')
#plt.plot(muscl[:,0],muscl[:,1],'go',fillstyle='none',label='MMOD')
plt.plot(lxf[:,0][::2],lxf[:,1][::2],'o',fillstyle='none',label='LF',c='b')
plt.plot(dflu[:,0],dflu[:,1],'-',fillstyle='none',label='DFLU',c='r')
plt.plot(nt[:,0][::2], nt[:,1][::2],'ko',fillstyle='none',label='NT')
#plt.plot(lxf[:,0][::2],lxf[:,1][::2],'--',fillstyle='none',label='NSLF', c='g')
#plt.plot(bex[:,0],bex[:,1],'-',fillstyle='none',label='Exact',c='k')
#plt.plot(ex[:,0],ex[:,1],'-',fillstyle='none',label='Exact',c='k')
#plt.xlim(-1,1)
plt.xlabel('x')
plt.ylabel('u')
#sfig.legend(loc=15)
plt.legend();
#plt.legend(bbox_to_anchor=(1,0), loc="lower right")
#plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
               # mode="expand", borderaxespad=0, ncol=2)
                
plt.grid(True, linestyle = '--')
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('u.pdf')

plt.xlim(0.2,0.75)
plt.savefig('uzoomed.pdf')