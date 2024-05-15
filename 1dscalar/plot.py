from mpl import *
import numpy as np
import matplotlib.pyplot as plt



lxf = np.loadtxt("lxf.txt")
dflu = np.loadtxt("dflu.txt")
nt   = np.loadtxt("nt.txt")


#ue  = np.loadtxt("exact.txt")


#--------------------------density--------------------

fig = plt.figure()
#plt.plot(ue[:,0],ue[:,1],label='Exact', c='k')
#plt.plot([:,0],initial[:,1],'-',fillstyle='none',label='I.C.')
plt.plot(lxf[:,0],lxf[:,1],'-',fillstyle='none',label='lxf',c='b')
plt.plot(dflu[:,0],dflu[:,1],'--',fillstyle='none',label='dflu',c='r')
plt.plot(nt[:,0], nt[:,1],'--',fillstyle='none',label='NT',c='c')
#plt.plot(fv[:,0],fv[:,1],'--',fillstyle='none',label='First-order', c='g')
#plt.plot(ex[:,0],ex[:,1],'-',fillstyle='none',label='Ref.',c='k')
plt.xlim(-4,4)
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