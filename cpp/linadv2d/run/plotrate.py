#------------------------------------------------------------------------------
# Requirements
# pip install latextable
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

s = "{:<15} {:<25} {:<25} {:<25}  {:<25}"

def convergence_rate(data):
   N1 = len(data[:,0])
   print("___________________________________________________________________________________________________________________")
   print(s.format("h","L1 error", "order","L2 errror", "order"))
   print("___________________________________________________________________________________________________________________")
   for i in range(N1):
      if i == 0:
         rate1, rate2 = '-', '-'
         row = [data[i,0], data[i,4], rate1, data[i,5], rate2]
      else:
         rate1, rate2 = ( np.log(data[i-1,4]/data[i,4])/np.log(2.0),
                          np.log(data[i-1,5]/data[i,5])/np.log(2.0) )
         row = [data[i,0], data[i,4], rate1, data[i,5], rate2]
      print (s.format(*row))

model1 = np.loadtxt("err.dat")
convergence_rate(model1)