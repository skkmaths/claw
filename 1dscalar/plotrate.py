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
         row = [data[i,0], data[i,1], rate1, data[i,2], rate2]
      else:
         rate1, rate2 = ( np.log(data[i-1,1]/data[i,1])/np.log(2.0),
                          np.log(data[i-1,2]/data[i,2])/np.log(2.0) )
         row = [data[i,0], data[i,1], rate1, data[i,2], rate2]
      print (s.format(*row))

model1 = np.loadtxt("error.txt")
convergence_rate(model1)