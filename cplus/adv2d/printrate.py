import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser()                                               

parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()

#use python plotrate.py -f filename.txt

s = "{:<15} {:<25} {:<25} "

def convergence_rate(data):
   N1 = len(data[:,0])
   print("___________________________________________________________________________________________________________________")
   print(s.format("dx","L1-error", "order"))
   print("___________________________________________________________________________________________________________________")
   for i in range(N1):
      if i == 0:
         rate1, rate2 = '-', '-'
         row = [data[i,1], data[i,2], rate1]
      else:
         rate1= ( np.log(data[i-1,2]/data[i,2])/np.log(2.0) )
         row = [data[i,1], data[i,2], rate1]
      print (s.format(*row))
    
model1 = np.loadtxt(args.file)
convergence_rate(model1)