#------------------------------------------------------------------------------
# Requirements
# pip install latextable
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Formatting string with 5-digit precision for errors and EOC
s = "{:<10} {:<15} {:<10} {:<15} {:<10} {:<15} {:<10}"

def convergence_rate(data):
    N1 = len(data[:, 0])
    print("____________________________________________________________________________________________________")
    print(s.format("h", "L1 error", "order", "L2 error", "order", "Linf error", "order"))
    print("____________________________________________________________________________________________________")
    for i in range(N1):
        h = f"{data[i, 0]:.5f}"
        L1 = f"{data[i, 2]:.5e}"
        L2 = f"{data[i, 3]:.5e}"
        Linf = f"{data[i, 4]:.5e}"
        
        if i == 0:
            rate1 = rate2 = rate3 = '-'
        else:
            rate1 = f"{np.log(data[i - 1, 2] / data[i, 2]) / np.log(2.0):.5f}"
            rate2 = f"{np.log(data[i - 1, 3] / data[i, 3]) / np.log(2.0):.5f}"
            rate3 = f"{np.log(data[i - 1, 4] / data[i, 4]) / np.log(2.0):.5f}"
        
        print(s.format(h, L1, rate1, L2, rate2, Linf, rate3))

# Assuming the columns are:
# [h, <optional>, L1 error, L2 error, Linf error]
model1 = np.loadtxt("error.txt")
convergence_rate(model1)