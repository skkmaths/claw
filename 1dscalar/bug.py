import numpy as np

dummy = 1.0

def compute_a():
    a = dummy**2 
    return a

def compute_b(dummy1):
    a = dummy1**2
    return a

def check(dummy):
    dummy = 3.0
    print( 'value a =',compute_a(), 'value b =', compute_b(dummy))
t = 0
while t <10:
    
    dummy = 2.0
    check(dummy)
    
    t +=1


