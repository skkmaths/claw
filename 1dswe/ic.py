import numpy as np

# Set intial condition for SWE equation


def density(x):
    if x < 0.0:
        return 0.2
    else: 
        return 0.5
density = np.vectorize(density)

def velocity(x):
    if x < 0.0:
        return 0.8
    else:
        return 0.4
velocity = np.vectorize(velocity)

def height(x):
    if x < 0.0:
        return 0.9
    else:
        return 0.2
height = np.vectorize(height)
    