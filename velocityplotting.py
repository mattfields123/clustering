import numpy as np
from matplotlib import pyplot as plt


with open('velocities.dat') as file_name:
    velocities = np.loadtxt(file_name)


tsteps = 100
meshsize = 10

velocities.reshape((tsteps,meshsize))





