from matplotlib import pyplot as plt
import numpy as np
from scipy import stats


with open("partavgs/partavg0.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")

b_array = np.zeros((7,np.shape(array)))
def opening_files(b_array):
    with open("partavgs/partavg0.dat") as file_name:
        array[0] = np.loadtxt(file_name, delimiter=",")

    with open("partavgs/partavg005.dat") as file_name:
        array[1] = np.loadtxt(file_name, delimiter=",")

    with open("partavgs/partavg01.dat") as file_name:
        array[2] = np.loadtxt(file_name, delimiter=",")

    with open("partavgs/partavg025.dat") as file_name:
        array[3] = np.loadtxt(file_name, delimiter=",")

    with open("partavgs/partavg05.dat") as file_name:
        array[4] = np.loadtxt(file_name, delimiter=",")

    with open("partavgs/partavg075.dat") as file_name:
        array[5] = np.loadtxt(file_name, delimiter=",")
    with open("partavgs/partavg1.dat") as file_name:
        array[0] = np.loadtxt(file_name, delimiter=",")
t_array = np.linspace(0,400*np.,np.shape(array)[0])
c_array = opening_files(b_array)

stats.linregress()



