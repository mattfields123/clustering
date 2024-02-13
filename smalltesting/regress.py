from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/mscigit/smalltesting/partavg.dat") as file_name:
    array = np.loadtxt(file_name)





t = np.linspace(0,10000*0.01,10001)

print(np.shape(array))
# array = array*400

y = linregress(t[1000:10000],array[1000:10000])
arr = y.slope
print(y)