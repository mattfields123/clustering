from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/mscigit/partavgfiles/g1d1vinfa1rd100.dat") as file_name:
    array = np.loadtxt(file_name)





t = np.linspace(0,1000*0.25,1001)

print(np.shape(array))
print(array)
# array = array*400

y = linregress(t[500:1000],array[500:1000])

print(y)


