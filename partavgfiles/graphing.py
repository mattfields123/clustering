from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/mscigit/partavgfiles/partavg.dat") as file_name:
    array = np.loadtxt(file_name)

print(np.shape(array))

t = np.linspace(0,1,101)
arr = np.zeros(10)
for x in range(10):
    y = linregress(t,array[x,:])
    arr[x] = y.slope



plt.figure()

plt.plot(arr)

plt.show()

# plt.figure()
# for x in range(10):
#     plt.plot(t,array[x,:])

# plt.show()
