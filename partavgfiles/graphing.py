from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/mscigit/partavgfiles/partavg.dat") as file_name:
    array = np.loadtxt(file_name)

print(np.shape(array))

t = np.linspace(0,1,4001)
arr = np.zeros(11)
for x in range(11):
    y = linregress(t,array[x,:])
    arr[x] = y.slope

t0 = [0,0.005,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1]

plt.figure()
plt.xlabel('Gamma')
plt.ylabel('Drift Velocity')
plt.plot(t0,arr,)

plt.show()

plt.figure()
plt.xlabel('Westward Drift')
plt.ylabel('Time')
for x in range(11):
    plt.plot(array[x,:],t)
plt.show()
