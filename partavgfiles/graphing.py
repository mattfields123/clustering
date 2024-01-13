from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/mscigit/partavgfiles/partavg.dat") as file_name:
    array = np.loadtxt(file_name)

print(np.shape(array))
array = array*400
t = np.linspace(0,4000*0.6,4001)
arr = np.zeros(11)
for x in range(11):
    y = linregress(t,array[x,:])
    arr[x] = y.slope

t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64]

fig1 = plt.figure()
plt.xlabel('Amplitude scaling ($alpha$)')
plt.ylabel('Drift Velocity / $\mathrm{kmday}^{-1}$')
plt.plot(t0,arr,)

plt.show()

fig2 = plt.figure()
plt.xlabel('Westward Drift / km')
plt.ylabel('Time / day')
for x in range(11):
    plt.plot(array[x,:],t,label=str(t0[x]))

plt.legend()
plt.show()
