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

# t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64] # amp
# t0 = [0.,0.005,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1] # gamma
# t0 = [0,0.005,0.01,0.05,0.075,0.1,0.125,0.25,0.5,0.75,1] #small delta
t0 = [10*np.pi,20*np.pi,50*np.pi,100*np.pi,200*np.pi,500*np.pi,1000*np.pi,2000*np.pi,2*np.pi,np.pi,0.2*np.pi] # delta
fig, ax = plt.subplots(1,2)
ax[0].set_xlabel('Delta')
ax[0].set_ylabel('Drift Velocity / $\mathrm{kmday}^{-1}$')
ax[0].plot(t0,-arr,)

ax[1].set_xlabel('Westward Drift / km')
ax[1].set_ylabel('Time / day')
for x in range(11):
    ax[1].plot(array[x,:],t,label=str(t0[x]))

ax[1].legend()
plt.show()
