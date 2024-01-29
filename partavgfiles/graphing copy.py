from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma small range4/partavg.dat") as file_name:
    array = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma 4/partavg.dat") as file_name:
    barray = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma small range5/partavg.dat") as file_name:
    array2 = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma 5/partavg.dat") as file_name:
    barray2 = np.loadtxt(file_name)






print(np.shape(array))
array = array*400
barray = barray*400
array2 = array2*400
barray2 = barray2*400


t = np.linspace(0,4000*0.6,4001)
arr = np.zeros(11)
barr = np.zeros(11)
arr2 = np.zeros(11)
barr2 = np.zeros(11)


for x in range(11):
    y = linregress(t[1000:4000],array[x,1000:4000])
    ylarge = linregress(t[1000:4000],barray[x,1000:4000])
    y2 = linregress(t[1000:4000],array2[x,1000:4000])
    ylarge2 = linregress(t[1000:4000],barray2[x,1000:4000])


    arr[x] = y.slope
    barr[x] = ylarge.slope
    arr2[x] = y2.slope
    barr2[x] = ylarge2.slope
t0b = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] # gamma large
# t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64] # amp
t0 = [0.,0.005,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1] # gamma
t02 = [0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.5,0.75,0.9]
t0b2 = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]

t0[0] = 10e-10
t0b[0] = 10e-10
t0b2[0] = 10e-10

logt = np.log10(t0)
logtb = np.log10(t0b)
logx = np.log10(-arr)
logb = np.log10(-barr)
logx2 = np.log10(-arr2)
logb2 = np.log10(-barr2)

logt02 = np.log10(t02)
logt0b2 = np.log10(t0b2)


print(t0,t0b)

combined = np.zeros(20)
combinedt = np.zeros(20)

combinedt[0:10] = logt[1:11]
combined[0:10] = logx[1:11]
combinedt[10:20] = logtb[1:11]
combined[10:20] = logb[1:11]

print(logt)
print(logtb)

gg = linregress(combinedt,combined)
print(gg)

fig3, ax3 = plt.subplots(2,2)
ax3[1,1].scatter(logt[1:11],logx[1:11],marker='x')
ax3[1,1].scatter(logtb[1:11],logb[1:11],marker='x')
ax3[1,1].scatter(logt02[1:11],logx2[1:11],marker='x')
ax3[1,1].scatter(logt0b2[1:11],logb2[1:11],marker='x')

ax3[1,1].set_xlabel('log $\gamma$')
ax3[1,1].set_ylabel('log C')

# ax3[1,0].scatter(t0[1:11],logx[1:11],marker='x')
# ax3[1,0].scatter(t0b[1:11],logb[1:11],marker='x')
# ax3[1,0].set_xlabel('$\gamma$')
# ax3[1,0].set_ylabel('log C')

# ax3[0,1].scatter(logt[1:11],-arr[1:11],marker='x')
# ax3[0,1].scatter(logtb[1:11],-barr[1:11],marker='x')
# ax3[0,1].set_xlabel('log $\gamma$')
# ax3[0,1].set_ylabel('C')

ax3[0,0].scatter(t0[1:11],-arr[1:11],marker='x')
ax3[0,0].scatter(t0b[1:11],-barr[1:11],marker='x')
ax3[0,0].scatter(t02[1:11],-arr2[1:11],marker='x')
ax3[0,0].scatter(t0b2[1:11],-barr2[1:11],marker='x')


ax3[0,0].set_xlabel('$\gamma$')
ax3[0,0].set_ylabel('C')

plt.show()





fig, ax = plt.subplots(1,2)
ax[0].set_xlabel('Gamma')
ax[0].set_ylabel('Drift Velocity / $\mathrm{kmday}^{-1}$')

#t1 = t0[0:3] + t0[4:11]
#arr2 = list(-arr[0:3]) + list(-arr[4:11])


#ax[0].plot(t0[1:3],-arr[1:3],)
#ax[0].plot(t0[4:11],-arr[4:11])

ax[0].scatter(t0,-arr,marker='x')
ax[0].scatter(t0b,-barr,marker='x')



#t0 = [10*np.pi,20*np.pi,50*np.pi,100*np.pi,200*np.pi,500*np.pi,1000*np.pi,2000*np.pi,2*np.pi,np.pi,0.2*np.pi] # delta


ax[1].set_xlabel('Westward Drift / km')
ax[1].set_ylabel('Time / day')
for x in range(11):
    ax[1].plot(array[x,:],t,label=str(t0b[x]))
    ax[1].plot(barray[x,:],t,label=str(t0[x]))

ax[1].legend()
plt.show()
