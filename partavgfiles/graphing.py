from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma/partavg.dat") as file_name:
    array = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma smaller range/partavg.dat") as file_name:
    barray = np.loadtxt(file_name)

with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma2/partavg.dat") as file_name:
    array2 = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma smaller range2/partavg.dat") as file_name:
    barray2 = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma3/partavg.dat") as file_name:
    array3 = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/PARTICLE AVERAGES/partavggamma smaller range3/partavg.dat") as file_name:
    barray3 = np.loadtxt(file_name)





print(np.shape(array))
array = array*400
barray = barray*400
array2 = array2*400
barray2 = barray2*400
array3 = array3*400
barray3 = barray3*400


t = np.linspace(0,4000*0.6,4001)
arr = np.zeros(11)
barr = np.zeros(11)
arr2 = np.zeros(11)
barr2 = np.zeros(11)
arr3 = np.zeros(11)
barr3 = np.zeros(11)

for x in range(11):
    y = linregress(t[1000:4000],array[x,1000:4000])
    ylarge = linregress(t[1000:4000],barray[x,1000:4000])
    y2 = linregress(t[1000:4000],array2[x,1000:4000])
    ylarge2 = linregress(t[1000:4000],barray2[x,1000:4000])
    y3 = linregress(t[1000:4000],array3[x,1000:4000])
    ylarge3 = linregress(t[1000:4000],barray3[x,1000:4000])

    arr[x] = y.slope
    barr[x] = ylarge.slope
    arr2[x] = y2.slope
    barr2[x] = ylarge2.slope
    arr3[x] = y3.slope
    barr3[x] = ylarge3.slope

t0b = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] # gamma large
# t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64] # amp
t0 = [0.,0.005,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1] # gamma
# t0 = [0,0.005,0.01,0.05,0.075,0.1,0.125,0.25,0.5,0.75,1] #small delta
#t0 = [10*np.pi,20*np.pi,50*np.pi,100*np.pi,200*np.pi,500*np.pi,1000*np.pi,2000*np.pi,2*np.pi,np.pi,0.2*np.pi] # delta

#t0 = [0.2*np.pi,np.pi,2*np.pi,10*np.pi,20*np.pi,50*np.pi,100*np.pi,200*np.pi,500*np.pi,1000*np.pi,2000*np.pi]

# arrb = np.copy(arr)
# arr[0] = arrb[10]
# arr[1] = arrb[9]
# arr[2] = arrb[8]
# arr[3] = arrb[0]
# arr[4] = arrb[1]
# arr[5] = arrb[2]
# arr[6] = arrb[3]
# arr[7] = arrb[4]
# arr[8] = arrb[5]
# arr[9] = arrb[6]
# arr[10] = arrb[7]

# arrayb = np.copy(array)
# array[0] = arrayb[10]
# array[1] = arrayb[9]
# array[2] = arrayb[8]
# array[3] = arrayb[0]
# array[4] = arrayb[1]
# array[5] = arrayb[2]
# array[6] = arrayb[3]
# array[7] = arrayb[4]
# array[8] = arrayb[5]
# array[9] = arrayb[6]
# array[10] = arrayb[7]

t0[0] = 10e-10
t0b[0] = 10e-10
logt = np.log(t0)
logtb = np.log(t0b)
logx = np.log(-arr)
logb = np.log(-barr)

print(t0,t0b)

combined = np.zeros(20)
combinedt = np.zeros(20)

combinedt[0:10] = logt[1:11]
combined[0:10] = -arr[1:11]
combinedt[10:20] = logtb[1:11]
combined[10:20] = -barr[1:11]

print(logt)
print(logtb)

gg = linregress(logt[1:11],-arr[1:11])
print(gg)

plt.figure()
plt.scatter(logt[1:11],-arr[1:11],marker='x')
plt.scatter(logtb[1:11],-barr[1:11])
plt.scatter(logt[1:11],-arr2[1:11],marker='x')
plt.scatter(logtb[1:11],-barr2[1:11])



fig, ax = plt.subplots(1,2)
ax[0].set_xlabel('Gamma')
ax[0].set_ylabel('Drift Velocity / $\mathrm{kmday}^{-1}$')

#t1 = t0[0:3] + t0[4:11]
#arr2 = list(-arr[0:3]) + list(-arr[4:11])


#ax[0].plot(t0[1:3],-arr[1:3],)
#ax[0].plot(t0[4:11],-arr[4:11])

ax[0].scatter(t0,-arr,marker='x')
ax[0].scatter(t0b,-barr,marker='x')
ax[0].scatter(t0,-arr2)
ax[0].scatter(t0b,-barr2)
ax[0].scatter(t0,-arr3,marker='^')
ax[0].scatter(t0b,-barr3,marker='^')


#t0 = [10*np.pi,20*np.pi,50*np.pi,100*np.pi,200*np.pi,500*np.pi,1000*np.pi,2000*np.pi,2*np.pi,np.pi,0.2*np.pi] # delta


ax[1].set_xlabel('Westward Drift / km')
ax[1].set_ylabel('Time / day')
for x in range(11):
    ax[1].plot(array[x,:],t,label=str(t0b[x]))
    ax[1].plot(barray[x,:],t,label=str(t0[x]))

ax[1].legend()
plt.show()
