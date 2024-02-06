from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress


with open("/Users/bunny/Documents/msci/mscigit/files/beta32avg.dat") as file_name:
    array = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/mscigit/files/beta64avg.dat") as file_name:
    barray = np.loadtxt(file_name)







print(np.shape(array))
print(np.shape(barray))
array = array*400
barray = barray*400
# array2 = array2*400
# barray2 = barray2*400
# array3 = array3*400
t = np.linspace(0,10000*0.05*(4*10**5/(3600*24)),10001)
arr = np.zeros(11)
barr = np.zeros(11)
arr2 = np.zeros(11)
barr2 = np.zeros(11)
arr3 = np.zeros(11)

for x in range(11):
    y = linregress(t[1000:10000],array[x,1000:10000])
    ylarge = linregress(t[1000:10000],barray[x,1000:10000])
    # y2 = linregress(t[1000:4000],array2[x,1000:4000])
    # ylarge2 = linregress(t[1000:4000],barray2[x,1000:4000])
    # y3 = linregress(t[1000:4000],array3[x,1000:4000])

    arr[x] = y.slope
    barr[x] = ylarge.slope
    # arr2[x] = y2.slope
    # barr2[x] = ylarge2.slope
    # arr3[x] = y3.slope
# t0b = [0.01,0.05,0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] # gamma large
# t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64] # amp
t0 = [0.01,0.1,0.5,1,2,3,4,8,16,32,64] # gamma
t0b = [0.01,0.1,0.5,1,2,3,4,8,16,32,64] # gamma


t02 = [0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.5,0.75,0.9]
t0b2 = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
t0c = [0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]


logt = np.log10(t0)
logtb = np.log10(t0b)
logx = np.log10(-arr)
logb = np.log10(-barr)
# logx2 = np.log10(-arr2)
# logb2 = np.log10(-barr2)

# logt0c = np.log10(t0c)
# logx3 = np.log10(-arr3)

# logt02 = np.log10(t02)
# logt0b2 = np.log10(t0b2)


print(t0,t0b)

combinedt = np.zeros(5*11-3)
combinedx = np.zeros(5*11-3)

combinedt[0:10] = logt[1:11]
combinedt[10:20] = logtb[1:11]
# combinedt[20:31] = logt02
# combinedt[31:41] = logt0b2[1:11]
# combinedt[41:52] = logt0c

combinedx[0:10] = logx[1:11]
combinedx[10:20] = logb[1:11]
# combinedx[20:31] = logx2
# combinedx[31:41] = logb2[1:11]
# combinedx[41:52] = logx3

combined = np.zeros((2,52))
combined[0] = combinedt
combined[1] = combinedx

smallgamma = []
largegamma = []
smallx = []
largex = []


# for r in range(52):
#     if  -2 < combined[0][r] < -1:
#         smallgamma.append(combined[0,r])
#         smallx.append(combined[1,r])
#     if combined[0][r] > -1:
#         largegamma.append(combined[0,r])
#         largex.append(combined[1,r])
# print(smallgamma)

# smalllin = linregress(smallgamma,smallx)
# largelin = linregress(largegamma,largex)

# x_array = np.linspace(-2,-1,100)
# x_array_large = np.linspace(-1,0,100)
# fig4 = plt.figure()
# plt.scatter(smallgamma,smallx,marker='x')
# plt.scatter(largegamma,largex,marker='x')
# plt.plot(x_array,smalllin.slope*x_array+smalllin.intercept)
# plt.plot(x_array_large,largelin.slope*x_array_large+largelin.intercept)
# # plt.show()

# print(smalllin)
# print(largelin)




fig3, ax3 = plt.subplots(2,2)
ax3[1,1].scatter(logt,logx,marker='x',label='$\gamma=1$')
ax3[1,1].scatter(logtb,logb,marker='x')

# ax3[0,1].scatter(t0,logx,marker='x',label='gamma=1')
# ax3[0,1].scatter(t0b,logb,marker='x')

# ax3[1,0].scatter(logt,-arr,marker='x',label='gamma=1')
# ax3[1,0].scatter(logtb,-barr,marker='x')



# ax3[1,1].scatter(logt02[0:11],logx2[0:11],marker='x')
# ax3[1,1].scatter(logt0b2[1:11],logb2[1:11],marker='x')
# ax3[1,1].scatter(logt0c[0:10], logx3[0:10],marker='x')


ax3[1,1].set_xlabel('log $Amplitude Scaling$')
ax3[1,1].set_ylabel('log C')

# ax3[1,0].scatter(t0[1:11],logx[1:11],marker='x')
# ax3[1,0].scatter(t0b[1:11],logb[1:11],marker='x')
# ax3[1,0].set_xlabel('$\gamma$')
# ax3[1,0].set_ylabel('log C')

# ax3[0,1].scatter(logt[1:11],-arr[1:11],marker='x')
# ax3[0,1].scatter(logtb[1:11],-barr[1:11],marker='x')
# ax3[0,1].set_xlabel('log $\gamma$')
# ax3[0,1].set_ylabel('C')

ax3[0,0].scatter(t0[1:11],-arr[1:11],marker='x',label='$\gamma=1$')
ax3[0,0].scatter(t0b[1:11],-barr[1:11],marker='x')
# ax3[0,0].scatter(t02[1:11],-arr2[1:11],marker='x')
# ax3[0,0].scatter(t0b2[1:11],-barr2[1:11],marker='x')


ax3[0,0].set_xlabel('$Amplitude Scaling$')
ax3[0,0].set_ylabel('C')
plt.legend()

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
