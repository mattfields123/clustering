from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress



with open("/Users/bunny/Documents/msci/mscigit/filesforplots/gammasmallvalues.dat") as file_name:
    array = np.loadtxt(file_name)






print(np.shape(array))
# array = array*400
# barray = barray*400
# array2 = array2*400
# barray2 = barray2*400

# array3 = array3*400

t = np.linspace(0,1000*0.25,1001)
arr = np.zeros(11)
barr = np.zeros(11)
arr2 = np.zeros(11)
barr2 = np.zeros(11)
arr3 = np.zeros(11)

for x in range(11):
    y = linregress(t[500:1000],array[x,500:1000])
   
    arr[x] = y.slope

t0b = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] # gamma large
# t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64] # amp
t0 = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15] # gamma


logt = np.log10(t0)
logx = np.log10(-arr)

print(t0,t0b)
combinedt = logt
combinedx = logx
combinedt[0:10] = logt[1:11]

combined = np.zeros((2,11))
combined[0] = combinedt
combined[1] = combinedx

smallgamma = []
largegamma = []
smallx = []
largex = []


for r in range(11):
    if  -2 < combined[0][r] < -1:
        smallgamma.append(combined[0,r])
        smallx.append(combined[1,r])
    if combined[0][r] > -1:
        largegamma.append(combined[0,r])
        largex.append(combined[1,r])
print(smallgamma)

smalllin = linregress(smallgamma,smallx)
#largelin = linregress(largegamma,largex) add in when got more data
largelin = linregress([0,1],[1,2])


x_array = np.linspace(-2,-1,100)
x_array_large = np.linspace(-1,0,100)
fig4 = plt.figure()
plt.scatter(smallgamma,smallx,marker='x')
plt.scatter(largegamma,largex,marker='x')
plt.plot(x_array,smalllin.slope*x_array+smalllin.intercept,label=str(smalllin.slope))
plt.plot(x_array_large,largelin.slope*x_array_large+largelin.intercept,label=str(largelin.slope))
plt.xlabel('$\log\gamma$')
plt.ylabel('$\log C$')
plt.legend()
plt.show()

print(smalllin)
print(largelin)




fig3, ax3 = plt.subplots(2,2)
ax3[1,1].scatter(logt[1:11],logx[1:11],marker='x')



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
ax[1].legend()
plt.show()

print(arr)
