from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress



with open("/Users/bunny/Documents/msci/mscigit/filesforplots/amplitudes.dat") as file_name:
    array = np.loadtxt(file_name)






# print(np.shape(barray))
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
 #   y1 = linregress(t[500:1000],barray[x,500:1000])
    arr[x] = y.slope
#    barr[x] = y1.slope
#t0b = [0.2,0.25,0.3,0.35,0.45,0.5,0.6,0.7,0.8,0.9,1] # gamma large
# t0 = [0.1,0.5,1,2,3,4,9,25,36,49,64] # amp
t0 = [0.01,0.05,0.1,0.25,0.5,1,2,4,8,32,64] # gamma


logt = np.log10(t0)
logx = np.log10(-arr)
#logt0b = np.log10(t0b)
#logxb = np.log10(-barr)
#print(t0,t0b)

combinedt = np.zeros(6)
combinedx = np.zeros(6)

combinedt[0:9] = logt[3:9]
combinedx[0:9] = logx[3:9]
# combinedt[11:22] = logtb
# combinedx[11:22] = logxb

combined = np.zeros((2,6))
combined[0] = combinedt
combined[1] = combinedx
print(combined,'combined')
smalllin = linregress(combinedt,combinedx)

x_array = np.linspace(-2,-1,100)
x_array_large = np.linspace(-1,0,100)
fig4 = plt.figure()
plt.scatter(combinedt,combinedx,marker='x')
plt.plot(x_array,smalllin.slope*x_array+smalllin.intercept,label=str(smalllin.slope))

plt.xlabel('$\log\gamma$')
plt.ylabel('$\log C$')
plt.legend()
plt.show()

print(smalllin.slope)
print(smalllin,'smalllin')




fig3, ax3 = plt.subplots(2,2)
ax3[1,1].scatter(logt[0:9],logx[0:9],marker='x')



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

ax3[0,0].plot(t0[0:9],-arr[0:9],marker='x')


ax3[0,0].set_xlabel('alpha')
ax3[0,0].set_ylabel('C')

plt.show()





fig, ax = plt.subplots(1,2)
ax[0].set_xlabel('Amplitude scaling')
ax[0].set_ylabel('Drift Velocity / $\mathrm{kmday}^{-1}$')

#t1 = t0[0:3] + t0[4:11]
#arr2 = list(-arr[0:3]) + list(-arr[4:11])


#ax[0].plot(t0[1:3],-arr[1:3],)
#ax[0].plot(t0[4:11],-arr[4:11])

ax[0].plot(t0,-arr,marker='x')


#t0 = [10*np.pi,20*np.pi,50*np.pi,100*np.pi,200*np.pi,500*np.pi,1000*np.pi,2000*np.pi,2*np.pi,np.pi,0.2*np.pi] # delta


ax[1].set_xlabel('Westward Drift / km')
ax[1].set_ylabel('Time / day')

ax[1].legend()
plt.show()

print(arr)
