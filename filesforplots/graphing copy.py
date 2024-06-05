from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress



with open("/Users/bunny/Documents/msci/mscigit/filesforplots/gammasmallvalues.dat") as file_name:
    array = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/mscigit/filesforplots/gammalarge.dat") as file_name:
    barray = np.loadtxt(file_name)

# with open("/Users/bunny/Documents/msci/mscigit/amplitudevary/gamma/small.dat") as file_name:
#     array = np.loadtxt(file_name)
# with open("/Users/bunny/Documents/msci/mscigit/amplitudevary/gamma/large.dat") as file_name:
#     barray = np.loadtxt(file_name)

# with open("/Users/bunny/Documents/msci/mscigit/filesforplots/gammavarydelta1small.dat") as file_name:
#     array = np.loadtxt(file_name)
# with open("/Users/bunny/Documents/msci/mscigit/filesforplots/delta1gammavarylarge.dat") as file_name:
#     barray = np.loadtxt(file_name)



print(np.shape(barray))
array = array * 400 
barray = barray * 400


t = np.linspace(0,1000*0.25*(4*10**6)/(3600*24),1001)
arr = np.zeros(11)
barr = np.zeros(11)

for x in range(11):
    y = linregress(t[500:1000],array[x,500:1000])
    y1 = linregress(t[500:1000],barray[x,500:1000])
    arr[x] = y.slope
    barr[x] = y1.slope
t0b = [0.2,0.25,0.3,0.35,0.45,0.5,0.6,0.7,0.8,0.9,1] # gamma large

t0 = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15] # gamma small

logt = np.log10(t0)
logx = np.log10(-arr)
logt0b = np.log10(t0b)
logxb = np.log10(-barr)
print(t0,t0b)

combinedt = np.zeros(22)
combinedx = np.zeros(22)

combinedt[0:11] = logt
combinedx[0:11] = logx
combinedt[11:22] = logt0b
combinedx[11:22] = logxb

combined = np.zeros((2,22))
combined[0] = combinedt
combined[1] = combinedx

smallgamma = []
largegamma = []
smallx = []
largex = []


for r in range(22):
    if  -1.6 < combined[0][r] < -1:
        smallgamma.append(combined[0,r])
        smallx.append(combined[1,r])
    if combined[0][r] > -1:
        largegamma.append(combined[0,r])
        largex.append(combined[1,r])
print(smallgamma)
print(largex,'largex')

print(largegamma)
smalllin = linregress(smallgamma,smallx)
largelin = linregress(largegamma,largex)


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


fig, ax = plt.subplots(1,2)
ax[0].set_xlabel('Gamma')
ax[0].set_ylabel('Drift Velocity / $\mathrm{kmday}^{-1}$')

ax[0].scatter(t0,-arr,marker='x')
ax[0].scatter(t0b,-barr,marker='x')

ax[1].set_xlabel('Westward Drift / km')
ax[1].set_ylabel('Time / years')
for x in range(11):
    ax[1].plot(array[x,:],t/365,label=str(t0b[x]))
ax[1].legend()
plt.show()

