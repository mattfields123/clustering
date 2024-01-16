from matplotlib import pyplot as plt
import numpy as np
with open("001.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("01.dat") as file_name:
    darray = np.loadtxt(file_name, delimiter=",")
with open("05.dat") as file_name:
    carray = np.loadtxt(file_name, delimiter=",")

array001 = array[:,0] * 400/5
array01 = darray[:,0] * 400/5
array05 = carray[:,0] * 400/5
print(np.shape(array05))

t = np.linspace(0,10000,10001)

fig = plt.figure()
plt.plot(array001,t,color='green')
plt.plot(array01,t,color='red')
plt.plot(array05,t,color='orange')
plt.savefig('lawrence.png')