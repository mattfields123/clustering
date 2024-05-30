import numpy as np
from matplotlib import pyplot as plt

with open("/Users/bunny/Documents/msci/mscigit/streamdata/str.dat") as file_name:
    stream = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/mscigit/streamdata/p.dat") as file_name:
    potential = np.loadtxt(file_name)

print(np.shape(stream))
print(np.shape(potential))

fig, ax = plt.subplots(2)


c = ax[0].imshow(stream[3].reshape(250,250),extent=[-5,5,-5,5],cmap='rainbow')
d = ax[1].imshow(potential[3].reshape(250,250),extent = [-5,5,-5,5],cmap='rainbow')
plt.colorbar(c)
plt.colorbar(d)
plt.show()
print(250*5*250)