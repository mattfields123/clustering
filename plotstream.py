import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
with open("/Users/bunny/Documents/msci/mscigit/streamdata/potentialfunction.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("/Users/bunny/Documents/msci/mscigit/streamdata/streamfunction.dat") as file_name:
    darray = np.loadtxt(file_name, delimiter=",")
print(np.shape(array))
points = np.zeros((10,10,10))
points = array.reshape((10,10,10))

points2 = darray.reshape((10,10,10))


print(points[0])
x_coords = np.linspace(-5,5,10)
y_coords = np.linspace(5,5,10)

print(points[1][3][1:10])
print(points[2][3][1:10])

fig, ax = plt.subplots(nrows=2,ncols=1)

c = ax[0].imshow(points[0],cmap='jet',extent=[-5,5,-5,5])

d = ax[1].imshow(points2[0],cmap='jet',extent=[-5,5,-5,5])

plt.show()
