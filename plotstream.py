import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
with open("/Users/bunny/Documents/msci/mscigit/streamdata/potentialfunction.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("/Users/bunny/Documents/msci/mscigit/streamdata/streamfunction.dat") as file_name:
    darray = np.loadtxt(file_name, delimiter=",")
print(np.shape(array))
grid = 100
timesteps = 20

points = np.zeros((timesteps,grid,grid))
points = array.reshape((timesteps,grid,grid))

points2 = darray.reshape((timesteps,grid,grid))


print(points[0])
x_coords = np.linspace(-5,5,grid)
y_coords = np.linspace(5,5,grid)

print(points[1][3][1:10])
print(points[2][3][1:10])

fig, ax = plt.subplots(nrows=2,ncols=1)

c = ax[0].imshow(points[10],cmap='jet',extent=[-5,5,-5,5])

d = ax[1].imshow(points2[10],cmap='jet',extent=[-5,5,-5,5])

plt.show()
