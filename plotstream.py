import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
with open("/Users/bunny/Documents/msci/mscigit/streamdata/potentialfunction.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("/Users/bunny/Documents/msci/mscigit/streamdata/streamfunction.dat") as file_name:
    darray = np.loadtxt(file_name, delimiter=",")
print(np.shape(array))
grid = 10
timesteps = 5

points = np.zeros((timesteps,grid,grid))
points = array.reshape((timesteps,grid,grid))

points2 = darray.reshape((timesteps,grid,grid))


x_coords = np.linspace(-5,5,grid)
y_coords = np.linspace(5,5,grid)

print(points[1][3][1:10])
print(points[1][3][1:10])

fig, ax = plt.subplots(nrows=2,ncols=1)

c = ax[0].imshow(points[2],cmap='jet',extent=[-5,5,-5,5])

d = ax[1].imshow(points2[1],cmap='jet',extent=[-5,5,-5,5])

plt.show()

def init_func():
    plt.cla()

fig,axes = plt.subplots(1,2,figsize=(12,12))



def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    axes[0].imshow(points[ii,:,:])
    axes[1].imshow(points2[ii,:,:])

anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, timesteps),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=1)
anim.save('stream.mp4', writer=writervideo)
