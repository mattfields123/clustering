import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
with open("/Users/bunny/Documents/msci/mscigit/streamdata/potentialfunction.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("/Users/bunny/Documents/msci/mscigit/streamdata/streamfunction.dat") as file_name:
    darray = np.loadtxt(file_name, delimiter=",")
print(np.shape(array))
grid = 100
timesteps = 50

points = np.zeros((timesteps,grid,grid))
points2 = np.zeros((timesteps,grid,grid))
for t in range(timesteps):
    for g in range(grid):
        for h in range(grid):
            points[t,g,h] = array[t*timesteps+g*grid+h]
            points2[t,g,h] = darray[t*timesteps+g*grid+h]
            

plt.imshow(darray[0:10000].reshape((100,100)))
for x in range(timesteps):
    print(darray[10000*x],darray[10000*x+1],darray[10000*x+5])


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

writervideo = FFMpegWriter(fps=5)
anim.save('stream.mp4', writer=writervideo)
