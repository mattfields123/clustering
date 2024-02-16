import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

with open('/Users/bunny/Documents/msci/mscigit/veltest/velocities.dat') as file_name:
    velocities = np.loadtxt(file_name)


tsteps = 100
meshsize = 10

velocities = velocities.reshape(tsteps,meshsize,meshsize)


print(velocities.shape)



fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.imshow(velocities[ii])


anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, 100),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=20)
anim.save('mscigit/velocityxprofile.mp4', writer=writervideo)

