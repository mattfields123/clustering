import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# with open('/Users/bunny/Documents/msci/mscigit/veltest/fixed.dat') as file_name:
#   velocities = np.loadtxt(file_name)

with open('fixed.dat') as file_name:
    velocities = np.loadtxt(file_name)




tsteps = 500
meshsize = 1000




velocities = velocities.reshape(tsteps,meshsize,meshsize)






        


print(velocities.shape)



fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()
c = plt.imshow(velocities[0])
fig.colorbar(c)



def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    c = plt.imshow(velocities[ii].transpose(),origin='lower',extent=[-5,5,-5,5])
    


anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, tsteps),
                     init_func=init_func)





writervideo = FFMpegWriter(fps=20)
anim.save('fixedxprofile.mp4', writer=writervideo)

