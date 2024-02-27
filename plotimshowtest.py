import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

N_part = 2
tsteps = 5
meshsize = 5


# with open('fixed.dat') as file_name:
#     velocities = np.loadtxt(file_name)

# with open("partx.dat") as file_name:
#     partx = np.loadtxt(file_name)

# with open("party.dat") as file_name:
#     party = np.loadtxt(file_name)

# with open("fixedpoint.dat") as file_name:
#     fixedpoint = np.loadtxt(file_name)
# with open("counters.dat") as file_name:
#     lengthvelocity = np.loadtxt(file_name)

velocities = np.zeros((tsteps,meshsize,meshsize))
velocities[0,3,2] = 1
velocities[1,2,1] = 2
velocities[2,1,4] = 3

partx = np.array([[1,2,3,4],[2,3,4,1],[2,1,3,4],[4,3,2,1],[2,3,1,2]])
party = np.array([[1,2,3,4],[2,3,4,1],[2,1,3,4],[4,3,2,1],[2,3,1,2]])




fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()
c = plt.imshow(velocities[0],zorder=1)
fig.colorbar(c)

def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    c = plt.imshow(velocities[ii].transpose(),origin='lower',zorder=1,extent=[-2000,2000,-2000,2000])
    a = plt.scatter(400*partx[ii+1,:],400*party[ii+1,:],c='red',s=10,zorder=2)
    b = plt.scatter(400*fixedpoint[ii,0:int(lengthvelocity[ii]),0],400*fixedpoint[ii,0:int(lengthvelocity[ii]),1],c='red',zorder=3)
    plt.xlim(-2000, 2000)
    plt.ylim(-2000, 2000)

anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, 4),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=20)
anim.save('anim.mp4', writer=writervideo)
