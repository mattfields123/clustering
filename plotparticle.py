import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

N_part = 200
tsteps = 1000

with open("partx.dat") as file_name:
    partx = np.loadtxt(file_name)

with open("party.dat") as file_name:
    party = np.loadtxt(file_name)

partx.reshape(tsteps+1,N_part**2)
party.reshape(tsteps+1,N_part**2)

partx = np.fmod(partx,10)-5
party = np.fmod(party,10)-5


fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.scatter(partx[ii,:],party[ii,:])
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)


anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, tsteps+1),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=10)
anim.save('inertial.mp4', writer=writervideo)
