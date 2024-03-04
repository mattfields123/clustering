import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter
import pickle
N_part = 100
tsteps = 500

with open("partx.dat") as file_name:
    partx = np.loadtxt(file_name)

with open("party.dat") as file_name:
    party = np.loadtxt(file_name)

with open("fpoint_x.dat",'rb') as file_name:
    fpoint_x = pickle.load(file_name)
with open("fpoint_y.dat",'rb') as file_name:
    fpoint_y = pickle.load(file_name)








partx = np.fmod(partx+1000,10)-5
party = np.fmod(party+1000,10)-5


print(partx.shape)

print(party.shape)

fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.scatter(400*partx[ii+1,:],400*party[ii+1,:],c='black',s=0.1)
    plt.scatter(400*fpoint_x[ii],400*fpoint_y[ii],c='red')
    plt.xlim(-2000, 2000)
    plt.ylim(-2000, 2000)
  

anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, 500),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=20)
anim.save('anim.mp4', writer=writervideo)
