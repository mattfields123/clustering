import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

N_part = 500
tsteps = 50

with open("partx.dat") as file_name:
    partx = np.loadtxt(file_name)

with open("party.dat") as file_name:
    party = np.loadtxt(file_name)

with open("fixedpoint.dat") as file_name:
    fixedpoint = np.loadtxt(file_name)
with open("lengthvelocity.dat") as file_name:
    lengthvelocity = np.loadtxt(file_name)

print(fixedpoint.shape)



print(lengthvelocity)
lenvel = np.zeros(tsteps+2)
lenvel[1:tsteps+2] = lengthvelocity

partx.reshape(tsteps+1,N_part**2)
party.reshape(tsteps+1,N_part**2)

partx = np.fmod(partx+1000,10)-5
party = np.fmod(party+1000,10)-5


fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.scatter(400*partx[ii,:],400*party[ii,:],c='black',s=0.1)
    plt.scatter(400*fixedpoint[int(lengthvelocity[ii]):int(lengthvelocity[ii+1]),0],400*fixedpoint[int(lengthvelocity[ii]):int(lengthvelocity[ii+1]),1])
    plt.xlim(-2000, 2000)
    plt.ylim(-2000, 2000)


anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, 50),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=20)
anim.save('anim.mp4', writer=writervideo)
