import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

N_part = 100
tsteps = 500

with open("partx.dat") as file_name:
    partx = np.loadtxt(file_name)

with open("party.dat") as file_name:
    party = np.loadtxt(file_name)

with open("fixedpoint.dat") as file_name:
    fixedpoint = np.loadtxt(file_name)
with open("counters.dat") as file_name:
    lengthvelocity = np.loadtxt(file_name)


print(fixedpoint.shape)

print(partx.shape)

print(lengthvelocity)



fixedpoint = fixedpoint.reshape(tsteps,5000,2,order='F')

#partx.reshape(tsteps+1,N_part**2)
#party.reshape(tsteps+1,N_part**2)

print(fixedpoint[1,int(lengthvelocity[1]),1])


partx = np.fmod(partx+1000,10)-5
party = np.fmod(party+1000,10)-5
fixedpoint = np.fmod(fixedpoint+1000,10)-5


print(partx.shape)

print(party.shape)

print(fixedpoint.shape)


print(fixedpoint[1,int(lengthvelocity[1]),1])

fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.scatter(400*partx[ii+1,:],400*party[ii+1,:],c='black',s=0.1)
    plt.scatter(400*fixedpoint[ii,0:int(lengthvelocity[ii]),0],400*fixedpoint[ii,0:int(lengthvelocity[ii]),1],c='red')
    plt.xlim(-2000, 2000)
    plt.ylim(-2000, 2000)
    print(fixedpoint[ii,0:int(lengthvelocity[ii]),0])

anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, 500),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=20)
anim.save('anim.mp4', writer=writervideo)
