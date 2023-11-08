import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter
time = 100
N_part = 10


with open("velocity.dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("partavg.dat") as file_name:
    arrayP2 = np.loadtxt(file_name, delimiter=",")
print('Other files: tick!')
with open("particles.dat") as file_name:
    print('Opened')
    arrayP = np.loadtxt(file_name, delimiter=",", usecols=(0, 1), max_rows=time*N_part**2) # noqa E501
print('I am done reading')

array1 = np.zeros((128, 128, 2*time))
# array1 = np.zeros((128, 128, 2))
array2 = np.zeros((N_part**2, 2*time))
array3 = np.zeros(2*time)
# array1[:, :, 0] = array[0:4*4096, 0].reshape(128, 128)
# array1[:, :, 1] = array[0:4*4096, 1].reshape(128, 128)
for ii in range(time):
    array1[:, :, 2*ii] = array[4*4096*ii:4*4096*(ii+1), 0].reshape(128, 128) # noqa E501
    array1[:, :, 2*ii+1] = array[4*4096*ii:4*4096*(ii+1), 1].reshape(128, 128) # noqa E501
    array2[:, 2*ii] = arrayP[(N_part**2)*ii:(N_part**2)*(ii+1), 0]
    array2[:, 2*ii+1] = arrayP[(N_part**2)*ii:(N_part**2)*(ii+1), 1]
    array3[2*ii] = arrayP2[ii, 0]
    array3[2*ii+1] = arrayP2[ii, 1]
    if ii % 10 == 0:
        print(ii)


x = np.linspace(-5, 5, 128)
y = np.linspace(-5, 5, 128)
xx, yy = np.meshgrid(x, y)
# x = np.linspace(-0.2, 0.2, 128)
# y = np.linspace(-0.2, 0.2, 128)
# xx, yy = np.meshgrid(x, y)
fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.title(f'{N_part**2}, inertial tracers, t={round(ii*0.4/86.4, 3)} days, x_avg={round(array3[2*ii]*400, 3)} km, y_avg={round(array3[2*ii+1]*400, 3)} km') # noqa E501
    # plt.title(f'{N_part**2} passive tracers in an eastward-propagating field, t={ii*0.25}') # noqa E501
    plt.quiver(xx*400, yy*400, array1[:, :, 2*ii], array1[:, :, 2*ii+1], np.arctan2(array1[:, :, 2*ii], array1[:, :, 2*ii+1])) # noqa E501
    # plt.quiver(xx*400, yy*400, array1[:, :, 0], array1[:, :, 1], np.arctan2(array1[:, :, 0], array1[:, :, 1])) # noqa E501
    # plt.quiver(xx*400, yy*400, array1[:, :, 2*ii], array1[:, :, 2*ii+1]) # noqa E501
    # plt.scatter(array2[:, 4*ii]*400, array2[:, 4*ii+1]*400, c=np.multiply(array2[:, 4*ii+2]-0.05, array2[:, 4*ii+2]-0.05), marker='.', s=15, cmap='plasma') # noqa E501
    plt.scatter(array2[:, 2*ii]*400, array2[:, 2*ii+1]*400, c=['#000000'], s=1) # noqa E501
    plt.xlim(-2000, 2000)
    plt.ylim(-2000, 2000)
    print(ii)


anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, time),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=100)
anim.save('inertial.mp4', writer=writervideo)

# x = np.linspace(-5,5,64)
# y = np.linspace(-5,5,64)
# X, Y = np.meshgrid(x, y)
# plt.figure(figsize=(12, 12))
# plt.quiver(X, Y, array1[:,:,0], array1[:,:,1], np.arctan2(array1[:,:,0],array1[:,:,1])) # noqa E501
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.savefig('Velocity0.png')
