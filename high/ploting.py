
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
with open("potentialfunction(1).dat") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
with open("streamfunction(2).dat") as file_name:
    darray = np.loadtxt(file_name, delimiter=",")
# print(np.shape(array))

points = np.zeros((3,500,500))
points = array.reshape((3,500,500))

points2 = darray.reshape((3,500,500))


# print(points[0])

x_coords = np.linspace(-5,5,500)
y_coords = np.linspace(5,5,500)

print(points[1][3][1:10])
print(points[2][3][1:10])

fig, ax = plt.subplots(nrows=2,ncols=1)

c = ax[0].imshow(points[0],cmap='jet',extent=[-5,5,-5,5])

d = ax[1].imshow(points2[0],cmap='jet',extent=[-5,5,-5,5])

plt.savefig("streampotential.png")




# fig = plt.figure(figsize=(12, 12))
# axes = plt.subplot()
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')


# def init_func():
#     plt.cla()


# def update_plot(ii):
#     plt.cla()
#     plt.xlabel('X (km)')
#     plt.ylabel('Y (km)')
#     plt.title('Streamfunction') # noqa E501
#     plt.imshow(points[ii],cmap='jet',extent=[-5,5,-5,5])
#     print(ii)
#     print(points[ii][10][5])


# anim = FuncAnimation(fig,
#                      update_plot,
#                      frames=np.arange(0, 3),
#                      init_func=init_func)

# writervideo = FFMpegWriter(fps=2)
# anim.save('inertial.mp4', writer=writervideo)
