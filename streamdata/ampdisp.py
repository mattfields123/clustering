import numpy as np
from matplotlib import pyplot as plt

with open("/Users/bunny/Documents/msci/mscigit/streamdata/amplitudes.dat") as file_name:
    amp = np.loadtxt(file_name)
with open("/Users/bunny/Documents/msci/mscigit/streamdata/dispersions.dat") as file_name:
    disp = np.loadtxt(file_name)

disp = np.reshape(disp,(65,65))

plt.imshow(disp)
plt.show()