from matplotlib import pyplot as plt
import numpy as np

with open('stream.dat') as f:
    psi = np.loadtxt(f)
with open('pot.dat') as f:
    phi = np.loadtxt(f)

im = plt.imshow(psi[100].reshape((200,200)))
plt.colorbar()
plt.savefig('figure.png')
