from matplotlib import pyplot as plt
import numpy as np

with open('stream.dat') as f:
    psi = np.loadtxt(f)
with open('pot.dat') as f:
    phi = np.loadtxt(f)

plt.imshow(psi[100])
plt.savefig('figure.png')
