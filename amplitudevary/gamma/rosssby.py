import numpy as np
from matplotlib import pyplot as plt

x = np.linspace(0,10,100)
l = 1
R_d = 0.1

y = -2*np.pi*x*32/(4*((np.pi)**2)*(x**2+l**2)+1/(R_d**2))
z = -32/(4*np.pi**2*(x**2)+100)

plt.plot(x,z,)
plt.xlabel('K')
plt.ylabel('$\omega / 2\pi$')

plt.show()