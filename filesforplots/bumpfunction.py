import numpy as np
from matplotlib import pyplot as plt
n =1000
x = np.linspace(0.0001,9.9999,n)

y = np.zeros(n)

for t in range(n):
    if x[t] < 0.1:
        y[t] = np.exp(1-1/(1-(x[t]-0.1)**2))
    elif  9.9 > x[t] >= 0.1:
        y[t] = 1
    else:
        y[t] = np.exp(1-1/(1-(x[t]-10+0.1)**2))

plt.plot(x,y)
plt.show()