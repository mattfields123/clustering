print('hello world')
import numpy as np
print(np.zeros((5,2)))
from matplotlib import pyplot as plt

x = np.linspace(0,10,100)
y = x**2

plt.plot(x,y)
plt.show()