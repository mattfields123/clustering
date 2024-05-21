
import numpy as np
from matplotlib import pyplot as plt
from time import time
# with open('uvel.dat') as file_name:
    # uvel = np.loadtxt(file_name)
# 
# with open('vvel.dat') as file_name:
    # vvel = np.loadtxt(file_name)

vel_domain = 250

uvel = np.zeros((vel_domain,vel_domain))
vvel = np.zeros((vel_domain,vel_domain))

x = y = np.linspace(-5,5,vel_domain)

X, Y = np.meshgrid(x,y)


# PHI = np.sin(2*np.pi*X**2+2*np.pi*Y+3*Y**3)
# PHI = np.sin(X) + np.cos(Y)
# PHI = np.cos(2*np.pi*X+2*np.pi*Y)
PHI = X**4 + Y**2

def parabola(a,b):
    return [a**2,b**2,a*b,a,b,1]
    
t0 = time()

fpoints_x = []
fpoints_y = []
fpoints_stable_x = []
fpoints_stable_y = []
fpoints_saddle_x = []
fpoints_saddle_y = []
fpoints_unstable_x = []
fpoints_unstable_y = []


for i in range(vel_domain-2):
    for j in range(vel_domain-2):
        
        t_i = i + 1
        t_j = j + 1

        matrix = np.zeros((6,6))
        psi = np.zeros(6)

        psi[0] = PHI[t_i,t_j]
        psi[1] = PHI[t_i+1,t_j]
        psi[2] = PHI[t_i,t_j+1]
        psi[3] = PHI[t_i-1,t_j]
        psi[4] = PHI[t_i,t_j-1]
        psi[5] = PHI[t_i-1,t_j+1]


        matrix[0] = parabola(x[t_i],y[t_j])
        matrix[1] = parabola(x[t_i+1],y[t_j])
        matrix[2] = parabola(x[t_i],y[t_j+1])
        matrix[3] = parabola(x[t_i-1],y[t_j])
        matrix[4] = parabola(x[t_i],y[t_j-1])
        matrix[5] = parabola(x[t_i-1],y[t_j+1])

        if np.linalg.det(matrix) != 0:

            matrix_inv = np.linalg.inv(matrix)
        
            coefficients = np.matmul(matrix_inv,psi)

            A,B,C,D,E,F = coefficients

            matcoef = np.zeros((2,2))
            matcoef[0,0] = 2*A
            matcoef[0,1] = C
            matcoef[1,0] = C
            matcoef[1,1] = 2*B

            nonhom = np.zeros(2)
            nonhom[0] = -D
            nonhom[1] = -E
            matcoefinv = np.linalg.inv(matcoef)
            xy = np.matmul(matcoefinv,nonhom)
            x_fix = xy[0]
            y_fix = xy[1]

            # x_fix = (2*B*D-E*C)/(C**2-4*A*B)
            # y_fix = (2*A*E-D*C)/(C**2-4*A*B)
            if abs(x[t_i]-x_fix) < 10/vel_domain and abs(y[t_j]-y_fix) < 10/vel_domain:
                DET = np.linalg.det(matcoef)
                TR = np.trace(matcoef)

                if DET > 0:
                    if TR < 0:
                        fpoints_stable_x.append(x_fix)
                        fpoints_stable_y.append(y_fix)
                    if TR > 0:
                        fpoints_unstable_x.append(x_fix)
                        fpoints_unstable_y.append(y_fix)
                elif DET < 0:
                    fpoints_saddle_x.append(x_fix)
                    fpoints_saddle_y.append(y_fix)
                fpoints_x.append(x_fix)
                fpoints_y.append(y_fix)
        

print(fpoints_x,fpoints_y)
print((time()-t0))
plt.scatter(fpoints_x,fpoints_y)
plt.imshow(PHI.transpose(),extent=[-5,5,-5,5],origin='lower')
plt.show()