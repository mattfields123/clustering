
import numpy as np
from matplotlib import pyplot as plt
from time import time
# with open('uvel.dat') as file_name:
    # uvel = np.loadtxt(file_name)
# 
# with open('vvel.dat') as file_name:
    # vvel = np.loadtxt(file_name)

vel_domain = 51

uvel = np.zeros((vel_domain,vel_domain))
vvel = np.zeros((vel_domain,vel_domain))
y_domain = int(vel_domain*2/np.sqrt(3))
x = np.linspace(-5,5,vel_domain)
y = np.linspace(-5,5,y_domain)




# PHI = np.sin(2*np.pi*X**2+2*np.pi*Y+3*Y**3)
# PHI = np.sin(X) + np.cos(Y)
# PHI = np.cos(2*np.pi*X+2*np.pi*Y)
PHI = np.zeros((vel_domain,y_domain))

for i in range(vel_domain):
    for j in range(y_domain):
        if (np.mod(i,2) > 0):
            PHI[i,j] = np.cos(np.pi/4*x[i])+np.sin(np.pi/4*y[j]) 
        else:
            PHI[i,j] = np.cos(np.pi/4*(x[i]+5/vel_domain)) + np.sin(np.pi/4*y[j])

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

total = 0

for i in range(1,vel_domain-1):
    for j in range(1,y_domain-1):
        
        if ((np.mod(i,2) == 1) and (np.mod(j,3) == 1)) or ((np.mod(i,2) == 0) and (np.mod(j,3)==2)):

            
            t_i = i
            t_j = j
            matrix = np.zeros((6,6))
            psi = np.zeros(6)

            psi[0] = PHI[t_i-1,t_j-1]
            psi[1] = PHI[t_i-1,t_j]
            psi[2] = PHI[t_i,t_j-1]
            psi[3] = PHI[t_i,t_j+1]
            psi[4] = PHI[t_i+1,t_j-1]
            psi[5] = PHI[t_i+1,t_j]

            if np.mod(t_i,2) == 0:

                matrix[0] = parabola(x[t_i-1],y[t_j-1])
                matrix[1] = parabola(x[t_i-1],y[t_j])
                matrix[2] = parabola(x[t_i]+5/vel_domain,y[t_j-1])
                matrix[3] = parabola(x[t_i]+5/vel_domain,y[t_j+1])
                matrix[4] = parabola(x[t_i+1],y[t_j-1])
                matrix[5] = parabola(x[t_i+1],y[t_j])

            else:
                matrix[0] = parabola(x[t_i-1]+5/vel_domain,y[t_j-1])
                matrix[1] = parabola(x[t_i-1]+5/vel_domain,y[t_j])
                matrix[2] = parabola(x[t_i],y[t_j-1])
                matrix[3] = parabola(x[t_i],y[t_j+1])
                matrix[4] = parabola(x[t_i+1]+5/vel_domain,y[t_j-1])
                matrix[5] = parabola(x[t_i+1]+5/vel_domain,y[t_j])


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
                r = (np.sqrt(3)/2)* 10 / vel_domain
                R = 10 / vel_domain
                if np.mod(t_i,2)  == 0:
                    x_loc = x[t_i] + 5 / vel_domain
                else:
                    x_loc = x[t_i]

                if ((x_loc-x_fix)**2 + (y[t_j]-y_fix)**2) < r**2:
                    DET = np.linalg.det(matcoef)
                    TR = np.trace(matcoef)

                    if DET > 0:
                        if TR < 0:
                            fpoints_stable_x.append(x_fix)
                            fpoints_stable_y.append(y_fix)
                            print('stable')
                        if TR > 0:
                            fpoints_unstable_x.append(x_fix)
                            fpoints_unstable_y.append(y_fix)
                            print('unstable')
                    elif DET < 0:
                        fpoints_saddle_x.append(x_fix)
                        fpoints_saddle_y.append(y_fix)
                        print('saddle')
                    fpoints_x.append(x_fix)
                    fpoints_y.append(y_fix)
                    print(t_i,t_j)
                    total = total + 1 
                elif ((x_loc-x_fix)**2 + (y[t_j]-y_fix)**2) < R**2:
                    print(x_fix,y_fix,'lmao')

print(fpoints_x,fpoints_y,'fpoints')
print((time()-t0))
print(total)
plt.scatter(fpoints_x,fpoints_y)
plt.imshow(PHI.transpose(),extent=[-5,5,-5,5],origin='lower')
plt.show()