
import numpy as np
from matplotlib import pyplot as plt
from time import time



with open('stream.dat') as file_name:
    PSI = np.loadtxt(file_name)

with open('pot.dat') as file_name:
    PHI = np.loadtxt(file_name)

vel_domain = 1000
tsteps = 500


PSI = PSI.reshape((tsteps,vel_domain,vel_domain))
PHI = PHI.reshape((tsteps,vel_domain,vel_domain))


# uvel = np.zeros((vel_domain,vel_domain))
# vvel = np.zeros((vel_domain,vel_domain))



x = y = np.linspace(-5,5,vel_domain)
# t = np.linspace(0,5,tsteps)
# X, Y, T = np.meshgrid(x,y,t)


# PHI = np.sin(2*np.pi*X**2+2*np.pi*Y+3*Y**3)
# PHI = np.sin(X) + np.cos(Y)
# PHI = np.cos(2*np.pi*X+2*np.pi*Y)
# PSI = np.cos(0*X) + np.sin(0*Y)

# PHI = np.cos(X+Y-T)


def parabola(a,b):
    return [a**2,b**2,a*b,a,b,1]
    
t0 = time()

overall_array_x = []
overall_array_y = []




for k in range(500):
    fpoints_x = []
    fpoints_y = []
    for i in range(vel_domain-2):
        for j in range(vel_domain-2):
            
            t_i = i + 1
            t_j = j + 1

            matrix = np.zeros((6,6))
            psi = np.zeros(6)

            psi[0] = PSI[k,t_i,t_j]
            psi[1] = PSI[k,t_i+1,t_j]
            psi[2] = PSI[k,t_i,t_j+1]
            psi[3] = PSI[k,t_i-1,t_j]
            psi[4] = PSI[k,t_i,t_j-1]
            psi[5] = PSI[k,t_i-1,t_j-1]


            matrix[0] = parabola(x[t_i],y[t_j])
            matrix[1] = parabola(x[t_i+1],y[t_j])
            matrix[2] = parabola(x[t_i],y[t_j+1])
            matrix[3] = parabola(x[t_i-1],y[t_j])
            matrix[4] = parabola(x[t_i],y[t_j-1])
            matrix[5] = parabola(x[t_i-1],y[t_j-1])

            phi = np.zeros(6)

            phi[0] = PHI[k,t_i,t_j]
            phi[1] = PHI[k,t_i+1,t_j]
            phi[2] = PHI[k,t_i,t_j+1]
            phi[3] = PHI[k,t_i-1,t_j]
            phi[4] = PHI[k,t_i,t_j-1]
            phi[5] = PHI[k,t_i-1,t_j-1]

            if np.linalg.det(matrix) != 0:

                matrix_inv = np.linalg.inv(matrix)
            


                coefficients = np.matmul(matrix_inv,psi)
                coefficientspot = np.matmul(matrix_inv,phi)
                A1,B1,C1,D1,E1,F1 = coefficients
                A2,B2,C2,D2,E2,F2 = coefficientspot
                
                matcoef = np.zeros((2,2))
                nonhom = np.zeros(2)
                nonhom[0] = -D1-E2
                nonhom[1] = E1-D2

                matcoef[0,0] = 2*A1+C2
                matcoef[0,1] = C1+2*B2
                matcoef[1,0] = -C1 + 2*A2
                matcoef[1,1] = -2*B1 + C2

                
                matcoefinv = np.linalg.inv(matcoef)
                xy = np.matmul(matcoefinv,nonhom)
                x_fix = xy[0]
                y_fix = xy[1]
                
                if x[t_i-1] <= x_fix <= x[t_i+1]:
                    if y[t_j-1] <= y_fix <= y[t_j+1]:
                        fpoints_x.append(x_fix)
                        fpoints_y.append(y_fix)
    overall_array_x.append(fpoints_x)            
    overall_array_y.append(fpoints_y)
    print(len(fpoints_x))





print((time()-t0))



