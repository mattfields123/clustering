
import numpy as np
from matplotlib import pyplot as plt
from time import time
import multiprocessing as mp
from copy import copy

time1 = time()

with open('stream.dat') as file_name:
    PSI = np.loadtxt(file_name)

with open('pot.dat') as file_name:
    PHI = np.loadtxt(file_name)

vel_domain = 200
tsteps = 1000

# t = np.linspace(0,5,tsteps)
x = np.linspace(-5,5,vel_domain)
y = np.linspace(-5,5,vel_domain)


# T,X,Y = np.meshgrid(t,x,y,indexing='ij')
# print(T.shape)


#PSI = np.cos(X) + np.sin(Y)
#PHI = np.sin(0*X+0*Y+0*T)

print(PSI.shape)
print(PHI.shape)
PSI = PSI.reshape((tsteps,vel_domain,vel_domain))
PHI = PHI.reshape((tsteps,vel_domain,vel_domain))


uvel = np.zeros((vel_domain,vel_domain))
vvel = np.zeros((vel_domain,vel_domain))


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

overall_array_x_stable = []
overall_array_y_stable = []
overall_array_x_unstable = []
overall_array_y_unstable = []
overall_array_x_saddle = []
overall_array_y_saddle = []

#for k in range(tsteps):

def compute_fixed(k):
    fpoints_x = []
    fpoints_y = []

    fpoints_stable_x = []
    fpoints_stable_y = []
    fpoints_unstable_x = []
    fpoints_unstable_y = []
    fpoints_saddle_x = []
    fpoints_saddle_y = []


    for i in range(1,vel_domain-1):
        for j in range(1,vel_domain-1):
            
            if ((np.mod(i,2) == 1) and (np.mod(j,3) == 1)) or ((np.mod(i,2) == 0) and (np.mod(j,3)==2)):


                t_i = i
                t_j = j
                matrix = np.zeros((6,6))
                psi = np.zeros(6)

                psi[0] = PSI[k,t_i-1,t_j-1]
                psi[1] = PSI[k,t_i-1,t_j]
                psi[2] = PSI[k,t_i,t_j-1]
                psi[3] = PSI[k,t_i,t_j+1]
                psi[4] = PSI[k,t_i+1,t_j-1]
                psi[5] = PSI[k,t_i+1,t_j]

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


                phi = np.zeros(6)

                phi[0] = PHI[k,t_i-1,t_j-1]
                phi[1] = PHI[k,t_i-1,t_j]
                phi[2] = PHI[k,t_i,t_j-1]
                phi[3] = PHI[k,t_i,t_j+1]
                phi[4] = PHI[k,t_i+1,t_j-1]
                phi[5] = PHI[k,t_i+1,t_j]

                if np.linalg.det(matrix) != 0:

                    matrix_inv = np.linalg.inv(matrix)
                


                    coefficients = np.matmul(matrix_inv,psi)
                    coefficientspot = np.matmul(matrix_inv,phi)
                    A1,B1,C1,D1,E1,F1 = coefficients
                    A2,B2,C2,D2,E2,F2 = coefficientspot
                    
                    matcoef = np.zeros((2,2))
                    nonhom = np.zeros(2)
                    nonhom[0] = -D2 - E1
                    nonhom[1] = D1 - E2

                    matcoef[0,0] = C1 + 2*A2
                    matcoef[0,1] = C2 + 2*B1
                    matcoef[1,0] = C2 - 2*A1
                    matcoef[1,1] = 2*B2 - C1

                    

                    
                    matcoefinv = np.linalg.inv(matcoef)
                    xy = np.matmul(matcoefinv,nonhom)
                    x_fix = xy[0]
                    y_fix = xy[1]

                    r = (np.sqrt(3)/2)*10/vel_domain

                    if (x[t_i]-x_fix)**2 + (y[t_j]-y_fix)**2 < r**2:
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
        
    fpoints_x_array = np.array(fpoints_x)
    fpoints_y_array = np.array(fpoints_y)

    fpoints_stable_x_array = np.array(fpoints_stable_x)
    fpoints_stable_y_array = np.array(fpoints_stable_y)

    fpoints_unstable_x_array = np.array(fpoints_unstable_x)
    fpoints_unstable_y_array = np.array(fpoints_unstable_y)


    fpoints_saddle_x_array = np.array(fpoints_saddle_x)
    fpoints_saddle_y_array = np.array(fpoints_saddle_y)


    print('stable numer = ',len(fpoints_stable_x))
    print('unstable numer = ',len(fpoints_unstable_x))
    print('saddle numer = ',len(fpoints_saddle_x))
    return fpoints_stable_x_array,fpoints_stable_y_array,fpoints_unstable_x_array,fpoints_unstable_y_array,fpoints_saddle_x_array,fpoints_saddle_y_array

pool = mp.Pool(30)

overall_array_x_stable,overall_array_y_stable,overall_array_x_unstable,overall_array_y_unstable,overall_array_x_saddle,overall_array_y_saddle = zip(*pool.map(compute_fixed,range(0,tsteps)))

combined_fixed_points = [overall_array_x_stable,overall_array_y_stable,overall_array_x_unstable,overall_array_y_unstable,overall_array_x_saddle,overall_array_y_saddle]


print('time to compute: ',time()-time1)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter

N_part = 500

with open("partx.dat") as file_name:
    partx = np.loadtxt(file_name)

with open("party.dat") as file_name:
    party = np.loadtxt(file_name)

partx = np.fmod(partx+1000,10)-5
party = np.fmod(party+1000,10)-5


print(partx.shape)

print(party.shape)



def field_radius(tstep):
    particles_in_vicinity = []
    particle_total = []
    new_fixed_points = copy(combined_fixed_points)
    for b in range(3):
        particles_in_vicinity = []
        a = 2*b
        for l in range(len(combined_fixed_points[a][tstep])):
            counter_l = 0
            for i in range(N_part):
                for j in range(N_part):
                    if (partx[tstep, i] - combined_fixed_points[a][tstep][l])**2 + (party[tstep,j]-combined_fixed_points[a+1][tstep][l])**2 < 0.01:
                        counter_l = counter_l + 1
            particles_in_vicinity.append(counter_l)
            if counter_l < 10:
                new_fixed_points[a][tstep] = np.delete(new_fixed_points[a][tstep], combined_fixed_points[a][tstep][l])
                new_fixed_points[a+1][tstep] = np.delete(new_fixed_points[a+1][tstep], combined_fixed_points[a+1][tstep][l])
        particle_total.append(particles_in_vicinity)


    return particle_total, new_fixed_points



fields_metric,new_fixed = field_radius(400)

print(len(new_fixed[1][400]))
print(len(combined_fixed_points[1][400]))


print(fields_metric)
print('sums',sum(fields_metric[0])/len(fields_metric[0]),sum(fields_metric[1])/len(fields_metric[1]),sum(fields_metric[2])/len(fields_metric[2]))

fig = plt.figure(figsize=(12, 12))
axes = plt.subplot()


def init_func():
    plt.cla()


def update_plot(ii):
    plt.cla()
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.scatter(400*partx[ii+1,:],400*party[ii+1,:],c='black',s=0.1)
    plt.scatter(400*overall_array_y_stable[ii],400*overall_array_x_stable[ii],c='red',s=10,alpha=0.5)
    plt.scatter(400*overall_array_y_unstable[ii],400*overall_array_x_unstable[ii],c='green',s=10,alpha=0.5)
    plt.scatter(400*overall_array_y_saddle[ii],400*overall_array_x_saddle[ii],c='blue',s=10,alpha=0.5)


    plt.xlim(-2000, 2000)
    plt.ylim(-2000, 2000)
  

anim = FuncAnimation(fig,
                     update_plot,
                     frames=np.arange(0, tsteps),
                     init_func=init_func)

writervideo = FFMpegWriter(fps=20)
anim.save('anim.mp4', writer=writervideo)




# plt.cla()
# plt.xlabel('X (km)')
# plt.ylabel('Y (km)')
# # plt.scatter(400*partx[ii+1,:],400*party[ii+1,:],c='black',s=0.1)
# plt.scatter(400*overall_array_x_stable[0],400*overall_array_y_stable[0],c='red',s=0.5)
# plt.scatter(400*overall_array_x_unstable[0],400*overall_array_y_unstable[0],c='green',s=0.5)
# plt.scatter(400*overall_array_x_saddle[0],400*overall_array_y_saddle[0],c='blue',s=0.5)
# plt.imshow(PSI[0].transpose(),extent=[-2000,2000,-2000,2000],origin='lower')

# plt.xlim(-2000, 2000)
# plt.ylim(-2000, 2000)
  
# plt.show()




print((time()-t0))



