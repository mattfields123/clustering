program streampotential
use constants
use admin
use random
use rossby_wave
use velocity
implicit none

! Actual program

real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
real(dp) :: psi_result(10,10), x, y, x_array(10), y_array(10)
real(dp) :: phi_result(10,10)
real(dp) :: dispersions(65,65), amplitudes(65,65)
integer :: c_x, c_y, c_t, timesteps, grid
real(dp) :: t_array(10)
real(dp) :: t
real :: dt
real(dp) ::  g=0.1
real(dp) :: nu=1.
phase1 = tau*random_matrix(65,65)
phase1 = tau*random_matrix(65,65)
time1 = nu*random_matrix(65,65)
time2 = nu*random_matrix(65,65)
        
timesteps = 50
grid = 100

dt = 0.25


x_array = linspace(-5.0,5.0,grid)        
y_array = linspace(-5.0,5.0,grid)
t_array = linspace(0.,(timesteps-1.)*dt,timesteps)


call amplitudes_array(amplitudes)
call dispersion_relation_array(dispersions)


open(1, file = 'streamfunction.dat')
open(2, file = 'potentialfunction.dat')
!open(3, file= 'amplitudes.dat')
!open(4, file = 'dispersions.dat')

do c_t = 1 , timesteps
t = t_array(c_t)

call MaduLawrence_loop(time1,phase1,t,nu)
call MaduLawrence_loop(time2,phase2,t,nu)
! Where there is 1 we have vu (testing)
do c_x = 1 , grid
do c_y = 1, grid
do hex = 1,6


x = x_array(c_x)
y = y_array(c_y)
psi_result(c_x,c_y) = streamfunction(x,y,t,time1,phase1)
phi_result(c_x,c_y) = potentialfunction(x,y,t,time1,time2,phase1,phase2)


end do
end do 
end do 
print*, phi_result(2,2)
print*, psi_result(2,2)


do c_x = 1,grid
do c_y = 1,grid

write(1,*) psi_result(c_x,c_y)
write(2,*) phi_result(c_x,c_y)
end do
end do

end do

!write(3,*) amplitudes
!write(4,*) dispersions

close(1)
close(2)



contains 

function streamfunction(x,y,t,time,phase1) result(psi)
implicit none
integer :: c_k,c_l
real(dp) :: t
real(dp) :: time(65,65), phase1(65,65)
real(dp) :: psi, k, l, x, y, scaling
real(dp) :: k_array(65), l_array(65)
real(dp) :: intermediary
k_array = linspace(-3.2,3.2,65)
l_array = linspace(-3.2,3.2,65)
psi = 0.
scaling = 513.5 
do c_k = 1,65
do c_l = 1,65
k = k_array(c_k)
l = l_array(c_l)
intermediary = amplitudes(c_k,c_l)
psi = psi + intermediary * cos(tau*k*x+tau*l*y-dispersions(c_k,c_l)*t + phase1(c_k,c_l)) 
end do
end do
psi = psi/scaling
end function streamfunction

function potentialfunction(x,y,t,time1,time2,phase1,phase2) result(phi)
implicit none
integer :: c_k, c_l
real(dp) :: t
real(dp) :: time1(65,65),time2(65,65), phase1(65,65), phase2(65,65)
real(dp) :: phi, k, l , x , y, scaling
real(dp) :: k_array(65), l_array(65)
real(dp) :: int1,int2,int11,int22,inter1,inter2

k_array = linspace(-3.2,3.2,65) 
l_array = linspace(-3.2,3.2,65) 
phi = 0.
scaling = 513.5

do c_k = 1,65
do c_l = 1,65
k = k_array(c_k)
l = l_array(c_l)
int1 = delta 
int2 = (1-delta)
int11 = cos(tau*k*x+tau*l*y-dispersions(c_k,c_l)*t+phase1(c_k,c_l))
int22 = cos(tau*k*x+tau*l*y-dispersions(c_k,c_l)*t+phase2(c_k,c_l))
inter1 = int1*int11
inter2 = int2*int22
phi = phi + inter1 + inter2


end do
end do
phi = phi/scaling


end function potentialfunction


end program streampotential
