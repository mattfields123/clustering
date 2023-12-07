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
real(dp) :: dispersions(65,65), amplitudes(65,65)
integer :: c_x, c_y, c_t
real(dp) :: t
phase1 = tau*random_matrix(65,65)
phase1 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)
        
        
        
x_array = linspace(-5.0,5.0,10)        
y_array = linspace(-5.0,5.0,10)
call amplitudes_array(amplitudes)
call dispersion_relation_array(dispersions)

open(1, file = 'streamfunction.dat')
open(2, file = 'potentialfunction.dat')
do c_t = 1 , 10
call MaduLawrence_loop(time1,phase1,t)
call MaduLawrence_loop(time2,phase2,t)

do c_x = 1 , 10
do c_y = 1, 10
x = x_array(c_x)
y = y_array(c_y)
psi_result(c_x,c_y) = streamfunction(x,y,t,time1,phase1)

end do 
end do 

do c_x = 1,10
do c_y = 1,10

write(1,*) psi_result(c_x,c_y)

end do
end do

end do

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
intermediary = (1-gamma) * Madulawrence(time(c_k,c_l)-t,vu) * amplitudes(c_k,c_l)
psi = psi + intermediary * cos(tau*k*x+tau*l*y-dispersions(c_k,c_l)*t + phase1(c_k,c_l)) 
end do
end do
psi = psi/scaling
end function streamfunction



end program streampotential
