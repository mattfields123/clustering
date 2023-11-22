program streampotential
use constants
use admin
implicit none

! Actual program

real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
real(dp) :: psi_result(10,10)

phase1 = tau*random_matrix(65,65)
phase1 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)
        
        
        
x_array = linspace(-5.0,5.0,10)        
y_array = linspace(-5.0,5.0,10)
amplitudes = amplitudes_array

open(1, file = 'streamfunction.dat')
open(2, file = 'potentialfunction.dat')
do t = 1 , 10
MaduLawrence_loop(time1,phase1,t)
MaduLawrence_loop(time2,phase2,t)

do c_x = 1 , 10
do c_y = 1, 10
x = x_array(c_x)
y = y_array(c_y)
psi_result(c_x,c_y) = streamfunction(x,y,t,time,phase1)

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
integer :: a, N, time(65,65), phase1(65,65)
real(dp) :: psi, k, l, x, y

k_array = linspace(-3.2,3.2,65)
l_array = linspace(-3.2,3.2,65)
psi = 0.
 
do c_k = 1,65
do c_l = 1,65
k = k_array(c_k)
l = l_array(c_l)
psi = psi + (1-gamma) * Madulawrence(time(c_k,c_l)-t,vu) amplitudes(k,l) * cos(tau*k*x+tau*l*y-dispersions(k,l)*t + phase1(c_k,c_l)) 

end do
end do
psi = psi/scaling
end function streamfunction



end program streampotential
