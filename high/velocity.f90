module velocity
use rossby_wave
use rossby_wave_attributes
use admin
use parameters
use random

implicit none

real(dp) :: mu = 0.3_dp
real(dp) :: scaling = 513.5_dp

!real(dp) :: scaling = 65._dp**2
real(dp) :: tausol = 8*atan(1.)*(1-0.15) !relies on gamma
real(dp) :: taupot = 8*atan(1.)*0.15 !relies on gamma

contains
    function velocity_pointML(x, y, t, phase1, phase2, time1, time2, dispersions, amplitudes) result(velocity)
        real(dp) :: x, y, t, k, l, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: velocity(2), k_array(65), l_array(65), dispersions(65,65), amplitudes(65,65)
        integer :: c_k, c_l
        
        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        velocity = (/0.,0./)
        
        
        do c_k = 1, 65
        k = k_array(c_k)
        do c_l = 1, 65
        l = l_array(c_l)
        psi1 = amplitudes(c_k,c_l)*sin(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase1(c_k,c_l))
        psi2 = amplitudes(c_k,c_l)*sin(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase2(c_k,c_l))
        velocity(1) = velocity(1) + (tausol*l - taupot*k*delta)*MaduLawrence(time1(c_k,c_l)-t,vu)*psi1 
        velocity(1) = velocity(1) - taupot*(1-delta)*k*MaduLawrence(time2(c_k,c_l)-t,vu)*psi2
        velocity(2) = velocity(2) - (tausol*k + taupot*l*delta)*MaduLawrence(time1(c_k,c_l)-t,vu)*psi1 
        velocity(2) = velocity(2) - taupot*l*(1-delta)*MaduLawrence(time2(c_k,c_l)-t,vu)*psi2
        end do
        end do
        
        velocity = velocity/scaling
    
    end function velocity_pointML

subroutine MaduLawrence_loop(time,phase,t)
        real(dp) :: t, phase(65,65), time(65,65), div, ML
        integer :: c_k, c_l

do c_k = 1,65
do c_l = 1,65

if (time(c_k,c_l)-t<0) then
time(c_k,c_l) = time(c_k,c_l) + vu
phase(c_k,c_l) = random_number1()*tau
end if

end do
end do
end subroutine MaduLawrence_loop 

end module velocity
