module velocity
use rossby_wave_attributes
use admin
use parameters
use random
implicit none

real(dp) :: mu = 0.3_dp
real(dp) :: scaling = 513.5_dp

!real(dp) :: scaling = 65._dp**2


contains
    function velocity_pointML(x, y, t, phase1, phase2, time1, time2, amplitudes, dispersions, g) result(velocity)
        real(dp) :: x, y, t, k, l, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: velocity(2), k_array(65), l_array(65)
        real(dp) :: amplitudes(65,65), dispersions(65,65)
        real(dp) :: g, tausol, taupot
        integer :: c_k, c_l
!     tausol = 8*atan(1.)*(1-gamma) !relies on gamma
 !    taupot = 8*atan(1.)*gamma !relies on gamma
  
        tausol = (1-gamma)*8*atan(1.)
        taupot = gamma*8*atan(1.)
      
        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        velocity = (/0.,0./)
        
        !$OMP PARALLEL DO private(psi1,psi2,k,l,c_k,c_l) reduction(+:velocity)
        do c_k = 1, 65
        k = k_array(c_k)
        do c_l = 1, 65
        l = l_array(c_l)
        psi1 = amp_scaling*amplitudes(c_k,c_l)*sin(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase1(c_k,c_l))
        psi2 = amp_scaling*amplitudes(c_k,c_l)*sin(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase2(c_k,c_l))
        end do
        end do
        !$OMP END PARALLEL DO
        
        velocity(1) = velocity(1) + (tausol*l - taupot*k*delta)*psi1 
        velocity(1) = velocity(1) - taupot*(1-delta)*k*psi2
        velocity(2) = velocity(2) - (tausol*k + taupot*l*delta)*psi1 
        velocity(2) = velocity(2) - taupot*l*(1-delta)*psi2

        velocity = velocity/scaling
    
    end function velocity_pointML
    
function streamfunction(x, y, t, phase1, phase2, time1, time2, amplitudes, dispersions, g) result(velocity)
        real(dp) :: x, y, t, k, l, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: velocity, k_array(65), l_array(65)
        real(dp) :: amplitudes(65,65), dispersions(65,65)
        real(dp) :: g, tausol, taupot
        integer :: c_k, c_l

        tausol = (1-gamma)*8*atan(1.)
        taupot = gamma*8*atan(1.)
      
        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        velocity = (/0.,0./)
        
        !$OMP PARALLEL DO private(psi1,psi2,k,l,c_k,c_l) reduction(+:velocity)
        do c_k = 1, 65
        k = k_array(c_k)
        do c_l = 1, 65
        l = l_array(c_l)
        psi1 = amp_scaling*amplitudes(c_k,c_l)*cos(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase1(c_k,c_l))

        velocity = velocity + psi1
        end do
        end do
        !$OMP END PARALLEL DO
        
        velocity = (1-gamma)*velocity/scaling
    
    end function streamfunction

function potentialfunction(x, y, t, phase1, phase2, time1, time2, amplitudes, dispersions, g) result(velocity)
        real(dp) :: x, y, t, k, l, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: velocity, k_array(65), l_array(65)
        real(dp) :: amplitudes(65,65), dispersions(65,65)
        real(dp) :: g, tausol, taupot
        integer :: c_k, c_l

        tausol = (1-gamma)*8*atan(1.)
        taupot = gamma*8*atan(1.)
      
        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        velocity = (/0.,0./)
        
        !$OMP PARALLEL DO private(psi1,psi2,k,l,c_k,c_l) reduction(+:velocity)
        do c_k = 1, 65
        k = k_array(c_k)
        do c_l = 1, 65
        l = l_array(c_l)
        psi1 = amp_scaling*amplitudes(c_k,c_l)*cos(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase1(c_k,c_l))
        psi2 = amp_scaling*amplitudes(c_k,c_l)*cos(tau*k*x + tau*l*y - dispersions(c_k,c_l)*t + phase2(c_k,c_l))
        
        velocity = delta * psi1 + (1 - delta) * psi2
        end do
        end do
        !$OMP END PARALLEL DO
        


        velocity = gamma*velocity/scaling
    
    end function potentialfunction


    
    subroutine MaduLawrence_loop(time, phase, t, vu)
        real(dp) :: vu,  t, phase(65,65), time(65,65), div, ML
        integer :: counter_k, counter_l
        do counter_k = 1, 65
            do counter_l = 1, 65
                if (time(counter_k,counter_l)-t<0) then
                    time(counter_k,counter_l) = time(counter_k,counter_l) + vu
                    phase(counter_k,counter_l) = random_number1()*tau
                end if
            end do
        end do
    end subroutine MaduLawrence_loop

end module velocity
