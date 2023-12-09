module velocity
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
    function streamfunction(x,y,t,phase1,time1) result(psi)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), psi, ML
        real(dp) :: k_array(65), l_array(65)
        integer :: counter_l, counter_k

        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        psi = 0.

        !$OMP PARALLEL DO private(ML,k,l) reduction(+:psi)
        do counter_k = 1, 65
            k = k_array(counter_k)
            do counter_l = 1, 65
                l = l_array(counter_l)
                ML = (1-gamma)*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)
                psi = psi + ML*amplitude(k,l)*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l))
            end do
        end do
        !$OMP END PARALLEL DO
        psi = psi/scaling
    end function streamfunction
    
    function potential(x,y,t,phase1,phase2,time1,time2) result(phi)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), phase2(65,65), time2(65,65), phi, ML1, ML2
        real(dp) :: k_array(65), l_array(65)
        integer :: counter_l, counter_k

        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        phi = 0.

        !$OMP PARALLEL DO private(ML1,ML2,k,l) reduction(+:phi)
        do counter_k = 1, 65
            k = k_array(counter_k)
            do counter_l = 1, 65
                l = l_array(counter_l)
                ML1 = gamma*delta*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)
                ML2 = gamma*(1-delta)*MaduLawrence(time2(counter_k,counter_l)-t,vu,mu)
                phi = phi + ML2*amplitude(k,l)*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase2(counter_k,counter_l))
                phi = phi + ML1*amplitude(k,l)*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l))
            end do
        end do
        !$OMP END PARALLEL DO
        phi = phi/scaling
    end function potential
    
    function velocity_pointML(x, y, t, phase1, phase2, time1, time2, amplitudes, dispersions) result(velocity)
        real(dp) :: x, y, t, k, l, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: velocity(2), k_array(65), l_array(65)
        real(dp) :: amplitudes(65,65), dispersions(65,65)
        integer :: counter_k, counter_l
        
        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        velocity = (/0.,0./)
        
        !$OMP PARALLEL DO private(psi1,psi2,k,l) reduction(+:velocity)
        do counter_k = 1, 65
        k = k_array(counter_k)
        do counter_l = 1, 65
        l = l_array(counter_l)
        psi1 = amplitudes(k,l)*sin(tau*k*x + tau*l*y - dispersions(k,l)*t + phase1(counter_k,counter_l))
        psi2 = amplitudes(k,l)*sin(tau*k*x + tau*l*y - dispersions(k,l)*t + phase2(counter_k,counter_l))
        velocity(1) = velocity(1) + (tausol*l - taupot*k*delta)*psi1 
        velocity(1) = velocity(1) - taupot*(1-delta)*k*psi2
        velocity(2) = velocity(2) - (tausol*k + taupot*l*delta)*psi1 
        velocity(2) = velocity(2) - taupot*l*(1-delta)*psi2
        end do
        end do
        !$OMP END PARALLEL DO
        velocity = velocity/scaling
    
    end function velocity_pointML
    
    function velocity_divergence(x,y,t,phase1,phase2,time1,time2) result(div)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), phase2(65,65), time2(65,65), div, ML1, ML2
        real(dp) :: k_array(65), l_array(65), tausqpot
        integer :: counter_l, counter_k

        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        div = 0.
        tausqpot = tau*taupot

        !$OMP PARALLEL DO private(ML1,ML2,k,l,counter_k,counter_l) reduction(+:div)
        do counter_k = 1, 65
            k = k_array(counter_k)
            do counter_l = 1, 65
                l = l_array(counter_l)
                ML1 = (k**2 + l**2)*tausqpot*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)*amplitude(k,l)
                ML2 = (k**2 + l**2)*tausqpot*MaduLawrence(time2(counter_k,counter_l)-t,vu,mu)*amplitude(k,l)
                div = div - delta*ML1*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l)) &
                & - (1-delta)*ML2*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase2(counter_k,counter_l))
            end do
        end do
        !$OMP END PARALLEL DO
        div = div/scaling
    end function velocity_divergence
    
    function vorticity(x,y,t,phase1,time1) result(zeta)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), zeta, ML
        real(dp) :: k_array(65), l_array(65), tausqsol
        integer :: counter_l, counter_k

        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        zeta = 0.
        tausqsol = tau*tausol

        !$OMP PARALLEL DO private(ML,k,l) reduction(+:zeta)
        do counter_k = 1, 65
            k = k_array(counter_k)
            do counter_l = 1, 65
                l = l_array(counter_l)
                ML = tausqsol*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)*amplitude(k,l)
                zeta = zeta - (k**2 + l**2)*ML*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l))
            end do
        end do
        !$OMP END PARALLEL DO
        zeta = zeta/scaling
    end function vorticity
    
    function normal_strain(x,y,t,phase1,phase2,time1,time2) result(ns)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), phase2(65,65), time2(65,65), ns, ML1, ML2
        real(dp) :: k_array(65), l_array(65), tausqsol, tausqpot
        integer :: counter_l, counter_k

        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        ns = 0.
        tausqsol = tau*tausol
        tausqpot = tau*taupot

        !$OMP PARALLEL DO private(ML1,ML2,k,l) reduction(+:ns)
        do counter_k = 1, 65
            k = k_array(counter_k)
            do counter_l = 1, 65
                l = l_array(counter_l)
                ML1 = (2*tausqsol*l*k + tausqpot*delta*(l**2-k**2))*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)
                ML2 = tausqpot*(l**2 - k**2)*(1-delta)*MaduLawrence(time2(counter_k,counter_l)-t,vu,mu)
                ns = ns + ML1*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l))*amplitude(k,l) &
                & + ML2*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase2(counter_k,counter_l))*amplitude(k,l)
            end do
        end do
        !$OMP END PARALLEL DO
        ns = ns/scaling
    end function normal_strain
    
    function shear_strain(x,y,t,phase1,phase2,time1,time2) result(ss)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), phase2(65,65), time2(65,65), ss, ML1, ML2
        real(dp) :: k_array(65), l_array(65), tausqsol, tausqpot
        integer :: counter_l, counter_k

        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        ss = 0.
        tausqsol = tau*tausol
        tausqpot = tau*taupot

        !$OMP PARALLEL DO private(ML1,ML2,k,l) reduction(+:ss)
        do counter_k = 1, 65
            k = k_array(counter_k)
            do counter_l = 1, 65
                l = l_array(counter_l)
                ML1 = (tausqsol*(l**2 - k**2)-2*tausqpot*k*l*delta)*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)
                ML2 = 2*tausqpot*l*k*(1-delta)*MaduLawrence(time2(counter_k,counter_l)-t,vu,mu)
                ss = ss + ML1*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l))*amplitude(k,l) &
                & - ML2*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase2(counter_k,counter_l))*amplitude(k,l)
            end do
        end do
        !$OMP END PARALLEL DO
        ss = ss/scaling
    end function shear_strain
    
    function Okubo_Weiss(x,y,t,phase1,phase2,time1,time2) result(OkW)
        real(dp) :: x, y, t, k, l, phase1(65,65), time1(65,65), phase2(65,65), time2(65,65), OkW, zeta, ns, ss

        ns = normal_strain(x,y,t,phase1,phase2,time1,time2)
        ss = shear_strain(x,y,t,phase1,phase2,time1,time2)
        zeta = vorticity(x,y,t,phase1,time1)
        OkW = ns**2 + ss**2 - zeta**2
    end function Okubo_Weiss

    function speed_point(x, y, t, phase1, phase2, time1, time2) result(spd)
        real(dp) :: x, y, t, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), spd

        velocity = velocity_pointML(x, y, t, phase1, phase2, time1, time2)
        spd = (velocity(1)**2 + velocity(2)**2)**0.5
    end function speed_point

    subroutine MaduLawrence_loop(time, phase, t)
        real(dp) :: t, phase(65,65), time(65,65), div, ML
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

    function vortex_vel(x,y) result(velocity)
        real(dp) :: x, y, velocity(2), r
        
        r = (x**2 + y**2)/10
        if (r > 0.001) then
            velocity = (/-y/r,x/r/)
        else
            velocity = (/0.,0./)
        end if
    end function vortex_vel
    
    function acceleration_ML(x, y, t, phase1, phase2, time1, time2) result(acceleration)
        real(dp) :: x, y, t, k, l, psi1, psi2, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: acceleration(2), velocity(2), k_array(65), l_array(65), psi(4)
        integer :: counter_k, counter_l
        
        k_array = linspace(-3.2,3.2,65)
        l_array = linspace(-3.2,3.2,65)
        acceleration = (/0.,0./)
        velocity = velocity_pointML(x, y, t, phase1, phase2, time1, time2)
        
        do counter_k = 1, 65
        k = k_array(counter_k)
        do counter_l = 1, 65
        l = l_array(counter_l)
        psi1 = tausol*amplitude2(k,l)*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase1(counter_k,counter_l))
        psi2 = taupot*amplitude2(k,l)*cos(tau*k*x + tau*l*y - dispersion(k,l)*t + phase2(counter_k,counter_l))
        psi(1) = l*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)*psi1 
        psi(2) = -k*MaduLawrence(time2(counter_k,counter_l)-t,vu,mu)*psi2
        psi(3) = -k*MaduLawrence(time1(counter_k,counter_l)-t,vu,mu)*psi1 
        psi(4) = -l*MaduLawrence(time2(counter_k,counter_l)-t,vu,mu)*psi2
        acceleration(1) = acceleration(1) - dispersion(k,l)*(psi(1) + psi(2))
        acceleration(1) = acceleration(1) + tau*k*(psi(1) + psi(2))*velocity(1)
        acceleration(1) = acceleration(1) + tau*l*(psi(1) + psi(2))*velocity(2)
        acceleration(2) = acceleration(2) - dispersion(k,l)*(psi(3) + psi(4))
        acceleration(2) = acceleration(2) + tau*k*(psi(3) + psi(4))*velocity(1)
        acceleration(2) = acceleration(2) + tau*l*(psi(3) + psi(4))*velocity(2)
        end do
        end do
        acceleration = acceleration/scaling
    
    end function acceleration_ML

    function vortex_acc(x,y) result(acc)
        real(dp) :: x, y, acc(2), r

        r = (x**2 + y**2)/10
        if (r > 0.001) then
            acc = (/-x/r,-y/r/)
        else
            acc = (/0.,0./)
        end if
    end function vortex_acc

end module velocity
