module particle
use velocity
implicit none

contains
    function passive(partx, party, t, phase1, phase2, time1, time2, dt, dispersions, amplitudes) result(parts)
        real(dp) :: partx, party, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), parts(2)
        real(dp) :: k1(2), k2(2), k3(2), k4(2), dispersions(65,65), amplitudes(65,65)
        real :: dt

        parts = (/partx,party/)
        k1 = velocity_pointML(partx,party,t,phase1,phase2,time1,time2, dispersions, amplitudes)
        k2 = velocity_pointML(partx+k1(1)*dt/2,party+k1(2)*dt/2,t+dt/2,phase1,phase2,time1,time2,dispersions, amplitudes)
        k3 = velocity_pointML(partx+k2(1)*dt/2,party+k2(2)*dt/2,t+dt/2,phase1,phase2,time1,time2, dispersions, amplitudes)
        k4 = velocity_pointML(partx+k3(1)*dt,party+k3(2)*dt,t+dt,phase1,phase2,time1,time2, dispersions, amplitudes)
        parts = parts + dt*(k1+2*k2+2*k3+k4)/6
    end function passive


end module particle
