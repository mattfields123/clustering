module particle
use velocity
implicit none

contains
    function passive(partx, party, t, phase1, phase2, time1, time2, dt) result(parts)
        real(dp) :: partx, party, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), parts(2)
        real(dp) :: k1(2), k2(2), k3(2), k4(2)
        real :: dt

        parts = (/partx,party/)
        k1 = velocity_pointML(partx,party,t,phase1,phase2,time1,time2)
        k2 = velocity_pointML(partx+k1(1)*dt/2,party+k1(2)*dt/2,t+dt/2,phase1,phase2,time1,time2)
        k3 = velocity_pointML(partx+k2(1)*dt/2,party+k2(2)*dt/2,t+dt/2,phase1,phase2,time1,time2)
        k4 = velocity_pointML(partx+k3(1)*dt,party+k3(2)*dt,t+dt,phase1,phase2,time1,time2)
        parts = parts + dt*(k1+2*k2+2*k3+k4)/6
    end function passive

    function soulsby(x, y, t, u, v, phase1, phase2, time1, time2) result(souls)
        real(dp) :: x, y, t, u, v, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), souls(2)
        real :: R, St
        R = 1./2.
        St = 946.

        souls = 1.5*R*acceleration_ML(x,y,t,phase1,phase2,time1,time2)
        souls = souls + (velocity_pointML(x,y,t,phase1,phase2,time1,time2) - (/u, v/))*St
    end function soulsby

    function inertial_soulsby(x, y, t, u, v, phase1, phase2, time1, time2, dt) result(posvel)
        real(dp) :: x, y, t, u, v, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        real(dp) :: k1(4), k2(4), k3(4), k4(4), posvel(4)
        real :: dt

        posvel = (/x, y, u, v/)

        k1 = (/u, v, soulsby(x,y,t,u,v,phase1,phase2,time1,time2)/)
        k2 = (/u+dt*k1(3)*0.5, v+dt*k1(4)*0.5,&
            &soulsby(x+dt*k1(1)*0.5,y+dt*k1(2)*0.5,t+dt*0.5,u+dt*k1(3)*0.5,v+dt*k1(4)*0.5,phase1,phase2,time1,time2)/)
        k3 = (/u+dt*k2(3)*0.5, v+dt*k2(4)*0.5,&
            &soulsby(x+dt*k2(1)*0.5,y+dt*k2(2)*0.5,t+dt*0.5,u+dt*k2(3)*0.5,v+dt*k2(4)*0.5,phase1,phase2,time1,time2)/)
        k4 = (/u+dt*k3(3), v+dt*k3(4), soulsby(x+dt*k3(1),y+dt*k3(2),t+dt,u+dt*k3(3),v+dt*k3(4),phase1,phase2,time1,time2)/)

        posvel = posvel + dt*(k1+2*k2+2*k3+k4)/6
    end function inertial_soulsby

end module particle