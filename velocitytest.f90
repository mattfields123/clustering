program velocitytest
use constants
use parameterslarge
use admin
use rossby_wave
use velocity

integer :: tsteps, meshsize









contains 

function velocitycomp(meshsize,tsteps) return done
integer :: done
integer :: tsteps, meshsize
real(dp), allocatable :: t_array(:)
real(dp), allocatable :: x_array(:), y_array(:)
real(dp) :: vel(2), t, x, y, time1(65,65), time2(65,65), phase1(65,65), phase2(65,65)
integer :: t_c, x_c, y_c
real(dp) :: dt = 0.25



t_array = linspace(0.,tsteps*dt,tsteps)
x_array = linspace(-5.,5.,meshsize)
y_array = linspace(-5.,5.,meshsize)






do t_c = 1, tsteps

t = t_array(t_c)


do x_c = 1, meshsize

x = x_array(x_c)

do y_c = 1, meshsize

y = y_array(y_c)

vel = velocitynoML(x,y,t,time1,time2,phase1,phase2)





end do
end do
end do







end function velocitycomp


end program velocity test

