program velocitytest
use constants
use parameters
use admin
use rossby_wave
use velocity

integer :: tsteps, meshsize

meshsize = 100
tsteps = 100

call velocitycomp(meshsize,tsteps)








contains 

subroutine velocitycomp(meshsize,tsteps)
integer :: done
integer :: tsteps, meshsize
real(dp), allocatable :: t_array(:)
real(dp), allocatable :: x_array(:), y_array(:)
real(dp) :: vel(2), t, x, y, time1(65,65), time2(65,65), phase1(65,65), phase2(65,65), dispersions(65,65), amplitudes(65,65)
integer :: t_c, x_c, y_c
real(dp) :: g=0.
real :: dt = 0.25
real(dp), allocatable :: vel_array(:,:)

allocate(t_array(tsteps))
allocate(x_array(meshsize))
allocate(y_array(meshsize))
allocate(vel_array(meshsize,meshsize))


t_array = linspace(0.,(tsteps-1)*dt,tsteps)
x_array = linspace(-5.,5.,meshsize)
y_array = linspace(-5.,5.,meshsize)


call dispersion_relation_array(dispersions)
call amplitudes_array(amplitudes)
    

open(1,file='velocities.dat')


do t_c = 1, tsteps

t = t_array(t_c)

call MaduLawrence_loop(time1, phase1, t, vu)
call MaduLawrence_loop(time2, phase2, t, vu)
    
!$OMP PARALLEL DO private(vel)
do x_c = 1, meshsize

x = x_array(x_c)

do y_c = 1, meshsize

y = y_array(y_c)

!vel = velocity_pointML(x_array(x_c),y_array(y_c),t,phase1,phase2,time1,time2,amplitudes,dispersions,g)

! vel = velocity_pointML(x,y,t,phase1,phase2,time1,time2,amplitudes,dispersions,g)
vel(1) = cos(tau*x+tau*y-t)
        
vel_array(x_c,y_c) = vel(1)






end do
end do
!$OMP END PARALLEL DO

write(1,*) vel_array


end do



end subroutine velocitycomp


end program velocitytest

