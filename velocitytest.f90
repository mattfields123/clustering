program velocitytest
use constants
use parameters
use admin
use rossby_wave
use velocity

integer :: tsteps, meshsize
real(dp) :: thresh
meshsize = 500
tsteps = 5
thresh = 0.1


call velocitycomp(meshsize,tsteps,thresh)








contains 

subroutine velocitycomp(meshsize,tsteps,thresh)
integer :: done
integer :: tsteps, meshsize
real(dp), allocatable :: t_array(:)
real(dp), allocatable :: x_array(:), y_array(:)
real(dp) :: vel(2), t, x, y, time1(65,65), time2(65,65), phase1(65,65), phase2(65,65), dispersions(65,65), amplitudes(65,65)
integer :: t_c, x_c, y_c
real(dp) :: g=0.
real :: dt = 0.01
real(dp), allocatable :: vel_array(:,:)
real(dp), allocatable :: fixed_array(:,:)
real(dp) :: thresh

allocate(t_array(tsteps))
allocate(x_array(meshsize))
allocate(y_array(meshsize))
allocate(vel_array(meshsize,meshsize))
allocate(fixed_array(meshsize,meshsize))

t_array = linspace(0.,(tsteps-1)*dt,tsteps)
x_array = linspace(-2.,2.,meshsize)
y_array = linspace(-2.,2.,meshsize)


call dispersion_relation_array(dispersions)
call amplitudes_array(amplitudes)
    

open(1,file='velocities.dat')
open(2,file='fixed.dat')

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
vel(1) = sin(tau*x+tau*y-t)
       
vel_array(x_c,y_c) = vel(1)



IF (abs(vel(1)) < thresh) THEN
fixed_array(x_c,y_c) = 1.
ELSE
fixed_array(x_c,y_c) = 0.
END IF




end do
end do
!$OMP END PARALLEL DO

write(1,*) vel_array

write(2,*) fixed_array
end do



end subroutine velocitycomp


end program velocitytest

