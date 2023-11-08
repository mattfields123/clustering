program inertialtracers
use particle
implicit none
real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)

phase1 = tau*random_matrix(65,65)
phase2 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)
!print*,time1(1,:)

!call vortex_inertial(1000000,3)
!call inertial_test(2000, phase1, phase2, time1, time2, 30, 0.0001_dp, '00010.csv')
!call inertial_test(1000, phase1, phase2, time1, time2, 30, 0.0002_dp, '00020.csv')
!call inertial_test(800, phase1, phase2, time1, time2, 30, 0.00025_dp, '00025.csv')
!call inertial_test(400, phase1, phase2, time1, time2, 30, 0.0005_dp, '00050.csv')
!call inertial_test(200, phase1, phase2, time1, time2, 30, 0.001_dp, '00100.csv')
!call inertial_test(100, phase1, phase2, time1, time2, 30, 0.002_dp, '00200.csv')
call super_inertial(1000, phase1, phase2, time1, time2, 200)

contains

function vortex_soulsby(x, y, u, v) result(soulsby)
    real(dp) :: x, y, u, v, soulsby(2)
    real :: R = 2.
    real :: St = 3800.

    soulsby = 1.5*R*vortex_acc(x,y)
    soulsby = soulsby + (vortex_vel(x,y) - (/u, v/))*St
end function vortex_soulsby

function inertial_vortex_soulsby(x, y, u, v, dt) result(posvel)
    real(dp) :: x, y, dt, u, v
    real(dp) :: k1(4), k2(4), k3(4), k4(4), posvel(4)

    posvel = (/x, y, u, v/)

    k1 = (/u, v, vortex_soulsby(x,y,u,v)/)
    k2 = (/u+dt*k1(3)*0.5, v+dt*k1(4)*0.5, vortex_soulsby(x+dt*k1(1)*0.5,y+dt*k1(2)*0.5,u+dt*k1(3)*0.5,v+dt*k1(4)*0.5)/)
    k3 = (/u+dt*k2(3)*0.5, v+dt*k2(4)*0.5, vortex_soulsby(x+dt*k2(1)*0.5,y+dt*k2(2)*0.5,u+dt*k2(3)*0.5,v+dt*k2(4)*0.5)/)
    k4 = (/u+dt*k3(3), v+dt*k3(4), vortex_soulsby(x+dt*k3(1),y+dt*k3(2),u+dt*k3(3),v+dt*k3(4))/)

    posvel = posvel + dt*(k1+2*k2+2*k3+k4)/6
end function inertial_vortex_soulsby

subroutine vortex_inertial(timesteps, N_part)
real(dp) :: x, y, t, velocity(2), t_array(timesteps), x_array(128), y_array(128)
real(dp) :: linx(N_part), liny(N_part), k1(2), k2(2), k3(2), k4(2), dt, posvel(4), vvel(2)
!real(dp) :: partx_array(N_part,N_part), party_array(N_part,N_part), partu_array(N_part,N_part), partv_array(N_part,N_part)
real(dp) :: partx_array, party_array, partu_array, partv_array
integer :: counter_x, counter_y, counter_t, timesteps, N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

x_array = linspace(-1.0,1.0,128)
y_array = linspace(-1.0,1.0,128)
t_array = linspace(0.,timesteps-1.,timesteps)

velocity = (/0.,0./)

linx = linspace(0.3,0.7,N_part)
liny = linspace(-1.0,1.0,N_part)
dt = 0.0005_dp

open(1, file = 'inertialvelocity.dat')
open(2, file = 'inertialparticles.dat')
!do counter_n1x = 1, N_part
    !do counter_n1y = 1, N_part
        !vvel = vortex_vel(linx(counter_n1x),liny(2))
        !partx_array(counter_n1x,:) = linx(counter_n1x)
        !party_array(counter_n1x,:) = liny(2)
        !partu_array(counter_n1x,:) = vvel(1)
        !partv_array(counter_n1x,:) = vvel(2)
        !write(2,*) partx_array(counter_n1x,2), ',', party_array(counter_n1x,2), ',', &
        !    &partu_array(counter_n1x,2), ',', partv_array(counter_n1x,2)
        vvel = vortex_vel(0.3_dp,0._dp)
        partx_array = 0.3_dp
        party_array = 0._dp
        partu_array = vvel(1)
        partv_array = vvel(2)
        write(2,*) partx_array, ',', party_array, ',', partu_array, ',', partv_array
    !enddo
!enddo
!print*,partx_array
!print*,party_array

do counter_y = 1,128
    y = y_array(counter_y)
    do counter_x = 1,128
    x = x_array(counter_x)
        velocity = vortex_vel(x,y)
        write(1,*) velocity(1), ',', velocity(2)
    end do
end do

do counter_t = 1,timesteps
!print*,counter_t
t = t_array(counter_t)
!do counter_n2x = 1, N_part
!do counter_n2y = 1, N_part
    !posvel = inertial_vortex_soulsby(partx_array(counter_n2x,2),party_array(counter_n2x,2),&
    !    &partu_array(counter_n2x,2),partv_array(counter_n2x,2),dt)
    !partx_array(counter_n2x,2) = posvel(1)
    !party_array(counter_n2x,2) = posvel(2)
    !partu_array(counter_n2x,2) = posvel(3)
    !partv_array(counter_n2x,2) = posvel(4)
    !write(2,*) posvel(1), ',', posvel(2), ',', posvel(3), ',', posvel(4)
    posvel = inertial_vortex_soulsby(partx_array,party_array,partu_array,partv_array,dt)
    partx_array = posvel(1)
    party_array = posvel(2)
    partu_array = posvel(3)
    partv_array = posvel(4)
    write(2,*) posvel(1), ',', posvel(2), ',', posvel(3), ',', posvel(4)
!end do
!end do
end do
close(2)
close(1) 
print*,'Test'
end subroutine vortex_inertial

subroutine inertial_test(timesteps, phase1, phase2, time1, time2, N_part, dt, str)
real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2)
real(dp) :: x_array(128), y_array(128), posvel(4), partu_array(N_part,N_part), partv_array(N_part,N_part)
real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
real :: dt
integer :: counter_x, counter_y, counter_t, timesteps
integer :: N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y
character(len = 9) :: str

x_array = linspace(-5.0,5.0,128)
y_array = linspace(-5.0,5.0,128)

linx = linspace(-3.0,3.0,N_part)
liny = linspace(-3.0,3.0,N_part)

!open(1, file = 'velocity.csv')
open(2, file = str)
do counter_n1x = 1, N_part
    do counter_n1y = 1, N_part
        velocity = velocity_pointML(linx(counter_n1x),liny(counter_n1y), 0._dp, phase1, phase2, time1, time2)
        partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
        party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
        partu_array(counter_n1x,counter_n1y) = velocity(1)
        partv_array(counter_n1x,counter_n1y) = velocity(2)
        !write(2,*) partx_array(counter_n1x,counter_n1y), ',', party_array(counter_n1x,counter_n1y)
    enddo
enddo

do counter_t = 1,timesteps
print*,dt, counter_t
!$OMP PARALLEL DO private(posvel)
do counter_n2x = 1, N_part
do counter_n2y = 1, N_part
    posvel = inertial_soulsby(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y),&
    & partu_array(counter_n2x,counter_n2y),partv_array(counter_n2x,counter_n2y),t,phase1,phase2,time1,time2,dt)
    partx_array(counter_n2x,counter_n2y) = modulo(posvel(1) + 5, 10.) - 5
    party_array(counter_n2x,counter_n2y) = modulo(posvel(2) + 5, 10.) - 5
    partu_array(counter_n2x,counter_n2y) = posvel(3)
    partv_array(counter_n2x,counter_n2y) = posvel(4)
end do
end do
!$OMP END PARALLEL DO
if (counter_t == timesteps) then
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part
        write(2,*) partx_array(counter_n2x,counter_n2y), ',', party_array(counter_n2x,counter_n2y)
    end do
    end do
endif
end do
close(2)
!print*,'128 grid points w/ ', N_part**2, ' particles'
print*,dt
print*,partx_array
print*,party_array
end subroutine inertial_test

subroutine super_inertial(timesteps, phase1, phase2, time1, time2, N_part)
real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2)
real(dp) :: t_array(timesteps), x_array(256), y_array(256), modx, mody, phi, psi, div, OkW
!real(dp) :: t_array(timesteps), x_array(128), y_array(128), modx, mody, phi, psi, div, OkW
real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
real(dp) :: partdiv, partpot, partu_array(N_part,N_part), partv_array(N_part,N_part), partOkw, posvel(4)
real :: dt = 0.001
integer :: counter_x, counter_y, counter_t, timesteps, Okw_int
integer :: N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

!x_array = linspace(-5.0,5.0,256)
!y_array = linspace(-5.0,5.0,256)
x_array = linspace(-5.0,5.0,128)
y_array = linspace(-5.0,5.0,128)
t_array = linspace(dt,timesteps*dt,timesteps)

linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

open(1, file = 'velocity.dat')
open(2, file = 'particles.dat')
open(3, file = 'partavg.dat')
!open(4, file = 'streamfunction.dat')
!open(7, file = 'potential.dat')
open(8, file = 'divergence.dat')
open(9, file = 'okuboweiss.dat')
open(10, file = 'okuboclass.dat')
!$OMP PARALLEL DO private(velocity)
do counter_n1x = 1, N_part
    do counter_n1y = 1, N_part
        velocity = velocity_pointML(linx(counter_n1x),liny(counter_n1y), 0._dp, phase1, phase2, time1, time2)
        partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
        party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
        partu_array(counter_n1x,counter_n1y) = velocity(1)
        partv_array(counter_n1x,counter_n1y) = velocity(2)
    enddo
enddo
!$OMP END PARALLEL DO
write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2), ',', &
    & sum(partu_array)/(N_part**2), ',', sum(partv_array)/(N_part**2)

do counter_t = 1,timesteps
t = t_array(counter_t)
call MaduLawrence_loop(time1, phase1, t)
call MaduLawrence_loop(time2, phase2, t)
print*,counter_t, 'Vel'
do counter_y = 1,256
!do counter_y = 1,128
y = y_array(counter_y)
do counter_x = 1,256
!do counter_x = 1,128
x = x_array(counter_x)
    if ((mod(counter_y, 2) == 0) .and. (mod(counter_x, 2) == 0)) then
        velocity = velocity_pointML(x, y, t, phase1, phase2, time1, time2)
        write(1,*) velocity(1), ',', velocity(2)
    end if
    !psi = streamfunction(x, y, t, phase1, time1)
    !phi = potential(x, y, t, phase1, phase2, time1, time2)
    div = velocity_divergence(x, y, t, phase1, phase2, time1, time2)
    OkW = Okubo_Weiss(x, y, t, phase1, phase2, time1, time2)
    Okw_int = okubo_class(div, Okw)
    !write(4,*) psi
    !write(7,*) phi
    write(8,*) div
    write(9,*) OkW
    write(10,*) OkW_int
end do
end do
print*,counter_t, 'Part'
!$OMP PARALLEL DO private(posvel)
do counter_n2x = 1, N_part
do counter_n2y = 1, N_part
    posvel = inertial_soulsby(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y),&
    & partu_array(counter_n2x,counter_n2y),partv_array(counter_n2x,counter_n2y),t,phase1,phase2,time1,time2,dt)
    partx_array(counter_n2x,counter_n2y) = modulo(posvel(1) + 5, 10.) - 5
    party_array(counter_n2x,counter_n2y) = modulo(posvel(2) + 5, 10.) - 5
    partu_array(counter_n2x,counter_n2y) = posvel(3)
    partv_array(counter_n2x,counter_n2y) = posvel(4)
    write(2,*) posvel(1),',',posvel(2),',',posvel(3),',',posvel(4)
end do
end do
!$OMP END PARALLEL DO
write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2), ',', &
    & sum(partu_array)/(N_part**2), ',', sum(partv_array)/(N_part**2)
end do
close(10)
close(9)
close(8)
!close(7)
!close(4)
close(3)
close(2)
close(1) 
print*, 'Inertial experiment with ', N_part**2, ' particles'
end subroutine super_inertial

end program inertialtracers