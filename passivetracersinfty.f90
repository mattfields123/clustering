program passivetracers
use particle
implicit none
real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
integer :: i

phase1 = tau*random_matrix(65,65)
phase2 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)

call super_passive(2000, phase1, phase2, time1, time2, 500)
!print *, char(7)  
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of external_prog.exe was ", i

contains

subroutine vel_passive(timesteps, phase1, phase2, time1, time2, N_part)
    real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), parts(2)
    real(dp) :: t_array(timesteps), x_array(128), y_array(128), modx, mody
    real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
    real :: dt = 0.25
    integer :: counter_x, counter_y, counter_t, timesteps
    integer :: N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

    x_array = linspace(-5.0,5.0,128)
    y_array = linspace(-5.0,5.0,128)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
    liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

    open(1, file = 'velocity.dat')
    open(2, file = 'particles.dat')
    open(3, file = 'partavg.dat')
    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
            party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
            modx = partx_array(counter_n1x,counter_n1y)
            mody = party_array(counter_n1x,counter_n1y)
            write(2,*) modx, ',', mody, ',', partx_array(counter_n1x,counter_n1y), ',', party_array(counter_n1x,counter_n1y)
        enddo
    enddo
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2)

    do counter_t = 1,timesteps
    t = t_array(counter_t)
    print*,counter_t, 'Vel'
    do counter_y = 1,128
    y = y_array(counter_y)
    do counter_x = 1,128
    x = x_array(counter_x)
        velocity = velocity_pointML(x, y, t, phase1, phase2, time1, time2)
        write(1,*) velocity(1), ',', velocity(2)
    end do
    end do
    print*,counter_t, 'Part'
    !$OMP PARALLEL DO private(parts,modx,mody)
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part
    parts = passive(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y),t,phase1,phase2,time1,time2,dt)
        partx_array(counter_n2x,counter_n2y) = parts(1)
        party_array(counter_n2x,counter_n2y) = parts(2)
        !modx = modulo(parts(1) + 5, 10.) - 5
        !mody = modulo(parts(2) + 5, 10.) - 5
        !write(2,*) modx, ',', mody, ',', parts(1), ',', parts(2)
    end do
    end do
    !$OMP END PARALLEL DO
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part        
        modx = modulo(partx_array(counter_n2x,counter_n2y) + 5, 10.) - 5
        mody = modulo(party_array(counter_n2x,counter_n2y) + 5, 10.) - 5
        write(2,*) modx, ',', mody, ',', partx_array(counter_n2x,counter_n2y), ',', party_array(counter_n2x,counter_n2y)
    end do
    end do
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2)
    end do
    close(3)
    close(2)
    close(1) 
    print*, N_part**2, ' particles: cluster big data'
end subroutine vel_passive

subroutine vel_passive2(timesteps, N_part)
real(dp) :: x, y, t, velocity(2), t_array(timesteps), x_array(128), y_array(128)
real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part), k1(2), k2(2), k3(2), k4(2), dt
integer :: counter_x, counter_y, counter_t, timesteps, N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

x_array = linspace(-5.0,5.0,128)
y_array = linspace(-5.0,5.0,128)
t_array = linspace(0.,timesteps-1.,timesteps)

velocity = (/0.,0./)

linx = linspace(-3.0,3.0,N_part)
liny = linspace(-3.0,3.0,N_part)
dt = 0.5

open(1, file = 'velocity.csv')
open(2, file = 'particles.csv')
do counter_n1x = 1, N_part
    do counter_n1y = 1, N_part
        partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
        party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
        write(2,*) partx_array(counter_n1x,counter_n1y), ',', party_array(counter_n1x,counter_n1y)
    enddo
enddo
print*,partx_array
print*,party_array

do counter_y = 1,128
    y = y_array(counter_y)
    do counter_x = 1,128
    x = x_array(counter_x)
        velocity = vortex_vel(x,y)
        write(1,*) velocity(1), ',', velocity(2)
    end do
end do

do counter_t = 1,timesteps
print*,counter_t
t = t_array(counter_t)
do counter_n2x = 1, N_part
do counter_n2y = 1, N_part
    k1 = vortex_vel(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y))
    k2 = vortex_vel(partx_array(counter_n2x,counter_n2y)+k1(1)*dt/2,party_array(counter_n2x,counter_n2y)+k1(2)*dt/2)
    k3 = vortex_vel(partx_array(counter_n2x,counter_n2y)+k2(1)*dt/2,party_array(counter_n2x,counter_n2y)+k2(2)*dt/2)
    k4 = vortex_vel(partx_array(counter_n2x,counter_n2y)+k3(1)*dt,party_array(counter_n2x,counter_n2y)+k3(2)*dt)
    partx_array(counter_n2x,counter_n2y) = partx_array(counter_n2x,counter_n2y) + dt*(k1(1)+2*k2(1)+2*k3(1)+k4(1))/6
    party_array(counter_n2x,counter_n2y) = party_array(counter_n2x,counter_n2y) + dt*(k1(2)+2*k2(2)+2*k3(2)+k4(2))/6
    write(2,*) partx_array(counter_n2x,counter_n2y), ',', party_array(counter_n2x,counter_n2y)
end do
end do
end do
close(2)
close(1) 
print*,'Test'
end subroutine vel_passive2

subroutine vel_passive4(timesteps, phase1, phase2, time1, time2, N_part, dt, str)
    real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2)
    real(dp) :: x_array(128), y_array(128), k1(2), k2(2), k3(2), k4(2), dt
    real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
    integer :: counter_x, counter_y, counter_t, timesteps
    integer :: N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y
    character(len = 11) :: str

    x_array = linspace(-5.0,5.0,128)
    y_array = linspace(-5.0,5.0,128)
 
    linx = linspace(-3.0,3.0,N_part)
    liny = linspace(-3.0,3.0,N_part)

    !open(1, file = 'velocity.csv')
    open(2, file = str)
    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
            party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
            !write(2,*) partx_array(counter_n1x,counter_n1y), ',', party_array(counter_n1x,counter_n1y)
        enddo
    enddo

    !do counter_y = 1,128
    !y = y_array(counter_y)
    !do counter_x = 1,128
    !x = x_array(counter_x)
    !    velocity = velocity_pointML(x, y, 0._dp, phase1, phase2, time1, time2)
    !    write(1,*) velocity(1), ',', velocity(2)
    !end do
    !end do
    !close(1) 

    do counter_t = 1,timesteps
    print*,dt, counter_t
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part
        k1 = velocity_pointML(partx_array(counter_n2x,counter_n2y),&
            &party_array(counter_n2x,counter_n2y),0._dp,phase1,phase2,time1,time2)
        k2 = velocity_pointML(partx_array(counter_n2x,counter_n2y)+k1(1)*dt/2,&
            &party_array(counter_n2x,counter_n2y)+k1(2)*dt/2,0._dp,phase1,phase2,time1,time2)
        k3 = velocity_pointML(partx_array(counter_n2x,counter_n2y)+k2(1)*dt/2,&
            &party_array(counter_n2x,counter_n2y)+k2(2)*dt/2,0._dp,phase1,phase2,time1,time2)
        k4 = velocity_pointML(partx_array(counter_n2x,counter_n2y)+k3(1)*dt,&
            &party_array(counter_n2x,counter_n2y)+k3(2)*dt,0._dp,phase1,phase2,time1,time2)
        partx_array(counter_n2x,counter_n2y) = modulo(partx_array(counter_n2x,counter_n2y) +&
            & 5 + dt*(k1(1)+2*k2(1)+2*k3(1)+k4(1))/6, 10.) - 5
        party_array(counter_n2x,counter_n2y) = modulo(party_array(counter_n2x,counter_n2y) + &
            &5 + dt*(k1(2)+2*k2(2)+2*k3(2)+k4(2))/6, 10.) - 5
        if (counter_t == timesteps) then
            write(2,*) partx_array(counter_n2x,counter_n2y), ',', party_array(counter_n2x,counter_n2y)
        endif
    end do
    end do
    end do
    close(2)
    !print*,'128 grid points w/ ', N_part**2, ' particles'
    print*,dt
    print*,partx_array
    print*,party_array
end subroutine vel_passive4

subroutine super_passive(timesteps, phase1, phase2, time1, time2, N_part)
    real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), parts(2)
    real(dp) :: t_array(timesteps), x_array(256), y_array(256), modx, mody, phi, psi, div, OkW
    !real(dp) :: t_array(timesteps), x_array(128), y_array(128), modx, mody, phi, psi, div, OkW
    real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
    real(dp) :: partdiv, partpot, velx_array(N_part,N_part), vely_array(N_part,N_part), partOkw, partvel(2)
    real :: dt = 0.25
    integer :: counter_x, counter_y, counter_t, timesteps, Okw_int
    integer :: N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

    x_array = linspace(-5.0,5.0,256)
    y_array = linspace(-5.0,5.0,256)
    !x_array = linspace(-5.0,5.0,128)
    !y_array = linspace(-5.0,5.0,128)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
    liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

    open(1, file = 'velocity.dat')
    open(2, file = 'particles.dat')
    open(3, file = 'partavg.dat')
    open(4, file = 'streamfunction.dat')
    open(7, file = 'potential.dat')
    open(8, file = 'divergence.dat')
    !open(9, file = 'betaokuboweiss.dat')
    !open(10, file = 'betaokuboclass.dat')
    !$OMP PARALLEL DO private(parts,modx,mody,partvel)
    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
            party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
            parts = (/linx(counter_n1x),liny(counter_n1y)/)
            modx = partx_array(counter_n1x,counter_n1y)
            mody = party_array(counter_n1x,counter_n1y)
            partvel = velocity_pointML(modx, mody, t, phase1, phase2, time1, time2)
            velx_array(counter_n1x,counter_n1y) = partvel(1)
            vely_array(counter_n1x,counter_n1y) = partvel(2)
            write(2,*) modx,',',mody,',',parts(1),',',parts(2),',',partvel(1),',',partvel(2)
        enddo
    enddo
    !$OMP END PARALLEL DO
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2), ',', &
        & sum(velx_array)/(N_part**2), ',', sum(vely_array)/(N_part**2)
    !print*, partx_array(:,1)
    !print*, velx_array(:,1)

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
        psi = streamfunction(x, y, t, phase1, time1)
        phi = potential(x, y, t, phase1, phase2, time1, time2)
        div = velocity_divergence(x, y, t, phase1, phase2, time1, time2)
        !OkW = Okubo_Weiss(x, y, t, phase1, phase2, time1, time2)
        !Okw_int = okubo_class(div, Okw)
        write(4,*) psi
        write(7,*) phi
        write(8,*) div
        !write(9,*) OkW
        !write(10,*) OkW_int
    end do
    end do
    print*,counter_t, 'Part'
    !$OMP PARALLEL DO private(parts,modx,mody,partvel)
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part
    parts = passive(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y),t,phase1,phase2,time1,time2,dt)
        partx_array(counter_n2x,counter_n2y) = parts(1)
        party_array(counter_n2x,counter_n2y) = parts(2)
        modx = modulo(parts(1) + 5, 10.) - 5
        mody = modulo(parts(2) + 5, 10.) - 5
        partvel = velocity_pointML(modx, mody, t, phase1, phase2, time1, time2)
        velx_array(counter_n2x,counter_n2y) = partvel(1)
        vely_array(counter_n2x,counter_n2y) = partvel(2)
        write(2,*) modx,',',mody,',',parts(1),',',parts(2),',',partvel(1),',',partvel(2)
    end do
    end do
    !$OMP END PARALLEL DO
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2), ',', &
        & sum(velx_array)/(N_part**2), ',', sum(vely_array)/(N_part**2)
    end do
    !close(10)
    !close(9)
    close(8)
    close(7)
    close(4)
    close(3)
    close(2)
    close(1) 
    print*, 'Super-sized experiment with ', N_part**2, ' particles'
end subroutine super_passive

end program passivetracers