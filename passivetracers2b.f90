program passivetracers2
use velocity2
implicit none
real(dp) :: phase1(65,65), time1(65,65)
integer :: i

phase1 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)

!call basin_speed2(500, phase1, time1)
!call passive_snapshot(3000, phase1, time1, 500)
call super_passive2(3000, phase1, time1, 500)
!print *, char(7)
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of email.ps1 was ", i

contains

subroutine super_passive2(timesteps, phase1, time1, N_part)
    real(dp) :: x, y, t, phase1(65,65), time1(65,65), velocity(2), parts(2), psi_array(256,256), div_array(256,256)
    real(dp) :: t_array(timesteps), x_array(256), y_array(256), psi, div, zeta, modx, mody
    !real(dp) :: t_array(timesteps), x_array(128), y_array(128), psi, div, zeta
    real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
    real :: dt = 0.025
    integer :: counter_x, counter_y, counter_t, timesteps, N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

    x_array = linspace(-0.5,0.5,256)
    y_array = linspace(-0.5,0.5,256)
    !x_array = linspace(-0.5,0.5,128)
    !y_array = linspace(-0.5,0.5,128)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-0.5,(N_part-1.)/(N_part)-0.5,N_part)
    liny = linspace(-0.5,(N_part-1.)/(N_part)-0.5,N_part)

    open(1, file = 'betavelocity2.dat')
    open(2, file = 'betaparticles2.dat')
    open(3, file = 'betapartavg2.dat')
    open(4, file = 'betastreamfunction2.dat')
    open(7, file = 'betadivergence2.dat')
    !open(8, file = 'betavorticity2.dat')

    print*,'0 Part'
    !$OMP PARALLEL DO private(parts)
    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
            party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
            parts = (/linx(counter_n1x),liny(counter_n1y)/)
            write(2,*) parts(1),',',parts(2),',',parts(1),',',parts(2)
        enddo
    enddo
    !$OMP END PARALLEL DO
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2), ',', sum(psi_array)/sum(div_array)

    do counter_t = 1,timesteps
    t = t_array(counter_t)
    call MaduLawrence_loop2(time1, phase1, t)
    print*,counter_t, 'Vel'
    do counter_y = 1,256
    !do counter_y = 1,128
    y = y_array(counter_y)
    do counter_x = 1,256
    !do counter_x = 1,128
    x = x_array(counter_x)
        if ((mod(counter_y, 2) == 0) .and. (mod(counter_x, 2) == 0)) then
            velocity = velocity_pointML2(x, y, t, phase1, time1)
            write(1,*) velocity(1), ',', velocity(2)
        end if
        psi = streamfunction2(x, y, t, phase1, time1)
        psi_array(counter_x,counter_y) = psi
        div = velocity_divergence2(x, y, t, phase1, time1)
        div_array(counter_x,counter_y) = div
        !zeta = vorticity2(x, y, t, phase1, time1)
        write(4,*) psi
        write(7,*) div
        !write(8,*) zeta
    end do
    end do
    print*,counter_t, 'Part'
    !$OMP PARALLEL DO private(parts)
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part
    parts = passive2(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y),t,phase1,time1,dt)
        partx_array(counter_n2x,counter_n2y) = parts(1)
        party_array(counter_n2x,counter_n2y) = parts(2)
        modx = modulo(parts(1) + 0.5, 1.) - 0.5
        mody = modulo(parts(2) + 0.5, 1.) - 0.5
        write(2,*) modx,',',mody,',',parts(1),',',parts(2)
    end do
    end do
    !$OMP END PARALLEL DO
        write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2), ',', sum(psi_array)/sum(div_array)
    end do
    !close(8)
    close(7)
    close(4)
    close(3)
    close(2)
    close(1) 
    print*, 'Beta plane experiment ', N_part**2, ' particles on a medium domain'
end subroutine super_passive2

subroutine basin_speed2(timesteps, phase1, time1)
    real(dp) :: x, y, t, spd, phase1(65,65), time1(65,65)
    real(dp) :: t_array(timesteps), x_array(256), y_array(256)
    integer :: counter_x, counter_y, counter_t, timesteps

    x_array = linspace(-0.5,0.5,256)
    y_array = linspace(-0.5,0.5,256)
    t_array = linspace(0.,timesteps-1.,timesteps)

    do counter_t = 1,timesteps
    print*,counter_t
    t = t_array(counter_t)
    call MaduLawrence_loop2(time1, phase1, t)
    do counter_y = 1,256
    y = y_array(counter_y)
    do counter_x = 1,256
    x = x_array(counter_x)
        spd = spd + speed_point2(x, y, t, phase1, time1)
    end do
    end do
    end do
    print*,spd/(256*256*timesteps)
end subroutine basin_speed2

subroutine passive_snapshot(timesteps, phase1, time1, N_part)
    real(dp) :: x, y, t, phase1(65,65), time1(65,65), velocity(2), parts(2), modx, mody
    real(dp) :: t_array(timesteps), x_array(256), y_array(256), psi, div, zeta
    real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
    real :: dt = 0.025
    integer :: counter_x, counter_y, counter_t, timesteps, N_part, counter_n1x, counter_n2x, counter_n1y, counter_n2y

    x_array = linspace(-0.5,0.5,256)
    y_array = linspace(-0.5,0.5,256)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-0.5,(N_part-1.)/(N_part)-0.5,N_part)
    liny = linspace(-0.5,(N_part-1.)/(N_part)-0.5,N_part)
    t = 0._dp

    open(1, file = 'betavelocitysnap2.dat')
    open(2, file = 'betaparticlessnap2.dat')
    open(3, file = 'betapartavgsnap2.dat')
    open(4, file = 'betastreamfunctionsnap2.dat')
    open(7, file = 'betadivergencesnap2.dat')
    !open(8, file = 'betavorticitysnap2.dat')
    
    print*,'0 Part'
    !$OMP PARALLEL DO private(parts)
    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx_array(counter_n1x,counter_n1y) = linx(counter_n1x)
            party_array(counter_n1x,counter_n1y) = liny(counter_n1y)
            parts = (/linx(counter_n1x),liny(counter_n1y)/)
            write(2,*) parts(1),',',parts(2),',',parts(1),',',parts(2)
        enddo
    enddo
    !$OMP END PARALLEL DO
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2)

    print*,'0 Vel'
    do counter_y = 1,256
    y = y_array(counter_y)
    do counter_x = 1,256
    x = x_array(counter_x)
        if ((mod(counter_y, 2) == 0) .and. (mod(counter_x, 2) == 0)) then
            velocity = velocity_pointML2(x, y, t, phase1, time1)
                write(1,*) velocity(1), ',', velocity(2)
        end if
        psi = streamfunction2(x, y, t, phase1, time1)
        div = velocity_divergence2(x, y, t, phase1, time1)
        !zeta = vorticity2(x, y, t, phase1, time1)
        write(4,*) psi
        write(7,*) div
        !write(8,*) zeta
    end do
    end do

    do counter_t = 1,timesteps
    print*,counter_t, 'Part'
    !$OMP PARALLEL DO private(parts)
    do counter_n2x = 1, N_part
    do counter_n2y = 1, N_part
    parts = passive2(partx_array(counter_n2x,counter_n2y),party_array(counter_n2x,counter_n2y),t,phase1,time1,dt)
        partx_array(counter_n2x,counter_n2y) = parts(1)
        party_array(counter_n2x,counter_n2y) = parts(2)
        modx = modulo(parts(1) + 0.5, 1.) - 0.5
        mody = modulo(parts(2) + 0.5, 1.) - 0.5
        write(2,*) modx,',',mody,',',parts(1),',',parts(2)
    end do
    end do
    !$OMP END PARALLEL DO
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2)
    end do
    !close(8)
    close(7)
    close(4)
    close(3)
    close(2)
    close(1) 
    print*, 'Beta snapshot with ', N_part**2, ' particles on a medium scale domain'
end subroutine passive_snapshot


end program passivetracers2