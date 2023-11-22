program passivetracers
use particle
use parameters
implicit none
real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), amplitudes(65,65), dispersions(65,65)
integer :: i

phase1 = tau*random_matrix(65,65)
phase2 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)

call dispersion_relation_array(dispersions)
call amplitudes_array(amplitudes)

call vel_passive(time_steps, phase1, phase2, time1, time2, N_particles,dispersions, amplitudes)
!print *, char(7)  
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of external_prog.exe was ", i

contains

subroutine vel_passive(timesteps, phase1, phase2, time1, time2, N_part,dispersions,amplitudes)
    real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), parts(2)
    real(dp) :: dispersions(65,65), amplitudes(65,65)
    real(dp) :: t_array(timesteps), x_array(128), y_array(128), modx, mody
    real(dp) :: linx(N_part), liny(N_part), partx_array(N_part,N_part), party_array(N_part,N_part)
    real :: dt = 0.25
    integer :: counter_x, counter_y, counter_t, timesteps
    integer :: N_part, c_n1x, c_n2x, c_n1y, c_n2y
    integer :: beginning,end
    real :: rate
    integer :: beg1,end1
    x_array = linspace(-5.0,5.0,128)
    y_array = linspace(-5.0,5.0,128)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
    liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

    open(1, file = 'velocity.dat')
    open(2, file = 'particles.dat')
    open(3, file = 'partavg.dat')
    call system_clock(beginning, rate)
    do c_n1x = 1, N_part
        do c_n1y = 1, N_part
            partx_array(c_n1x,c_n1y) = linx(c_n1x)
            party_array(c_n1x,c_n1y) = liny(c_n1y)
            modx = partx_array(c_n1x,c_n1y)
            mody = party_array(c_n1x,c_n1y)
            write(2,*) modx, ',', mody, ',', partx_array(c_n1x,c_n1y), ',', party_array(c_n1x,c_n1y)
        enddo
    enddo
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2)

    do counter_t = 1,timesteps
    t = t_array(counter_t)

        call MaduLawrence_loop(time1,phase1,t)
        call MaduLawrence_loop(time2,phase2,t)


    print*,counter_t, 'Vel'
   
    call system_clock(beg1,rate)

    do counter_y = 1,128
    y = y_array(counter_y)
    do counter_x = 1,128
    x = x_array(counter_x)
        velocity = velocity_pointML(x, y, t, phase1, phase2, time1, time2,dispersions,amplitudes)
        write(1,*) velocity(1), ',', velocity(2)
    end do
    end do
    
    call system_clock(end1)
    print*, 'time to vel : ', real(end1-beg1) / real(rate)

    print*,counter_t, 'Part'
    !$OMP PARALLEL DO private(parts,modx,mody)
    do c_n2x = 1, N_part
    do c_n2y = 1, N_part
    parts = passive(partx_array(c_n2x,c_n2y),party_array(c_n2x,c_n2y),t,phase1,phase2,time1,time2,dt,dispersions,amplitudes)
        partx_array(c_n2x,c_n2y) = parts(1)
        party_array(c_n2x,c_n2y) = parts(2)
        !modx = modulo(parts(1) + 5, 10.) - 5
        !mody = modulo(parts(2) + 5, 10.) - 5
        !write(2,*) modx, ',', mody, ',', parts(1), ',', parts(2)
    end do
    end do
    !$OMP END PARALLEL DO
    do c_n2x = 1, N_part
    do c_n2y = 1, N_part        
        modx = modulo(partx_array(c_n2x,c_n2y) + 5, 10.) - 5
        mody = modulo(party_array(c_n2x,c_n2y) + 5, 10.) - 5
        write(2,*) modx, ',', mody, ',', partx_array(c_n2x,c_n2y), ',', party_array(c_n2x,c_n2y)
    end do
    end do
    write(3,*) sum(partx_array)/(N_part**2), ',', sum(party_array)/(N_part**2)
    end do
    close(3)
    close(2)
    close(1) 
    print*, N_part**2, ' particles: cluster big data'
    call system_clock(end)

    print *, "time : ", real(end-beginning) / real(rate)

end subroutine vel_passive

end program passivetracers
