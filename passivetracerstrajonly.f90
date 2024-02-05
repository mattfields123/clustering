module passivetracers
use particle
use parameters
use rossby_wave
implicit none 

contains

function vel_passive(timesteps, phase1, phase2, time1, time2, N_part, g) result(partavg)

real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), parts(2)
    real(dp) :: t_array(timesteps), x_array(128), y_array(128), modx, mody
    real(dp) :: linx(N_part), liny(N_part), partx(N_part,N_part), party(N_part,N_part)
    real(dp) :: dispersions(65,65), amplitudes(65,65)
    real(dp) :: g
    real(dp), allocatable :: partavg(:)

    real :: dt = 0.05
    integer :: counter_x, counter_y, counter_t, timesteps
    integer :: N_part, counter_n1x, c_n2x, counter_n1y, c_n2y
    integer :: beginning,end
    real :: rate
    integer :: beg1,end1

    allocate(partavg(timesteps+1))

    x_array = linspace(-5.0,5.0,128)
    y_array = linspace(-5.0,5.0,128)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
    liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

    
    call system_clock(beginning, rate)
    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx(counter_n1x,counter_n1y) = linx(counter_n1x)
            party(counter_n1x,counter_n1y) = liny(counter_n1y)
        enddo
    enddo

    
    
    partavg(1) = sum(partx)/(N_part**2)
    
    
    call dispersion_relation_array(dispersions)
    call amplitudes_array(amplitudes)
    
    do counter_t = 1,timesteps
    
    t = t_array(counter_t)
    
    call MaduLawrence_loop(time1, phase1, t, vu)
    call MaduLawrence_loop(time2, phase2, t, vu)
    
    !$OMP PARALLEL DO private(parts,modx,mody)
    do c_n2x = 1, N_part
        do c_n2y = 1, N_part

            parts = passive(partx(c_n2x,c_n2y),party(c_n2x,c_n2y),t,phase1,phase2,time1,time2,dt,amplitudes,dispersions,g)
            partx(c_n2x,c_n2y) = parts(1)
            party(c_n2x,c_n2y) = parts(2)

        end do
    end do
    !$OMP END PARALLEL DO 
    partavg(counter_t+1) = sum(partx)/(N_part**2)

    end do
   
    print*, N_part**2, ' particles: cluster big data'
    call system_clock(end)

    print *, "time : ", real(end-beginning) / real(rate)

end function vel_passive

end module passivetracers


