module passivetracers
use particle
use parameters
use rossby_wave
implicit none 

contains

function vel_passive(timesteps, phase1, phase2, time1, time2, N_part, g) result(partavg)

real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), parts(2)
    real(dp) :: t_array(timesteps), x_array(1000), y_array(1000), modx, mody
    real(dp) :: linx(N_part), liny(N_part), partx(N_part,N_part), party(N_part,N_part)
    real(dp) :: dispersions(65,65), amplitudes(65,65)
    real(dp) :: g, vel(2)
    real(dp), allocatable :: partavg(:)
    integer :: c_v1, c_v2
    real :: dt = 0.25
    integer :: counter_x, counter_y, counter_t, timesteps
    integer :: N_part, counter_n1x, c_n2x, counter_n1y, c_n2y
    integer :: beginning,end
    real :: rate
    integer :: beg1,end1
    real(dp) :: t0 = 0.
    integer :: counter_i
    

    allocate(partavg(timesteps+1))
    

    x_array = linspace(-5.0,5.0,1000)
    y_array = linspace(-5.0,5.0,1000)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
    liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

        print*, 'gamma= ', gamma    
    
    call system_clock(beginning, rate)
    open(2,file='partx.dat')
    open(3,file='party.dat')
    open(4,file='fixedpoint.dat')
    open(5,file='lengthvelocity.dat')



    do counter_n1x = 1, N_part
        do counter_n1y = 1, N_part
            partx(counter_n1x,counter_n1y) = linx(counter_n1x)
            party(counter_n1x,counter_n1y) = liny(counter_n1y)
        enddo
    enddo

    write(2,*) partx
    write(3,*) party
    
    partavg(1) = sum(partx)/(N_part**2)
    
    
    call dispersion_relation_array(dispersions)
    call amplitudes_array(amplitudes)
    
    counter_i = 0

    do c_v1 = 1, 1000
    do c_v2 = 1, 1000

        vel = velocity_pointML(x_array(c_v1),y_array(c_v2),t0,phase1,phase2,time1,time2,amplitudes,dispersions,g)
        if (vel(1) < 10e-15 .AND. vel(2) < 10e-15) then
            write(4,*) x_array(c_v1), y_array(c_v2)
            counter_i = counter_i + 1
        end if
    end do
    end do

    write(5,*) counter_i

    do counter_t = 1,timesteps
    
    t = t_array(counter_t)
    
    call MaduLawrence_loop(time1, phase1, t, vu)
    call MaduLawrence_loop(time2, phase2, t, vu)
    

    counter_i = 0 
    do c_v1 = 1, 1000
    do c_v2 = 1, 1000

        vel = velocity_pointML(x_array(c_v1),y_array(c_v2),t,phase1,phase2,time1,time2,amplitudes,dispersions,g)
        if (vel(1) < 10e-15 .AND. vel(2) < 10e-15) then
            write(4,*) x_array(c_v1), y_array(c_v2)
            counter_i = counter_i + 1
        end if
    end do
    end do
    
    write(5,*) counter_i
    



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
    write(2,*) partx
    write(3,*) party
    end do
   
    print*, N_part**2, ' particles: cluster big data'
    call system_clock(end)

    print *, "time : ", real(end-beginning) / real(rate)

end function vel_passive

end module passivetracers


