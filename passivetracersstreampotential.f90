module passivetracers
use particle
use parameters
use rossby_wave
implicit none 

contains

function vel_passive(timesteps, phase1, phase2, time1, time2, N_part, g) result(partavg)

real(dp) :: x, y, t, phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), velocity(2), parts(2)
    real(dp) :: t_array(timesteps), modx, mody
    real(dp) :: linx(N_part), liny(N_part), partx(N_part,N_part), party(N_part,N_part)
    real(dp) :: dispersions(65,65), amplitudes(65,65)
    real(dp) :: g, vel(2)
    real(dp) :: threshold = 0.01

    real(dp), allocatable :: partavg(:)
    real(dp), allocatable :: x_array(:)
    real(dp), allocatable :: y_array(:)
    integer, allocatable :: fixed(:,:)
    integer :: c_v1, c_v2
    real :: dt = 0.25
    integer :: counter_x, counter_y, counter_t, timesteps
    integer :: N_part, counter_n1x, c_n2x, counter_n1y, c_n2y
    integer :: beginning,end
    real :: rate
    integer :: beg1,end1
    integer :: vel_domain = 100
    integer :: counter_i=0
    


    real(dp), allocatable :: fixedpointlocation(:,:,:)
    integer, allocatable :: fixedpointnumbers(:) 
    real(dp), allocatable :: stream(:,:)
    real(dp), allocatable :: pot(:,:)
    integer, allocatable :: bothvel(:,:)

    allocate(partavg(timesteps+1))
    allocate(x_array(vel_domain))
    allocate(y_array(vel_domain))
    allocate(fixed(vel_domain,vel_domain))
    allocate(fixedpointlocation(timesteps,8000,2))
    allocate(fixedpointnumbers(timesteps))
    allocate(stream(vel_domain,vel_domain))
    allocate(pot(vel_domain,vel_domain))
    allocate(bothvel(vel_domain,vel_domain))    


    x_array = linspace(-5.0,(vel_domain-1.)/(0.1*vel_domain)-5.0,vel_domain)
    y_array = linspace(-5.0,(vel_domain-1.)/(0.1*vel_domain)-5.0,vel_domain)
    t_array = linspace(0.,(timesteps-1)*dt,timesteps)
 
    linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
    liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)

        print*, 'gamma= ', gamma    
        print*, 'tsteps =', timesteps
        print*, 'N_part = ', N_part
        print*, 'threshold = ',threshold


 
    call system_clock(beginning, rate)
    open(2,file='partx.dat')
    open(3,file='party.dat')
    open(8,file='stream.dat')
    open(9,file='pot.dat')
      


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

    do counter_t = 1,timesteps
   
        print*, counter_t
 
    t = t_array(counter_t)
    
    call MaduLawrence_loop(time1, phase1, t, vu)
    call MaduLawrence_loop(time2, phase2, t, vu)
    
    !$OMP PARALLEL DO private(c_v1,c_v2,vel)
    do c_v1 = 1, vel_domain
    do c_v2 = 1, vel_domain

        stream(c_v1,c_v2) = streamfunction(x_array(c_v1),y_array(c_v2),t,phase1,phase2,time1,time2,amplitudes,dispersions,g)
        pot(c_v1,c_v2) = potentialfunction(x_array(c_v1),y_array(c_v2),t,phase1,phase2,time1,time2,amplitudes,dispersions,g)
    end do
    end do
    !$OMP END PARALLEL DO
    
 
    !$OMP PARALLEL DO private(parts,modx,mody,c_n2x,c_n2y)
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
    write(8,*) stream
    write(9,*) pot

    print*, 'timestep : ', counter_t
        
    end do

    print*, N_part**2, ' particles: cluster big data'
    call system_clock(end)

    print *, "time : ", real(end-beginning) / real(rate)

end function vel_passive

end module passivetracers


