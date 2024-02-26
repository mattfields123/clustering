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
    integer :: vel_domain = 1000
    integer :: counter_i=0
    
    real(dp), allocatable :: fixedpointlocation(:,:,:)
    integer, allocatable :: fixedpointnumbers(:) 
    real(dp), allocatable :: uvel(:,:)
    real(dp), allocatable :: vvel(:,:)
    allocate(partavg(timesteps+1))
    allocate(x_array(vel_domain))
    allocate(y_array(vel_domain))
    allocate(fixed(vel_domain,vel_domain))
    allocate(fixedpointlocation(timesteps,5000,2))
    allocate(fixedpointnumbers(timesteps))
    allocate(uvel(vel_domain,vel_domain))
    allocate(vvel(vel_domain,vel_domain))

    x_array = linspace(-5.0,5.0,vel_domain)
    y_array = linspace(-5.0,5.0,vel_domain)
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
    open(4,file='fixedpoint.dat')
    open(5,file='counters.dat')
    open(6,file='uvel.dat')
    open(7,file='vvel.dat')

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
    
    t = t_array(counter_t)
    
    call MaduLawrence_loop(time1, phase1, t, vu)
    call MaduLawrence_loop(time2, phase2, t, vu)
    
    counter_i = 0

    !$OMP PARALLEL DO private(c_v1,c_v2,vel) reduction(+:counter_i)
    do c_v1 = 1, vel_domain
    do c_v2 = 1, vel_domain

        vel = velocity_pointML(x_array(c_v1),y_array(c_v2),t,phase1,phase2,time1,time2,amplitudes,dispersions,g)
        
        IF ((abs(vel(1)) < threshold) .AND. (abs(vel(2)) < threshold)) then
            counter_i = counter_i + 1 
            fixedpointlocation(counter_t,counter_i,1) = x_array(c_v1)
            fixedpointlocation(counter_t,counter_i,2) = y_array(c_v2)
        END IF
        uvel = vel(1)
        vvel = vel(2)

    end do
    end do
    !$OMP END PARALLEL DO
    
    fixedpointnumbers(counter_t+1) = counter_i

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
    write(6,*) uvel
    write(7,*) vvel
    
        print*, 'timestep : ', counter_t
        
    end do


    write(4,*) fixedpointlocation
    write(5,*) fixedpointnumbers

    print*, N_part**2, ' particles: cluster big data'
    call system_clock(end)

    print *, "time : ", real(end-beginning) / real(rate)

end function vel_passive

end module passivetracers


