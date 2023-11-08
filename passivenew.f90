program passivenew
        use parameters
        use particle
        implicit none
        real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
        integer :: i
        
        ! Initialisation of phases and lifetimes

        phase1 = tau*random_matrix(65,65)
        phase2 = tau*random_matrix(65,65)
        time1 = vu*random_matrix(65,65)
        time2 = vu*random_matrix(65,65)

        ! Runs program

        call tracers(time_steps,phase1,phase2,time1,time2,N_particles)



        ! Subroutines
        contains

                subroutine tracers(timesteps, phase1, phase2, time1, time2, N_part)

                        real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65), parts(2), partsalpha(2)
                        real(dp) :: t_array(timesteps)
                        real(dp) :: linx(N_part), liny(N_part)
                        ! real(dp) :: particledata(2,timesteps,N_part)
                        real(dp), allocatable :: particledata(:,:,:,:) 
                        real :: dt = 0.25
                        integer :: counter_x,counter_y,counter_t,timesteps
                        integer :: N_part
                        integer :: start, end
                        real :: rate
                        
                        t_array = linspace(0., (timesteps-1)*dt,timesteps)

                        linx = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
                   
                        liny = linspace(-5.0,(N_part-1.)/(0.1*N_part)-5.0,N_part)
                        allocate(particledata(2,timesteps,N_part,N_part))                       
print*, timesteps, N_part




                        call system_clock(start,rate)
                        

                        !$OMP PARALLEL DO private(particledata)
                
                        do counter_x = 1, N_part
                        do counter_y = 1, N_part 
                                
                                particledata(1,1,counter_x,counter_y) = linx(counter_x)
                                particledata(2,1,counter_x,counter_y) = liny(counter_y)
                                
                                    do counter_t = 2,timesteps
                                        
                                        partsalpha(1) = particledata(1,counter_t-1,counter_x,counter_y)
                                        partsalpha(2) = particledata(2,counter_t-1,counter_x,counter_y)
                                        parts = passive(partsalpha(1),partsalpha(2),t_array(counter_t),phase1,phase2,time1,time2,dt)
                                        particledata(1,counter_t,counter_x,counter_y) = parts(1)
                                        particledata(2,counter_t,counter_x,counter_y) = parts(2)
                                        
                                    end do
                                    end do
                                    end do
                        !$OMP END PARALLEL DO
                        call system_clock(end)
                        print*,'successful', real(end-start) / real(rate)
                                        
                end subroutine tracers
                                        
                        

        end program passivenew

