module parameters
use constants
        implicit none

real, parameter :: gamma = 0.1
real :: delta = 0.1
real(dp) :: vu = 10000*atan(1.)
real :: amp_scaling = 1.
integer :: time_steps = 5
integer :: N_particles = 2

end module parameters
