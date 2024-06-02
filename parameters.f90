module parameters
use constants
        implicit none

real, parameter :: gamma = 0.1
real :: delta = 0.1
real(dp) :: vu = 10*atan(1.)
real :: amp_scaling = 1.
integer :: time_steps = 500
integer :: N_particles = 500

end module parameters
