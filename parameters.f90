module parameters
use constants
        implicit none

real, parameter :: gamma = 1.0
real :: delta = 1.0
real(dp) :: vu = 10000*atan(1.)
real :: amp_scaling = 1.
integer :: time_steps = 500
integer :: N_particles = 200

end module parameters
