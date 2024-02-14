module parameters
use constants
        implicit none

real, parameter :: gamma = 1.
real :: delta = 1.
real(dp) :: vu = 100000*atan(1.)
real :: amp_scaling = 1.
integer :: time_steps = 10000
integer :: N_particles = 200

end module parameters
