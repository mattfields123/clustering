module parameters
use constants
        implicit none


real, parameter :: gamma = 0.15
real :: delta = 0.
real(dp) :: vu = 100*atan(1.)
real :: amp_scaling = 1
integer :: time_steps = 1000
integer :: N_particles = 50

end module parameters
