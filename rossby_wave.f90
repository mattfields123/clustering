module rossby_wave 

use constants
use parameters
use admin
implicit none
contains

subroutine amplitudes_array(amplitudes)
implicit none
real(dp) :: amplitudes(65,65)
integer :: a,b
real(dp) :: k_array(65), l_array(65), k, l
k_array = linspace(-3.2,3.2,65)
l_array = linspace(-3.2,3.2,65)

do a = 1,65
do b = 1,65
k = k_array(a)
l = l_array(b)
amplitudes(a,b) = amp_scaling * exp(-(25.*k/32.)**2-(25.*l/32.)**2)*((25.*k/32.)**2 + (25.*l/32.)**2)
end do 
end do
end subroutine amplitudes_array

subroutine dispersion_relation_array(dispersions)
implicit none
real(dp) :: dispersions(65,65)
real(dp) :: k_array(65), l_array(65)
integer :: a,b
real(dp) :: k,l

k_array = linspace(-3.2,3.2,65)
l_array = linspace(-3.2,3.2,65)

do a = 1,65
do b = 1,65
k = k_array(a)
l = l_array(b)
dispersions(a,b) = (-beta * k ) / (k**2+l**2+Rd**(-2))

end do
end do

end subroutine dispersion_relation_array

end module rossby_wave


