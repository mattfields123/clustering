module rossby_wave_attributes
use constants
use parameters
implicit none

contains
    function amplitude(k, l) result(amp)
        real(dp) :: k, l, amp

        amp = exp(-k**2 - l**2) * (k**2 + l**2)
    
    end function amplitude

    function amplitude2(k, l) result(amp)
        real(dp) :: k, l, amp

        amp = amp_scaling * exp(-(25.*k/32.)**2 - (25.*l/32.)**2) * ((25.*k/32.)**2 + (25.*l/32.)**2)
    
    end function amplitude2

    function dispersion(k, l) result(omega)
        real(dp) :: k, l, omega

        omega = -beta * k/(k**2 + l**2 + Rd**(-2))
        !omega = 0.
 
    end function dispersion
end module rossby_wave_attributes