program executearray
use passivetracers
use constants
use parameters
implicit none
real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
integer :: c_g
real(dp) :: gammas(11)

phase1 = tau*random_matrix(65,65)
phase2 = tau*random_matrix(65,65)
! time1 = vu*random_matrix(65,65)
! time2 = vu*random_matrix(65,65)

gammas(1) = 5**tau 
gammas(2) = 10*tau
gammas(3) = 25*tau
gammas(4) = 50*tau
gammas(5) = 100*tau
gammas(6) = 250*tau
gammas(7) = 500*tau
gammas(8) = 1000*tau
gammas(9) = 1*tau
gammas(10) = 0.5*tau
gammas(11) = 0.01*tau


print*, gammas
open(1, file="partavg.dat")
do c_g=1,11

time1 = gammas(c_g) * random_matrix(65,65)
time2 = gammas(c_g) * random_matrix(65,65)

write(1,*) vel_passive(time_steps, phase1, phase2, time1, time2, N_particles ,gammas(c_g))

end do

!print *, char(7)  
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of external_prog.exe was ", i

end program executearray
