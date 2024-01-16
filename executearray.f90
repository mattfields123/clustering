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
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)

gammas(1) = 0.
gammas(2) = 0.005
gammas(3) = 0.01
gammas(4) = 0.05
gammas(5) = 0.075
gammas(6) = 0.1
gammas(7) = 0.125
gammas(8) = 0.25
gammas(9) = 0.5
gammas(10) = 0.75
gammas(11) = 1.


print*, gammas
open(1, file="partavg.dat")
do c_g=1,11

write(1,*) vel_passive(time_steps, phase1, phase2, time1, time2, N_particles ,gammas(c_g))

end do

!print *, char(7)  
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of external_prog.exe was ", i

end program executearray
