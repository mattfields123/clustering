program executearray
use passivetracers
use constants
use parameters
implicit none
real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
integer :: c_g
real(dp) :: gammas(11)

gammas(1) = 0.2
gammas(2) = 0.25
gammas(3) = 0.3
gammas(4) = 0.35
gammas(5) = 0.45
gammas(6) = 0.5
gammas(7) = 0.6
gammas(8) = 0.7
gammas(9) = 0.8
gammas(10) = 0.9
gammas(11) = 1.0


print*, gammas
print*, vu
print*, delta

open(1, file="partavg.dat")
do c_g=1,11

phase1 = tau*random_matrix(65,65)
phase2 = tau*random_matrix(65,65)
time1 = vu*random_matrix(65,65)
time2 = vu*random_matrix(65,65)

write(1,*) vel_passive(time_steps, phase1, phase2, time1, time2, N_particles ,gammas(c_g))

end do

!print *, char(7)  
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of external_prog.exe was ", i

end program executearray
