program executearray
use passivetracers
use constants
use parameters
implicit none
real(dp) :: phase1(65,65), phase2(65,65), time1(65,65), time2(65,65)
integer :: c_g
real(dp) :: gammas(11)

gammas(1) = 10
gammas(2) = 20
gammas(3) = 30
gammas(4) = 40
gammas(5) = 50
gammas(6) = 60
gammas(7) = 70
gammas(8) = 80
gammas(9) = 90
gammas(10) = 100
gammas(11) = 500


print*, gammas
print*, vu
print*, delta

open(1, file="partavg.dat")
do c_g=1,11

phase1 = tau*random_matrix(65,65)
phase2 = tau*random_matrix(65,65)
time1 = gammas(c_g)*random_matrix(65,65)
time2 = gammas(c_g)*random_matrix(65,65)

write(1,*) vel_passive(time_steps, phase1, phase2, time1, time2, N_particles ,gammas(c_g))

end do

!print *, char(7)  
!call execute_command_line ("powershell.exe -File ../email.ps1", exitstat=i)
!print *, "Exit status of external_prog.exe was ", i

end program executearray
