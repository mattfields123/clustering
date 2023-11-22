module phase
use constants
use rossby_wave
use admin
MaduLawrence(ktime ,Nvu)

contains

subroutine generatephases(2,array,Tsteps,vu,amplitudes)
implicit none
integer :: Tsteps
real(dp),allocatable :: array(:,:,:,:)
real(dp) :: phase1(65,65),phase2(65,65)
integer :: t,k,l

allocate(array(2,65,65,Tsteps))

do t = 1,Tsteps

if (mod(t,vu) == 1) then
phase1 = random_number(phase1)
phase2 = random_number(phase2)
end if
do k = 1,65
do l = 1,65
        MLtime = mod(t,vu)
        array(1,k,l,t) = MaduLawrence(MLtime,vu) * amplitudes(k,l)

      


array(k,l,t) = random

end do
end do
end do




end subroutine


end module phase

