program streamdata


use constants
use parameters
use admin
use streampotential

implicit none
real(dp) :: x_array(250), y_array(250), stream(250,250), pot(250,250)
integer :: i,j
real(dp) :: phase(65,65),phase2(65,65)
real(dp) :: time1(65,65),time2(65,65)

call random_number(phase)
call random_number(phase2)


phase = tau*phase
phase2 = tau*phase2
x_array = linspace(-3.2,3.2,250)
y_array = linspace(-3.2,3.2,250)

do i=1,250
do j=1,250
stream(i,j) = streamfunction(x_array(i),y_array(j),7)
pot(i,j) = potentialfunction(x_array(i),y_array(j),7,time1,time2,phase,phase2)

end do
end do

open(1,file='streamimage.dat')
open(2,file='potimage.dat')

write(1,*) stream
write(2,*) pot



end program streamdata