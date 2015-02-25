program main
use config
use linear_system
implicit none

real,dimension(3,3)	::	A
real,dimension(:,:),allocatable	::	L,U
integer	::	i

A = transpose(reshape((/ 1,3,5,2,4,7,1,1,0 /), shape(A)))
call LU(A,L,U)

do i=1,3
	write(*,*)A(i,:)
end do
write(*,*)
do i=1,3
	write(*,*)L(i,:)
end do
write(*,*)
do i=1,3
	write(*,*)U(i,:)
end do






end program main