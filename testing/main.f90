program main
use config
use linear_system
implicit none

real,dimension(3,3)	::	A,A_prime
real,dimension(3)	::	b
integer	::	i

A = transpose(reshape((/ 1,-3,1,2,-8,8,-6,3,-15 /), shape(A)))
b = (/ 4,-2,9 /)
write(*,*)solve_system(A,b)






end program main