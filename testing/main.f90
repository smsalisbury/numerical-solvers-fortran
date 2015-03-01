program main
use config
use integrate1D
implicit none

real(wp),dimension(:),allocatable	::	x,w

call gauss_quadrature_points(3,x,w)
write(*,*)x
write(*,*)w





end program main