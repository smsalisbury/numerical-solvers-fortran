program problem4
use config
use integrate1D
use functions
use file_io
implicit none

integer		::	i,j
real(wp)	::	g
real(wp)	::	a1=2.0_wp*pi,b1=1.5_wp*pi
real(wp)	::	a2=0.0_wp,b2=3.0_wp


do i=1,5
	do j=1,20
		g = gauss_quadrature(f1,a1,b1,i,j)
		write(*,*)i,j,g,abs(g-f1_actual)
	enddo
enddo


end program problem4