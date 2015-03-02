program problem6
use config
use global
use integrate2D
use functions
use file_io
implicit none

integer		::	i,n
real(wp)	::	m

do i=1,15
	n = 2**i
	m = monte_carlo(f,0.0_wp,3.0_wp,2.0_wp*pi,1.5_wp*pi,n)
	write(*,*)n,m,abs(m - f_actual)
enddo
	

end program problem6