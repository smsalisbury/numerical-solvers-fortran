program problem5
use config
use integrate1D
use functions
use file_io
implicit none

real(wp)	::	a1,a2

a1 = adaptive_quadrature(f1,2.0_wp*pi,1.5_wp*pi,'gauss')
a2 = adaptive_quadrature(f2,0.0_wp,3.0_wp,'gauss')

write(*,*)a1,abs(a1-f1_actual)
write(*,*)a2,abs(a2-f2_actual)


end program problem5