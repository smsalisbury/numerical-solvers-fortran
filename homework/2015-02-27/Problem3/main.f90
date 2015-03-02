program problem3
use config
use integrate1D
use functions
use file_io
implicit none

real(wp)	::	a1,a2

a1 = adaptive_quadrature(f1,2.0_wp*pi,1.5_wp*pi,'trapezoid',1.0E-6_wp)
a2 = adaptive_quadrature(f2,0.0_wp,3.0_wp,'trapezoid',1.0E-6_wp)

write(*,*)a1,(a1-f1_actual)
write(*,*)a2,(a2-f2_actual)


end program problem3