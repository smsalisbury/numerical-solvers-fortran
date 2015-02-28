program problem3
use config
use integrate1D
use functions
use file_io
implicit none

write(*,*)adaptive_quadrature(f1,2.0_wp*pi,1.5_wp*pi,'trapezoid',1.0E-6_wp)


end program problem3