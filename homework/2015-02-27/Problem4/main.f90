program problem4
use config
use integrate1D
use functions
use file_io
implicit none

write(*,*)f2_actual
write(*,*)gauss_quadrature(f2,0.0_wp,3.0_wp,10,3)


end program problem4