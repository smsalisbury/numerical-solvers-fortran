program problem5
use config
use integrate1D
use functions
use file_io
implicit none

write(*,*)adaptive_quadrature(f2,0.0_wp,3.0_wp,'gauss')


end program problem5