program main
use config
use polynomial_interpolation
use file_io
implicit none

! DATA DICTIONARY
real(wp),dimension(:),allocatable             ::  t,y
character(20),parameter                       ::  filename='points.xy'

integer                                       ::  i
real(wp),dimension(:),allocatable             ::  w

! READ IN t AND y
call read_xy(filename,t,y)

w = lagrange_interpolation_weights(t,y)
do i=-11,22
  write(*,*)real(i,wp)*0.1_wp,lagrange_interpolation(t,y,real(i,wp)*0.1_wp,w)
end do

! DEALLOCATE ARRAYS
deallocate(t)
deallocate(y)

end program main