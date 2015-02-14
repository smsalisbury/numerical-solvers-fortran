program main
use config
implicit none

! DATA DICTIONARY
real(wp)                  ::  a,h
integer                   ::  i

h = 0.0_wp
do i=1,9
  h = 10.0_wp**(-i)
  a = approx(0.0_wp,h)
  write(*,*)i,h,a,abs(1.0_wp-a)
end do

contains
  function approx(x,h)
  ! INPUT VARIABLES
  real(wp)                  ::  h,x
  
  ! OUTPUT VARIABLES
  real(wp)                  ::  approx
  
  ! CALCULATE
  approx = (-1.0_wp/(2.0_wp*h**3))*exp(x-2.0_wp*h) + (1.0_wp/(h**3))*exp(x-h) &
            + (-1.0_wp/(h**3))*exp(x+h) + (1.0_wp/(2.0_wp*h**3))*exp(x+2.0_wp*h)
  end function approx

end program main