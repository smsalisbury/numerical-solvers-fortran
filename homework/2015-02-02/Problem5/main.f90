program main
use config
implicit none

! DATA DICTIONARY
real(wp),dimension(-3,3)      ::  c
real(wp),dimension(3)         ::  h = (/ 0.1_wp, 0.05_wp, 0.025_wp /)
real(wp)                      ::  x=0.4_wp,actual=-4.0180396_wp,a

integer                       ::  i,j

write(*,*)f(x)

do i=1,size(h)
  a = approx(x,h(i))
  write(*,*)i,h(i),a,abs(a-actual)
end do

contains
  function approx(x,h)
    ! INPUTS
    real(wp)          :: x,h
    
    ! OUTPUTS
    real(wp)          ::  approx
    
    ! CALCULATE
    approx = 1.0_wp/h**2*((1.0_wp/90.0_wp)*f(x-3.0_wp*h) &
              - (3.0_wp/20.0_wp)*f(x-2.0_wp*h) &
              + (3.0_wp/2.0_wp)*f(x-h) &
              - (49.0_wp/18.0_wp)*f(x) &
              + (3.0_wp/2.0_wp)*f(x+h) &
              - (3.0_wp/20.0_wp)*f(x+2.0_wp*h) &
              + (1.0_wp/90.0_wp)*f(x+3.0_wp*h))
  end function approx
  
  function f(x)
    ! INPUTS
    real(wp)        ::  x
    
    ! OUTPUTS
    real(wp)        ::  f
    
    ! CALCULATE
    f = exp(x)*sin(3.0_wp*x)
  end function
  
end program main