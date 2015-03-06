program problem4
use config
use functions
use initial_value_1D
implicit none

! DATA DICTIONARY
real(wp)			::	t0,tf,y0,h
integer				::	n

write(*,*)'Enter the size of h:'
read(*,*)h

y0 = 30.0_wp
t0 = 0.0_wp
tf = 20.0_wp
n = int((tf-t0)/h)

write(*,*)taylor2(func,ft,fy,t0,tf,y0,n)

end program problem4