program problem1
use config
use initial_value_1D
use functions
implicit none

integer					::	i
real(wp),dimension(3)	::	h1 = (/0.01_wp,0.005_wp,0.0025_wp/)
real(wp)				::	h2,h_hi,h_low
real(wp)				::	E

do i=1,size(h1)
	E = linear_multistep(f,0.0_wp,5.0_wp,10.0_wp,int(5.0_wp/h1(i)))
	write(*,*)h1(i),E
enddo
write(*,*)'------------------'
h_hi = 0.51_wp
h_low = 0.49_wp
do i=1,11
	h2 = h_low + real(i-1,wp)*(h_hi-h_low)/10.0_wp
	E = linear_multistep(f,0.0_wp,100.0_wp,10.0_wp,int(100.0_wp/h2))
	write(*,*)h2,E
enddo











end program problem1