program main
use config
use integrate1D
use functions
use file_io
implicit none

! DATA DICTIONARY
real(wp),parameter				::	pi = 4.0_wp*atan(1.0_wp)
real(wp),parameter				::	f1_actual = 0.0000164883_wp
real(wp),parameter				::	f2_actual = -0.54734_wp

integer							::	k,n
real(wp)						::	mc1,mc2,diff

do k=1,20
	n = 2**k
	!mc1 = monte_carlo(f2,2.0_wp*pi,1.5_wp*pi,n)
	mc2 = monte_carlo(f2,0.0_wp,3.0_wp,n)
	diff = abs(mc2 - f2_actual)
	write(*,*)n,mc2,diff
end do

end program main