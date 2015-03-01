program problem2
use config
use integrate1D
use functions
use file_io
implicit none

! DATA DICTIONARY
real(wp),parameter				::	pi = 4.0_wp*atan(1.0_wp)
real(wp),parameter				::	f1_actual = 0.0000164883_wp
real(wp),parameter				::	f2_actual = -0.54734_wp

integer							::	k,n,i,j
real(wp)						::	mc1,mc2,diff,sum
real(wp)						::	a = 0.0_wp
real(wp)						::	b = 3.0_wp

integer,dimension(4)			::	intervals = (/ 10,20,40,80 /)

do i=1,size(intervals)
	write(*,*)intervals(i)," intervals"
	do k=1,10
		n = 2**k
		sum = 0.0_wp
		do j=1,intervals(i)
			sum = sum + monte_carlo(f2, a+real(j-1,wp)*(b-a)/real(intervals(i),wp), a+real(j,wp)*(b-a)/real(intervals(i),wp), n)
		enddo
		mc2 = sum
		diff = abs(mc2 - f2_actual)
		write(*,*)n,mc2,abs(diff/f2_actual)*100.0_wp
	enddo
	write(*,*)
enddo

end program problem2