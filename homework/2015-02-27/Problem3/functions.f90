module functions
use config
implicit none

! Just some functions to integrate

contains
	function f1(x)
		real(wp)		::	f1,x
		f1=sin(x)/(1.0_wp + 2.0_wp*exp(2.0_wp*x))
	end function f1
	
	function f2(x)
		real(wp)		::	f2,x
		f2=(x**2 - 3.0_wp*x)/((x**2+1.0_wp)*(x**2+2.0_wp))
	end function f2
	
end module functions
