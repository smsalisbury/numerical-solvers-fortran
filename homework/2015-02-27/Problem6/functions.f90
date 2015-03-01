module functions
use config
implicit none

! Just some functions to integrate
real(wp)			::	f_actual = -9.02462E-6_wp
real(wp),parameter	::	pi = 4.0_wp*atan(1.0_wp)

contains
	function f1(x)
		real(wp)		::	f1,x
		f1=sin(x)/(1.0_wp + 2.0_wp*exp(2.0_wp*x))
	end function f1
	
	function f2(x)
		real(wp)		::	f2,x
		f2=(x**2 - 3.0_wp*x)/((x**2+1.0_wp)*(x**2+2.0_wp))
	end function f2
	
	function f(x,y)
		real(wp)	::	f,x,y
		f = f2(x) * f1(y)
	end function f
	
end module functions
