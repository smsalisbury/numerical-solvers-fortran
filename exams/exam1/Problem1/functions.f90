module functions
use config
implicit none

contains
	function f(t,P)
		real(wp)	::	f,P,t
		f = 2.0_wp*P - 0.005_wp*P**2	
	end function
	
	function P(t)
		real(wp)	::	P,t
		P = 10.256_wp*exp(2.0_wp*t)/(1.0_wp+0.02564*exp(2.0_wp*t))
	end function P






end module functions