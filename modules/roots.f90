module roots
use config
use derivate1D
implicit none

contains
	function secant(f,g,e_in)
		!------------------------------------------
		!	This function finds a root of a function
		!	using a specified guess g and the secant
		! 	variation of Newton's method.
		!------------------------------------------
		
		! INPUTS
		real(wp),external		::	f
		real(wp)				::	g
		real(wp),optional		::	e_in
		
		! OUTPUTS
		real(wp)				::	secant
		
		! INTERNAL
		real(wp)				::	e
		integer					::	k,max_iter=1000000
		
		! DEFAULT VALUES
		if (present(e_in)) then
			e = e_in
		else
			e = 1.0E-5_wp
		endif
		
		! CALCULATE
		secant = g
		do k=1,max_iter
			secant = secant - f(secant)/c_diff_first(f,secant)
		enddo
	end function secant
	
end module roots