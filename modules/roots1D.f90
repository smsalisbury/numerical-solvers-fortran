module roots1D
use config
use derivate1D
implicit none

contains
	function newton(f,fy,g,e_in)
		!------------------------------------------
		!	This function finds a root of a function
		!	using a specified guess g and Newton's
		!	method.	Note that the derivative of f (fy)
		!	must be specified.
		!------------------------------------------
		
		! INPUTS
		real(wp),external		::	f
		real(wp),external		::	fy
		real(wp)				::	g
		real(wp),optional		::	e_in
		
		! OUTPUTS
		real(wp)				::	newton
		
		! INTERNAL
		real(wp)				::	e
		real(wp)				::	der
		real(wp)				::	x
		real(wp)				::	newton_old
		integer					::	k,max_iter=1000000
		
		! DEFAULT VALUES
		if (present(e_in)) then
			e = e_in
		else
			e = 1.0E-5_wp
		endif
		
		! CALCULATE
		newton = g
		newton_old = newton
		do k=1,max_iter
			der = fy(newton)
			if (abs(der) < e) then
				x = x + 0.1_wp
				cycle
			endif
			newton = newton - f(newton)/der
			if (abs(newton-newton_old) < e) exit
			newton_old = newton
		enddo
	endfunction newton
	
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
		real(wp)				::	der
		real(wp)				::	x
		real(wp)				::	secant_old
		integer					::	k,max_iter=1000000
		
		! DEFAULT VALUES
		if (present(e_in)) then
			e = e_in
		else
			e = 1.0E-5_wp
		endif
		
		! CALCULATE
		secant = g
		secant_old = secant
		do k=1,max_iter
			der = c_diff_first(f,secant)
			if (abs(der) < e) then
				x = x + 0.1_wp
				cycle
			endif
			secant = secant - f(secant)/der
			if (abs(secant-secant_old) < e) exit
			secant_old = secant
		enddo
	endfunction secant
	
end module roots1D