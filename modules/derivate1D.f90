module derivate1D
use config
implicit none

contains
	function c_diff_first(f,x,h_in)
		!------------------------------------------
		!	This function finds the derivative of f
		!	at point x using a Central Difference
		!	approximation (o(h^2)). If h is not
		!	specified, the function iterates until
		!	a converged value is found.
		!------------------------------------------
		
		! INPUTS
		real(wp),external		::	f
		real(wp)				::	x
		real(wp),optional		::	h_in
		
		! OUTPUTS
		real(wp)				::	c_diff_first
		
		! INTERNAL VARIABLES
		real(wp)				::	c_diff_first_old
		real(wp)				::	h
		integer					::	k,max_iter = 10
		real(wp)				::	e = 1.0E-5_wp
		
		! CALCULATE
		if (present(h_in)) then
			h = h_in
			c_diff_first = (f(x+h)-f(x-h))/(2.0_wp*h)
		else
			c_diff_first_old = 0.0_wp
			do k=1,max_iter
				h = 2.0_wp**(-k)
				c_diff_first = (f(x+h)-f(x-h))/(2.0_wp*h)
				if (k .EQ. 1) then
					c_diff_first_old = c_diff_first
					cycle
				endif
				
				if (abs(c_diff_first_old - c_diff_first) < e) exit
				c_diff_first_old = c_diff_first
			enddo
		endif
	end function c_diff_first
	
	function f(x)
		real(wp) :: f,x
		f = x**2 - 4*x - 5
	end function f
	
end module derivate1D