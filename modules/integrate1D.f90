module integrate1D
use config
implicit none

contains
	function monte_carlo(f,a,b,n)
		!------------------------------------------
		!	This function integrates the function f
		!	using a Monte Carlo method. If n is not
		!	specified, then the method runs until it
		!	converges.
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),external				::	f
		real(wp)						::	a,b
		integer,optional				:: 	n
		
		! OUTPUT VARIABLES
		real(wp)						::	monte_carlo
		
		! INTERNAL VARIABLES
		real(wp),parameter				::	e=1.0E-5_wp
		logical							::	to_convergence
		real(wp)						:: 	r,sum,x
		integer							::	k
		
		! DEFAULT VALUES
		if (present(n)) then
			to_convergence = .FALSE.
		else
			to_convergence = .TRUE.
		end if
		
		! CALCULATE INTEGRAL
		!if (to_convergence) then
		
	!	else
			sum = 0.0_wp
			do k=1,n
				call random_number(r)
				x = r*(b-a)+a
				sum = sum + f(x)
			end do
			monte_carlo = (b-a)*sum/real(n,wp)
	!	end if
	end function monte_carlo

end module integrate1D