module integrate2D
use config
use global
implicit none

contains
	function monte_carlo(f,x1,x2,y1,y2,n)
		!------------------------------------------
		!	This function integrates the function f
		!	using a Monte Carlo method.
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),external				::	f
		real(wp)						::	x1,x2,y1,y2
		integer							:: 	n
		
		! OUTPUT VARIABLES
		real(wp)						::	monte_carlo
		
		! INTERNAL VARIABLES
		real(wp),parameter				::	e=1.0E-5_wp
		logical							::	to_convergence
		real(wp)						:: 	rx,ry,sum,x,y
		integer							::	k
		
		! CALCULATE INTEGRAL
		sum = 0.0_wp
		do k=1,n
			call seed_random()
			call random_number(rx)
			call random_number(ry)
			x = rx*(x2-x1)+x1
			y = ry*(y2-y1)+y1
			sum = sum + f(x,y)
		end do
		monte_carlo = (x2-x1)*(y2-y1)*sum/real(n,wp)
	end function monte_carlo
	
end module integrate2D