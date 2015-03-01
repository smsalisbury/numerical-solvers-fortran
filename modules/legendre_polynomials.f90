module legendre_polynomials
use config
implicit none

contains
	function legendre_poly_root(n,g,e_in)
		!------------------------------------------
		!	This function finds a root of the Legendre
		!	polynomials. There is a separate function for this
		!	because the function that generates the polynomial
		!	cannot be used in the secant function.
		!------------------------------------------
		
		! INPUTS
		real(wp)				::	g
		integer					::	n
		real(wp),optional		::	e_in
		
		! OUTPUTS
		real(wp)				::	legendre_poly_root
		
		! INTERNAL
		real(wp)				::	e
		real(wp)				::	x
		real(wp)				::	x_old
		real(wp)				::	der
		integer					::	k,max_iter=1000
		
		! DEFAULT VALUES
		if (present(e_in)) then
			e = e_in
		else
			e = 1.0E-5_wp
		endif
		
		! CALCULATE
		x = g
		x_old = x
		do k=1,max_iter
			der = legendre_poly_1der(n,x)
			if (abs(der) < e) then
				x = x + 0.1_wp
				cycle
			endif
			x = x - legendre_poly(n,x)/legendre_poly_1der(n,x)
			if (abs(x - x_old) < e) exit
			x_old = x
		enddo
		legendre_poly_root = x
	end function legendre_poly_root
	
	function legendre_poly_1der(n,x,h_in)
		!------------------------------------------
		!	This function finds the derivative of
		!	a Legendre Polynomial of degree n.
		!------------------------------------------
		
		! INPUTS
		real(wp)				::	x
		integer					::	n
		real(wp),optional		::	h_in
		
		! OUTPUTS
		real(wp)				::	legendre_poly_1der
		
		! INTERNAL
		real(wp)				::	legendre_poly_1der_old
		real(wp)				::	h
		integer					::	k,max_iter = 10
		real(wp)				::	e = 1.0E-5_wp
		
		! CALCULATE
		if (present(h_in)) then
			h = h_in
			legendre_poly_1der = (legendre_poly(n,x+h)-legendre_poly(n,x-h))/(2.0_wp*h)
		else
			legendre_poly_1der_old = 0.0_wp
			do k=1,max_iter
				h = 2.0_wp**(-k)
				legendre_poly_1der = (legendre_poly(n,x+h)-legendre_poly(n,x-h))/(2.0_wp*h)
				if (k .EQ. 1) then
					legendre_poly_1der_old = legendre_poly_1der
					cycle
				endif
				
				if (abs(legendre_poly_1der_old - legendre_poly_1der) < e) exit
				legendre_poly_1der_old = legendre_poly_1der
			enddo
		endif
	end function legendre_poly_1der
	
	recursive function legendre_poly(n,x) result(res)
		!------------------------------------------
		!	Returns a Legendre Polynomial of degree
		!	n.
		!------------------------------------------
		
		! INPUTS
		real(wp)				::	x
		integer					::	n
		
		! OUTPUTS
		real(wp)				::	res
		
		! INTERNAL VARIABLES
		real(wp)				::	j
		
		! CALCULATE
		if (n < 0) then
			write(*,*)"Cannot create a Legendre Polynomial of a negative degree."
			write(*,*)"Aborting..."
		else if (n .EQ. 0) then
			res = 1.0_wp
		else if (n .EQ. 1) then
			res = x
		else
			j = real(n-1,wp)
			res = ((2.0_wp*j+1.0_wp)/(j+1.0_wp))*x*legendre_poly(n-1,x) - (j/(j+1.0_wp))*legendre_poly(n-2,x)
		endif
	end function legendre_poly


end module legendre_polynomials