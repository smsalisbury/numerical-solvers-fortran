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
		endif
		
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
	
	subroutine gauss_quadrature_points(n,x,w)
		!------------------------------------------
		!	This function uses the Legendre Polynomials
		!	to calculate the Gauss Points and Weights
		!	for a specified number of n points
		!------------------------------------------
		
		! INPUT VARIABLES
		integer,intent(in)								:: 	n
		
		! OUTPUT VARIABLES
		real(wp),dimension(:),allocatable,intent(out)	::	x,w
		
		! INTERNAL VARIABLES
		
		! ALLOCATE OUTPUTS
		allocate(x(n),w(n))
		
		! CALCULATE INTEGRAL
	end subroutine gauss_quadrature_points
  
  function trapezoid(f,a,b)
		!------------------------------------------
		!	Evaluate an integral approximation from
		!	a to b using the trapezoid rule.
		!------------------------------------------
		
		! INPUTS
		real(wp),external			::	f
		real(wp)					::	a,b
		
		! OUTPUTS
		real(wp)					::	trapezoid
		
		! CALCULATE
		trapezoid = trapezoid_values(f(a),f(b),a,b)
  end function trapezoid
  
  function trapezoid_values(fa,fb,a,b)
		!------------------------------------------
		!	Evaluate an integral approximation from
		!	a to b using the trapezoid rule, where
		!	explicit function values are provided
		!	rather than calculated. This is useful
		!	if the function is called recursively on
		!	the same interval.
		!------------------------------------------
		
		! INPUTS
		real(wp)					::	fa,fb
		real(wp)					::	a,b
		
		! OUTPUTS
		real(wp)					::	trapezoid_values
		
		! CALCULATE
		trapezoid_values = 0.5_wp*(fa + fb)*(b-a)
  end function trapezoid_values
  
  function simpson(f,a,b)
		!------------------------------------------
		!	Evaluate an integral approximation from
		!	a to b using the Simpson rule.
		!------------------------------------------
		
		! INPUTS
		real(wp),external			::	f
		real(wp)					::	a,b
		
		! OUTPUTS
		real(wp)					::	simpson
		
		! CALCULATE
		simpson = simpson_values(f(a),f(real(a+b,wp)/2.0_wp),f(b),a,b)
  end function simpson
  
  function simpson_values(fa,fab,fb,a,b)
		!------------------------------------------
		!	Evaluate an integral approximation from
		!	a to b using the Simpson's rule, where
		!	explicit function values are provided
		!	rather than calculated. This is useful
		!	if the function is called recursively on
		!	the same interval.
		!------------------------------------------
		
		! INPUTS
		real(wp)					::	fa,fab,fb
		real(wp)					::	a,b
		
		! OUTPUTS
		real(wp)					::	simpson_values
		
		! CALCULATE
		simpson_values = (fa + 4.0_wp*fab + fb)*(b-a)/6.0_wp
  end function simpson_values
  
  function adaptive_quadrature(f,a,b,method_in,e_in)
		!------------------------------------------
		!	Use adaptive quadrature to evaluate an
		!	integral
		!------------------------------------------
		
		! INPUTS
		real(wp),external			::	f
		real(wp)					::	a,b
		character(*),optional		::	method_in
		real(wp),optional			::	e_in
		
		! OUTPUTS
		real(wp)					::	adaptive_quadrature
		
		! INTERNAL VARIABLES
		real(wp)					::	e
		character(40)				::	method
		real(wp)					::	fa,fb,fab
		real(wp)					::	first_try
		
		! DEFAULT VALUES
		if (present(method_in)) then
			method = method_in
		else
			method = 'simpson'
		endif
		if (present(e_in)) then
			e = e_in
		else
			e = 1.0E-5_wp
		endif
		
		! CALCULATE INITIAL VALUES
		if (method .EQ. 'trapezoid') then
			fa = f(a)
			fb = f(b)
			first_try = trapezoid_values(fa,fb,a,b)
		else if ((method .EQ. 'simpson') .OR. (method .EQ. 'simpsons')) then
			fa = f(a)
			fb = f(b)
			fab = f(real(a+b,wp)/2.0_wp)
			first_try = simpson_values(fa,fab,fb,a,b)
		endif
		
		! RECURSIVELY CALCULATE THE INTEGRAL
		! Use the helper function
		if (method .EQ. 'trapezoid') then
			adaptive_quadrature = adaptive_quadrature_helper(f,(/ fa,fb /),a,b,method,first_try,10,e)
		else if ((method .EQ. 'simpson') .OR. (method .EQ. 'simpsons')) then
			adaptive_quadrature = adaptive_quadrature_helper(f,(/ fa,fab,fb /),a,b,method,first_try,10,e)
		endif
		
  end function adaptive_quadrature
  
  recursive function adaptive_quadrature_helper(f,f_values,a,b,method,first,level,e) result(res)
		!------------------------------------------
		!	This is a helper function explicitly
		!	used in the adaptive_quadrature function.
		!------------------------------------------
		
		! INPUTS
		real(wp),dimension(:)			::	f_values
		real(wp),external				::	f
		real(wp)						::	a,b
		character(*)					::	method
		real(wp)						::	first
		integer							::	level
		real(wp)						::	e
		
		! OUTPUT
		real(wp)						::	res
		
		! INTERNAL
		real(wp)						::	sub1,sub2,next_sub1,next_sub2
		real(wp)						::	ab,aab,abb
		real(wp)						::	fab,faab,fabb
		
		! VALIDATE
		if (method .EQ. 'trapezoid' .AND. size(f_values) /= 2) then
			write(*,*)"The number of function values and the method are not consistent in adaptive_quadrature_helper."
			write(*,*)"The method is ",method," and there are ",size(f_values)," function values."
			write(*,*)"Aborting..."
			stop
		else if (((method .EQ. 'simpson') .OR. (method .EQ. 'simpsons')) .AND. size(f_values) /= 3) then
			write(*,*)"The number of function values and the method are not consistent in adaptive_quadrature_helper."
			write(*,*)"The method is ",method," and there are ",size(f_values)," function values."
			write(*,*)"Aborting..."
			stop
		endif
		
		! CALCULATE SUB INTEGRALS
		if (method .EQ. 'trapezoid') then
			ab = real(a+b,wp)/2.0_wp
			fab = f(ab)
			sub1 = trapezoid_values(f_values(1),fab,a,ab)
			sub2 = trapezoid_values(fab,f_values(2),ab,b)
		else if ((method .EQ. 'simpson') .OR. (method .EQ. 'simpsons')) then
			ab = real(a+b,wp)/2.0_wp
			aab = (3.0_wp*real(a,wp)+real(b,wp))/4.0_wp
			abb = (real(a,wp)+3.0_wp*real(b,wp))/4.0_wp
			faab = f(aab)
			fabb = f(abb)
			sub1 = simpson_values(f_values(1),faab,f_values(2),a,ab)
			sub2 = simpson_values(f_values(2),fabb,f_values(3),ab,b)
		endif
		
		if (abs(first - (sub1+sub2)) < e) then
			! If the subdivision is close enough to the original, return the value
			res = sub1+sub2
		else if (level == 0) then
			! Maximum subdivision levels reached
			write(*,*)"Maximum subdivition levels reached."
			res = sub1 + sub2
		else
			! Subdivide some more
			if (method .EQ. 'trapezoid') then
				next_sub1 = adaptive_quadrature_helper(f,(/ f_values(1),fab /),a,ab,method,sub1,level-1,e)
				next_sub2 = adaptive_quadrature_helper(f,(/ fab,f_values(2) /),ab,b,method,sub2,level-1,e)
			else if ((method .EQ. 'simpson') .OR. (method .EQ. 'simpsons')) then
				next_sub1 = adaptive_quadrature_helper(f,(/ f_values(1),faab,fab /),a,ab,method,sub1,level-1,e)
				next_sub2 = adaptive_quadrature_helper(f,(/ fab,fabb,f_values(2) /),ab,b,method,sub2,level-1,e)
			endif
			res = next_sub1 + next_sub2
		endif
  end function adaptive_quadrature_helper
  
  
  
  

end module integrate1D