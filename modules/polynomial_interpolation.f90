module polynomial_interpolation
use config
implicit none

public nested_evaluation, polynomial_evaluation

contains
  ! This function interpolates a polynomial interpolation function
  ! using the nested evaluation method:
  !
  !   polynomial evaluation:
  !     f(x) = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
  !   
  !   nested evaluation:
  !     f(x) = c_0 + x*(c_1 + x*(c_2 + ... x(c_n-1 + c_n*x)))
  !
  ! The nested evaluation method is generally quicker than the traditional
  ! term-by-term method.
  !
  ! real nested_evaluation(array(:) c, real x)
  !   c   The coefficient array. This is a 1-D array of any size.
  !   x   The value at which to evaluate the polynomial.
  function nested_evaluation(c,x)
    !INPUT VARIABLES
    real(wp),dimension(:)::c
    real(wp)::x
    
    !OUTPUT VARIABLES
    real(wp)::nested_evaluation
    
    !INTERNAL VARIABLES
    real(wp)::s
    integer::j,n
    
    !CONSTRAINTS
    n=size(c)
    if (n < 1) then
      write(*,*)'The constants array must be bigger than size 0.'
      stop
    end if
    
    !CALCULATE
    s=c(n)
    do j=n-1,1,-1
      s=c(j) + s*x
    end do
    
    !RETURN VALUE
    nested_evaluation = s
  end function nested_evaluation
  
  ! This function interpolates a polynomial interpolation function
  ! using the traditional term-by-term evaluation method:
  !
  !   polynomial evaluation:
  !     f(x) = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
  !
  ! The nested evaluation method is generally quicker than the traditional
  ! term-by-term method.
  !
  ! real polynomial_evaluation(array(:) c, real x)
  !   c   The coefficient array. This is a 1-D array of any size.
  !   x   The value at which to evaluate the polynomial.
  function polynomial_evaluation(c,x)
    !INPUT VARIABLES
    real(wp),dimension(:)::c
    real(wp)::x
    
    !OUTPUT VARIABLES
    real(wp)::polynomial_evaluation
    
    !INTERNAL VARIABLES
    real(wp)::s
    integer::j,n
    
    !CONSTRAINTS
    n=size(c)
    if (n < 1) then
      write(*,*)'The constants array must be bigger than size 0.'
      stop
    end if
    
    !CALCULATE
    s=c(1)
    do j=2,n
      s=s + c(j)*(x**(j-1))
    end do
    
    !RETURN VALUE
    polynomial_evaluation = s
  end function polynomial_evaluation

end module polynomial_interpolation