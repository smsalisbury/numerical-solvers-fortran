module polynomial_interpolation
use config
use array
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
  
  function lagrange_interpolation(x_values,y_values,x,weights)
    ! INPUT VARIABLES
    real(wp),dimension(:)           ::  x_values,y_values
    real(wp),dimension(:),optional  ::  weights
    real(wp)                        ::  x

    ! OUTPUT VARIABLES      
    real(wp)                        ::  lagrange_interpolation
    
    ! INTERNAL VARIABLES
    real(wp)                        ::  numerator,denominator
    integer                         ::  m,n
    integer                         ::  j
    integer                         ::  test
    real(wp),dimension(:),allocatable ::  w
    
    ! CONSTRAINTS
    n=size(x_values)
    m=size(y_values)
    if (n < 1 .OR. m < 1) then
      write(*,*)'Both point arrays must be bigger than size 0.'
      stop
    else if (n /= m) then
      write(*,*)'The point arrays are not the same size.'
      stop
    end if
    
    ! CHECK TO SEE IF VALUE IS GIVEN
    test = array_search(x,x_values)
    if (test /= -1) then
      lagrange_interpolation = y_values(test)
      return
    end if
    
    ! CALCULATE
    ! -- Get weights if not present
    if (.NOT. present(weights)) then
      w = lagrange_interpolation_weights(x_values,y_values)
    else
      allocate(w(size(weights)))
      w = weights
    end if
    
    ! -- Calculate numerator and denominator
    numerator = 0.0_wp
    denominator = 0.0_wp
    do j=1,n
      numerator = numerator + w(j)*y_values(j)/(x - x_values(j))
      denominator = denominator + w(j)/(x - x_values(j))
    end do
    
    ! RETURN VALUE
    lagrange_interpolation = numerator/denominator
  end function lagrange_interpolation
  
  function lagrange_interpolation_weights(x_values,y_values)
    ! INPUT VARIABLES
    real(wp),dimension(:)           ::  x_values,y_values

    ! OUTPUT VARIABLES      
    real(wp),allocatable,dimension(:) ::  lagrange_interpolation_weights
    
    ! INTERNAL VARIABLES
    real(wp)                        ::  p
    integer                         ::  m,n
    integer                         ::  i,j
    
    ! CONSTRAINTS
    n=size(x_values)
    m=size(y_values)
    if (n < 1 .OR. m < 1) then
      write(*,*)'Both point arrays must be bigger than size 0.'
      stop
    else if (n /= m) then
      write(*,*)'The point arrays are not the same size.'
      stop
    end if
    
    ! CALCULATE
    allocate(lagrange_interpolation_weights(n))
    
    ! -- Barycentric weights
    do j=1,n
      p = 1.0_wp
      do i=1,n
        if (i == j) cycle
        p = p * (x_values(j) - x_values(i))
      end do
      lagrange_interpolation_weights(j) = 1.0_wp/p
    end do
    
  end function lagrange_interpolation_weights

end module polynomial_interpolation