program main
use config
use file_io
implicit none

! DATA DICTIONARY
integer                 ::  output_file
integer                 ::  i
real(wp)                ::  k,h
real(wp)                ::  actual
real(wp)                ::  approx

! SET UP FILE
call create_dir('output')
output_file = file_open('output/data.out')

! ITERATE AND SOLVE
do i=0,16
  k = 0.5_wp*real(i,wp)
  h = 10.0_wp**(-k)
  actual = -sin(1.2_wp)
  approx = second_derivative(f,1.2_wp,h)
  write(output_file,*)i,k,h,approx,abs(approx-actual)
  write(*,*)i,k,h,approx,abs(approx-actual)
end do

! CLOSE FILE
close(output_file)

contains
  function second_derivative(f,x,h)
    ! INPUTS
    real(wp)            ::  x,h
    real(wp),external   ::  f
    
    ! OUTPUTS
    real(wp)            ::  second_derivative
    
    ! CALCULATION
    second_derivative = (f(x-h) - 2.0_wp*f(x) + f(x+h))/h**2
  end function second_derivative

  function f(x)
    ! INPUTS
    real(wp)            ::  x
    
    ! OUTPUTS
    real(wp)            ::  f
    
    ! CALCULATION
    f = sin(x)  
  end function f




end program main