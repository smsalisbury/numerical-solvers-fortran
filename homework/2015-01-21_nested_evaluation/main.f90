program main
use config
use polynomial_interpolation
implicit none

! DATA DICTIONARY
integer::n,i
real(wp),allocatable::c(:)
real(wp)::x,y,z
real(wp)::nested_time,normal_time
integer::nested_start,nested_finish,normal_start,normal_finish,clock_resolution
integer::k,k_max=100000

! Retrieve number of terms from the user
write(*,*)"Enter the number of terms:"
read(*,*)n

! Retrieve x-vale from the user
write(*,*)"Enter the x-value:"
read(*,*)x

! Allocate the constants array
allocate(c(n))

! Constants for e^x
do i=1,n
  c(i) = 1.0_wp/real(factorial_real(i-1),wp)
end do

! Uncomment this section for ln(1+x)
! c(1) = 0.0_wp
! do i=2,n
  ! c(i) = (-1.0_wp)**(i)/real((i-1),wp)
! end do

! Run the nested evaluation
! This is repeated in a loop because my computer's
! timing resolution is only a millisecond.
call system_clock(nested_start,clock_resolution)
  do k=1,k_max
    y=nested_evaluation(c,x)
  end do
call system_clock(nested_finish)
nested_time = real(nested_finish-nested_start,wp)/(real(clock_resolution,wp)*real(k_max,wp))
write(*,*)'Nested Time: ',nested_time,'Nested Value: ',y

write(*,*)'----------'

! Run the normal polynomial evaluation
call system_clock(normal_start,clock_resolution)
  do k=1,k_max
    z=polynomial_evaluation(c,x)
  end do
call system_clock(normal_finish)
normal_time = real(normal_finish-normal_start,wp)/(real(clock_resolution,wp)*real(k_max,wp))
write(*,*)'Poly Time: ',normal_time,'Poly Value: ',z

! Deallocate the constants array
deallocate(c)

contains
  ! This function calculates the factorial of a number
  ! It returns a real number rather than an integer because
  ! it extends the range of the function (overflow problems)
  !
  ! real factorial_real(int x)
  !   x   The number to calculate the factorial of
  function factorial_real(x)
    ! Dictionary
    integer::x
    real(wp)::ans,factorial_real
    integer::i
    
    ans = 1.0_wp
    do i=1,x
      ans = ans*real(i)
    end do
    
    factorial_real = ans
  end function factorial_real
end program main