module linear_system
use config
use matrix
implicit none

contains
	function solve_system(A,b,method_in)
		!------------------------------------------
		!	This function finds a solution to the
		!	system of equations Ax=b
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),dimension(:,:)				::	A
		real(wp),dimension(:)				::	b
		character(20),optional				::	method_in
		
		! OUTPUT VARIABLES
		real(wp),dimension(:),allocatable	::	solve_system
		
		! INTERNAL VARIABLES
		character(20)						::	method
		integer								::	i_size,j_size,i,j,k
		real(wp),dimension(:,:),allocatable	::	lu				! FOR LU DECOMPOSITION METHOD
		real(wp),dimension(:),allocatable	::	y,x				! FOR LU DECOMPOSITION METHOD
		
		! DEFAULT VALUES
		if (present(method_in)) then
			method = method_in
		else
			method = 'lu'
		end if
		
		! VALIDATION
		i_size = size(A,1)
		j_size = size(A,2)
		if (i_size /= j_size) then
			write(*,*)"The coefficient matrix in solve_system() is not a square matrix."
			write(*,*)"Aborting..."
			stop
		else if ((det(A) - 0.0_wp) < 1.0E-5_wp) then
			write(*,*)"The coefficient matrix in solve_system() is singular or nearly singular."
			write(*,*)"det = ",det(A)
		end if
		
		! ALLOCATE RETURN
		allocate(solve_system(i_size))
		
		! COMPUTE SOLUTION
		select case (method)
			case default
				! LU decomposition method
				allocate(y(i_size),x(i_size))
				allocate(lu(i_size,j_size))
				lu = lu_decomp(A)
				y = 0.0_wp
				x = 0.0_wp
				do i=1,i_size
					y(i) = b(i)
					do k=1,i-1
						y(i) = y(i) - lu(i,k)*y(k)
					end do
				end do
				do i=i_size,1,-1
					x(i) = y(i)
					do k=i+1,i_size
						x(i) = x(i) - lu(i,k)*x(k)
					end do
					x(i) = x(i)/lu(i,i)
				end do
				solve_system = x
		end select
	end function solve_system



end module linear_system