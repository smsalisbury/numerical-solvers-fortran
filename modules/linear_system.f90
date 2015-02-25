module linear_system
use config
implicit none

contains
	function det(A)
		!------------------------------------------
		!	This function computes the determinant
		! 	of a 2-dimensional matrix
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),dimension(:,:)				::	A
		
		! OUTPUT VARIABLES
		real(wp)							::	det
		
		! INTERNAL VARIABLES
		integer								::	i
		
		! COMPUTE DETERMINANT
		det = 0.0_wp
	end function det

	subroutine LU(A,L,U)
		!------------------------------------------
		!	This subroutine computes the L and U
		! 	matricies of an LU decomposition of A
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),dimension(:,:),intent(in)					:: 	A
		
		! OUTPUT VARIABLES
		real(wp),dimension(:,:),allocatable,intent(out)		::	L,U ! ** should they be allocatable?
		
		! INTERNAL VARIABLES
		integer												::	i,j,k,i_size,j_size,a_size
		
		! VALIDATION
		i_size = size(A,1)
		j_size = size(A,2)
		! **** Need to validate the size and dimension of A
		
		! SET UP L and U
		a_size = i_size
		allocate(L(a_size,a_size))
		allocate(U(a_size,a_size))
		L = 0.0_wp
		U = 0.0_wp
		
		! COMPUTE L and U
		do i=1,a_size
			do j=1,a_size
				U(i,j) = A(i,j)
				L(i,j) = A(i,j)
				do k=1,(i-1)
					U(i,j) = U(i,j) - U(k,j)*L(i,k)
				end do
				do k=1,(j-1)
					L(i,j) = U(i,j) - U(k,j)*L(i,k)
				end do
				L(i,j) = L(i,j)/U(j,j)
			end do		
		end do
		
	end subroutine LU



end module linear_system