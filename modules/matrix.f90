module matrix
use config
implicit none

contains

	function mat_mult(A,B)
		!------------------------------------------
		!	This function multiplies 2 matricies
		!	It is an alias of the intrinsic matmul(A,B)
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),dimension(:,:)								:: 	A,B
		
		! OUTPUT VARIABLES
		real(wp),dimension(:,:),allocatable					::	mat_mult
		
		! INTERNAL VARIABLES
		integer												::	i,j,k,a_i_size,a_j_size,b_i_size,b_j_size
		
		! VALIDATION
		a_i_size = size(A,1)
		a_j_size = size(A,2)
		b_i_size = size(B,1)
		b_j_size = size(B,2)
		if (a_j_size /= b_i_size) then
			write(*,*)"Matrix sizes aren't compatible for multiplication in mat_mult."
			write(*,*)"A is ",a_i_size," by ",a_j_size
			write(*,*)"B is ",b_i_size," by ",b_j_size
			write(*,*)"Aborting..."
			stop
		end if
		
		! SET UP OUTPUT
		allocate(mat_mult(a_i_size,b_j_size))
		mat_mult = matmul(A,B)
	end function mat_mult
	
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
		integer								::	i,i_size,j_size
		real(wp),dimension(:,:),allocatable	::	lu
		logical								::  test
		
		! VALIDATION
		i_size = size(A,1)
		j_size = size(A,2)
		if (i_size /= j_size) then
			write(*,*)"The matrix in det() is not a square matrix."
			write(*,*)"Aborting..."
			stop
		end if
		
		! COMPUTE DETERMINANT
		allocate(lu(i_size,j_size))
		lu = lu_decomp(A)
		
		det = 1.0_wp
		do i=1,i_size
			det = det * lu(i,i)
		end do
	end function det

	function lu_decomp(A)
		!------------------------------------------
		!	This subroutine computes the L and U
		! 	matricies of an LU decomposition of A
		!	This stores in the combined format:
		!	
		!		U	U	U	U
		!		L	U	U	U
		!		L	L	U	U
		!		L	L	L	U
		!
		!	The matricies can then be split later 
		!	using the lu_split() subroutine
		!------------------------------------------
		
		! INPUT VARIABLES
		real(wp),dimension(:,:)								:: 	A
		
		! OUTPUT VARIABLES
		real(wp),dimension(:,:)	,allocatable				:: 	lu_decomp
		
		! INTERNAL VARIABLES
		real(wp),dimension(:,:),allocatable					::	L,U
		integer												::	i,j,k,i_size,j_size,a_size
		
		! VALIDATION
		i_size = size(A,1)
		j_size = size(A,2)
		! **** Need to validate the size and dimension of A
		! **** Also check fo singularity
		
		! SET UP L and U and lu_decomp
		a_size = i_size
		allocate(L(a_size,a_size))
		allocate(U(a_size,a_size))
		allocate(lu_decomp(a_size,a_size))
		L = 0.0_wp
		U = 0.0_wp
		
		! SET DIAGONALS OF L
		do i=1,a_size
			L(i,i) = 1.0_wp
		end do
		
		! COMPUTE L and U and store back in A
		do j=1,a_size
			do i=1,j
				U(i,j) = A(i,j)
				do k=1,(i-1)
					U(i,j) = U(i,j) - U(k,j)*L(i,k)
				end do
				lu_decomp(i,j) = U(i,j)
			end do
			do i=j+1,a_size
				L(i,j) = A(i,j)
				do k=1,(j-1)
					L(i,j) = L(i,j) - U(k,j)*L(i,k)
				end do
				L(i,j) = L(i,j)/U(j,j)
				lu_decomp(i,j) = L(i,j)
			end do
		end do
	end function lu_decomp

end module matrix