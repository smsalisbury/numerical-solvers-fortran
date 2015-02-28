module matrix
implicit none

contains
  function solve_system(A,b)
    ! INPUTS
    real(wp),dimension(*,*)             ::  A
    real(wp),dimension(*)               ::  b
    
    ! OUTPUTS
    real(wp),dimension(:),allocatable   ::  solve_system
    
    ! ALLOCATE ARRAY
    allocate(solve_system(size(b)))
    
    
  end function







end module matrix