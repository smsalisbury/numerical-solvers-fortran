module array
use config
implicit none

!DATA DICTIONARY

public in_array

contains
  function in_array(needle,haystack,e)
    ! INPUT VARIABLES
    real(wp)                    ::  needle
    real(wp),dimension(:)       ::  haystack
    real(wp),optional           ::  e
    
    ! OUTPUT VARIABLES
    logical                     ::  in_array
    
    ! INTERNAL VARIABLES
    integer                     ::  i
    real(wp)                    ::  error
    
    ! SET ERROR
    if (present(e)) then
      error = e
    else
      error = 1.0E-5_wp
    end if
    
    ! SEARCH
    in_array = .FALSE.
    do i=1,size(haystack)
      if (abs(needle - haystack(i)) <= error) then
        in_array = .TRUE.
        return
      end if
    end do
    
  end function in_array
  
  function array_search(needle,haystack,e)
    ! INPUT VARIABLES
    real(wp)                    ::  needle
    real(wp),dimension(:)       ::  haystack
    real(wp),optional           ::  e
    
    ! OUTPUT VARIABLES
    integer                     ::  array_search
    
    ! INTERNAL VARIABLES
    integer                     ::  i
    real(wp)                    ::  error
    
    ! SET ERROR
    if (present(e)) then
      error = e
    else
      error = 1.0E-5_wp
    end if
    
    ! SEARCH
    array_search = -1
    do i=1,size(haystack)
      if (abs(needle - haystack(i)) <= error) then
        array_search = i
        return
      end if
    end do
    
  end function array_search

end module array
