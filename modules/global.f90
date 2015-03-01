module global
use config
implicit none

contains
	subroutine seed_random()
		integer	::	clock,i,n
		integer,dimension(:),allocatable :: seed
		
		call random_seed(size = n)
		allocate(seed(n))
		
		call system_clock(clock)
		
		seed = clock + 37*(/ (i-1, i=1, n) /)
		call random_seed(PUT=seed)
		
		deallocate(seed)
	end subroutine seed_random


end module global