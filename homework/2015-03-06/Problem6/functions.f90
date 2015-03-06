module functions
use config
implicit none

! Just some functions

contains
	function f(t,y)
		real(wp)		::	f,y,t
		f = 0.3_wp*y - 0.009*y**2
	end function f
	
end module functions
