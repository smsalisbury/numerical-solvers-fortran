module initial_value_1D
use config
use file_io
use functions
implicit none

contains
	function explicit_euler(f,x_0,x_f,y_0,n)
		!	--------------------------------
		!	This function evaluates an initial
		!	value problem from x_0 to x_f with
		!	initial value y_0 and taking steps
		!	of size (x_f-x_0)/n.
		!	--------------------------------
		
		! 	INPUTS
		real(wp),external			::	f					! The derivative of the solution, as a function of y and x
		real(wp)					::	x_0,x_f				! The domain limits
		real(wp)					::	y_0					! The initial value
		integer						::	n					! The number of steps to take
		
		!	OUTPUTS
		real(wp)					::	explicit_euler		! Returns the final estimate at x_f
		
		!	INTERNAL
		integer						::	output_file
		integer						::	i
		real(wp)					::	h,y,x
		
		!	OPEN OUTPUT FILE
		call create_dir('output')
		output_file = file_open('output/explicit_euler.out')
		
		!	TIME STEPS
		y = y_0
		x = x_0
		h = (x_f-x_0)/real(n,wp)
		write(output_file,*)x,y
		do i=1,n
			x = x + h
			y = y + h*f(x,y)
			write(output_file,*)x,y
		enddo
		explicit_euler = y
		
		!	CLOSE OUTPUT FILE
		close(output_file)	
	end function
	
	function linear_multistep(f,x_0,x_f,y_0,n,order_in)
		!	--------------------------------
		!	This function performs a linear
		!	multistep method to approximate y.
		! 	This uses another method to start
		!	the approximation.
		!	--------------------------------
		
		!	INPUTS
		real(wp),external			::	f
		real(wp)					::	x_0,x_f
		real(wp)					::	y_0
		integer						::	n
		integer,optional			::	order_in
		
		!	OUTPUTS
		real(wp)					::	linear_multistep
		
		!	INTERNAL
		integer								::	output_file
		integer								::	rk_file
		integer								::	order
		integer								::	i
		real(wp),dimension(:),allocatable	::	Y_values,f_values
		real(wp)							::	h
		real(wp)							::	x,y
		real(wp)							::	tmp,error
		
		!	DEFAULTS
		if (present(order_in)) then
			order = order_in
		else
			order = 4
		endif
		
		!	OPEN OUTPUT FILE
		call create_dir('output')
		output_file = file_open('output/linear_multistep.out')
		
		!	ALLOCATE ARRAYS
		allocate(Y_values(n+1),f_values(n+1))
		
		! 	CALCULATE STEP SIZE
		h = (x_f-x_0)/real(n,wp)
		
		!	METHOD INITIALIZATION
		!	Use a RK method of proper order to initialize
		if (order <= 2) then
			! Use RK2
		else if (order <= 4) then
			tmp = rk4(f,x_0,x_0+h,y_0,1)
			rk_file = file_open('output/rk4.out')
				read(rk_file,*)x,Y_values(1)
				f_values(1) = f(x,Y_values(1))
				read(rk_file,*)x,Y_values(2)
				f_values(2) = f(x,Y_values(2))
			close(rk_file)
		endif
		
		!	USE LINEAR MULTISTEP TO COMPLETE
		error = 0.0_wp
		select case (order)
			case default ! order = 4
				do i=1,2
					write(output_file,*)(x_0 + h*real(i-1,wp)),Y_values(i)
				enddo
				do i=2,n
					x = x_0 + h*real(i,wp)
					y = Y_values(i)
					f_values(i) = f(x,y)
					!Y_values(i+1) = -4.0_wp*Y_values(i) + 5.0_wp*Y_values(i-1) + h*4.0_wp*f_values(i) + h*2.0_wp*f_values(i-1)
					Y_values(i+1) = Y_values(i) + 1.5_wp*h*f_values(i) - 0.5_wp*h*f_values(i-1)
					tmp = abs(Y_values(i+1)-P(x))
					if (tmp > error) error = tmp
					write(output_file,*)x,Y_values(i+1),tmp
				enddo
		end select
		linear_multistep = error
		
		!	CLOSE OUTPUT FILE
		close(output_file)
	end function linear_multistep
	
	function rk4(f,x_0,x_f,y_0,n)
		!	--------------------------------
		!	This function performs a RK4
		!	method to approximate y.
		!	--------------------------------
		
		!	INPUTS
		real(wp),external			::	f
		real(wp)					::	x_0,x_f
		real(wp)					::	y_0
		integer						::	n
		
		!	OUTPUTS
		real(wp)					::	rk4
		
		!	INTERNAL
		integer						::	output_file
		integer						::	i
		real(wp)					::	h
		real(wp)					::	x,y
		real(wp)					::	k1,k2,k3,k4
		
		!	OPEN OUTPUT FILE
		call create_dir('output')
		output_file = file_open('output/rk4.out')
		
		!	CALCULATE
		x = x_0
		y = y_0
		h = (x_f-x_0)/real(n,wp)
		write(output_file,*)x,y
		do i=1,n
			k1 = f(x,y)
			k2 = f(x+0.5_wp*h,y+0.5_wp*h*k1)
			k3 = f(x+0.5_wp*h,y+0.5_wp*h*k2)
			k4 = f(x+h,y+h*k3)
			x = x_0 + h*real(i,wp)
			y = y + (h/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
			write(output_file,*)x,y
		enddo
		rk4 = y
		
		!	CLOSE OUTPUT FILE
		close(output_file)
	end function rk4





end module