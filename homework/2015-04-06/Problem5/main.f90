program problem5
use config
use gnuplot
implicit none

! DATA DICTIONARY
real(wp)									::	y1,y2,t,y1old,y2old,told
real(wp),dimension(:),allocatable			::	f1,f2 ! Most recent first
real(wp)									::	k1,k2,k3,k4
real(wp)									::	y10,y20,t0,tf
real(wp),dimension(:),allocatable			::	y1_values,y2_values,t_values
real(wp)									::	alpha,beta,beta_c
real(wp)									::	h
integer										::	n,i,steps,j

alpha = 10.0_wp
beta = 3.6_wp
beta_c = (3.0_wp*alpha)/(5.0_wp) - (25.0_wp/alpha)

y10 = 0.0_wp
y20 = 2.0_wp
t0 = 0.0_wp
tf = 100.0_wp
h = 0.01_wp

n = int((tf - t0)/h)

allocate(y1_values(n))
allocate(y2_values(n))
allocate(t_values(n))

!	TWO STEP
steps = 2
y1 = y10
y2 = y20
t = t0
y1old = y1
y2old = y2
told = t
allocate(f1(steps))
allocate(f2(steps))
f1 = 0.0_wp
f2 = 0.0_wp
do i=1,steps
	do j=steps,2,-1
		f1(j) = f1(j-1)
		f2(j) = f2(j-1)
	enddo
	f1(1) = calcf1(y1old,y2old)
	f2(1) = calcf2(y1old,y2old)
	
	t = t + h
	!	RK-4
	!	--	y1
	k1 = h*calcf1(y1old,y2old)
	k2 = h*calcf1(y1old+0.5_wp*k1,y2old)
	k3 = h*calcf1(y1old+0.5_wp*k2,y2old)
	k4 = h*calcf1(y1old+k3,y2old)
	y1 = y1old + (1.0_wp/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
	!	--	y1
	k1 = h*calcf2(y1old,y2old)
	k2 = h*calcf2(y1old,y2old+0.5_wp*k1)
	k3 = h*calcf2(y1old,y2old+0.5_wp*k2)
	k4 = h*calcf2(y1old,y2old+k3)
	y2 = y2old + (1.0_wp/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
	
	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
do i=steps+1,n
	do j=steps,2,-1
		f1(j) = f1(j-1)
		f2(j) = f2(j-1)
	enddo
	f1(1) = calcf1(y1old,y2old)
	f2(1) = calcf2(y1old,y2old)
	
	t = t + h
	y1 = y1old + h/2.0_wp*(3.0_wp*f1(1) - f1(2))
	y2 = y2old + h/2.0_wp*(3.0_wp*f2(1) - f2(2))

	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
deallocate(f1)
deallocate(f2)
call scatter_plot2D(t_values,y1_values,1.0_wp,'adams2_1.png','t','y1')
call scatter_plot2D(y1_values,y2_values,1.0_wp,'adams2_2.png','y1','y2')

!	FOUR STEP
steps = 4
y1 = y10
y2 = y20
t = t0
y1old = y1
y2old = y2
told = t
allocate(f1(steps))
allocate(f2(steps))
f1 = 0.0_wp
f2 = 0.0_wp
do i=1,steps
	do j=steps,2,-1
		f1(j) = f1(j-1)
		f2(j) = f2(j-1)
	enddo
	f1(1) = calcf1(y1old,y2old)
	f2(1) = calcf2(y1old,y2old)
	
	t = t + h
	!	RK-4
	!	--	y1
	k1 = h*calcf1(y1old,y2old)
	k2 = h*calcf1(y1old+0.5_wp*k1,y2old)
	k3 = h*calcf1(y1old+0.5_wp*k2,y2old)
	k4 = h*calcf1(y1old+k3,y2old)
	y1 = y1old + (1.0_wp/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
	!	--	y1
	k1 = h*calcf2(y1old,y2old)
	k2 = h*calcf2(y1old,y2old+0.5_wp*k1)
	k3 = h*calcf2(y1old,y2old+0.5_wp*k2)
	k4 = h*calcf2(y1old,y2old+k3)
	y2 = y2old + (1.0_wp/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
	
	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
do i=steps+1,n
	do j=steps,2,-1
		f1(j) = f1(j-1)
		f2(j) = f2(j-1)
	enddo
	f1(1) = calcf1(y1old,y2old)
	f2(1) = calcf2(y1old,y2old)
	
	t = t + h
	y1 = y1old + h/24.0_wp*(55.0_wp*f1(1) - 59.0_wp*f1(2) + 37.0_wp*f1(3) - 9.0_wp*f1(4))
	y2 = y2old + h/24.0_wp*(55.0_wp*f2(1) - 59.0_wp*f2(2) + 37.0_wp*f2(3) - 9.0_wp*f2(4))

	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
deallocate(f1)
deallocate(f2)
call scatter_plot2D(t_values,y1_values,1.0_wp,'adams4_1.png','t','y1')
call scatter_plot2D(y1_values,y2_values,1.0_wp,'adams4_2.png','y1','y2')

!	FIVE STEP
steps = 5
y1 = y10
y2 = y20
t = t0
y1old = y1
y2old = y2
told = t
allocate(f1(steps))
allocate(f2(steps))
f1 = 0.0_wp
f2 = 0.0_wp
do i=1,steps
	do j=steps,2,-1
		f1(j) = f1(j-1)
		f2(j) = f2(j-1)
	enddo
	f1(1) = calcf1(y1old,y2old)
	f2(1) = calcf2(y1old,y2old)
	
	t = t + h
	!	RK-4
	!	--	y1
	k1 = h*calcf1(y1old,y2old)
	k2 = h*calcf1(y1old+0.5_wp*k1,y2old)
	k3 = h*calcf1(y1old+0.5_wp*k2,y2old)
	k4 = h*calcf1(y1old+k3,y2old)
	y1 = y1old + (1.0_wp/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
	!	--	y1
	k1 = h*calcf2(y1old,y2old)
	k2 = h*calcf2(y1old,y2old+0.5_wp*k1)
	k3 = h*calcf2(y1old,y2old+0.5_wp*k2)
	k4 = h*calcf2(y1old,y2old+k3)
	y2 = y2old + (1.0_wp/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
	
	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
do i=steps+1,n
	do j=steps,2,-1
		f1(j) = f1(j-1)
		f2(j) = f2(j-1)
	enddo
	f1(1) = calcf1(y1old,y2old)
	f2(1) = calcf2(y1old,y2old)
	
	t = t + h
	y1 = y1old + h/720.0_wp*(1901.0_wp*f1(1) - 2774.0_wp*f1(2) + 2616.0_wp*f1(3) - 1274.0_wp*f1(4) + 251.0_wp*f1(5))
	y2 = y2old + h/720.0_wp*(1901.0_wp*f2(1) - 2774.0_wp*f2(2) + 2616.0_wp*f2(3) - 1274.0_wp*f2(4) + 251.0_wp*f2(5))

	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
deallocate(f1)
deallocate(f2)
call scatter_plot2D(t_values,y1_values,1.0_wp,'adams5_1.png','t','y1')
call scatter_plot2D(y1_values,y2_values,1.0_wp,'adams5_2.png','y1','y2')

deallocate(y1_values)
deallocate(y2_values)
deallocate(t_values)

contains
	function calcf1(y1,y2)
		real(wp)	::	y1,y2
		real(wp)	:: 	calcf1
		
		calcf1 = alpha - y1 - (4.0_wp*y1*y2)/(1.0_wp+y1*y1)
	end function
	
	function calcf2(y1,y2)
		real(wp)	::	y1,y2
		real(wp)	:: 	calcf2
		
		calcf2 = beta*y1*(1.0_wp - y2/(1.0_wp+y1*y1))
	end function
end program problem5