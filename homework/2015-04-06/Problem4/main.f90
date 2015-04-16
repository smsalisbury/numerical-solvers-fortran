program problem4
use config
use gnuplot
implicit none

! DATA DICTIONARY
real(wp)									::	y1,y2,t,y1old,y2old,told
real(wp)									::	y1tmp,y2tmp
real(wp)									::	y10,y20,t0,tf
real(wp),dimension(:),allocatable			::	y1_values,y2_values,t_values
real(wp)									::	alpha,beta,beta_c
real(wp)									::	h
integer										::	n,i

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

!	EXPLICIT EULER
y1 = y10
y2 = y20
t = t0
y1old = y1
y2old = y2
told = t
do i=1,n
	t = t + h
	y1 = y1old + h*(alpha - y1old - (4.0_wp*y1old*y2old/(1.0_wp+y1old*y1old)))
	y2 = y2old + h*(beta*y1old*(1.0_wp - y2old/(1.0_wp + y1old*y1old)))

	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
call scatter_plot2D(t_values,y1_values,1.0_wp,'euler_1.png','t','y1')
call scatter_plot2D(y1_values,y2_values,1.0_wp,'euler_2.png','y1','y2')

!	RK-2
y1 = y10
y2 = y20
t = t0
y1old = y1
y2old = y2
told = t
do i=1,n
	t = t + h
	y1tmp = y1old + h*(alpha - y1old - (4.0_wp*y1old*y2old/(1.0_wp+y1old*y1old)))
	y2tmp = y2old + h*(beta*y1old*(1.0_wp - y2old/(1.0_wp + y1old*y1old)))
	y1 = y1old + h*(alpha - y1tmp - (4.0_wp*y1tmp*y2tmp/(1.0_wp+y1tmp*y1tmp)))
	y2 = y2old + h*(beta*y1tmp*(1.0_wp - y2tmp/(1.0_wp + y1tmp*y1tmp)))
	
	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
call scatter_plot2D(t_values,y1_values,1.0_wp,'rk2_1.png','t','y1')
call scatter_plot2D(y1_values,y2_values,1.0_wp,'rk2_2.png','y1','y2')

!	Trapezoidal
y1 = y10
y2 = y20
t = t0
y1old = y1
y2old = y2
told = t
do i=1,n
	t = t + h
	y1tmp = y1old + h*(alpha - y1old - (4.0_wp*y1old*y2old/(1.0_wp+y1old*y1old)))
	y2tmp = y2old + h*(beta*y1old*(1.0_wp - y2old/(1.0_wp + y1old*y1old)))
	y1 = y1old + 0.5_wp*h*((alpha - y1tmp - (4.0_wp*y1tmp*y2tmp/(1.0_wp+y1tmp*y1tmp)))+ &
		(alpha - y1old - (4.0_wp*y1old*y2old/(1.0_wp+y1old*y1old))))
	y2 = y2old + 0.5_wp*h*((beta*y1tmp*(1.0_wp - y2tmp/(1.0_wp + y1tmp*y1tmp)))+ &
		(beta*y1old*(1.0_wp - y2old/(1.0_wp + y1old*y1old))))
	
	t_values(i) = t
	y1_values(i) = y1
	y2_values(i) = y2
	
	y1old = y1
	y2old = y2
	told = t
enddo
call scatter_plot2D(t_values,y1_values,1.0_wp,'trap_1.png','t','y1')
call scatter_plot2D(y1_values,y2_values,1.0_wp,'trap_2.png','y1','y2')

deallocate(y1_values)
deallocate(y2_values)
deallocate(t_values)
end program problem4