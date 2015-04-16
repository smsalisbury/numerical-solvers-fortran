program problem1
use config
use gnuplot
implicit none

! DATA DICTIONARY
real(wp)							::	theta0,thetaf,h
real(wp)							::	x0,y0
real(wp)							::	x,y,xold,yold
real(wp),dimension(:),allocatable	::	x_values,y_values	
integer								::	n,i

h = 0.02_wp
theta0 = 0.0_wp
thetaf = 120.0_wp
x0 = 1.0_wp
y0 = 0.0_wp

n = int((thetaf-theta0)/h)

allocate(x_values(n+1))
allocate(y_values(n+1))

x = x0
y = y0
xold = x
yold = y
!	FORWARD EULER
do i = 1,n
	x = xold - h*yold
	y = yold + h*xold
	x_values(i) = x
	y_values(i) = y
	xold = x
	yold = y
enddo
call scatter_plot2D(x_values,y_values,1.0_wp,'forward.png')

x = x0
y = y0
xold = x
yold = y
!	BACKWARD EULER
do i = 1,n
	x = xold - h*yold
	y = yold + h*x
	x = xold - h*y
	x_values(i) = x
	y_values(i) = y
	xold = x
	yold = y
enddo
call scatter_plot2D(x_values,y_values,1.0_wp,'backward.png')

x = x0
y = y0
xold = x
yold = y
!	TRAPEZOIDAL
do i = 1,n
	y = yold + h*xold
	x = xold - 0.5_wp*h*(yold + y)
	y = yold + 0.5_wp*h*(xold + x)
	x_values(i) = x
	y_values(i) = y
	xold = x
	yold = y
enddo
call scatter_plot2D(x_values,y_values,1.0_wp,'trapezoid.png')

deallocate(x_values)
deallocate(y_values)
end program problem1