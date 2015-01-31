program main
implicit none

! DATA DICTIONARY
real::a,b,c,x,v,u
real,dimension(2)::p,q,r
integer::i

! POINTS
p = (/ 0.0, 1.0 /)
q = (/ 1.0, 0.9 /)
r = (/ 3.0, 0.5 /)

! COEFFICIENTS
! ln(gamma_0)
a = q(1)*r(1)*log(p(2))/((p(1)-q(1))*(p(1)-r(1))) - p(1)*r(1)*log(q(2))/((p(1)-q(1))*(q(1)-r(1))) &
      + p(1)*q(1)*log(r(2))/((p(1)-r(1))*(q(1)-r(1)))
! gamma_1
b = -(q(1)+r(1))*log(p(2))/((p(1)-q(1))*(p(1)-r(1))) + (p(1)+r(1))*log(q(2))/((p(1)-q(1))*(q(1)-r(1))) &
      + (p(1)+q(1))*log(r(2))/((p(1)-r(1))*(r(1)-q(1)))
! gamma_2
c = log(p(2))/((p(1)-q(1))*(p(1)-r(1))) - log(q(2))/((p(1)-q(1))*(q(1)-r(1))) + log(r(2))/((p(1)-r(1))*(q(1)-r(1)))

write(*,*)exp(a),b,c

end program main