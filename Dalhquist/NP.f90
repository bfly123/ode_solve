subroutine pertur(x,y,f,dx,per_k,n_method,n_fun)
implicit none
integer i,j,k
integer n_fun  !Function number
integer n_order  !method order 
integer n_method !method  number 1-7 NP method 11-71 transformed NP
double precision  x,fy,dx,fymax,dx1 !fy=y''/y' 
double precision y(n_fun),f(n_fun),per_k(n_fun)
double precision yp(7,n_fun) !derivatives
double precision a1,a2,a3,a4,a5,a6,a7 !perturbation coefficients
double precision b1,b2,b3,b4,b5,b6,b7 !transformed coefficients
double precision c1,c2,c3,c4,c5       !transformed coefficients

double precision tmp
!*****************F******************
!write(*,*) G

select case(n_method)
case(1)
n_order=1
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
per_k(:)=1
case(2)
n_order=2
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
	per_k(i)=1+dx*a1
enddo
case(3)
n_order=3
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
	per_k(i)=1+dx*a1+dx**2*a2
enddo
case(4)
n_order=4
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
	per_k(i)=1+dx*a1+dx**2*a2+dx**3*a3
enddo
case(5)
n_order=5
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
    a4=-yp(5,i)/120/yp(1,i)-a1*yp(4,i)/24/yp(1,i)-a2*yp(3,i)/6/yp(1,i)-a3*yp(2,i)/2/yp(1,i)
	per_k(i)=1+dx*a1+dx**2*a2+dx**3*a3+dx**4*a4
enddo

case(6)
n_order=6
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
    a4=-yp(5,i)/120/yp(1,i)-a1*yp(4,i)/24/yp(1,i)-a2*yp(3,i)/6/yp(1,i)-a3*yp(2,i)/2/yp(1,i)
    a5=-yp(6,i)/720/yp(1,i)-a1*yp(5,i)/120/yp(1,i)-a2*yp(4,i)/24/yp(1,i)-a3*yp(3,i)/6/yp(1,i)-a4*yp(2,i)/2/yp(1,i)
	per_k(i)=1+dx*a1+dx**2*a2+dx**3*a3+dx**4*a4+dx**5*a5
enddo

case(7)
n_order=7
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
    a4=-yp(5,i)/120/yp(1,i)-a1*yp(4,i)/24/yp(1,i)-a2*yp(3,i)/6/yp(1,i)-a3*yp(2,i)/2/yp(1,i)
    a5=-yp(6,i)/720/yp(1,i)-a1*yp(5,i)/120/yp(1,i)-a2*yp(4,i)/24/yp(1,i)-a3*yp(3,i)/6/yp(1,i)-a4*yp(2,i)/2/yp(1,i)
    a6=-yp(7,i)/5040/yp(1,i)-a1*yp(6,i)/720/yp(1,i)-a2*yp(5,i)/120/yp(1,i)-a3*yp(4,i)/24/yp(1,i)&
		-a4*yp(3,i)/6/yp(1,i)-a5*yp(2,i)/2/yp(1,i)
	per_k(i)=1+dx*a1+dx**2*a2+dx**3*a3+dx**4*a4+dx**5*a5+dx**6*a6
enddo

! a3=-y4/24/y1-a11*y3/6/y1-a12*y2/2/y1
! a4=-y5/120/y1-a11*y4/24/y1-a12*y3/6/y1-a13*y2/2/y1
! a5=-y6/720/y1-a11*y5/120/y1-a12*y4/24/y1-a13*y3/6/y1-a14*y2/2/y1
! a6=-y7/5040/y1-a11*y6/720/y1-a12*y5/120/y1-a13*y4/24/y1-a14*y3/6/y1-a15*y2/2/y1
!
case(31)
!!**************3rdTransform******************
!
n_order=3
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
fymax=0
do i=1,n_fun
	fy=yp(2,i)/yp(1,i)
	enddo

do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
	fy=yp(2,i)/yp(1,i)
	b2=a2/(1+a1/fy)
	b1=a1-a2/(a1+fy)
	per_k(i)=(1+b1*dx +b2*dx**2)/(1-b2/fy*dx)
enddo
case(41)
!***********4thTransform********************
n_order=4
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
	fy=yp(2,i)/yp(1,i)
	b1=-5*a1
	b2=a2-6*a1**2-(a3-6*a1*a2)/(fy+a1)
	b3=(a3-6*a1*a2)/(1+a1/fy)
	c1=-6*a1
	c2=-b3/fy
	per_k(i)=(1+b1*dx +b2*dx**2+b3*dx**3)/(1+c1*dx+c2*dx**2) 
enddo
case(51)
!*********5thorderTransform******************************
n_order=5
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
    a4=-yp(5,i)/120/yp(1,i)-a1*yp(4,i)/24/yp(1,i)-a2*yp(3,i)/6/yp(1,i)-a3*yp(2,i)/2/yp(1,i)
	fy=yp(2,i)/yp(1,i)
	c1=-6*a1
	c2=-144*a2
	b1=a1+c1
	b2=a2+a1*c1+c2
	b4=(a4+a3*c1+a2*c2)/(1+a1/fy)
	c3=-b4/fy
	b3=a3+a2*c1+a1*c2+c3

	per_k(i)=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4)/(1+c1*dx+c2*dx**2+c3*dx**3) 
	enddo

case(61)
n_order=6
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
    a4=-yp(5,i)/120/yp(1,i)-a1*yp(4,i)/24/yp(1,i)-a2*yp(3,i)/6/yp(1,i)-a3*yp(2,i)/2/yp(1,i)
    a5=-yp(6,i)/720/yp(1,i)-a1*yp(5,i)/120/yp(1,i)-a2*yp(4,i)/24/yp(1,i)-a3*yp(3,i)/6/yp(1,i)-a4*yp(2,i)/2/yp(1,i)
	fy=yp(2,i)/yp(1,i)

c1=-6*a1
c2=144*a2
c3=600*a3

b1=a1+c1
b2=a2+c1*a1+c2
b3=a3+a2*c1+a1*c2+c3
b5=(a5+a4*c1+a3*c2+a2*c3)/(1+a1/fy)
c4=-b5/fy
b4=a4+a3*c1+a2*c2+a1*c3+c4
per_k(i)=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4) 
enddo

case(71)
n_order=7
call f_derivative(n_order,n_fun,y,yp)
f(:)=yp(1,:)
do i = 1,n_fun
	a1=-yp(2,i)/2/yp(1,i)
    a2=-yp(3,i)/6/yp(1,i)-a1*yp(2,i)/2/yp(1,i)
    a3=-yp(4,i)/24/yp(1,i)-a1*yp(3,i)/6/yp(1,i)-a2*yp(2,i)/2/yp(1,i)
    a4=-yp(5,i)/120/yp(1,i)-a1*yp(4,i)/24/yp(1,i)-a2*yp(3,i)/6/yp(1,i)-a3*yp(2,i)/2/yp(1,i)
    a5=-yp(6,i)/720/yp(1,i)-a1*yp(5,i)/120/yp(1,i)-a2*yp(4,i)/24/yp(1,i)-a3*yp(3,i)/6/yp(1,i)-a4*yp(2,i)/2/yp(1,i)
    a6=-yp(7,i)/5040/yp(1,i)-a1*yp(6,i)/720/yp(1,i)-a2*yp(5,i)/120/yp(1,i)-a3*yp(4,i)/24/yp(1,i)&
		-a4*yp(3,i)/6/yp(1,i)-a5*yp(2,i)/2/yp(1,i)
	fy=yp(2,i)/yp(1,i)

c1=-6*a1
c2=144*a2
c3=-144*a3
c4=1440*a4

b1=a1+c1
b2=a2+c1*a1+c2
b3=a3+a2*c1+a1*c2+c3
b4=a4+a3*c1+a2*c2+a1*c3+c4
b6=(a6+a5*c1+a4*c2+a3*c3+a2*c4)/(1+a1/fy)
c5=-b6/fy
b5=a5+a4*c1+a3*c2+a2*c3+a1*c4+c5

per_k(i)=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5+b6*dx**6)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4+c5*dx**5) 
enddo

endselect

do i=1,n_fun
if (abs(yp(1,i)).le.1.d-16) then
	per_k(i)=1
endif
enddo

end


