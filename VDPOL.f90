program main
!_________________________________________________________________
!|****************************************************************|
!|*******This is to solve  y'=F(y) by pertrubation method*******  |
!|**************2016年 01月 15日 星期五 15:50:32 CST**************|
!|****************************************************************|
!|________________________________________________________________|

implicit none
integer i,j,k,ix,nx
double precision x,dx,xx
double precision y,z,f1,f2,per_k1,per_k2
double precision ld,exact
double precision err,err1


y=2.d0
z=0.d0
xx=11.d0
dx=0.1d-2

nx=0
x=0

      open(1,file='result.dat',status='unknown')
	  do ix=1,int(xx/dx)+10
	     if (x>=xx) then
			 exit
		 else if (x+dx>xx)then
			 dx=xx-x
		 endif
         call pertur(y,z,f1,f2,ld,dx,x,per_k1,per_k2)
		    y=y+dx*f1/per_k1 
		    z=z+dx*f2/per_k2 
			x=x+dx

        write(*,*)x, y,z
        write(1,*)x, y,z
	enddo
	close(1)

	!err=dabs((exact-y)/exact)
!
!	dx=0.1d0
!      nx=0
!	  x=0
!	y=0.d0
!
!	  do ix=1,int(xx/dx)+10
!	     if (x>=xx) then
!			 exit
!		 else if (x+dx>xx)then
!			 dx=xx-x
!		 endif
!         call pertur(y,f,ld,dx,x,per_K)
!		    y=y+dx*f /per_k 
!			x=x+dx
!		!	exact=10-(10+x)*exp(-x)
!	enddo
!
!      !open(1,file='result.dat',status='unknown')
!      !  write(*,*)x, dabs(exact-y)/exact
!		err1=dabs(exact-y)/exact
!
!		!write(*,*)lde-2, dlog(err)/ dlog(10.d0)  !, dlog(err/err1)/dlog( 2.d0 )
!		!write(*,*)lde+log(0.2), err,err1
!		write(*,*)lde+log(0.2), dlog(err)/dlog(10.d0), dlog(err/err1)/dlog( 2.d0 )
!		!write(*,*)lde-2, dlog( err/err1 )/dlog( 2.d0 )
!		lde=lde+0.1d0
!	enddo
!     !close(1)
  end

subroutine pertur(y,z,f1,f2,ld,dx,x,per_k1,per_k2)
implicit none
integer i,j,k
	  integer nf
double precision y,z,f1,f2,per_k1,per_k2,fy1,fy2
double precision a11,a12,a13,a14,a15,a16,a17
double precision a21,a22,a23,a24,a25,a26,a27
double precision b1,b2,b3,b4,b5,b6,b7
double precision g,g1,g2,g3,g4,g5,g6,g7
double precision y1,y2,y3,y4,y5,y6,y7
double precision z1,z2,z3,z4,z5,z6,z7
double precision c1,c2,c3,c4,c5
double precision tmp
double precision ss,dq,dx,x,ld
integer method
!*****************F******************
!write(*,*) G

method=61
ss=1.d-2

y1=z
z1=((1-y**2)*z-y)/ss

f1=y1
f2=z1

y2=z1
z2=((1-y**2)*z1-2*y*z*y1-y1)/ss

fy1=y2/y1
fy2=z2/z1

y3=z2
z3=-((y**2 - 1)*z2 + 2*y*z*y2 + 4*y*y1*z1 + 2*z*y1**2 + y2)/ss

y4=z3
z4=-((y**2 - 1)*z3 + 2*y*z*y3 + 6*y*y1*z2 + 6*y*y2*z1 + 6*z*y1*y2&
   	+ 6*y1**2*z1 + y3)/ss

y5=z4
z5=-((y**2 - 1)*z4 + 2*y*z*y4 + 8*y*y1*z3 + 12*y*y2*z2 + 8*y*y3*z1&
   	+ 8*z*y1*y3 + 6*z*y2**2 + 12*y1**2*z2 + 24*y1*y2*z1 + y4)/ss
y6=z5
z6=-((y**2 - 1)*z5 + 2*y*z*y5 + 10*y*y1*z4 + 20*y*y2*z3 + 20*y*y3*z2&
   	+ 10*y*y4*z1 + 10*z*y1*y4 + 20*z*y2*y3 + 20*y1**2*z3 + 60*y1*y2*z2&
   	+ 40*y1*y3*z1 + 30*y2**2*z1 + y5)/ss
y7=z6
z7=-((y**2 - 1)*z6 + 2*y*z*y6 + 12*y*y1*z5 + 30*y*y2*z4 + 40*y*y3*z3&
   	+ 30*y*y4*z2 + 12*y*y5*z1 + 12*z*y1*y5 + 30*z*y2*y4 + 20*z*y3**2 &
	+ 30*y1**2*z4 + 120*y1*y2*z3 + 120*y1*y3*z2 + 60*y1*y4*z1 + 90*y2**2*z2 &
   	+ 120*y2*y3*z1 + y6)/ss

 a11=-y2/2/y1
 a12=-y3/6/y1-a11*y2/2/y1
 a13=-y4/24/y1-a11*y3/6/y1-a12*y2/2/y1
 a14=-y5/120/y1-a11*y4/24/y1-a12*y3/6/y1-a13*y2/2/y1
 a15=-y6/720/y1-a11*y5/120/y1-a12*y4/24/y1-a13*y3/6/y1-a14*y2/2/y1
 a16=-y7/5040/y1-a11*y6/720/y1-a12*y5/120/y1-a13*y4/24/y1-a14*y3/6/y1-a15*y2/2/y1

 a21=-z2/2/z1
 a22=-z3/6/z1-a21*z2/2/z1
 a23=-z4/24/z1-a21*z3/6/z1-a22*z2/2/z1
 a24=-z5/120/z1-a21*z4/24/z1-a22*z3/6/z1-a23*z2/2/z1
 a25=-z6/720/z1-a21*z5/120/z1-a22*z4/24/z1-a23*z3/6/z1-a24*z2/2/z1
 a26=-z7/5040/z1-a21*z6/720/z1-a22*z5/120/z1-a23*z4/24/z1-a24*z3/6/z1-a25*z2/2/z1

select case(method)
case(2)
per_k1=1+a11*dx  
per_k2=1+a21*dx  
case(3)
per_k1=1+a11*dx+a12*dx**2
per_k2=1+a21*dx+a22*dx**2
case(4)
per_k1=1+a11*dx+a12*dx**2+a13*dx**3 
per_k2=1+a21*dx+a22*dx**2+a23*dx**3 
case(5)
per_k1=1+a11*dx+a12*dx**2+a13*dx**3+a14*dx**4 
per_k2=1+a21*dx+a22*dx**2+a23*dx**3+a24*dx**4 
case(6)
per_k1=1+a11*dx+a12*dx**2+a13*dx**3+a14*dx**4+a15*dx**5 
per_k2=1+a21*dx+a22*dx**2+a23*dx**3+a24*dx**4+a25*dx**5 
case(7)
per_k1=1+a11*dx+a12*dx**2+a13*dx**3+a14*dx**4+a15*dx**5+a16*dx**6
per_k2=1+a21*dx+a22*dx**2+a23*dx**3+a24*dx**4+a25*dx**5+a26*dx**6

case(31)
!!**************3rdTransform******************
!
b2=a12/(1+a11/fy1)
b1=a11-a12/(a11+fy1)
per_k1=(1+b1*dx +b2*dx**2)/(1-b2/fy1*dx)

b2=a22/(1+a21/fy2)
b1=a21-a22/(a21+fy2)
per_k2=(1+b1*dx +b2*dx**2)/(1-b2/fy2*dx)

case(41)
!***********4thTransform********************
b1=-5*a11
b2=a12-6*a11**2-(a13-6*a11*a12)/(fy1+a11)
b3=(a13-6*a11*a12)/(1+a11/fy1)
c1=-6*a11
c2=-b3/fy1
per_k1=(1+b1*dx +b2*dx**2+b3*dx**3)/(1+c1*dx+c2*dx**2) 

b1=-5*a21
b2=a22-6*a21**2-(a23-6*a21*a22)/(fy2+a21)
b3=(a23-6*a21*a22)/(1+a21/fy2)
c1=-6*a21
c2=-b3/fy2
per_k2=(1+b1*dx +b2*dx**2+b3*dx**3)/(1+c1*dx+c2*dx**2) 
!
case(51)
!*********5thorderTransform******************************
c1=-6*a11
c2=-144*a12

b1=a11+c1
b2=a12+a11*c1+c2
b4=(a14+a13*c1+a12*c2)/(1+a11/fy1)
c3=-b4/fy1
b3=a13+a12*c1+a11*c2+c3

per_k1=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4)/(1+c1*dx+c2*dx**2+c3*dx**3) 

c1=-6*a21
c2=-144*a22

b1=a21+c1
b2=a22+a21*c1+c2
b4=(a24+a23*c1+a22*c2)/(1+a21/fy2)
c3=-b4/fy2
b3=a23+a22*c1+a21*c2+c3

per_k2=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4)/(1+c1*dx+c2*dx**2+c3*dx**3) 

case(61)
c1=-6*a11
c2=144*a12
c3=600*a13

b1=a11+c1
b2=a12+c1*a11+c2
b3=a13+a12*c1+a11*c2+c3
b5=(a15+a14*c1+a13*c2+a12*c3)/(1+a11/fy1)
c4=-b5/fy1
b4=a14+a13*c1+a12*c2+a11*c3+c4
per_k1=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4) 

c1=-6*a21
c2=144*a22
c3=600*a23

b1=a21+c1
b2=a22+c1*a21+c2
b3=a23+a22*c1+a21*c2+c3
b5=(a25+a24*c1+a23*c2+a22*c3)/(1+a21/fy2)
c4=-b5/fy2
b4=a24+a23*c1+a22*c2+a21*c3+c4

per_k2=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4) 

case(71)
c1=-6*a11
c2=144*a12
c3=-144*a13
c4=1440*a14

b1=a11+c1
b2=a12+c1*a11+c2
b3=a13+a12*c1+a11*c2+c3
b4=a14+a13*c1+a12*c2+a11*c3+c4
b6=(a16+a15*c1+a14*c2+a13*c3+a12*c4)/(1+a11/fy1)
c5=-b6/fy1
b5=a15+a14*c1+a13*c2+a12*c3+a11*c4+c5

per_k1=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5+b6*dx**6)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4+c5*dx**5) 

c1=-6*a21
c2=144*a22
c3=-144*a23
c4=1440*a24

b1=a21+c1
b2=a22+c1*a21+c2
b3=a23+a22*c1+a21*c2+c3
b4=a24+a23*c1+a22*c2+a21*c3+c4
b6=(a26+a25*c1+a24*c2+a23*c3+a22*c4)/(1+a21/fy2)
c5=-b6/fy2
b5=a25+a24*c1+a23*c2+a22*c3+a21*c4+c5

per_k2=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5+b6*dx**6)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4+c5*dx**5) 

endselect

if (abs(y1).le.1.d-16)then
	per_k1=1
endif
if (abs(z1).le.1.d-16)then
	per_k2=1
endif

end



