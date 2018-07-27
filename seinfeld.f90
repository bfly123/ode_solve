program main
!_________________________________________________________________
!|****************************************************************|
!|*******This is to solve  y'=F(y) by pertrubation method*******  |
!|**************2016年 01月 15日 星期五 15:50:32 CST**************|
!|****************************************************************|
!|________________________________________________________________|

implicit none
integer i,j,k,ix,nx,nt2
double precision x,dx,xx
double precision y,y0,f,g,per_k
double precision ld,exact,lde
double precision err,err1


lde=-3
!lde=-.2-log(0.4)

do while (lde.le.4.6)


y=0.d0
xx=5.d0
dx=0.2d0
ld=-10**lde

      nx=0
	  x=0

	  do ix=1,int(xx/dx)+10
	     if (x>=xx) then
			 exit
		 else if (x+dx>xx)then
			 dx=xx-x
		 endif
         call pertur(y,f,ld,dx,x,per_K)
		    y=y+dx*f/per_k 
			x=x+dx
			exact=10-(10+x)*exp(-x)
	enddo

	err=dabs((exact-y)/exact)

        !write(*,*)x, err

	dx=0.1d0
      nx=0
	  x=0
	y=0.d0

	  do ix=1,int(xx/dx)+10
	     if (x>=xx) then
			 exit
		 else if (x+dx>xx)then
			 dx=xx-x
		 endif
         call pertur(y,f,ld,dx,x,per_K)
		    y=y+dx*f /per_k 
			x=x+dx
		!	exact=10-(10+x)*exp(-x)
	enddo

      !open(1,file='result.dat',status='unknown')
      !  write(*,*)x, dabs(exact-y)/exact
		err1=dabs(exact-y)/exact

		!write(*,*)lde-2, dlog(err)/ dlog(10.d0)  !, dlog(err/err1)/dlog( 2.d0 )
		!write(*,*)lde+log(0.2), err,err1
		write(*,*)lde+log(0.2), dlog(err)/dlog(10.d0), dlog(err/err1)/dlog( 2.d0 )
		!write(*,*)lde-2, dlog( err/err1 )/dlog( 2.d0 )
		lde=lde+0.1d0
	enddo
     !close(1)
  end


subroutine pertur(y,f,ld,dx,x,per_k)
implicit none
integer i,j,k
	  integer nf
double precision y,f,per_k
double precision a1,a2,a3,a4,a5,a6,a7
double precision b1,b2,b3,b4,b5,b6,b7
double precision g,g1,g2,g3,g4,g5,g6,g7
double precision y1,y2,y3,y4,y5,y6,y7
double precision c1,c2,c3,c4,c5
double precision tmp
double precision ss,dq,dx,x,ld
integer method
!*****************F******************
!write(*,*) G

method=71

g=10-(10+x)*exp(-x)
g1=(9+x)*exp(-x)
g2=-(8+x)*exp(-x)
g3=(7+x)*exp(-x)
g4=-(6+x)*exp(-x)
g5=(5+x)*exp(-x)
g6=-(4+x)*exp(-x)
g7=(3+x)*exp(-x)

 tmp=(y-g)
 f=g1+ld*tmp
 y1=f
 y2=g2+ld**2*tmp
 y3=g3+ld**3*tmp
 y4=g4+ld**4*tmp
 y5=g5+ld**5*tmp
 y6=g6+ld**6*tmp
 y7=g7+ld**7*tmp

 a1=-y2/2/y1
 a2=-y3/6/y1-a1*y2/2/y1
 a3=-y4/24/y1-a1*y3/6/y1-a2*y2/2/y1
 a4=-y5/120/y1-a1*y4/24/y1-a2*y3/6/y1-a3*y2/2/y1
 a5=-y6/720/y1-a1*y5/120/y1-a2*y4/24/y1-a3*y3/6/y1-a4*y2/2/y1
 a6=-y7/5040/y1-a1*y6/720/y1-a2*y5/120/y1-a3*y4/24/y1-a4*y3/6/y1-a5*y2/2/y1

select case(method)
case(2)
per_k=1+a1*dx  
case(3)
per_k=1+a1*dx+a2*dx**2
case(4)
per_k=1+a1*dx+a2*dx**2+a3*dx**3 
case(5)
per_k=1+a1*dx+a2*dx**2+a3*dx**3+a4*dx**4 
case(6)
per_k=1+a1*dx+a2*dx**2+a3*dx**3+a4*dx**4+a5*dx**5 
case(7)
per_k=1+a1*dx+a2*dx**2+a3*dx**3+a4*dx**4+a5*dx**5+a6*dx**6

case(31)
!**************3rdTransform******************

b2=a2/(1+a1/ld)
b1=a1-a2/(a1+ld)
per_k=(1+b1*dx +b2*dx**2)/(1-b2/ld*dx)

case(41)
!***********4thTransform********************
b1=-5*a1
b2=a2-6*a1**2-(a3-6*a1*a2)/(ld+a1)
b3=(a3-6*a1*a2)/(1+a1/ld)

c1=-6*a1
c2=-b3/ld

per_k=(1+b1*dx +b2*dx**2+b3*dx**3)/(1+c1*dx+c2*dx**2) 

case(51)
!*********5thorderTransform******************************
c1=-6*a1
c2=-144*a2

b1=a1+c1
b2=a2+a1*c1+c2
b4=(a4+a3*c1+a2*c2)/(1+a1/ld)
c3=-b4/ld
b3=a3+a2*c1+a1*c2+c3

per_k=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4)/(1+c1*dx+c2*dx**2+c3*dx**3) 

case(61)
c1=-6*a1
c2=144*a2
c3=600*a3

b1=a1+c1
b2=a2+c1*a1+c2
b3=a3+a2*c1+a1*c2+c3
b5=(a5+a4*c1+a3*c2+a2*c3)/(1+a1/ld)
c4=-b5/ld
b4=a4+a3*c1+a2*c2+a1*c3+c4

per_k=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4) 

case(71)
c1=-6*a1
c2=144*a2
c3=-144*a3
c4=1440*a4

b1=a1+c1
b2=a2+c1*a1+c2
b3=a3+a2*c1+a1*c2+c3
b4=a4+a3*c1+a2*c2+a1*c3+c4
b6=(a6+a5*c1+a4*c2+a3*c3+a2*c4)/(1+a1/ld)
c5=-b6/ld
b5=a5+a4*c1+a3*c2+a2*c3+a1*c4+c5

per_k=(1+b1*dx +b2*dx**2+b3*dx**3+b4*dx**4+b5*dx**5+b6*dx**6)/(1+c1*dx+c2*dx**2+c3*dx**3+c4*dx**4+c5*dx**5) 
endselect

end



