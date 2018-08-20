program main
!_________________________________________________________________
!|****************************************************************|
!|*******This is to solve  y'=F(y) by pertrubation method*******  |
!|**************2016年 01月 15日 星期五 15:50:32 CST**************|
!|****************************************************************|
!|________________________________________________________________|

implicit none
integer(kind=8) i,j,k,ix,nx
integer n_fun,n_method
parameter(n_fun=1)
double precision x,dx,xx
double precision y(n_fun),per_k(n_fun),f(n_fun)
double precision exact
double precision err,err1
double precision fymax


!init
y(1)=1.d-16
xx=0.1d0
dx=1.0d-3
nx=10000000000

n_method=2

	open(1,file='result.dat',status='unknown')
	x=0
	j=1
	!write(*,*) int(xx/dx)
	do j=1,nx
	call pertur(x,y,f,dx,per_k,n_method,n_fun)
     if (x>=xx) then
		 exit
	 else if (x+dx>xx)then
		 dx=xx-x
	 endif
         y(:)=y(:)+dx*f(:)/per_k(:)
		x=x+dx
		exact=1.d-16 *dexp(1.d3*x)
        write(*,*)x, (exact- y(:))/exact
        write(1,*)x, exact, y(:)
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


