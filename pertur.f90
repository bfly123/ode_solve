subroutine pertur(y,f,ld,dx,x,per_k)
implicit none
integer i,j,k
	  integer nf
double precision y,f,per_k
double precision a1,a2,a3,a4,a5,a6
double precision g,g1,g2,g3,g4,g5,g6
double precision y1,y2,y3,y4,y5,y6
double precision tmp
double precision ss,dq,dx,x
!*****************F******************
!write(*,*) G
g=10-(10+x)*exp(-x)
g1=(9+x)*exp(-x)
g2=-(8+x)exp(-x)
g3=(7+x)exp(-x)
g4=-(6+x)exp(-x)
g5=(5+x)exp(-x)
g6=-(4+x)exp(-x)

 tmp=(y-g)
 f=g1+ld*tmp
 y1=f
 y2=g2+ld**2*tmp
 y3=g3+ld**3*tmp
 y4=g4+ld**4*tmp
 y5=g5+ld**5*tmp
 y6=g6+ld**6*tmp

 a1=-y2/2/y1
 a2=-y3/6/y1-a1*y2/2/y1
 a3=-y4/24/y1-a1*y3/6/y1-a2*y2/2/y1
 a4=-y5/120/y1-a1*y4/24/y1-a2*y3/6/y1-a3*y2/2/y1

	per_k=1+a1*dx+a2*dx**2+a3*dx**3

  enddo
end

