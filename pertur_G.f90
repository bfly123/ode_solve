subroutine pertur_G(y,g,nf)
implicit none
integer i,j,k
double precision ss,dq
integer nf
double precision y(nf)
double precision g(nf)
double precision f(nf)
double precision y_tmp(nf)
double precision f_tmp(nf)
	

call sub_f(y,f,nf)
y_tmp=y
	g=0
ss=1.d-3

do j=1,nf
	  y_tmp(j)=y_tmp(j)+ss
	  call sub_f(y_tmp,f_tmp,nf)
	  dq=ss
!	  if (abs(y(j)).le.1.d-12) dq=1
	  do i=1,nf
		G(i)=G(i)+(f_tmp(i)-f(i))/dq*f(j)
	enddo
	  y_tmp(j)=y(j)
enddo
end
