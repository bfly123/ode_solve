subroutine f_derivative(n,n_fun,y,yp)
	  integer n, n_fun
	  double precision yp(7,n_fun)
	  double precision y(n_fun)
	  double precision q
		q=1.d3
	  select case(n)
	  case(1)
		! Begin f_y 3
		yp(1,1)=q*y(1)
		! End f_y
	  case(2)
		yp(1,1)=q*y(1)
		yp(2,1)=q**2*y(1)
	  case(3)
		yp(1,1)=q*y(1)
		yp(2,1)=q**2*y(1)
		yp(3,1)=q**3*y(1)
	  case(4)
		yp(1,1)=q*y(1)
		yp(2,1)=q**2*y(1)
		yp(3,1)=q**3*y(1)
		yp(4,1)=q**4*y(1)
	  case(5)
		yp(1,1)=q*y(1)
		yp(2,1)=q**2*y(1)
		yp(3,1)=q**3*y(1)
		yp(4,1)=q**4*y(1)
		yp(5,1)=q**5*y(1)
	  case(6)

		yp(1,1)=q*y(1)
		yp(2,1)=q**2*y(1)
		yp(3,1)=q**3*y(1)
		yp(4,1)=q**4*y(1)
		yp(5,1)=q**5*y(1)
		yp(6,1)=q**6*y(1)
		  case(7)

		yp(1,1)=q*y(1)
		yp(2,1)=q**2*y(1)
		yp(3,1)=q**3*y(1)
		yp(4,1)=q**4*y(1)
		yp(5,1)=q**5*y(1)
		yp(6,1)=q**6*y(1)
		yp(7,1)=q**7*y(1)
		
	
	  endselect
	  end



