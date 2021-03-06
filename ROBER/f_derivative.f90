subroutine f_derivative(n,n_fun,y,yp)
	  integer n, n_fun
	  double precision yp(7,n_fun)
	  double precision y(n_fun)

	  select case(n)
	  case(1)
		! Begin f_y 3
		yp(1,1)=-0.04*y(1)+1.0e4*y(2)*y(3)
		yp(1,2)=0.04*y(1)-1.0e4*y(2)*y(3)-3.0e7*y(2)**2
		yp(1,3)=3.0e7*y(2)**2
		! End f_y
	  case(2)
		yp(1,1)=-0.04*y(1)+1.0e4*y(2)*y(3)
		yp(1,2)=0.04*y(1)-1.0e4*y(2)*y(3)-3.0e7*y(2)**2
		yp(1,3)=3.0e7*y(2)**2

		yp(2,1)=10000.0*y(2)*yp(1,3) + 10000.0*y(3)*yp(1,2) - 0.04*yp(1,1) 
		yp(2,2)=-60000000.0*y(2)*yp(1,2) - 10000.0*y(2)*yp(1,3) - 10000.0*y(3)*yp(1,2) + & 
		 0.04*yp(1,1) 
		yp(2,3)=60000000.0*y(2)*yp(1,2) 
		

	  case(3)
		yp(1,1)=-0.04*y(1)+1.0e4*y(2)*y(3)
		yp(1,2)=0.04*y(1)-1.0e4*y(2)*y(3)-3.0e7*y(2)**2
		yp(1,3)=3.0e7*y(2)**2

		yp(2,1)=10000.0*y(2)*yp(1,3) + 10000.0*y(3)*yp(1,2) - 0.04*yp(1,1) 
		yp(2,2)=-60000000.0*y(2)*yp(1,2) - 10000.0*y(2)*yp(1,3) - 10000.0*y(3)*yp(1,2) + & 
		 0.04*yp(1,1) 
		yp(2,3)=60000000.0*y(2)*yp(1,2) 
		
		yp(3,1)=10000.0*y(2)*yp(2,3) + 10000.0*y(3)*yp(2,2) - 0.04*yp(2,1) + 20000.0*yp( & 
		1,2)*yp(1,3) 
		yp(3,2)=-60000000.0*y(2)*yp(2,2) - 10000.0*y(2)*yp(2,3) - 10000.0*y(3)*yp(2,2) + & 
		 0.04*yp(2,1) - 60000000.0*yp(1,2)**2 - 20000.0*yp(1,2)*yp(1,3) 
		yp(3,3)=60000000.0*(y(2)*yp(2,2) + yp(1,2)**2) 

	  case(4)

		yp(1,1)=-0.04*y(1)+1.0e4*y(2)*y(3)
		yp(1,2)=0.04*y(1)-1.0e4*y(2)*y(3)-3.0e7*y(2)**2
		yp(1,3)=3.0e7*y(2)**2

		yp(2,1)=10000.0*y(2)*yp(1,3) + 10000.0*y(3)*yp(1,2) - 0.04*yp(1,1) 
		yp(2,2)=-60000000.0*y(2)*yp(1,2) - 10000.0*y(2)*yp(1,3) - 10000.0*y(3)*yp(1,2) + & 
		 0.04*yp(1,1) 
		yp(2,3)=60000000.0*y(2)*yp(1,2) 
		
		yp(3,1)=10000.0*y(2)*yp(2,3) + 10000.0*y(3)*yp(2,2) - 0.04*yp(2,1) + 20000.0*yp( & 
		1,2)*yp(1,3) 
		yp(3,2)=-60000000.0*y(2)*yp(2,2) - 10000.0*y(2)*yp(2,3) - 10000.0*y(3)*yp(2,2) + & 
		 0.04*yp(2,1) - 60000000.0*yp(1,2)**2 - 20000.0*yp(1,2)*yp(1,3) 
		yp(3,3)=60000000.0*(y(2)*yp(2,2) + yp(1,2)**2) 

		yp(4,1)=10000.0*y(2)*yp(3,3) + 10000.0*y(3)*yp(3,2) - 0.04*yp(3,1) + 30000.0*yp( & 
		1,2)*yp(2,3) + 30000.0*yp(2,2)*yp(1,3) 
		yp(4,2)=-60000000.0*y(2)*yp(3,2) - 10000.0*y(2)*yp(3,3) - 10000.0*y(3)*yp(3,2) + & 
		 0.04*yp(3,1) - 180000000.0*yp(1,2)*yp(2,2) - 30000.0*yp(1,2)*yp(2,3) - 30000.0* & 
		yp(2,2)*yp(1,3) 
		yp(4,3)=60000000.0*y(2)*yp(3,2) + 180000000.0*yp(1,2)*yp(2,2) 
		
	  case(5)
		yp(1,1)=-0.04*y(1)+1.0e4*y(2)*y(3)
		yp(1,2)=0.04*y(1)-1.0e4*y(2)*y(3)-3.0e7*y(2)**2
		yp(1,3)=3.0e7*y(2)**2

		yp(2,1)=10000.0*y(2)*yp(1,3) + 10000.0*y(3)*yp(1,2) - 0.04*yp(1,1) 
		yp(2,2)=-60000000.0*y(2)*yp(1,2) - 10000.0*y(2)*yp(1,3) - 10000.0*y(3)*yp(1,2) + & 
		 0.04*yp(1,1) 
		yp(2,3)=60000000.0*y(2)*yp(1,2) 
		
		yp(3,1)=10000.0*y(2)*yp(2,3) + 10000.0*y(3)*yp(2,2) - 0.04*yp(2,1) + 20000.0*yp( & 
		1,2)*yp(1,3) 
		yp(3,2)=-60000000.0*y(2)*yp(2,2) - 10000.0*y(2)*yp(2,3) - 10000.0*y(3)*yp(2,2) + & 
		 0.04*yp(2,1) - 60000000.0*yp(1,2)**2 - 20000.0*yp(1,2)*yp(1,3) 
		yp(3,3)=60000000.0*(y(2)*yp(2,2) + yp(1,2)**2) 

		yp(4,1)=10000.0*y(2)*yp(3,3) + 10000.0*y(3)*yp(3,2) - 0.04*yp(3,1) + 30000.0*yp( & 
		1,2)*yp(2,3) + 30000.0*yp(2,2)*yp(1,3) 
		yp(4,2)=-60000000.0*y(2)*yp(3,2) - 10000.0*y(2)*yp(3,3) - 10000.0*y(3)*yp(3,2) + & 
		 0.04*yp(3,1) - 180000000.0*yp(1,2)*yp(2,2) - 30000.0*yp(1,2)*yp(2,3) - 30000.0* & 
		yp(2,2)*yp(1,3) 
		yp(4,3)=60000000.0*y(2)*yp(3,2) + 180000000.0*yp(1,2)*yp(2,2) 

		yp(5,1)=10000.0*y(2)*yp(4,3) + 10000.0*y(3)*yp(4,2) - 0.04*yp(4,1) + 40000.0*yp( & 
		1,2)*yp(3,3) + 60000.0*yp(2,2)*yp(2,3) + 40000.0*yp(3,2)*yp(1,3) 
		yp(5,2)=-60000000.0*y(2)*yp(4,2) - 10000.0*y(2)*yp(4,3) - 10000.0*y(3)*yp(4,2) + & 
		 0.04*yp(4,1) - 240000000.0*yp(1,2)*yp(3,2) - 40000.0*yp(1,2)*yp(3,3) - 180000000.0&
		 *yp(2,2)**2 - 60000.0*yp(2,2)*yp(2,3) - 40000.0*yp(3,2)*yp(1,3) 
		yp(5,3)=60000000.0*y(2)*yp(4,2) + 240000000.0*yp(1,2)*yp(3,2) + 180000000.0*yp(2 & 
		,2)**2 
	  case(6)
		yp(1,1)=-0.04*y(1)+1.0e4*y(2)*y(3)
		yp(1,2)=0.04*y(1)-1.0e4*y(2)*y(3)-3.0e7*y(2)**2
		yp(1,3)=3.0e7*y(2)**2

		yp(2,1)=10000.0*y(2)*yp(1,3) + 10000.0*y(3)*yp(1,2) - 0.04*yp(1,1) 
		yp(2,2)=-60000000.0*y(2)*yp(1,2) - 10000.0*y(2)*yp(1,3) - 10000.0*y(3)*yp(1,2) + & 
		 0.04*yp(1,1) 
		yp(2,3)=60000000.0*y(2)*yp(1,2) 
		
		yp(3,1)=10000.0*y(2)*yp(2,3) + 10000.0*y(3)*yp(2,2) - 0.04*yp(2,1) + 20000.0*yp( & 
		1,2)*yp(1,3) 
		yp(3,2)=-60000000.0*y(2)*yp(2,2) - 10000.0*y(2)*yp(2,3) - 10000.0*y(3)*yp(2,2) + & 
		 0.04*yp(2,1) - 60000000.0*yp(1,2)**2 - 20000.0*yp(1,2)*yp(1,3) 
		yp(3,3)=60000000.0*(y(2)*yp(2,2) + yp(1,2)**2) 

		yp(4,1)=10000.0*y(2)*yp(3,3) + 10000.0*y(3)*yp(3,2) - 0.04*yp(3,1) + 30000.0*yp( & 
		1,2)*yp(2,3) + 30000.0*yp(2,2)*yp(1,3) 
		yp(4,2)=-60000000.0*y(2)*yp(3,2) - 10000.0*y(2)*yp(3,3) - 10000.0*y(3)*yp(3,2) + & 
		 0.04*yp(3,1) - 180000000.0*yp(1,2)*yp(2,2) - 30000.0*yp(1,2)*yp(2,3) - 30000.0* & 
		yp(2,2)*yp(1,3) 
		yp(4,3)=60000000.0*y(2)*yp(3,2) + 180000000.0*yp(1,2)*yp(2,2) 

		yp(5,1)=10000.0*y(2)*yp(4,3) + 10000.0*y(3)*yp(4,2) - 0.04*yp(4,1) + 40000.0*yp( & 
		1,2)*yp(3,3) + 60000.0*yp(2,2)*yp(2,3) + 40000.0*yp(3,2)*yp(1,3) 
	
		yp(5,2)=-60000000.0*y(2)*yp(4,2) - 10000.0*y(2)*yp(4,3) - 10000.0*y(3)*yp(4,2) + & 
		 0.04*yp(4,1) - 240000000.0*yp(1,2)*yp(3,2) - 40000.0*yp(1,2)*yp(3,3) - 180000000.0&
		 *yp(2,2)**2 - 60000.0*yp(2,2)*yp(2,3) - 40000.0*yp(3,2)*yp(1,3)
		yp(5,3)=60000000.0*y(2)*yp(4,2) + 240000000.0*yp(1,2)*yp(3,2) + 180000000.0*yp(2 & 
		,2)**2 

		yp(6,1)=10000.0*y(2)*yp(5,3) + 10000.0*y(3)*yp(5,2) - 0.04*yp(5,1) + 50000.0*yp( & 
		1,2)*yp(4,3) + 100000.0*yp(2,2)*yp(3,3) + 100000.0*yp(3,2)*yp(2,3) + 50000.0*yp( & 
		4,2)*yp(1,3) 
		yp(6,2)=-60000000.0*y(2)*yp(5,2) - 10000.0*y(2)*yp(5,3) - 10000.0*y(3)*yp(5,2) + & 
		 0.04*yp(5,1) - 300000000.0*yp(1,2)*yp(4,2) - 50000.0*yp(1,2)*yp(4,3) - 600000000.0&
		 *yp(2,2)*yp(3,2) - 100000.0*yp(2,2)*yp(3,3) - 100000.0*yp(3,2)*yp(2,3) - 50000.0*&
		 yp(4,2)*yp(1,3) 
		yp(6,3)=60000000.0*y(2)*yp(5,2) + 300000000.0*yp(1,2)*yp(4,2) + 600000000.0*yp(2 & 
		,2)*yp(3,2) 
	  endselect
	  end



