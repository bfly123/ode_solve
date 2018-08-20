#This is a script from f_derivative.sh
from sympy import *
#from scipy import *
x =symbols("x",real=True)
yp11 =symbols("yp11",cls=Function)
yp12 =symbols("yp12",cls=Function)
yp13 =symbols("yp13",cls=Function)
y1 =symbols("y1",cls=Function)
y2 =symbols("y2",cls=Function)
y3 =symbols("y3",cls=Function)
f=open(r'b.dat','w')
for i in range(1,7+1):
     z1=diff( -4.0e-2*y1(x)+1.0e4*y2(x)*y3(x) ,x,i)
     s='yp('+str(i+1)+',1)='+str(z1)  
     print(s,file=f)  
     z2=diff( 4.0e-2*y1(x)-1.0e4*y2(x)*y3(x)-3.0e7*y2(x)**2 ,x,i)
     s='yp('+str(i+1)+',2)='+str(z2)  
     print(s,file=f)  
     z3=diff( 3.0e7*y2(x)**2 ,x,i)
     s='yp('+str(i+1)+',3)='+str(z3)  
     print(s,file=f)  
