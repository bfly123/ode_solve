#from __future__ import print_function
from sympy import *

x=symbols('x')
y=symbols('y',cls=Function)
z=symbols('z',cls=Function)

f=open(r'a.txt','w')
for a in range(1,7+1):
    zn=diff((1-y(x)**2)*z(x)-y(x),x,a)
    s='y'+str(a+1)+'=z'+str(a)
    #s=str(s).replace(' ','')
    print(s, file = f)
    s='z'+str(a+1)+'='+str(zn)+'/ss'
    #s=str(s).replace(' ','')
    print(s,file=f)
f.close()
