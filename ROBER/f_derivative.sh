#!/bin/bash
cat /dev/null > coefficient.py
echo "#This is a script from f_derivative.sh" >> coefficient.py
echo "from sympy import *" >> coefficient.py

ns=`awk '/Begin f_y/{print $4}' f_derivative.f90`
awk '/yp\(1\,/{print $0}' f_derivative.f90 > a.dat
for ((i=1;i<ns+2;i++));do
	sed -i ''  's/yp(1,'$i')/yp1'$i'/g;s/y('$i')/y'$i'\\(x\\)/g' a.dat
done
sed -i ''  's/=/= /' a.dat
echo "x =symbols(\"x\",real=True)" >> coefficient.py
for ((i=1;i<ns+1;i++));do
	echo "yp1$i =symbols(\"yp1$i\",cls=Function)" >> coefficient.py
done
for ((i=1;i<ns+1;i++));do
	echo "y$i =symbols(\"y$i\",cls=Function)" >> coefficient.py
done
echo "f=open(r'b.dat','w')" >> coefficient.py
echo "for i in range(1,7+1):">>coefficient.py

for ((i=1;i<ns+1;i++));do
	num=$i
eval	N$num=`awk 'NR=='$num' {
print $2
}' a.dat`
#eval echo "\$N$num"
done
for ((i=1;i<ns+1;i++));do
	eval echo  "\ \ \ \    z$i=diff\( \$N$i ,x,i\)"  >> coefficient.py
	echo  "     s='yp('+str(i+1)+',$i)='+str(z$i)  " >>coefficient.py
	echo  "     print(s,file=f)  " >>coefficient.py
done

python3 coefficient.py

for ((i=1;i<ns+2;i++));do
	sed  -i '' 's/Derivative(y'$i'(x), x)/yp(1,'$i')/g' b.dat
	for ((j=1;j<20;j++));do
		sed -i '' 's/Derivative(y'$i'(x), (x, '$j'))/yp('$j','$i')/g' b.dat
		#sed  's/Derivative(z(x), (x, '$i'))/z'$i'/g' b.txt > a.txt
	done
done
for ((i=1;i<ns+2;i++));do
	sed  -i '' 's/y'$i'(x)/y('$i')/g' b.dat
done

awk 'BEGIN{FS=""}{for (j=1;j<5;j++){
if(NF>j*80){
	for(i=(j-1)*80+1;i<=j*80;i++) printf $i
		printf " &"
		 print " "}
else {
	for(i=(j-1)*80+1;i<NF+1;i++) printf $i
	print " "}
}}' b.dat > a.dat
	#cat b.dat|cut -b -80
	sed -i '' '/^ $/d' a.dat
	#sed -i '' '/&/!G' a.dat
	sed  -i '' '/1)=/{x;p;x;}' a.dat
	sed  -i '' '1i \
	The derevatives are automaticly corrected by f_derivative.sh. Need yn from f_derivative.f90' a.dat
