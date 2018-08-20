#!/bin/bash 
python3 coeffient.py
for ((i=1;i<20;i++));do
	echo $i
	sed  's/Derivative(y(x), x)/y1/g' a.txt > b.txt
	sed  's/Derivative(z(x), x)/z1/g' b.txt > a.txt
	sed  's/Derivative(y(x), (x, '$i'))/y'$i'/g' a.txt > b.txt
	sed  's/Derivative(z(x), (x, '$i'))/z'$i'/g' b.txt > a.txt
done
sed  's/y(x)/y/g' a.txt > b.txt
sed  's/z(x)/z/g' b.txt > a.txt
#
