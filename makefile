Objects = global.o van_der_pol_per.o pertur.o pertur_G.o
source = global.f90 van_der_pol_per.f90 pertur.f90 pertur_G.f90
Bin = main.exe
F90 = gfortran -openmp

${Bin}:${Objects}
	${F90} -o ${Bin} ${Objects}

%.o :%.f90
	${F90} -c -o $@ $<


.PHONY:clean
clean:
	rm -f main.exe $(Objects) *.mod


