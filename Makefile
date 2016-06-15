all: kscf

kscf: kscf_init.o kscf_energx.o 
	cp *.o ../object/
	mv *.o objects

kscf_init.o: kscf_init.c
	icc -I../ -w3 -diag-disable:remark -mkl -DMKL_ILP64 -c kscf_init.c 

kscf_energx.o: kscf_energx.c
	icc -I../ -w3 -diag-disable:remark -mkl -DMKL_ILP64 -c kscf_energx.c 

clean:
	rm objects/*.o
