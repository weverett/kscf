all: kscf

kscf: kscf_init.o kscf_energx.o 
	cp *.o ../object/
	mv *.o objects

kscf_init.o: kscf_init.c
	icc -std=c99 -I../ -L../simint/libsimint.a -c kscf_init.c 

kscf_energx.o: kscf_energx.c
	icc -std=c99 -I../ -L../simint/libsimint.a -c kscf_energx.c 

clean:
	rm objects/*.o
