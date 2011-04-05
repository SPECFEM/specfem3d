# Makefile

#############################################################
## modify to match your compiler defaults 

# compilers
F90 = gfortran
CC = gcc

#############################################################


LIB_OBJS = create_adjsrc_amplitude.o rw_ascfile_c.o 

# targets
all: xcreate_adjsrc_amplitude

xcreate_adjsrc_amplitude: $(LIB_OBJS)
	${F90} -Wall -o xcreate_adjsrc_amplitude $(LIB_OBJS)

create_adjsrc_amplitude.o: create_adjsrc_amplitude.f90
	${F90} -Wall -c create_adjsrc_amplitude.f90

rw_ascfile_c.o: ../../lib/rw_ascfile_c.c
	${CC} -c -o rw_ascfile_c.o ../../lib/rw_ascfile_c.c


clean:
	rm -f *.o xcreate_adjsrc_amplitude
