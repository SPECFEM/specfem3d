# Makefile

#############################################################
## modify to match your compiler defaults 

# compilers
F90 = gfortran
CC = gcc

#############################################################


LIB_OBJS = create_adjsrc_traveltime.o rw_ascfile_c.o 

# targets
all: xcreate_adjsrc_traveltime

xcreate_adjsrc_traveltime: $(LIB_OBJS)
	${F90} -Wall -o xcreate_adjsrc_traveltime $(LIB_OBJS)

create_adjsrc_traveltime.o: create_adjsrc_traveltime.f90
	${F90} -Wall -c create_adjsrc_traveltime.f90

rw_ascfile_c.o: ../../lib/rw_ascfile_c.c
	${CC} -c -o rw_ascfile_c.o ../../lib/rw_ascfile_c.c


clean:
	rm -f *.o xcreate_adjsrc_traveltime
