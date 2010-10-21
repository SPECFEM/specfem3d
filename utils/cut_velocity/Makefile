# Makefile

#############################################################
## modify to match your compiler defaults 

# compilers
F90 = gfortran
CC = gcc

#############################################################


LIB_OBJS = cut_velocity.o rw_ascfile_c.o 

# targets
all: xcut_velocity

xcut_velocity: $(LIB_OBJS)
	${F90} -Wall -o xcut_velocity $(LIB_OBJS)

cut_velocity.o: cut_velocity.f90
	${F90} -Wall -c cut_velocity.f90

rw_ascfile_c.o: ../lib/rw_ascfile_c.c
	${CC} -c -o rw_ascfile_c.o ../lib/rw_ascfile_c.c


clean:
	rm -f *.o xcut_velocity
