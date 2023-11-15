# Makefile

#############################################################
## modify to match your compiler defaults

# compilers
F90 = gfortran
CC = gcc

FLAGS = #-Wall

#############################################################


LIB_OBJS = create_adjsrc_waveform.o rw_ascfile_c.o

# targets
all: xcreate_adjsrc_waveform

xcreate_adjsrc_waveform: $(LIB_OBJS)
	${F90} $(FLAGS) -o xcreate_adjsrc_waveform $(LIB_OBJS)

create_adjsrc_waveform.o: create_adjsrc_waveform.f90
	${F90} $(FLAGS) -c create_adjsrc_waveform.f90

rw_ascfile_c.o: ../../scripts/lib/rw_ascfile_c.c
	${CC} -c -o rw_ascfile_c.o ../../scripts/lib/rw_ascfile_c.c


clean:
	rm -f *.o xcreate_adjsrc_waveform
