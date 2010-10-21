# Makefile

#############################################################
## modify to match your compiler defaults 

F90 = gfortran 
CC = gcc

#############################################################


F90FLAGS = -O2
CCFLAGS = -O2

CSRC = rw_ascfile_c
F90SRC = rw_ascfile_f

F90OBJ = $(patsubst %,%.o,$(F90SRC))
COBJ = $(patsubst %,%.o,$(CSRC))
OBJ = $(F90OBJ) $(COBJ)


all :  lib

lib : $(OBJ)


$(F90OBJ) : %.o : %.f90
	$(F90)  $(F90FLAGS) -c $*.f90

$(COBJ): %.o : %.c
	$(CC) $(CCFLAGS) -c $*.c 

clean: 
	\rm -f $(OBJS) *~ *.o 

cleanall:
	\rm -f $(OBJS)  *.o *~
