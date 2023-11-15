# Makefile

#############################################################
## modify to match your compiler defaults

F90 = gfortran
CC = gcc

# modify to point to your own SAC installation path
SAC_DIR = /opt/sac
SAC_INC = ${SAC_DIR}/include
SAC_LIB = ${SAC_DIR}/lib


#############################################################

F90FLAGS = -O2
CCFLAGS = -O2

CSRC = rw_ascfile_c
F90SRC = rw_ascfile_f

F90OBJ = $(patsubst %,%.o,$(F90SRC))
COBJ = $(patsubst %,%.o,$(CSRC))
OBJ = $(F90OBJ) $(COBJ)


all :  lib convolve_stf

lib : $(OBJ)


#
# convolve_stf
#
CFLAGS = -I$(SAC_INC)
LIBS = -L$(SAC_LIB) -lsac -lsacio

convolve_stf: convolve_stf.c
	$(CC) $(CFLAGS) -o convolve_stf convolve_stf.c $(LIBS)



#
# library routines
#
$(F90OBJ) : %.o : %.f90
	$(F90)  $(F90FLAGS) -c $*.f90

$(COBJ): %.o : %.c
	$(CC) $(CCFLAGS) -c $*.c

clean:
	\rm -f $(OBJS) *~ *.o convolve_stf

cleanall:
	\rm -f $(OBJS)  *.o *~ convolve_stf
