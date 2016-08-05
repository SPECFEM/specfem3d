## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##                    Dummy program Makefile to test FTI library

##   PLEASE SET THIS VARIABLE BEFORE COMPILING
FTIPATH		= /path/to/fti/install/directory

##   COMPILERS
MPICC 		= mpicc
MPIRUN 		= mpirun

##   FLAGS
# Compiling using shared library
FTIFLAG 	= -I$(FTIPATH)/include -L$(FTIPATH)/lib -lfti
# Compiling using static library
#FTIFLAG 	= -I$(FTIPATH)/include $(FTIPATH)/lib/libfti.a

##   TARGETS
all: 		prog
prog:		prog.c
		$(MPICC) -o prog prog.c $(FTIFLAG)
# 1st arg: Memory/process in MB & 2nd arg: Number of iterations
test:
		$(MPIRUN) -np 8 ./prog 128 1500 config.fti
clean:
		rm -f prog
.PHONY:		prog test clean

