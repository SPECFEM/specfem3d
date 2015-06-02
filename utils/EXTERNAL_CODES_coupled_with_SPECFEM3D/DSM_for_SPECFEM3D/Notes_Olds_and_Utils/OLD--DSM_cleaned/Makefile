
# Fortran90 + MPI and C + MPI compilers
F90 = mpif90
CC  = mpicc

# obj directory
OBJDIR = obj

# Intel ifort with full checking options to debug (slow but very useful to check everything)
# option "-assume buffered_io" is important especially on
# parallel file systems like SFS 3. / Lustre 1.8. If omitted
# I/O throughput lingers at .5 MB/s, with it it can increase to ~44 MB/s
# However it does not make much of a difference on NFS mounted volumes or with SFS 3.1.1 / Lustre 1.6.7.1
# Intel ifort optimized for speed for production runs (add -mcmodel=medium -shared-intel to use more than 2 GB of memory)
FFLAGS = -O3 -check nobounds -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -diag-disable 7165 -gen-interfaces -warn all -module $(OBJDIR) # -DUSE_SERIAL
        # useful for debugging...
        # for debugging (the code will be slow but will check many things), change to this:
#FFLAGS = -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -diag-disable 7165 -gen-interfaces -warn all -module $(OBJDIR) # -DUSE_SERIAL

# GNU gfortran
# add -frange-check for range checking
#FFLAGS = -std=f2003 -fimplicit-none -O3 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -J $(OBJDIR) # -DUSE_SERIAL

# Cray Fortran with full checking
#FFLAGS = -eC -eD -ec -en -eI -ea -g -G0

# Portland pgf90 flags
#FFLAGS = -O1 -Mdclchk -Mdepchk -Minfo -Mstandard # -DUSE_SERIAL

CFLAGS = -O2

default: | $(OBJDIR) xdsmti 

all: clean default

# EXECUTABLES:

LIBS_FULL =

$(OBJDIR):
	@test -d $(OBJDIR) || (rm -f $(OBJDIR); mkdir -p $(OBJDIR))

xdsmti: $(OBJDIR)/compwave.o \
              $(OBJDIR)/dsmti.o \
              $(OBJDIR)/matrix.o \
              $(OBJDIR)/msub.o \
              $(OBJDIR)/define_Gauss_points_and_weights.o \
              $(OBJDIR)/param.o \
              $(OBJDIR)/sacbw.o \
              $(OBJDIR)/sacconv.o
	 ${F90} $(FFLAGS) -o xdsmti \
	                  $(OBJDIR)/compwave.o \
	                  $(OBJDIR)/dsmti.o \
	                  $(OBJDIR)/matrix.o \
	                  $(OBJDIR)/msub.o \
	                  $(OBJDIR)/define_Gauss_points_and_weights.o \
	                  $(OBJDIR)/param.o \
	                  $(OBJDIR)/sacbw.o \
	                  $(OBJDIR)/sacconv.o \
	                  $(LIBS_FULL) 

clean:
	rm -rf *.o *.mod xdsmti $(OBJDIR) 

# targets to clean up all object directories
purge: clean
realclean: purge
deepclean: purge
mrproper:  purge

# OBJECTS:  

$(OBJDIR)/define_Gauss_points_and_weights.o:  define_Gauss_points_and_weights.f90 
	 ${F90} $(FFLAGS) -c -o $(OBJDIR)/define_Gauss_points_and_weights.o define_Gauss_points_and_weights.f90

$(OBJDIR)/param.o:  param.f90 
	 ${F90} $(FFLAGS) -c -o $(OBJDIR)/param.o param.f90

$(OBJDIR)/compwave.o: compwave.f90 
	${F90} $(FFLAGS) -c -o $(OBJDIR)/compwave.o compwave.f90

$(OBJDIR)/dsmti.o: dsmti.F90 
	${F90} $(FFLAGS) -c -o $(OBJDIR)/dsmti.o dsmti.F90

$(OBJDIR)/sacconv.o: sacconv.f90 
	${F90} $(FFLAGS) -c -o $(OBJDIR)/sacconv.o sacconv.f90

$(OBJDIR)/matrix.o: matrix.f90
	${F90} $(FFLAGS) -c -o $(OBJDIR)/matrix.o matrix.f90

$(OBJDIR)/msub.o: msub.F90
	${F90} $(FFLAGS) -c -o $(OBJDIR)/msub.o msub.F90

############# the routine below is in C rather than Fortran90 #############
	
$(OBJDIR)/sacbw.o: sacbw.c 
	${CC} $(CFLAGS) -c -o $(OBJDIR)/sacbw.o sacbw.c

