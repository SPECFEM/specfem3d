FC = ifort
MPIFC = mpif90
FLAGS_CHECK = -O3 -xHost -ftz -assume byterecl -funroll-loops -unroll5 -ftz -align sequence -vec-report0 -std95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds

####
#### targets
####

# default targets for the pure Fortran version
DEFAULT = \
	TraPSV_write_coef_mpi \
	TraPSV_read_coef_mpi_vert\
	TraPSV_MPI_read_Zmin\
	FFT_MPI_FACES_ZMIN\
	FFT_MPI_FACES_VERT\
	ChangeFormat_zmin_disp\
	ChangeFormat_zmin\
	ChangeFormat_disp\
	ChangeFormat

CLEAN   = \
	cleanTraPSV_read_coef_mpi_vert\
	cleanTraPSV_MPI_read_Zmin\
	cleanFFT_MPI_FACES_ZMIN\
	cleanFFT_MPI_FACES_VERT\
	cleanChangeFormat_zmin_disp\
	cleanChangeFormat_zmin\
	cleanChangeFormat_disp\
	cleanChangeFormat

default: $(DEFAULT)

all: default

clean : $(CLEAN)

TraPSV_write_coef_mpi :
	(cd TraPSV_MPI_write_computes_coefs_with_FEM; make)

TraPSV_read_coef_mpi_vert :
	(cd TraPSV_MPI_vertical_faces; make)

TraPSV_MPI_read_Zmin :
	(cd TraPSV_MPI_read_Zmin; make)

FFT_MPI_FACES_ZMIN :
	(cd FFT_MPI_FACES_ZMIN; make)

FFT_MPI_FACES_VERT :
	(cd FFT_MPI_FACES_VERT; make)

ChangeFormat_zmin_disp :
	(cd ChangeFormat_zmin_disp; make)

ChangeFormat_zmin :
	(cd ChangeFormat_zmin; make)

ChangeFormat_disp :
	(cd ChangeFormat_disp; make)

ChangeFormat :
	(cd ChangeFormat; make)

cleanTraPSV_read_coef_mpi_vert :
	(cd TraPSV_MPI_vertical_faces; make clean)

cleanTraPSV_MPI_read_Zmin :
	(cd TraPSV_MPI_read_Zmin; make clean)

cleanFFT_MPI_FACES_ZMIN :
	(cd FFT_MPI_FACES_ZMIN; make clean)

cleanFFT_MPI_FACES_VERT :
	(cd FFT_MPI_FACES_VERT; make clean)

cleanChangeFormat_zmin_disp :
	(cd ChangeFormat_zmin_disp; make clean)

cleanChangeFormat_zmin :
	(cd ChangeFormat_zmin; make clean)

cleanChangeFormat_disp :
	(cd ChangeFormat_disp; make clean)

cleanChangeFormat :
	(cd ChangeFormat; make clean)

