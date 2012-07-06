
GPU_MODE Notes (by Max Rietmann)

src/shared/constants.h:

Change the following settings (for NOISE simulations):

! sources and receivers Z coordinates given directly instead of with depth
  logical, parameter :: USE_SOURCES_RECVS_Z = .true.

! the seismograms are normal to surface
! Z record corresponds to the normal, while E and N are two tangent vectors
! that completes an orthonormal.
  logical, parameter :: EXT_MESH_RECV_NORMAL = .true.

Settings for Eiger
./configure --with-cuda FC=mpif90 MPIFC=mpif90 CUDA_LIB=-L/apps/eiger/Cuda-4.0/cuda/lib64 MPI_INC=-I/apps/eiger/mvapich2/1.5.1p1/mvapich2-gnu/include
./configure --with-cuda FC=mpif90 MPIFC=mpif90 FLAGS_CHECK="-O3 -ffree-line-length-none -fopenmp" FLAGS_NO_CHECK="-Ofast -mfpmath=sse -funroll-loops -ffree-line-length-none -fopenmp" CUDA_LIB=-L/apps/eiger/Cuda-4.0/cuda/lib64 MPI_INC=-I/apps/eiger/mvapich2/1.5.1p1/mvapich2-gnu/include

./configure --with-cuda FC=mpif90 MPIFC=mpif90 FLAGS_CHECK="`echo $ICC_CHECK`" FLAGS_NO_CHECK='`echo $ICC_FCFLAGS`' CUDA_LIB=-L/apps/eiger/Cuda-4.0/cuda/lib64 MPI_INC=-I/apps/eiger/mvapich2/1.5.1p1/mvapich2-gnu/include

NOTE: Need to add below to get nvcc to use the older version of gcc, as it doesn't support gcc4.5 or above.
--compiler-bindir /usr/bin/gcc-4.3

 CUDA_LIB=-L$CUDA_HOME/lib64 MPI_INC=-I$MPICH_DIR/include


Settings for todi.cscs.ch
export GNU_FCFLAGS="-Ofast -mfpmath=sse -funroll-loops -ffree-line-length-none"
export GNU_CHECK="-O3 -ffree-line-length-none"
export DEBUG_FLAGS="-g -fbacktrace"
export ICC_FCFLAGS="-O3 -fp-model fast=2 -x SSE4.2 -ftz -funroll-loops -unroll5"
export ICC_FCFLAGS="-O3 -fp-model fast=2 -funroll-loops -unroll5 -msse3 -ftree-vectorize"
export ICC_CHECK="-O2"
export CRAY_FCFLAGS="-eF -em -rm -O3,fp3"
export CRAY_CHECK="-eF -em -rm"

# Cray configure
./configure --with-cuda CC=cc FC=ftn MPIFC=ftn MPICC=cc FLAGS_CHECK="`echo $CRAY_CHECK`" FLAGS_NO_CHECK='`echo $CRAY_FCFLAGS`' CUDA_LIB=-L$CUDA_HOME/lib64 MPI_INC=-I$MPICH_DIR/include
# GNU configure
./configure --with-cuda CC=cc FC=ftn MPIFC=ftn MPICC=cc FLAGS_CHECK="`echo $GNU_CHECK`" FLAGS_NO_CHECK='`echo $GNU_FCFLAGS`' CUDA_LIB=-L$CUDA_HOME/lib64 MPI_INC=-I$MPICH_DIR/include
# Intel convifugre
./configure --with-cuda CC=cc FC=ftn MPIFC=ftn MPICC=cc FLAGS_CHECK="`echo $ICC_CHECK`" FLAGS_NO_CHECK='`echo $ICC_FCFLAGS`' CUDA_LIB=-L$CUDA_HOME/lib64 MPI_INC=-I$MPICH_DIR/include

./configure CC=cc FC=ftn MPIFC=ftn MPICC=cc FLAGS_CHECK="`echo $CRAY_CHECK`" FLAGS_NO_CHECK='`echo $CRAY_FCFLAGS`' MPI_INC=-I$MPICH_DIR/include

suggested compiler options:
intel: -O3 -fp-model fast=2 -x SSE4.2 -ftz -funroll-loops -unroll5
gnu: -Ofast -mfpmath=sse -funroll-loops
pgi: 
####################################
Notes from Jeff Poznanovic:

Configure and build with:

module load cuda cudatools cudasdk
module load PrgEnv-cray
./configure --with-cuda --with-mpi MPIFC=ftn MPICC=cc FC=ftn CC=cc FCFLAGS="-eF -em -rm" CFLAGS="-h list=m" MPI_INC=-I/opt/cray/mpt/default/xt/gemini/mpich2-cray/73/include

 
make


