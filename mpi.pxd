# Pyrex module for MPI

cdef extern from "mpi.h":
    enum MPI_Error:
        MPI_SUCCESS
    int MPI_Init(int *, char ***)
    int MPI_Finalize()

# end of file
