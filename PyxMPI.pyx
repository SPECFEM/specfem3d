# Process this file with Pyrex to produce PyxMPI.c


cimport mpi


cdef extern from "stdlib.h":
    void *malloc(int)
    void free(void *)


class MPI_Error(EnvironmentError):
    pass


def MPI_Init(argv):
    cdef int error, cargc, i
    cdef char **cargv, **mycargv
    myargv = []

    # Construct a C char** argument vector from 'argv'.
    cargc = len(argv)
    cargv = <char **>malloc((cargc + 1) * sizeof(char *))
    for i from 0 <= i < cargc:
        arg = argv[i]
        myargv.append(arg) # add ref
        cargv[i] = arg
    cargv[cargc] = NULL

    # Call MPI_Init().
    mycargv = cargv; # MPI might allocate & return its own.
    error = mpi.MPI_Init(&cargc, &cargv)
    if error != mpi.MPI_SUCCESS:
        free(mycargv)
        raise MPI_Error, error

    # Reconstruct Python's 'argv' from the modified 'cargv'.
    del argv[:]
    for i from 0 <= i < cargc:
        argv.append(cargv[i])
    free(mycargv)
    
    return


def MPI_Finalize():
    cdef int error
    error = mpi.MPI_Finalize()
    if error != mpi.MPI_SUCCESS:
        raise MPI_Error, error
    return


# end of file
