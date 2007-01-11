# Process this file with Pyrex to produce PyxSpecfem.c


"""Python bindings for the SPECFEM3D Global Solver."""


# include 'config.h' in order to get the definitions of FC_FUNC and FC_FUNC_
cdef extern from "config.h":
    pass


# external Fortran functions

cdef extern void specfem3D_f "FC_FUNC(specfem3d, SPECFEM3D)" () except *
def specfem3D(arg):
    """Run the SPECFEM3D Global Solver."""
    import PyxParameters
    PyxParameters.component = arg
    specfem3D_f()


# end of file
