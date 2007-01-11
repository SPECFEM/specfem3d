# Process this file with Pyrex to produce PyxMeshfem.c


"""Python bindings for the SPECFEM3D Global Solver."""


# include 'config.h' in order to get the definitions of FC_FUNC and FC_FUNC_
cdef extern from "config.h":
    pass


# external Fortran functions

cdef extern void meshfem3D_f "FC_FUNC(meshfem3d, MESHFEM3D)" () except *
def meshfem3D(arg):
    """Run the SPECFEM3D Global Mesher."""
    import PyxParameters
    PyxParameters.component = arg
    meshfem3D_f()


# end of file
