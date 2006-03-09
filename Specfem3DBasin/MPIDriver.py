#!/usr/bin/env python


from Driver import Driver
from mpi.Application import Application as MPIApplication


class MPIDriver(Driver, MPIApplication):

    """Base class for scripts which run SPECFEM3D Fortran code in
    parallel under MPI."""

    class Inventory(Driver.Inventory, MPIApplication.Inventory): pass

    def rank(self):
        # NYI: MPI_Comm_rank()
        return self.inventory.mode == "server" and 0 or 1
    

# end of file
