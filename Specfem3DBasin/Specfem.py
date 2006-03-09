#!/usr/bin/env python


from MPIDriver import MPIDriver
from Driver import Driver


try:
    import PyxMPI
except ImportError:
    Base = Driver
else:
    Base = MPIDriver


class Specfem(Base):
    
    def __init__(self):
        super(Specfem, self).__init__("output_solver.txt")

    def main(self, *args, **kwds):
        #self.code.specfem3D(self)
        self.code.xxxxfem3D(self)


# end of file
