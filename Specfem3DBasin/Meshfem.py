#!/usr/bin/env python


from MPIDriver import MPIDriver
from Driver import Driver


try:
    import PyxMPI
except ImportError:
    Base = Driver
else:
    Base = MPIDriver


class Meshfem(Base):
    
    def __init__(self):
        super(Meshfem, self).__init__("output_mesher.txt")

    def main(self, *args, **kwds):
        #self.code.meshfem3D(self)
        self.code.xxxxfem3D(self)


# end of file
