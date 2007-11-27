#!/usr/bin/env python


from pyre.components import Component
from pyre.units.length import km


class Mesher(Component):


    # name by which the user refers to this component
    name = "mesher"


    #
    #--- parameters
    #

    import pyre.inventory as pyre

    SAVE_MESH_FILES               = pyre.bool("save-files")
    SUPPRESS_UTM_PROJECTION       = pyre.bool("suppress-utm-projection")
    dry                           = pyre.bool("dry")

    depth_block                   = pyre.dimensional('depth-block', default=0.0*km)

    LATITUDE_MAX                  = pyre.float("latitude-max")
    LATITUDE_MIN                  = pyre.float("latitude-min")
    LONGITUDE_MAX                 = pyre.float("longitude-max")
    LONGITUDE_MIN                 = pyre.float("longitude-min")

    NEX_ETA                       = pyre.int("nex-eta", default=64)
    NEX_XI                        = pyre.int("nex-xi", default=64)
    NPROC_ETA                     = pyre.int("nproc-eta", validator=pyre.greaterEqual(1), default=1)
    NPROC_XI                      = pyre.int("nproc-xi", validator=pyre.greaterEqual(1), default=1)
    UTM_PROJECTION_ZONE           = pyre.int("utm-projection-zone")


    #
    #--- configuration
    #
    
    def _configure(self):
        Component._configure(self)

        # convert to kilometers
        self.DEPTH_BLOCK_KM = self.depth_block / km
        
        return


    def nproc(self):
        """Return the total number of processors needed."""
        return self.NPROC_XI * self.NPROC_ETA


    #
    #--- execution
    #
    
    def execute(self, script):
        """Execute the mesher."""
        from PyxMeshfem import meshfem3D
        if self.dry:
            print "execute", meshfem3D
        else:
            meshfem3D(script) # call into Fortran
        return


# end of file
