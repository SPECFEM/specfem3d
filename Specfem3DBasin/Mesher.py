#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.units.angle import deg
from pyre.units.length import km


class Mesher(Component):
    
    class Inventory(Component.Inventory):

        from pyre.inventory import bool, dimensional, float, greaterEqual, int
    
        SAVE_MESH_FILES               = bool("save-files")
        SUPPRESS_UTM_PROJECTION       = bool("suppress-utm-projection")
        
        depth_block                   = dimensional('depth-block', default=0.0*km)

        LATITUDE_MAX                  = float("latitude-max")
        LATITUDE_MIN                  = float("latitude-min")
        LONGITUDE_MAX                 = float("longitude-max")
        LONGITUDE_MIN                 = float("longitude-min")
        
        NEX_ETA                       = int("nex-eta", default=64)
        NEX_XI                        = int("nex-xi", default=64)
        NPROC_ETA                     = int("nproc-eta", validator=greaterEqual(1), default=1)
        NPROC_XI                      = int("nproc-xi", validator=greaterEqual(1), default=1)
        UTM_PROJECTION_ZONE           = int("utm-projection-zone")

    def __init__(self, name):
        Component.__init__(self, name, "mesher")

    def _init(self):
        Component._init(self)

        # convert to kilometers
        self.DEPTH_BLOCK_KM = self.inventory.depth_block / km
        
        return

    def nproc(self):
        """Return the total number of processors needed."""
        return (self.inventory.NPROC_XI *
                self.inventory.NPROC_ETA)


# end of file
