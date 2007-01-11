#!/usr/bin/env python


#
# base class for all models
#


from pyre.components import Component


class Model(Component):
    
    import pyre.inventory as pyre
    
    # parameters common to all models
    ATTENUATION                   = pyre.bool("attenuation")
    OCEANS                        = pyre.bool("oceans")
    TOPOGRAPHY                    = pyre.bool("topography")
    USE_OLSEN_ATTENUATION         = pyre.bool("use-olsen-attenuation")


# end of file
