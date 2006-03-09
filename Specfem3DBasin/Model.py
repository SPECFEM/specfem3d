#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.inventory.Facility import Facility


class Model(Component):

    class Inventory(Component.Inventory):

        from pyre.inventory import bool, inputFile
        
        ATTENUATION                   = bool("attenuation")
        OCEANS                        = bool("oceans")
        TOPOGRAPHY                    = bool("topography")
        USE_OLSEN_ATTENUATION         = bool("use-olsen-attenuation")

        topoBathyFile                 = inputFile("topo-bathy-file", default="DATA/la_topography/topo_bathy_final.dat")
        
    def __init__(self, name):
        Component.__init__(self, name, "model")
        self.aliases.extend(self.classAliases)
        self.PATHNAME_TOPO_FILE = None

    def _init(self):
        Component._init(self)
        
        # Access our InputFile inventory items to make sure they're
        # readable.  (They will be reopened by the Fortran code.)
        if self.inventory.TOPOGRAPHY or self.inventory.OCEANS:
            f = self.inventory.topoBathyFile
            self.PATHNAME_TOPO_FILE = f.name
            f.close()


# built-in models

def BuiltInModel(name, *aliases):
    return type(Model)(
        'BuiltInModel(%s)' % name,
        (Model,),
        {'className': name,
        'classAliases': [name] + list(aliases)}
        )


class ModelHarvardLA(Model):
    className = "Harvard_LA"
    classAliases = ["Harvard_LA"]
    class Inventory(Model.Inventory):
        from pyre.inventory import inputFile
        basin_3d_high_res            = inputFile("basin-3d-high-res",
                                                 default="DATA/la_3D_block_harvard/la_3D_high_res/LA_HR_voxet_extracted.txt")
        basin_3d_medium_res          = inputFile("basin-3d-medium-res",
                                                 default="DATA/la_3D_block_harvard/la_3D_medium_res/LA_MR_voxet_extracted.txt")
        hauksson_regional            = inputFile("hauksson-regional",
                                                 default="DATA/hauksson_model/hauksson_final_grid_smooth.dat")
        moho_map                     = inputFile("moho-map",
                                                 default="DATA/moho_map/moho_lupei_zhu.dat")
        basement_map                 = inputFile("basement-map",
                                                 default="DATA/la_basement/reggridbase2_filtered_ascii.dat")
        salton_sea                   = inputFile("salton-sea",
                                                 default="DATA/st_3D_block_harvard/regrid3_vel_p.bin")
    def __init__(self, name):
        Model.__init__(self, name)
        self.BASIN_MODEL_3D_HIGH_RES_FILE    = None
        self.BASIN_MODEL_3D_MEDIUM_RES_FILE  = None
        self.HAUKSSON_REGIONAL_MODEL_FILE    = None
        self.MOHO_MAP_FILE                   = None
        self.BASEMENT_MAP_FILE               = None
        self.SALTON_SEA_MODEL_FILE           = None
    def _init(self):
        Model._init(self)
        # Access our InputFile inventory items to make sure they're
        # readable.  (They will be reopened by the Fortran code.)
        f = self.inventory.basin_3d_high_res;    self.BASIN_MODEL_3D_HIGH_RES_FILE    = f.name;  f.close()
        f = self.inventory.basin_3d_medium_res;  self.BASIN_MODEL_3D_MEDIUM_RES_FILE  = f.name;  f.close()
        f = self.inventory.hauksson_regional;    self.HAUKSSON_REGIONAL_MODEL_FILE    = f.name;  f.close()
        f = self.inventory.moho_map;             self.MOHO_MAP_FILE                   = f.name;  f.close()
        f = self.inventory.basement_map;         self.BASEMENT_MAP_FILE               = f.name;  f.close()
        f = self.inventory.salton_sea;           self.SALTON_SEA_MODEL_FILE           = f.name;  f.close()


builtInModelClasses = [
    BuiltInModel("SoCal"),
    BuiltInModel("Lacq_gas_field_France"),
    BuiltInModel("Min_Chen_anisotropy"),
    ModelHarvardLA,
    ]

def retrieveBuiltInModelClass(componentName):
    for cls in builtInModelClasses:
        for name in cls.classAliases:
            if name == componentName:
                return cls
    return None


# model facility

class ModelFacility(Facility):

    def __init__(self, name):
        Facility.__init__(self, name, default="SoCal")
        return
    
    def _retrieveComponent(self, instance, componentName):
        cls = retrieveBuiltInModelClass(componentName)
        if cls is None:
            return Facility._retrieveComponent(self, instance, componentName)
        model = cls(componentName)
        model.aliases.append(self.name)
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('built-in')
        return model, locator


# end of file
