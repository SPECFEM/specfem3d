#!/usr/bin/env python


from CodecConfig import CodecConfig
import journal
from pyre.applications.Script import Script as PyreScript
import sys


# patch Pyre to our liking
import PyrePatches


def myexcepthook(type, value, traceback):
    sys.__excepthook__(type, value, traceback)
    import pdb
    pdb.post_mortem(traceback)


class Script(PyreScript):

    class Inventory(PyreScript.Inventory):

        from pyre.inventory import bool, facility, str
        from Launchers import LauncherFacility
        from Mesher import Mesher
        from Model import ModelFacility
        from Properties import OutputDir
        from Solver import Solver
        
        LOCAL_PATH                    = str("local-path")
        OUTPUT_FILES                  = OutputDir("output-dir")
        
        launcher                      = LauncherFacility("launcher")
        mesher                        = facility("mesher", factory=Mesher, args=["mesher"])
        model                         = ModelFacility("model")
        solver                        = facility("solver", factory=Solver, args=["solver"])
    
    def __init__(self):
        super(Script, self).__init__("Specfem3DBasin")
        self.error = journal.error(self.name)

    def _init(self):
        super(Script, self)._init()
        self.MODEL = self.inventory.model.className
        
        # compute the total number of processors needed
        nodes = self.inventory.mesher.nproc()
        self.inventory.launcher.inventory.nodes = nodes
        self.inventory.launcher.nodes = nodes
        
    
    def run(self, *args, **kwds):
        for i in xrange(0, len(sys.argv)):
            if sys.argv[i] == '-pdb':
                if sys.stdin.isatty():
                    sys.excepthook = myexcepthook
                del sys.argv[i]
                break
        try:
            super(Script, self).run(*args, **kwds)
        except PyrePatches.PropertyValueError, error:
            self.reportPropertyValueError(error)

    def raisePropertyValueError(self, traitName):
        desc = self.inventory.getTraitDescriptor(traitName)
        trait = self.inventory.getTrait(traitName)
        raise PyrePatches.PropertyValueError, (trait, desc.locator)

    def reportPropertyValueError(self, error):
        import linecache
        j = journal.journal()
        j.device.renderer = PyrePatches.PropertyValueError.Renderer()
        property = error.args[0]
        locator = error.args[1]
        entry = journal.diagnostics.Entry.Entry()
        meta = entry.meta
        meta["name"] = property.name
        meta["error"] = error
        if locator:
            meta["filename"] = locator.source
            try:
                line = locator.line
            except AttributeError:
                meta["src"] = ""
                meta["line"] = "???"
            else:
                src = linecache.getline(locator.source, locator.line)
                src = src.rstrip()
                meta["src"] = src
                meta["line"] = locator.line
        else:
            meta["filename"] = "???"
            meta["line"] = "???"
            meta["src"] = ""
        j.record(entry)
        
    def initializeCurator(self, curator, registry):
        cfg = CodecConfig()
        curator.registerCodecs(cfg)
        return super(Script, self).initializeCurator(curator, registry)
        
    def collectUserInput(self, registry):
        # read INI-style .cfg files
        curator = self.getCurator()
        configRegistry = curator.getTraits(self.name, extraDepositories=[], encoding='cfg')
        self.updateConfiguration(configRegistry)
        # read parameter files given on the command line
        from os.path import isfile, splitext
        for arg in self.argv:
            if isfile(arg):
                base, ext = splitext(arg)
                encoding = ext[1:] # NYI: not quite
                codec = self.getCurator().codecs.get(encoding)
                if codec:
                    shelf = codec.open(base)
                    paramRegistry = shelf['inventory'].getFacility(self.name)
                    if paramRegistry:
                        self.updateConfiguration(paramRegistry)
                else:
                    self.error.log("unknown encoding: %s" % ext)
            else:
                self.error.log("cannot open '%s'" % arg)
        return
        

# end of file
