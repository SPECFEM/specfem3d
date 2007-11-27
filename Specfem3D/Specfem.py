#!/usr/bin/env python


try:
    import _mpi
except ImportError:
    from pyre.applications import Script as Base
else:
    from mpi import Application as Base


class Specfem(Base):
    
    
    # name by which the user refers to this application
    name = "Specfem3DBasin"
    
    
    #
    #--- inventory
    #
    
    import pyre.inventory as pyre
    
    from Mesher import Mesher
    from Solver import Solver

    outputDir                     = pyre.str("output-dir", default="OUTPUT_FILES")
    scratchDir                    = pyre.str("scratch-dir", default="/scratch")
        
    model                         = pyre.facility("model", default="Harvard_LA")
    mesher                        = pyre.facility("mesher", factory=Mesher)
    solver                        = pyre.facility("solver", factory=Solver)
    
    command = pyre.str("command", validator=pyre.choice(['mesh', 'solve', 'go']), default="go")


    #
    #--- configuration
    #
    
    def _configure(self):
        super(Specfem, self)._configure()
        
        from os import makedirs
        from os.path import isdir
        if not isdir(self.outputDir):
            makedirs(self.outputDir)
        
        # declare the interpreter to be used on the compute nodes
        from os.path import join
        self.mpiExecutable = join(self.outputDir, "pyspecfem3D") # includes solver

        # compute the total number of processors needed
        self.nodes = self.mesher.nproc()
        
        return


    #
    #--- final initialization
    #
    
    def _init(self):
        super(Specfem, self)._init()
        self.OUTPUT_FILES = self.outputDir
        self.LOCAL_PATH = self.scratchDir
        self.solver.setOutputDirectories(
            LOCAL_PATH = self.LOCAL_PATH,
            OUTPUT_FILES = self.OUTPUT_FILES
            )


    #
    #--- .odb files
    #
    

    def _getPrivateDepositoryLocations(self):
        from os.path import dirname, isdir, join
        models = join(dirname(__file__), 'models')
        assert isdir(models)
        return [models]

    
    #
    #--- executed on the login node in response to the user's command
    #
    
    def onLoginNode(self, *args, **kwds):
        
        # build the solver
        import __main__
        from os.path import dirname
        srcdir = dirname(__main__.__file__)
        self.solver.build(self, srcdir)
        
        # schedule the job (bsub)
        super(Specfem, self).onLoginNode(*args, **kwds)
        
        return


    #
    #--- executed when the batch job is scheduled
    #

    def onLauncherNode(self, *args, **kwds):
        super(Specfem, self).onLauncherNode(*args, **kwds) # mpirun
        self.solver.collectOutputFiles()
        return

    
    #
    #--- executed in parallel on the compute nodes
    #
    
    def onComputeNodes(self, *args, **kwds):

        # execute mesher and/or solver
        
        if self.command == "mesh" or self.command == "go":
            self.mesher.execute(self)
        
        if self.command == "solve" or self.command == "go":
            self.solver.execute(self)
            self.solver.collectSeismograms()

        return


    #
    #--- executed by the serial version of the code
    #
    
    def main(self, *args, **kwds):
        self.onComputeNodes(self, *args, **kwds)


# entry points

def main(*args, **kwds):
    script = Specfem()
    script.run(*args, **kwds)


# end of file
