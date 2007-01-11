#!/usr/bin/env python


from pyre.components import Component
from pyre.units.time import second


class Solver(Component):

    
    # name by which the user refers to this component
    name = "solver"
    
    
    #
    #--- parameters
    #
    
    import pyre.inventory as pyre
    #from CMTSolution import cmtValidator

    ABSORBING_CONDITIONS          = pyre.bool("absorbing-conditions")
    CREATE_SHAKEMAP               = pyre.bool("create-shakemap")
    MOVIE_SURFACE                 = pyre.bool("movie-surface")
    MOVIE_VOLUME                  = pyre.bool("movie-volume")
    PRINT_SOURCE_TIME_FUNCTION    = pyre.bool("print-source-time-function")
    SAVE_DISPLACEMENT             = pyre.bool("save-displacement")
    SAVE_FORWARD                  = pyre.bool("save-forward")
    USE_HIGHRES_FOR_MOVIES        = pyre.bool("use-highres-for-movies")
    dry                           = pyre.bool("dry")

    record_length                 = pyre.dimensional("record-length", default=0.0*second)

    HDUR_MOVIE                    = pyre.float("hdur-movie")

    cmtSolution                   = pyre.inputFile("cmt-solution", default="DATA/CMTSOLUTION")
    stations                      = pyre.inputFile("stations", default="DATA/STATIONS")
    STATIONS_FILTERED             = pyre.str("stations-filtered", default="DATA/STATIONS_FILTERED")

    NTSTEP_BETWEEN_FRAMES         = pyre.int("ntstep-between-frames")
    NTSTEP_BETWEEN_OUTPUT_INFO    = pyre.int("ntstep-between-output-info")
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = pyre.int("ntstep-between-output-seismos")

    simulation_type               = pyre.str("simulation-type", validator=pyre.choice(['forward', 'adjoint', 'both']), default='forward')

    
    #
    #--- configuration
    #
    
    def _configure(self):
        Component._configure(self)
        
        # convert to seconds
        self.RECORD_LENGTH_IN_SECONDS = self.record_length / second

        st = { 'forward':1, 'adjoint':2, 'both':3 }
        self.SIMULATION_TYPE = st[self.simulation_type]

        return


    #
    #--- final initialization
    #
    
    def _init(self):
        Component._init(self)

        from os.path import abspath, join
        self.CMTSOLUTION = abspath(self.cmtSolution.name)
        self.STATIONS = abspath(self.stations.name)


    def setOutputDirectories(self, LOCAL_PATH, OUTPUT_FILES):
        from os.path import abspath, join
        
        self.LOCAL_PATH = LOCAL_PATH
        self.OUTPUT_FILES = OUTPUT_FILES

        # always written by the mesher
        self.HEADER_FILE = abspath(join(OUTPUT_FILES, 'values_from_mesher.h'))


    #
    #--- building
    #
    
    def build(self, script, srcdir):
        """Build the solver."""

        import os, os.path, sys
        from os.path import abspath, join

        outputDir = abspath(script.outputDir)
        pyspecfem3D = abspath(script.mpiExecutable)

        wd = os.getcwd()
        print "cd", srcdir
        os.chdir(srcdir)
        
        # create the include file for the solver
        self.createheader(script)
        
        # now finally build the solver
        argv = ['make', 'OUTPUT_DIR=' + outputDir, pyspecfem3D]
        print ' '.join(argv)
        status = os.spawnvp(os.P_WAIT, argv[0], argv)
        if status != 0:
            sys.exit("%s: %s: exit %d" % (sys.argv[0], argv[0], status))

        print "cd", wd
        os.chdir(wd)

        return
    
    
    def createheader(self, script):
        """Create the include file for the solver."""

        import os, shutil, sys
        from os.path import exists
        from PyxParameters import create_header_file
        
        # This path is hardwired into the Fortran source.
        oldHeader = 'OUTPUT_FILES/values_from_mesher.h'

        # If the header doesn't exist, simply create it.
        if not exists(oldHeader):
            self.HEADER_FILE = oldHeader
            create_header_file(script) # call into Fortran
            return
        
        # First generate the header into a temporary file.
        from tempfile import mktemp
        newHeader  = mktemp()
        self.HEADER_FILE = newHeader
        create_header_file(script) # call into Fortran
        
        # Did the header file change?
        argv = ['diff', oldHeader, newHeader]
        print ' '.join(argv)
        status = os.spawnvp(os.P_WAIT, argv[0], argv)
        if status == 0:
            # Nope!  Nothing to do here.
            os.remove(newHeader)
            return
        if status != 1:
            # diff countered a problem
            os.remove(newHeader)
            sys.exit("%s: %s: exit %d" % (sys.argv[0], argv[0], status))

        # Replace the old header with the new one.
        print "mv", newHeader, oldHeader
        shutil.move(newHeader, oldHeader)

        return
    
    
    #
    #--- execution
    #
    
    def execute(self, script):
        """Execute the solver."""
        from PyxSpecfem import specfem3D
        if self.dry:
            print "execute", specfem3D
        else:
            specfem3D(script) # call into Fortran
        return


    #
    #--- clean-up
    #
    
    def collectSeismograms(self):
        """collect seismograms"""

        if self.dry:
            return
        
        import os, shutil, tarfile
        from os.path import basename, join
        from glob import glob
        from cig.seismo.sac import asc2sac
        from mpi import MPI_Comm_rank, MPI_COMM_WORLD
        
        rank = MPI_Comm_rank(MPI_COMM_WORLD)
        scratchSeismogramArchive = join(self.LOCAL_PATH, "seismograms-%d.tar.gz" % rank)
        archive = self.scratchSeismogramArchive
        
        files = []
        for sem in glob(join(self.LOCAL_PATH, "*.sem")):
            files.append(sem)
            sac = sem + ".sac"
            asc2sac(sem, sac)
            files.append(sac)
        for semd in glob(join(self.LOCAL_PATH, "*.semd")):
            files.append(semd)
            sac = semd + ".sac"
            asc2sac(semd, sac)
            files.append(sac)
        if len(files) == 0:
            # no seismograms on this node
            archive.close()
            os.remove(archive.name)
            return
        
        # A compressed tar file is between 10-15% of the size of the
        # raw data files.  By taring and compressing it now -- in
        # parallel, on local filesystems -- we sharply reduce the
        # amount of data we have to shovel over the network.
        
        tgz = tarfile.open(archive.name, "w:gz", archive)
        for name in files:
            arcname = basename(name)
            tgz.add(name, arcname)
        tgz.close()
        archive.close()

        # Copy the archive to the shared filesystem.

        src = archive.name
        dst = join(self.OUTPUT_FILES, basename(src))
        shutil.copyfile(src, dst)

        return


    def collectOutputFiles(self):
        """collect output files"""
        
        if self.dry:
            return
        
        import os, tarfile
        from os.path import basename, join

        seismogramArchive = join(self.OUTPUT_FILES, "seismograms.tar.gz")
        archiveOut = open(seismogramArchive, "w")
        skipList = ['pyspecfem3D', basename(archiveOut.name)]

        # Archive output files -- including the intermediate seismogram
        # archives delivered from the compute nodes.

        filesIn = []
        archivesIn = []
        for name in os.listdir(self.OUTPUT_FILES):
            if name in skipList:
                continue
            pathname = join(self.OUTPUT_FILES, name)
            if name.startswith("seismograms-"):
                archivesIn.append(pathname)
            else:
                filesIn.append((pathname, name))
        if len(filesIn) == 0:
            self._warning.log("No output files!")
            archiveOut.close()
            os.remove(archiveOut.name)
            return

        tgzOut = tarfile.open(archiveOut.name, "w:gz", archiveOut)
        
        # Rearchive seismograms.

        for archiveIn in archivesIn:
            tgzIn = tarfile.open(archiveIn, "r:gz")
            for member in tgzIn.getmembers():
                seismogram = tgzIn.extractfile(member)
                tgzOut.addfile(member, seismogram)
            tgzIn.close()

        # Archive other output files.
        
        for name, arcname in filesIn:
            tgzOut.add(name, arcname)
        
        tgzOut.close()
        archiveOut.close()

        # Delete the intermediate seismogram archives.
        
        for archiveIn in archivesIn:
            os.remove(archiveIn)

        return


# end of file
