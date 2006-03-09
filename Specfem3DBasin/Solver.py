#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.units.time import second


class Solver(Component):
    
    class Inventory(Component.Inventory):

        from pyre.inventory import bool, choice, dimensional, float, inputFile, int, str
    
        ABSORBING_CONDITIONS          = bool("absorbing-conditions")
        CREATE_SHAKEMAP               = bool("create-shakemap")
        MOVIE_SURFACE                 = bool("movie-surface")
        MOVIE_VOLUME                  = bool("movie-volume")
        PRINT_SOURCE_TIME_FUNCTION    = bool("print-source-time-function")
        SAVE_DISPLACEMENT             = bool("save-displacement")
        SAVE_FORWARD                  = bool("save-forward")
        USE_HIGHRES_FOR_MOVIES        = bool("use-highres-for-movies")
        
        record_length                 = dimensional("record-length", default=0.0*second)
        
        HDUR_MOVIE                    = float("hdur-movie")
        
        cmtSolution                   = inputFile("cmt-solution", default="DATA/CMTSOLUTION")
        stations                      = inputFile("stations", default="DATA/STATIONS")
        stations_filtered             = inputFile("stations-filtered", default="DATA/STATIONS_FILTERED")

        NTSTEP_BETWEEN_FRAMES         = int("ntstep-between-frames")
        NTSTEP_BETWEEN_OUTPUT_INFO    = int("ntstep-between-output-info")
        NTSTEP_BETWEEN_OUTPUT_SEISMOS = int("ntstep-between-output-seismos")

        simulation_type               = str("simulation-type", validator=choice(['forward', 'adjoint', 'both']), default='forward')
        
    def __init__(self, name):
        Component.__init__(self, name, "solver")
        self.CMTSOLUTION = None
        self.STATIONS = None

    def _init(self):
        Component._init(self)

        # convert to seconds
        self.RECORD_LENGTH_IN_SECONDS = self.inventory.record_length / second

        st = { 'forward':1, 'adjoint':2, 'both':3 }
        self.SIMULATION_TYPE = st[self.inventory.simulation_type]

        # Access our InputFile inventory items to make sure they're
        # readable.  (They will be reopened by the Fortran code.)
        f = self.inventory.cmtSolution;  self.checkCMTSolution(f);  self.CMTSOLUTION = f.name;  f.close()
        f = self.inventory.stations;                                self.STATIONS    = f.name;  f.close()
        f = self.inventory.stations;                                self.STATIONS    = f.name;  f.close()

    def checkCMTSolution(self, f):
        NLINES_PER_CMTSOLUTION_SOURCE = 13 # constants.h
        lineTally = 0
        for line in f:
            lineTally = lineTally + 1
        if lineTally % NLINES_PER_CMTSOLUTION_SOURCE != 0:
            raise ValueError("total number of lines in 'cmt-solution' file '%s' should be a multiple of %d"
                             % (f.name, NLINES_PER_CMTSOLUTION_SOURCE))
        NSOURCES = lineTally / NLINES_PER_CMTSOLUTION_SOURCE
        if NSOURCES < 1:
            raise ValueError("need at least one source in 'cmt-solution' file '%s'" % f.name)
        return


# end of file
