#!/usr/bin/env python


from Script import Script


class Driver(Script):

    """Base class for scripts which call into SPECFEM3D Fortran code."""
    
    def __init__(self, outputFilename=None):
        super(Driver, self).__init__()
        self.outputFilename = outputFilename
        import Specfem3DBasinCode as code
        self.code = code
    
    def _init(self):
        super(Driver, self)._init()
        # make sure the output directory is writable
        if self.outputFilename and self.rank() == 0:
            self.checkOutputDir()
        return

    def rank(self):
        return 0
    
    def checkOutputDir(self):
        from os import remove
        from os.path import join
        outputDir = self.inventory.OUTPUT_FILES
        temp = join(outputDir, self.outputFilename)
        try:
            f = open(temp, 'w')
        except IOError:
            self.raisePropertyValueError('output-dir')
        f.close()
        remove(temp)
    
    def readValue(self, name):
        """Callback from Fortran into Python."""
        l = name.split('.')
        o = self
        for n in l:
            try:
                o = getattr(o, n)
            except AttributeError:
                o = getattr(o.inventory, n)
        return o


# end of file
