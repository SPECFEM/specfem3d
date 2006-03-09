#!/usr/bin/env python


# courtesy Bob Ippolito  bob@redivi.com 

class MetaCategory(type):
     def __new__(cls, name, bases, dct):
         if '__category_of__' in dct:
             return type.__new__(cls, name, bases, dct)
         if not len(bases) == 1 and isinstance(bases[0], cls):
             raise TypeError("Categories may only have a Category(...) as their base")
         cls = bases[0].__category_of__
         for k,v in dct.iteritems():
             if k == '__module__' or k == '__doc__':
                 continue
             setattr(cls, k, v)
         return cls

def Category(cls):
     return MetaCategory(
         'Category(%s)' % (cls.__name__),
         (object,),
         {'__category_of__': cls}
     )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Property


from pyre.inventory.Property import Property
from pyre.inventory.Trait import Trait
from journal.devices.Renderer import Renderer as PyreRenderer
import sys


class PropertyValueError(ValueError):
     

     class Renderer(PyreRenderer):
          def __init__(self, header=None, format=None, footer=None):
               PyreRenderer.__init__(
                    self,
                    header="%(filename)s:%(line)s: property '%(name)s': %(error)s\n >> %(src)s",
                    format=format, footer=footer
                    )


     def __init__(self, *args):
          ValueError.__init__(self, *args)
          self.origExc = sys.exc_info()


     def __str__(self):
          return self.origExc[0].__name__ + ": " + self.origExc[1].__str__()


class PropertyPatches(Category(Property)):

    
    """Add location info to ValueError exceptions.
    
    Instead of _cast, call _richCast;
    instead of validator, call _richValidator.
    
    """
    

    def _set(self, instance, value, locator):
        # None is a special value; it means that a property is not set
        if value is not None:
            # convert
            value = self._richCast(value, locator)
            # validate 
            value = self._richValidator(value, locator)

        # record
        return Trait._set(self, instance, value, locator)


    def _getDefaultValue(self, instance):
        """retrieve the default value and return it along with a locator"""

        value = self.default

        # None is a special value and shouldn't go through the _cast
        if value is not None:
            # convert
            value = self._richCast(value)
            # validate
            value = self._richValidator(value)
        
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('default')

        return value, locator

    
    def _richCast(self, value, locator=None):
        try:
            value = self._cast(value)
        except:
             raise PropertyValueError, (self, locator)
        else:
             return value
    

    def _richValidator(self, value, locator=None):
        if self.validator:
            try:
                value = self.validator(value)
            except ValueError:
                raise PropertyValueError, (self, locator)
        return value


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mpi.Application


from mpi.Application import Application as MPIApp


class MPIAppPatches(Category(MPIApp)):


    def run(self, *args, **kwds):
        
        # If we are a worker on a compute node -- i.e., we are a
        # process started by 'mpirun' -- call MPI_Init() now to clean
        # the MPI arguments from sys.argv before Pyre can see them.
        # If we are the launcher/"server" on the login node -- i.e.,
        # we are the process started directly by the user -- don't
        # call MPI_Init(), as MPICH-GM will die with SIGPIPE
        # ("<MPICH-GM> Error: Need to obtain the job magic number in
        # GMPI_MAGIC !").
        
        scratchRegistry = self.createRegistry()
        help, unprocessedArguments = self.processCommandline(scratchRegistry)
        mode = scratchRegistry.getProperty('mode', "server")
        
        if mode == "worker":
            from PyxMPI import MPI_Init, MPI_Finalize
            MPI_Init(sys.argv)
            super(MPIApp, self).run(self, *args, **kwds)
            MPI_Finalize()
        else:
            super(MPIApp, self).run(self, *args, **kwds)
        
        return

    
    def onServer(self, *args, **kwds):
        self._debug.log("%s: onServer" % self.name)

        launcher = self.inventory.launcher
        launched = launcher.launch()
        if not launched:
            # only one node -- nothing to launch
            from PyxMPI import MPI_Init, MPI_Finalize
            MPI_Init(sys.argv)
            self.onComputeNodes(*args, **kwds)
            MPI_Finalize()
        
        return



# end of file
