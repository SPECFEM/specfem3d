#!/usr/bin/env python


from mpi.Launcher import Launcher
from pyre.util import expandMacros
from pyre.inventory.Facility import Facility
import os, sys


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def which(filename):
    from os.path import abspath, exists, join
    from os import environ, pathsep
    dirs = environ['PATH'].split(pathsep)
    for dir in dirs:
       pathname = join(dir, filename)
       if exists(pathname):
           return abspath(pathname)
    return filename


def hms(t):
    return (int(t / 3600), int((t % 3600) / 60), int(t % 60))


defaultEnvironment = "[EXPORT_ROOT,LD_LIBRARY_PATH,PYTHONPATH]"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MPI Launchers:
#     (Replacement) Launcher for MPICH
#     Launcher for LAM/MPI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class LauncherMPI(Launcher):


    class Inventory(Launcher.Inventory):

        from pyre.inventory import bool, str

        dry = bool("dry", default=False)
        dry.meta['tip'] = "prints the command line and exits"
        
        debug = bool("debug", default=False)

        Launcher.Inventory.nodes.meta['tip'] = """number of machine nodes"""
        Launcher.Inventory.nodelist.meta['tip'] = """a comma-separated list of machine names in square brackets (e.g., [101-103,105,107])"""
        nodegen = str("nodegen")
        nodegen.meta['tip'] = """a printf-style format string, used in conjunction with 'nodelist' to generate the list of machine names (e.g., "n%03d")"""
        
        extra = str("extra")
        extra.meta['tip'] = "extra arguments to pass to mpirun"
        
        command = str("command", default="mpirun")
        re_exec = bool("re-exec", default=True)


    def launch(self):
        args = self._buildArgumentList()
        if not args:
            return self.inventory.dry
        
        command = " ".join(args)
        self._info.log("executing: {%s}" % command)

        if self.inventory.dry:
            print command
            return True
        
        os.system(command)
        return True

            
    def _buildArgumentList(self):
        if not self.nodes:
            self.nodes = len(self.nodelist)

        if self.nodes < 2:
            self.inventory.nodes = 1
            return []

        if False:
            import mpi
            if mpi.world().handle():
                self.inventory.nodes = 1
                return []
            elif self.inventory.re_exec:
                # re-exec under mpipython.exe
                args = list(sys.argv)
                args.append("--launcher.re-exec=False") # protect against infinite regress
                return args
            else:
                assert False
        
        # build the command
        args = []
        args.append(self.inventory.command)
        self._appendMpiRunArgs(args)

        args += sys.argv
        args.append("--mode=worker")

        return args

    
    def _appendMpiRunArgs(self, args):
        args.append(self.inventory.extra)
        args.append("-np %d" % self.nodes)
        
        # use only the specific nodes specified explicitly
        if self.nodelist:
            self._appendNodeListArgs(args)


class LauncherMPICH(LauncherMPI):


    class Inventory(LauncherMPI.Inventory):
 
        from pyre.inventory import str

        machinefile = str("machinefile", default="mpirun.nodes")
        machinefile.meta['tip'] = """filename of machine file"""


    def __init__(self):
        LauncherMPI.__init__(self, "mpich")


    def _appendNodeListArgs(self, args):
        machinefile = self.inventory.machinefile
        nodegen = self.inventory.nodegen
        file = open(machinefile, "w")
        for node in self.nodelist:
            file.write((nodegen + '\n') % node)
        file.close()
        args.append("-machinefile %s" % machinefile)


class LauncherLAMMPI(LauncherMPI):


    class Inventory(LauncherMPI.Inventory):

        from pyre.inventory import list

        environment = list("environment", default=defaultEnvironment)
        environment.meta['tip'] = """a comma-separated list of environment variables to export to the batch job"""


    def __init__(self):
        LauncherMPI.__init__(self, "lam-mpi")


    def _appendMpiRunArgs(self, args):
        args.append("-x %s" % ','.join(self.inventory.environment))
        super(LauncherLAMMPI, self)._appendMpiRunArgs(args)


    def _appendNodeListArgs(self, args):
        nodegen = self.inventory.nodegen
        args.append("n" + ",".join([(nodegen) % node for node in self.nodelist]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are Launchers for batch schedulers found on the TeraGrid.
# These should be incorporated into Pythia eventually.

# With something like StringTemplate by Terence Parr and Marq Kole,
# the batch scripts could be generated entirely from an
# inventory-data-driven template.

#     http://www.stringtemplate.org/doc/python-doc.html

# This code uses a hybrid approach, mixing Python logic with primitive
# templates powered by pyre.util.expandMacros().
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class LauncherBatch(Launcher):


    class Inventory(Launcher.Inventory):

        from pyre.inventory import bool, dimensional, list, str
        from pyre.units.time import minute

        dry = bool("dry", default=False)
        debug = bool("debug", default=False)

        task = str("task")

        # Ignore 'nodegen' so that the examples will work without modification.
        nodegen = str("nodegen")
        nodegen.meta['tip'] = """(ignored)"""

        walltime = dimensional("walltime", default=0*minute)
        mail = bool("mail", default=False)
        queue = str("queue")

        directory = str("directory", default="${cwd}")
        script = str("script", default="${directory}/script")
        stdout = str("stdout", default="${directory}/stdout")
        stderr = str("stderr", default="${directory}/stderr")
        environment = list("environment", default=defaultEnvironment)

        cwd = str("cwd", default=os.getcwd())
        argv = str("argv", default=(' '.join(sys.argv)))


    def launch(self):
        if self.inventory.dry:
            print self._buildScript()
            print "# submit with:"
            print "#", self._buildBatchCommand()
            return True
        
        # write the script
        scriptFile = open(expandMacros("${script}", self.inv), "w")
        scriptFile.write(self._buildScript())
        scriptFile.close()

        # build the batch command
        command = self._buildBatchCommand()
        self._info.log("executing: {%s}" % command)

        os.system(command)
        return True


    def __init__(self, name):
        Launcher.__init__(self, name)
        
        # Used to recursively expand ${macro) in format strings using my inventory.
        class InventoryAdapter(object):
            def __init__(self, launcher):
                self.launcher = launcher
            def __getitem__(self, key):
                return expandMacros(str(self.launcher.inventory.getTraitValue(key)), self)
        self.inv = InventoryAdapter(self)
        
        return


    def _buildBatchCommand(self):
        return expandMacros("${batch-command} ${script}", self.inv)


    def _buildScript(self):
        script = [
            "#!/bin/sh",
            ]
        self._buildScriptDirectives(script)
        script += [
            expandMacros('''\

cd ${directory}
${command} ${argv} --mode=worker
''', self.inv)
            ]
        script = "\n".join(script) + "\n"
        return script


# Note: mpi.LauncherPBS in Pythia-0.8 does not work!

class LauncherPBS(LauncherBatch):


    class Inventory(LauncherBatch.Inventory):
        
        from pyre.inventory import str
        
        command = str("command", default="mpirun -np ${nodes} -machinefile $PBS_NODEFILE") # Sub-launcher?
        batch_command = str("batch-command", default="qsub")


    def __init__(self):
        LauncherBatch.__init__(self, "pbs")


    def _buildScriptDirectives(self, script):
        
        queue = self.inventory.queue
        if queue:
            script.append("#PBS -q %s" % queue)

        task = self.inventory.task
        if task:
            script.append("#PBS -N %s" % task)

        if self.inventory.stdout:
            script.append(expandMacros("#PBS -o ${stdout}", self.inv))
        if self.inventory.stderr:
            script.append(expandMacros("#PBS -e ${stderr}", self.inv))

        resourceList = self._buildResourceList()

        script += [
            "#PBS -V", # export qsub command environment to the batch job
            "#PBS -l %s" % resourceList,
            ]

        return script



    def _buildResourceList(self):

        resourceList = [
            "nodes=%d" % self.nodes,
            ]

        walltime = self.inventory.walltime.value
        if walltime:
            resourceList.append("walltime=%d:%02d:%02d" % hms(walltime))

        resourceList = ",".join(resourceList)

        return resourceList


class LauncherLSF(LauncherBatch):


    class Inventory(LauncherBatch.Inventory):
        
        from pyre.inventory import list, str
        
        command = str("command", default="mpijob mpirun")
        batch_command = str("batch-command", default="bsub")
        bsub_options = list("bsub-options")


    def __init__(self):
        LauncherBatch.__init__(self, "lsf")


    def _buildBatchCommand(self):
        return expandMacros("${batch-command} < ${script}", self.inv)


    def _buildScriptDirectives(self, script):

        # LSF scripts must have a job name; otherwise strange "/bin/sh: Event not found"
        # errors occur (tested on TACC's Lonestar system).
        task = self.inventory.task
        if not task:
            task = "jobname"
        script.append("#BSUB -J %s" % task)
        
        queue = self.inventory.queue
        if queue:
            script.append("#BSUB -q %s" % queue)

        walltime = self.inventory.walltime.value
        if walltime:
            script.append("#BSUB -W %d:%02d" % hms(walltime)[0:2])
        
        if self.inventory.stdout:
            script.append(expandMacros("#BSUB -o ${stdout}", self.inv))
        if self.inventory.stderr:
            script.append(expandMacros("#BSUB -e ${stderr}", self.inv))
            
        script += [
            "#BSUB -n %d" % self.nodes,
            ]

        script += ["#BSUB " + option for option in self.inventory.bsub_options]

        return script


class LauncherGlobus(LauncherBatch):


    class Inventory(LauncherBatch.Inventory):

        from pyre.inventory import str
        
        batch_command = str("batch-command", default="globusrun")
        executable = str("executable", default=sys.argv[0])
        resource = str("resource", default="localhost")


    def _buildBatchCommand(self):
        return expandMacros("${batch-command} -b -r ${resource} -f ${script}", self.inv)


    def __init__(self):
        LauncherBatch.__init__(self, "globus")


    def _buildScript(self):
        script = [
            expandMacros('''\
&   (jobType=mpi)
    (executable="${executable}")
    (count=${nodes})
    (directory="${directory}")
    (stdout="${stdout}")
    (stderr="${stderr}")''', self.inv),
            ]
        
        script.append('    (environment = %s)' % self._buildEnvironment())

        # add the arguments
        args = sys.argv[1:]
        args.append("--mode=worker")
        command = '    (arguments= ' + ' '.join([('"%s"' % arg) for arg in args]) + ')'
        script.append(command)

        script = '\n'.join(script) + '\n'

        return script


    def _buildEnvironment(self):
        from os import environ
        #vars = environ.keys()
        vars = self.inventory.environment
        env = [('(%s "%s")' % (var, environ.get(var,""))) for var in vars]
        env = ' '.join(env)
        return env


# main
if __name__ == "__main__":

    
    from pyre.applications.Script import Script

    
    class TestApp(Script):

        
        class Inventory(Script.Inventory):
            
            from pyre.inventory import facility
            from mpi.Launcher import Launcher

            launcher = facility("launcher", default="mpich")


        def main(self, *args, **kwds):
            launcher = self.inventory.launcher
            if launcher:
                try:
                    # batch launcher
                    print launcher._buildScript()
                    print
                    print "# submit with", launcher._buildBatchCommand()
                except AttributeError:
                    # direct launcher
                    print ' '.join(launcher._buildArgumentList())
            return

        
        def __init__(self):
            Script.__init__(self, "CitcomS")

    
    app = TestApp()
    app.run()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LauncherFacility
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


builtInLauncherClasses = {
    'globus':  LauncherGlobus,
    'lam-mpi': LauncherLAMMPI,
    'lsf':     LauncherLSF,
    'mpich':   LauncherMPICH,
    'pbs':     LauncherPBS,
    }

class LauncherFacility(Facility):

    def __init__(self, name):
        Facility.__init__(self, name, default="lsf")
        return
    
    def _retrieveComponent(self, instance, componentName):
        cls = builtInLauncherClasses.get(componentName, None)
        if cls is None:
            return Facility._retrieveComponent(self, instance, componentName)
        launcher = cls()
        launcher.aliases.append(self.name)
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('built-in')
        return launcher, locator


# end of file 
