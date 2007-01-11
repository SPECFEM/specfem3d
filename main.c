
#include <Python.h>
#include <mpi.h>

#include "config.h"


extern void initPyxParameters(void);
extern void initPyxMeshfem(void);
#ifdef WITH_SOLVER
extern void initPyxSpecfem(void);
#endif

static int g_status;
int g_argc;
char **g_argv;

#define COMMAND \
"import sys; " \
"path = sys.argv[1]; " \
"requires = sys.argv[2]; " \
"entry = sys.argv[3]; " \
"path = path.split(':'); " \
"path.extend(sys.path); " \
"sys.path = path; " \
"from merlin import loadObject; " \
"entry = loadObject(entry); " \
"entry(sys.argv[3:], kwds={'requires': requires})"

/* include the implementation of _mpi */
#include "mpi/_mpi.c"

struct _inittab inittab[] = {
#ifdef WITH_MPI
    { "_mpi", init_mpi },
#endif
    { "PyxParameters", initPyxParameters },
    { "PyxMeshfem", initPyxMeshfem },
#ifdef WITH_SOLVER
    { "PyxSpecfem", initPyxSpecfem },
#endif
    { 0, 0 }
};


#define FC_PY_MAIN FC_FUNC_(fc_py_main, FC_PY_MAIN)
void FC_PY_MAIN()
{
    if (g_argc < 3 || strcmp(g_argv[1], "--pyre-start") != 0) {
        g_status = Py_Main(g_argc, g_argv);
        return;
    }
    
    /* make sure 'sys.executable' is set to the path of this program  */
    Py_SetProgramName(g_argv[0]);
    
    /* initialize Python */
    Py_Initialize();
    
    /* initialize sys.argv */
    PySys_SetArgv(g_argc - 1, g_argv + 1);
    
    /* run the Python command */
    g_status = PyRun_SimpleString(COMMAND) != 0;
    
    /* shut down Python */
    Py_Finalize();
}


int main(int argc, char **argv)
{
#if defined(WITH_MPI) && defined(USE_MPI)
    /* initialize MPI */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "%s: MPI_Init failed! Exiting ...", argv[0]);
        return 1;
    }
#endif
    
    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }
    
    g_argc = argc;
    g_argv = argv;
    
#define main 42
#if FC_MAIN == main
    /* start Python */
    FC_PY_MAIN();
#else
    /* call the Fortran trampoline (which, in turn, starts Python) */
    FC_MAIN();
#endif
    
#if defined(WITH_MPI) && defined(USE_MPI)
    /* shut down MPI */
    MPI_Finalize();
#endif
    
    return g_status;
}


/* end of file */
