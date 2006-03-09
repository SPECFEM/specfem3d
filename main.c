
#include <Python.h>
#include <stdio.h>
#include "config.h"


extern void initSpecfem3DBasinCode(void);
extern void initPyxMPI(void);

static int status;

struct _inittab inittab[] = {
    { "Specfem3DBasinCode", initSpecfem3DBasinCode },
#ifdef USE_MPI
    { "PyxMPI", initPyxMPI },
#endif
    { 0, 0 }
};


#define FC_RUN_PYTHON_SCRIPT FC_FUNC_(run_python_script, RUN_PYTHON_SCRIPT)
void FC_RUN_PYTHON_SCRIPT()
{
    /* run the Python script */
#ifndef SCRIPT
#define SCRIPT Specfem
#endif
#define STR(s) #s
#define COMMAND(s) "from Specfem3DBasin."STR(s)" import "STR(s)"; app = "STR(s)"(); app.run()"
    status = PyRun_SimpleString(COMMAND(SCRIPT)) != 0;
}


int main(int argc, char **argv)
{
    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }
    
    /* initialize Python */
    Py_Initialize();
    
    /* initialize sys.argv */
    PySys_SetArgv(argc, argv);
    
#define main 42
#if FC_MAIN == main
    /* run the Python script */
    FC_RUN_PYTHON_SCRIPT();
#else
    /* call the Fortran trampoline (which runs the Python script) */
    FC_MAIN();
#endif
    
    /* shut down Python */
    Py_Finalize();
    
    return status;
}


void xxxxfem3D_dispatch()
{
#define Meshfem 42
#define Specfem 24
#if SCRIPT == Meshfem
    FC_FUNC(meshfem3d, MESHFEM3D)();
#else
    FC_FUNC(specfem3d, SPECFEM3D)();
#endif
}


/* end of file */
