
#include <Python.h>
#include "config.h"

/* called from Fortran to propagate Python exceptions */
int FC_FUNC_(err_occurred, ERR_OCCURRED)()
{
    return PyErr_Occurred() ? 1 : 0;
}

/* end of file */
