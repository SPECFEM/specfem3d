# Process this file with Pyrex to produce PyxParameters.c


"""Python bindings for the SPECFEM3D Global Solver."""


# include 'config.h' in order to get the definitions of FC_FUNC and FC_FUNC_
cdef extern from "config.h":
    pass


cdef extern from "string.h":
    char *strncpy(char *, char *, int)


cdef extern from "Python.h":
    object PyString_FromStringAndSize(char *, int)


# In the future, this could be passed through the Fortran layer as an
# opaque context argument.
component = None


def getValue(o, name):
    """Get a value from the Python scripts."""
    l = name.split('.')
    for n in l:
        o = getattr(o, n)
    return o


# replacements for Fortran functions

cdef public void read_value_integer "FC_FUNC_(read_value_integer, READ_VALUE_INTEGER)" (int *value, char *name, int nameLen) except *:
    attrName = PyString_FromStringAndSize(name, nameLen)
    value[0] = getValue(component, attrName)

cdef public void read_value_double_precision "FC_FUNC_(read_value_double_precision, READ_VALUE_DOUBLE_PRECISION)" (double *value, char *name, int nameLen) except *:
    attrName = PyString_FromStringAndSize(name, nameLen)
    value[0] = getValue(component, attrName)

cdef public void read_value_logical "FC_FUNC_(read_value_logical, READ_VALUE_LOGICAL)" (int *value, char *name, int nameLen) except *:
    attrName = PyString_FromStringAndSize(name, nameLen)
    value[0] = getValue(component, attrName)

cdef public void read_value_string "FC_FUNC_(read_value_string, READ_VALUE_STRING)" (char *value, char *name, int valueLen, int nameLen) except *:
    cdef char *vp
    cdef int vl, i
    attrName = PyString_FromStringAndSize(name, nameLen)
    v = getValue(component, attrName)
    vl = len(v)
    if vl > valueLen:
        raise ValueError("%s value '%s' is too long (%d bytes) for destination Fortran buffer (%d bytes)" % (attrName, v, vl, valueLen))
    vp = v
    strncpy(value, vp, vl)
    for i from vl <= i < valueLen:
        value[i] = c' '
    return

cdef public void open_parameter_file "FC_FUNC_(open_parameter_file, OPEN_PARAMETER_FILE)" () except *:
    return

cdef public void close_parameter_file "FC_FUNC_(close_parameter_file, CLOSE_PARAMETER_FILE)" () except *:
    return

cdef public void get_value_integer "FC_FUNC_(get_value_integer, GET_VALUE_INTEGER)" (int *value, char *name, int *default, int nameLen):
    value[0] = default[0]
    read_value_integer(value, name, nameLen)

cdef public void get_value_double_precision "FC_FUNC_(get_value_double_precision, GET_VALUE_DOUBLE_PRECISION)" (double *value, char *name, double *default, int nameLen):
    value[0] = default[0]
    read_value_double_precision(value, name, nameLen)
    
cdef public void get_value_logical "FC_FUNC_(get_value_logical, GET_VALUE_LOGICAL)" (int *value, char *name, int *default, int nameLen):
    value[0] = default[0]
    read_value_logical(value, name, nameLen)
    
cdef public void get_value_string "FC_FUNC_(get_value_string, GET_VALUE_STRING)" (char *value, char *name, char *default, int valueLen, int nameLen, int defaultLen):
    if defaultLen > valueLen:
        strncpy(value, default, valueLen)
    else:
        strncpy(value, default, defaultLen)
        for i from defaultLen <= i < valueLen:
            value[i] = c' '
    read_value_string(value, name, valueLen, nameLen)


# external Fortran functions

cdef extern void create_header_file_f "FC_FUNC_(create_header_file, CREATE_HEADER_FILE)" () except *
def create_header_file(arg):
    """Create the include file for the solver."""
    global component
    component = arg
    create_header_file_f()


# end of file
