
#include <string.h>

/* this function is missing from libsac.a v101.3 */
void copy_float(float *src, float *dest, int n)
{
    memmove(dest, src, n * sizeof(float));
}
