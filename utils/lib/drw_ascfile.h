#ifndef _drw_ascfile_h
#define _drw_ascfile_h

void dread_ascfile(const char *ascfile,
                   double *t0, double *dt, int *n,
                   double *data);
void dwrite_ascfile(const char *ascfile,
                    double t0, double dt, int n,
                    const double *data);

void dread_ascfile_c(const char *ascfile,
                   double *t0, double *dt, int *n,
                   double *data);
void dwrite_ascfile_c(const char *ascfile,
                    const double *t0, const double *dt, const int *n,
                    const double *data);

void DREAD_ASCFILE_C(const char *ascfile,
                   double *t0, double *dt, int *n,
                   double *data);
void DWRITE_ASCFILE_C(const char *ascfile,
                    const double *t0, const double *dt, const int *n,
                    const double *data);





#endif
