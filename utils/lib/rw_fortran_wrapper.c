#include "drw_ascfile.h"

void dread_ascfile_c(const char *name, double *t0, double *dt, int * n, double *data) {
  
  dread_ascfile(name,t0,dt,n,data);}

void DREAD_ASCFILE_C(const char *name, double *t0, double *dt, int * n, double *data) {
  
  dread_ascfile(name,t0,dt,n,data);}

void dread_ascfile_c_(const char *name, double *t0, double *dt, int * n, double *data) {
  
  dread_ascfile(name,t0,dt,n,data);}

void DREAD_ASCFILE_C_(const char *name, double *t0, double *dt, int * n, double *data) {
  
  dread_ascfile(name,t0,dt,n,data);}

void dwrite_ascfile_c(const char *name, const double *t0, const double *dt, const int *n, const double *data) {
  
  dwrite_ascfile(name,*t0,*dt,*n,data);}

void DWRITE_ASCFILE_C(const char *name, const double *t0, const double *dt, const int *n, const double *data) {
  
  dwrite_ascfile(name,*t0,*dt,*n,data);}

void dwrite_ascfile_c_(const char *name, const double *t0, const double *dt, const int *n, const double *data) {
  
  dwrite_ascfile(name,*t0,*dt,*n,data);}

void DWRITE_ASCFILE_C_(const char *name, const double *t0, const double *dt, const int *n, const double *data) {
  
  dwrite_ascfile(name,*t0,*dt,*n,data);}
