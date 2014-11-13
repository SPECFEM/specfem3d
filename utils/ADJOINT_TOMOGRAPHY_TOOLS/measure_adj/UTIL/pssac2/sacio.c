/*******************************************************************
*			sacio.c
* SAC I/O functions:
*	read_sachead	read SAC header
*	read_sac	read SAC binary data
*	read_sac2	read SAC data with cut option
*	write_sac	write SAC binary data
*	wrtsac0		write 1D array as evenly-spaced SAC binary data
*	wrtsac0_	fortran wrap for wrtsac0
*	wrtsac2		write 2 1D arrays as XY SAC data
*	wrtsac2_	fortran wrap for wrtsac2
*	swab4		reverse byte order for integer/float
*	rtrend		remove trend
*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SAC_NULL
#include "sac.h"


/***********************************************************

  read_sachead

  Description:	read binary SAC header from file.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled

  Return:	0 if success, -1 if failed

  Modify history:
	05/29/97	Lupei Zhu	Initial coding
************************************************************/

int	read_sachead(const char	*name,
		SACHEAD		*hd
	)
{
  FILE		*strm;
  
  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return -1;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     fclose(strm);
     return -1;
  }

#ifdef i386
  swab4((char *) hd, HD_SIZE);
#endif

  fclose(strm);
  return 0;

}


/***********************************************************

  read_sac

  Description:	read binary data from file. If succeed, it will return
		a float pointer to the read data array. The SAC header
		is also filled. A NULL pointer is returned if failed.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled

  Return:	float pointer to the data array, NULL if failed

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
************************************************************/

float*	read_sac(const char	*name,
		SACHEAD		*hd
	)
{
  FILE		*strm;
  float		*ar;
  unsigned	sz;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

#ifdef i386
  swab4((char *) hd, HD_SIZE);
#endif

  sz = hd->npts*sizeof(float);
  if ((ar = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     return NULL;
  }

  if (fread((char *) ar, sz, 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }

  fclose(strm);

#ifdef i386
  swab4((char *) ar, sz);
#endif

  return ar;

}

/* sac_byte_order (hd)
   Determine if the header variable (internal4)
   which defines the header version is between 0 and 6
   do_swap  = 1 Try and byte swap
            = 0 Do not byte swap
   Returns:
              0 - No swapping needs to be done
	      1 - Swapping is needed
*/
int
sac_byte_order(SACHEAD *hd) {
  if( ( hd->internal4 > 0  ) &&
      ( hd->internal4 <= 6 ) ) {
    return(0);
  }
  return(1);
}


float*	read_sac_swap(const char	*name,
		      SACHEAD	*hd,
		      int           do_swap
	)
{
  FILE		*strm;
  float		*ar;
  unsigned	sz;
  int           swap;
  
  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     fclose(strm);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     fclose(strm);
     return NULL;
  }
  swap = sac_byte_order(hd);
  if(swap) {
    if(do_swap) {
      swab4((char *) hd, HD_SIZE);
      if(sac_byte_order(hd)) {
	fprintf(stderr, "Error deciphering SAC header %s\n", name);
	fclose(strm);
	return(NULL);
      }
    } else {
      fprintf(stderr, "Error deciphering SAC header %s\n", name);
      fclose(strm);
      return(NULL);
    }
  }

  
  sz = hd->npts*sizeof(float);
  if ((ar = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     fclose(strm);
     return NULL;
  }

  if (fread((char *) ar, sz, 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     fclose(strm);
     free(ar);
     return NULL;
  }

  fclose(strm);
  if(swap) {
    if(do_swap) {
      swab4((char *) ar, sz);
    } else{
      fprintf(stderr, "I am not really certain how I got here\n");
      fprintf(stderr, "%s needs the data swapped but I was told not to\n",name);
      fclose(strm);
      free(ar);
      return(NULL);
    }
  }

  return ar;
  
}


/***********************************************************

  write_sac

  Description:	write SAC binary data.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD hd		SAC header
		const float *ar		pointer to the data

  Return:	0 if succeed; -1 if failed

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
************************************************************/

int	write_sac(const char	*name,
		SACHEAD		hd,
		const float	*ar
	)
{
  FILE		*strm;
  unsigned	sz;
  float		*data;
  int		error = 0;

  strm = NULL;
  sz = hd.npts*sizeof(float);
  if (hd.iftype == IXY) sz *= 2;

  if ((data = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for writing %s\n",name);
     error = 1;
  }

  if ( !error && memcpy(data, ar, sz) == NULL) {
     fprintf(stderr, "Error in copying data for writing %s\n",name);
     error = 1;
  }

#ifdef i386
  swab4((char *) data, sz);
  swab4((char *) &hd, HD_SIZE);
#endif

  if ( !error && (strm = fopen(name, "w")) == NULL ) {
     fprintf(stderr,"Error in opening file for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(&hd, sizeof(SACHEAD), 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC header for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(data, sz, 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC data for writing %s\n",name);
     error = 1;
  }

  free(data);
  fclose(strm);

  return (error==0) ? 0 : -1;

}



/*****************************************************

  swab4

  Description:	reverse byte order for float/integer

  Author:	Lupei Zhu

  Arguments:	char *pt	pointer to byte array
		int    n	number of bytes

  Return:	none

  Modify history:
	12/03/96	Lupei Zhu	Initial coding

************************************************************/

void	swab4(	char	*pt,
		int	n
	)
{
  int i;
  char temp;
  for(i=0;i<n;i+=4) {
    temp = pt[i+3];
    pt[i+3] = pt[i];
    pt[i] = temp;
    temp = pt[i+2];
    pt[i+2] = pt[i+1];
    pt[i+1] = temp;
  }
}



/***********************************************************

  wrtsac0

  Description:	write 1D array into evenly spaced SAC data.

  Author:	Lupei Zhu

  Arguments:	const char *name	file name
		float	dt		sampling interval
		int	ns		number of points
		float	b0		starting time
		float	dist		distance range
		const float *ar		data array

  Return:	0 if succeed; -1 if failed

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swab byte-order on PC
************************************************************/

int    wrtsac0(const char	*name,
		float		dt,
		int		ns,
		float		b0,
		float		dist,
		const float	*ar
	)
{

  SACHEAD	hd = sac_null;

  hd.npts = ns;
  hd.delta = dt;
  hd.dist = dist;
  hd.b = b0;
  hd.o = 0.;
  hd.e = b0+(hd.npts-1)*hd.delta;
  hd.iztype = IO;
  hd.iftype = ITIME;
  hd.leven = TRUE;

  return write_sac(name, hd, ar);

}
 


/***********************************************************

  wrtsac2

  Description:	write 2 arrays into XY SAC data.

  Author:	Lupei Zhu

  Arguments:	const char *name	file name
		int	ns		number of points
		const float *x		x data array
		const float *y		y data array

  Return:	0 succeed, -1 fail

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
************************************************************/

int	wrtsac2(const char	*name,
		int		n,
		const float	*x,
		const float	*y
	)
{
  SACHEAD	hd = sac_null;
  float		*ar;
  unsigned	sz;
  int		exit_code;

  hd.npts = n;
  hd.iftype = IXY;
  hd.leven = FALSE;

  sz = n*sizeof(float);

  if ( (ar = (float *) malloc(2*sz)) == NULL ) {
     fprintf(stderr, "error in allocating memory%s\n",name);
     return -1;
  }

  if (memcpy(ar, x, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }
  if (memcpy(ar+sz, y, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }

  exit_code = write_sac(name, hd, ar);
  
  free(ar);
  
  return exit_code;

}


/* write user data t[n] into sac header starting at pos (0-9)*/
void sacUdata(float *hd, int pos, float *t, int n, int type)
{
  int i;
  hd += type;
  for(i=pos;i<n && i<10;i++) hd[i] = t[i];
}


/* for fortran--write evenly-spaced data */
void    wrtsac0_(const char *name, float dt, int ns, float b0, float dist, const float *ar) {
  wrtsac0(name, dt, ns, b0, dist, ar);
}

/* for fortran--write x-y data */
void    wrtsac2_(const char *name, int n, const float *x, const float *y) {
  wrtsac2(name, n, x, y);
}


/***********************************************************

  read_sac2

  Description:	read portion ofdata from file. If succeed, it will return
		a float pointer to the read data array. The SAC header
		is also filled. A NULL pointer is returned if failed.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled
		int	tmark,		time mark in sac header
		float	t1		begin time is tmark + t1
		int	npts		number of points

  Return:	float pointer to the data array, NULL if failed

  Modify history:
	11/08/00	Lupei Zhu	Initial coding
************************************************************/

float*	read_sac2(const char	*name,
		SACHEAD		*hd,
		int		tmark,
		float		t1,
		int		npts
	)
{
  FILE		*strm;
  int		i, nt1, nt2;
  float		*ar, *fpt;
  void		rtrend(float *, int);

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

#ifdef i386
  swab4((char *) hd, HD_SIZE);
#endif

  t1 += *( (float *) hd + 10 + tmark);
  nt1 = (int) rint((t1-hd->b)/hd->delta);
  nt2 = nt1+npts-1;
  if (nt1>=hd->npts-1 || nt2<=0) {
     fprintf(stderr,"data not in the specified window %s\n",name);
     return NULL;
  }
  if ((ar = (float *) malloc(npts*sizeof(float))) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     return NULL;
  }
  for(i=0;i<npts;i++) ar[i]=0.;

  if (nt1 > 0 ) {
     if (fseek(strm,nt1*sizeof(float),SEEK_CUR) < 0) {
	fprintf(stderr, "error in seek %s\n",name);
	return NULL;
     }
     fpt = ar;
  } else {
     fpt = ar-nt1;
     nt1 = 0;
  }
  if (nt2>hd->npts-1) nt2=hd->npts-1;
  i = nt2-nt1+1;
  if (fread((char *) fpt, sizeof(float), i, strm) != i) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }
  fclose(strm);

#ifdef i386
  swab4((char *) fpt, i*sizeof(float));
#endif

  /* remove trend */
  rtrend(fpt, i);

  hd->npts = npts;
  hd->b = t1;
  hd->e = hd->b+npts*hd->delta;
  return ar;

}

/* remove trend a*i + b */
void	rtrend(float *y, int n) {
     int i;
     double a, b, a11, a12, a22, y1, y2;
     y1 = y2 = 0.;
     for(i=0;i<n;i++) {
       y1 += i*y[i];
       y2 += y[i];
     }
     a12 = 0.5*n*(n-1);
     a11 = a12*(2*n-1)/3.;
     a22 = n;
     b = a11*a22-a12*a12;
     a = (a22*y1-a12*y2)/b; 
     b = (a11*y2-a12*y1)/b;
     for(i=0;i<n;i++) {
       y[i] = y[i] - a*i - b;
     }
}
