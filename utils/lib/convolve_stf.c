/*
 * convolve_stf.c by Qinya Liu, Eh Tan
 *
 * Copyright (C) 2003-2009, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "sacio.h"
#include "config.h"


/* defined in libsac */
void fft(float *xreal, float *ximag, sac_int_t n, sac_int_t idir);


static void fzero(float *dst, sac_int_t n) {
    while (n)
        dst[--n] = 0.0;
}

static void fcpy(float *dst, const float *src, sac_int_t n) {
    while (n) {
        --n;
        dst[n] = src[n];
    }
}

void convolve(float **pconv, sac_int_t *pnconv,
              const float *data, sac_int_t ndata,
              const float *stf,  sac_int_t nstf)
{
    sac_int_t nconv, ncorr, i;
    struct { float *xreal, *ximag; } cdata, cstf, ccorr;
    float *conv, *buffer;

    nconv = ndata + nstf - 1;
    conv = (float *)malloc(nconv * sizeof(float));

    for (ncorr = 2; ncorr < nconv; ncorr *= 2)
        ;

    buffer = (float *)malloc(2 * 3 * ncorr * sizeof(float));
    cdata.xreal = buffer + 0 * ncorr;
    cdata.ximag = buffer + 1 * ncorr;
    cstf.xreal  = buffer + 2 * ncorr;
    cstf.ximag  = buffer + 3 * ncorr;
    ccorr.xreal = buffer + 4 * ncorr;
    ccorr.ximag = buffer + 5 * ncorr;

    fcpy(cdata.xreal, data, ndata);
    fzero(cdata.xreal + ndata, ncorr - ndata);
    fzero(cdata.ximag, ncorr);

    fcpy(cstf.xreal, stf, nstf);
    fzero(cstf.xreal + nstf, ncorr - nstf);
    fzero(cstf.ximag, ncorr);

    fft(cdata.xreal, cdata.ximag, ncorr, 1);
    fft(cstf.xreal, cstf.ximag, ncorr, 1);

    for (i = 0; i < ncorr; ++i) {
        /* ccorr[i] = cdata[i] * cstf[i] */
        ccorr.xreal[i] = cdata.xreal[i] * cstf.xreal[i] - cdata.ximag[i] * cstf.ximag[i];
        ccorr.ximag[i] = cdata.xreal[i] * cstf.ximag[i] + cdata.ximag[i] * cstf.xreal[i];
    }

    fft(ccorr.xreal, ccorr.ximag, ncorr, -1);
    for (i = 0; i < nconv; ++i) {
        conv[i] = ccorr.xreal[i] * ncorr;
    }

    free(buffer);

    *pconv = conv;
    *pnconv = nconv;

    return;
}


int
main(int argc, char *argv[])
{
    char cstf, *endpt;
    float hdur, *data;
    int j;
    sac_int_t datasize;
    const int min_nhdur = 10;

    if(argc < 4) {
        fprintf(stderr,
                "Usage: %s g[|t] HDUR FILE...\n"
                "  This program convolves SAC files with a gaussian(|triangle)\n"
                "  source time function of given half duration\n",
                argv[0]);
        return 1;
    }

    if(strcmp(argv[1], "t") != 0  && strcmp(argv[1], "g") != 0) {
        fprintf(stderr,"The source time function type could only be triangle(t) or gaussian(g) \n");
        return 1;
    }
    cstf = argv[1][0];

    errno = 0;
    hdur = (float)strtod(argv[2], &endpt);
    if(errno || endpt == argv[2] || *endpt != 0) {
        fprintf(stderr,"No floating point number can be formed from %s\n",argv[2]);
        return 1;
    }


    /* loop over input sac files */
    datasize = 0;
    data = NULL;
    for (j=3; j<argc; j++) {
        sac_int_t max, npts, nlen, nerr;
        float beg, del, dt, origin, tmp[1];
        sac_int_t nstf, nconv, i;
        float hstf, *stf, *conv;
        char *outf;


        fprintf(stderr, "convolving sac file %s with half duration %7.3f\n",
                argv[j], hdur);

        if((outf = (char *) malloc(strlen(argv[j]) + strlen(".conv") + 1)) == NULL) {
            fprintf(stderr,"Out of memory\n");
            return 1;
        }
        strcpy(outf, argv[j]);
        strcat(outf, ".conv");

        /* read header to get the length of time series */
        max = 1;
        rsac1(argv[j], tmp, &nlen, &beg, &del, &max, &nerr, strlen(argv[j]));
        if(nerr != -803) {
            fprintf(stderr,"Error reading sac file %s\n",argv[j]);
            return 1;
        }

        getnhv("npts", &npts, &nerr, strlen("npts"));

        /* now we know how much memory we need to allocate */
        if(npts > datasize) {
            data = (float*)realloc(data, npts * sizeof(float));
            if(data == NULL) {
                fprintf(stderr, "out of memory\n");
                return 1;
            }
            datasize = npts;
        }

        /* read the data */
        max = npts;
        rsac1(argv[j], data, &nlen, &beg, &del, &max, &nerr, strlen(argv[j]));
        if(nerr) {
            fprintf(stderr,"Error reading sac file %s\n",argv[j]);
            return 1;
        }

        /* get additional info */
        getfhv("delta", &dt, &nerr, strlen("delta"));

        getfhv("o", &origin, &nerr, strlen("o"));
        if(nerr) {
            /* Assuming origin time is 0, per convention of SPECFEM */
      fprintf(stderr, " Ignoring the error, assuming origin time is 0\n");
            origin = 0.0;
            setfhv("o", &origin, &nerr, strlen("o"));
        }



        /* creat source time function time series */
        if(min_nhdur * dt / 2 > hdur) {
            fprintf(stderr,"The half duration %f is too small to convolve\n", hdur);
            return 1;
        }

        nstf = (int)ceil(2 * hdur/dt)+1;
        hstf = (nstf-1)*dt/2;

        if((stf = (float *) malloc(nstf*sizeof(float))) == NULL) {
            fprintf(stderr,"Error in allocating memory for source time function\n");
            return 1;
        }
        if(cstf == 't') {
            /* triangular */
            for (i=0; i<nstf/2; i++)
                stf[i] = i * dt / (hstf * hstf);
            for (i=nstf/2; i<nstf; i++)
                stf[i] = (2 * hstf - i * dt)/ (hstf*hstf);
        } else {
            /* gaussian */
            const float decay_rate = 1.628;
            float alpha = decay_rate / hdur;
            float divisor = sqrt(4*atan(1.0) /*pi*/);
            for (i=0; i<nstf; i++) {
                float tao_i = fabs(i*dt - hstf);
                stf[i] = alpha * exp(- alpha*alpha* tao_i*tao_i) / divisor;
            }
        }


        /* creat convolution time series */
        convolve(&conv,&nconv,data,npts,stf,nstf);
        for(i=0; i<nconv; i++)
            conv[i] *= dt;


        /* update sac file header */
        {
            int n;
            float bnew, enew, aminm, amaxm, amean;

            bnew = origin - hstf + beg;
            enew = bnew + del * nconv;

            aminm = amaxm = conv[0];
            amean = 0;
            for(n=0; n<nconv; n++) {
                aminm = (aminm < conv[n]) ? aminm : conv[n];
                amaxm = (amaxm > conv[n]) ? amaxm : conv[n];
                amean += conv[n];
            }
            amean /= nconv;

            setfhv("b", &bnew, &nerr, strlen("b"));
            setfhv("e", &enew, &nerr, strlen("e"));
            setnhv("npts", &nconv, &nerr, strlen("npts"));
            setfhv("depmin", &aminm, &nerr, strlen("depmin"));
            setfhv("depmax", &amaxm, &nerr, strlen("depmax"));
            setfhv("depmen", &amean, &nerr, strlen("depmen"));
        }

        /* output to .conv sac file */
        wsac0(outf, data, conv, &nerr, strlen(outf));

        if(nerr) {
            fprintf(stderr, "Not able to write to file %s\n", argv[j]);
            return 1;
        }

        free(conv);
        free(stf);
        free(outf);
    }

    return 0;
}
