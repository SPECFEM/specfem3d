/*
 * asc2sac.c by Eh Tan
 *
 * Copyright (C) 2009, California Institute of Technology.
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
#include <limits.h>
#include <errno.h>
#include "sacio.h"
#include "config.h"


void asc2sac(char *ascfn) {
    sac_int_t npts, max_npts, nerr, nconv, itmp;
    float ftmp;
    char *sacfn;
    float *time, *data, a, b;
    FILE *f;
    
    sacfn = (char *)malloc(strlen(ascfn) + strlen(".sac") + 1);
    if (!sacfn) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    strcpy(sacfn, ascfn);
    strcat(sacfn, ".sac");

    /* reading ascii file */
    f = fopen(ascfn, "r");
    if(f == NULL) {
	fprintf(stderr, "Cannot open file '%s' to read\n", ascfn);
	exit(1);
    }
    
    time = data = 0;
    max_npts = 0;
    for (npts = 0; (nconv = fscanf(f, "%f %f\n", &a, &b)) == 2; ++npts) {
        if (nconv != 2) {
            fprintf(stderr, "error while reading file '%s'\n", ascfn);
            exit(1);
        }
        if (npts >= max_npts) {
            max_npts = max_npts ? 2 * max_npts : 1024;
            time = (float*) realloc(time, max_npts * sizeof(float));
            data = (float*) realloc(data, max_npts * sizeof(float));
            if(time == NULL || data == NULL) {
                fprintf(stderr, "Out of memory\n");
                exit(1);
            }
        }
        time[npts] = a;
        data[npts] = b;
    }
    if (nconv != EOF || ferror(f)) {
        fprintf(stderr, "error while reading file '%s' (on or near line %ld)\n", ascfn, (long)npts);
	exit(1);
    }
    fclose(f);
    /* finished reading ascii file */

    /* write SAC data usng SAC IO library */
    nerr = 0;
    newhdr();
    setnhv("npts", &npts, &nerr, strlen("npts"));
    itmp = 1;
    setlhv("leven", &itmp, &nerr, strlen("leven"));
    ftmp = time[1] - time[0];
    setfhv("delta", &ftmp, &nerr, strlen("delta"));
    setfhv("b", &(time[0]), &nerr, strlen("b"));
    setfhv("e", &(time[npts-1]), &nerr, strlen("e"));
    setihv("iftype", "itime", &nerr, strlen("iftype"), strlen("itime"));
    setihv("idep", "idisp", &nerr, strlen("idep"), strlen("idisp"));

    if(nerr) {
	fprintf(stderr, "error when setting header for '%s'\n", sacfn);
        exit(1);
    }

    wsac0(sacfn, time, data, &nerr, strlen(sacfn));

    if(nerr) {
        fprintf(stderr, "error when writing '%s'\n", sacfn);
        exit(1);
    }
    
    free(time);
    free(data);
    free(sacfn);

    return;
}


int main(int argc, char *argv[]) {
    int i;
     
    if (argc < 2) {
        fprintf(stderr,
                "usage: %s FILE...\n"
                "\n"
                "    Converts ASCII time series files to SAC binary files.\n"
                "    Input files must have two columns:  t, f(t).\n"
                "\n",
                argv[0]);
        exit(1);
    }
    
    for (i = 1; i < argc; ++i) {
        asc2sac(argv[i]);
    }
    
    return 0;
}
