/*!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static FILE *fp = NULL;
static int useBinary = 0;
static int numInColumn = 0;

#define FC_FUNC(name,NAME) name ## _
#define FC_FUNC_(name,NAME) name ## _

static void end_line(void)
{
    if (!useBinary)
    {
        char str2[8] = "\n";
        fprintf(fp, "%s", str2);
        numInColumn = 0;
    }
}

static void open_file(const char *filename)
{
    char full_filename[1024];
    if (strstr(filename, ".vtk") != NULL)
    {
        strcpy(full_filename, filename);
    }
    else
    {
        sprintf(full_filename, "%s.vtk", filename);
    }

    fp = fopen(full_filename, "w+");
}


static void close_file(void)
{
    end_line();
    fclose(fp);
    fp = NULL;
}

static void force_big_endian(unsigned char *bytes)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap & useBinary)
    {
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
}

static void write_string(const char *str)
{
    fprintf(fp, "%s", str);
}


static void new_section(void)
{
    if (numInColumn != 0)
        end_line();
    numInColumn = 0;
}

static void write_int(int val)
{
    if (useBinary)
    {
        force_big_endian((unsigned char *) &val);
        fwrite(&val, sizeof(int), 1, fp);
    }
    else
    {
        char str[128];
        sprintf(str, "%d ", val);
        fprintf(fp, "%s", str);
        if (((numInColumn++) % 9) == 8)
        {
            char str2[8] = "\n";
            fprintf(fp, "%s", str2);
            numInColumn = 0;
        }
    }
}

static void write_float(float val)
{
    if (useBinary)
    {
        force_big_endian((unsigned char *) &val);
        fwrite(&val, sizeof(float), 1, fp);
    }
    else
    {
        char str[128];
        sprintf(str, "%20.12e ", val);
        fprintf(fp, "%s", str);
        if (((numInColumn++) % 9) == 8)
        {
            end_line();
        }
    }
}

static void write_header(void)
{
    fprintf(fp, "# vtk DataFile Version 3.1\n");
    fprintf(fp, "material model VTK file\n");
    if (useBinary)
        fprintf(fp, "BINARY\n");
    else
        fprintf(fp, "ASCII\n");
}

void FC_FUNC_(write_unstructured_mesh,
              WRITE_UNSTRUCTURED_MESH)( char *filename,int * filename_size, int *ub, int *npts,
                                        float *pts, int * ncells, int *celltypes, int *conn,
                                        char * varname,int * varname_size, float *vars)
{
    int   i,j;
    char  str[128];
    int   conn_size = 0;
    int  *curr_conn = conn;

    filename[*filename_size] = '\0';
    varname[*varname_size] = '\0';
    useBinary = *ub;
    open_file(filename);
    write_header();

    write_string("DATASET UNSTRUCTURED_GRID\n");
    sprintf(str, "POINTS %d float\n", *npts);
    write_string(str);
    for (j = 0 ; j < 3*(*npts) ; j++)  write_float(pts[j]);

    new_section();
    conn_size = 9*(*ncells);
    sprintf(str, "CELLS %d %d\n", *ncells, conn_size);
    write_string(str);
    for (i = 0 ; i < *ncells ; i++)
    {
        int npt = 8;
        write_int(npt);
        for (j = 0 ; j < npt ; j++)
            write_int(*curr_conn++);
        end_line();
    }

    new_section();
    sprintf(str, "CELL_TYPES %d\n", *ncells);
    write_string(str);
    for (j = 0 ; j < *ncells ; j++)
    {
        write_int(celltypes[j]);
        end_line();
    }

    new_section();
     sprintf(str, "POINT_DATA %d\n", *npts);
     write_string(str);
     sprintf(str, "SCALARS %s float\n",varname);
            write_string(str);
            write_string("LOOKUP_TABLE default\n");
            for (j = 0 ; j < *npts ; j++) write_float(vars[j]);
            end_line();

    close_file();
}
