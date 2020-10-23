//-----------------------------------------------------------------------------
// Program:     mesh2vtu
// Description: Converts x, y, z, s points to VTK XML Unstructured Grid
//              using definition of cells
// File:        mesh2vtu.cxx
// Author:      Nicholas Schwarz, schwarz@evl.uic.edu
//              Electronic Visualization Laboratory
//              University of Illinois at Chicago
//               - wrote ugrid and introduced vtk to Seismo Lab, Caltech
//
//              Brian Savage savage13@gps.caltech.edu
//              California Institute of Technology
//              Geologial and Planetary Sciences
//
// Input:       in binary
//              integer    number of points
//              4 floats   (x,y,z,s) point 0
//                ...
//              4 floats   (x,y,z,s) point n-1
//              integer    number of cells
//              8 integers (1-8) cell 0
//                ...      define a hexahedron of 8 points
//              8 integers (1-8) cell n-1
//
// Date:        4  June 2004 ver 1.0 (was ugrid)
//                 - original version, only read in x,y,z,s points
//              25 June 2004 ver 2.0 (mesh2vtu)
//                 - reads in cell definition
//                 - input is done in binary
//              26 February 2010
//                 - changes array allocation for reading in points
//
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#if VTK_MAJOR_VERSION <= 5
#include <vtkUnstructuredGridToPolyDataFilter.h>
#endif
#include <vtkXMLPolyDataWriter.h>
#include <vtkDelaunay3D.h>
#include <vtkCellArray.h>
#include <vtkPointSet.h>
#include <vtkHexahedron.h>

void usage(char *progname)
{
  printf("Usage: %s -i input-file -o output-file\n"
         "    Takes an input file (binary) with a number of points and a number of cells\n"
         "    and transforms them into an unstructured grid file\n"
         "\n"
         "    -i input-file (Binary file)\n"
         "    -o output-file (XML Unstructured Grid File)\n"
         "    -s Perform byte swapping on input file\n"
         "\n"
         "    Input Binary files have this structure:\n"
         "      number_of_points          integer (4 bytes)\n"
         "      x_1, y_1, z_1, scalar_1   4 reals (4 bytes each)\n"
         "      ...\n"
         "      x_n, y_n, z_n, scalar_n   4 reals (4 bytes each)\n"
         "      number_of_cells           integer (4 bytes)\n"
         "      cell_1 (eight points)     8 integers (4 bytes each)\n"
         "      ...\n"
         "      cell_n                    8 integers (4 bytes each)\n"
         "\n", progname);
}

bool parse_args(int argc, char **argv, char **input, char **output, bool *swap)
{
  int c;

  *input = *output = NULL;
  *swap = false;

  while ( (c = getopt(argc, argv, "i:o:s")) != -1) {
    switch (c) {
    case 'i':
      *input = optarg;
      break;
    case 'o':
      *output = optarg;
      break;
    case 's':
      *swap = true;
      break;
    case '?':
      usage(argv[0]);
      return false;
    default:
      printf("?? getopt returned character code 0%o ??\n", c);
      return false;
    }
  }

  if (*input == NULL) {
    printf("ERROR: Must specify input file -i input-file\n\n");
    usage(argv[0]);
    return false;
  }

  if (*output == NULL) {
    printf("ERROR: Must specify output file -o output-file\n\n");
    usage(argv[0]);
    return false;
  }

  return true;
}

bool read_int32_normal(int fd, int *val)
{
  return read(fd, val, sizeof(*val)) == sizeof(*val);
}

bool read_int32_swap(int fd, int *val)
{
  if (!read_int32_normal(fd, val))
    return false;
  *val = (*val<<24) | (*val<<8 & 0xff0000) | (*val>>8 & 0xff00) | (*val>>24);
  return true;
}

bool read_float32_normal(int fd, float *val)
{
  return read(fd, val, sizeof(*val)) == sizeof(*val);
}

bool read_float32_swap(int fd, float *val)
{
  int tmp;
  if (!read_int32_normal(fd, &tmp))
    return false;
  tmp = (tmp<<24) | (tmp<<8 & 0xff0000) | (tmp>>8 & 0xff00) | (tmp>>24);
  memcpy(val, &tmp, sizeof(tmp));
  return true;
}

int main(int argc, char** argv) {
  char *input, *output;
  bool swap;
  float xyz[3];
  float scalar;
  int cell[8];
  int i, j;
  int npts, ncells;
  int fd;
  bool (*read_int32)(int fd, int *val);
  bool (*read_float32)(int fd, float *val);

  if (!parse_args(argc, argv, &input, &output, &swap)) {
    return 1;
  }

  if (swap) {
    read_int32 = &read_int32_swap;
    read_float32 = &read_float32_swap;
  } else {
    read_int32 = &read_int32_normal;
    read_float32 = &read_float32_normal;
  }

  if ((fd = open(input, O_RDONLY)) == -1) {
    printf("Error opening file: %s.\n", input);
    return 1;
  }

  if (!read_int32(fd, &npts)) {
    printf("Bad read on file (in points): %s\n", input);
    return 1;
  }

  printf("mesh2vtu: Reading in points: %d\n", npts);

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();
  for (i = 0 ; i < npts ; i++)
  {
    read_float32(fd, &xyz[0]);
    read_float32(fd, &xyz[1]);
    read_float32(fd, &xyz[2]);
    read_float32(fd, &scalar);

    newPts -> InsertPoint(i, xyz);
    newScalars -> InsertValue(i, scalar);
  }

  vtkCellArray *cells = vtkCellArray::New();
  if (!read_int32(fd, &ncells)) {
    printf("Bad read on file (in cells): %s\n", input);
    return 1;
  }

  printf("mesh2vtu: Reading in cells: %d\n", ncells);

  int *cellTypes = new int[ncells];
  vtkHexahedron *hex = vtkHexahedron::New();
  hex->GetPointIds()->SetNumberOfIds(8);

  for(i = 0; i < ncells; i++) {
    for(j = 0; j < 8; j++) {
      read_int32(fd, &cell[j]);
      hex->GetPointIds()->SetId(j,cell[j]);
    }
    cells->InsertNextCell(hex);
    cellTypes[i] = hex->GetCellType();
  }

  close(fd);

  printf("Creating unstructured grid...\n");

  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
  dataSet -> SetPoints(newPts);
  dataSet -> GetPointData() -> SetScalars(newScalars);
  dataSet -> SetCells(cellTypes, cells);

  cells-> Delete();
  newPts -> Delete();
  newScalars -> Delete();

  printf("Writing out VTU file...\n");

  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION <= 5
  writer -> SetInput(dataSet);
#else
  writer -> SetInputData(dataSet);
#endif
  writer -> SetFileName(output);
  writer -> Write();

  writer -> Delete();

  printf("Done.\n");

  return 0;

}
