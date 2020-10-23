//-----------------------------------------------------------------------------
// Program:     surf2vtu
// Description: Converts x, y, z, s points to VTK XML Unstructured Grid
//              using definition of cells
// File:        surf2vtu.cxx
// Author:      Nicholas Schwarz, schwarz@evl.uic.edu
//              Electronic Visualization Laboratory
//              University of Illinois at Chicago
//               - wrote ugrid and introduced vtk to Seismo Lab, Caltech
//
//              Brian Savage savage13@gps.caltech.edu
//              California Institute of Technology
//              Geologial and Planetary Sciences
//
//              Qinya Liu, modification from mesh2vtu
//
// Input:       in binary
//              integer    number of points
//              4 floats   (x,y,z,s) point 0
//                ...
//              4 floats   (x,y,z,s) point n-1
//              integer    number of quad cells
//              4 integers (1-4) cell 0
//                ...      define a hexahedron of 8 points
//              4 integers (1-4) cell n-1
//
// Date:        4  June 2004 ver 1.0 (was ugrid)
//                 - original version, only read in x,y,z,s points
//              25 June 2004 ver 2.0 (mesh2vtu)
//                 - reads in cell definition
//                 - input is done in binary
//              29 Sep 2005 modified to quad2vtu
//                 - to deal with quad elements instead of Hexahedron
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
#include <vtkQuad.h>

void usage(char *progname)
{
  printf("Usage: %s -i input-file -o output-file\n"
         "    Takes an input file (binary) with a number of points and a number of cells\n"
         "    and transforms them into an unstructured grid file\n"
         "\n"
         "    -i input-file (Binary file)\n"
         "    -o output-file (XML Unstructured Grid File)\n"
         "\n"
         "    Input Binary files have this structure:\n"
         "      number_of_points          integer (4 bytes)\n"
         "      x_1, y_1, z_1, scalar_1   4 reals (4 bytes each)\n"
         "      ...\n"
         "      x_n, y_n, z_n, scalar_n   4 reals (4 bytes each)\n"
         "      number_of_cells           integer (4 bytes)\n"
         "      cell_1 (four points)      4 integers (4 bytes each)\n"
         "      ...\n"
         "      cell_n                    4 integers (4 bytes each)\n"
         "\n", progname);
}

bool parse_args(int argc, char **argv, char **input, char **output)
{
  int c;

  *input = *output = NULL;

  while ( (c = getopt(argc, argv, "i:o:")) != -1) {
    switch (c) {
    case 'i':
      *input = optarg;
      break;
    case 'o':
      *output = optarg;
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

int main(int argc, char** argv) {
  char *input, *output;
  float xyz[3];
  float scalar;
  int cell[8];
  int i, j;
  int npts, ncells;
  int fd;

  if (!parse_args(argc, argv, &input, &output)) {
    return 1;
  }

  if ((fd = open(input, O_RDONLY)) == -1) {
    printf("Error opening file: %s.\n", input);
    return 1;
  }

  if(read(fd, &npts, sizeof(int)) != sizeof(int)) {
    printf("Bad read on file (in points): %s\n", input);
    return 1;
  }

  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();
  printf("surf2vtu: Reading in points: %d\n", npts);
  for (i = 0 ; i < npts ; i++)
    {
      read(fd, &xyz[0], sizeof(float));
      read(fd, &xyz[1], sizeof(float));
      read(fd, &xyz[2], sizeof(float));
      read(fd, &scalar, sizeof(float));
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, scalar);
    }

  vtkCellArray *cells = vtkCellArray::New();
  if(read(fd, &ncells, sizeof(int)) != sizeof(int)) {
    printf("Bad read on file (in cells): %s\n", input);
    return 1;
  }
  printf("surf2vtu: Reading in cells: %d\n", ncells);
  int *cellTypes = new int[ncells];
  vtkQuad *quad = vtkQuad::New();
  quad->GetPointIds()->SetNumberOfIds(4);

  for(i = 0; i < ncells; i++) {
    for(j = 0; j < 4; j++) {
      read(fd, &cell[j], sizeof(int));
      quad->GetPointIds()->SetId(j,cell[j]);
    }
    cells->InsertNextCell(quad);
    cellTypes[i] = quad->GetCellType();
  }

  close(fd);

  dataSet -> SetPoints(newPts);
  dataSet -> GetPointData() -> SetScalars(newScalars);
  dataSet -> SetCells(cellTypes, cells);

  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION <= 5
  writer -> SetInput(dataSet);
#else
  writer -> SetInputData(dataSet);
#endif
  writer -> SetFileName(output);
  writer -> Write();

  writer -> Delete();
  newPts -> Delete();
  newScalars -> Delete();
  dataSet -> Delete();
  cells -> Delete();

  printf("Done.\n");

  return 0;

}
