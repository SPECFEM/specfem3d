//-----------------------------------------------------------------------------
// Program:     ugrid
// Description: Converts x, y, z, s data to VTK XML Unstructured Grid
// File:        ugrid.cxx
// Author:      Nicholas Schwarz, schwarz@evl.uic.edu
//              Electronic Visualization Laboratory
//              University of Illinois at Chicago
// Date:        4 June 2004
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <unistd.h>

#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

void usage(char *progname)
{
  printf("Usage: %s -i input-file -o output-file\n"
         "    Takes an input file (ascii) with a number of points\n"
         "    and transforms them into an unstructured grid file\n"
         "\n"
         "    -i input-file (ascii file)\n"
         "    -o output-file (XML Unstructured Grid File)\n"
         "    -L - input is in lon lat depth scalar\n"
         "\n"
         "    Input ascii files have this structure:\n"
         "      number_of_points\n"
         "      x_1 y_1 z_1 scalar_1\n"
         "      ...\n"
         "      x_n y_n z_n scalar_n\n"
         "\n", progname);
}

bool parse_args(int argc, char **argv, char **input, char **output)
{
  int c;

  *input = *output = NULL;

  while ( (c = getopt(argc, argv, "i:o:L")) != -1) {
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
  FILE *file;
  int i;
  int npts, ncells;
  int pid[8];

  if (!parse_args(argc, argv, &input, &output)) {
    return 1;
  }

  fprintf(stderr,"Welcome to ugrid\n");
  if ((file = fopen(input, "r")) == 0)
    {
      cerr << "ERROR: Can't read file: " << input << "\n";
      return 1;
    }

  // get points

  fscanf(file, "%d", &npts);
  fprintf(stderr, "npts: %d\n", npts);
  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();

  for (i = 0 ; i < npts ; i++)
    {
      fscanf(file, "%f", &xyz[0]);
      fscanf(file, "%f", &xyz[1]);
      fscanf(file, "%f", &xyz[2]);
      fscanf(file, "%f", &scalar);
      fprintf(stderr, "%f %f %f %f\n", xyz[0], xyz[1], xyz[2], scalar);
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, scalar);
    }

  // get cells
  fscanf(file, "%d", &ncells);

  vtkHexahedron* ahex;

  for (i = 0 ; i < ncells ; i++)
    {
      fscanf(file, "%d %d %d %d %d %d %d %d", 
	     &pid[0], &pid[1], &pid[2], &pid[3], 
	     &pid[4], &pid[5], &pid[6], &pid[7]);
      ahex = vtkHexahedron::New();
      ahex -> GetPointIds() -> InsertNextId(pid[0]);
      ahex -> GetPointIds() -> InsertNextId(pid[1]);
      ahex -> GetPointIds() -> InsertNextId(pid[2]);
      ahex -> GetPointIds() -> InsertNextId(pid[3]);
      ahex -> GetPointIds() -> InsertNextId(pid[4]);
      ahex -> GetPointIds() -> InsertNextId(pid[5]);
      ahex -> GetPointIds() -> InsertNextId(pid[6]);
      ahex -> GetPointIds() -> InsertNextId(pid[7]);
      dataSet -> InsertNextCell(ahex -> GetCellType(), ahex -> GetPointIds());
    }

  dataSet -> SetPoints(newPts);
  dataSet -> GetPointData() -> SetScalars(newScalars);

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

  return 0;

}
