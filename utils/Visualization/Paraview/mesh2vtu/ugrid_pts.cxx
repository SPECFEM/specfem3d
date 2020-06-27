//-----------------------------------------------------------------------------
// Program:     ugrid
// Description: Converts x, y, z, s data to VTK XML Unstructured Grid
// File:        ugrid.cxx
// Author:      Nicholas Schwarz, schwarz@evl.uic.edu
//              Electronic Visualization Laboratory
//              University of Illinois at Chicago
// Date:        3 June 2004
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <unistd.h>

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

bool parse_args(int argc, char **argv, char **input, char **output, bool *latlon)
{
  int c;

  *input = *output = NULL;
  *latlon = false;

  while ( (c = getopt(argc, argv, "i:o:L")) != -1) {
    switch (c) {
    case 'i':
      *input = optarg;
      break;
    case 'o':
      *output = optarg;
      break;
    case 'L':
      *latlon = true;
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

void lat_lon_depth_2_xyz(float *x, float *y, float *z)
{
  float lat, lon, depth;
  float theta, phi, r;

  const float R_EARTH_KM = 6371.0;
  const float PI = 3.141592653589793;
  const float D2R = PI/180.0;

  lat = *x; lon = *y; depth = *z;

  theta = (PI / 2.0) - atan(0.99329534 * tan(lat * D2R));
  phi = lon * D2R;

  r = (R_EARTH_KM - depth) / R_EARTH_KM;
  *x = r * sin(theta) * cos(phi);
  *y = r * sin(theta) * sin(phi);
  *z = r * cos(theta);
}

int main(int argc, char** argv) {
  char *input, *output;
  bool latlon;
  float xyz[3];
  float scalar;
  FILE *file;
  int i;
  int npts;

  if (!parse_args(argc, argv, &input, &output, &latlon)) {
    return 1;
  }

  if ((file = fopen(input, "r")) == 0)
    {
      cerr << "ERROR: Can't read file: " << input << "\n";
      return 1;
    }

  fscanf(file, "%d", &npts);

  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();

  for(i = 0; i < npts; i++)
    {
      fscanf(file, "%f", &xyz[0]);
      fscanf(file, "%f", &xyz[1]);
      fscanf(file, "%f", &xyz[2]);
      fscanf(file, "%f", &scalar);
      if (latlon)
        {
          lat_lon_depth_2_xyz(&xyz[0], &xyz[1], &xyz[2]);
        }
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, scalar);
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
