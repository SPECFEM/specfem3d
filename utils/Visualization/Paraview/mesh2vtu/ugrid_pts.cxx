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
#include <vtk/vtkFloatArray.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("Usage: ugrid input_file output_file\n");
    return 0;
  }

  float xyz[3];
  FILE *file;
  int i;
  int npts;
  
  if ((file = fopen(argv[1], "r")) == 0)
    {
      //cerr << "ERROR: Can't read file: " << filename << "\n";
      return 0;
    }
   
  fscanf(file, "%d", &npts);

  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
  float *xV = new float[npts];
  float *yV = new float[npts];
  float *zV = new float[npts];
  float *sV = new float[npts];

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();

  for(i = 0; i < npts; i++)
    {
      fscanf(file, "%f", &xV[i]);
      fscanf(file, "%f", &yV[i]);
      fscanf(file, "%f", &zV[i]);
      fscanf(file, "%f", &sV[i]);
      xyz[0] = xV[i]; 
      xyz[1] = yV[i]; 
      xyz[2] = zV[i];
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, sV[i]);
    }

  dataSet -> SetPoints(newPts);
  dataSet -> GetPointData() -> SetScalars(newScalars);

  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  writer -> SetInput(dataSet);
  writer -> SetFileName(argv[2]);
  writer -> Write();

  writer -> Delete();
  newPts -> Delete();
  newScalars -> Delete();
  dataSet -> Delete();
 
  return 0;

}
