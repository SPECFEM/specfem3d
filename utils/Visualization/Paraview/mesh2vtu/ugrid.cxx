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
#include <vtk/vtkHexahedron.h>
#include <vtk/vtkCellArray.h>
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
  int npts, ncells;
  int pid[8];
  fprintf(stderr,"Welcome to ugrid\n");
  if ((file = fopen(argv[1], "r")) == 0)
    {
      //cerr << "ERROR: Can't read file: " << filename << "\n";
      return 0;
    }


  // get points

  fscanf(file, "%d", &npts);
  fprintf(stderr, "npts: %d\n", npts);
  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
  float *xV = new float[npts];
  float *yV = new float[npts];
  float *zV = new float[npts];
  float *sV = new float[npts];

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();

  for (i = 0 ; i < npts ; i++)
    {
      fscanf(file, "%f", &xV[i]);
      fscanf(file, "%f", &yV[i]);
      fscanf(file, "%f", &zV[i]);
      fscanf(file, "%f", &sV[i]);
      fprintf(stderr, "%f %f %f %f\n", xV[i], yV[i], zV[i], sV[i]);
      xyz[0] = xV[i]; 
      xyz[1] = yV[i]; 
      xyz[2] = zV[i];
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, sV[i]);
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
  writer -> SetInput(dataSet);
  writer -> SetFileName(argv[2]);
  writer -> Write();

  writer -> Delete();
  newPts -> Delete();
  newScalars -> Delete();
  dataSet -> Delete();
 
  return 0;

}
