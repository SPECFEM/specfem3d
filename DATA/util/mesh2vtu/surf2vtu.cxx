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
//              integer    number of quad cells
//              4 integers (1-4) cell 0
//                ...      define a hexahedron of 8 points
//              4 integers (1-8) cell n-1
//              
// Date:        4  June 2004 ver 1.0 (was ugrid)
//                 - original version, only read in x,y,z,s points
//              25 June 2004 ver 2.0 (mesh2vtu)
//                 - reads in cell definition
//                 - input is done in binary
//              29 Sep 2005 modified to quad2vtu
//                 - to deal with quad elements instead of Hexhedron
// 
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <vtk/vtkFloatArray.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>
#include <vtk/vtkUnstructuredGridToPolyDataFilter.h>
#include <vtk/vtkXMLPolyDataWriter.h>
#include <vtk/vtkUnstructuredGridToPolyDataFilter.h>
#include <vtk/vtkDelaunay3D.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkPointSet.h>
#include <vtk/vtkQuad.h>


int main(int argc, char** argv) {

  if (argc < 3) {
    printf("Usage: ugrid input_file output_file\n");
    return 0;
  }

  float xyz[3];
  int cell[8];
  FILE *file;
  int i, j;
  int npts, ncells;
  int pid[8];

  int fd;
  
  if((fd = open(argv[1], O_RDONLY)) == -1) {
    printf("Error opening file: %s.\n", argv[1]);
    return 0;
  }

  if(read(fd, &npts, sizeof(int)) != sizeof(int)) {
    printf("Bad read on file (in points): %s\n", argv[1]);
  }
  
  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
  float *xV = new float[npts];
  float *yV = new float[npts];
  float *zV = new float[npts];
  float *sV = new float[npts];

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();
  printf("mesh2vtu: Reading in points: %d\n", npts);
  for (i = 0 ; i < npts ; i++)
    {
      read(fd, &xV[i], sizeof(float));
      read(fd, &yV[i], sizeof(float));
      read(fd, &zV[i], sizeof(float));
      read(fd, &sV[i], sizeof(float));
      xyz[0] = xV[i]; 
      xyz[1] = yV[i]; 
      xyz[2] = zV[i];
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, sV[i]);
    }

  vtkCellArray *cells = vtkCellArray::New();
  if(read(fd, &ncells, sizeof(int)) != sizeof(int)) {
    printf("Bad read on file (in cells): %s\n", argv[1]);
  }
  printf("mesh2vtu: Reading in cells: %d\n", ncells);  
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
  writer -> SetInput(dataSet);
  writer -> SetFileName(argv[2]);
  writer -> Write();

  writer -> Delete();
  newPts -> Delete();
  newScalars -> Delete();
  dataSet -> Delete();
  cells -> Delete();
  
  //  printf("Done.\n");
 
  return 0;

}
