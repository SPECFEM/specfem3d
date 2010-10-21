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
#include <vtk/vtkHexahedron.h>


int main(int argc, char** argv) {

  if (argc < 3) {
    printf("Usage: ugrid input_file output_file\n");
    return 0;
  }

  float xyz[3];
  float scalar;
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
  
  printf("mesh2vtu: Reading in points: %d\n", npts);

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();
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
    printf("Bad read on file (in cells): %s\n", argv[1]);
  }
  
  printf("mesh2vtu: Reading in cells: %d\n", ncells);  

  int *cellTypes = new int[ncells];
  vtkHexahedron *hex = vtkHexahedron::New();
  hex->GetPointIds()->SetNumberOfIds(8);

  for(i = 0; i < ncells; i++) {
    for(j = 0; j < 8; j++) {
      read(fd, &cell[j], sizeof(int));
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
  writer -> SetInput(dataSet);
  writer -> SetFileName(argv[2]);
  writer -> Write();

  writer -> Delete();
  
  printf("Done.\n");
 
  return 0;

}
