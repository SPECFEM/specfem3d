# Makefile

#############################################################
## modify to match your compiler defaults 
CXX=g++

## modify to match your library paths
# VTK libraries
VTK=/opt/seismo-util
LIBS=-L$(VTK)/lib/vtk -lvtkCommon -lvtkIO -lvtkRendering -lvtkFiltering -lvtkGraphics

INCLUDES=-I$(VTK)/include -I$(VTK)/include/vtk
CPPFLAGS= -O3 $(INCLUDES) -Wno-deprecated

#############################################################


all: mesh2vtu ugrid ugrid_pts

mesh2vtu: mesh2vtu.o
	$(CXX) $(CPPFLAGS) -o mesh2vtu mesh2vtu.o $(LIBS)

mesh2vtu.o: mesh2vtu.cxx
	$(CXX) $(CPPFLAGS) -c mesh2vtu.cxx -o mesh2vtu.o

ugrid: ugrid.o
	$(CXX) $(CPPFLAGS) -o ugrid ugrid.o $(LIBS)

ugrid.o: ugrid.cxx
	$(CXX) $(CPPFLAGS) -c ugrid.cxx -o ugrid.o

ugrid_pts: ugrid_pts.o
	$(CXX) $(CPPFLAGS) -o ugrid_pts ugrid_pts.o $(LIBS)

ugrid_pts.o: ugrid_pts.cxx
	$(CXX) $(CPPFLAGS) -c ugrid_pts.cxx -o ugrid_pts.o

clean:
	rm -f mesh2vtu ugrid ugrid_pts *.o *~ core
