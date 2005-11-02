
CXX=g++
VTK=/opt/seismo-util/source/vtk # please change according to the location of vtk in your system
INCLUDES=-I$(VTK)/include
CPPFLAGS=$(INCLUDES) -Wno-deprecated
LIBS=-L$(VTK)/lib/vtk -lvtkCommon -lvtkIO -lvtkRendering -lvtkFiltering -lvtkGraphics

all: mesh2vtu ugrid ugrid_pts cell2vtu  surf2vtu

mesh2vtu: mesh2vtu.o
	$(CXX) $(CPPFLAGS) -o mesh2vtu mesh2vtu.o $(LIBS)

mesh2vtu.o: mesh2vtu.cxx
	$(CXX) $(CPPFLAGS) -c mesh2vtu.cxx -o mesh2vtu.o

surf2vtu: surf2vtu.o
	$(CXX) $(CPPFLAGS) -o surf2vtu surf2vtu.o $(LIBS)

surf2vtu.o: surf2vtu.cxx
	$(CXX) $(CPPFLAGS) -c surf2vtu.cxx -o surf2vtu.o

cell2vtu: cell2vtu.o
	$(CXX) $(CPPFLAGS) -o cell2vtu cell2vtu.o $(LIBS)

cell2vtu.o: cell2vtu.cxx
	$(CXX) $(CPPFLAGS) -c cell2vtu.cxx -o cell2vtu.o

ugrid: ugrid.o
	$(CXX) $(CPPFLAGS) -o ugrid ugrid.o $(LIBS)

ugrid.o: ugrid.cxx
	$(CXX) $(CPPFLAGS) -c ugrid.cxx -o ugrid.o

ugrid_pts: ugrid_pts.o
	$(CXX) $(CPPFLAGS) -o ugrid_pts ugrid_pts.o $(LIBS)

ugrid_pts.o: ugrid_pts.cxx
	$(CXX) $(CPPFLAGS) -c ugrid_pts.cxx -o ugrid_pts.o

clean:
	rm -f mesh2vtu ugrid ugrid_pts *.o *~ core surf2vtu cell2vtu
