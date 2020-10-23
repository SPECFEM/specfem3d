# Makefile

#############################################################
## modify to match your compiler defaults
CXX=g++

## modify to match your library paths
# VTK libraries
VTK := /usr

EXISTS := $(shell cmake --find-package -DNAME:STRING=VTK -DCOMPILER_ID:STRING=$(CXX) -DLANGUAGE:STRING=CXX -DMODE:STRING=EXIST &> /dev/null ; echo $$?)
ifeq ($(EXISTS),0)
VTK_CPPFLAGS := $(shell cmake --find-package -DNAME:STRING=VTK -DCOMPILER_ID:STRING=$(CXX) -DLANGUAGE:STRING=CXX -DMODE:STRING=COMPILE)
VTK_LDFLAGS := $(shell cmake --find-package -DNAME:STRING=VTK -DCOMPILER_ID:STRING=$(CXX) -DLANGUAGE:STRING=CXX -DMODE:STRING=LINK)
else
$(warning CMake not found or could not find VTK development headers. Using $$(VTK)=$(VTK) instead.)
VTK_CPPFLAGS := -I$(VTK)/include -I$(VTK)/include/vtk
VTK_LDFLAGS := -L$(VTK)/lib64/vtk -L$(VTK)/lib/vtk -lvtkCommon -lvtkIO -lvtkFiltering -lvtkGraphics
endif

CPPFLAGS = -O3 $(VTK_CPPFLAGS) -Wall -Wno-deprecated
LIBS = $(VTK_LDFLAGS)

#############################################################

TARGETS = mesh2vtu cell2vtu surf2vtu ugrid ugrid_pts

all: $(TARGETS)

%: %.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LIBS)

%.o: %.cxx
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f $(TARGETS) $(addsuffix .o, $(TARGETS))
.PHONY: clean

