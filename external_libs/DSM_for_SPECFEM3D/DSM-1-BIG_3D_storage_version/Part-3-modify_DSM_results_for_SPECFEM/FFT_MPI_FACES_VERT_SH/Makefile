BIN=../../bin

include ../../config.h
program = $(BIN)/fft_face_vert_SH
objs =    rotations_matrix.o others.o TraFFT.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(objs) *.lst *.o

all: $(program)

