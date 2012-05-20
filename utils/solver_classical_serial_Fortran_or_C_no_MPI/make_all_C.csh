
rm -f xspecfem3D* *.o

# add -vec-report3 to get information about loops that are vectorized or not
# do NOT suppress -ftz, which is critical for performance

# on some machines -fast does NOT work when linking C with Fortran for some reason
# in that case you can switch back to -O3 -xSSE4.1 .
# Can add -ftrapuv -traceback to debug if needed.
#icc -c -O3 -xSSE4.1 -ftz -funroll-loops -unroll5 -vec-report1 -std=c99 -x c -Wcheck serial_specfem3D_single_no_Deville.c
#ifort -o xspecfem3D_C -O3 -xSSE4.1 -ftz -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds serial_specfem3D_single_no_Deville.o read_arrays_solver.f90 -nofor_main

icc -c -O3 -xSSE4.1 -ftz -funroll-loops -unroll5 -vec-report1 -std=c99 -x c -Wcheck serial_specfem3D_single_with_Deville.c
ifort -o xspecfem3D_C -O3 -xSSE4.1 -ftz -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds serial_specfem3D_single_with_Deville.o read_arrays_solver.f90 -nofor_main

# GNU gcc (but I have found NO way of really turning flush-to-zero (FTZ) on i.e. turning
# gradual underflow off, and as a result the code is twice slower than it should when using gcc!!!)
#
# gcc -c -fbounds-check -Wall -fno-trapping-math -fno-signaling-nans -std=gnu99 -O3 serial_specfem3D_single_with_Deville.c
# gfortran -o xspecfem3D_C -O3 -std=f95 -fimplicit-none -frange-check -O3 -pedantic -pedantic-errors -Waliasing -Wampersand -Wline-truncation -Wsurprising -Wunderflow -fno-trapping-math serial_specfem3D_single_with_Deville.o read_arrays_solver.f90

# g++ -fno-trapping-math -O3 -o xspecfem3D_C serial_specfem3D_single_with_Deville.c -lm 

# pgcc -fast -Mnobounds -Minline -Mneginfo -Knoieee -Minform=warn -fastsse -tp amd64e -o xspecfem3D_C serial_specfem3D_single_with_Deville.c -lm

rm -f *.o

