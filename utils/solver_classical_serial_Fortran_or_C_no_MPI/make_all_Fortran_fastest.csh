
rm -f xspecfem3D* *.o

# Intel ifort compiler
# add -Winline to get information about routines that are inlined
# add -vec-report3 to get information about loops that are vectorized or not
# do NOT suppress -ftz, which is critical for performance
ifort -O3 -xSSE4.2 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds -align sequence -assume byterecl -ftz -o xspecfem3D_F90 serial_specfem3D_inlined_v03_is_the_fastest_no_more_function_calls.f90 read_arrays_solver.f90

# other compilers

#gfortran -std=gnu -fimplicit-none -O3 -fno-trapping-math -Wunused-labels -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow -o xspecfem3D_F90 serial_specfem3D_inlined_v03_is_the_fastest_no_more_function_calls.f90 read_arrays_solver.f90

#pgf90 -fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -Mstandard -fastsse -tp amd64e -o xspecfem3D_F90 serial_specfem3D_inlined_v03_is_the_fastest_no_more_function_calls.f90 read_arrays_solver.f90

#xlf_r -O3 -qsave -qstrict -qtune=ppc970 -qarch=ppc64v -qcache=auto -qfree=f90 -Q -qflttrap=en:ov:zero:inv -o xspecfem3D_f90 serial_specfem3D_inlined_v03_is_the_fastest_no_more_function_calls.f90 read_arrays_solver.f90

