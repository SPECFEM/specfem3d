
rm -f xspecfem3D* xsmooth_a_gradient *.o

############
############  normal (i.e. simpler) version
############

# Intel ifort compiler
# add -Winline to get information about routines that are inlined
# add -vec-report3 to get information about loops that are vectorized or not
# do NOT suppress -ftz, which is critical for performance
#ifort -O3 -xHost -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds -align sequence -assume byterecl -ftz -o xspecfem3D_F90_normal specfem3D_normal_no_Deville.f90 read_arrays_solver.f90

# other compilers

gfortran -std=gnu -fimplicit-none -O3 -fno-trapping-math -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow -o xspecfem3D_F90_normal specfem3D_normal_no_Deville.f90 read_arrays_solver.f90

#pgf90 -fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -Mstandard -fastsse -tp amd64e -o xspecfem3D_F90_normal specfem3D_normal_no_Deville.f90 read_arrays_solver.f90

#xlf_r -O3 -qsave -qstrict -qtune=ppc970 -qarch=ppc64v -qcache=auto -qfree=f90 -Q -qflttrap=en:ov:zero:inv -o xspecfem3D_F90_normal specfem3D_normal_no_Deville.f90 read_arrays_solver.f90

############
############  faster (i.e. optimized) version
############

# Intel ifort compiler
# add -Winline to get information about routines that are inlined
# add -vec-report3 to get information about loops that are vectorized or not
# do NOT suppress -ftz, which is critical for performance
#ifort -O3 -xHost -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds -align sequence -assume byterecl -ftz -o xspecfem3D_F90_faster specfem3D_fastest_version_with_Deville_and_inlining.f90 read_arrays_solver.f90

# other compilers

gfortran -std=gnu -fimplicit-none -O3 -fno-trapping-math -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow -o xspecfem3D_F90_faster specfem3D_fastest_version_with_Deville_and_inlining.f90 read_arrays_solver.f90

#pgf90 -fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -Mstandard -fastsse -tp amd64e -o xspecfem3D_F90_faster specfem3D_fastest_version_with_Deville_and_inlining.f90 read_arrays_solver.f90

#xlf_r -O3 -qsave -qstrict -qtune=ppc970 -qarch=ppc64v -qcache=auto -qfree=f90 -Q -qflttrap=en:ov:zero:inv -o xspecfem3D_F90_faster specfem3D_fastest_version_with_Deville_and_inlining.f90 read_arrays_solver.f90

############
############  demo code to show how to smooth a gradient
############

gfortran -std=gnu -fimplicit-none -O3 -fno-trapping-math -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow -o xsmooth_a_gradient specfem3D_how_to_consistently_average_a_gradient.f90 read_arrays_solver.f90

