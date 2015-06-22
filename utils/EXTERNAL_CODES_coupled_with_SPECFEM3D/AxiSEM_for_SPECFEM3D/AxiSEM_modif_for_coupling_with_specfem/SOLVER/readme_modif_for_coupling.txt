list of files modified for the coupling AxiSEM with an other simulation code
(eg FD, DG , SEM):

main.f90
data_io.f90
seismograms.f90
def_precomp_terms.f90
parameters.F90
time_evol_wave.F90

file added :
coupling_mod.f90


  short description of modification:

main.f90:
 add read_boundary_coordinates calling routine and finalize_coupling. In order
 to read and write parameters for the coupling stuff.

 data_io.f90 :
   add logical flag : coupling

 seismograms.f90 :
   add test when define dumping files, in the coupluing case we use an other way to
   dump outputs files in disk.

 def_precomp_terms.f90 :
   add arrays to store the physical properties of the model. Thoses arrays
   will use when commuting stress from strain.

parameters.F90 :
  hardcoded dump_type='coupling' and coupling=.true. to  be changed in next
  version

time_evol_wave.F90 :

 added a subroutine call for writting the solution in the edge of regional box
 boundary.


  All the new routines needed for the coupling are in coupling_mod.f90 that is
  a new file in AxiSEM code

TO DO : for now tthe choice to run axisem in coupling mode is hardcoded, we
must change it and leave the choice in input parameter.

