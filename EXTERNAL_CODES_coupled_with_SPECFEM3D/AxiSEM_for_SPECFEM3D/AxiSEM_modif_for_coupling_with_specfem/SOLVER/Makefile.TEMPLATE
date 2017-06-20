PROG = axisem

include ../make_axisem.macros

SRCS =	background_models.F90 commpi.F90 commun.F90 get_model.F90 \
	lateral_heterogeneities.F90 meshes_io.F90 nc_helpers.F90 \
	nc_routines.F90 nc_snapshots.F90 parameters.F90 time_evol_wave.F90 \
	analytic_mapping.f90 analytic_semi_mapping.f90 \
	analytic_spheroid_mapping.f90 apply_masks.f90 attenuation.f90 \
	clocks.f90 clocks_wrapper_solver.f90 data_comm.f90 \
	data_heterogeneous.f90 data_io.f90 data_matr.f90 data_mesh.f90 \
	data_pointwise.f90 data_proc.f90 data_source.f90 data_spec.f90 \
	data_time.f90 def_grid.f90 def_precomp_terms.f90 get_mesh.f90 \
	global_parameters.f90 interpolation.f90 kdtree2.f90 list.f90 main.f90 \
	pointwise_derivatives.f90 rotations.f90 seismograms.f90 source.f90 \
	stiffness_di.f90 stiffness_fluid.f90 stiffness_mono.f90 \
	stiffness_quad.f90 subpar_mapping.f90 unrolled_loops.f90 utlity.f90 \
	coupling_mod.f90 wavefields_io.f90 ftz.c pthread.c

OBJS =	background_models.o commpi.o commun.o get_model.o \
	lateral_heterogeneities.o meshes_io.o nc_helpers.o nc_routines.o \
	nc_snapshots.o parameters.o time_evol_wave.o analytic_mapping.o \
	analytic_semi_mapping.o analytic_spheroid_mapping.o apply_masks.o \
	attenuation.o clocks.o clocks_wrapper_solver.o data_comm.o \
	data_heterogeneous.o data_io.o data_matr.o data_mesh.o \
	data_pointwise.o data_proc.o data_source.o data_spec.o data_time.o \
	def_grid.o def_precomp_terms.o get_mesh.o global_parameters.o \
	interpolation.o kdtree2.o list.o main.o pointwise_derivatives.o \
	rotations.o seismograms.o source.o stiffness_di.o stiffness_fluid.o \
	stiffness_mono.o stiffness_quad.o subpar_mapping.o unrolled_loops.o \
	utlity.o coupling_mod.o wavefields_io.o ftz.o pthread.o

ifeq ($(strip $(USE_NETCDF)),true)
   FFLAGS += -Denable_netcdf
   ifeq ($(strip $(USE_PAR_NETCDF)),true)
   	   FFLAGS += -Denable_parallel_netcdf
   endif

   ifdef NETCDF_PATH
       LIBS = -L $(strip $(NETCDF_PATH))/lib -lnetcdff -Wl,-rpath,$(strip $(NETCDF_PATH))/lib
       INCLUDE = -I $(strip $(NETCDF_PATH))/include
   else
       LIBS = -lnetcdff
       INCLUDE = -I /usr/include
   endif
else
   LIBS = 
   INCLUDE = 
endif

ifeq ($(strip $(SERIAL)),true)
   FFLAGS += -Dserial
   LDFLAGS += -pthread
endif
ifeq ($(strip $(INCLUDE_MPI)),true)
    FFLAGS += -Dinclude_mpi
endif


# cancel m2c implicit rule 
%.o : %.mod 
 

# SOLVER specific code
FFLAGS += -Dsolver
 
all: $(PROG) utils 

$(PROG): $(OBJS)
    ifeq ($(strip $(SERIAL)),true)
        ifeq ($(strip $(USE_PAR_NETCDF)),true)
        	$(error SERIAL and USE_PAR_NETCDF cannot be true at the same time)
        endif
    endif
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.M *.mod *.d *.il core 

utils:
	cd UTILS; make

.SUFFIXES: $(SUFFIXES) .f90 .F90

.f90.o:
	$(FC) $(FFLAGS) -c $(INCLUDE) $<

.F90.o:
	$(FC) $(FFLAGS) -c $(INCLUDE) $<

analytic_mapping.o: analytic_semi_mapping.o analytic_spheroid_mapping.o \
	data_mesh.o global_parameters.o subpar_mapping.o Makefile ../make_axisem.macros 
analytic_semi_mapping.o: global_parameters.o Makefile ../make_axisem.macros 
analytic_spheroid_mapping.o: data_mesh.o global_parameters.o Makefile ../make_axisem.macros 
apply_masks.o: global_parameters.o Makefile ../make_axisem.macros 
attenuation.o: analytic_mapping.o commun.o data_io.o data_matr.o data_mesh.o \
	data_pointwise.o data_proc.o data_source.o data_spec.o data_time.o \
	get_mesh.o global_parameters.o pointwise_derivatives.o utlity.o Makefile ../make_axisem.macros 
clocks.o:  Makefile ../make_axisem.macros 
clocks_wrapper_solver.o: clocks.o data_io.o data_proc.o data_time.o Makefile ../make_axisem.macros 
data_comm.o: global_parameters.o list.o Makefile ../make_axisem.macros 
data_heterogeneous.o: global_parameters.o Makefile ../make_axisem.macros 
data_io.o: global_parameters.o Makefile ../make_axisem.macros 
data_matr.o: global_parameters.o Makefile ../make_axisem.macros 
data_mesh.o: data_io.o data_spec.o global_parameters.o Makefile ../make_axisem.macros 
data_pointwise.o: global_parameters.o Makefile ../make_axisem.macros 
data_proc.o:  Makefile ../make_axisem.macros 
data_source.o: global_parameters.o Makefile ../make_axisem.macros 
data_spec.o: global_parameters.o Makefile ../make_axisem.macros 
data_time.o: global_parameters.o Makefile ../make_axisem.macros 
def_grid.o: analytic_mapping.o apply_masks.o commun.o data_comm.o data_io.o \
	data_mesh.o data_pointwise.o data_proc.o data_spec.o get_mesh.o \
	global_parameters.o utlity.o Makefile ../make_axisem.macros 
def_precomp_terms.o: analytic_mapping.o attenuation.o background_models.o \
	commun.o data_io.o data_matr.o data_mesh.o data_pointwise.o \
	data_proc.o data_source.o data_spec.o def_grid.o get_mesh.o \
	get_model.o global_parameters.o pointwise_derivatives.o utlity.o coupling_mod.o Makefile ../make_axisem.macros 
get_mesh.o: background_models.o commun.o data_comm.o data_io.o data_mesh.o \
	data_proc.o data_time.o global_parameters.o Makefile ../make_axisem.macros 
global_parameters.o:  Makefile ../make_axisem.macros 
interpolation.o: global_parameters.o Makefile ../make_axisem.macros 
list.o: global_parameters.o Makefile ../make_axisem.macros 
main.o: clocks_wrapper_solver.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o def_grid.o get_mesh.o meshes_io.o nc_routines.o \
	nc_snapshots.o parameters.o time_evol_wave.o Makefile ../make_axisem.macros 
pointwise_derivatives.o: data_mesh.o data_pointwise.o data_spec.o \
	global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
rotations.o: data_io.o data_mesh.o data_proc.o data_source.o \
	global_parameters.o Makefile ../make_axisem.macros 
seismograms.o: commun.o data_io.o data_mesh.o data_proc.o data_source.o \
	data_time.o global_parameters.o nc_routines.o rotations.o utlity.o Makefile ../make_axisem.macros 
source.o: apply_masks.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o data_time.o global_parameters.o nc_routines.o \
	pointwise_derivatives.o utlity.o Makefile ../make_axisem.macros 
stiffness_di.o: attenuation.o data_matr.o data_mesh.o data_source.o \
	data_spec.o global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
stiffness_fluid.o: data_matr.o data_mesh.o data_source.o data_spec.o \
	global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
stiffness_mono.o: attenuation.o data_matr.o data_mesh.o data_source.o \
	data_spec.o global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
stiffness_quad.o: attenuation.o data_matr.o data_mesh.o data_source.o \
	data_spec.o global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
subpar_mapping.o: global_parameters.o Makefile ../make_axisem.macros 
unrolled_loops.o: data_mesh.o global_parameters.o Makefile ../make_axisem.macros 
utlity.o: analytic_mapping.o data_mesh.o data_spec.o global_parameters.o Makefile ../make_axisem.macros 
wavefields_io.o: data_io.o data_mesh.o data_pointwise.o data_proc.o \
	data_source.o global_parameters.o nc_routines.o nc_snapshots.o \
	pointwise_derivatives.o coupling_mod.o Makefile ../make_axisem.macros 
background_models.o: data_proc.o global_parameters.o interpolation.o Makefile ../make_axisem.macros 
commpi.o: clocks.o data_comm.o data_mesh.o data_proc.o data_time.o \
	global_parameters.o list.o Makefile ../make_axisem.macros 
commun.o: clocks.o commpi.o data_io.o data_mesh.o data_proc.o data_time.o \
	global_parameters.o Makefile ../make_axisem.macros 
get_model.o: background_models.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o data_time.o def_grid.o global_parameters.o \
	lateral_heterogeneities.o nc_routines.o utlity.o Makefile ../make_axisem.macros 
lateral_heterogeneities.o: background_models.o commun.o data_heterogeneous.o \
	data_io.o data_mesh.o data_proc.o data_source.o global_parameters.o \
	kdtree2.o utlity.o Makefile ../make_axisem.macros 
meshes_io.o: data_io.o data_matr.o data_mesh.o data_pointwise.o data_proc.o \
	data_source.o data_spec.o global_parameters.o nc_routines.o \
	nc_snapshots.o utlity.o Makefile ../make_axisem.macros 
nc_helpers.o: data_io.o data_proc.o global_parameters.o Makefile ../make_axisem.macros 
nc_routines.o: clocks.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o data_spec.o data_time.o global_parameters.o \
	nc_helpers.o Makefile ../make_axisem.macros 
nc_snapshots.o: data_io.o data_mesh.o data_proc.o data_source.o \
	global_parameters.o nc_helpers.o Makefile ../make_axisem.macros 
parameters.o: attenuation.o commun.o data_comm.o data_io.o data_mesh.o \
	data_proc.o data_source.o data_time.o global_parameters.o \
	nc_helpers.o nc_routines.o utlity.o Makefile ../make_axisem.macros 
time_evol_wave.o: apply_masks.o attenuation.o clocks.o commun.o data_io.o \
	data_matr.o data_mesh.o data_pointwise.o data_proc.o data_source.o \
	data_time.o def_precomp_terms.o global_parameters.o meshes_io.o \
	nc_routines.o parameters.o pointwise_derivatives.o rotations.o \
	seismograms.o source.o stiffness_di.o stiffness_fluid.o \
	stiffness_mono.o stiffness_quad.o wavefields_io.o coupling_mod.o Makefile ../make_axisem.macros 
ftz.o:  Makefile ../make_axisem.macros 
pthread.o:  Makefile ../make_axisem.macros 
kdtree2.o:  Makefile ../make_axisem.macros
coupling_mod.o: global_parameters.o utlity.o commun.o Makefile ../make_axisem.macros

