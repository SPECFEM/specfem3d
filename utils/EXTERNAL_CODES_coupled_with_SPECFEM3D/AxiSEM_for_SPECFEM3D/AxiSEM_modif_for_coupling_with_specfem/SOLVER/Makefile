PROG = axisem

include ../make_axisem.macros

SRCS =	commpi.F90 commun.F90 get_model.F90 lateral_heterogeneities.F90 \
	meshes_io.F90 nc_routines.F90 parameters.F90 time_evol_wave.F90 \
	analytic_mapping.f90 analytic_semi_mapping.f90 \
	analytic_spheroid_mapping.f90 apply_masks.f90 attenuation.f90 \
	background_models.f90 clocks.f90 data_comm.f90 data_heterogeneous.f90 \
	data_io.f90 data_matr.f90 data_mesh.f90 data_pointwise.f90 \
	data_proc.f90 data_source.f90 data_spec.f90 data_time.f90 \
	def_grid.f90 def_precomp_terms.f90 geom_transf.f90 get_mesh.f90 \
	global_parameters.f90 interpolation.f90 kdtree2.f90 list.f90 main.f90 \
	pointwise_derivatives.f90 rotations.f90 seismograms.f90 source.f90 \
	splib.f90 stiffness.f90 subpar_mapping.f90 unrolled_loops.f90 \
	utlity.f90 wavefields_io.f90 coupling_mod.f90 ftz.c pthread.c

OBJS =	commpi.o commun.o get_model.o lateral_heterogeneities.o meshes_io.o \
	nc_routines.o parameters.o time_evol_wave.o analytic_mapping.o \
	analytic_semi_mapping.o analytic_spheroid_mapping.o apply_masks.o \
	attenuation.o background_models.o clocks.o data_comm.o \
	data_heterogeneous.o data_io.o data_matr.o data_mesh.o \
	data_pointwise.o data_proc.o data_source.o data_spec.o data_time.o \
	def_grid.o def_precomp_terms.o geom_transf.o get_mesh.o \
	global_parameters.o interpolation.o kdtree2.o list.o main.o \
	pointwise_derivatives.o rotations.o seismograms.o source.o splib.o \
	stiffness.o subpar_mapping.o unrolled_loops.o utlity.o \
	wavefields_io.o coupling_mod.o ftz.o pthread.o

ifeq ($(strip $(USE_NETCDF)),true)
   FFLAGS += -Dunc
   LIBS = -L $(strip $(NETCDF_PATH))/lib -lnetcdff -Wl,-rpath,$(strip $(NETCDF_PATH))/lib
   INCLUDE = -I $(strip $(NETCDF_PATH))/include
else
   LIBS = 
   INCLUDE =
endif

ifeq ($(strip $(SERIAL)),true)
    FFLAGS += -Dserial
    LDFLAGS += -pthread
endif


# cancel m2c implicit rule 
%.o : %.mod 
 

all: $(PROG) 

$(PROG): $(OBJS)
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
	geom_transf.o get_mesh.o global_parameters.o pointwise_derivatives.o \
	utlity.o Makefile ../make_axisem.macros 
background_models.o: global_parameters.o interpolation.o Makefile ../make_axisem.macros 
clocks.o:  Makefile ../make_axisem.macros 
data_comm.o: global_parameters.o list.o Makefile ../make_axisem.macros 
data_heterogeneous.o: global_parameters.o Makefile ../make_axisem.macros 
data_io.o: global_parameters.o Makefile ../make_axisem.macros 
data_matr.o: global_parameters.o Makefile ../make_axisem.macros 
data_mesh.o: data_io.o global_parameters.o Makefile ../make_axisem.macros 
data_pointwise.o: global_parameters.o Makefile ../make_axisem.macros 
data_proc.o:  Makefile ../make_axisem.macros 
data_source.o: global_parameters.o Makefile ../make_axisem.macros 
data_spec.o: global_parameters.o Makefile ../make_axisem.macros 
data_time.o: global_parameters.o Makefile ../make_axisem.macros 
def_grid.o: apply_masks.o commun.o data_comm.o data_io.o data_mesh.o \
	data_pointwise.o data_proc.o data_spec.o geom_transf.o get_mesh.o \
	global_parameters.o splib.o utlity.o Makefile ../make_axisem.macros 
def_precomp_terms.o: analytic_mapping.o attenuation.o background_models.o \
	commun.o data_io.o data_matr.o data_mesh.o data_pointwise.o \
	data_proc.o data_source.o data_spec.o def_grid.o geom_transf.o \
	get_mesh.o get_model.o global_parameters.o pointwise_derivatives.o \
	splib.o utlity.o Makefile ../make_axisem.macros 
geom_transf.o: analytic_mapping.o data_mesh.o subpar_mapping.o Makefile ../make_axisem.macros 
get_mesh.o: background_models.o commun.o data_comm.o data_io.o data_mesh.o \
	data_proc.o data_time.o global_parameters.o Makefile ../make_axisem.macros 
global_parameters.o:  Makefile ../make_axisem.macros 
interpolation.o: global_parameters.o Makefile ../make_axisem.macros 
list.o: global_parameters.o Makefile ../make_axisem.macros 
main.o: clocks.o commun.o data_io.o data_mesh.o data_proc.o data_source.o \
	data_time.o def_grid.o get_mesh.o meshes_io.o nc_routines.o \
	parameters.o time_evol_wave.o Makefile ../make_axisem.macros 
pointwise_derivatives.o: data_mesh.o data_pointwise.o data_spec.o \
	global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
rotations.o: data_io.o data_mesh.o data_proc.o data_source.o \
	global_parameters.o Makefile ../make_axisem.macros 
seismograms.o: commun.o data_io.o data_mesh.o data_pointwise.o data_proc.o \
	data_source.o data_time.o global_parameters.o nc_routines.o \
	pointwise_derivatives.o rotations.o utlity.o Makefile ../make_axisem.macros 
source.o: apply_masks.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o data_spec.o data_time.o global_parameters.o \
	nc_routines.o pointwise_derivatives.o utlity.o Makefile ../make_axisem.macros 
splib.o: global_parameters.o Makefile ../make_axisem.macros 
stiffness.o: attenuation.o data_matr.o data_mesh.o data_source.o data_spec.o \
	global_parameters.o unrolled_loops.o Makefile ../make_axisem.macros 
subpar_mapping.o: global_parameters.o Makefile ../make_axisem.macros 
unrolled_loops.o: data_mesh.o global_parameters.o Makefile ../make_axisem.macros 
utlity.o: data_mesh.o data_spec.o geom_transf.o global_parameters.o Makefile ../make_axisem.macros 
wavefields_io.o: data_io.o data_mesh.o data_pointwise.o data_proc.o \
	data_source.o data_time.o global_parameters.o nc_routines.o \
	pointwise_derivatives.o Makefile ../make_axisem.macros 
commpi.o: clocks.o data_comm.o data_io.o data_mesh.o data_proc.o data_time.o \
	global_parameters.o list.o Makefile ../make_axisem.macros 
commun.o: clocks.o commpi.o data_mesh.o data_proc.o data_time.o \
	global_parameters.o Makefile ../make_axisem.macros 
get_model.o: background_models.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o data_time.o def_grid.o global_parameters.o \
	lateral_heterogeneities.o nc_routines.o utlity.o Makefile ../make_axisem.macros 
lateral_heterogeneities.o: background_models.o commun.o data_heterogeneous.o \
	data_io.o data_mesh.o data_proc.o data_source.o global_parameters.o \
	kdtree2.o utlity.o Makefile ../make_axisem.macros 
meshes_io.o: data_io.o data_matr.o data_mesh.o data_pointwise.o data_proc.o \
	data_source.o data_spec.o global_parameters.o nc_routines.o utlity.o Makefile ../make_axisem.macros 
nc_routines.o: clocks.o commun.o data_io.o data_mesh.o data_proc.o \
	data_source.o data_time.o global_parameters.o Makefile ../make_axisem.macros 
parameters.o: attenuation.o commun.o data_comm.o data_io.o data_mesh.o \
	data_proc.o data_source.o data_time.o global_parameters.o \
	nc_routines.o utlity.o Makefile ../make_axisem.macros 
time_evol_wave.o: apply_masks.o attenuation.o clocks.o commun.o data_io.o \
	data_matr.o data_mesh.o data_pointwise.o data_proc.o data_source.o \
	data_time.o def_precomp_terms.o global_parameters.o meshes_io.o \
	nc_routines.o parameters.o pointwise_derivatives.o rotations.o \
	seismograms.o source.o stiffness.o wavefields_io.o Makefile ../make_axisem.macros 
ftz.o:  Makefile ../make_axisem.macros 
pthread.o:  Makefile ../make_axisem.macros 
kdtree2.o:  Makefile ../make_axisem.macros
