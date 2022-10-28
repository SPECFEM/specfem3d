PROG = xmesh

include ../make_axisem.macros

SRCS =	background_models.F90 test_bkgrdmodel.F90 analytic_mapping.f90 \
	analytic_semi_mapping.f90 analytic_spheroid_mapping.f90 clocks.f90 \
	clocks_wrapper_mesher.f90 data_bkgrdmodel.f90 data_coarse.f90 \
	data_diag.f90 data_gllmesh.f90 data_grid.f90 data_mesh.f90 \
	data_numbering.f90 data_pdb.f90 data_spec.f90 data_time.f90 \
	discont_meshing.f90 gllmeshgen.f90 global_parameters.f90 input.f90 \
	interpolation.f90 main.f90 mesh_info.f90 meshgen.f90 \
	model_discontinuities.f90 numbering.f90 parallelization.f90 pdb.f90 \
	sorting.f90 splib.f90 subpar_mapping.f90

OBJS =	background_models.o test_bkgrdmodel.o analytic_mapping.o \
	analytic_semi_mapping.o analytic_spheroid_mapping.o clocks.o \
	clocks_wrapper_mesher.o data_bkgrdmodel.o data_coarse.o data_diag.o \
	data_gllmesh.o data_grid.o data_mesh.o data_numbering.o data_pdb.o \
	data_spec.o data_time.o discont_meshing.o gllmeshgen.o \
	global_parameters.o input.o interpolation.o main.o mesh_info.o \
	meshgen.o model_discontinuities.o numbering.o parallelization.o pdb.o \
	sorting.o splib.o subpar_mapping.o

ifeq ($(strip $(USE_NETCDF)),true)
   FFLAGS += -Dunc
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



# cancel m2c implicit rule
%.o : %.mod


all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.M *.mod *.d *.il core

.SUFFIXES: $(SUFFIXES) .f90 .F90

.f90.o:
	$(FC) $(FFLAGS) -c $(INCLUDE) $<

.F90.o:
	$(FC) $(FFLAGS) -c $(INCLUDE) $<

analytic_mapping.o: analytic_semi_mapping.o analytic_spheroid_mapping.o \
	data_mesh.o subpar_mapping.o Makefile ../make_axisem.macros
analytic_semi_mapping.o: global_parameters.o Makefile ../make_axisem.macros
analytic_spheroid_mapping.o: global_parameters.o Makefile ../make_axisem.macros
clocks.o:  Makefile ../make_axisem.macros
clocks_wrapper_mesher.o: clocks.o data_time.o Makefile ../make_axisem.macros
data_bkgrdmodel.o: global_parameters.o Makefile ../make_axisem.macros
data_coarse.o:  Makefile ../make_axisem.macros
data_diag.o:  Makefile ../make_axisem.macros
data_gllmesh.o: global_parameters.o Makefile ../make_axisem.macros
data_grid.o: global_parameters.o Makefile ../make_axisem.macros
data_mesh.o: global_parameters.o Makefile ../make_axisem.macros
data_numbering.o:  Makefile ../make_axisem.macros
data_pdb.o: global_parameters.o Makefile ../make_axisem.macros
data_spec.o: global_parameters.o Makefile ../make_axisem.macros
data_time.o: global_parameters.o Makefile ../make_axisem.macros
discont_meshing.o: background_models.o data_bkgrdmodel.o data_diag.o \
	data_grid.o data_spec.o global_parameters.o model_discontinuities.o \
	splib.o Makefile ../make_axisem.macros
gllmeshgen.o: analytic_mapping.o data_diag.o data_gllmesh.o data_mesh.o \
	data_spec.o splib.o Makefile ../make_axisem.macros
global_parameters.o:  Makefile ../make_axisem.macros
input.o: background_models.o data_bkgrdmodel.o data_coarse.o data_diag.o \
	data_grid.o data_mesh.o data_spec.o global_parameters.o Makefile ../make_axisem.macros
interpolation.o: global_parameters.o Makefile ../make_axisem.macros
main.o: clocks.o clocks_wrapper_mesher.o data_bkgrdmodel.o data_grid.o \
	data_mesh.o data_time.o discont_meshing.o gllmeshgen.o input.o \
	mesh_info.o meshgen.o numbering.o parallelization.o pdb.o \
	test_bkgrdmodel.o Makefile ../make_axisem.macros
mesh_info.o: data_bkgrdmodel.o data_diag.o data_gllmesh.o data_grid.o \
	data_mesh.o data_numbering.o data_spec.o global_parameters.o Makefile ../make_axisem.macros
meshgen.o: analytic_spheroid_mapping.o data_bkgrdmodel.o data_coarse.o \
	data_diag.o data_grid.o data_mesh.o data_numbering.o data_pdb.o \
	data_spec.o global_parameters.o numbering.o Makefile ../make_axisem.macros
model_discontinuities.o: background_models.o data_bkgrdmodel.o data_diag.o \
	data_grid.o global_parameters.o Makefile ../make_axisem.macros
numbering.o: clocks.o data_diag.o data_gllmesh.o data_grid.o data_mesh.o \
	data_numbering.o data_spec.o data_time.o sorting.o Makefile ../make_axisem.macros
parallelization.o: data_bkgrdmodel.o data_coarse.o data_diag.o data_gllmesh.o \
	data_grid.o data_mesh.o data_pdb.o data_spec.o sorting.o \
	test_bkgrdmodel.o Makefile ../make_axisem.macros
pdb.o: background_models.o clocks.o data_bkgrdmodel.o data_diag.o \
	data_gllmesh.o data_grid.o data_mesh.o data_numbering.o data_pdb.o \
	data_spec.o data_time.o global_parameters.o numbering.o sorting.o Makefile ../make_axisem.macros
sorting.o: clocks.o data_time.o global_parameters.o Makefile ../make_axisem.macros
splib.o: global_parameters.o Makefile ../make_axisem.macros
subpar_mapping.o: global_parameters.o Makefile ../make_axisem.macros
background_models.o: global_parameters.o interpolation.o Makefile ../make_axisem.macros
test_bkgrdmodel.o: background_models.o data_bkgrdmodel.o data_diag.o \
	data_gllmesh.o data_grid.o data_mesh.o data_spec.o \
	global_parameters.o Makefile ../make_axisem.macros
