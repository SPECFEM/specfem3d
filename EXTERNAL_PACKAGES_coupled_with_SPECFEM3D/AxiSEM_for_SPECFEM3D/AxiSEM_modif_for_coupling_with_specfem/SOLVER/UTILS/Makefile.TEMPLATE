include ../../make_axisem.macros

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

all: xpost_processing xfield_transform

xpost_processing: post_processing.o nc_postroutines.o Makefile ../../make_axisem.macros 
	$(FC) $(LDFLAGS) post_processing.o nc_postroutines.o -o xpost_processing $(LIBS) 

post_processing.o: post_processing.F90 nc_postroutines.o Makefile ../../make_axisem.macros 
	$(FC) $(FFLAGS) $(INCLUDE) -c post_processing.F90 

nc_postroutines.o: nc_postroutines.F90 Makefile ../../make_axisem.macros
	$(FC) $(FFLAGS) $(INCLUDE) -c nc_postroutines.F90

xfield_transform: field_transform.o  Makefile ../../make_axisem.macros 
	$(FC) $(LDFLAGS) field_transform.o -o xfield_transform $(LIBS)

field_transform.o: field_transform.F90 Makefile ../../make_axisem.macros 
	$(FC) $(FFLAGS) $(INCLUDE) -c field_transform.F90

clean: 
	rm -rf *.o *.mod xpost_processing xfield_transform


