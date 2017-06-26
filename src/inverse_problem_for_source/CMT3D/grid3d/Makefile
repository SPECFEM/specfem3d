F90 = gfortran
F90_FLAGS = -O3
F77 = gfortran
F77_FLAGS = -ffixed-line-length-132

LIB = /opt/sac/lib/libsacio.a

GRD_F90_SRC = get_cmt \
    grid3d_constants grid3d_variables grid3d_sub3 grid3d_sub2 grid3d_sub
GRD_F77_SRC = distaz
GRD_SRC = grid3d_flexwin

MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = .


GRD_F90_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(GRD_F90_SRC))
GRD_F77_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(GRD_F77_SRC))
GRD_OBJ = $(GRD_F77_OBJ) $(GRD_F90_OBJ)

$(GRD_SRC) : $(GRD_SRC).f90 $(GRD_OBJ)
	$(F90) -o $(BIN_DIR)/$@ $(F90_FLAGS) $@.f90 -I $(MOD_DIR) $(GRD_OBJ) $(LIB)

# notice that theortically you need to put the mod dependence here
$(GRD_F90_OBJ): $(OBJ_DIR)/%.o : %.f90
	$(F90) -o $@ $(F90_FLAGS) -c $*.f90 -I $(MOD_DIR) -J$(MOD_DIR)
#	@ \mv -f *.mod $(MOD_DIR)

$(GRD_F77_OBJ): $(OBJ_DIR)/%.o : %.f
	$(F77) -o $@ $(F77_FLAGS) -c $*.f




.PHONY : clean

clean:
	\rm -f *.o *.mod *~ $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod $(GRD_SRC)

