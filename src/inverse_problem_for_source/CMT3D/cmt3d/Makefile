F90 = gfortran
F90_FLAGS = -O3
F77 = gfortran
F77_FLAGS = -ffixed-line-length-132

LIB = /opt/sac/lib/libsacio.a

CMT_F90_SRC = get_cmt  \
    cmt3d_constants cmt3d_variables cmt3d_sub4 cmt3d_sub3 cmt3d_sub2 cmt3d_sub
CMT_F77_SRC = cmt3d_utils distaz
CMT_SRC = cmt3d_flexwin

MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = .


CMT_F90_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(CMT_F90_SRC))
CMT_F77_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(CMT_F77_SRC))
CMT_OBJ = $(CMT_F77_OBJ) $(CMT_F90_OBJ)

$(CMT_SRC) : $(CMT_SRC).f90 $(CMT_OBJ)
	$(F90) -o $(BIN_DIR)/$@ $(F90_FLAGS) $@.f90 -I $(MOD_DIR) $(CMT_OBJ) $(LIB)

$(CMT_F90_OBJ): $(OBJ_DIR)/%.o : %.f90
	$(F90) -o $@ $(F90_FLAGS) -c $*.f90 -I $(MOD_DIR) -J$(MOD_DIR)

$(CMT_F77_OBJ): $(OBJ_DIR)/%.o : %.f
	$(F77) -o $@ $(F77_FLAGS) -c $*.f




.PHONY : clean

clean:
	\rm -f *.o *.mod *~ $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod $(CMT_SRC)

