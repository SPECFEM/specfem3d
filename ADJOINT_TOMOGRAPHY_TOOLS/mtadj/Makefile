F90 = gfortran
F90_FLAGS = -fbounds-check
F77 = gfortran
F77_FLAGS = -ffixed-line-length-132

LIB = -L/opt/sac/lib/ -lsacio -lsac

MT_F90_SRC = mtadj_constants mtadj_variables mtadj_sub3 mtadj_sub2 mtadj_sub
MT_SRC = mtadj

MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = .


MT_F90_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(MT_F90_SRC))
MT_OBJ = $(MT_F77_OBJ) $(MT_F90_OBJ)

$(MT_SRC) : $(MT_SRC).f90 $(MT_OBJ)
	$(F90) -o $(BIN_DIR)/$@ $(F90_FLAGS) $@.f90 -I $(MOD_DIR) $(MT_OBJ) $(LIB)

$(MT_F90_OBJ): $(OBJ_DIR)/%.o : %.f90
	$(F90) -o $@ $(F90_FLAGS) -c $*.f90 -I $(MOD_DIR) -J$(MOD_DIR)

.PHONY : clean

clean:
	\rm -f *.o *.mod *~ $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod $(MT_SRC) OUTPUT/*

