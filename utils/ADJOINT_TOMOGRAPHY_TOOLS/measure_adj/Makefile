F90 = gfortran
F90_FLAGS = -Wall -g -fbounds-check

#64-bit
F90_FLAGS = -O2 -fbounds-check     #-Wall

#32-bit
#F90_FLAGS = -O2 -fbounds-check -m32  #-Wall

# new version with default sac libraries
SACLIBDIR = ${SACHOME}/lib
LIB = -lsacio -lsac
#LIB = -L/opt/seismo/lib -lDRWFiles -lf90recipes -lDSacio -lDSacLib -lSacTools -lm



# NOTE: order matters if modules depend on other modules
MOD = ma_constants ma_variables ascii_rw ma_sub2 ma_sub

SRC_DIR = .
MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = .
MAIN = measure_adj rotate_adj_src
MOD_FLAG = module

MOD_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(MOD))
F90_OBJ = #$(patsubst %,$(OBJ_DIR)/%.o,$(SUB))
OBJ = $(F90_OBJ) $(MOD_OBJ)

all : measure_adj rotate_adj_src

$(MAIN) : % : $(SRC_DIR)/%.f90 $(F90_OBJ) $(MOD_OBJ)
	$(F90) -o $(BIN_DIR)/$* $(F90_FLAGS) $(SRC_DIR)/$*.f90 -I$(MOD_DIR) -J$(MOD_DIR) $(OBJ) -L${SACLIBDIR} $(LIB)

$(F90_OBJ): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(F90) -o $@ $(F90_FLAGS) -c $(SRC_DIR)/$*.f90

$(MOD_OBJ): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(F90) -o $@ $(F90_FLAGS) -c $(SRC_DIR)/$*.f90 -I$(MOD_DIR) -J$(MOD_DIR)

#rotate_adj_src: rotate_adj_src.f90
#	$(F90) -o rotate_adj_src rotate_adj_src.f90 -L/opt/seismo/lib -lDRWFiles

.PHONY : clean

clean:
	\rm -f *.o *.mod $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod  *.txt* STA.*  OUTPUT_FILES/*  *.sac rotate_adj_src measure_adj


