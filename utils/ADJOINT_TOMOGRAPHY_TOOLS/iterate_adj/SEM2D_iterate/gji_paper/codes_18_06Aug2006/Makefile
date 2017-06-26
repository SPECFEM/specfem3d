F90 = ifort
#F90_FLAGS = -O2
F90_FLAGS = -O2
# first is for FFT, second is for bandpass filter
LIB = -L/opt/seismo-util/lib -lfftw3 -lDSacLib

# libraries for the original filter.f90
#LIB = -L/opt/seismo-util/lib -lDRWFiles -lDSacLib -lDSacio -lSacTools -lm

MOD = wave2d_constants wave2d_variables wave2d_define_der_matrices wave2d_solver wave2d_sub wave2d_sub2 wave2d_sub3
SUB = gll_library lagrange_poly  numerical_recipes

SRC_DIR = src
MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = .
MAIN = wave2d wave2d_2 wave2d_3 wave2d_4
MOD_FLAG = module

MOD_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(MOD))
F90_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(SUB))
OBJ = $(F90_OBJ) $(MOD_OBJ)

all : wave2d wave2d_2 wave2d_3 wave2d_4

$(MAIN) : % : $(SRC_DIR)/%.f90 $(F90_OBJ) $(MOD_OBJ)
	$(F90) -o $(BIN_DIR)/$* $(F90_FLAGS) $(SRC_DIR)/$*.f90 -$(MOD_FLAG) $(MOD_DIR) $(OBJ) $(LIB)

$(F90_OBJ): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(F90) -o $@ $(F90_FLAGS) -c $(SRC_DIR)/$*.f90

$(MOD_OBJ): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(F90) -o $@ $(F90_FLAGS) -c $(SRC_DIR)/$*.f90 -$(MOD_FLAG) $(MOD_DIR)

.PHONY : clean

clean:
	\rm -f *.o *.mod *~ $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod  *.gif *.ps $(MAIN)  plot*.csh *.eps *.cpt src/*~
cleanall:
	\rm -f OUTPUT/* *.o *.mod *~ $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod  *.gif *.jpg *.ps $(MAIN)  plot*.csh *.eps *.cpt


