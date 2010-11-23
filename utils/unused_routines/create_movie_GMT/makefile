F90 = gfortran
FLAGS =  #-O3 
LIB =  
OBJ = obj/
BIN = .

SRC_FILES = create_movie_GMT \
            read_parameter_file compute_parameters \
            utm_geo  read_value_parameters get_value_parameters

OBJ_FILES = $(patsubst %, $(OBJ)%.o, $(SRC_FILES))

all : movie convert_nonlinear

movie : $(OBJ_FILES)
	$(F90) $(FLAGS)  -o $(BIN)/xcreate_movie_GMT $(OBJ_FILES) $(LIB)

$(OBJ)%.o : %.f90
	$(F90) $(FLAGS)  -o $(OBJ)$*.o -c $*.f90

convert_nonlinear: convert_nonlinear.f90
	$(F90) $(FLAGS)  -o $(BIN)/xconvert_nonlinear convert_nonlinear.f90 $(LIB)

all : movie convert_nonlinear

clean:
	\rm -f $(OBJ)*.o *~ xcreate_movie_GMT xconvert_nonlinear

