rm -rf xcombine_vol_data *.o
mpif90 -O3 -c -o combine_vol_data_mod.o combine_vol_data_mod.f90
cc -c -o write_c_binary.o write_c_binary.c
mpif90 -O3 -o xcombine_vol_data combine_vol_data_mod.o write_c_binary.o
