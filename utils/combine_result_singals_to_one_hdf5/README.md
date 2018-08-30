# combine result singals

A set of scripts and samples to combine seismo files into single binary file in HDF5 format.  

## Description of included files:  
* combine_signal_files_into_a_single_one.f90: fortran code which do the combine process.  
* test_signal_files.txt: a list of signal pathes that will be combined.
* test_script_convert_signals_to_one_single_binalize.py: sample script for the process.  
* test_script_open_binarized_file.py: sample script for opening the combined binary file.  
* compile_on_occigen.sh: useful script when one comile the fortran code on OCCIGEN/CINES. 
* test_signals: folder where the sample signals are put.  
  

## To use in your environment
1. compile combine_signal_files_into_a_single_one.f90 with compile flags for HDF5  
`h5pfc (or h5fc) combine_signal_files_into_a_single_one.f90 -o combine`  
(just do `./compile_on_occigen.sh` if you use it on OCCIGEN)  
2. prepare list_signal_files.txt which is a list of pathes for seismo signal files.  
3. adjust the values `total_timesteps` and `total_num_sigs` in `test_script_convert_signals_to_single_binary.py` to your simulation.  
4. `python test_script_convert_signals_to_single_binary.py` then the result will be saved as test_out.h5  
   
