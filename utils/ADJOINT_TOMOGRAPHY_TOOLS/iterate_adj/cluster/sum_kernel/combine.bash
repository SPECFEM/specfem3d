#for file in dm0; do for slice in `seq 1 3`; do xcombine_vol_data ../topo_input/slice_file${slice} $file topo OUTPUT_SUM . 0; mv ${file}.mesh ${file}_${slice}.mesh; done; done

#for file in mu_kernel kappa_kernel; do xcombine_vol_data ../topo_input/slice_file $file topo OUTPUT_SUM . 0; done

for file in mu_kernel_smooth kappa_kernel_smooth; do xcombine_vol_data ../topo_input/slice_file $file topo OUTPUT_SUM . 0; done
