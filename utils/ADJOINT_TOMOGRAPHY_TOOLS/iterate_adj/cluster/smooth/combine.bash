#for file in dm0; do for slice in `seq 1 3`; do xcombine_vol_data slice_file${slice} $file topo inout_smooth . 0; mv ${file}.mesh ${file}_${slice}.mesh; done; done

for file in kappa_kernel_smooth mu_kernel_smooth; do xcombine_vol_data slice_file $file topo inout_smooth . 0; done
