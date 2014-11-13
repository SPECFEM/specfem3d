#for file in Beta Bulk; do xcombine_vol_data slice_file $file topo OUTPUT_MODEL . 0; done

for file in vp_new vs_new; do xcombine_vol_data slice_file $file topo OUTPUT_MODEL . 0; done

#for file in vp vs; do xcombine_vol_data slice_file $file topo INPUT_MODEL . 0; done
