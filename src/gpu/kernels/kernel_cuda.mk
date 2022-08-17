cuda_kernels_OBJS := \
	$O/add_source_main_rec_noise_cuda_kernel.cuda-kernel.o \
	$O/add_sources_ac_SIM_TYPE_2_OR_3_kernel.cuda-kernel.o \
	$O/add_sources_el_SIM_TYPE_2_OR_3_kernel.cuda-kernel.o \
	$O/assemble_boundary_accel_on_device.cuda-kernel.o \
	$O/assemble_boundary_potential_on_device.cuda-kernel.o \
	$O/compute_acoustic_seismogram_kernel.cuda-kernel.o \
	$O/compute_add_sources_acoustic_kernel.cuda-kernel.o \
	$O/compute_add_sources_kernel.cuda-kernel.o \
	$O/compute_coupling_acoustic_el_kernel.cuda-kernel.o \
	$O/compute_coupling_elastic_ac_kernel.cuda-kernel.o \
	$O/compute_coupling_ocean_cuda_kernel.cuda-kernel.o \
	$O/compute_dynamic_fault_cuda.cuda-kernel.o \
	$O/compute_elastic_seismogram_kernel.cuda-kernel.o \
	$O/compute_element_strain_cudakernel.cuda-kernel.o \
	$O/compute_kernels_acoustic_kernel.cuda-kernel.o \
	$O/compute_kernels_ani_cudakernel.cuda-kernel.o \
	$O/compute_kernels_cudakernel.cuda-kernel.o \
	$O/compute_kernels_hess_ac_cudakernel.cuda-kernel.o \
	$O/compute_kernels_hess_el_cudakernel.cuda-kernel.o \
	$O/compute_kernels_strength_noise_cuda_kernel.cuda-kernel.o \
	$O/compute_stacey_acoustic_kernel.cuda-kernel.o \
	$O/compute_stacey_elastic_kernel.cuda-kernel.o \
	$O/enforce_free_surface_cuda_kernel.cuda-kernel.o \
	$O/get_maximum_field_kernel.cuda-kernel.o \
	$O/get_maximum_vector_kernel.cuda-kernel.o \
	$O/Kernel_2_acoustic_impl.cuda-kernel.o \
	$O/Kernel_2_viscoelastic_impl.cuda-kernel.o \
	$O/kernel_3_accel_cuda_device.cuda-kernel.o \
	$O/kernel_3_acoustic_cuda_device.cuda-kernel.o \
	$O/kernel_3_cuda_device.cuda-kernel.o \
	$O/kernel_3_veloc_cuda_device.cuda-kernel.o \
	$O/noise_read_add_surface_movie_cuda_kernel.cuda-kernel.o \
	$O/prepare_boundary_accel_on_device.cuda-kernel.o \
	$O/prepare_boundary_potential_on_device.cuda-kernel.o \
	$O/process_smooth.cuda-kernel.o \
	$O/synchronize_boundary_accel_on_device.cuda-kernel.o \
	$O/transfer_surface_to_host_kernel.cuda-kernel.o \
	$O/UpdateDispVeloc_kernel.cuda-kernel.o \
	$O/UpdatePotential_kernel.cuda-kernel.o \
	$(EMPTY_MACRO)


