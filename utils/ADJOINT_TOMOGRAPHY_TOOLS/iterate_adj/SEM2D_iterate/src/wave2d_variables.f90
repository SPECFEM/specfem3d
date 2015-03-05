module wave2d_variables

  use wave2d_constants

! material properties and kernels
! For safety reasons, we keep the models for synthetics and data separate,
! and then assign one to rho-kappa-mu when we run the solver.
  double precision kappal,mul,lambdal,lambdalplus2mul
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: rho,kappa,mu
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: rho_dat, kappa_dat, mu_dat, alpha_dat, beta_dat, bulk_dat
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: rho_syn, kappa_syn, mu_syn, alpha_syn, beta_syn, bulk_syn
!  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: bulk_initial_syn, beta_initial_syn
!  double precision, dimension(NGLOB) :: mu_global, kappa_global, rho_global
!  double precision, dimension(NGLOB) :: alpha_global, beta_global

! temporary arrays
  double precision, dimension(NGLOB) :: temp_global1, temp_global2
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: temp_local1, temp_local2

! area associated with each GLL point
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: da_local
  double precision, dimension(NLOCAL) :: da_local_vec
  double precision, dimension(NGLOB) :: da_global

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz

! anchors (x1,x2,z1,z2 --> left,right,bottom,top
  double precision, dimension(NSPEC) :: z1,z2
  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(NGLOB) :: x, z
  integer, dimension(NGLOB) :: valence

! corners of the elements, in a global vector
  integer, dimension(NSPEC_CORNER) :: ielement_corner

! phase velocity map (CHT)
!  double precision, dimension(NGLOB) :: c_glob, c_glob_dat, c_glob_syn
  double precision, dimension(NGLOB) :: x_lon, z_lat, x_plot, z_plot
  double precision :: alpha0,beta0,rho0,kappa0,mu0,bulk0
!  double precision :: c,cmin,cmax
  integer :: Nfac
  double precision :: w_scale, afac

! 1D model
  integer :: nlayer
  double precision, dimension(10) :: z_breaks, r_layers, a_layers, b_layers

! gridpoints per wavelength
  double precision :: alpha_min, alpha_max, beta_min, beta_max

! Jacobian matrix and Jacobian
  double precision dxidxl,dxidzl,dgammadxl,dgammadzl,jacobianl
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: dxidx,dxidz,dgammadx,dgammadz,jacobian

! boundary elements: 1-xmin, 2-xmax, 3-zmin, 4-zmax (1,2,3,4 --> left, right, bottom, top)
! ibelm contains ispec indices for elements on the boundary of the grid
  integer, dimension(4,NELE) :: ibelm
  double precision, dimension(4,NGLL,NELE) :: jacobianb
  integer nspecb(4)

! absorbing boundary conditions
  integer j1, j2, ib, ibb
  double precision nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,weight

! local to global numbering
  integer, dimension(NGLLX,NGLLZ,NSPEC) :: ibool

! Mass matrix
  double precision mass_local
  double precision, dimension(NGLOB) :: mass_global

! time marching
  double precision deltat,deltatover2,deltatsqover2
  double precision b_deltat,b_deltatover2,b_deltatsqover2
  double precision dh,time_step

! displacement, velocity and acceleration
  double precision, dimension(NCOMP,NGLOB) :: displ,veloc,accel
  double precision, dimension(NCOMP,NGLOB) :: b_displ,b_veloc,b_accel

! plotting
!  double precision, dimension(NGLOB) :: norm

! space derivatives
  double precision tempx1l,tempx2l,tempy1l,tempy2l,tempz1l,tempz2l
  double precision fac1,fac2,hp1,hp2
  double precision dsxdxl,dszdxl,dsydxl,dsydzl,dsxdzl,dszdzl
  double precision sigma_xx,sigma_xy,sigma_xz,sigma_zx,sigma_zy,sigma_zz
  double precision, dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempy1,tempy2,tempz1,tempz2
  double precision, dimension(3,3) :: ds

  double precision b_tempx1l,b_tempx2l,b_tempy1l,b_tempy2l,b_tempz1l,b_tempz2l
  double precision b_dsxdxl,b_dszdxl,b_dsydxl,b_dsydzl,b_dsxdzl,b_dszdzl
  double precision b_sigma_xx,b_sigma_xy,b_sigma_xz,b_sigma_zx,b_sigma_zy,b_sigma_zz
  double precision, dimension(NGLLX,NGLLZ) :: b_tempx1,b_tempx2,b_tempy1,b_tempy2,b_tempz1,b_tempz2
  double precision, dimension(3,3) :: b_ds

  double precision :: utm_xmin, utm_zmin

! measurements (data)
!  integer :: imeasure
  double precision :: meas_pert
  double precision, dimension(:,:), allocatable :: measure_vec
  double precision, dimension(:), allocatable :: measure_pert_vec, cov_data

! indexing arrays
  integer, dimension(:,:,:), allocatable :: index_data
  integer, dimension(:,:), allocatable :: index_source
  integer, dimension(NVAR,2) :: m_inds

! misfit function
  double precision :: chi_data(MAX_EVENT, MAX_SR, MAX_COMP, MAX_PHASE), chi_val, chi_val_0
  double precision :: data_norm, chi_data_stop, chi_model_stop
  double precision :: model_norm_target, model_norm, model_norm_diff
  double precision :: model_norm_target_struct, model_norm_struct, model_norm_diff_struct
  double precision :: model_norm_target_source, model_norm_source, model_norm_diff_source
  double precision :: gradient_norm, gradient_norm_model, gradient_norm_data
  double precision :: gradient_norm_struct, gradient_norm_model_struct, gradient_norm_data_struct
  double precision :: gradient_norm_source, gradient_norm_model_source, gradient_norm_data_source
  double precision, dimension(NVAR) :: covm_weight_parts, covg_weight_parts
  double precision, dimension(NVAR) :: model_norm_target_parts, model_norm_parts, model_norm_diff_parts
  double precision, dimension(NVAR) :: gradient_norm_parts, gradient_norm_model_parts, &
      gradient_norm_data_parts, gbalance_parts
!  double precision :: chi_model_norm, chi_model_norm_struct, chi_model_norm_source
!  double precision :: chi_model_norm_target, chi_model_norm_target_struct, chi_model_norm_target_source
!  double precision :: chi_gradient_norm, chi_gradient_norm_struct, chi_gradient_norm_source
!  double precision :: chi_gradient_norm_data, chi_gradient_norm_data_struct, chi_gradient_norm_data_source
!  double precision :: chi_gradient_norm_model, chi_gradient_norm_model_struct, chi_gradient_norm_model_source
  double precision :: var_red_val

! kernels
!  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: kernel_basis_sum, kernel_basis
!  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: alpha_basis_sum, beta_basis_sum

! half duration of the source (put into wave2d_constants.f90)
!  double precision hdur

! file open or write status variable
  integer ios

! number of time steps to store wavefield
  integer NINT

! input/output directory
  character(len=100) :: in_dir,out_dir

end module wave2d_variables
