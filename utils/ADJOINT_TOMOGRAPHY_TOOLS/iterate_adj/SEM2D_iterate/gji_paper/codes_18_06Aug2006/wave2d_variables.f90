module wave2d_variables

  use wave2d_constants

! sound velocity
  double precision c

! material properties
  double precision kappal,mul,lambdal,lambdalplus2mul
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: rho,kappa,mu

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

! anchors
  double precision, dimension(NSPEC) :: z1,z2
  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(NGLOB) :: x, z, x_lon, z_lat, da
  integer, dimension(NGLOB) :: valence

! phase velocity map (CHT)
  double precision, dimension(NGLOB) :: c_glob
  double precision c0,cmin,cmax
  integer ihomo

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
  double precision, dimension(NGLOB) :: mass_global, mu_global, kappa_global, rho_global

! time marching
  double precision deltat,deltatover2,deltatsqover2
  double precision b_deltat,b_deltatover2,b_deltatsqover2
  double precision dh,time_step

! displacement, velocity and acceleration
  double precision, dimension(NCOMP,NGLOB) :: displ,veloc,accel
  double precision, dimension(NCOMP,NGLOB) :: b_displ,b_veloc,b_accel
  
! plotting
  double precision, dimension(NGLOB) :: norm

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

  integer :: imeasure
  double precision, dimension(:), allocatable :: measure_vec

! misfit function
  double precision :: chi(MAX_EVENT, MAX_SR, MAX_COMP, MAX_PHASE)

! kernel
  double precision, dimension(NGLOB) :: kernel_basis, kernel_basis_sum

! half duration of the source (put into wave2d_constants.f90)
!  double precision hdur

! file open or write status variable
  integer ios

! number of time steps to store wavefield
  integer NINT 

! input/output directory
  character(len=100) :: in_dir,out_dir
end module wave2d_variables
