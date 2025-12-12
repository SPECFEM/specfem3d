!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


! von Karman perturbation distribution
!
! apply to velocity models for investigating scattering effects

  module model_scattering_par

  use constants, only: myrank,CUSTOM_REAL,IMAIN,MAX_STRING_LEN

  use shared_parameters, only: SCATTERING_STRENGTH,SCATTERING_CORRELATION,SCATTERING_MATERIAL_IDS

  implicit none

  ! apply perturbation to different mesh materials
  integer, dimension(:), allocatable :: target_material_ids
  integer :: num_target_materials = 0
  logical :: USE_ALL_MATERIALS = .false.

  ! von Karman distribution
  ! perturbation array (regular x/y/z grid) for each process
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: perturbation_grid
  integer :: pert_Nx,pert_Ny,pert_Nz

  ! normalized fft cube mesh extends in range [0,1]
  double precision, parameter :: grid_length = 1.d0   ! range [0,1]

  ! flag to use small grid cell sizes
  logical, parameter :: USE_SMALLER_FFT_GRID_CELLS = .false.

  ! grid origin location, grid spacing
  double precision :: grid_origin_x,grid_origin_y,grid_origin_z
  double precision :: grid_dim_x,grid_dim_y,grid_dim_z
  double precision :: grid_dx

  double precision :: grid_cell_size                  ! fft grid cell size
  double precision :: grid_normalization_factor       ! normalization factor
  double precision :: grid_lambda_min                 ! normalized estimated minimum wavelength

  ! local (partial) grid dimensions for each process
  double precision :: grid_xmin_proc,grid_ymin_proc,grid_zmin_proc
  double precision :: grid_xmax_proc,grid_ymax_proc,grid_zmax_proc

  ! fft size
  integer :: grid_N

  ! note: not working yet, mesh point positions are determined on the fly and not yet known at the onset
  !       of velocity determination get_model() routine - todo for future use: shrink fft grid size per process

  ! special functions
  ! log2 interface
  interface log2
    module procedure rlog2,dlog2
  end interface

  contains

    !---------------------------------------------------
    ! logarithm base 2 (single precision)
    real function rlog2(x)

    implicit none
    real, intent(in) :: x

    rlog2 = log(x) / log(2.0)

    end function rlog2

    !---------------------------------------------------
    ! logarithm base 2 (double precision)
    double precision function dlog2(x)

    implicit none
    double precision, intent(in) :: x

    dlog2 = log(x) / log(2.d0)

    end function dlog2

  end module model_scattering_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_scattering_broadcast()

  ! standard routine to setup model

  use model_scattering_par

  implicit none

  ! user info
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'adding scattering perturbations:'
    call flush_IMAIN()
  endif

  ! initialize grid
  grid_normalization_factor = 0.d0
  grid_origin_x = 0.d0
  grid_origin_y = 0.d0
  grid_origin_z = 0.d0

  ! get target material ids from Par_file setting
  call read_target_material_ids()

  ! setup grid for ffts
  call setup_grid()

  ! generate perturbation distribution
  call generate_perturbations()

  end subroutine model_scattering_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_target_material_ids()

! get material ids from Par_file string

  use model_scattering_par

  implicit none
  ! local parameters
  integer :: i,ntokens,ier,istart,pos
  logical :: previous_is_delim
  character(len=1), parameter :: delimiter = ','     ! token string like "1,2"
  character(len=MAX_STRING_LEN) :: tmpstr

  ! parse tokens

  ! initializes
  num_target_materials = 0
  ntokens = 0
  tmpstr = SCATTERING_MATERIAL_IDS

  ! suppress leading white space and trim
  tmpstr = trim(adjustl(tmpstr))

  ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
  if (index(tmpstr,achar(13)) > 0) tmpstr = tmpstr(1:index(tmpstr,achar(13))-1)

  ! suppress trailing comments '#..'
  if (index(tmpstr,'#') > 0) tmpstr = tmpstr(1:index(tmpstr,'#')-1)

  ! checks if anything to do
  if (len_trim(tmpstr) > 0) then
    ! counts tokens
    ntokens = 1
    previous_is_delim = .true.
    do i = 1, len_trim(tmpstr)
      ! finds next delimiter (space or tabular)
      if (tmpstr(i:i) == delimiter) then
        if (.not. previous_is_delim) then
          ntokens = ntokens + 1
          previous_is_delim = .true.
        endif
      else
        if (previous_is_delim) previous_is_delim = .false.
      endif
    enddo
  endif

  ! set number of targets
  num_target_materials = ntokens

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  number of target materials: ',num_target_materials
    if (num_target_materials == 0) then
      write(IMAIN,*) '  no targets recognized, will proceed without scattering perturbations... '
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! check if anything left do to
  if (num_target_materials == 0) return

  ! setup array with target ids
  allocate(target_material_ids(num_target_materials),stat=ier)
  if (ier /= 0) stop 'Error allocating scattering target_material_ids'
  target_material_ids(:) = 0

  ! parse tokens
  istart = 1
  do i = 1,num_target_materials
    pos = index(tmpstr(istart:), ',')
    if (pos > 0) then
      read(tmpstr(istart:istart+pos-2), *) target_material_ids(i)
      istart = istart + pos
    else
      read(tmpstr(istart:len_trim(tmpstr)), *) target_material_ids(i)
    endif
  enddo

  !debug
  !print *, 'debug: scattering: target_material_ids =', target_material_ids(:),'tmpstr:',tmpstr

  ! target id 0 means to apply to all materials
  if (any(target_material_ids(:) == 0)) USE_ALL_MATERIALS = .true.

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  target material ids       : ',target_material_ids
    if (USE_ALL_MATERIALS) then
      write(IMAIN,*) '  using scattering perturbations for all materials'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_target_material_ids

!
!--------------------------------------------------------------------------------------------------
!

  subroutine setup_grid()

  use constants, only: NGLLX,HUGEVAL

  ! mesh element size
  use generate_databases_par, only: ibool,NSPEC_AB,NGLOB_AB,mat_ext_mesh

  use create_regions_mesh_ext_par, only: &
    xstore => xstore_unique, &
    ystore => ystore_unique, &
    zstore => zstore_unique

  use model_scattering_par

  implicit none
  real(kind=CUSTOM_REAL) :: x_min,x_min_all,y_min,y_min_all,z_min,z_min_all, &
                            x_max,x_max_all,y_max,y_max_all,z_max,z_max_all
  real(kind=CUSTOM_REAL) :: x_min_elem,y_min_elem,z_min_elem,x_max_elem,y_max_elem,z_max_elem
  real(kind=CUSTOM_REAL) :: dim_x,dim_y,dim_z,dim_max
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_slice,elemsize_min_norm
  integer :: ispec,imaterial_id

  ! checks if anything to do
  if (num_target_materials == 0) return

  ! get mesh properties
  if (USE_ALL_MATERIALS) then
    ! Calculation of whole computational domain
    x_min = minval(xstore(:))
    x_max = maxval(xstore(:))
    y_min = minval(ystore(:))
    y_max = maxval(ystore(:))
    z_min = minval(zstore(:))
    z_max = maxval(zstore(:))
  else
    ! only domains of target materials
    x_min = HUGEVAL
    y_min = HUGEVAL
    z_min = HUGEVAL

    x_max = - HUGEVAL
    y_max = - HUGEVAL
    z_max = - HUGEVAL

    do ispec = 1,NSPEC_AB
      ! material id: index 1 has associated material number
      imaterial_id = mat_ext_mesh(1,ispec)

      ! check only elements for target materials (target id 0 means to apply to all materials)
      if (any(target_material_ids(:) == imaterial_id)) then
        ! determine min/max extent of element
        call get_elem_minmax_xyz(x_min_elem,x_max_elem,y_min_elem,y_max_elem,z_min_elem,z_max_elem, &
                                 ispec,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

        ! domain min/max
        x_min = min(x_min,x_min_elem)
        y_min = min(y_min,y_min_elem)
        z_min = min(z_min,z_min_elem)

        x_max = max(x_max,x_max_elem)
        y_max = max(y_max,y_max_elem)
        z_max = max(z_max,z_max_elem)
      endif
    enddo
  endif

  ! store local grid dimensions
  grid_xmin_proc = x_min
  grid_xmax_proc = x_max
  grid_ymin_proc = y_min
  grid_ymax_proc = y_max
  grid_zmin_proc = z_min
  grid_zmax_proc = z_max

  ! determine min/max on all processes
  x_min_all = HUGEVAL
  y_min_all = HUGEVAL
  z_min_all = HUGEVAL

  x_max_all = - HUGEVAL
  y_max_all = - HUGEVAL
  z_max_all = - HUGEVAL

  call min_all_all_cr(x_min,x_min_all)
  call min_all_all_cr(y_min,y_min_all)
  call min_all_all_cr(z_min,z_min_all)

  call max_all_all_cr(x_max,x_max_all)
  call max_all_all_cr(y_max,y_max_all)
  call max_all_all_cr(z_max,z_max_all)

  ! dimensions
  dim_x = x_max_all - x_min_all
  dim_y = y_max_all - y_min_all
  dim_z = z_max_all - z_min_all

  if (dim_x < 0.0_CUSTOM_REAL) dim_x = 0.0_CUSTOM_REAL
  if (dim_y < 0.0_CUSTOM_REAL) dim_y = 0.0_CUSTOM_REAL
  if (dim_z < 0.0_CUSTOM_REAL) dim_z = 0.0_CUSTOM_REAL

  ! maximum mesh length
  dim_max = max(dim_x,dim_y,dim_z)

  ! store grid origin at minimum values
  grid_origin_x = x_min_all
  grid_origin_y = y_min_all
  grid_origin_z = z_min_all

  grid_dim_x = dim_x
  grid_dim_y = dim_y
  grid_dim_z = dim_z

  ! store normalization factor to get fft grid in range [0,1]
  grid_normalization_factor = dim_max

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  domain dimensions               : ',dim_x,'/',dim_y,'/',dim_z
    write(IMAIN,*) '  domain origins                  : ',x_min_all,'/',y_min_all,'/',z_min_all
    write(IMAIN,*) '  domain maximum extend           : ',dim_max
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check if domain volume found, otherwise target material ids won't match any in the mesh
  if (dim_x == 0.0_CUSTOM_REAL .or. dim_y == 0.0_CUSTOM_REAL .or. dim_z == 0.0_CUSTOM_REAL) then
    ! user info
    if (myrank == 0) then
      write(IMAIN,*) '  no valid domain found (target material ids not found in mesh)'
      write(IMAIN,*) '  no perturbations will be applied...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! set zero targets
    num_target_materials = 0
    ! all done
    return
  endif

  ! determine cell size for FFT grid
  !
  ! get minimum element size
  elemsize_min_slice = HUGEVAL
  do ispec = 1, NSPEC_AB
    ! material id: index 1 has associated material number
    imaterial_id = mat_ext_mesh(1,ispec)

    ! check only element sizes for target materials (target id 0 means to apply to all materials)
    if (any(target_material_ids(:) == imaterial_id) .or. USE_ALL_MATERIALS) then
      ! computes minimum and maximum size of this grid cell
      call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                               NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

      ! overall minimum element size in this slice
      elemsize_min_slice = min(elemsize_min_slice,elemsize_min)
    endif
  enddo
  ! fallback: if no target material was found this uses first element for size estimate
  if (elemsize_min_slice == HUGEVAL) then
    call get_elem_minmaxsize(elemsize_min,elemsize_max,1, &
                             NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)
    elemsize_min_slice = min(elemsize_min_slice,elemsize_min)
  endif

  ! minimum element size over all processes/slices
  call min_all_all_cr(elemsize_min_slice,elemsize_min)

  ! normalized element size
  elemsize_min_norm = real(elemsize_min / grid_normalization_factor, kind=CUSTOM_REAL)

  ! minimum grid cell size for FFT grid
  if (USE_SMALLER_FFT_GRID_CELLS) then
    ! using smaller cell sizes (makes FFTs much more expensive)
    ! at least half of the element size
    grid_cell_size = 0.5d0 * elemsize_min_norm
    ! or
    ! average GLL point distance for grid estimate
    !grid_cell_size = elemsize_min_norm / (NGLLX-1)
  else
    ! cell size ~ minimum element size
    grid_cell_size = elemsize_min_norm
  endif

  ! estimated minimum wavelength (assuming 5 grid points per wavelength), normalized
  grid_lambda_min = elemsize_min_norm / (NGLLX-1) * 5

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  domain minimum element size     : ',elemsize_min,'(m)'
    write(IMAIN,*) '  estimated minimum wavelength    : ',sngl(grid_lambda_min * grid_normalization_factor),'(m)'
    write(IMAIN,*) '  normalization factor            : ',sngl(grid_normalization_factor)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine setup_grid

!
!--------------------------------------------------------------------------------------------------
!


  real function psd_vonKarman_3D(a_in,kx,ky,kz)

  ! von Karman distribution
  !
  ! 3D power spectral density: Imperatori & Mai (2012)
  !
  ! P(k) = sigma**2 ( 2 sqrt(pi * a) )**E gamma(H + E/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + E/2) )
  !      = sigma**2 ( 2 sqrt(pi * a) )**3 gamma(H + 3/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + 3/2) )
  !
  ! E = euclidian dimension -> E = 3 for 3D medium
  !
  !from scipy.special import gamma

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: a_in,kx,ky,kz

  ! local parameters
  real(kind=CUSTOM_REAL) :: a,H,sigma
  real(kind=CUSTOM_REAL) :: g1,g2,amp,k,psd

  real(kind=CUSTOM_REAL), parameter :: PI = real(acos(-1.d0),kind=CUSTOM_REAL)
  real(kind=CUSTOM_REAL), parameter :: CONST_3HALF = real(3.d0 / 2.d0,kind=CUSTOM_REAL)

  a = a_in    ! correlation length
  H = 0.3     ! Hurst exponent: Imperatori & Mai use H = 0.1 and 0.3
  sigma = 0.1 ! standard deviation: 10%

  ! gamma function - a standard intrinsic function for Fortran 2008 and later
  !
  ! coefficients
  ! see: https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function
  ! with gamma(1/2) = sqrt(pi)
  !      gamma(1)   = 1
  !      gamma(3/2) = sqrt(pi)/2
  !
  g1 = gamma(H + CONST_3HALF)                            ! gamma function: gamma(H+3/2)
  g2 = gamma(H)                                          ! gamma function: gamma(H)
  amp = sigma * sigma * ( 2.0 * sqrt(PI * a) )**3        ! coefficient   : sigma**2 * ( 2 * sqrt(pi * a) )**3

  ! wavenumber k = sqrt(kx**2 + ky**2 + kz**2)
  k = sqrt(kx*kx + ky*ky + kz*kz)

  ! 3D von Karman power spectral density
  ! P(k) = sigma**2 ( 2 sqrt(pi * a) )**3 *  gamma(H + 3/2) / (gamma(H) * (1 + k**2 a**2 )**(H + 3/2)
  !      =             amp                *        g1       / (   g2    * (1 + k**2 a**2 )**(H + 3/2)
  psd = amp * g1 / ( g2 * (1.0 + k**2 * a**2)**(H + CONST_3HALF) )

  ! return value
  psd_vonKarman_3D = psd
  return

  end function psd_vonKarman_3D

!
!--------------------------------------------------------------------------------------------------
!

  subroutine generate_perturbations()

  use constants, only: PI,TINYVAL,itag
  use shared_parameters, only: SAVE_MESH_FILES,NPROC

  use model_scattering_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: a_corr
  real(kind=CUSTOM_REAL) :: dk,k_min,k_max,k_lambda,kx,ky,kz
  real(kind=CUSTOM_REAL) :: psd,rand_phase
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: freqs
  integer :: N,npower_of_2,index_k_lambda
  integer :: i,j,k,ik,ier
  double precision :: mb_size
  double precision :: grid_dx_in_m

  ! complex array in double precision (must match with fft.f90 precision)
  integer, parameter :: CUSTOM_CMPLX = 8
  complex(kind=CUSTOM_CMPLX),dimension(:,:,:),allocatable :: kxyz_dist
  complex(kind=CUSTOM_CMPLX) :: k_random

  ! helper array for FFT
  real(kind=CUSTOM_REAL) :: mpow(30)
  real(kind=CUSTOM_REAL) :: val_avg,val_max

  ! random seed
  integer :: num_seed
  integer,dimension(:),allocatable :: myseed

  ! external function
  real,external :: psd_vonKarman_3D

  ! white noise (otherwise von Karman noise distribution)
  logical :: USE_WHITE_NOISE = .false.
  real(kind=CUSTOM_REAL),external :: get_random_perturbation_value

  ! for mesh interpolations
  ! global perturbation array (for rank 0 only)
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: perturbation_grid_global
  integer :: pert_Nx_global,pert_Ny_global,pert_Nz_global

  ! local grids
  integer :: pert_imin, pert_imax, pert_jmin, pert_jmax, pert_kmin, pert_kmax
  ! local gathers
  integer :: iproc,imin,imax,jmin,jmax,kmin,kmax,Nx,Ny,Nz
  integer, dimension(9) :: pert_index_minmax
  integer, dimension(9,0:NPROC-1) :: pert_index_minmax_all
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: tmp_slice

  ! timing
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! checks if anything to do
  if (num_target_materials == 0) return

  ! white noise
  ! in case correlation length is smaller than GLL point distance, we could use white noise
  !
  ! correlation length
  ! such that k * a ~ 1 (for correlation factor == 1)
  !   a_corr = real(grid_lambda_min / (2.d0 * PI) * SCATTERING_CORRELATION,kind=CUSTOM_REAL)
  !
  ! therefore
  !   grid_lambda_min / (2.d0 * PI) * SCATTERING_CORRELATION < elemsize_min_norm / (NGLLX-1)
  ! using the estimates for grid_lambda_min = ( elemsize_min_norm / (NGLLX-1) * 5), this is equal to
  !   5 / (2.d0 * PI) * SCATTERING_CORRELATION < 1
  ! and
  !   SCATTERING_CORRELATION < 1 / 5 * 2 * PI ~ 1.25
  !
  ! for now, let's use white noise for factors == 0
  if (SCATTERING_CORRELATION < TINYVAL) then
    ! set flag
    USE_WHITE_NOISE = .true.
  endif

  ! user output
  if (myrank == 0) then
    if (USE_WHITE_NOISE) then
      write(IMAIN,*) '  perturbation distribution       : ','white noise (correlation factor == 0.0)'
    else
      write(IMAIN,*) '  perturbation distribution       : ','von Karman'
    endif
    write(IMAIN,*) '  perturbation correlation factor : ',sngl(SCATTERING_CORRELATION)
    write(IMAIN,*) '  perturbation maximum amplitude  : ',sngl(SCATTERING_STRENGTH)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

  ! user output
  if (myrank == 0) then
    ! FFT
    if (.not. USE_WHITE_NOISE) then
      write(IMAIN,*) '  FFT grid size (normalized)      : ',sngl(grid_length)
      write(IMAIN,*) '  FFT grid cell size (normalized) : ',sngl(grid_cell_size)
      call flush_IMAIN()
    endif
  endif

  ! next power of 2 for FFT
  npower_of_2 = int(log2(grid_length / grid_cell_size)) + 1
  N = 2**npower_of_2

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  estimated number of grid points : ',N
    call flush_IMAIN()
  endif

  ! limits minimum/maximum number of points for FFT
  !if (N < 1024) N = 1024  ! number of minimum points along x-direction (power of 2 for fft) 2**(10)
  !if (N < 512) N = 512    ! number of minimum points along x-direction (power of 2 for fft) 2**(9)
  if (N < 256) N = 256    ! number of minimum points along x-direction (power of 2 for fft) 2**(8)

  !if (N > 8192) N = 8192  ! number of maximum points along x-direction (power of 2 for fft) 2**(13)
  !if (N > 4096) N = 4096  ! number of maximum points along x-direction (power of 2 for fft) 2**(12)
  if (N > 2048) N = 2048  ! number of minimum points along x-direction (power of 2 for fft) 2**(11)

  ! power of 2
  npower_of_2 = ceiling(log2(real(N)))

  ! for mesh interpolations (uses double precision values)
  grid_dx       = grid_length / (N-1)  ! grid space increment, dx == dy == dz
  grid_N        = N

  ! grid space in m
  grid_dx_in_m = grid_dx * grid_normalization_factor

  ! user output
  if (myrank == 0) then
    ! FFT
    if (.not. USE_WHITE_NOISE) then
      write(IMAIN,*) '  FFT using number of grid points : N           = ',N
      write(IMAIN,*) '  FFT using power of 2            : npower_of_2 = ',npower_of_2
      ! memory required
      mb_size = dble(N) * dble(N) * dble(N) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  FFT using memory size           : ',sngl(mb_size),'MB'
    endif
    write(IMAIN,*) '  grid spacing                    : dx          = ',sngl(grid_dx_in_m),'(m)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! for global perturbation grid dimension (needed on rank 0)
  pert_Nx_global = ceiling(grid_dim_x / grid_dx_in_m) + 1
  pert_Ny_global = ceiling(grid_dim_y / grid_dx_in_m) + 1
  pert_Nz_global = ceiling(grid_dim_z / grid_dx_in_m) + 1

  ! make sure perturbation grid is not bigger than wavenumber distribution array
  if (pert_Nx_global > N) pert_Nx_global = N
  if (pert_Ny_global > N) pert_Ny_global = N
  if (pert_Nz_global > N) pert_Nz_global = N
  ! ensure minimum size for valid grid if domain is very small (1x1x1)
  if (pert_Nx_global < 1) pert_Nx_global = 1
  if (pert_Ny_global < 1) pert_Ny_global = 1
  if (pert_Nz_global < 1) pert_Nz_global = 1

  ! for local (partial) perturbation grids on each process
  pert_imin = int((grid_xmin_proc - grid_origin_x) / grid_dx_in_m) + 1
  pert_imax = int((grid_xmax_proc - grid_origin_x) / grid_dx_in_m) + 1
  pert_jmin = int((grid_ymin_proc - grid_origin_y) / grid_dx_in_m) + 1
  pert_jmax = int((grid_ymax_proc - grid_origin_y) / grid_dx_in_m) + 1
  pert_kmin = int((grid_zmin_proc - grid_origin_z) / grid_dx_in_m) + 1
  pert_kmax = int((grid_zmax_proc - grid_origin_z) / grid_dx_in_m) + 1

  ! add margin elements to extend local grids (such that they overlap on shared nodes to avoid interpolation artefacts)
  pert_imin = pert_imin - 1
  pert_jmin = pert_jmin - 1
  pert_kmin = pert_kmin - 1

  pert_imax = pert_imax + 1
  pert_jmax = pert_jmax + 1
  pert_kmax = pert_kmax + 1

  ! handle egde cases
  if (pert_imin < 1) pert_imin = 1
  if (pert_jmin < 1) pert_jmin = 1
  if (pert_kmin < 1) pert_kmin = 1

  if (pert_imax < 1) pert_imax = 1
  if (pert_jmax < 1) pert_jmax = 1
  if (pert_kmax < 1) pert_kmax = 1

  if (pert_imin > N) pert_imin = N
  if (pert_jmin > N) pert_jmin = N
  if (pert_kmin > N) pert_kmin = N

  if (pert_imax > N) pert_imax = N
  if (pert_jmax > N) pert_jmax = N
  if (pert_kmax > N) pert_kmax = N

  ! stay within global mesh extend
  if (pert_imax > pert_Nx_global) pert_imax = pert_Nx_global
  if (pert_jmax > pert_Ny_global) pert_jmax = pert_Ny_global
  if (pert_kmax > pert_Nz_global) pert_kmax = pert_Nz_global

  ! update local perturbation grid xmin/xmax/.. positions
  grid_xmin_proc = grid_origin_x + (pert_imin-1) * grid_dx_in_m
  grid_xmax_proc = grid_origin_x + (pert_imax-1) * grid_dx_in_m
  grid_ymin_proc = grid_origin_y + (pert_jmin-1) * grid_dx_in_m
  grid_ymax_proc = grid_origin_y + (pert_jmax-1) * grid_dx_in_m
  grid_zmin_proc = grid_origin_z + (pert_kmin-1) * grid_dx_in_m
  grid_zmax_proc = grid_origin_z + (pert_kmax-1) * grid_dx_in_m

  ! total number of grid element for this slice
  pert_Nx = pert_imax - pert_imin + 1
  pert_Ny = pert_jmax - pert_jmin + 1
  pert_Nz = pert_kmax - pert_kmin + 1

  ! in case where a process might not own any target material regions
  if (pert_Nx < 1) pert_Nx = 0
  if (pert_Ny < 1) pert_Ny = 0
  if (pert_Nz < 1) pert_Nz = 0

  if (pert_Nx > N) pert_Nx = N
  if (pert_Ny > N) pert_Ny = N
  if (pert_Nz > N) pert_Nz = N

  ! wavenumbers
  ! min/max: k_max = 2 pi / (2 dx)
  k_min = real(2.d0 * PI / (N * grid_dx),kind=CUSTOM_REAL)
  k_max = real(2.d0 * PI / (2.0 * grid_dx),kind=CUSTOM_REAL)

  ! wavenumber increment
  dk = real(2.d0 * PI / (N * grid_dx),kind=CUSTOM_REAL)

  ! maximum index for k_max
  !index_k_max = N / 2 # k_max / dk = (2 pi / (2 dx) ) / (2 pi / (N dx) ) = (N * dx) / (2 * dx) = N / 2

  ! index for k_lambda_min
  k_lambda = real(2.d0 * PI / (grid_lambda_min),kind=CUSTOM_REAL)
  index_k_lambda = int(k_lambda / dk)

  ! correlation length
  ! such that k * a ~ 1 (for correlation factor == 1)
  a_corr = real(grid_lambda_min / (2.d0 * PI) * SCATTERING_CORRELATION,kind=CUSTOM_REAL)

  ! user output
  if (myrank == 0) then
    ! FFT
    if (.not. USE_WHITE_NOISE) then
      write(IMAIN,*) '  wavenumbers: k min/max          = ',k_min,'/',k_max
      write(IMAIN,*) '               dk increment       = ',dk
      write(IMAIN,*) '               k lambda_min       = ',k_lambda
      write(IMAIN,*) '               index k lambda_min = ',index_k_lambda
      write(IMAIN,*)
    endif
    write(IMAIN,*) '  correlation length              : ',a_corr
    write(IMAIN,*) '                       (in meters): ',sngl(a_corr * grid_normalization_factor),'(m)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! setup wavenumber arrays
  if (myrank == 0) then
    ! only main process creates pertubation grid, then distributes it to all others
    ! to have the same perturbation grid for all processes.

    ! initializes random number generator
    ! seed with fixed value to make successive runs repeatable
    call random_seed(size=num_seed)
    allocate(myseed(num_seed))
    myseed(1:num_seed) = 12345
    call random_seed(put=myseed)

    ! debug
    !call random_number(rand_phase)
    !print *,'debug: rank ',myrank,'random number 1: ',rand_phase

    ! FFT
    if (.not. USE_WHITE_NOISE) then
      ! 3D arrays
      !allocate(freqs(N), wavenumbers_kx(N,N,N), wavenumbers_ky(N,N,N), wavenumbers_kz(N,N,N),stat=ier)
      allocate(freqs(N),stat=ier)
      if (ier /= 0) stop 'Error allocating wavenumbers arrays'

      ! fft indexing
      freqs(:) = 0.0
      do i = 1,N
        ! numpy-like
        !if (i <= N/2+1) then
        !  ! indices: 0,1,..,N/2
        !  ik = i - 1
        !else
        !  ! indices: -N/2-1,-N/2-2,..,-1
        !  ik = - (N - (i - 1))
        !endif

        ! here positive values only; will take conjugates later for second half of array (when calling to apply symmetries)
        if (i <= N/2+1) then
          ! indices: 0,1,..,N/2
          ik = i - 1
        else
          ! indices: N/2-1,N/2-2,..,1
          ik = N - (i - 1)
        endif
        freqs(i) = real(ik,kind=CUSTOM_REAL)     ! for d = 1.0/N: freqs = ik / (d*N) = ik / ((1.0 / N) * N) = ik / ( 1.0 ) = ik
      enddo
      ! debug
      !print *,'debug: freqs ',freqs(:); print *
      ! same as:
      !freqs = (/0.0,(real(i,kind=CUSTOM_REAL),i=1,N/2)/)
      !freqs = (/freqs(1:N/2+1),(real(i,kind=CUSTOM_REAL),i=N/2-1,1,-1)/)
      !print *,'debug: new freqs ',freqs(:); print *
      ! debug
      !print *,'maximum wavenumber = ',maxval(freqs(:)) * dk ! maximum wavenumber

      ! sets up helper array for FFTs
      do i = 1,npower_of_2
        mpow(i) = 2**(npower_of_2-i)
      enddo

      ! wavenumber distribution (for FFTs)
      allocate(kxyz_dist(N,N,N),stat=ier)
      if (ier /= 0) stop 'Error allocating kxyz_dist array'
      kxyz_dist(:,:,:) = cmplx(0.0,0.0)

      ! applies amplitudes, which follow defined power spectral density (psd), to random phases
      do k = 1,N
        do j = 1,N
          do i = 1,N
            ! wavenumbers
            kx = freqs(i) * dk  ! (2.0 * np.pi * freqs[icol]) / L
            ky = freqs(j) * dk  ! (2.0 * np.pi * freqs[irow]) / L
            kz = freqs(k) * dk  ! (2.0 * np.pi * freqs[iz]) / L

            ! amplitudes w/ von Karman distribution
            psd = psd_vonKarman_3D(a_corr,kx,ky,kz)

            ! random phase
            call random_number(rand_phase)
            ! range [0,2pi]
            rand_phase = real(rand_phase * 2.d0 * PI,kind=CUSTOM_REAL)
            k_random = cmplx( cos(rand_phase), sin(rand_phase) )

            ! stores wavenumber distribution
            kxyz_dist(i,j,k) = k_random * sqrt(psd)
          enddo
        enddo
      enddo
      ! free memory
      deallocate(freqs)

      ! user output
      write(IMAIN,*) '  starting 3D FFTs'
      call flush_IMAIN()

      ! define symmetry conditions for 3D FFT
      call fft_apply_3D_symmetry(kxyz_dist,N)

      ! FFT arrays
      ! example: 1D fft
      !allocate(k_line(N),x_FFT(N),stat=ier)
      !if (ier /= 0) stop 'Error allocating x_FFT array'
      !k_line(:) = cmplx(0.0,0.0)
      !x_FFT(:) = 0.0_CUSTOM_REAL
      ! inverse Fourier transform
      ! w/ 1D FFTs
      !do k = 1,N
      !  do j = 1,N
      !    ! takes 1D line
      !    k_line(:) = kxyz_dist(:,j,k)
      !    ! pad negative k
      !    do ii = 2, N/2
      !      ! fills from N,N-1,..,N/2+2
      !      k_line(N+2-ii) = conjg(k_line(ii))
      !    enddo
      !    ! 1D inverse FFT
      !    call FFTinv(npower_of_2, k_line(:), 1.0_CUSTOM_REAL, dk, x_FFT(:), mpow) ! inverse FFT, outputs real array x_FFT
      !
      !    !call rspec(k_line(:),N/2)                                    ! restructuring
      !    !call FFT(npower_of_2, k_line(:), -1.0_CUSTOM_REAL, dk, mpow) ! inverse FFT, outputs complex array k_line
      !    !x_FFT(1:N) = real(k_line(1:N))                               ! takes the real part
      !
      !    ! stores perturbations array
      !    perturbation_grid(:,j,k) = x_FFT(:)
      !  enddo
      !enddo

      ! 3D FFT
#ifdef USE_FFTW
      call FFT_3D_FFTW(N, kxyz_dist, -1.0_CUSTOM_REAL)
#else
      call FFT_3D(N, npower_of_2, kxyz_dist, -1.0_CUSTOM_REAL, dk, mpow) ! inverse 3D FFT (zign == -1)
#endif
    endif

    ! perturbations
    ! user output
    write(IMAIN,*)
    write(IMAIN,*) '  perturbations: global grid Nx/Ny/Nz = ',pert_Nx_global,'/',pert_Ny_global,'/',pert_Nz_global
    ! memory required
    mb_size = dble(pert_Nx_global) * dble(pert_Ny_global) * dble(pert_Nz_global) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
    write(IMAIN,*) '                 memory size   = ',sngl(mb_size),'MB'
    write(IMAIN,*)
    call flush_IMAIN()

    ! perturbation grid array (for model perturbations)
    allocate(perturbation_grid_global(pert_Nx_global,pert_Ny_global,pert_Nz_global),stat=ier)
    if (ier /= 0) stop 'Error allocating global perturbation_grid array'
    perturbation_grid_global(:,:,:) = 0.0_CUSTOM_REAL

    ! stores real part
    do k = 1,pert_Nz_global
      do j = 1,pert_Ny_global
        do i = 1,pert_Nx_global
          ! stores perturbations array values (in range [-1,1])
          if (USE_WHITE_NOISE) then
            ! white noise
            perturbation_grid_global(i,j,k) = get_random_perturbation_value()
          else
            ! from FFTs (real part)
            perturbation_grid_global(i,j,k) = real(kxyz_dist(i,j,k),kind=CUSTOM_REAL)
          endif
        enddo
      enddo
    enddo

    ! user output
    write(IMAIN,*) '  applying scaling: maximum amplitude = ',SCATTERING_STRENGTH
    call flush_IMAIN()

    ! debug
    !print *,'before scattering perturbation: min/max = ',minval(perturbation_grid),'/',maxval(perturbation_grid)
    !print *,'                                average = ',sum(perturbation_grid)/(N*N*N)

    ! applies scaling
    !
    ! note: let's get the average and absolute maximum value from the actual perturbation_grid() array.
    !       this might differ a bit from the full kxyz_dist() array, depending on the chosen (partial)
    !       mesh domain to cover with the grid perturbations. Still, we want the applied perturbations to have
    !       a zero mean and scaled such that the maximum perturbation becomes the scattering strength defined by the user.
    !
    call get_grid_average_max(perturbation_grid_global,pert_Nx_global,pert_Ny_global,pert_Nz_global,val_avg,val_max)

    ! makes sure it has a zero average
    perturbation_grid_global(:,:,:) = perturbation_grid_global(:,:,:) - val_avg

    ! normalizes to range [-1,1]
    if (val_max > 0.0_CUSTOM_REAL) then
      perturbation_grid_global(:,:,:) = perturbation_grid_global(:,:,:) / val_max
    endif

    ! scales with maximum strength
    perturbation_grid_global(:,:,:) = perturbation_grid_global(:,:,:) * real(SCATTERING_STRENGTH,kind=CUSTOM_REAL)

    ! debug
    !print *,'scattering perturbation: min/max = ',minval(perturbation_grid),'/',maxval(perturbation_grid)
    !print *,'                         average = ',sum(perturbation_grid)/(N*N*N)
    !print *,'debug: loop done'

    ! FFT
    if (.not. USE_WHITE_NOISE) then
      ! free memory
      deallocate(kxyz_dist)
    endif
  endif

  ! synchronize all processes
  call synchronize_all()

  ! distribute perturbation grid to all processes
  ! allocates local perturbation grid on each process
  allocate(perturbation_grid(pert_Nx,pert_Ny,pert_Nz),stat=ier)
  if (ier /= 0) stop 'Error allocating perturbation_grid array'
  perturbation_grid(:,:,:) = 0.0_CUSTOM_REAL

  ! get dimensions from all processes
  pert_index_minmax(1) = pert_imin
  pert_index_minmax(2) = pert_imax
  pert_index_minmax(3) = pert_jmin
  pert_index_minmax(4) = pert_jmax
  pert_index_minmax(5) = pert_kmin
  pert_index_minmax(6) = pert_kmax
  pert_index_minmax(7) = pert_Nx
  pert_index_minmax(8) = pert_Ny
  pert_index_minmax(9) = pert_Nz

  ! gather all index bounds (on main process only)
  pert_index_minmax_all(:,:) = 0
  call gather_all_i(pert_index_minmax,9,pert_index_minmax_all,9,NPROC)

  ! main process distributes grid to all other arrays
  if (myrank == 0) then
    ! main process
    ! setup own grid first
    if (pert_Nx > 0 .and. pert_Ny > 0 .and. pert_Nz > 0) then
      ! slice array portion
      perturbation_grid(1:pert_Nx,1:pert_Ny,1:pert_Nz) = perturbation_grid_global(pert_imin:pert_imax, &
                                                                                  pert_jmin:pert_jmax, &
                                                                                  pert_kmin:pert_kmax)
    endif

    ! send partial grids to all other processes
    do iproc = 1,NPROC-1
      ! grid dimension of process
      imin = pert_index_minmax_all(1,iproc)
      imax = pert_index_minmax_all(2,iproc)
      jmin = pert_index_minmax_all(3,iproc)
      jmax = pert_index_minmax_all(4,iproc)
      kmin = pert_index_minmax_all(5,iproc)
      kmax = pert_index_minmax_all(6,iproc)
      Nx = pert_index_minmax_all(7,iproc)
      Ny = pert_index_minmax_all(8,iproc)
      Nz = pert_index_minmax_all(9,iproc)
      ! send slice
      if (Nx > 0 .and. Ny > 0 .and. Nz > 0) then
        ! temporary slice for sending
        allocate(tmp_slice(Nx,Ny,Nz),stat=ier)
        if (ier /= 0) stop 'Error allocating temporary slice array'
        ! setup slice
        tmp_slice(:,:,:) = perturbation_grid_global(imin:imax,jmin:jmax,kmin:kmax)
        ! or explicit
        !tmp_slice(:,:,:) = 0.0_CUSTOM_REAL
        !do k = 1,Nz
        !  do j = 1,Ny
        !    do i = 1,Nx
        !      tmp_slice(i,j,k) = perturbation_grid_global(imin+i-1,jmin+j-1,kmin+k-1)
        !    enddo
        !  enddo
        !enddo

        ! send to process iproc
        call sendv_cr(tmp_slice,Nx*Ny*Nz,iproc,itag)

        ! free memory
        deallocate(tmp_slice)
      endif
    enddo

    ! deallocate global grid on rank 0 after distribution
    deallocate(perturbation_grid_global)
  else
    ! secondary process
    ! receive slice from rank 0
    if (pert_Nx > 0 .and. pert_Ny > 0 .and. pert_Nz > 0) then
      call recvv_cr(perturbation_grid,pert_Nx*pert_Ny*pert_Nz,0,itag)
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  perturbations (slice 0): min/max = ',minval(perturbation_grid),'/',maxval(perturbation_grid)
    write(IMAIN,*) '                           average = ', &
                   sum(perturbation_grid) / real(pert_Nx * pert_Ny * pert_Nz,kind=CUSTOM_REAL)
    write(IMAIN,*)
    ! timing
    tCPU = wtime() - time_start
    write(IMAIN,*) '  Elapsed time for perturbation setup in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! visualization output
  if (SAVE_MESH_FILES) then
    ! output VTU file for visual inspection
    call plot_grid_data()
  endif

  end subroutine generate_perturbations

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_grid_average_max(array,Nx,Ny,Nz,val_avg,val_max)

  use constants, only: CUSTOM_REAL,SIZE_DOUBLE

  implicit none

  integer, intent(in) :: Nx,Ny,Nz
  real(kind=CUSTOM_REAL),dimension(Nx,Ny,Nz), intent(in) :: array
  real(kind=CUSTOM_REAL), intent(out) :: val_avg,val_max

  ! local parameters
  integer :: i,j,k
  double precision :: davg,dmax,dval

  ! initializes
  val_avg = 0.0_CUSTOM_REAL
  val_max = 0.0_CUSTOM_REAL

  ! checks if anything to do
  if (Nx == 0 .or. Ny == 0 .or. Nz == 0) return

  ! in custom real
  !val_avg = sum(perturbation_grid) / real(pert_Nx * pert_Ny * pert_Nz, kind=CUSTOM_REAL)
  !val_max = maxval(abs(perturbation_grid))

  ! determine average and maximum value (in double precision)
  davg = 0.d0
  dmax = 0.d0
  do k = 1,Nz
    do j = 1,Ny
      do i = 1,Nx
        ! in double precision
        dval = real(array(i,j,k),kind=SIZE_DOUBLE)
        ! for avg
        davg = davg + dval
        ! for max
        dmax = max(dmax,abs(dval))
      enddo
    enddo
  enddo
  ! takes average
  davg = davg / dble(Nx * Ny * Nz)

  ! return as custom real values
  val_avg = real(davg, kind=CUSTOM_REAL)
  val_max = real(dmax, kind=CUSTOM_REAL)

  end subroutine get_grid_average_max

!
!-------------------------------------------------------------------------------------------------
!

  function get_random_perturbation_value() result(pert_val)

  use constants, only: CUSTOM_REAL

  implicit none
  real(kind=CUSTOM_REAL) :: pert_val
  real(kind=CUSTOM_REAL) :: rand_val

  ! random value between [0,1]
  call random_number(rand_val)

  ! white noise between [-1,1]
  pert_val = 2.0_CUSTOM_REAL * rand_val - 1.0_CUSTOM_REAL

  end function get_random_perturbation_value

!
!-------------------------------------------------------------------------------------------------
!

  subroutine plot_grid_data()

  use constants, only: myrank,IMAIN,CUSTOM_REAL,MAX_STRING_LEN
  use shared_parameters, only: LOCAL_PATH

  use model_scattering_par, only: grid_xmin_proc,grid_ymin_proc,grid_zmin_proc, &
    grid_dx,grid_normalization_factor

  use model_scattering_par, only: &
    array => perturbation_grid, &
    Nx => pert_Nx, &
    Ny => pert_Ny, &
    Nz => pert_Nz

  implicit none
  ! local parameters
  integer :: ne,np,ixyz,n1,n2,n3,n4,n5,n6,n7,n8
  integer :: i,j,k,ier

  ! global point data
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: total_dat
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: total_dat_xyz
  real(kind=CUSTOM_REAL) :: cell_size
  integer,dimension(:,:),allocatable :: total_dat_con
  ! file output
  character(len=MAX_STRING_LEN) :: mesh_file,var_name
  character(len=MAX_STRING_LEN) :: prname

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  saving perturbation grid...'
    call flush_IMAIN()
  endif

  ! each process plots out its own perturbation grid
  ! checks if anything to do
  if (Nx == 0 .or. Ny == 0 .or. Nz == 0) return

  ! regular grid
  np = Nx * Ny * Nz             ! total number of points
  ne = (Nx-1) * (Ny-1) * (Nz-1) ! total number of elements

  ! creates array to hold point data
  allocate(total_dat(np),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat array'
  total_dat(:) = 0.0_CUSTOM_REAL
  allocate(total_dat_xyz(3,np),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_xyz array'
  total_dat_xyz(:,:) = 0.0_CUSTOM_REAL
  allocate(total_dat_con(8,ne),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_con array'
  total_dat_con(:,:) = 0

  ! cell size (same in x/y/z direction)
  cell_size = grid_dx * grid_normalization_factor

  ! regular grid points values
  np = 0
  do k = 1,Nz
    do j = 1,Ny
      do i = 1,Nx
        ! number counter
        np = np + 1

        ! checks point index
        ixyz = i + (j-1) * (Nx) + (k-1) * (Nx * Ny)
        if (np /= ixyz) stop 'Invalid grid point'

        ! array value
        total_dat(np) = array(i,j,k)

        ! grid point position [-L/2,L/2]
        total_dat_xyz(1,np) = real(grid_xmin_proc + (i-1) * cell_size,kind=CUSTOM_REAL)
        total_dat_xyz(2,np) = real(grid_ymin_proc + (j-1) * cell_size,kind=CUSTOM_REAL)
        total_dat_xyz(3,np) = real(grid_zmin_proc + (k-1) * cell_size,kind=CUSTOM_REAL)
      enddo
    enddo
  enddo
  ! element connections
  ne = 0
  do k = 1,Nz-1
    do j = 1,Ny-1
      do i = 1,Nx-1
        ! element counter
        ne = ne + 1

        ! point index
        !ixyz = i + (j-1) * (N) + (k-1) * (N * N)

        ! element corners
        n1 = i     + (j-1) * (Nx) + (k-1) * (Nx * Ny)  ! i  ,j  ,k
        n2 = (i+1) + (j-1) * (Nx) + (k-1) * (Nx * Ny)  ! i+1,j  ,k
        n3 = (i+1) + (j  ) * (Nx) + (k-1) * (Nx * Ny)  ! i+1,j+1,k
        n4 = i     + (j  ) * (Nx) + (k-1) * (Nx * Ny)  ! i  ,j+1,k
        n5 = i     + (j-1) * (Nx) + (k  ) * (Nx * Ny)  ! i  ,j  ,k+1
        n6 = (i+1) + (j-1) * (Nx) + (k  ) * (Nx * Ny)  ! i+1,j  ,k+1
        n7 = (i+1) + (j  ) * (Nx) + (k  ) * (Nx * Ny)  ! i+1,j+1,k+1
        n8 = i     + (j  ) * (Nx) + (k  ) * (Nx * Ny)  ! i  ,j+1,k+1

        ! note: indices for VTK start at 0
        total_dat_con(1,ne) = n1 - 1
        total_dat_con(2,ne) = n2 - 1
        total_dat_con(3,ne) = n3 - 1
        total_dat_con(4,ne) = n4 - 1
        total_dat_con(5,ne) = n5 - 1
        total_dat_con(6,ne) = n6 - 1
        total_dat_con(7,ne) = n7 - 1
        total_dat_con(8,ne) = n8 - 1
      enddo
    enddo
  enddo

  ! VTU binary format
  ! stores arrays in databases folder
  call create_name_database(prname,myrank,LOCAL_PATH)
  mesh_file = trim(prname) // 'perturbation_grid.vtu'
  var_name = 'val'
  call write_VTU_movie_data_binary(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  VTU file written: ',trim(mesh_file)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! free memory
  deallocate(total_dat,total_dat_xyz,total_dat_con)

  end subroutine plot_grid_data


!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_scattering_add_perturbations(imaterial_id,xmesh,ymesh,zmesh, &
                                                rho,vp,vs, &
                                                c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                                ANISOTROPY)

  ! overwrites values with updated model values (from iteration step) here, given at all GLL points

  use model_scattering_par

  implicit none

  integer,intent(in) :: imaterial_id
  double precision,intent(in) :: xmesh,ymesh,zmesh

  real(kind=CUSTOM_REAL), intent(inout) :: vp,vs,rho
  real(kind=CUSTOM_REAL), intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  logical, intent(in) :: ANISOTROPY

  ! local parameters
  double precision :: pert_val
  integer :: ix,iy,iz,Nx,Ny,Nz
  ! positioning
  double precision :: offset_x,offset_y,offset_z
  double precision :: spac_x,spac_y,spac_z
  double precision :: gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision :: val1,val2,val3,val4,val5,val6,val7,val8

  ! checks if anything to do
  if (num_target_materials == 0) return

  ! checks if we apply perturbation (to target materials only; target id 0 means to apply to all materials)
  if (.not. (any(target_material_ids(:) == imaterial_id) .or. USE_ALL_MATERIALS)) return

  ! number of grid points (in x,y,z direction)
  Nx = pert_Nx
  Ny = pert_Ny
  Nz = pert_Nz

  ! checks if anything to do in this slice
  if (Nx == 0 .or. Ny == 0 .or. Nz == 0) return

  ! determine spacing and cell for linear interpolation
  offset_x = (xmesh - grid_xmin_proc) / grid_normalization_factor     ! position offset in local mesh (normalized)
  offset_y = (ymesh - grid_ymin_proc) / grid_normalization_factor
  offset_z = (zmesh - grid_zmin_proc) / grid_normalization_factor

  spac_x = offset_x / grid_dx
  spac_y = offset_y / grid_dx
  spac_z = offset_z / grid_dx

  ix = int(spac_x)
  iy = int(spac_y)
  iz = int(spac_z)

  gamma_interp_x = spac_x - dble(ix)
  gamma_interp_y = spac_y - dble(iy)
  gamma_interp_z = spac_z - dble(iz)

  ! shifting indices to start at 1
  ix = ix + 1
  iy = iy + 1
  iz = iz + 1

  ! suppress edge effects for points outside of the model SPOSTARE DOPO
  if (ix < 1) then
     ix = 1
     gamma_interp_x = 0.d0
  endif
  if (ix > Nx-1) then
     ix = Nx-1
     gamma_interp_x = 1.d0
  endif

  if (iy < 1) then
     iy = 1
     gamma_interp_y = 0.d0
  endif
  if (iy > Ny-1) then
     iy = Ny-1
     gamma_interp_y = 1.d0
  endif

  if (iz < 1) then
     iz = 1
     gamma_interp_z = 0.d0
  endif
  if (iz > Nz-1) then
     iz = Nz-1
     gamma_interp_z = 1.d0
  endif

  ! checks
  if (ix < 0 .or. iy < 0 .or. iz < 0) then
    print *,'Error: scattering position has invalid index in local perturbation grid:'
    print *,'  rank        : ',myrank
    print *,'  corner index: ',ix,iy,iz
    print *,'  location    : ',sngl(xmesh),sngl(ymesh),sngl(zmesh)
    print *,'  origin      : x/y/z = ',sngl(grid_origin_x),'/',sngl(grid_origin_y),'/',sngl(grid_origin_z)
    print *,'  pertubation grid: min x/y/z  = ',sngl(grid_xmin_proc),'/',sngl(grid_ymin_proc),'/',sngl(grid_zmin_proc)
    print *,'                    max x/y/z  = ',sngl(grid_xmax_proc),'/',sngl(grid_ymax_proc),'/',sngl(grid_zmax_proc)
    print *,'                    Nx/Ny/Nz  = ',Nx,'/',Ny,'/',Nz
    call exit_MPI(myrank,'Error corner index in model_scattering routine')
  endif

  ! perturbation values
  val1 = perturbation_grid(ix  ,iy  ,iz  )
  val2 = perturbation_grid(ix+1,iy  ,iz  )
  val3 = perturbation_grid(ix+1,iy+1,iz  )
  val4 = perturbation_grid(ix  ,iy+1,iz  )
  val5 = perturbation_grid(ix  ,iy  ,iz+1)
  val6 = perturbation_grid(ix+1,iy  ,iz+1)
  val7 = perturbation_grid(ix+1,iy+1,iz+1)
  val8 = perturbation_grid(ix  ,iy+1,iz+1)

  ! interpolation rule
  ! use trilinear interpolation in cell to define perturbation value
  pert_val =  &
       val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val5 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val6 * gamma_interp_x        * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val7 * gamma_interp_x        * gamma_interp_y        * gamma_interp_z + &
       val8 * (1.d0-gamma_interp_x) * gamma_interp_y        * gamma_interp_z

  !debug
  !if (myrank == 0) &
  !  print *,'debug: loc ',sngl(xmesh),sngl(ymesh),sngl(zmesh),'ixyz',ix,iy,iz, &
  !          'perturbation ',pert_val,'corners',val1,val2,val3,val4,val5,val6,val7,val8
  !debug shared node location
  !if (abs(xmesh-295257) < 1.d1 .and. abs(ymesh-(-152295.d0)) < 1.d1 .and. abs(zmesh-(-1346.64)) < 1.d1) then
  !  print *,'debug: rank ',myrank,'x/y/z',xmesh,ymesh,zmesh,'ix/iy/iz',ix,iy,iz,'pert',pert_val, &
  !        'spac',spac_x,spac_y,spac_z,'gamma',gamma_interp_x,gamma_interp_y,gamma_interp_z, &
  !        'val',val1,val2,val3,val4,val5,val6,val7,val8
  !endif

  ! apply perturbation to model parameters
  call apply_perturbation(pert_val)

contains

  subroutine apply_perturbation(pert_val)
  implicit none
  double precision, intent(in) :: pert_val
  ! local parameters
  real(kind=CUSTOM_REAL) :: scaling

  ! adds perturbation
  scaling = real(1.d0 + pert_val,kind=CUSTOM_REAL)

  ! apply perturbation to model
  ! define new elastic parameters in the model
  rho = scaling * rho
  vp  = scaling * vp
  vs  = scaling * vs
  ! anisotropy
  if (ANISOTROPY) then
    c11 = scaling * c11
    c12 = scaling * c12
    c13 = scaling * c13
    c14 = scaling * c14
    c15 = scaling * c15
    c16 = scaling * c16
    c22 = scaling * c22
    c23 = scaling * c23
    c24 = scaling * c24
    c25 = scaling * c25
    c26 = scaling * c26
    c33 = scaling * c33
    c34 = scaling * c34
    c35 = scaling * c35
    c36 = scaling * c36
    c44 = scaling * c44
    c45 = scaling * c45
    c46 = scaling * c46
    c55 = scaling * c55
    c56 = scaling * c56
    c66 = scaling * c66
  endif

  end subroutine apply_perturbation

  end subroutine model_scattering_add_perturbations

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_FFTW

  subroutine FFT_3D_FFTW(N, array, zign)

! uses FFTW version 3 library
!
! for compilation, modify Makefile:
! ```
!   FLAGS_CHECK = -DUSE_FFTW -O3 ..
!   ..
!   MPILIBS = -lfftw3 -L/opt/homebrew/lib ..
!   ..
! ```

  use, intrinsic :: iso_c_binding
  use constants, only: CUSTOM_REAL

  implicit none

  include 'fftw3.f03'

  integer, parameter :: CUSTOM_CMPLX = 8

  integer, intent(in) :: N
  complex(kind=CUSTOM_CMPLX), dimension(N,N,N), intent(inout) :: array
  real(kind=CUSTOM_REAL), intent(in) :: zign  ! +1 for forward, -1 for inverse

  ! local parameters
  integer :: ier
  complex(kind=CUSTOM_CMPLX), dimension(:,:,:), allocatable :: array_in
  ! FFTW plan
  type(C_PTR) :: plan
  integer(C_INT) :: fftw_direction

  ! user info
  print *,'***'
  print *,'using FFTW library calls'
  print *,'***'

  ! Determine direction
  if (zign > 0.0_CUSTOM_REAL) then
    fftw_direction = FFTW_FORWARD   ! Forward FFT
  else
    fftw_direction = FFTW_BACKWARD  ! Inverse FFT
  endif

  ! Copy to avoid aliasing warning
  allocate(array_in(N,N,N),stat=ier)
  if (ier /= 0) stop 'Error allocating array_in'
  array_in(:,:,:) = array(:,:,:)

  ! Create plan for 3D complex-to-complex FFT
  ! FFTW_ESTIMATE: quick planning, reasonable performance
  ! Alternative: FFTW_MEASURE for better performance (takes longer to plan)
  plan = fftw_plan_dft_3d(N, N, N, array_in, array, fftw_direction, FFTW_ESTIMATE)

  ! Execute the FFT
  call fftw_execute_dft(plan, array_in, array)

  ! Cleanup
  call fftw_destroy_plan(plan)

  ! free temporary array
  deallocate(array_in)

  ! Note: FFTW inverse transform is unnormalized
  ! Divide by N^3 for inverse transform if needed
  if (zign < 0.0_CUSTOM_REAL) then
    array = array / real(N*N*N, kind=CUSTOM_CMPLX)
  endif

  end subroutine FFT_3D_FFTW

#endif
