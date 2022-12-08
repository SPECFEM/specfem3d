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

! Base module for kinematic and dynamic fault solvers
!
! Authors:
! Percy Galvez, Surendra Somala, Jean-Paul Ampuero

module fault_solver_common

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IMAIN,IN_DATA_FILES

  implicit none

  type fault_type
    integer :: nspec=0, nglob=0
    real(kind=CUSTOM_REAL), dimension(:,:),   pointer :: T => null(),V => null(),D => null(),coord => null()
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: R => null()
    real(kind=CUSTOM_REAL), dimension(:),     pointer :: B => null(),invM1 => null(),invM2 => null(),Z => null()
    real(kind=CUSTOM_REAL) :: dt
    integer, dimension(:), pointer :: ibulk1 => null(), ibulk2 => null()  ! global nodes id
    integer, dimension(:), pointer :: ispec1 => null(), ispec2 => null()  ! fault elements id
  end type fault_type

  ! outputs(dyn) /inputs (kind) at selected times for all fault nodes:
  ! strength, state, slip, slip velocity, fault stresses, rupture time, process zone time
  ! rupture time = first time when slip velocity = threshold V_RUPT (defined below)
  ! process zone time = first time when slip = Dc
  type dataXZ_type
    real(kind=CUSTOM_REAL), dimension(:), pointer :: stg => null(), sta => null(), d1 => null(), d2 => null(), &
                                                     v1 => null(), v2 => null(), &
                                                     t1 => null(), t2 => null(), t3 => null(), tRUP => null(), tPZ => null()
    real(kind=CUSTOM_REAL), dimension(:), pointer :: xcoord => null(), ycoord => null(), zcoord => null()
    integer                                       :: npoin=0
  end type dataXZ_type

  type swf_type
    integer :: kind
    logical :: healing = .false.
    real(kind=CUSTOM_REAL), dimension(:), pointer :: Dc => null(), mus => null(), mud => null(), &
                                                     theta => null(), T => null(), C => null()
  end type swf_type

  type twf_type
    real(kind=CUSTOM_REAL) ::  nuc_x, nuc_y, nuc_z, nuc_r, nuc_t0, nuc_v
  end type twf_type


  type rsf_type
    integer :: StateLaw = 1 ! 1=ageing law, 2=slip law
    real(kind=CUSTOM_REAL), dimension(:), pointer :: V0 => null(), f0 => null(), L => null(), &
                                                     V_init => null(), &
                                                     a => null(), b => null(), theta => null(), &
                                                     T => null(), C => null(), &
                                                     fw => null(), Vw => null()
  end type rsf_type

 ! outputs on selected fault nodes at every time step:
  type dataT_type
    integer :: npoin=0, ndat=0, nt=0
    real(kind=CUSTOM_REAL) :: dt
    real(kind=CUSTOM_REAL) :: element_size
    integer, dimension(:), pointer :: iglob => null()   ! on-fault global index of output nodes
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: dat => null()
    character(len=MAX_STRING_LEN), dimension(:), pointer :: name => null(),longFieldNames => null()
    character(len=MAX_STRING_LEN) :: shortFieldNames
  end type dataT_type

  type, extends (fault_type) :: bc_dynandkinflt_type
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: T0 => null()
    real(kind=CUSTOM_REAL), dimension(:),   pointer :: mu => null(), Fload => null()
    integer, dimension(:),   pointer :: npoin_perproc => null(), poin_offset => null()
    type(dataT_type)        :: dataT
    type(dataXZ_type)       :: dataXZ,dataXZ_all
    type(swf_type), pointer :: swf => null()
    type(rsf_type), pointer :: rsf => null()
    type(twf_type), pointer :: twf => null()
    logical                 :: allow_opening = .false. ! default : do not allow opening

!! DK DK added this in order to be able to use the type for both dynamic and kinematic faults
    real(kind=CUSTOM_REAL) :: kin_dt
    integer :: kin_it
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: v_kin_t1,v_kin_t2
  end type bc_dynandkinflt_type

  ! fault array
  type(bc_dynandkinflt_type),dimension(:), pointer :: faults => null()
  ! Number of faults
  integer, save                                    :: Nfaults = 0

  ! Kelvin-Voigt damping
  real(kind=CUSTOM_REAL), allocatable, save :: Kelvin_Voigt_eta(:)
  logical :: USE_KELVIN_VOIGT_DAMPING = .false.

  public :: fault_type, &
            initialize_fault, get_jump, get_weighted_jump, rotate, add_BT, &
            dataT_type, init_dataT, store_dataT, SCEC_write_dataT, &
            Kelvin_Voigt_eta, USE_KELVIN_VOIGT_DAMPING, &
            fault_check_mesh_resolution

contains

!---------------------------------------------------------------------

  subroutine initialize_fault (bc,IIN_BIN)

  use constants, only: PARALLEL_FAULT,NDIM,NGLLSQUARE

  use specfem_par, only: NPROC,DT,NGLOB_AB, &
    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
    my_neighbors_ext_mesh

  use specfem_par_elastic, only: rmassx

  implicit none
!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(inout) :: bc
  integer, intent(in) :: IIN_BIN

  real(kind=CUSTOM_REAL) :: tmp_vec(3,NGLOB_AB)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: jacobian2Dw
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nxyz
  integer, dimension(:,:), allocatable :: ibool1
  integer :: ij,k,e,ier

  ! number of elements and global points on fault
  read(IIN_BIN) bc%nspec,bc%nglob

  if (bc%nspec > 0) then
    ! checks nglob
    if (bc%nglob < 1) then
      print *,'Error: invalid fault surface with nspec ',bc%nspec,' and nglob ',bc%nglob
      stop 'Invalid fault with zero global points'
    endif

    ! array allocations
    ! fault points (upper surface)
    allocate(bc%ibulk1(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2159')
    ! fault points (lower surface)
    allocate(bc%ibulk2(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2160')
    ! rotation matrix (rotates from x/y/z <-> strike/dip/normal)
    allocate(bc%R(NDIM,NDIM,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2161')
    ! fault point coordinates
    allocate(bc%coord(NDIM,(bc%nglob)),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2162')
    ! inverse of mass matrix (upper surface)
    allocate(bc%invM1(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2163')
    ! inverse of mass matrix (lower surface)
    allocate(bc%invM2(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2164')
    ! fault boundary matrix
    allocate(bc%B(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2165')
    ! fault impedance
    allocate(bc%Z(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2166')
    ! fault element iglob
    allocate(ibool1(NGLLSQUARE,bc%nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2167')
    ! fault element normal
    allocate(normal(NDIM,NGLLSQUARE,bc%nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2168')
    ! fault element jacobian
    allocate(jacobian2Dw(NGLLSQUARE,bc%nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2169')
    ! fault elements ispec (touching upper surface)
    allocate(bc%ispec1(bc%nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2170')
    ! fault elements ispec (touching lower surface)
    allocate(bc%ispec2(bc%nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2171')

    read(IIN_BIN) ibool1
    read(IIN_BIN) jacobian2Dw
    read(IIN_BIN) normal
    read(IIN_BIN) bc%ibulk1
    read(IIN_BIN) bc%ibulk2
    read(IIN_BIN) bc%coord(1,:)
    read(IIN_BIN) bc%coord(2,:)
    read(IIN_BIN) bc%coord(3,:)

    ! optional, checks if ispec arrays for elements on fault are saved
    ! (older fault_db.bin files won't have these records)
    read(IIN_BIN,iostat=ier) bc%ispec1
    ! dummy value if not stored
    if (ier /= 0) bc%ispec1(:) = 0

    read(IIN_BIN,iostat=ier) bc%ispec2
    ! dummy value if not stored
    if (ier /= 0) bc%ispec2(:) = 0

    ! done reading faults_db.bin file, no more records in read(IIN_BIN).. from here on

    ! sets simulation time step size
    bc%dt = real(DT,kind=CUSTOM_REAL)

    ! normal vector
    allocate(nxyz(NDIM,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2170')

    nxyz(:,:) = 0.0_CUSTOM_REAL
    bc%B(:) = 0.0_CUSTOM_REAL

    do e = 1,bc%nspec
      do ij = 1,NGLLSQUARE
        k = ibool1(ij,e)
        nxyz(:,k) = nxyz(:,k) + normal(:,ij,e)
        bc%B(k) = bc%B(k) + jacobian2Dw(ij,e)
      enddo
    enddo
  else
    ! dummy allocations (for subroutine arguments)
    allocate(bc%coord(NDIM,1), &
             bc%B(1), &
             bc%Z(1), &
             bc%invM1(1), &
             bc%invM2(1), &
             bc%ibulk1(1), &
             bc%ibulk2(1), &
             bc%R(NDIM,NDIM,1))
  endif

  ! fault parallelization across multiple MPI processes
  if (PARALLEL_FAULT) then
    ! assembles fault boundary matrix B
    tmp_vec(:,:) = 0.0_CUSTOM_REAL
    if (bc%nspec > 0) tmp_vec(1,bc%ibulk1(:)) = bc%B(:)

    ! assembles with other MPI processes
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,tmp_vec, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    if (bc%nspec > 0) bc%B(:) = tmp_vec(1,bc%ibulk1(:))

    ! assembles fault normal n
    tmp_vec(:,:) = 0.0_CUSTOM_REAL
    if (bc%nspec > 0) tmp_vec(:,bc%ibulk1(:)) = nxyz(:,:)

    ! assembles with other MPI processes
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,tmp_vec, &
                                     num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     my_neighbors_ext_mesh)

    if (bc%nspec > 0) nxyz(:,:) = tmp_vec(:,bc%ibulk1(:))
  endif

  if (bc%nspec > 0) then
    ! normalizes fault normal vector
    call normalize_3d_vector(nxyz)
    ! sets up rotation matrix R
    call compute_R(bc%R,bc%nglob,nxyz)

    ! inverse mass matrix
    !SURENDRA : WARNING! Assuming rmassx=rmassy=rmassz
    ! Needed in dA_Free = -K2*d2/M2 + K1*d1/M1
    bc%invM1(:) = rmassx(bc%ibulk1(:))
    bc%invM2(:) = rmassx(bc%ibulk2(:))

    ! Fault impedance
    !   Z in :  Trac = T_Stick - Z*dV
    !   Z = 1/( B1/M1 + B2/M2 ) / (0.5*dt)
    ! T_stick = Z*Vfree traction as if the fault was stuck (no displ discontinuity)
    ! NOTE: same Bi on both sides, see note above

    bc%Z(:) = 1.0_CUSTOM_REAL/(0.5_CUSTOM_REAL * bc%dt * bc%B(:) *( bc%invM1(:) + bc%invM2(:) ))

    ! WARNING: In non-split nodes at fault edges M is assembled across the fault.
    ! hence invM1+invM2=2/(M1+M2) instead of 1/M1+1/M2
    ! In a symmetric mesh (M1=M2) Z will be twice its intended value
  endif

  end subroutine initialize_fault

!---------------------------------------------------------------------

  subroutine normalize_3d_vector(v)

  implicit none
  real(kind=CUSTOM_REAL), intent(inout) :: v(:,:)

  ! local parameters
  real(kind=CUSTOM_REAL) :: norm
  integer :: k

 ! assume size(v) = [3,N]
  do k = 1,size(v,2)
    ! vector length
    norm = sqrt( v(1,k)*v(1,k) + v(2,k)*v(2,k) + v(3,k)*v(3,k) )

    ! normalizes vector
    ! checks norm to avoid division by zero
    if (norm > 1.e-24) then
      v(:,k) = v(:,k) / norm
    else
      v(:,k) = 0.0_CUSTOM_REAL
    endif
  enddo

  end subroutine normalize_3d_vector

!---------------------------------------------------------------------

! Percy: define fault directions according to SCEC conventions
! Fault coordinates (s,d,n) = (1,2,3)
!   s = strike , d = dip , n = normal
!   1 = strike , 2 = dip , 3 = normal
! with dip pointing downwards
!
  subroutine compute_R(R,nglob,n)

  implicit none
  integer :: nglob
  real(kind=CUSTOM_REAL), intent(out) :: R(3,3,nglob)
  real(kind=CUSTOM_REAL), intent(in) :: n(3,nglob)

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: s,d

  ! strike direction
  s(1,:) =  n(2,:)   ! sx = ny
  s(2,:) = -n(1,:)   ! sy =-nx
  s(3,:) = 0.0_CUSTOM_REAL

  ! set the along strike direction when the fault is a horizontal plane.
  where(abs(s(1,:))+abs(s(2,:)) < 1e-6)
      s(1,:) = 1.0_CUSTOM_REAL
      s(2,:) = 0.0_CUSTOM_REAL
  endwhere

  call normalize_3d_vector(s)

  ! dip direction
  d(1,:) = -s(2,:)*n(3,:)                 ! dx = -sy*nz
  d(2,:) =  s(1,:)*n(3,:)                 ! dy = sx*nz
  d(3,:) =  s(2,:)*n(1,:) - s(1,:)*n(2,:) ! dz = sy*nx-ny*sx

  call normalize_3d_vector(d)

  ! dz is always dipwards (negative), because
  ! (nx*sy-ny*sx) = -(nx^2+ny^2)/sqrt(nx^2+ny^2)
  !               = -sqrt(nx^2+ny^2) < 0

  R(1,:,:) = s(:,:)  ! strike
  R(2,:,:) = d(:,:)  ! dip
  R(3,:,:) = n(:,:)  ! normal

  end subroutine compute_R


!===============================================================

  function get_jump (bc,v) result(dv)

  implicit none
!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: v(:,:)
  real(kind=CUSTOM_REAL) :: dv(3,bc%nglob)

  ! difference between side 2 and side 1 of fault nodes. dv
  dv(1,:) = v(1,bc%ibulk2(:)) - v(1,bc%ibulk1(:))
  dv(2,:) = v(2,bc%ibulk2(:)) - v(2,bc%ibulk1(:))
  dv(3,:) = v(3,bc%ibulk2(:)) - v(3,bc%ibulk1(:))

  end function get_jump

!---------------------------------------------------------------------

  function get_weighted_jump (bc,f) result(da)

  implicit none
!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: f(:,:)

  real(kind=CUSTOM_REAL) :: da(3,bc%nglob)

  ! difference between side 2 and side 1 of fault nodes. M-1 * F
  da(1,:) = bc%invM2(:)*f(1,bc%ibulk2(:)) - bc%invM1(:)*f(1,bc%ibulk1(:))
  da(2,:) = bc%invM2(:)*f(2,bc%ibulk2(:)) - bc%invM1(:)*f(2,bc%ibulk1(:))
  da(3,:) = bc%invM2(:)*f(3,bc%ibulk2(:)) - bc%invM1(:)*f(3,bc%ibulk1(:))

  ! NOTE: In non-split nodes at fault edges M and f are assembled across the fault.
  ! Hence, f1=f2, invM1=invM2=1/(M1+M2) instead of invMi=1/Mi, and da=0.

  end function get_weighted_jump

!----------------------------------------------------------------------

  function rotate(bc,v,fb) result(vr)

  implicit none
!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: v(3,bc%nglob)
  integer, intent(in) :: fb

  real(kind=CUSTOM_REAL) :: vr(3,bc%nglob)

  ! Percy, tangential direction Vt, equation 7 of Pablo's notes in agreement with SPECFEM3D

  if (fb == 1) then
    ! forward rotation: (v * R)_i = v_j R_ij
    vr(1,:) = v(1,:)*bc%R(1,1,:)+v(2,:)*bc%R(1,2,:)+v(3,:)*bc%R(1,3,:) ! vs  strike
    vr(2,:) = v(1,:)*bc%R(2,1,:)+v(2,:)*bc%R(2,2,:)+v(3,:)*bc%R(2,3,:) ! vd  dip
    vr(3,:) = v(1,:)*bc%R(3,1,:)+v(2,:)*bc%R(3,2,:)+v(3,:)*bc%R(3,3,:) ! vn  normal direction
  else
    ! backward rotation: (v * R^T)_i = v_j R_ji
    vr(1,:) = v(1,:)*bc%R(1,1,:)+v(2,:)*bc%R(2,1,:)+v(3,:)*bc%R(3,1,:)  !vx
    vr(2,:) = v(1,:)*bc%R(1,2,:)+v(2,:)*bc%R(2,2,:)+v(3,:)*bc%R(3,2,:)  !vy
    vr(3,:) = v(1,:)*bc%R(1,3,:)+v(2,:)*bc%R(2,3,:)+v(3,:)*bc%R(3,3,:)  !vz
  endif

  end function rotate

!----------------------------------------------------------------------

  subroutine add_BT(bc,MxA,T)

  implicit none
!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(inout) :: MxA(:,:)
  real(kind=CUSTOM_REAL), dimension(3,bc%nglob), intent(in) :: T

  MxA(1,bc%ibulk1(:)) = MxA(1,bc%ibulk1(:)) + bc%B(:)*T(1,:)
  MxA(2,bc%ibulk1(:)) = MxA(2,bc%ibulk1(:)) + bc%B(:)*T(2,:)
  MxA(3,bc%ibulk1(:)) = MxA(3,bc%ibulk1(:)) + bc%B(:)*T(3,:)

  MxA(1,bc%ibulk2(:)) = MxA(1,bc%ibulk2(:)) - bc%B(:)*T(1,:)
  MxA(2,bc%ibulk2(:)) = MxA(2,bc%ibulk2(:)) - bc%B(:)*T(2,:)
  MxA(3,bc%ibulk2(:)) = MxA(3,bc%ibulk2(:)) - bc%B(:)*T(3,:)

  end subroutine add_BT


!===============================================================
! dataT outputs

  subroutine init_dataT(dataT,coord,nglob,NT,DT,ndat,iflt)

  use constants, only: PARALLEL_FAULT
  use specfem_par, only: NPROC,myrank

  implicit none
  integer, intent(in) :: nglob,NT,iflt,ndat
  real(kind=CUSTOM_REAL), intent(in) :: coord(3,nglob),DT

!! DK DK use type(dataT_type) instead of class(dataT_type) for compatibility with some current compilers
  type(dataT_type), intent(inout) :: dataT

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: dist_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dist_loc
  integer, dimension(:,:), allocatable :: iglob_all
  integer, dimension(:), allocatable :: iproc,iglob_tmp,glob_indx
  real(kind=CUSTOM_REAL) :: xtarget,ytarget,ztarget,dist,distkeep
  integer :: i, iglob , IIN, ier, jflt, np, k
  character(len=MAX_STRING_LEN) :: tmpname
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: name_tmp
  integer :: ipoin, ipoin_local, npoin_local

  !  1. read fault output coordinates from user file,
  !  2. define iglob: the fault global index of the node nearest to user
  !     requested coordinate

  IIN = 251 ! WARNING: not safe, should check that unit is not aleady opened

  ! initializes dataT
  dataT%ndat = ndat
  dataT%nt = NT
  dataT%dt = DT
  dataT%element_size = 0.0

  ! count the number of output points on the current fault (#iflt)
  open(IIN,file=IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FAULT_STATIONS',status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file ',IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FAULT_STATIONS'
    if (myrank == 0) write(IMAIN,*) 'Fatal error opening FAULT_STATIONS file. Abort.'
    stop 'Error opening file FAULT_STATIONS'
  endif

  ! number of points
  read(IIN,*) np

  ! counts fault points on specified fault
  dataT%npoin = 0
  do i = 1,np
    ! format : #x  #y  #z  #name  #fault_id
    ! example: -4500.0  0.0  0.0  faultst-045dp000 1
    read(IIN,*) xtarget,ytarget,ztarget,tmpname,jflt
    ! only points on this fault
    if (jflt == iflt) dataT%npoin = dataT%npoin + 1
  enddo
  close(IIN)

  ! allocates fault point arrays
  if (dataT%npoin > 0) then
    allocate(dataT%iglob(dataT%npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2171')
    allocate(dataT%name(dataT%npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2172')
    ! for parallel fault
    allocate(dist_loc(dataT%npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2173')
  else
    ! dummy arrays
    allocate(dataT%iglob(1), &
             dataT%name(1), &
             dist_loc(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2174')
  endif

  ! initializes arrays
  dataT%iglob(:) = 0
  dataT%name(:) = ''
  dist_loc(:) = huge(distkeep)

  ! checks if anything left to do
  if (dataT%npoin == 0) return

  ! opens in fault stations
  open(IIN,file=IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FAULT_STATIONS',status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file ',IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FAULT_STATIONS'
    stop 'Error opening file FAULT_STATIONS'
  endif

  ! number of points
  read(IIN,*) np

  ! reads in fault point positions
  k = 0
  do i = 1,np
    ! format : #x  #y  #z  #name  #fault_id
    ! example: -4500.0  0.0  0.0  faultst-045dp000 1
    read(IIN,*) xtarget,ytarget,ztarget,tmpname,jflt

    ! only points on this fault
    if (jflt /= iflt) cycle

    k = k+1
    dataT%name(k) = tmpname

    ! search nearest node
    distkeep = huge(distkeep)
    do iglob = 1,nglob
      ! todo: if this takes too long - this could also be done without taking sqrt(..)
      dist = sqrt( (coord(1,iglob)-xtarget)**2 &
                 + (coord(2,iglob)-ytarget)**2 &
                 + (coord(3,iglob)-ztarget)**2)
      if (dist < distkeep) then
        distkeep = dist
        dataT%iglob(k) = iglob
      endif
    enddo
    dist_loc(k) = distkeep
  enddo
  close(IIN)

  if (PARALLEL_FAULT) then
    ! For each output point, find the processor that contains the nearest node
    ! temporary arrays
    allocate(iproc(dataT%npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2174')
    allocate(iglob_all(dataT%npoin,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2175')
    allocate(dist_all(dataT%npoin,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2176')

    call gather_all_i(dataT%iglob,dataT%npoin,iglob_all,dataT%npoin,NPROC)
    call gather_all_cr(dist_loc,dataT%npoin,dist_all,dataT%npoin,NPROC)

    if (myrank == 0) then
      ! NOTE: output points lying at an interface between procs are assigned to a unique proc
      iproc = minloc(dist_all,2) - 1
      do ipoin = 1,dataT%npoin
        dataT%iglob(ipoin) = iglob_all(ipoin,iproc(ipoin))
      enddo
    endif
    call bcast_all_i(iproc,dataT%npoin)
    call bcast_all_i(dataT%iglob,dataT%npoin)

    ! Number of output points contained in the current processor
    npoin_local = count( iproc == myrank )

    if (npoin_local > 0) then
      ! Make a list of output points contained in the current processor
      allocate(glob_indx(npoin_local),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2177')

      ipoin_local = 0
      do ipoin = 1,dataT%npoin
        if (myrank == iproc(ipoin)) then
          ipoin_local = ipoin_local + 1
          glob_indx(ipoin_local) = ipoin
        endif
      enddo

      ! Consolidate the output information (remove output points outside current proc)
      allocate(iglob_tmp(dataT%npoin),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2178')
      allocate(name_tmp(dataT%npoin),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2179')

      iglob_tmp(:) = dataT%iglob(:)
      name_tmp(:) = dataT%name(:)

      deallocate(dataT%iglob)
      deallocate(dataT%name)

      dataT%npoin = npoin_local

      allocate(dataT%iglob(dataT%npoin),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2180')
      allocate(dataT%name(dataT%npoin),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2181')

      dataT%iglob = iglob_tmp(glob_indx)
      dataT%name = name_tmp(glob_indx)

      deallocate(glob_indx,iglob_tmp,name_tmp)
    else
      ! no local points
      dataT%npoin = 0
    endif

    ! free temporary arrays
    deallocate(iproc,iglob_all,dist_all)
  else
    ! fault in single slice
    ! no parallel fault means that all fault element are in a single process, that is, in rank 0 process
    ! checks if this process has any fault points at all
    if (nglob == 0) then
      ! cannot have local station points
      dataT%npoin = 0
    endif
  endif

  !  3. initialize arrays
  if (dataT%npoin > 0) then
    allocate(dataT%dat(dataT%ndat,dataT%npoin,dataT%nt),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2182')
    dataT%dat(:,:,:) = 0.0_CUSTOM_REAL

    allocate(dataT%longFieldNames(dataT%ndat),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2183')
    dataT%longFieldNames(1) = "horizontal right-lateral slip (m)"
    dataT%longFieldNames(2) = "horizontal right-lateral slip rate (m/s)"
    dataT%longFieldNames(3) = "horizontal right-lateral shear stress (MPa)"
    dataT%longFieldNames(4) = "vertical up-dip slip (m)"
    dataT%longFieldNames(5) = "vertical up-dip slip rate (m/s)"
    dataT%longFieldNames(6) = "vertical up-dip shear stress (MPa)"
    dataT%longFieldNames(7) = "normal stress (MPa)"
    dataT%shortFieldNames = "h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress"
  else
    ! dummy allocations (for subroutine arguments)
    allocate(dataT%dat(1,1,1))
  endif

  end subroutine init_dataT

!---------------------------------------------------------------

  subroutine store_dataT(dataT,D,V,T,itime)

  implicit none
!! DK DK use type() instead of class() for compatibility with some current compilers
  type(dataT_type), intent(inout) :: dataT
  real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: D,V,T
  integer, intent(in) :: itime

  integer :: ipoin,iglob

  do ipoin = 1,dataT%npoin
    iglob = dataT%iglob(ipoin)
    dataT%dat(1,ipoin,itime) = D(1,iglob)                     ! horizontal right-lateral slip (m)
    dataT%dat(2,ipoin,itime) = V(1,iglob)                     ! horizontal right-lateral slip rate (m/s)
    dataT%dat(3,ipoin,itime) = T(1,iglob)/1.0e6_CUSTOM_REAL   ! horizontal right-lateral shear stress (MPa)
    dataT%dat(4,ipoin,itime) = D(2,iglob)                     ! vertical up-dip slip (m)
    dataT%dat(5,ipoin,itime) = V(2,iglob)                     ! vertical up-dip slip rate (m/s)
    dataT%dat(6,ipoin,itime) = T(2,iglob)/1.0e6_CUSTOM_REAL   ! vertical up-dip shear stress (MPa)
    dataT%dat(7,ipoin,itime) = T(3,iglob)/1.0e6_CUSTOM_REAL   ! normal stress (MPa)
  enddo

  end subroutine store_dataT

!------------------------------------------------------------------------

  subroutine SCEC_write_dataT(dataT)

  use constants, only: NGLLX
  use specfem_par, only: OUTPUT_FILES

  implicit none

!! DK DK use type() instead of class() for compatibility with some current compilers
  type(dataT_type), intent(in) :: dataT

  ! local parameters
  integer :: ipoin,k,elem_size
  character(len=10) :: my_fmt
  character(len=24) :: my_fmt_size
  integer, dimension(8) :: time_values
  integer, parameter :: IOUT_SC = 121 !WARNING: not very robust. Could instead look for an available ID

  call date_and_time(VALUES=time_values)

  ! element size in m
  elem_size = int(dataT%element_size)
  if (elem_size == 0) then
    write(my_fmt_size,'(a,i1,a,i1,a)') '(a,i',1,',a,i',int(log10(1.0*NGLLX))+1,',a)'
  else
    write(my_fmt_size,'(a,i1,a,i1,a)') '(a,i',int(log10(1.0*elem_size))+1,',a,i',int(log10(1.0*NGLLX))+1,',a)'
  endif

  write(my_fmt,'(a,i1,a)') '(',dataT%ndat+1,'(E15.7))'

  do ipoin = 1,dataT%npoin
    open(IOUT_SC,file=trim(OUTPUT_FILES)//trim(dataT%name(ipoin))//'.dat',status='replace')

    write(IOUT_SC,'(a)') "# problem=TPV104" ! WARNING: this should be a user input
    write(IOUT_SC,'(a)') "# author=Surendra Nadh Somala" ! WARNING: this should be a user input
    write(IOUT_SC,1000) time_values(2), time_values(3), time_values(1), time_values(5), time_values(6), time_values(7)
    write(IOUT_SC,'(a)') "# code=SPECFEM3D_Cartesian (split nodes)"
    write(IOUT_SC,'(a)') "# code_version=1.1"
    write(IOUT_SC,my_fmt_size) "# element_size=",elem_size," m  (w/ ",NGLLX," GLL nodes)"
    write(IOUT_SC,'(a,E15.7)') "# time_step=",dataT%dt
    write(IOUT_SC,'(a)') "# location=" // trim(dataT%name(ipoin))
    write(IOUT_SC,'(a)') "# Column #1 = Time (s)"

    do k = 1,dataT%ndat
      write(IOUT_SC,1100) k+1,trim(dataT%longFieldNames(k))
    enddo

    write(IOUT_SC,'(a)') "#"
    write(IOUT_SC,'(a)') "# The line below lists the names of the data fields:"
    write(IOUT_SC,'(a)') "# t " // trim(dataT%shortFieldNames)
    write(IOUT_SC,'(a)') "#"

    do k = 1,dataT%nt
      write(IOUT_SC,my_fmt) k*dataT%dt, dataT%dat(:,ipoin,k)
    enddo

    close(IOUT_SC)
  enddo

1000 format ( '# Date = ', i2.2, '/', i2.2, '/', i4.4, '; time = ',i2.2, ':', i2.2, ':', i2.2 )
1100 format ( '# Column #', i1, ' = ',a )

  end subroutine SCEC_write_dataT

!---------------------------------------------------------------------

  subroutine fault_check_mesh_resolution()

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,CUSTOM_REAL,HUGEVAL,COURANT_SUGGESTED,myrank
  use specfem_par, only: NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore,kappastore,mustore,rhostore,DT
  use specfem_par_elastic, only: ispec_is_elastic

  implicit none

  ! local parameters
  type(bc_dynandkinflt_type),pointer :: bc
  real(kind=CUSTOM_REAL) :: eta_max,eta_max_glob
  real(kind=CUSTOM_REAL) :: h_fault,csp_max,DIM
  real(kind=CUSTOM_REAL) :: dtc_fault,dtc_kv
  integer :: NGLL,ispec,iside,elem,ifault
  logical :: USE_KELVIN_VOIGT_DAMPING_all

  ! elements info
  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax,vpmin_glob,vpmax_glob
  real(kind=CUSTOM_REAL) :: poissonmin,poissonmax
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: avg_distance,dt_suggested
  real(kind=CUSTOM_REAL) :: eta_suggested_min,eta_suggested_max
  logical :: has_vs_zero

  ! tabulated values of critical frequency (non-dimensional, 1D)
  ! (taken from utils/critical_timestep.m)
  real(kind=CUSTOM_REAL),dimension(19),parameter :: &
    Omega_max = (/2.0000000e+00, 4.8989795e+00, 8.6203822e+00, 1.3540623e+01, 1.9797952e+01, 2.7378050e+01, &
                  3.6256848e+01, 4.6421894e+01, 5.7867306e+01, 7.0590158e+01, 8.4588883e+01, 9.9862585e+01, &
                  1.1641072e+02, 1.3423295e+02, 1.5332903e+02, 1.7369883e+02, 1.9534221e+02, 2.1825912e+02, &
                  2.4244948e+02 /)

  ! stability factor for leapfrog timescheme
  real(kind=CUSTOM_REAL),parameter :: C = 2.0_CUSTOM_REAL

  ! gets fault element size and maximum Vp velocity
  vpmin_glob = HUGEVAL
  vpmax_glob = -HUGEVAL

  elemsize_min_glob = HUGEVAL
  elemsize_max_glob = -HUGEVAL

  ! loop over all faults
  do ifault = 1,Nfaults
    bc => faults(ifault)
    ! loops over fault elements
    do elem = 1,bc%nspec
      ! loops over both sides of fault
      do iside = 1,2
        ! gets element id
        if (iside == 1) then
          ispec = bc%ispec1(elem)
        else
          ispec = bc%ispec2(elem)
        endif

        ! checks if valid entry
        if (ispec < 1) cycle

        ! checks if element is elastic
        if (.not. ispec_is_elastic(ispec)) then
          print *,'Error fault elements must be elastic: element ',ispec,' is in different acoustic/poroelastic domain'
          print *,'Please check mesh setup for dynamic rupture simulation...'
          stop 'Invalid fault element, elements in acoustic/poroelastic domains not supported yet'
        endif

        ! determines minimum/maximum velocities within this element
        ! elastic/acoustic
        call get_vpvs_minmax(vpmin,vpmax,vsmin,vsmax,poissonmin,poissonmax, &
                             ispec,has_vs_zero, &
                             NSPEC_AB,kappastore,mustore,rhostore)

        ! min/max for whole cpu partition
        vpmin_glob = min(vpmin_glob, vpmin)
        vpmax_glob = max(vpmax_glob, vpmax)

        ! computes minimum and maximum size of this grid cell
        call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                                 NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

        elemsize_min_glob = min(elemsize_min_glob, elemsize_min)
        elemsize_max_glob = max(elemsize_max_glob, elemsize_max)
      enddo ! iside
    enddo ! elem
  enddo

  ! main process gathers overall values
  ! Vp velocity
  vpmin = vpmin_glob
  vpmax = vpmax_glob
  call min_all_cr(vpmin,vpmin_glob)
  call max_all_cr(vpmax,vpmax_glob)

  ! element size
  elemsize_min = elemsize_min_glob
  elemsize_max = elemsize_max_glob
  call min_all_cr(elemsize_min,elemsize_min_glob)
  call max_all_cr(elemsize_max,elemsize_max_glob)

  ! sets element size in dataT structure
  ! broadcast to all processes
  call bcast_all_singlecr(elemsize_max_glob)
  ! all faults will have same element size set
  do ifault = 1,Nfaults
    bc => faults(ifault)
    bc%dataT%element_size = elemsize_max_glob
  enddo

  ! gets maximum damping value set by Par_file_faults
  if (USE_KELVIN_VOIGT_DAMPING) then
    eta_max_glob = maxval(Kelvin_Voigt_eta(:))
  else
    eta_max_glob = 0.0_CUSTOM_REAL
  endif
  eta_max = eta_max_glob
  call max_all_cr(eta_max,eta_max_glob)

  ! check damping flag; see if any damping will be used, for user output
  call any_all_l( USE_KELVIN_VOIGT_DAMPING, USE_KELVIN_VOIGT_DAMPING_all )

  ! only main process has valid values
  if (myrank == 0) then
    ! average distance between GLL points within this element
    avg_distance = elemsize_max_glob / ( NGLLX - 1)  ! since NGLLX = NGLLY = NGLLZ
    ! rough estimate of time step
    dt_suggested = COURANT_SUGGESTED * avg_distance / vpmax_glob

    ! critical time step estimation
    ! see: utils/critical_timestep.m
    NGLL = max(NGLLX,NGLLY,NGLLZ)
    ! checks
    if (NGLL < 2 .or. NGLL > 20) stop 'Error critical time step estimation: NGLL must be from 2 to 20'

    ! critical time step,
    ! assumes a cube element
    h_fault = elemsize_min_glob
    csp_max = vpmax_glob
    DIM = NDIM

    ! critical time step (for fault elements)
    dtc_fault = C * h_fault / csp_max / sqrt(DIM) / Omega_max(NGLL-1)

    ! proposed eta value range (see manual about explanation) is between [0.1,0.3] * dt_c
    ! (0.1 to 0.3)$\ensuremath{\times}\mathtt{dtc\_fault}$
    eta_suggested_min = min(0.1 * dtc_fault,0.1 * dt_suggested)
    eta_suggested_max = min(0.3 * dtc_fault,0.3 * dt_suggested)

    ! critical time step with given Kelvin-Voigt damping
    ! (see manual for explanation)
    if (eta_max_glob > 0.0_CUSTOM_REAL) then
      dtc_kv = eta_max_glob * (sqrt(1.0_CUSTOM_REAL + dtc_fault**2/eta_max_glob**2) - 1.0_CUSTOM_REAL)
    else
      dtc_kv = dtc_fault
    endif

    ! user output
    write(IMAIN,*) "  Fault resolution info:"
    write(IMAIN,*) "    Element size   min/max                      = ",elemsize_min_glob,"/",elemsize_max_glob
    write(IMAIN,*) "    P velocity     min/max                      = ",vpmin_glob,"/",vpmax_glob
    write(IMAIN,*) "    estimated maximum suggested time step       = ",dt_suggested
    write(IMAIN,*) "    estimated critical time step                = ",dtc_fault
    write(IMAIN,*) "    suggested eta value range min/max           = ",eta_suggested_min,"/",eta_suggested_max
    if (USE_KELVIN_VOIGT_DAMPING_all) then
      write(IMAIN,*) "    Kelvin Voigt eta(damping):"
      write(IMAIN,*) "      maximum eta value                         = ",eta_max_glob
      write(IMAIN,*) "      estimated critical time step (w/ damping) = ",dtc_kv
    endif
    write(IMAIN,*) "    current time step size                      = ",sngl(DT)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine fault_check_mesh_resolution


end module fault_solver_common
