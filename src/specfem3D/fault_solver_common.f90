! Base module for kinematic and dynamic fault solvers
!
! Authors:
! Percy Galvez, Surendra Somala, Jean-Paul Ampuero

module fault_solver_common

  use constants

  implicit none

!!!!! DK DK  private

  type fault_type
    integer :: nspec=0, nglob=0
    real(kind=CUSTOM_REAL), dimension(:,:),   pointer :: T=>null(),V=>null(),D=>null(),coord=>null()
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: R=>null()
    real(kind=CUSTOM_REAL), dimension(:),     pointer :: B=>null(),invM1=>null(),invM2=>null(),Z=>null()
    real(kind=CUSTOM_REAL) :: dt
    integer, dimension(:), pointer :: ibulk1=>null(), ibulk2=>null()
  end type fault_type

  ! outputs(dyn) /inputs (kind) at selected times for all fault nodes:
  ! strength, state, slip, slip velocity, fault stresses, rupture time, process zone time
  ! rupture time = first time when slip velocity = threshold V_RUPT (defined below)
  ! process zone time = first time when slip = Dc
  type dataXZ_type
    real(kind=CUSTOM_REAL), dimension(:), pointer :: stg=>null(), sta=>null(), d1=>null(), d2=>null(), v1=>null(), v2=>null(), &
                                                     t1=>null(), t2=>null(), t3=>null(), tRUP=>null(), tPZ=>null()
    real(kind=CUSTOM_REAL), dimension(:), pointer :: xcoord=>null(), ycoord=>null(), zcoord=>null()
    integer                                       :: npoin=0
  end type dataXZ_type

  type swf_type
!! DK DK    private
    integer :: kind
    logical :: healing = .false.
    real(kind=CUSTOM_REAL), dimension(:), pointer :: Dc=>null(), mus=>null(), mud=>null(), &
                                                     theta=>null(), T=>null(), C=>null()
  end type swf_type

  type rsf_type
!! DK DK    private
    integer :: StateLaw = 1 ! 1=ageing law, 2=slip law
    real(kind=CUSTOM_REAL), dimension(:), pointer :: V0=>null(), f0=>null(), L=>null(), &
                                                     V_init=>null(), &
                                                     a=>null(), b=>null(), theta=>null(), &
                                                     T=>null(), C=>null(), &
                                                     fw=>null(), Vw=>null()
  end type rsf_type

 ! outputs on selected fault nodes at every time step:
  type dataT_type
    integer :: npoin=0, ndat=0, nt=0
    real(kind=CUSTOM_REAL) :: dt
    integer, dimension(:), pointer :: iglob=>null()   ! on-fault global index of output nodes
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: dat=>null()
    character(len=70), dimension(:), pointer :: name=>null(),longFieldNames=>null()
    character(len=100) :: shortFieldNames
  end type dataT_type

  type, extends (fault_type) :: bc_dynandkinflt_type
!!!!!!!! DK DK      private
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: T0=>null()
    real(kind=CUSTOM_REAL), dimension(:),   pointer :: MU=>null(), Fload=>null()
    integer, dimension(:),   pointer :: npoin_perproc=>null(), poin_offset=>null()
    type(dataT_type)        :: dataT
    type(dataXZ_type)       :: dataXZ,dataXZ_all
    type(swf_type), pointer :: swf => null()
    type(rsf_type), pointer :: rsf => null()
    logical                 :: allow_opening = .false. ! default : do not allow opening

!! DK DK added this in order to be able to use the type for both dynamic and kinematic faults
    real(kind=CUSTOM_REAL) :: kin_dt
    integer :: kin_it
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: v_kin_t1,v_kin_t2

  end type bc_dynandkinflt_type

  logical, parameter :: PARALLEL_FAULT = .true.
 ! NOTE: PARALLEL_FAULT has to be the same
 !       in fault_solver_common.f90, fault_generate_databases.f90 and fault_scotch.f90

  public :: fault_type, PARALLEL_FAULT, &
            initialize_fault, get_jump, get_weighted_jump, rotate, add_BT, &
            dataT_type, init_dataT, store_dataT, SCEC_write_dataT

contains

!---------------------------------------------------------------------

subroutine initialize_fault (bc,IIN_BIN)

  use specfem_par
  use specfem_par_elastic, only : rmassx !,rmassy,rmassz

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(inout) :: bc
  integer, intent(in)                 :: IIN_BIN

  real(kind=CUSTOM_REAL) :: tmp_vec(3,NGLOB_AB)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable   :: jacobian2Dw
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nxyz
  integer, dimension(:,:), allocatable :: ibool1
  integer :: ij,k,e

  read(IIN_BIN) bc%nspec,bc%nglob
  if (.NOT.PARALLEL_FAULT .and. bc%nspec==0) return
  if (bc%nspec>0) then

    allocate(bc%ibulk1(bc%nglob))
    allocate(bc%ibulk2(bc%nglob))
    allocate(bc%R(3,3,bc%nglob))
    allocate(bc%coord(3,(bc%nglob)))
    allocate(bc%invM1(bc%nglob))
    allocate(bc%invM2(bc%nglob))
    allocate(bc%B(bc%nglob))
    allocate(bc%Z(bc%nglob))

    allocate(ibool1(NGLLSQUARE,bc%nspec))
    allocate(normal(NDIM,NGLLSQUARE,bc%nspec))
    allocate(jacobian2Dw(NGLLSQUARE,bc%nspec))

    read(IIN_BIN) ibool1
    read(IIN_BIN) jacobian2Dw
    read(IIN_BIN) normal
    read(IIN_BIN) bc%ibulk1
    read(IIN_BIN) bc%ibulk2
    read(IIN_BIN) bc%coord(1,:)
    read(IIN_BIN) bc%coord(2,:)
    read(IIN_BIN) bc%coord(3,:)

    bc%dt = dt

    bc%B = 0e0_CUSTOM_REAL
    allocate(nxyz(3,bc%nglob))
    nxyz = 0e0_CUSTOM_REAL
    do e=1,bc%nspec
      do ij = 1,NGLLSQUARE
        k = ibool1(ij,e)
        nxyz(:,k) = nxyz(:,k) + normal(:,ij,e)
        bc%B(k) = bc%B(k) + jacobian2Dw(ij,e)
      enddo
    enddo
  endif

  if (PARALLEL_FAULT) then

    tmp_vec = 0._CUSTOM_REAL
    if (bc%nspec>0) tmp_vec(1,bc%ibulk1) = bc%B
    ! assembles with other MPI processes
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,tmp_vec, &
                                     num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     my_neighbours_ext_mesh)
    if (bc%nspec>0) bc%B = tmp_vec(1,bc%ibulk1)

    tmp_vec = 0._CUSTOM_REAL
    if (bc%nspec>0) tmp_vec(:,bc%ibulk1) = nxyz
    ! assembles with other MPI processes
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,tmp_vec, &
                                     num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     my_neighbours_ext_mesh)
    if (bc%nspec>0) nxyz = tmp_vec(:,bc%ibulk1)

  endif

  if (bc%nspec>0) then
    call normalize_3d_vector(nxyz)
    call compute_R(bc%R,bc%nglob,nxyz)

    !SURENDRA : WARNING! Assuming rmassx=rmassy=rmassz
    ! Needed in dA_Free = -K2*d2/M2 + K1*d1/M1
    bc%invM1 = rmassx(bc%ibulk1)
    bc%invM2 = rmassx(bc%ibulk2)

    ! Fault impedance, Z in :  Trac=T_Stick-Z*dV
    !   Z = 1/( B1/M1 + B2/M2 ) / (0.5*dt)
    ! T_stick = Z*Vfree traction as if the fault was stuck (no displ discontinuity)
    ! NOTE: same Bi on both sides, see note above
    bc%Z = 1.e0_CUSTOM_REAL/(0.5e0_CUSTOM_REAL*bc%dt * bc%B *( bc%invM1 + bc%invM2 ))
    ! WARNING: In non-split nodes at fault edges M is assembled across the fault.
    ! hence invM1+invM2=2/(M1+M2) instead of 1/M1+1/M2
    ! In a symmetric mesh (M1=M2) Z will be twice its intended value

  endif

end subroutine initialize_fault

!---------------------------------------------------------------------
subroutine normalize_3d_vector(v)

  real(kind=CUSTOM_REAL), intent(inout) :: v(:,:)

  real(kind=CUSTOM_REAL) :: norm
  integer :: k

 ! assume size(v) = [3,N]
  do k=1,size(v,2)
    norm = sqrt( v(1,k)*v(1,k) + v(2,k)*v(2,k) + v(3,k)*v(3,k) )
    v(:,k) = v(:,k) / norm
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

  integer :: nglob
  real(kind=CUSTOM_REAL), intent(out) :: R(3,3,nglob)
  real(kind=CUSTOM_REAL), intent(in) :: n(3,nglob)

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: s,d

  s(1,:) =  n(2,:)   ! sx = ny
  s(2,:) = -n(1,:)   ! sy =-nx
  s(3,:) = 0.e0_CUSTOM_REAL
  call normalize_3d_vector(s)

  d(1,:) = -s(2,:)*n(3,:) ! dx = -sy*nz
  d(2,:) =  s(1,:)*n(3,:) ! dy = sx*nz
  d(3,:) =  s(2,:)*n(1,:) - s(1,:)*n(2,:) ! dz = sy*nx-ny*sx
  call normalize_3d_vector(d)
  ! dz is always dipwards (negative), because
  ! (nx*sy-ny*sx) = -(nx^2+ny^2)/sqrt(nx^2+ny^2)
  !               = -sqrt(nx^2+ny^2) < 0

  R(1,:,:) = s
  R(2,:,:) = d
  R(3,:,:) = n

end subroutine compute_R


!===============================================================
function get_jump (bc,v) result(dv)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: v(:,:)
  real(kind=CUSTOM_REAL) :: dv(3,bc%nglob)

  ! difference between side 2 and side 1 of fault nodes. dv
  dv(1,:) = v(1,bc%ibulk2)-v(1,bc%ibulk1)
  dv(2,:) = v(2,bc%ibulk2)-v(2,bc%ibulk1)
  dv(3,:) = v(3,bc%ibulk2)-v(3,bc%ibulk1)

end function get_jump

!---------------------------------------------------------------------
function get_weighted_jump (bc,f) result(da)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: f(:,:)

  real(kind=CUSTOM_REAL) :: da(3,bc%nglob)

  ! difference between side 2 and side 1 of fault nodes. M-1 * F
  da(1,:) = bc%invM2*f(1,bc%ibulk2)-bc%invM1*f(1,bc%ibulk1)
  da(2,:) = bc%invM2*f(2,bc%ibulk2)-bc%invM1*f(2,bc%ibulk1)
  da(3,:) = bc%invM2*f(3,bc%ibulk2)-bc%invM1*f(3,bc%ibulk1)

  ! NOTE: In non-split nodes at fault edges M and f are assembled across the fault.
  ! Hence, f1=f2, invM1=invM2=1/(M1+M2) instead of invMi=1/Mi, and da=0.

end function get_weighted_jump

!----------------------------------------------------------------------
function rotate(bc,v,fb) result(vr)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: v(3,bc%nglob)
  integer, intent(in) :: fb
  real(kind=CUSTOM_REAL) :: vr(3,bc%nglob)

  ! Percy, tangential direction Vt, equation 7 of Pablo's notes in agreement with SPECFEM3D

  ! forward rotation
  if (fb==1) then
    vr(1,:) = v(1,:)*bc%R(1,1,:)+v(2,:)*bc%R(1,2,:)+v(3,:)*bc%R(1,3,:) ! vs
    vr(2,:) = v(1,:)*bc%R(2,1,:)+v(2,:)*bc%R(2,2,:)+v(3,:)*bc%R(2,3,:) ! vd
    vr(3,:) = v(1,:)*bc%R(3,1,:)+v(2,:)*bc%R(3,2,:)+v(3,:)*bc%R(3,3,:) ! vn

    !  backward rotation
  else
    vr(1,:) = v(1,:)*bc%R(1,1,:)+v(2,:)*bc%R(2,1,:)+v(3,:)*bc%R(3,1,:)  !vx
    vr(2,:) = v(1,:)*bc%R(1,2,:)+v(2,:)*bc%R(2,2,:)+v(3,:)*bc%R(3,2,:)  !vy
    vr(3,:) = v(1,:)*bc%R(1,3,:)+v(2,:)*bc%R(2,3,:)+v(3,:)*bc%R(3,3,:)  !vz

  endif

end function rotate

!----------------------------------------------------------------------

subroutine add_BT(bc,MxA,T)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  type(bc_dynandkinflt_type), intent(in) :: bc
  real(kind=CUSTOM_REAL), intent(inout) :: MxA(:,:)
  real(kind=CUSTOM_REAL), dimension(3,bc%nglob), intent(in) :: T

  MxA(1,bc%ibulk1) = MxA(1,bc%ibulk1) + bc%B*T(1,:)
  MxA(2,bc%ibulk1) = MxA(2,bc%ibulk1) + bc%B*T(2,:)
  MxA(3,bc%ibulk1) = MxA(3,bc%ibulk1) + bc%B*T(3,:)

  MxA(1,bc%ibulk2) = MxA(1,bc%ibulk2) - bc%B*T(1,:)
  MxA(2,bc%ibulk2) = MxA(2,bc%ibulk2) - bc%B*T(2,:)
  MxA(3,bc%ibulk2) = MxA(3,bc%ibulk2) - bc%B*T(3,:)

end subroutine add_BT


!===============================================================
! dataT outputs

subroutine init_dataT(dataT,coord,nglob,NT,DT,ndat,iflt)

  use specfem_par, only : NPROC,myrank

  integer, intent(in) :: nglob,NT,iflt,ndat
  real(kind=CUSTOM_REAL), intent(in) :: coord(3,nglob),DT
!! DK DK use type(dataT_type) instead of class(dataT_type) for compatibility with some current compilers
  type(dataT_type), intent(out) :: dataT

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: dist_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dist_loc
  integer, dimension(:,:), allocatable :: iglob_all
  integer, dimension(:), allocatable :: iproc,iglob_tmp,glob_indx
  real(kind=CUSTOM_REAL) :: xtarget,ytarget,ztarget,dist,distkeep
  integer :: i, iglob , IIN, ier, jflt, np, k
  character(len=70) :: tmpname
  character(len=70), dimension(:), allocatable :: name_tmp
  integer :: ipoin, ipoin_local, npoin_local

  !  1. read fault output coordinates from user file,
  !  2. define iglob: the fault global index of the node nearest to user
  !     requested coordinate

  IIN = 251 ! WARNING: not safe, should check that unit is not aleady opened

 ! count the number of output points on the current fault (#iflt)
  open(IIN,file='../DATA/FAULT_STATIONS',status='old',action='read',iostat=ier)
  if (ier /= 0) then
    if (myrank==0) write(IMAIN,*) 'Fatal error opening FAULT_STATIONS file. Abort.'
    stop
  endif
  read(IIN,*) np
  dataT%npoin =0
  do i=1,np
    read(IIN,*) xtarget,ytarget,ztarget,tmpname,jflt
    if (jflt==iflt) dataT%npoin = dataT%npoin +1
  enddo
  close(IIN)

  if (dataT%npoin == 0) return

  allocate(dataT%iglob(dataT%npoin))
  allocate(dataT%name(dataT%npoin))
  allocate(dist_loc(dataT%npoin)) !Surendra : for parallel fault

  open(IIN,file='../DATA/FAULT_STATIONS',status='old',action='read')
  read(IIN,*) np
  k = 0
  do i=1,np
    read(IIN,*) xtarget,ytarget,ztarget,tmpname,jflt
    if (jflt/=iflt) cycle
    k = k+1
    dataT%name(k) = tmpname

   ! search nearest node
    distkeep = huge(distkeep)
    do iglob=1,nglob
      dist = sqrt( (coord(1,iglob)-xtarget)**2 &
                 + (coord(2,iglob)-ytarget)**2 &
                 + (coord(3,iglob)-ztarget)**2 )
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
    allocate(iproc(dataT%npoin))
    allocate(iglob_all(dataT%npoin,0:NPROC-1))
    allocate(dist_all(dataT%npoin,0:NPROC-1))
    call gather_all_i(dataT%iglob,dataT%npoin,iglob_all,dataT%npoin,NPROC)
    call gather_all_cr(dist_loc,dataT%npoin,dist_all,dataT%npoin,NPROC)
    if (myrank==0) then
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

    if (npoin_local>0) then
     ! Make a list of output points contained in the current processor
      allocate(glob_indx(npoin_local))
      ipoin_local = 0
      do ipoin = 1,dataT%npoin
        if (myrank == iproc(ipoin)) then
          ipoin_local = ipoin_local + 1
          glob_indx(ipoin_local) = ipoin
        endif
      enddo
     ! Consolidate the output information (remove output points outside current proc)
      allocate(iglob_tmp(dataT%npoin))
      allocate(name_tmp(dataT%npoin))
      iglob_tmp = dataT%iglob
      name_tmp = dataT%name
      deallocate(dataT%iglob)
      deallocate(dataT%name)
      dataT%npoin = npoin_local
      allocate(dataT%iglob(dataT%npoin))
      allocate(dataT%name(dataT%npoin))
      dataT%iglob = iglob_tmp(glob_indx)
      dataT%name = name_tmp(glob_indx)
      deallocate(glob_indx,iglob_tmp,name_tmp)

    else
      dataT%npoin = 0
      deallocate(dataT%iglob)
      deallocate(dataT%name)
    endif

    deallocate(iproc,iglob_all,dist_all)
  endif

 !  3. initialize arrays
  if (dataT%npoin>0) then
    dataT%ndat = ndat
    dataT%nt = NT
    dataT%dt = DT
    allocate(dataT%dat(dataT%ndat,dataT%nt,dataT%npoin))
    dataT%dat = 0e0_CUSTOM_REAL
    allocate(dataT%longFieldNames(dataT%ndat))
    dataT%longFieldNames(1) = "horizontal right-lateral slip (m)"
    dataT%longFieldNames(2) = "horizontal right-lateral slip rate (m/s)"
    dataT%longFieldNames(3) = "horizontal right-lateral shear stress (MPa)"
    dataT%longFieldNames(4) = "vertical up-dip slip (m)"
    dataT%longFieldNames(5) = "vertical up-dip slip rate (m/s)"
    dataT%longFieldNames(6) = "vertical up-dip shear stress (MPa)"
    dataT%longFieldNames(7) = "normal stress (MPa)"
    dataT%shortFieldNames = "h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress"
  endif

end subroutine init_dataT

!---------------------------------------------------------------
subroutine store_dataT(dataT,d,v,t,itime)

  !use specfem_par, only : myrank
!! DK DK use type() instead of class() for compatibility with some current compilers
  type(dataT_type), intent(inout) :: dataT
  real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: d,v,t
  integer, intent(in) :: itime

  integer :: i,k

  do i=1,dataT%npoin
    k = dataT%iglob(i)
    dataT%dat(1,itime,i) = d(1,k)
    dataT%dat(2,itime,i) = v(1,k)
    dataT%dat(3,itime,i) = t(1,k)/1.0e6_CUSTOM_REAL
    dataT%dat(4,itime,i) = -d(2,k)
    dataT%dat(5,itime,i) = -v(2,k)
    dataT%dat(6,itime,i) = -t(2,k)/1.0e6_CUSTOM_REAL
    dataT%dat(7,itime,i) = t(3,k)/1.0e6_CUSTOM_REAL
  enddo

end subroutine store_dataT

!------------------------------------------------------------------------
subroutine SCEC_write_dataT(dataT)

!! DK DK use type() instead of class() for compatibility with some current compilers
  type(dataT_type), intent(in) :: dataT

  integer   :: i,k,IOUT
  character(len=10) :: my_fmt

  integer, dimension(8) :: time_values

  call date_and_time(VALUES=time_values)

  IOUT = 121 !WARNING: not very robust. Could instead look for an available ID

  write(my_fmt,'(a,i1,a)') '(',dataT%ndat+1,'(E15.7))'

  do i=1,dataT%npoin
    open(IOUT,file='../OUTPUT_FILES/'//trim(dataT%name(i))//'.dat',status='replace')
    write(IOUT,*) "# problem=TPV104" ! WARNING: this should be a user input
    write(IOUT,*) "# author=Surendra Nadh Somala" ! WARNING: this should be a user input
    write(IOUT,1000) time_values(2), time_values(3), time_values(1), time_values(5), time_values(6), time_values(7)
    write(IOUT,*) "# code=SPECFEM3D_Cartesian (split nodes)"
    write(IOUT,*) "# code_version=1.1"
    write(IOUT,*) "# element_size=100 m  (*5 GLL nodes)" ! WARNING: this should be a user input
    write(IOUT,*) "# time_step=",dataT%dt
    write(IOUT,*) "# location=",trim(dataT%name(i))
    write(IOUT,*) "# Column #1 = Time (s)"
    do k=1,dataT%ndat
      write(IOUT,1100) k+1,dataT%longFieldNames(k)
    enddo
    write(IOUT,*) "#"
    write(IOUT,*) "# The line below lists the names of the data fields:"
    write(IOUT,*) "# t " // trim(dataT%shortFieldNames)
    write(IOUT,*) "#"
    do k=1,dataT%nt
      write(IOUT,my_fmt) k*dataT%dt, dataT%dat(:,k,i)
    enddo
    close(IOUT)
  enddo

1000 format ( ' # Date = ', i2.2, '/', i2.2, '/', i4.4, '; time = ',i2.2, ':', i2.2, ':', i2.2 )
1100 format ( ' # Column #', i1, ' = ',a )

end subroutine SCEC_write_dataT

end module fault_solver_common
