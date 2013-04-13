! Generates database for faults (dynamic or kinematic)
!
! Splitting fault nodes (with opening) is done in CUBIT.
! See sections "Mesh generation with split nodes"
! and "Cubit-python scripts for faults" in the README_SPECFEM3D_FAULT file.
!
! Authors:
! Percy Galvez, Jean-Paul Ampuero and Tarje Nissen-Meyer

module fault_generate_databases

  use create_regions_mesh_ext_par, only: NGLLX,NGLLY,NGLLZ,NGLLSQUARE,NDIM,CUSTOM_REAL,IMAIN
  use generate_databases_par, only : NGNOD2D

  implicit none
  private

  type fault_db_type
    private
    integer :: nspec=0,nglob=0
    real(kind=CUSTOM_REAL) :: eta
    integer, dimension(:), pointer:: ispec1, ispec2, ibulk1, ibulk2, iface1, iface2
    real(kind=CUSTOM_REAL), dimension(:), pointer :: xcoordbulk1,ycoordbulk1,zcoordbulk1,xcoordbulk2,ycoordbulk2,zcoordbulk2
    integer, dimension(:,:), pointer :: ibool1, ibool2, inodes1, inodes2
    integer, dimension(:,:,:), pointer :: ijk1, ijk2
    real(kind=CUSTOM_REAL), dimension(:,:), pointer:: jacobian2Dw
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer:: normal
  end type fault_db_type

  type(fault_db_type), allocatable, save :: fault_db(:)
  ! fault_db(i) is the database of the i-th fault in the mesh
  real(kind=CUSTOM_REAL), allocatable, save :: Kelvin_Voigt_eta(:)
  double precision, allocatable, save :: nodes_coords_open(:,:)
  integer, save :: nnodes_coords_open
  logical, save :: ANY_FAULT_IN_THIS_PROC = .false.
  logical, save :: ANY_FAULT = .false.

  logical, parameter :: PARALLEL_FAULT = .true.
 ! NOTE: PARALLEL_FAULT has to be the same
 !       in fault_solver_common.f90, fault_generate_databases.f90 and fault_scotch.f90

  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/))   ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
             reshape( (/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/))   ! xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
             reshape( (/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/))   ! ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
             reshape( (/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/))   ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/))  ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
             reshape( (/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/))   ! top
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
             reshape( (/ iface1_corner_ijk,iface2_corner_ijk, &
                 iface3_corner_ijk,iface4_corner_ijk, &
                 iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/))   ! all faces

  public :: fault_read_input, fault_setup, fault_save_arrays_test, fault_save_arrays, &
            nodes_coords_open, nnodes_coords_open, ANY_FAULT_IN_THIS_PROC, ANY_FAULT, PARALLEL_FAULT

contains

!=================================================================================================================
subroutine fault_read_input(prname,myrank)

  character(len=256), intent(in) :: prname
  integer, intent(in) :: myrank

  integer :: nb,i,iflt,ier,nspec,dummy_node
  integer, parameter :: IIN_PAR = 100

 ! read fault input file
  nb = 0
  open(unit=IIN_PAR,file='../DATA/Par_file_faults',status='old',action='read',iostat=ier)
  if (ier==0) then
    read(IIN_PAR,*) nb
    if (myrank==0) write(IMAIN,*) '  ... reading ', nb,' faults from file DATA/Par_file_faults'
  else
    if (myrank==0) write(IMAIN,*) 'File DATA/Par_file_faults not found: assume no faults'
    close(IIN_PAR)
  endif

  ANY_FAULT = (nb>0)
  if (.not. ANY_FAULT)  return

  allocate(fault_db(nb))
  do i=1,nb
    read(IIN_PAR,*) fault_db(i)%eta
  enddo
  close(IIN_PAR)

 ! read fault database file
  open(unit=IIN_PAR,file=prname(1:len_trim(prname))//'Database_fault', &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    write(IMAIN,*) 'error opening file: ',prname(1:len_trim(prname))//'Database_fault'
    write(IMAIN,*) 'make sure file exists'
    stop
  endif

  do iflt=1,size(fault_db)

    read(IIN_PAR) nspec
    fault_db(iflt)%nspec  = nspec

    if (nspec == 0) cycle

    ANY_FAULT_IN_THIS_PROC = .true.

    allocate(fault_db(iflt)%ispec1(nspec))
    allocate(fault_db(iflt)%inodes1(4,nspec))
    allocate(fault_db(iflt)%ispec2(nspec))
    allocate(fault_db(iflt)%inodes2(4,nspec))

    do i=1,nspec
      read(IIN_PAR) fault_db(iflt)%ispec1(i), fault_db(iflt)%inodes1(:,i)
    enddo
    do i=1,nspec
      read(IIN_PAR) fault_db(iflt)%ispec2(i), fault_db(iflt)%inodes2(:,i)
    enddo

   ! loading ispec1 ispec2 iface1 iface2 of fault elements.
!    allocate(fault_db(iflt)%iface1(nspec))
!    allocate(fault_db(iflt)%iface2(nspec))
!    do i=1,fault_db(iflt)%nspec
!      read(IIN_PAR,*) fault_db(iflt)%ispec1(i), fault_db(iflt)%ispec2(i), &
!                  fault_db(iflt)%iface1(i), fault_db(iflt)%iface2(i)
!    enddo

  enddo

 ! read nodes coordinates of the original version of the mesh, in which faults are open
  read(IIN_PAR) nnodes_coords_open
  allocate(nodes_coords_open(NDIM,nnodes_coords_open))
  do i = 1, nnodes_coords_open
     read(IIN_PAR) dummy_node, nodes_coords_open(:,i)
  enddo

  close(IIN_PAR)

end subroutine fault_read_input

!==================================================================================================================
subroutine fault_setup(ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                       xstore,ystore,zstore,nspec,nglob,myrank)

  integer, intent(in) :: nspec
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(inout) :: xstore,ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  integer, intent(in) :: myrank
  integer, intent(in) :: nglob
  integer, intent(in) :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh),intent(in) :: nodes_coords_ext_mesh

  integer :: iflt

  if (.not. ANY_FAULT_IN_THIS_PROC ) return

  do iflt=1,size(fault_db)

    if (fault_db(iflt)%nspec == 0) cycle
    !NOTE: the small fault opening in *_dummy does not affect this subroutine (see get_element_face_id)
    call setup_iface(fault_db(iflt),nnodes_ext_mesh,nodes_coords_ext_mesh,nspec,nglob,ibool)

    ! saving gll indices for each fault face, needed for ibulks
    call setup_ijk(fault_db(iflt))

    ! ibools = mapping from local indices on the fault (GLL index, element index)
    !          to global indices on the fault
    call setup_ibools(fault_db(iflt),xstore,ystore,zstore,nspec,fault_db(iflt)%nspec*NGLLSQUARE)

    ! ibulks = mapping global indices of fault nodes
    !          from global indices on the fault to global indices on the bulk
    call setup_ibulks(fault_db(iflt),ibool,nspec)

    ! close the fault in (xyz)store_dummy
    call close_fault(fault_db(iflt))

    call setup_Kelvin_Voigt_eta(fault_db(iflt),nspec)

    call save_fault_xyzcoord_ibulk(fault_db(iflt))

    call setup_normal_jacobian(fault_db(iflt),ibool,nspec,nglob,myrank)

  enddo

end subroutine fault_setup


!=============================================================================================================
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    ! sets face id of reference element associated with this face
subroutine setup_iface(fdb,nnodes_ext_mesh,nodes_coords_ext_mesh,nspec,nglob,ibool)

  use create_regions_mesh_ext_par, only: xstore_dummy,ystore_dummy,zstore_dummy

  type(fault_db_type), intent(inout) :: fdb
  integer, intent(in) :: nnodes_ext_mesh,nspec,nglob
  double precision, dimension(NDIM,nnodes_ext_mesh), intent(in) :: nodes_coords_ext_mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGNOD2D) :: xcoord,ycoord,zcoord
  integer :: icorner,e

  allocate(fdb%iface1(fdb%nspec))
  allocate(fdb%iface2(fdb%nspec))
  do e=1,fdb%nspec
   ! side 1
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,fdb%inodes1(icorner,e))
      ycoord(icorner) = nodes_coords_ext_mesh(2,fdb%inodes1(icorner,e))
      zcoord(icorner) = nodes_coords_ext_mesh(3,fdb%inodes1(icorner,e))
    enddo
    call get_element_face_id(fdb%ispec1(e),xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            fdb%iface1(e))
   ! side 2
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,fdb%inodes2(icorner,e))
      ycoord(icorner) = nodes_coords_ext_mesh(2,fdb%inodes2(icorner,e))
      zcoord(icorner) = nodes_coords_ext_mesh(3,fdb%inodes2(icorner,e))
    enddo
    call get_element_face_id(fdb%ispec2(e),xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            fdb%iface2(e))
  enddo

end subroutine setup_iface

!=============================================================================================================
subroutine setup_ijk(fdb)

  type(fault_db_type), intent(inout) :: fdb

  integer :: e,i,j,igll
  integer :: ijk_face1(3,NGLLX,NGLLY), ijk_face2(3,NGLLX,NGLLY)

  allocate(fdb%ijk1(3,NGLLX*NGLLY,fdb%nspec))
  allocate(fdb%ijk2(3,NGLLX*NGLLY,fdb%nspec))

  do e=1,fdb%nspec
    call get_element_face_gll_indices(fdb%iface1(e),ijk_face1,NGLLX,NGLLY)
    call get_element_face_gll_indices(fdb%iface2(e),ijk_face2,NGLLX,NGLLY)
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll + 1
        fdb%ijk1(:,igll,e)=ijk_face1(:,i,j)
        fdb%ijk2(:,igll,e)=ijk_face2(:,i,j)
      enddo
    enddo
  enddo

end subroutine setup_ijk

!=============================================================================================================
 subroutine setup_Kelvin_Voigt_eta(fdb,nspec)

  type(fault_db_type), intent(in) :: fdb
  integer, intent(in) :: nspec ! number of spectral elements in each block

  if (fdb%eta > 0.0_CUSTOM_REAL) then
    if (.not.allocated(Kelvin_Voigt_eta)) then
      allocate(Kelvin_Voigt_eta(nspec))
      Kelvin_Voigt_eta(:) = 0.0_CUSTOM_REAL
    endif
    Kelvin_Voigt_eta(fdb%ispec1) = fdb%eta
    Kelvin_Voigt_eta(fdb%ispec2) = fdb%eta
  endif

 end subroutine

!===============================================================================================================
! The lexicographic ordering of node coordinates
! guarantees that the fault nodes are
! consistently ordered on both sides of the fault,
! such that the K-th node of side 1 is facing the K-th node of side 2

subroutine setup_ibools(fdb,xstore,ystore,zstore,nspec,npointot)

  use generate_databases_par, only: nodes_coords_ext_mesh

  type(fault_db_type), intent(inout) :: fdb
  integer, intent(in) :: nspec,npointot
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: xstore,ystore,zstore

  double precision :: xp(npointot),yp(npointot),zp(npointot),xmin,xmax
  integer :: loc(npointot)
  logical :: ifseg(npointot)
  integer :: ispec,k,igll,ie,je,ke,e

  xmin = minval(nodes_coords_ext_mesh(1,:))
  xmax = maxval(nodes_coords_ext_mesh(1,:))

  k = 0
  do e = 1,fdb%nspec
    ispec = fdb%ispec1(e)
    do igll=1,NGLLSQUARE
      ie=fdb%ijk1(1,igll,e)
      je=fdb%ijk1(2,igll,e)
      ke=fdb%ijk1(3,igll,e)
      k = k+1
      xp(k) = xstore(ie,je,ke,ispec)
      yp(k) = ystore(ie,je,ke,ispec)
      zp(k) = zstore(ie,je,ke,ispec)
    enddo
  enddo
  allocate( fdb%ibool1(NGLLSQUARE,fdb%nspec) )
  call get_global(fdb%nspec,xp,yp,zp,fdb%ibool1,loc,ifseg,fdb%nglob,npointot,xmin,xmax)

! xp,yp,zp need to be recomputed on side 2
! because they are generally not in the same order as on side 1,
! because ispec1(e) is not necessarily facing ispec2(e).

  k = 0
  do e = 1,fdb%nspec
    ispec = fdb%ispec2(e)
    do igll=1,NGLLSQUARE
      ie=fdb%ijk2(1,igll,e)
      je=fdb%ijk2(2,igll,e)
      ke=fdb%ijk2(3,igll,e)
      k = k+1
      xp(k) = xstore(ie,je,ke,ispec)
      yp(k) = ystore(ie,je,ke,ispec)
      zp(k) = zstore(ie,je,ke,ispec)
    enddo
  enddo
  allocate( fdb%ibool2(NGLLSQUARE,fdb%nspec) )
  call get_global(fdb%nspec,xp,yp,zp,fdb%ibool2,loc,ifseg,fdb%nglob,npointot,xmin,xmax)

end subroutine setup_ibools


!=================================================================================

subroutine setup_ibulks(fdb,ibool,nspec)

  type(fault_db_type), intent(inout) :: fdb
  integer, intent(in) :: nspec, ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer :: e,k, K1, K2, ie,je,ke

  allocate( fdb%ibulk1(fdb%nglob) )
  allocate( fdb%ibulk2(fdb%nglob) )

  do e=1, fdb%nspec
    do k=1, NGLLSQUARE

      ie=fdb%ijk1(1,k,e)
      je=fdb%ijk1(2,k,e)
      ke=fdb%ijk1(3,k,e)
      K1= fdb%ibool1(k,e)
      fdb%ibulk1(K1)=ibool(ie,je,ke,fdb%ispec1(e))

      ie=fdb%ijk2(1,k,e)
      je=fdb%ijk2(2,k,e)
      ke=fdb%ijk2(3,k,e)
      K2= fdb%ibool2(k,e)
      fdb%ibulk2(K2)=ibool(ie,je,ke,fdb%ispec2(e))

    enddo
  enddo

end subroutine setup_ibulks

!=============================================================================================================
! We only close *store_dummy.  *store is close already in create_regions_mesh.f90
! Fortunately only *store_dummy is needed to compute jacobians and normals

subroutine close_fault(fdb)

  use create_regions_mesh_ext_par, only: xstore_dummy,ystore_dummy,zstore_dummy

  type(fault_db_type), intent(inout) :: fdb

  integer :: i,K1,K2

  do i=1,fdb%nglob
    K1 = fdb%ibulk1(i)
    K2 = fdb%ibulk2(i)
    xstore_dummy(K1) = 0.5d0*( xstore_dummy(K1) + xstore_dummy(K2) )
    xstore_dummy(K2) = xstore_dummy(K1)
    ystore_dummy(K1) = 0.5d0*( ystore_dummy(K1) + ystore_dummy(K2) )
    ystore_dummy(K2) = ystore_dummy(K1)
    zstore_dummy(K1) = 0.5d0*( zstore_dummy(K1) + zstore_dummy(K2) )
    zstore_dummy(K2) = zstore_dummy(K1)
  enddo

end subroutine close_fault

!=================================================================================

 subroutine save_fault_xyzcoord_ibulk(fdb)

  use create_regions_mesh_ext_par, only: xstore_dummy,ystore_dummy,zstore_dummy

  type(fault_db_type), intent(inout) :: fdb

  integer :: K1, K2, i

  allocate( fdb%xcoordbulk1(fdb%nglob) )
  allocate( fdb%ycoordbulk1(fdb%nglob) )
  allocate( fdb%zcoordbulk1(fdb%nglob) )
  allocate( fdb%xcoordbulk2(fdb%nglob) )
  allocate( fdb%ycoordbulk2(fdb%nglob) )
  allocate( fdb%zcoordbulk2(fdb%nglob) )

  do i=1, fdb%nglob
      K1 =fdb%ibulk1(i)
      K2 =fdb%ibulk2(i)
      fdb%xcoordbulk1(i) = xstore_dummy(K1)
      fdb%ycoordbulk1(i) = ystore_dummy(K1)
      fdb%zcoordbulk1(i) = zstore_dummy(K1)

      fdb%xcoordbulk2(i) = xstore_dummy(K2)
      fdb%ycoordbulk2(i) = ystore_dummy(K2)
      fdb%zcoordbulk2(i) = zstore_dummy(K2)
  enddo

 end subroutine save_fault_xyzcoord_ibulk


!=================================================================================

 subroutine setup_normal_jacobian(fdb,ibool,nspec,nglob,myrank)

  use create_regions_mesh_ext_par, only: xstore_dummy,ystore_dummy,zstore_dummy, &
                                         dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz
  use generate_databases_par, only : NGNOD2D


  type(fault_db_type), intent(inout) :: fdb
  integer, intent(in) :: nspec,nglob, ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer, intent(in) :: myrank

  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  integer,dimension(NGNOD2D) :: iglob_corners_ref
  integer :: ispec_flt,ispec,i,j,k,igll
  integer :: iface_ref,icorner

  allocate(fdb%normal(NDIM,NGLLSQUARE,fdb%nspec))
  allocate(fdb%jacobian2Dw(NGLLSQUARE,fdb%nspec))

  do ispec_flt=1,fdb%nspec

    iface_ref= fdb%iface1(ispec_flt)
    ispec = fdb%ispec1(ispec_flt)

    ! takes indices of corners of reference face
    do icorner = 1,NGNOD2D
      i = iface_all_corner_ijk(1,icorner,iface_ref)
      j = iface_all_corner_ijk(2,icorner,iface_ref)
      k = iface_all_corner_ijk(3,icorner,iface_ref)

      ! global reference indices
      iglob_corners_ref(icorner) = ibool(i,j,k,ispec)

      ! reference corner coordinates
      xcoord(icorner) = xstore_dummy(iglob_corners_ref(icorner))
      ycoord(icorner) = ystore_dummy(iglob_corners_ref(icorner))
      zcoord(icorner) = zstore_dummy(iglob_corners_ref(icorner))
    enddo

    ! gets face GLL 2Djacobian, weighted from element face
    call get_jacobian_boundary_face(myrank,nspec, &
           xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
           dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
           wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
           ispec,iface_ref,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

    ! normal convention: points away from domain1, reference element.
    do j=1,NGLLY
      do i=1,NGLLX
        ! directs normals such that they point outwards of element
        call get_element_face_normal(ispec,iface_ref,xcoord,ycoord,zcoord, &
          ibool,nspec,nglob, &
          xstore_dummy,ystore_dummy,zstore_dummy, &
          normal_face(:,i,j) )
      enddo
    enddo

    ! stores informations about this face
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        ! adds all gll points on that face
        igll = igll + 1
        ! stores weighted jacobian and normals
        fdb%jacobian2Dw(igll,ispec_flt) = jacobian2Dw_face(i,j)
        fdb%normal(:,igll,ispec_flt) = normal_face(:,i,j)
      enddo
    enddo

  enddo ! ispec_flt

end subroutine setup_normal_jacobian

!====================================================================================
! saves all fault data in ASCII files for verification
subroutine fault_save_arrays_test(prname)

  character(len=256), intent(in) :: prname ! 'proc***'

  integer, parameter :: IOUT = 121 !WARNING: not very robust. Could instead look for an available ID
  integer :: nbfaults,iflt,ier
  character(len=256) :: filename

  if (.not.ANY_FAULT) return

! saves mesh file proc***_fault_db.txt
  filename = prname(1:len_trim(prname))//'fault_db.txt'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',iostat=ier)
  if( ier /= 0 ) stop 'error opening database proc######_external_mesh.bin'

  nbfaults = size(fault_db)
  write(IOUT,*) 'NBFAULTS = ',nbfaults
  do iflt=1,nbfaults
    write(IOUT,*) 'BEGIN FAULT # ',iflt
    call save_one_fault_test(fault_db(iflt),IOUT)
    write(IOUT,*) 'END FAULT # ',iflt
  enddo
  close(IOUT)

end subroutine fault_save_arrays_test

!-------------------------------------------------------------------------------------

subroutine save_one_fault_test(f,IOUT)

  type(fault_db_type), intent(in) :: f
  integer, intent(in) :: IOUT

  integer :: e,k
  character(15) :: fmt1,fmt2

  write(fmt1,'("(a,",I0,"(x,I7))")') NGLLSQUARE+1   ! fmt = (a,(NGLL^2+1)(x,I7))
  write(fmt2,'("(a,",I0,"(x,F0.4))")') NGLLSQUARE+1   ! fmt = (a,(NGLL^2+1)(x,F0.16))

  write(IOUT,*) 'NSPEC NGLOB NGLL = ',f%nspec,f%nglob,NGLLX
  if (f%nspec==0) return
  do e=1,f%nspec
    write(IOUT,*) 'FLT_ELEM = ',e
    write(IOUT,*) 'ISPEC1 ISPEC2 = ',f%ispec1(e),f%ispec2(e)
    write(IOUT,fmt1) 'IBOOL1 = ',f%ibool1(:,e)
    write(IOUT,fmt1) 'IBOOL2 = ',f%ibool2(:,e)
    write(IOUT,fmt1) 'I1 = ',f%ijk1(1,:,e)
    write(IOUT,fmt1) 'J1 = ',f%ijk1(2,:,e)
    write(IOUT,fmt1) 'K1 = ',f%ijk1(3,:,e)
    write(IOUT,fmt1) 'I2 = ',f%ijk2(1,:,e)
    write(IOUT,fmt1) 'J2 = ',f%ijk2(2,:,e)
    write(IOUT,fmt1) 'K2 = ',f%ijk2(3,:,e)
    write(IOUT,fmt2) 'JAC2DW = ',f%jacobian2Dw(:,e)
    write(IOUT,fmt2) 'N1 = ',f%normal(1,:,e)
    write(IOUT,fmt2) 'N2 = ',f%normal(2,:,e)
    write(IOUT,fmt2) 'N3 = ',f%normal(3,:,e)
  enddo

  write(IOUT,*) 'FLT_NODE IBULK1 IBULK2'
  do k=1,f%nglob
    write(IOUT,*) k,f%ibulk1(k),f%ibulk2(k)
  enddo

  write(IOUT,*) 'FLT_NODE xcoordbulk ycoordbulk zcoordbulk'
  do k=1,f%nglob
    write(IOUT,*) f%ibulk1(k),f%xcoordbulk1(k),f%ycoordbulk1(k),f%zcoordbulk1(k)
    write(IOUT,*) f%ibulk2(k),f%xcoordbulk2(k),f%ycoordbulk2(k),f%zcoordbulk2(k)
  enddo

end subroutine save_one_fault_test

!=================================================================================
! saves fault data needed by the solver in binary files
subroutine fault_save_arrays(prname)

  character(len=256), intent(in) :: prname ! 'proc***'

  integer, parameter :: IOUT = 121 !WARNING: not very robust. Could instead look for an available ID
  integer :: nbfaults,iflt,ier
  character(len=256) :: filename
  integer :: size_Kelvin_Voigt

  if (.not.ANY_FAULT) return

! opening Kelvin_voig_eta.bin for each processor
! if number of fault elements = 0 then the file is empty
  filename = prname(1:len_trim(prname))//'Kelvin_voigt_eta.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    write(IMAIN,*) 'error opening file ',trim(filename)
    stop
  endif

! saves mesh file proc***_Kelvin_voigt_eta.bin
  if (allocated(Kelvin_Voigt_eta)) then
    size_Kelvin_Voigt = size(Kelvin_Voigt_eta)
  else
    size_Kelvin_Voigt = 0
  endif
  write(IOUT) size_Kelvin_Voigt
  if (size_Kelvin_Voigt /= 0) Write(IOUT) Kelvin_Voigt_eta
  close(IOUT)

! saves mesh file proc***_fault_db.bin
  filename = prname(1:len_trim(prname))//'fault_db.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    write(IMAIN,*) 'error opening file ',trim(filename)
    stop
  endif

  nbfaults = size(fault_db)
  write(IOUT) nbfaults
  do iflt=1,nbfaults
    call save_one_fault_bin(fault_db(iflt),IOUT)
  enddo
  close(IOUT)


end subroutine fault_save_arrays

!----------------------------------------------

subroutine save_one_fault_bin(f,IOUT)

  type(fault_db_type), intent(in) :: f
  integer, intent(in) :: IOUT

  write(IOUT) f%nspec,f%nglob
  if (f%nspec==0) return
  write(IOUT) f%ibool1
  write(IOUT) f%jacobian2Dw
  write(IOUT) f%normal
  write(IOUT) f%ibulk1
  write(IOUT) f%ibulk2
  write(IOUT) f%xcoordbulk1
  write(IOUT) f%ycoordbulk1
  write(IOUT) f%zcoordbulk1

end subroutine save_one_fault_bin

!------------------------------------------------


end module fault_generate_databases
