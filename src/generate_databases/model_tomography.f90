!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

!--------------------------------------------------------------------------------------------------
!
! generic tomography file
!
! note: the idea is to use external tomography velocity models
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------

  module tomography_par

  use constants,only: CUSTOM_REAL,IMODEL_TOMO

  implicit none

  ! for external tomography:
  ! (regular spaced, xyz-block file in ascii)

  ! number of external tomographic models
  integer :: NFILES_TOMO

  ! models dimensions
  double precision  :: END_X,END_Y,END_Z

  double precision, dimension(:), allocatable :: ORIG_X,ORIG_Y,ORIG_Z
  double precision, dimension(:), allocatable :: SPACING_X,SPACING_Y,SPACING_Z

  ! models parameter records
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vp_tomography,vs_tomography,rho_tomography,z_tomography

  ! models entries
  integer, dimension(:), allocatable :: NX,NY,NZ
  integer, dimension(:), allocatable :: nrecord

  ! min/max statistics
  double precision, dimension(:), allocatable :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX

  ! process rank
  integer :: myrank_tomo

  end module tomography_par

!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_tomography_broadcast(myrank)

  use tomography_par, only: myrank_tomo

  implicit none

  integer :: myrank

  ! stores rank
  myrank_tomo = myrank

  ! all processes read in same file
  ! note: for a high number of processes this might lead to a bottleneck

  ! determines the number of external tomographic models and allocates tomographic arrays
  call init_tomography_files()

  ! makes sure all processes are initialized
  call synchronize_all()

  ! reads Vp Vs and rho from extracted text file
  call read_model_tomography()

  ! otherwise:

  ! only master reads in model file
  !integer, dimension(1) :: nrecord
  !if(myrank == 0) call read_model_tomography()
  ! broadcast the information read on the master to the nodes, e.g.
  !call bcast_all_i(nrecord,1)

  !if( myrank /= 0 ) then
  ! allocate( vp_tomography(1:nrecord) ,stat=ier)
  ! if( ier /= 0 ) stop 'error allocating array vp_tomography'
  !endif
  !call bcast_all_cr(vp_tomography,size(vp_tomography))

  ! synchronizes processes
  call synchronize_all()

  end subroutine model_tomography_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine init_tomography_files()

! determines the number of external tomographic models and sets a total maximum number of element records

  use constants, only: MAX_STRING_LEN,IIN,IMAIN

  use generate_databases_par, only: TOMOGRAPHY_PATH,undef_mat_prop,nundefMat_ext_mesh,IMODEL,COUPLE_WITH_EXTERNAL_CODE

  use tomography_par

  implicit none

  ! local parameters
  double precision :: dummy,temp_x,temp_y,temp_z
  integer :: ier,iundef,nrecord_max,ifiles_tomo,nrec,nlines
  character(len=MAX_STRING_LEN*2) :: tomo_filename
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: string_read
  integer :: nmaterials

  ! sets number of materials to loop over
  nmaterials = nundefMat_ext_mesh

  NFILES_TOMO = 0
  nrecord_max = 0
  ifiles_tomo = 0

  ! checks if we over-impose a tomography model by Par_file setting: MODEL = tomo
  if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
    nmaterials = 1
  endif

  ! loops over number of undefined materials
  do iundef = 1, nmaterials

    ! sets filename
    if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
      ! note: since we have no undefined materials, we cannot access undef_mat_prop(:,:) to read in values
      ! uses default name
      filename = 'tomography_model.xyz'
    else
      ! checks if associated material is a tomography model
      if(trim(undef_mat_prop(2,iundef)) /= 'tomography') cycle

      ! gets filename from material property
      read(undef_mat_prop(4,iundef),*) filename
    endif

    ! counter
    ifiles_tomo = ifiles_tomo + 1

    ! sets filename with path (e.g. "DATA/tomo_files/" + "tomo.xyz")
    ! corrects the path and filename of tomography model
    if ( TOMOGRAPHY_PATH(len_trim(TOMOGRAPHY_PATH):len_trim(TOMOGRAPHY_PATH)) == "/"  )then
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // trim(filename)
    else
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // '/' // trim(filename)
    endif

    !! WANGYI test for the benchmark of hybrid DSM-SPECFEM3D coupling
    if (COUPLE_WITH_EXTERNAL_CODE) then
       write(*,*) 'iundef',iundef  ! add by WANGYI for test
       write(*,*) 'tomo_filename',tomo_filename ! add by WANGYI for test
    endif

    ! opens file for reading
    open(unit=IIN,file=trim(tomo_filename),status='old',action='read',iostat=ier)
    if(ier /= 0) then
      print*,'Error: could not open tomography file: ',trim(tomo_filename)
      print*,'Please check your settings in Par_file ...'
      call exit_MPI(myrank_tomo,'error reading tomography file')
    endif

    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) dummy,dummy,dummy,dummy,dummy,dummy

    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) dummy,dummy,dummy

    ! reads in models entries
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) temp_x,temp_y,temp_z

    ! determines total maximum number of element records
    nrec = int(temp_x*temp_y*temp_z)
    nrecord_max = max(nrecord_max,nrec)

    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) dummy,dummy,dummy,dummy,dummy,dummy

    ! checks the number of records for points definition
    nlines = 0
    do while(ier == 0)
      read(IIN,*,iostat=ier)
      if (ier == 0) nlines = nlines + 1
    enddo

    if( nlines /= nrec .and. myrank_tomo == 0 ) then
       print*, 'Error: ',trim(tomo_filename),' has invalid number of records'
       print*, '     number of grid points specified (= NX*NY*NZ):',nrec
       print*, '     number of file lines for grid points        :',nlines
       stop 'Error in tomography data file for the grid points definition'
    endif

    ! closes file
    close(IIN)

  enddo

  ! number of external tomographic models
  NFILES_TOMO = ifiles_tomo

  ! user output
  if( myrank_tomo == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) '     number of tomographic models       = ',NFILES_TOMO
    write(IMAIN,*) '     maximum number of data records     = ',nrecord_max
    write(IMAIN,*) '     size of required tomography arrays = ', &
      sngl( 4.d0 * NFILES_TOMO * nrecord_max * CUSTOM_REAL / 1024.d0 /1024.d0),'MB per process'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks if we found a tomography model
  if (NFILES_TOMO == 0 ) call exit_MPI(myrank_tomo,'Error no tomography model was read in')

  ! allocates models dimensions
  allocate(ORIG_X(NFILES_TOMO),ORIG_Y(NFILES_TOMO),ORIG_Z(NFILES_TOMO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')
  allocate(SPACING_X(NFILES_TOMO),SPACING_Y(NFILES_TOMO),SPACING_Z(NFILES_TOMO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! allocate models parameter records
  allocate(vp_tomography(NFILES_TOMO,nrecord_max), &
           vs_tomography(NFILES_TOMO,nrecord_max), &
           rho_tomography(NFILES_TOMO,nrecord_max), &
           z_tomography(NFILES_TOMO,nrecord_max),stat=ier)
  if(ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! allocate models entries
  allocate(NX(NFILES_TOMO),NY(NFILES_TOMO),NZ(NFILES_TOMO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')
  allocate(nrecord(NFILES_TOMO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! allocate models min/max statistics
  allocate(VP_MIN(NFILES_TOMO),VS_MIN(NFILES_TOMO),RHO_MIN(NFILES_TOMO), &
           VP_MAX(NFILES_TOMO),VS_MAX(NFILES_TOMO),RHO_MAX(NFILES_TOMO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

end subroutine init_tomography_files

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_tomography()

! read Vp Vs and rho from extracted text file

  use constants, only: MAX_STRING_LEN,IIN,IMAIN

  use generate_databases_par, only: TOMOGRAPHY_PATH,undef_mat_prop,nundefMat_ext_mesh,IMODEL

  use tomography_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo
  integer :: irecord,ier,iundef,imat
  character(len=MAX_STRING_LEN*2) :: tomo_filename
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: string_read
  integer :: nmaterials

  ! sets number of materials to loop over
  nmaterials = nundefMat_ext_mesh

  ! checks if we over-impose a tomography model by Par_file setting: MODEL = tomo
  if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
    nmaterials = 1
  endif

  imat = 0
  do iundef = 1, nmaterials

    ! sets filename
    if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
      ! note: since we have no undefined materials, we cannot access undef_mat_prop(:,:) to read in values
      ! uses default name
      filename = 'tomography_model.xyz'
    else
      ! checks if associated material is a tomography model
      if(trim(undef_mat_prop(2,iundef)) /= 'tomography') cycle

      ! gets filename from material property
      read(undef_mat_prop(4,iundef),*) filename
    endif

    ! counter
    imat = imat + 1

    ! sets filename with path (e.g. "DATA/tomo_files/" + "tomo.xyz")
    ! corrects the path and filename of tomography model
    if ( TOMOGRAPHY_PATH(len_trim(TOMOGRAPHY_PATH):len_trim(TOMOGRAPHY_PATH)) == "/"  )then
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // trim(filename)
    else
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // '/' // trim(filename)
    endif

    ! user output
    if( myrank_tomo == 0 ) then
       write(IMAIN,*) '     reading: ',trim(tomo_filename)
       call flush_IMAIN()
    endif

    ! opens file for reading
    open(unit=IIN,file=trim(tomo_filename),status='old',action='read',iostat=ier)
    if(ier /= 0) call exit_MPI(myrank_tomo,'error reading tomography file')

    rewind(unit=IIN,iostat=ier)
    if(ier /= 0) call exit_MPI(myrank_tomo,'error rewinding tomography file')

    !--------------------------------------
    ! header infos
    !--------------------------------------
    ! reads in model dimensions
    ! format: #origin_x #origin_y #origin_z #end_x #end_y #end_z
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) ORIG_X(imat), ORIG_Y(imat), ORIG_Z(imat), END_X, END_Y, END_Z

    ! model increments
    ! format: #dx #dy #dz
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) SPACING_X(imat), SPACING_Y(imat), SPACING_Z(imat)

    ! reads in models entries
    ! format: #nx #ny #nz
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) NX(imat), NY(imat), NZ(imat)

    ! reads in models min/max statistics
    ! format: #vp_min #vp_max #vs_min #vs_max #density_min #density_max
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) VP_MIN(imat), VP_MAX(imat), &
                        VS_MIN(imat), VS_MAX(imat), &
                        RHO_MIN(imat), RHO_MAX(imat)

    ! total number of element records
    nrecord(imat) = NX(imat)*NY(imat)*NZ(imat)

    !--------------------------------------
    ! data records
    !--------------------------------------
    ! first record
    ! format: #x #y #z #vp #vs #density
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo

    ! stores record values
    vp_tomography(imat,1) = vp_tomo
    vs_tomography(imat,1) = vs_tomo
    rho_tomography(imat,1) = rho_tomo
    z_tomography(imat,1) = z_tomo

    ! reads in record sections
    do irecord = 2,nrecord(imat)
       read(IIN,*) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo

       ! stores record values
       vp_tomography(imat,irecord) = vp_tomo
       vs_tomography(imat,irecord) = vs_tomo
       rho_tomography(imat,irecord) = rho_tomo
       z_tomography(imat,irecord) = z_tomo
    enddo

    close(IIN)

    ! user output
    if( myrank_tomo == 0 ) then
      write(IMAIN,*) '     number of grid points = NX*NY*NZ:',nrecord(imat)
      write(IMAIN,*)
    endif
  enddo

  end subroutine read_model_tomography

!
!-------------------------------------------------------------------------------------------------
!

  subroutine tomo_read_next_line(unit_in,string_read)

  use constants, only: MAX_STRING_LEN
  implicit none

  integer :: unit_in
  character(len=MAX_STRING_LEN) :: string_read

  integer :: ios

  do
     read(unit=unit_in,fmt="(a)",iostat=ios) string_read
     if(ios /= 0) stop 'error while reading tomography file'

     ! suppress leading white spaces, if any
     string_read = adjustl(string_read)

     ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
     if(index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

     ! exit loop when we find the first line that is not a comment or a white line
     if(len_trim(string_read) == 0) cycle
     if(string_read(1:1) /= '#') exit
  enddo

  ! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

  ! suppress trailing comments, if any
  if(index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

  ! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine tomo_read_next_line

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_tomography(xmesh,ymesh,zmesh,rho_model,vp_model,vs_model,qkappa_atten,qmu_atten,imaterial_id, &
                              has_tomo_value)

  use generate_databases_par, only: undef_mat_prop,nundefMat_ext_mesh,IMODEL,ATTENUATION_COMP_MAXIMUM

  use tomography_par

  implicit none

  double precision, intent(in) :: xmesh,ymesh,zmesh

  real(kind=CUSTOM_REAL), intent(out) :: qkappa_atten,qmu_atten

  real(kind=CUSTOM_REAL), intent(out) :: vp_model,vs_model,rho_model

  integer, intent(in) :: imaterial_id
  logical,intent(out) :: has_tomo_value

  ! local parameters
  integer :: ix,iy,iz,imat,iundef
  integer :: p0,p1,p2,p3,p4,p5,p6,p7

  double precision :: spac_x,spac_y,spac_z
  double precision :: gamma_interp_x,gamma_interp_y
  double precision :: gamma_interp_z1,gamma_interp_z2,gamma_interp_z3, &
    gamma_interp_z4,gamma_interp_z5,gamma_interp_z6,gamma_interp_z7,gamma_interp_z8

  real(kind=CUSTOM_REAL) :: vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8, &
    vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8

  real(kind=CUSTOM_REAL), dimension(NFILES_TOMO) :: vp_final,vs_final,rho_final
  integer :: nmaterials,imaterial

  ! initializes flag
  has_tomo_value = .false.

  ! sets number of materials to loop over
  nmaterials = nundefMat_ext_mesh

  ! sets material number
  imaterial = abs(imaterial_id)

  ! checks if we over-impose a tomography model by Par_file setting: MODEL = tomo
  if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
    nmaterials = 1
    imaterial = 1
  else
    ! checks if material is a tomographic material (negative id)
    if( imaterial_id >= 0 ) return

    ! checks if associated type is a tomography model
    if( trim(undef_mat_prop(2,abs(imaterial_id))) /= 'tomography' ) return
  endif

  ! loops over all undefined materials
  imat = 0
  do iundef = 1, nmaterials

    ! counter
    imat = imat + 1

    if( imat /= imaterial ) cycle

    ! determine spacing and cell for linear interpolation
    spac_x = (xmesh - ORIG_X(imat)) / SPACING_X(imat)
    spac_y = (ymesh - ORIG_Y(imat)) / SPACING_Y(imat)
    spac_z = (zmesh - ORIG_Z(imat)) / SPACING_Z(imat)

    ix = int(spac_x)
    iy = int(spac_y)
    iz = int(spac_z)

    gamma_interp_x = spac_x - dble(ix)
    gamma_interp_y = spac_y - dble(iy)

    ! suppress edge effects for points outside of the model SPOSTARE DOPO
    if(ix < 0) then
       ix = 0
       gamma_interp_x = 0.d0
    endif
    if(ix > NX(imat)-2) then
       ix = NX(imat)-2
       gamma_interp_x = 1.d0
    endif

    if(iy < 0) then
       iy = 0
       gamma_interp_y = 0.d0
    endif
    if(iy > NY(imat)-2) then
       iy = NY(imat)-2
       gamma_interp_y = 1.d0
    endif

    if(iz < 0) then
       iz = 0
       !   gamma_interp_z = 0.d0
    endif
    if(iz > NZ(imat)-2) then
       iz = NZ(imat)-2
       !  gamma_interp_z = 1.d0
    endif


    ! define 8 corners of interpolation element
    p0 = ix+iy*NX(imat)+iz*(NX(imat)*NY(imat))
    p1 = (ix+1)+iy*NX(imat)+iz*(NX(imat)*NY(imat))
    p2 = (ix+1)+(iy+1)*NX(imat)+iz*(NX(imat)*NY(imat))
    p3 = ix+(iy+1)*NX(imat)+iz*(NX(imat)*NY(imat))
    p4 = ix+iy*NX(imat)+(iz+1)*(NX(imat)*NY(imat))
    p5 = (ix+1)+iy*NX(imat)+(iz+1)*(NX(imat)*NY(imat))
    p6 = (ix+1)+(iy+1)*NX(imat)+(iz+1)*(NX(imat)*NY(imat))
    p7 = ix+(iy+1)*NX(imat)+(iz+1)*(NX(imat)*NY(imat))

    if(p0 < 0 .or. p1 < 0 .or. p2 < 0 .or. p3 < 0 .or. p4 < 0 .or. p5 < 0 .or. p6 < 0 .or. p7 < 0) then
       print*,'model: ',imat
       print*,'error rank: ',myrank_tomo
       print*,'corner index: ',p0,p1,p2,p3,p4,p5,p6,p7
       print*,'location: ',sngl(xmesh),sngl(ymesh),sngl(zmesh)
       print*,'origin: ',sngl(ORIG_X(imat)),sngl(ORIG_Y(imat)),sngl(ORIG_Z(imat))
       call exit_MPI(myrank_tomo,'error corner index in tomography routine')
    endif

    if(z_tomography(imat,p4+1) == z_tomography(imat,p0+1)) then
       gamma_interp_z1 = 1.d0
    else
       gamma_interp_z1 = (zmesh-z_tomography(imat,p0+1))/(z_tomography(imat,p4+1)-z_tomography(imat,p0+1))
    endif
    if(gamma_interp_z1 > 1.d0) then
       gamma_interp_z1 = 1.d0
    endif
    if(gamma_interp_z1 < 0.d0) then
       gamma_interp_z1 = 0.d0
    endif


    if(z_tomography(imat,p5+1) == z_tomography(imat,p1+1)) then
       gamma_interp_z2 = 1.d0
    else
       gamma_interp_z2 = (zmesh-z_tomography(imat,p1+1))/(z_tomography(imat,p5+1)-z_tomography(imat,p1+1))
    endif
    if(gamma_interp_z2 > 1.d0) then
       gamma_interp_z2 = 1.d0
    endif
    if(gamma_interp_z2 < 0.d0) then
       gamma_interp_z2 = 0.d0
    endif


    if(z_tomography(imat,p6+1) == z_tomography(imat,p2+1)) then
       gamma_interp_z3 = 1.d0
    else
       gamma_interp_z3 = (zmesh-z_tomography(imat,p2+1))/(z_tomography(imat,p6+1)-z_tomography(imat,p2+1))
    endif
    if(gamma_interp_z3 > 1.d0) then
       gamma_interp_z3 = 1.d0
    endif
    if(gamma_interp_z3 < 0.d0) then
       gamma_interp_z3 = 0.d0
    endif


    if(z_tomography(imat,p7+1) == z_tomography(imat,p3+1)) then
       gamma_interp_z4 = 1.d0
    else
       gamma_interp_z4 = (zmesh-z_tomography(imat,p3+1))/(z_tomography(imat,p7+1)-z_tomography(imat,p3+1))
    endif
    if(gamma_interp_z4 > 1.d0) then
       gamma_interp_z4 = 1.d0
    endif
    if(gamma_interp_z4 < 0.d0) then
       gamma_interp_z4 = 0.d0
    endif

    gamma_interp_z5 = 1.d0 - gamma_interp_z1
    gamma_interp_z6 = 1.d0 - gamma_interp_z2
    gamma_interp_z7 = 1.d0 - gamma_interp_z3
    gamma_interp_z8 = 1.d0 - gamma_interp_z4

    vp1 = vp_tomography(imat,p0+1)
    vp2 = vp_tomography(imat,p1+1)
    vp3 = vp_tomography(imat,p2+1)
    vp4 = vp_tomography(imat,p3+1)
    vp5 = vp_tomography(imat,p4+1)
    vp6 = vp_tomography(imat,p5+1)
    vp7 = vp_tomography(imat,p6+1)
    vp8 = vp_tomography(imat,p7+1)

    vs1 = vs_tomography(imat,p0+1)
    vs2 = vs_tomography(imat,p1+1)
    vs3 = vs_tomography(imat,p2+1)
    vs4 = vs_tomography(imat,p3+1)
    vs5 = vs_tomography(imat,p4+1)
    vs6 = vs_tomography(imat,p5+1)
    vs7 = vs_tomography(imat,p6+1)
    vs8 = vs_tomography(imat,p7+1)

    rho1 = rho_tomography(imat,p0+1)
    rho2 = rho_tomography(imat,p1+1)
    rho3 = rho_tomography(imat,p2+1)
    rho4 = rho_tomography(imat,p3+1)
    rho5 = rho_tomography(imat,p4+1)
    rho6 = rho_tomography(imat,p5+1)
    rho7 = rho_tomography(imat,p6+1)
    rho8 = rho_tomography(imat,p7+1)

    ! use trilinear interpolation in cell to define Vp Vs and rho
    vp_final(imat) = &
         vp1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
         vp2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
         vp3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
         vp4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
         vp5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
         vp6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
         vp7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
         vp8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4

    vs_final(imat) = &
         vs1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
         vs2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
         vs3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
         vs4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
         vs5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
         vs6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
         vs7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
         vs8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4

    rho_final(imat) = &
         rho1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
         rho2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
         rho3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
         rho4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
         rho5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
         rho6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
         rho7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
         rho8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4

    ! impose minimum and maximum velocity and density if needed
    if(vp_final(imat) < VP_MIN(imat)) vp_final(imat) = VP_MIN(imat)
    if(vp_final(imat) > VP_MAX(imat)) vp_final(imat) = VP_MAX(imat)

    if(vs_final(imat) < VS_MIN(imat)) vs_final(imat) = VS_MIN(imat)
    if(vs_final(imat) > VS_MAX(imat)) vs_final(imat) = VS_MAX(imat)

    if(rho_final(imat) > RHO_MAX(imat)) rho_final(imat) = RHO_MAX(imat)
    if(rho_final(imat) < RHO_MIN(imat)) rho_final(imat) = RHO_MIN(imat)

  enddo

  ! checks that material was found
  if (imat /= imaterial) stop 'Error specified tomographic material id was not found'

  ! model parameters for the associated negative imaterial_id index in materials file
  rho_model = rho_final(imaterial)
  vp_model = vp_final(imaterial)
  vs_model = vs_final(imaterial)

  ! attenuation: arbitrary value, see maximum in constants.h
  qmu_atten = ATTENUATION_COMP_MAXIMUM

!! DK DK Q_kappa is not implemented in this model_tomography routine yet, thus set it to dummy value
  qkappa_atten = 9999.

  ! value found
  has_tomo_value = .true.

  end subroutine model_tomography

!
!-------------------------------------------------------------------------------------------------
!

  subroutine deallocate_tomography_files()

    use tomography_par

    implicit none

    ! deallocates models dimensions
    deallocate(ORIG_X,ORIG_Y,ORIG_Z)
    deallocate(SPACING_X,SPACING_Y,SPACING_Z)

    ! deallocates models parameter records
    deallocate(vp_tomography)
    deallocate(vs_tomography)
    deallocate(rho_tomography)
    deallocate(z_tomography)

    ! deallocates models entries
    deallocate(NX,NY,NZ)
    deallocate(nrecord)

    ! deallocates models min/max statistics
    deallocate(VP_MIN,VS_MIN,RHO_MIN)
    deallocate(VP_MAX,VS_MAX,RHO_MAX)

  end subroutine deallocate_tomography_files
