!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

!--------------------------------------------------------------------------------------------------
!
! generic tomography file
!
! note: the idea is to use external tomography velocity models
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------

  module model_tomography_par

  use constants, only: CUSTOM_REAL,IMODEL_TOMO
  use generate_databases_par, only: IMODEL,ANISOTROPY

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
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: qp_tomography,qs_tomography
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: c_tomography

  ! models entries
  integer, dimension(:), allocatable :: NX,NY,NZ
  integer, dimension(:), allocatable :: nrecord

  ! min/max statistics
  double precision, dimension(:), allocatable :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX

  ! process rank
  integer :: myrank_tomo

  ! data format
  logical,dimension(:),allocatable :: tomo_has_q_values

  end module model_tomography_par

!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_tomography_broadcast(myrank)

  use model_tomography_par, only: myrank_tomo

  implicit none

  integer,intent(in) :: myrank

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

  ! only main reads in model file
  !integer, dimension(1) :: nrecord
  !if (myrank == 0) call read_model_tomography()
  ! broadcast the information read on the main to the nodes, e.g.
  !call bcast_all_i(nrecord,1)

  !if (myrank /= 0) then
  ! allocate( vp_tomography(1:nrecord) ,stat=ier)
  ! if (ier /= 0) stop 'error allocating array vp_tomography'
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

  use generate_databases_par, only: TOMOGRAPHY_PATH,undef_mat_prop,nundefMat_ext_mesh

  use model_tomography_par

  implicit none

  ! local parameters
  double precision :: dummy,temp_x,temp_y,temp_z
  integer :: ier,iundef,nrecord_max,ifiles_tomo,nrec,nlines
  character(len=MAX_STRING_LEN*2) :: tomo_filename
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: string_read
  integer :: nmaterials
  ! data format
  logical :: has_q_values
  integer :: ntokens
  logical,dimension(:),allocatable :: materials_with_q

  ! sets number of materials to loop over
  nmaterials = nundefMat_ext_mesh

  NFILES_TOMO = 0
  nrecord_max = 0
  ifiles_tomo = 0

  ! checks if we over-impose a tomography model by Par_file setting: MODEL = tomo
  if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
    nmaterials = 1
  endif

  ! data format flag
  allocate(materials_with_q(nmaterials),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 853')
  if (ier /= 0) stop 'Error allocating array materials_with_q'
  materials_with_q(:) = .false.

  ! loops over number of undefined materials
  do iundef = 1, nmaterials

    ! sets filename
    if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
      ! note: since we have no undefined materials, we cannot access undef_mat_prop(:,:) to read in values
      ! uses default name
      filename = 'tomography_model.xyz'
    else
      ! checks if associated material is a tomography model
      if (trim(undef_mat_prop(2,iundef)) /= 'tomography') cycle

      ! gets filename from material property
      read(undef_mat_prop(4,iundef),*) filename
    endif

    ! counter
    ifiles_tomo = ifiles_tomo + 1

    ! sets filename with path (e.g. "DATA/tomo_files/" + "tomo.xyz")
    ! corrects the path and filename of tomography model
    if (TOMOGRAPHY_PATH(len_trim(TOMOGRAPHY_PATH):len_trim(TOMOGRAPHY_PATH)) == "/") then
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // trim(filename)
    else
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // '/' // trim(filename)
    endif

    ! opens file for reading
    open(unit=IIN,file=trim(tomo_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open tomography file: ',trim(tomo_filename)
      print *,'Please check your settings in Par_file ...'
      call exit_MPI(myrank_tomo,'error reading tomography file')
    endif

    ! header infos
    ! format: #origin_x #origin_y #origin_z #end_x #end_y #end_z
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) dummy,dummy,dummy,dummy,dummy,dummy

    ! format: #dx #dy #dz
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) dummy,dummy,dummy

    ! reads in models entries
    ! format: #nx #ny #nz
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) temp_x,temp_y,temp_z

    ! determines total maximum number of element records
    nrec = int(temp_x*temp_y*temp_z)
    nrecord_max = max(nrecord_max,nrec)

    ! format: #vp_min #vp_max #vs_min #vs_max #density_min #density_max
    call tomo_read_next_line(IIN,string_read)
    read(string_read,*) dummy,dummy,dummy,dummy,dummy,dummy

    ! data records
    ! checks the number of records for points definition
    ! first record
    nlines = 1
    call tomo_read_next_line(IIN,string_read) ! allows to skip comment lines before data section

    ! checks number of entries of first data line
    call tomo_get_number_of_tokens(string_read,ntokens)
    !print *,'tomography file: number of tokens on first data line: ',ntokens,' line: ',trim(string_read)
    if (ntokens /= 6 .and. ntokens /= 8 .and. ntokens /= 25 .and. ntokens /= 27) then
      print *,'Error reading tomography file, data line has wrong number of entries: ',trim(string_read)
      stop 'Error reading tomography file'
    endif

    if (ANISOTROPY .and. IMODEL == IMODEL_TOMO) then
      if (ntokens == 6 .or. ntokens == 8) then
        print *,'Error reading tomography file, data line has wrong number of entries: ',trim(string_read)
        stop 'Error reading tomography file'
      endif
    endif

    ! determines data format
    if (ntokens == 8 .or. ntokens == 27) then
      has_q_values = .true.
    else
      has_q_values = .false.
    endif
    materials_with_q(ifiles_tomo) = has_q_values

    ! counts remaining records
    do while (ier == 0)
      read(IIN,*,iostat=ier)
      if (ier == 0) nlines = nlines + 1
    enddo

    if (nlines /= nrec .and. myrank_tomo == 0) then
       print *, 'Error: ',trim(tomo_filename),' has invalid number of records'
       print *, '     number of grid points specified (= NX*NY*NZ):',nrec
       print *, '     number of file lines for grid points        :',nlines
       stop 'Error in tomography data file for the grid points definition'
    endif

    ! closes file
    close(IIN)

  enddo

  ! number of external tomographic models
  NFILES_TOMO = ifiles_tomo

  ! user output
  if (myrank_tomo == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '     number of tomographic models       = ',NFILES_TOMO
    write(IMAIN,*) '     maximum number of data records     = ',nrecord_max
    write(IMAIN,*) '     size of required tomography arrays = ', &
      sngl( 4.d0 * NFILES_TOMO * nrecord_max * CUSTOM_REAL / 1024.d0 /1024.d0),'MB per process'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks if we found a tomography model
  if (NFILES_TOMO == 0) call exit_MPI(myrank_tomo,'Error no tomography model was read in')

  ! allocates models dimensions
  allocate(ORIG_X(NFILES_TOMO),ORIG_Y(NFILES_TOMO),ORIG_Z(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 854')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')
  allocate(SPACING_X(NFILES_TOMO),SPACING_Y(NFILES_TOMO),SPACING_Z(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 855')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! allocate models parameter records
  ! only allocate anisotropy arrays if needed
  if (ANISOTROPY .and. IMODEL == IMODEL_TOMO) then
    allocate(c_tomography(NFILES_TOMO,nrecord_max,21),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 904X')
    if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo anisotropy arrays')
  else
    allocate(vp_tomography(NFILES_TOMO,nrecord_max),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 856')
    allocate(vs_tomography(NFILES_TOMO,nrecord_max),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 857')
  endif

  allocate(rho_tomography(NFILES_TOMO,nrecord_max),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 858')
  allocate(z_tomography(NFILES_TOMO,nrecord_max),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 859')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! allocate models entries
  allocate(NX(NFILES_TOMO),NY(NFILES_TOMO),NZ(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 860')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')
  allocate(nrecord(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 861')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! allocate models min/max statistics
  allocate(VP_MIN(NFILES_TOMO),VS_MIN(NFILES_TOMO),RHO_MIN(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 862')
  allocate(VP_MAX(NFILES_TOMO),VS_MAX(NFILES_TOMO),RHO_MAX(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 863')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo arrays')

  ! q values
  allocate(tomo_has_q_values(NFILES_TOMO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 864')
  if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo q-flag array')
  tomo_has_q_values(:) = .false.
  ! stores data format flag
  do ifiles_tomo = 1,NFILES_TOMO
    tomo_has_q_values(ifiles_tomo) = materials_with_q(ifiles_tomo)
  enddo
  deallocate(materials_with_q)

  ! only allocate q arrays if needed
  if (any(tomo_has_q_values)) then
    allocate(qp_tomography(NFILES_TOMO,nrecord_max),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 865')
    allocate(qs_tomography(NFILES_TOMO,nrecord_max),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 866')
    if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate tomo q-value arrays')
  endif

end subroutine init_tomography_files

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_tomography()

! read Vp Vs and rho from extracted text file
! also read Qp Qs if needed
! also read c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66 if needed

  use constants, only: MAX_STRING_LEN,IIN,IMAIN

  use generate_databases_par, only: TOMOGRAPHY_PATH,undef_mat_prop,nundefMat_ext_mesh

  use model_tomography_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo
  real(kind=CUSTOM_REAL) :: qp_tomo,qs_tomo
  real(kind=CUSTOM_REAL) :: c11_tomo,c12_tomo,c13_tomo,c14_tomo,c15_tomo,c16_tomo,c22_tomo,c23_tomo,c24_tomo,c25_tomo,c26_tomo, &
                            c33_tomo,c34_tomo,c35_tomo,c36_tomo,c44_tomo,c45_tomo,c46_tomo,c55_tomo,c56_tomo,c66_tomo
  integer :: irecord,ier,iundef,imat
  character(len=MAX_STRING_LEN*2) :: tomo_filename
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: string_read
  integer :: nmaterials
  logical :: has_q_values

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
      if (trim(undef_mat_prop(2,iundef)) /= 'tomography') cycle

      ! gets filename from material property
      read(undef_mat_prop(4,iundef),*) filename
    endif

    ! sets filename with path (e.g. "DATA/tomo_files/" + "tomo.xyz")
    ! corrects the path and filename of tomography model
    if (TOMOGRAPHY_PATH(len_trim(TOMOGRAPHY_PATH):len_trim(TOMOGRAPHY_PATH)) == "/") then
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // trim(filename)
    else
      tomo_filename = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH)) // '/' // trim(filename)
    endif

    ! counter
    imat = imat + 1

    ! user output
    if (myrank_tomo == 0) then
       write(IMAIN,*) '     material id: ',-imat
       write(IMAIN,*) '     file       : ',trim(tomo_filename)
       call flush_IMAIN()
    endif

    ! opens file for reading
    open(unit=IIN,file=trim(tomo_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank_tomo,'Error opening tomography file')

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

    ! determines data format
    has_q_values = tomo_has_q_values(imat)

    if (ANISOTROPY .and. IMODEL == IMODEL_TOMO) then
      ! user output
      if (myrank_tomo == 0) then
        if (has_q_values) then
          write(IMAIN,*) '     data format: #x #y #z #c11 #c12 .... #c55 #c56 #c66 #density #Q_p #Q_s'
          ! #c22 #c23 #c24 #c25 #c26 #c33 #c34 #c35 #c36 #c44 #c45 #c46
        else
          write(IMAIN,*) '     data format: #x #y #z #c11 #c12 .... #c55 #c56 #c66 #density'
          ! #c22 #c23 #c24 #c25 #c26 #c33 #c34 #c35 #c36 #c44 #c45 #c46
        endif
        call flush_IMAIN()
      endif

      ! reads in first data values
      if (has_q_values) then
        ! format: #x #y #z #c11 #c12 #c13 #c14 #c15 #c16 #c22 #c23 #c24 #c25 #c26 #c33 #c34 #c35 #c36 #c44 #c45 #c46 #c55 #c56 #c66
        !         #density #Q_p #Q_s
        read(string_read,*) x_tomo,y_tomo,z_tomo,c11_tomo,c12_tomo,c13_tomo,c14_tomo,c15_tomo,c16_tomo,c22_tomo,c23_tomo,c24_tomo, &
                                                 c25_tomo,c26_tomo,c33_tomo,c34_tomo,c35_tomo,c36_tomo,c44_tomo,c45_tomo,c46_tomo, &
                                                 c55_tomo,c56_tomo,c66_tomo,rho_tomo,qp_tomo,qs_tomo
        qp_tomography(imat,1) = qp_tomo
        qs_tomography(imat,1) = qs_tomo
      else
        ! format: #x #y #z #c11 #c12 #c13 #c14 #c15 #c16 #c22 #c23 #c24 #c25 #c26 #c33 #c34 #c35 #c36 #c44 #c45 #c46 #c55 #c56 #c66
        !         #density
        read(string_read,*) x_tomo,y_tomo,z_tomo,c11_tomo,c12_tomo,c13_tomo,c14_tomo,c15_tomo,c16_tomo,c22_tomo,c23_tomo,c24_tomo, &
                                                 c25_tomo,c26_tomo,c33_tomo,c34_tomo,c35_tomo,c36_tomo,c44_tomo,c45_tomo,c46_tomo, &
                                                 c55_tomo,c56_tomo,c66_tomo,rho_tomo
      endif

      ! stores record values
      c_tomography(imat,1,1) = c11_tomo
      c_tomography(imat,1,2) = c12_tomo
      c_tomography(imat,1,3) = c13_tomo
      c_tomography(imat,1,4) = c14_tomo
      c_tomography(imat,1,5) = c15_tomo
      c_tomography(imat,1,6) = c16_tomo
      c_tomography(imat,1,7) = c22_tomo
      c_tomography(imat,1,8) = c23_tomo
      c_tomography(imat,1,9) = c24_tomo
      c_tomography(imat,1,10) = c25_tomo
      c_tomography(imat,1,11) = c26_tomo
      c_tomography(imat,1,12) = c33_tomo
      c_tomography(imat,1,13) = c34_tomo
      c_tomography(imat,1,14) = c35_tomo
      c_tomography(imat,1,15) = c36_tomo
      c_tomography(imat,1,16) = c44_tomo
      c_tomography(imat,1,17) = c45_tomo
      c_tomography(imat,1,18) = c46_tomo
      c_tomography(imat,1,19) = c55_tomo
      c_tomography(imat,1,20) = c56_tomo
      c_tomography(imat,1,21) = c66_tomo

      rho_tomography(imat,1) = rho_tomo
      z_tomography(imat,1) = z_tomo

      ! reads in record sections
      if (has_q_values) then
        do irecord = 2,nrecord(imat)
          !format: #x #y #z #c11 #c12 #c13 #c14 #c15 #c16 #c22 #c23 #c24 #c25 #c26 #c33 #c34 #c35 #c36 #c44 #c45 #c46 #c55 #c56 #c66
          !        #density #Q_p #Q_s
          read(IIN,*,iostat=ier) x_tomo,y_tomo,z_tomo,c11_tomo,c12_tomo,c13_tomo,c14_tomo,c15_tomo,c16_tomo,c22_tomo,c23_tomo, &
                                                      c24_tomo,c25_tomo,c26_tomo,c33_tomo,c34_tomo,c35_tomo,c36_tomo,c44_tomo, &
                                                      c45_tomo,c46_tomo,c55_tomo,c56_tomo,c66_tomo,rho_tomo,qp_tomo,qs_tomo
          if (ier /= 0) stop 'Error reading tomo file line format with q values'

          ! stores record values
          c_tomography(imat,irecord,1) = c11_tomo
          c_tomography(imat,irecord,2) = c12_tomo
          c_tomography(imat,irecord,3) = c13_tomo
          c_tomography(imat,irecord,4) = c14_tomo
          c_tomography(imat,irecord,5) = c15_tomo
          c_tomography(imat,irecord,6) = c16_tomo
          c_tomography(imat,irecord,7) = c22_tomo
          c_tomography(imat,irecord,8) = c23_tomo
          c_tomography(imat,irecord,9) = c24_tomo
          c_tomography(imat,irecord,10) = c25_tomo
          c_tomography(imat,irecord,11) = c26_tomo
          c_tomography(imat,irecord,12) = c33_tomo
          c_tomography(imat,irecord,13) = c34_tomo
          c_tomography(imat,irecord,14) = c35_tomo
          c_tomography(imat,irecord,15) = c36_tomo
          c_tomography(imat,irecord,16) = c44_tomo
          c_tomography(imat,irecord,17) = c45_tomo
          c_tomography(imat,irecord,18) = c46_tomo
          c_tomography(imat,irecord,19) = c55_tomo
          c_tomography(imat,irecord,20) = c56_tomo
          c_tomography(imat,irecord,21) = c66_tomo

          rho_tomography(imat,irecord) = rho_tomo
          z_tomography(imat,irecord) = z_tomo
          qp_tomography(imat,irecord) = qp_tomo
          qs_tomography(imat,irecord) = qs_tomo
        enddo
      else
        do irecord = 2,nrecord(imat)
          !format: #x #y #z #c11 #c12 #c13 #c14 #c15 #c16 #c22 #c23 #c24 #c25 #c26 #c33 #c34 #c35 #c36 #c44 #c45 #c46 #c55 #c56 #c66
          !        #density
          read(IIN,*,iostat=ier) x_tomo,y_tomo,z_tomo,c11_tomo,c12_tomo,c13_tomo,c14_tomo,c15_tomo,c16_tomo,c22_tomo,c23_tomo, &
                                                      c24_tomo,c25_tomo,c26_tomo,c33_tomo,c34_tomo,c35_tomo,c36_tomo,c44_tomo, &
                                                      c45_tomo,c46_tomo,c55_tomo,c56_tomo,c66_tomo,rho_tomo
          if (ier /= 0) stop 'Error reading tomo file line format'

          ! stores record values
          c_tomography(imat,irecord,1) = c11_tomo
          c_tomography(imat,irecord,2) = c12_tomo
          c_tomography(imat,irecord,3) = c13_tomo
          c_tomography(imat,irecord,4) = c14_tomo
          c_tomography(imat,irecord,5) = c15_tomo
          c_tomography(imat,irecord,6) = c16_tomo
          c_tomography(imat,irecord,7) = c22_tomo
          c_tomography(imat,irecord,8) = c23_tomo
          c_tomography(imat,irecord,9) = c24_tomo
          c_tomography(imat,irecord,10) = c25_tomo
          c_tomography(imat,irecord,11) = c26_tomo
          c_tomography(imat,irecord,12) = c33_tomo
          c_tomography(imat,irecord,13) = c34_tomo
          c_tomography(imat,irecord,14) = c35_tomo
          c_tomography(imat,irecord,15) = c36_tomo
          c_tomography(imat,irecord,16) = c44_tomo
          c_tomography(imat,irecord,17) = c45_tomo
          c_tomography(imat,irecord,18) = c46_tomo
          c_tomography(imat,irecord,19) = c55_tomo
          c_tomography(imat,irecord,20) = c56_tomo
          c_tomography(imat,irecord,21) = c66_tomo

          rho_tomography(imat,irecord) = rho_tomo
          z_tomography(imat,irecord) = z_tomo
        enddo
      endif
      close(IIN)

      ! user output
      if (myrank_tomo == 0) then
        write(IMAIN,*) '     number of grid points = NX*NY*NZ:',nrecord(imat)
        write(IMAIN,*)
        call flush_IMAIN()
      endif

    else
      ! user output
      if (myrank_tomo == 0) then
        if (has_q_values) then
          write(IMAIN,*) '     data format: #x #y #z #vp #vs #density #Q_p #Q_s'
        else
          write(IMAIN,*) '     data format: #x #y #z #vp #vs #density'
        endif
        call flush_IMAIN()
      endif

      ! reads in first data values
      if (has_q_values) then
        ! format: #x #y #z #vp #vs #density #Q_p #Q_s
        read(string_read,*) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo,qp_tomo,qs_tomo
        qp_tomography(imat,1) = qp_tomo
        qs_tomography(imat,1) = qs_tomo
      else
        ! format: #x #y #z #vp #vs #density
        read(string_read,*) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo
      endif

      ! stores record values
      vp_tomography(imat,1) = vp_tomo
      vs_tomography(imat,1) = vs_tomo
      rho_tomography(imat,1) = rho_tomo
      z_tomography(imat,1) = z_tomo

      ! reads in record sections
      if (has_q_values) then
        do irecord = 2,nrecord(imat)
          ! format: #x #y #z #vp #vs #density #Q_p #Q_s
          read(IIN,*,iostat=ier) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo,qp_tomo,qs_tomo
          if (ier /= 0) stop 'Error reading tomo file line format with q values'

          ! stores record values
          vp_tomography(imat,irecord) = vp_tomo
          vs_tomography(imat,irecord) = vs_tomo
          rho_tomography(imat,irecord) = rho_tomo
          z_tomography(imat,irecord) = z_tomo
          qp_tomography(imat,irecord) = qp_tomo
          qs_tomography(imat,irecord) = qs_tomo
        enddo
      else
        do irecord = 2,nrecord(imat)
          ! format: #x #y #z #vp #vs #density
          read(IIN,*,iostat=ier) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo
          if (ier /= 0) stop 'Error reading tomo file line format'

          ! stores record values
          vp_tomography(imat,irecord) = vp_tomo
          vs_tomography(imat,irecord) = vs_tomo
          rho_tomography(imat,irecord) = rho_tomo
          z_tomography(imat,irecord) = z_tomo
        enddo
      endif
      close(IIN)

      ! user output
      if (myrank_tomo == 0) then
        write(IMAIN,*) '     number of grid points = NX*NY*NZ:',nrecord(imat)
        write(IMAIN,*)
        call flush_IMAIN()
      endif
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

  integer :: ier

  do
     read(unit=unit_in,fmt="(a)",iostat=ier) string_read
     if (ier /= 0) stop 'error while reading tomography file'

     ! suppress leading white spaces, if any
     string_read = adjustl(string_read)

     ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
     if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

     ! reads next line if empty
     if (len_trim(string_read) == 0) cycle

     ! exit loop when we find the first line that is not a comment or a white line
     if (string_read(1:1) /= '#') exit
  enddo

  ! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

  ! suppress trailing comments, if any
  if (index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

  ! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine tomo_read_next_line

!
!-------------------------------------------------------------------------------------------------
!

  subroutine tomo_get_number_of_tokens(string_read,ntokens)

  use constants, only: MAX_STRING_LEN
  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: string_read
  integer,intent(out) :: ntokens

  ! local parameters
  integer :: i
  logical :: previous_is_delim
  character(len=1), parameter :: delim_space = ' '
  character(len=1), parameter :: delim_tab = achar(9) ! tab delimiter

  ! initializes
  ntokens = 0

  ! checks if anything to do
  if (len_trim(string_read) == 0) return

  ! counts tokens
  ntokens = 1
  previous_is_delim = .true.
  do i = 1, len_trim(string_read)
    ! finds next delimiter (space or tabular)
    if (string_read(i:i) == delim_space .or. string_read(i:i) == delim_tab) then
      if (.not. previous_is_delim) then
        ntokens = ntokens + 1
        previous_is_delim = .true.
      endif
    else
      if (previous_is_delim) previous_is_delim = .false.
    endif
  enddo

  end subroutine tomo_get_number_of_tokens


!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_tomography(xmesh,ymesh,zmesh,rho_model,vp_model,vs_model,qkappa_atten,qmu_atten, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              imaterial_id,has_tomo_value)

  use generate_databases_par, only: undef_mat_prop,nundefMat_ext_mesh,ATTENUATION_COMP_MAXIMUM

  use model_tomography_par

  implicit none

  double precision, intent(in) :: xmesh,ymesh,zmesh

  real(kind=CUSTOM_REAL), intent(out) :: qkappa_atten,qmu_atten

  real(kind=CUSTOM_REAL), intent(out) :: vp_model,vs_model,rho_model

  real(kind=CUSTOM_REAL), intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  integer, intent(in) :: imaterial_id
  logical,intent(out) :: has_tomo_value

  ! local parameters
  integer :: ier,i
  integer :: ix,iy,iz,imat
  integer :: p0,p1,p2,p3,p4,p5,p6,p7

  double precision :: spac_x,spac_y,spac_z
  double precision :: gamma_interp_x,gamma_interp_y
  double precision :: gamma_interp_z1,gamma_interp_z2,gamma_interp_z3,gamma_interp_z4

  real(kind=CUSTOM_REAL) :: vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8
  real(kind=CUSTOM_REAL) :: vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8
  real(kind=CUSTOM_REAL) :: rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8

  real(kind=CUSTOM_REAL) :: vp_final,vs_final,rho_final

  ! attenuation
  real(kind=CUSTOM_REAL) :: qp1,qp2,qp3,qp4,qp5,qp6,qp7,qp8
  real(kind=CUSTOM_REAL) :: qs1,qs2,qs3,qs4,qs5,qs6,qs7,qs8
  real(kind=CUSTOM_REAL) :: qp_final,qs_final
  real(kind=CUSTOM_REAL) :: L_val

  ! anisotropy
  real(kind=CUSTOM_REAL) :: c1,c2,c3,c4,c5,c6,c7,c8
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: c_final

  ! initializes flag
  has_tomo_value = .false.

  ! checks if we over-impose a tomography model by Par_file setting: MODEL = tomo
  if (nundefMat_ext_mesh == 0 .and. IMODEL == IMODEL_TOMO) then
    ! sets material number
    imat = 1
  else
    ! checks if material is a tomographic material (negative id)
    if (imaterial_id >= 0) return

    ! sets material number
    imat = abs(imaterial_id)

    ! checks if associated type is a tomography model
    if (trim(undef_mat_prop(2,imat)) /= 'tomography') return

    ! checks material
    if (imat < 1 .or. imat > nundefMat_ext_mesh) then
      print *,'Error tomography model: unknown material id ',imaterial_id,' for ',nundefMat_ext_mesh,' undefined materials'
      stop 'Error unknown material id in tomography model'
    endif
  endif

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
  if (ix < 0) then
     ix = 0
     gamma_interp_x = 0.d0
  endif
  if (ix > NX(imat)-2) then
     ix = NX(imat)-2
     gamma_interp_x = 1.d0
  endif

  if (iy < 0) then
     iy = 0
     gamma_interp_y = 0.d0
  endif
  if (iy > NY(imat)-2) then
     iy = NY(imat)-2
     gamma_interp_y = 1.d0
  endif

  if (iz < 0) then
     iz = 0
     !   gamma_interp_z = 0.d0
  endif
  if (iz > NZ(imat)-2) then
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

  if (p0 < 0 .or. p1 < 0 .or. p2 < 0 .or. p3 < 0 .or. p4 < 0 .or. p5 < 0 .or. p6 < 0 .or. p7 < 0) then
     print *,'model: ',imat
     print *,'error rank: ',myrank_tomo
     print *,'corner index: ',p0,p1,p2,p3,p4,p5,p6,p7
     print *,'location: ',sngl(xmesh),sngl(ymesh),sngl(zmesh)
     print *,'origin: ',sngl(ORIG_X(imat)),sngl(ORIG_Y(imat)),sngl(ORIG_Z(imat))
     call exit_MPI(myrank_tomo,'error corner index in tomography routine')
  endif

  ! interpolation gamma factors
  if (z_tomography(imat,p4+1) == z_tomography(imat,p0+1)) then
     gamma_interp_z1 = 1.d0
  else
     gamma_interp_z1 = (zmesh-z_tomography(imat,p0+1))/(z_tomography(imat,p4+1)-z_tomography(imat,p0+1))
  endif
  if (gamma_interp_z1 > 1.d0) then
     gamma_interp_z1 = 1.d0
  endif
  if (gamma_interp_z1 < 0.d0) then
     gamma_interp_z1 = 0.d0
  endif

  if (z_tomography(imat,p5+1) == z_tomography(imat,p1+1)) then
     gamma_interp_z2 = 1.d0
  else
     gamma_interp_z2 = (zmesh-z_tomography(imat,p1+1))/(z_tomography(imat,p5+1)-z_tomography(imat,p1+1))
  endif
  if (gamma_interp_z2 > 1.d0) then
     gamma_interp_z2 = 1.d0
  endif
  if (gamma_interp_z2 < 0.d0) then
     gamma_interp_z2 = 0.d0
  endif

  if (z_tomography(imat,p6+1) == z_tomography(imat,p2+1)) then
     gamma_interp_z3 = 1.d0
  else
     gamma_interp_z3 = (zmesh-z_tomography(imat,p2+1))/(z_tomography(imat,p6+1)-z_tomography(imat,p2+1))
  endif
  if (gamma_interp_z3 > 1.d0) then
     gamma_interp_z3 = 1.d0
  endif
  if (gamma_interp_z3 < 0.d0) then
     gamma_interp_z3 = 0.d0
  endif

  if (z_tomography(imat,p7+1) == z_tomography(imat,p3+1)) then
     gamma_interp_z4 = 1.d0
  else
     gamma_interp_z4 = (zmesh-z_tomography(imat,p3+1))/(z_tomography(imat,p7+1)-z_tomography(imat,p3+1))
  endif
  if (gamma_interp_z4 > 1.d0) then
     gamma_interp_z4 = 1.d0
  endif
  if (gamma_interp_z4 < 0.d0) then
     gamma_interp_z4 = 0.d0
  endif

  ! density
  rho1 = rho_tomography(imat,p0+1)
  rho2 = rho_tomography(imat,p1+1)
  rho3 = rho_tomography(imat,p2+1)
  rho4 = rho_tomography(imat,p3+1)
  rho5 = rho_tomography(imat,p4+1)
  rho6 = rho_tomography(imat,p5+1)
  rho7 = rho_tomography(imat,p6+1)
  rho8 = rho_tomography(imat,p7+1)
  ! use trilinear interpolation in cell to define rho
  rho_final = interpolate_trilinear(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8)

  ! impose minimum and maximum density if needed
  if (rho_final > RHO_MAX(imat)) rho_final = RHO_MAX(imat)
  if (rho_final < RHO_MIN(imat)) rho_final = RHO_MIN(imat)

  ! model parameters for the associated negative imaterial_id index in materials file
  rho_model = rho_final

  if (ANISOTROPY .and. IMODEL == IMODEL_TOMO) then

    ! anisotropy
    allocate(c_final(21),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 905X')
    if (ier /= 0) call exit_MPI(myrank_tomo,'not enough memory to allocate interpolated anisotropy parameters array')

    do i = 1,21
      c1 = c_tomography(imat,p0+1,i)
      c2 = c_tomography(imat,p1+1,i)
      c3 = c_tomography(imat,p2+1,i)
      c4 = c_tomography(imat,p3+1,i)
      c5 = c_tomography(imat,p4+1,i)
      c6 = c_tomography(imat,p5+1,i)
      c7 = c_tomography(imat,p6+1,i)
      c8 = c_tomography(imat,p7+1,i)
      ! use trilinear interpolation in cell to define Vp
      c_final(i) = interpolate_trilinear(c1,c2,c3,c4,c5,c6,c7,c8)
    enddo

    c11 = c_final(1)
    c12 = c_final(2)
    c13 = c_final(3)
    c14 = c_final(4)
    c15 = c_final(5)
    c16 = c_final(6)
    c22 = c_final(7)
    c23 = c_final(8)
    c24 = c_final(9)
    c25 = c_final(10)
    c26 = c_final(11)
    c33 = c_final(12)
    c34 = c_final(13)
    c35 = c_final(14)
    c36 = c_final(15)
    c44 = c_final(16)
    c45 = c_final(17)
    c46 = c_final(18)
    c55 = c_final(19)
    c56 = c_final(20)
    c66 = c_final(21)

    vp_model = sqrt(c11)/sqrt(rho_model) ! a better estimate of equivalent vp is needed for anisotropic models
    vs_model = sqrt(c66)/sqrt(rho_model) ! a better estimate of equivalent vs is needed for anisotropic models

    deallocate(c_final)

  else

    ! Vp
    vp1 = vp_tomography(imat,p0+1)
    vp2 = vp_tomography(imat,p1+1)
    vp3 = vp_tomography(imat,p2+1)
    vp4 = vp_tomography(imat,p3+1)
    vp5 = vp_tomography(imat,p4+1)
    vp6 = vp_tomography(imat,p5+1)
    vp7 = vp_tomography(imat,p6+1)
    vp8 = vp_tomography(imat,p7+1)
    ! use trilinear interpolation in cell to define Vp
    vp_final = interpolate_trilinear(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8)

    ! Vs
    vs1 = vs_tomography(imat,p0+1)
    vs2 = vs_tomography(imat,p1+1)
    vs3 = vs_tomography(imat,p2+1)
    vs4 = vs_tomography(imat,p3+1)
    vs5 = vs_tomography(imat,p4+1)
    vs6 = vs_tomography(imat,p5+1)
    vs7 = vs_tomography(imat,p6+1)
    vs8 = vs_tomography(imat,p7+1)
    ! use trilinear interpolation in cell to define Vs
    vs_final = interpolate_trilinear(vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8)

    ! impose minimum and maximum velocity if needed

    if (vp_final < VP_MIN(imat)) vp_final = VP_MIN(imat)
    if (vp_final > VP_MAX(imat)) vp_final = VP_MAX(imat)

    if (vs_final < VS_MIN(imat)) vs_final = VS_MIN(imat)
    if (vs_final > VS_MAX(imat)) vs_final = VS_MAX(imat)

    ! model parameters for the associated negative imaterial_id index in materials file
    vp_model = vp_final
    vs_model = vs_final

  endif

  ! attenuation
  if (tomo_has_q_values(imat)) then
    qp1 = qp_tomography(imat,p0+1)
    qp2 = qp_tomography(imat,p1+1)
    qp3 = qp_tomography(imat,p2+1)
    qp4 = qp_tomography(imat,p3+1)
    qp5 = qp_tomography(imat,p4+1)
    qp6 = qp_tomography(imat,p5+1)
    qp7 = qp_tomography(imat,p6+1)
    qp8 = qp_tomography(imat,p7+1)
    ! use trilinear interpolation in cell
    qp_final = interpolate_trilinear(qp1,qp2,qp3,qp4,qp5,qp6,qp7,qp8)

    qs1 = qs_tomography(imat,p0+1)
    qs2 = qs_tomography(imat,p1+1)
    qs3 = qs_tomography(imat,p2+1)
    qs4 = qs_tomography(imat,p3+1)
    qs5 = qs_tomography(imat,p4+1)
    qs6 = qs_tomography(imat,p5+1)
    qs7 = qs_tomography(imat,p6+1)
    qs8 = qs_tomography(imat,p7+1)
    ! use trilinear interpolation in cell
    qs_final = interpolate_trilinear(qs1,qs2,qs3,qs4,qs5,qs6,qs7,qs8)

    ! attenuation zero (means negligible attenuation)
    if (qs_final <= 1.e-5) qs_final = ATTENUATION_COMP_MAXIMUM
    if (qp_final <= 1.e-5) qp_final = ATTENUATION_COMP_MAXIMUM

    ! Anderson & Hart (1978) conversion between (Qp,Qs) and (Qkappa,Qmu)
    ! factor L
    L_val = 4.0/3.0 * (vs_model/vp_model)**2

    ! shear attenuation
    qmu_atten = qs_final

    ! converts to bulk attenuation
    if (abs(qs_final - L_val * qp_final) <= 1.e-5) then
      ! negligible bulk attenuation
      qkappa_atten = ATTENUATION_COMP_MAXIMUM
    else
      qkappa_atten = (1.0 - L_val) * qp_final * qs_final / (qs_final - L_val * qp_final)
    endif

    ! attenuation zero (means negligible attenuation)
    if (qmu_atten <= 1.e-5) qmu_atten = ATTENUATION_COMP_MAXIMUM
    if (qkappa_atten <= 1.e-5) qkappa_atten = ATTENUATION_COMP_MAXIMUM

    ! limits Q values
    if (qmu_atten < 1.0d0) qmu_atten = 1.0d0
    if (qmu_atten > ATTENUATION_COMP_MAXIMUM) qmu_atten = ATTENUATION_COMP_MAXIMUM
    if (qkappa_atten < 1.0d0) qkappa_atten = 1.0d0
    if (qkappa_atten > ATTENUATION_COMP_MAXIMUM) qkappa_atten = ATTENUATION_COMP_MAXIMUM

  else
    ! attenuation: arbitrary value, see maximum in constants.h
    qmu_atten = ATTENUATION_COMP_MAXIMUM
    ! Q_Kappa is not implemented in this model_tomography routine yet, thus set it to dummy value
    qkappa_atten = ATTENUATION_COMP_MAXIMUM
  endif

  ! value found
  has_tomo_value = .true.

  contains

  function interpolate_trilinear(val1,val2,val3,val4,val5,val6,val7,val8)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL) :: interpolate_trilinear
  real(kind=CUSTOM_REAL),intent(in) :: val1,val2,val3,val4,val5,val6,val7,val8

  ! note: we use gamma factors from parent routine (with 'contains' we can still use the scope of the parent routine).
  !       gamma parameters are global entities here, and used for briefty of the calling routine command,
  !       just to be aware...

  ! interpolation rule
  interpolate_trilinear =  &
         val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z1) + &
         val2 * gamma_interp_x        * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z2) + &
         val3 * gamma_interp_x        * gamma_interp_y        * (1.d0-gamma_interp_z3) + &
         val4 * (1.d0-gamma_interp_x) * gamma_interp_y        * (1.d0-gamma_interp_z4) + &
         val5 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * gamma_interp_z1 + &
         val6 * gamma_interp_x        * (1.d0-gamma_interp_y) * gamma_interp_z2 + &
         val7 * gamma_interp_x        * gamma_interp_y        * gamma_interp_z3 + &
         val8 * (1.d0-gamma_interp_x) * gamma_interp_y        * gamma_interp_z4

  end function interpolate_trilinear

  end subroutine model_tomography

!
!-------------------------------------------------------------------------------------------------
!

  subroutine deallocate_tomography_files()

    use model_tomography_par

    implicit none

    ! deallocates models dimensions
    deallocate(ORIG_X,ORIG_Y,ORIG_Z)
    deallocate(SPACING_X,SPACING_Y,SPACING_Z)

    ! deallocates models parameter records
    if (ANISOTROPY .and. IMODEL == IMODEL_TOMO) then
      !anisotropy
      deallocate(c_tomography)
    else
      deallocate(vp_tomography)
      deallocate(vs_tomography)
    endif

    deallocate(rho_tomography)
    deallocate(z_tomography)

    ! deallocates models entries
    deallocate(NX,NY,NZ)
    deallocate(nrecord)

    ! deallocates models min/max statistics
    deallocate(VP_MIN,VS_MIN,RHO_MIN)
    deallocate(VP_MAX,VS_MAX,RHO_MAX)

    ! q values
    if (any(tomo_has_q_values)) then
      deallocate(qp_tomography)
      deallocate(qs_tomography)
    endif
    deallocate(tomo_has_q_values)

  end subroutine deallocate_tomography_files
