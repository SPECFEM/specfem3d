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


  subroutine cmm_determine_cavity(nglob)

  use constants, only: MF_IN_DATA_FILES,MAX_STRING_LEN,IMAIN,HUGEVAL,TINYVAL,NDIM,myrank
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use create_meshfem_par, only: nodes_coords,ispec_material_id,iboun,iMPIcut_xi,iMPIcut_eta

  use meshfem_par, only: ibool,xstore,ystore,zstore,nspec,NPROC_XI,NPROC_ETA,CAVITY_FILE

  implicit none

  integer,intent(inout) :: nglob

  ! local parameters
  !-------------------------------cavity----------------------------------------
  integer :: ncavity
  integer :: i,j,k,ier,iproc
  integer :: i_cavity,i_node,inode_new_mesh
  integer :: ispec,ispec_new_mesh
  integer :: nspec_old,nglob_old

  double precision :: x0,x1,y0,y1,z0,z1,xmid,ymid,zmid

  double precision,allocatable,dimension(:) :: cavity_x0,cavity_x1,cavity_y0,cavity_y1,cavity_z0,cavity_z1

  logical,allocatable :: is_elmt(:),is_node(:)
  logical,allocatable :: iboun_old(:,:)
  logical, dimension(:,:), allocatable :: iMPIcut_xi_old,iMPIcut_eta_old
  logical :: cavity_file_exists,in_cavity

  integer,allocatable :: ibool_old(:,:,:,:),ispec_material_id_old(:)
  integer,allocatable :: ispec_new(:),inode_new(:)

  double precision,allocatable :: nodes_coords_old(:,:)

  integer :: cavity_num_elements(1),cavity_num_elements_all(0:NPROC_ETA*NPROC_XI-1)
  integer :: num_cav_total,icav,icav_glob
  double precision,allocatable,dimension(:,:) :: cavity_boundary
  double precision,allocatable,dimension(:,:) :: tmp_all

  character(len=MAX_STRING_LEN) :: filename

  !-------------------------------cavity----------------------------------------

  ! begin cavity
  ! default
  ncavity = 0

  ! read cavity file
  filename = trim(MF_IN_DATA_FILES)//trim(CAVITY_FILE)
  open(111,file=filename,action='read',status='old',iostat=ier)
  if (ier /= 0) then
    cavity_file_exists = .false.
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*)'File "'//trim(filename)//'" not found: assume no cavity'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  else
    cavity_file_exists = .true.
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Using cavity file: ',trim(filename)
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! check if the file is blank
  if (cavity_file_exists) then
    ! skip one comment line
    read(111,*)
    read(111,*) ncavity

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'cavity:'
      write(IMAIN,*) '  number of cavities = ',ncavity
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    !! checks that only 1 entry, more not supported yet...
    !if (ncavity > 1) then
    !  stop 'Error: only 1 cavity supported so far! Please check your cavity file...'
    !endif

    ! reads in cavity dimensions
    if (ncavity > 0) then
      allocate(cavity_x0(ncavity),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1321')
      allocate(cavity_x1(ncavity),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1322')
      allocate(cavity_y0(ncavity),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1323')
      allocate(cavity_y1(ncavity),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1324')
      allocate(cavity_z0(ncavity),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1325')
      allocate(cavity_z1(ncavity),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1326')
      cavity_x0=HUGEVAL; cavity_x1=HUGEVAL
      cavity_y0=HUGEVAL; cavity_y1=HUGEVAL
      cavity_z0=HUGEVAL; cavity_z1=HUGEVAL
      ! skip one comment line
      read(111,*)
      !read cavity range
      do i_cavity = 1,ncavity
        read(111,*) cavity_x0(i_cavity),cavity_x1(i_cavity), &
                    cavity_y0(i_cavity),cavity_y1(i_cavity), &
                    cavity_z0(i_cavity),cavity_z1(i_cavity)

        ! user output
        if (myrank == 0) then
          write(IMAIN,*)'  cavity range: x min / max = ',cavity_x0(i_cavity),cavity_x1(i_cavity),'(m)'
          write(IMAIN,*)'                y min / max = ',cavity_y0(i_cavity),cavity_y1(i_cavity),'(m)'
          write(IMAIN,*)'                z min / max = ',cavity_z0(i_cavity),cavity_z1(i_cavity),'(m)'
          write(IMAIN,*)
          call flush_IMAIN()
        endif
      enddo
    endif
    ! close cavity file
    close(111)
  endif

  ! add cavity if necessary
  if (ncavity > 0) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  creating cavity'
      call flush_IMAIN()
    endif

    allocate(is_elmt(nspec),is_node(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1327')
    if (ier /= 0) stop 'Error allocating is_elmt, is_node arrays'

    is_elmt(:) = .true.
    is_node(:) = .false.

    do ispec = 1,nspec
      ! find mid point of the spectral element
      x0 = xstore(1,1,1,ispec)
      y0 = ystore(1,1,1,ispec)
      z0 = zstore(1,1,1,ispec)

      x1 = xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      y1 = ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      z1 = zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)

      xmid = 0.5d0*(x0+x1)
      ymid = 0.5d0*(y0+y1)
      zmid = 0.5d0*(z0+z1)

      in_cavity = .false.
cavity: do i_cavity = 1,ncavity
        ! checks if midpoint within cavity range
        if ((xmid >= cavity_x0(i_cavity) .and. xmid <= cavity_x1(i_cavity)) .and. &
            (ymid >= cavity_y0(i_cavity) .and. ymid <= cavity_y1(i_cavity)) .and. &
            (zmid >= cavity_z0(i_cavity) .and. zmid <= cavity_z1(i_cavity))) then
          ! deactivate spectral element
          is_elmt(ispec) = .false.
          in_cavity = .true.
          exit cavity
        endif
      enddo cavity

      if (.not. in_cavity) then
        ! intact
        ! activate nodes
        do k = 1,NGLLZ_M
          do j = 1,NGLLY_M
            do i = 1,NGLLX_M
              is_node(ibool(i,j,k,ispec)) = .true.
            enddo
          enddo
        enddo
      endif
    enddo ! ispec=1,nspec

    deallocate(cavity_x0,cavity_x1,cavity_y0,cavity_y1,cavity_z0,cavity_z1)

    nspec_old = nspec
    nglob_old = nglob
    nspec = count(is_elmt(:))
    nglob = count(is_node(:))

    ! determines number of cavity elements found
    cavity_num_elements(1) = nspec_old - nspec

    ! collects on all processes
    call gather_all_all_i(cavity_num_elements, 1, cavity_num_elements_all, 1, NPROC_ETA*NPROC_XI)
    num_cav_total = sum(cavity_num_elements_all(:))

    ! user output
    if (myrank == 0) then
      do iproc = 0,NPROC_ETA*NPROC_XI - 1
        write(IMAIN,*)'    found ',cavity_num_elements_all(iproc),'cavity elements in slice ',iproc
      enddo
      write(IMAIN,*)
      write(IMAIN,*) '    total cavity elements found = ',num_cav_total
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! handles effect on MPI boundary:
    !   this will exchange all cavity surfaces (determined by midpoints) to all processes and
    !   check with their MPI boundary surfaces, if an MPI surface will have to be removed from the list
    if (NPROC_ETA*NPROC_XI > 1) then
      ! note: index (0,*) == 1 indicates a boundary point
      !       and there can be 4 boundaries maximum for each element: xi-min, xi-max, eta-min, eta-max side
      allocate(cavity_boundary(0:3,4*num_cav_total),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1328')
      if (ier /= 0) stop 'Error allocating cavity_boundary arrays'
      cavity_boundary(:,:) = 0.0

      ! checks if cavity elements are on a MPI boundary
      icav = 0
      do ispec = 1,nspec_old
        if (.not. is_elmt(ispec)) then
          ! local cavity element, will be removed
          icav = icav + 1

          ! checks if on an MPI boundary
          ! note: mesher only stored element corner points NGLLX_M == NGLLY_M == NGLLZ_M == 2
          ! xi-min side
          if (iMPIcut_xi(1,ispec) .eqv. .true.) then
            ! surface corners: (see save_databases.F90)
            ! ibool(1,1,1,ispec) , ibool(1,2,1,ispec) , ibool(1,1,2,ispec) , ibool(1,2,2,ispec)

            ! surface mid-point
            x1 = 0.25*(xstore(1,1,1,ispec) + xstore(1,NGLLY_M,1,ispec)+ &
                       xstore(1,1,NGLLZ_M,ispec) + xstore(1,NGLLY_M,NGLLZ_M,ispec))
            y1 = 0.25*(ystore(1,1,1,ispec) + ystore(1,NGLLY_M,1,ispec)+ &
                       ystore(1,1,NGLLZ_M,ispec) + ystore(1,NGLLY_M,NGLLZ_M,ispec))
            z1 = 0.25*(zstore(1,1,1,ispec) + zstore(1,NGLLY_M,1,ispec)+ &
                       zstore(1,1,NGLLZ_M,ispec) + zstore(1,NGLLY_M,NGLLZ_M,ispec))

            ! index in total array (counts from rank 0 up)
            icav_glob = icav
            if (NPROC_XI*NPROC_ETA > 1 .and. myrank > 0) icav_glob = icav_glob + sum(cavity_num_elements_all(0:myrank-1)) * 4
            cavity_boundary(0,icav_glob) = 1.0 ! flag to indicate being on a boundary
            cavity_boundary(1,icav_glob) = x1
            cavity_boundary(2,icav_glob) = y1
            cavity_boundary(3,icav_glob) = z1

            ! removes element from MPI interface
            iMPIcut_xi(1,ispec) = .false.
          endif

          ! xi-max side
          if (iMPIcut_xi(2,ispec) .eqv. .true.) then
            ! surface corners: (see save_databases.F90, line 332)
            ! ibool(2,1,1,ispec) , ibool(2,2,1,ispec) , ibool(2,1,2,ispec) , ibool(2,2,2,ispec)

            ! surface mid-point
            x1 = 0.25*(xstore(NGLLX_M,1,1,ispec) + xstore(NGLLX_M,NGLLY_M,1,ispec)+ &
                       xstore(NGLLX_M,1,NGLLZ_M,ispec) + xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec))
            y1 = 0.25*(ystore(NGLLX_M,1,1,ispec) + ystore(NGLLX_M,NGLLY_M,1,ispec)+ &
                       ystore(NGLLX_M,1,NGLLZ_M,ispec) + ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec))
            z1 = 0.25*(zstore(NGLLX_M,1,1,ispec) + zstore(NGLLX_M,NGLLY_M,1,ispec)+ &
                       zstore(NGLLX_M,1,NGLLZ_M,ispec) + zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec))

            ! index in total array (counts from rank 0 up)
            icav_glob = icav + 1
            if (NPROC_XI*NPROC_ETA > 1 .and. myrank > 0) icav_glob = icav_glob + sum(cavity_num_elements_all(0:myrank-1)) * 4
            cavity_boundary(0,icav_glob) = 1.0 ! flag to indicate being on a boundary
            cavity_boundary(1,icav_glob) = x1
            cavity_boundary(2,icav_glob) = y1
            cavity_boundary(3,icav_glob) = z1

            ! removes element from MPI interface
            iMPIcut_xi(2,ispec) = .false.
          endif

          ! eta-min side
          if (iMPIcut_eta(1,ispec) .eqv. .true.) then
            ! surface corners: (see save_databases.F90)
            ! ibool(1,1,1,ispec),ibool(2,1,1,ispec),ibool(1,1,2,ispec),ibool(2,1,2,ispec)

            ! surface mid-point
            x1 = 0.25*(xstore(1,1,1,ispec) + xstore(NGLLX_M,1,1,ispec) + &
                       xstore(1,1,NGLLZ_M,ispec) + xstore(NGLLX_M,1,NGLLZ_M,ispec))
            y1 = 0.25*(ystore(1,1,1,ispec) + ystore(NGLLX_M,1,1,ispec) + &
                       ystore(1,1,NGLLZ_M,ispec) + ystore(NGLLX_M,1,NGLLZ_M,ispec))
            z1 = 0.25*(zstore(1,1,1,ispec) + zstore(NGLLX_M,1,1,ispec) + &
                       zstore(1,1,NGLLZ_M,ispec) + zstore(NGLLX_M,1,NGLLZ_M,ispec))

            ! index in total array (counts from rank 0 up)
            icav_glob = icav + 2
            if (NPROC_XI*NPROC_ETA > 1 .and. myrank > 0) icav_glob = icav_glob + sum(cavity_num_elements_all(0:myrank-1)) * 4
            cavity_boundary(0,icav_glob) = 1.0 ! flag to indicate being on a boundary
            cavity_boundary(1,icav_glob) = x1
            cavity_boundary(2,icav_glob) = y1
            cavity_boundary(3,icav_glob) = z1

            ! removes element from MPI interface
            iMPIcut_eta(1,ispec) = .false.
          endif

          ! eta-max side
          if (iMPIcut_eta(2,ispec) .eqv. .true.) then
            ! surface corners: (see save_databases.F90)
            ! ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)

            ! surface mid-point
            x1 = 0.25*(xstore(NGLLX_M,NGLLY_M,1,ispec) + xstore(1,NGLLY_M,1,ispec) + &
                       xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) + xstore(1,NGLLY_M,NGLLZ_M,ispec))
            y1 = 0.25*(ystore(NGLLX_M,NGLLY_M,1,ispec) + ystore(1,NGLLY_M,1,ispec) + &
                       ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) + ystore(1,NGLLY_M,NGLLZ_M,ispec))
            z1 = 0.25*(zstore(NGLLX_M,NGLLY_M,1,ispec) + zstore(1,NGLLY_M,1,ispec) + &
                       zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) + zstore(1,NGLLY_M,NGLLZ_M,ispec))

            ! index in total array (counts from rank 0 up)
            icav_glob = icav + 3
            if (NPROC_XI*NPROC_ETA > 1 .and. myrank > 0) icav_glob = icav_glob + sum(cavity_num_elements_all(0:myrank-1)) * 4
            cavity_boundary(0,icav_glob) = 1.0 ! flag to indicate being on a boundary
            cavity_boundary(1,icav_glob) = x1
            cavity_boundary(2,icav_glob) = y1
            cavity_boundary(3,icav_glob) = z1

            ! removes element from MPI interface
            iMPIcut_eta(2,ispec) = .false.
          endif
        endif
      enddo
      if (icav /= cavity_num_elements(1)) stop 'Invalid cavity elements in loop'

      !print *,'cavity boundary:',myrank,'array:',cavity_boundary(:,:)

      ! collects on main processes
      if (myrank == 0) then
        allocate(tmp_all(4,num_cav_total*4),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1329')
      else
        allocate(tmp_all(1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1330')
      endif
      call sum_all_1Darray_dp(cavity_boundary,tmp_all,size(cavity_boundary))
      if (myrank == 0) then
        cavity_boundary(:,:) = tmp_all(:,:)
      endif
      deallocate(tmp_all)
      ! broadcasts to all others
      call bcast_all_dp(cavity_boundary,size(cavity_boundary))

      !print *,'cavity boundary after:',myrank,'array:',cavity_boundary(:,:)

      ! note: this is not well optimized, but a basic routine open for improvement...
      ! checks if any cavity elements are on an MPI boundary
      do ispec = 1,nspec_old
        ! checks if additional points from cavity-boundary are on an MPI boundary
        ! xi-min side
        if (iMPIcut_xi(1,ispec) .eqv. .true.) then
          ! surface corners: (see save_databases.F90)
          ! ibool(1,1,1,ispec) , ibool(1,2,1,ispec) , ibool(1,1,2,ispec) , ibool(1,2,2,ispec)

          ! surface mid-point
          x1 = 0.25*(xstore(1,1,1,ispec) + xstore(1,NGLLY_M,1,ispec) + &
                     xstore(1,1,NGLLZ_M,ispec) + xstore(1,NGLLY_M,NGLLZ_M,ispec))
          y1 = 0.25*(ystore(1,1,1,ispec) + ystore(1,NGLLY_M,1,ispec) + &
                     ystore(1,1,NGLLZ_M,ispec) + ystore(1,NGLLY_M,NGLLZ_M,ispec))
          z1 = 0.25*(zstore(1,1,1,ispec) + zstore(1,NGLLY_M,1,ispec) + &
                     zstore(1,1,NGLLZ_M,ispec) + zstore(1,NGLLY_M,NGLLZ_M,ispec))

          ! loops over all cavity surface points
          do icav = 1, 4*num_cav_total
            if (cavity_boundary(0,icav) > 1.e-5) then ! flag to indicate being on a MPI boundary
              if (abs(cavity_boundary(1,icav) - x1) < TINYVAL .and. &
                  abs(cavity_boundary(2,icav) - y1) < TINYVAL .and. &
                  abs(cavity_boundary(3,icav) - z1) < TINYVAL) then
                    ! removes element from MPI interface
                    iMPIcut_xi(1,ispec) = .false.
              endif
            endif
          enddo
        endif

        ! xi-max side
        if (iMPIcut_xi(2,ispec) .eqv. .true.) then
          ! surface corners: (see save_databases.F90, line 332)
          ! ibool(2,1,1,ispec) , ibool(2,2,1,ispec) , ibool(2,1,2,ispec) , ibool(2,2,2,ispec)

          ! surface mid-point
          x1 = 0.25*(xstore(NGLLX_M,1,1,ispec) + xstore(NGLLX_M,NGLLY_M,1,ispec) + &
                     xstore(NGLLX_M,1,NGLLZ_M,ispec) + xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec))
          y1 = 0.25*(ystore(NGLLX_M,1,1,ispec) + ystore(NGLLX_M,NGLLY_M,1,ispec) + &
                     ystore(NGLLX_M,1,NGLLZ_M,ispec) + ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec))
          z1 = 0.25*(zstore(NGLLX_M,1,1,ispec) + zstore(NGLLX_M,NGLLY_M,1,ispec) + &
                     zstore(NGLLX_M,1,NGLLZ_M,ispec) + zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec))

          ! loops over all cavity surface points
          do icav = 1, 4*num_cav_total
            if (cavity_boundary(0,icav) > 1.e-5) then ! flag to indicate being on a MPI boundary
              if (abs(cavity_boundary(1,icav) - x1) < TINYVAL .and. &
                  abs(cavity_boundary(2,icav) - y1) < TINYVAL .and. &
                  abs(cavity_boundary(3,icav) - z1) < TINYVAL) then
                    ! removes element from MPI interface
                    iMPIcut_xi(2,ispec) = .false.
              endif
            endif
          enddo
        endif

        ! eta-min side
        if (iMPIcut_eta(1,ispec) .eqv. .true.) then
          ! surface corners: (see save_databases.F90)
          ! ibool(1,1,1,ispec),ibool(2,1,1,ispec),ibool(1,1,2,ispec),ibool(2,1,2,ispec)

          ! surface mid-point
          x1 = 0.25*(xstore(1,1,1,ispec) + xstore(NGLLX_M,1,1,ispec) + &
                     xstore(1,1,NGLLZ_M,ispec) + xstore(NGLLX_M,1,NGLLZ_M,ispec))
          y1 = 0.25*(ystore(1,1,1,ispec) + ystore(NGLLX_M,1,1,ispec) + &
                     ystore(1,1,NGLLZ_M,ispec) + ystore(NGLLX_M,1,NGLLZ_M,ispec))
          z1 = 0.25*(zstore(1,1,1,ispec) + zstore(NGLLX_M,1,1,ispec) + &
                     zstore(1,1,NGLLZ_M,ispec) + zstore(NGLLX_M,1,NGLLZ_M,ispec))

          ! loops over all cavity surface points
          do icav = 1, 4*num_cav_total
            if (cavity_boundary(0,icav) > 1.e-5) then ! flag to indicate being on a MPI boundary
              if (abs(cavity_boundary(1,icav) - x1) < TINYVAL .and. &
                  abs(cavity_boundary(2,icav) - y1) < TINYVAL .and. &
                  abs(cavity_boundary(3,icav) - z1) < TINYVAL) then
                    ! removes element from MPI interface
                    iMPIcut_eta(1,ispec) = .false.
              endif
            endif
          enddo
        endif

        ! eta-max side
        if (iMPIcut_eta(2,ispec) .eqv. .true.) then
          ! surface corners: (see save_databases.F90)
          ! ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)

          ! surface mid-point
          x1 = 0.25*(xstore(NGLLX_M,NGLLY_M,1,ispec) + xstore(1,NGLLY_M,1,ispec) + &
                     xstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) + xstore(1,NGLLY_M,NGLLZ_M,ispec))
          y1 = 0.25*(ystore(NGLLX_M,NGLLY_M,1,ispec) + ystore(1,NGLLY_M,1,ispec) + &
                     ystore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) + ystore(1,NGLLY_M,NGLLZ_M,ispec))
          z1 = 0.25*(zstore(NGLLX_M,NGLLY_M,1,ispec) + zstore(1,NGLLY_M,1,ispec) + &
                     zstore(NGLLX_M,NGLLY_M,NGLLZ_M,ispec) + zstore(1,NGLLY_M,NGLLZ_M,ispec))

          ! loops over all cavity surface points
          do icav = 1, 4*num_cav_total
            if (cavity_boundary(0,icav) > 1.e-5) then ! flag to indicate being on a MPI boundary
              if (abs(cavity_boundary(1,icav) - x1) < TINYVAL .and. &
                  abs(cavity_boundary(2,icav) - y1) < TINYVAL .and. &
                  abs(cavity_boundary(3,icav) - z1) < TINYVAL) then
                    ! removes element from MPI interface
                    iMPIcut_eta(2,ispec) = .false.
              endif
            endif
          enddo
        endif
      enddo
      call synchronize_all()
    endif ! NPROC > 1

    ! checks if anything to do further
    if (cavity_num_elements(1) == 0) return

    ! allocates new mesh arrays
    allocate(ispec_new(nspec_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1331')
    if (ier /= 0) stop 'Error allocating ispec_new array'

    ispec_new(:) = -1

    ispec_new_mesh = 0
    do ispec = 1,nspec_old
      if (is_elmt(ispec)) then
        ispec_new_mesh = ispec_new_mesh + 1
        ispec_new(ispec) = ispec_new_mesh
      endif
    enddo
    if (ispec_new_mesh /= nspec) call exit_MPI(myrank,'ERROR: new number of spectral elements mismatch!')

    allocate(inode_new(nglob_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1332')
    inode_new(:) = -1

    inode_new_mesh = 0
    do i_node = 1,nglob_old
      if (is_node(i_node)) then
        inode_new_mesh = inode_new_mesh + 1
        inode_new(i_node) = inode_new_mesh
      endif
    enddo
    if (inode_new_mesh /= nglob) call exit_MPI(myrank,'ERROR: new number of spectral elements mismatch!')

    ! old mesh arrays
    allocate(nodes_coords_old(nglob_old,NDIM),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1333')
    allocate(ispec_material_id_old(nspec_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1334')
    allocate(ibool_old(NGLLX_M,NGLLY_M,NGLLZ_M,nspec_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1335')
    allocate(iboun_old(6,nspec_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1336')
    allocate(iMPIcut_xi_old(2,nspec_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1337')
    allocate(iMPIcut_eta_old(2,nspec_old),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1338')
    if (ier /= 0 ) stop 'Error allocating old mesh arrays for cavity'

    nodes_coords_old(:,:) = nodes_coords(:,:)
    ispec_material_id_old(:) = ispec_material_id(:)
    ibool_old(:,:,:,:) = ibool(:,:,:,:)
    iboun_old (:,:)= iboun(:,:)
    iMPIcut_xi_old(:,:) = iMPIcut_xi(:,:)
    iMPIcut_eta_old(:,:) = iMPIcut_eta(:,:)

    deallocate(nodes_coords)
    deallocate(ispec_material_id)
    deallocate(ibool)
    deallocate(iboun)
    deallocate(iMPIcut_xi,iMPIcut_eta)

    ! re-allocates new mesh arrays
    allocate(nodes_coords(nglob,NDIM),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1339')
    allocate(ispec_material_id(nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1340')
    allocate(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1341')
    allocate(iboun(6,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1342')
    allocate(iMPIcut_xi(2,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1343')
    allocate(iMPIcut_eta(2,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1344')
    if (ier /= 0 ) stop 'Error allocating updated mesh arrays for cavity'

    ! new specs
    do ispec = 1,nspec_old
      if (is_elmt(ispec)) then
        ! new mesh element
        ispec_material_id(ispec_new(ispec)) = ispec_material_id_old(ispec)

        ! absorbing boundary
        iboun(:,ispec_new(ispec)) = iboun_old(:,ispec)

        ! MPI boundary
        iMPIcut_xi(:,ispec_new(ispec)) = iMPIcut_xi_old(:,ispec)
        iMPIcut_eta(:,ispec_new(ispec)) = iMPIcut_eta_old(:,ispec)

        ! activate nodes
        do k = 1,NGLLZ_M
          do j = 1,NGLLY_M
            do i = 1,NGLLX_M
              ibool(i,j,k,ispec_new(ispec)) = inode_new(ibool_old(i,j,k,ispec))
            enddo
          enddo
        enddo
      endif
    enddo

    deallocate(is_elmt)
    deallocate(ispec_new)
    deallocate(ispec_material_id_old)
    deallocate(ibool_old)
    deallocate(iboun_old)
    deallocate(iMPIcut_xi_old,iMPIcut_eta_old)

    ! new coordinates
    do i_node = 1,nglob_old
      if (is_node(i_node)) then
        nodes_coords(inode_new(i_node),:) = nodes_coords_old(i_node,:)
      endif
    enddo

    deallocate(is_node)
    deallocate(inode_new)
    deallocate(nodes_coords_old)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  cavity setup done'
      call flush_IMAIN()
    endif

  endif ! of if (ncavity > 0)

  end subroutine cmm_determine_cavity
