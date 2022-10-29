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

!> Module dealing with SEP model files.
!! Constraints:
!! * Only acoustic and elastic elements
!! * Requires VP, VS and RHO models
module model_sep_mod

  implicit none

contains

!==============================================================================
!> Reads a SEP elastic model, with vp, vs and rho files.

  subroutine model_sep()

  use generate_databases_par, only: &
    SEP_MODEL_DIRECTORY, FOUR_THIRDS, myrank, IMAIN

  use create_regions_mesh_ext_par, only: &
    rhostore, rho_vp, rho_vs, &
    kappastore, mustore, &
    rho_vpI, rho_vsI, rhoarraystore

  implicit none
  real(kind=4), allocatable, dimension(:,:,:) :: vp_sep, vs_sep, rho_sep
  integer :: NX, NY, NZ
  real :: OX, OY, OZ, DX, DY, DZ
  integer :: NX_alt, NY_alt, NZ_alt,ier
  real :: OX_alt, OY_alt, OZ_alt, DX_alt, DY_alt, DZ_alt
  character(len=512) :: sep_header_name_vp, sep_header_name_vs, &
                        sep_header_name_rho
  character(len=512) :: sep_bin_vp, sep_bin_vs, sep_bin_rho
  real :: xmin, ymin
  integer :: imin, jmin, kmin, ni, nj, nk
  logical :: vp_exists, vs_exists, rho_exists

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using SEP model from directory: ',trim(SEP_MODEL_DIRECTORY)
  endif

  ! Get files from SEP_MODEL_DIRECTORY
  sep_header_name_vp = trim(SEP_MODEL_DIRECTORY) // "/vp.H"
  sep_header_name_vs = trim(SEP_MODEL_DIRECTORY) // "/vs.H"
  sep_header_name_rho = trim(SEP_MODEL_DIRECTORY) // "/rho.H"

  inquire(file=trim(sep_header_name_vp), exist=vp_exists)
  if (.not. vp_exists) stop "SEP vp model should exist"
  inquire(file=trim(sep_header_name_vs), exist=vs_exists)
  if (.not. vp_exists) stop "SEP vp model should exist"
  inquire(file=trim(sep_header_name_rho), exist=rho_exists)
  if (.not. vp_exists) stop "SEP vp model should exist"

  ! Read SEP header for each SEP file
  !    -> n1 (NX), n2 (NY), n3 (NZ)
  !    -> o1 (O1), o2 (O2), o3 (O3)
  !    -> d1 (D1), d2 (D2), d3 (D3)
  !    -> in (filename)
  !    -> Check for unit* to be meter
  ! Parse only one of the header, vp is the most likely to be present
  ! It might be useful to make sure that all SEP files have coherent dimensions.
  call parse_sep_header(trim(sep_header_name_vp) // char(0), &
                        NX, NY, NZ, OX, OY, OZ, DX, DY, DZ, &
                        sep_bin_vp)
  if (vs_exists) then
    call parse_sep_header(trim(sep_header_name_vs) // char(0), &
                          NX_alt, NY_alt, NZ_alt, &
                          OX_alt, OY_alt, OZ_alt, &
                          DX_alt, DY_alt, DZ_alt, &
                          sep_bin_vs)
    if ((NX /= NX_alt) .and. (NY /= NX_alt) .and. (NZ /= NZ_alt) .and. &
        (OX /= OX_alt) .and. (OY /= OX_alt) .and. (OZ /= OZ_alt) .and. &
        (DX /= DX_alt) .and. (DY /= DX_alt) .and. (DZ /= DZ_alt)) then
      stop "SEP headers should be consistant."
    endif
  endif
  if (rho_exists) then
    call parse_sep_header(trim(sep_header_name_rho) // char(0), &
                          NX_alt, NY_alt, NZ_alt, &
                          OX_alt, OY_alt, OZ_alt, &
                          DX_alt, DY_alt, DZ_alt, &
                          sep_bin_rho)
    if ((NX /= NX_alt) .and. (NY /= NX_alt) .and. (NZ /= NZ_alt) .and. &
        (OX /= OX_alt) .and. (OY /= OX_alt) .and. (OZ /= OZ_alt) .and. &
        (DX /= DX_alt) .and. (DY /= DX_alt) .and. (DZ /= DZ_alt)) then
      stop "SEP headers should be consistant."
    endif
  endif

  ! Match slice topography and parts of the SEP file.
  call find_slice_bounds_sep(NX, NY, NZ, OX, OY, OZ, DX, DY, DZ, &
                             xmin, ymin, imin, jmin, kmin, ni, nj, nk)

  !---------.
  ! Read VP |
  !---------'
  ! Read available SEP files, assign default values for unfound files.
  allocate(vp_sep(ni, nj, NZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 625')
  vp_sep(:,:,:) = 0.0

  call read_sep_binary_mpiio(trim(SEP_MODEL_DIRECTORY) // "/" // sep_bin_vp, &
                             NX, NY, NZ, ni, nj, NZ, &
                             imin, jmin, kmin, vp_sep)
  ! Interpolate SEP values on meshfem mesh.
  rho_vp = 0.0
  call interpolate_sep_on_mesh(vp_sep, xmin, ymin, ni, nj, NZ, &
                               DX, DY, DZ, rho_vp)

  !---------.
  ! Read VS |
  !---------'
  if (vs_exists) then
    allocate(vs_sep(ni, nj, NZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 626')
    vs_sep(:,:,:) = 0.0

    call read_sep_binary_mpiio(trim(SEP_MODEL_DIRECTORY) // "/" // sep_bin_vs, &
                               NX, NY, NZ, ni, nj, NZ, &
                               imin, jmin, kmin, vs_sep)
    call interpolate_sep_on_mesh(vs_sep, xmin, ymin, ni, nj, NZ, &
                                 DX, DY, DZ, rho_vs)
  endif

  !----------.
  ! Read RHO |
  !----------'
  if (rho_exists) then
    allocate(rho_sep(ni, nj, NZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 627')
    rho_sep(:,:,:) = 0.0

    call read_sep_binary_mpiio(trim(SEP_MODEL_DIRECTORY) // "/" // sep_bin_rho, &
                               NX, NY, NZ, ni, nj, NZ, &
                               imin, jmin, kmin, rho_sep)
    call interpolate_sep_on_mesh(rho_sep, xmin, ymin, ni, nj, NZ, &
                                 DX, DY, DZ, rhostore)
  endif

  ! Check that values around the acoustic / elastic interface is correct.
  ! Apply eventual correction.
  ! Note that an element should be entirely elastic or acoustic
  if (vs_exists) then
    call correct_sep_interface()
  endif ! vs_exists

  kappastore = rhostore *( rho_vp * rho_vp - FOUR_THIRDS*rho_vs*rho_vs )
  mustore = rhostore*rho_vs*rho_vs

  ! SPECFEM expects rho*velocity
  rho_vp = rho_vp * rhostore
  if (vs_exists) then
    rho_vs = rho_vs * rhostore
  endif
  ! Stacey, a completer par la suite
  !
  rho_vpI = rho_vp
  rho_vsI = rho_vs
  rhoarraystore(1, :, :, :, :) = rhostore
  rhoarraystore(2, :, :, :, :) = rhostore

  if (allocated(vp_sep)) deallocate(vp_sep)
  if (allocated(vs_sep)) deallocate(vs_sep)
  if (allocated(rho_sep)) deallocate(rho_sep)

  end subroutine model_sep


!==============================================================================
!> \brief Parallel read of a SEP file.
!!
!! \param filename Name of the SEP file to read
!! \param NX Dimension in x for the SEP file (n1 in SEP header)
!! \param NY Dimension in y for the SEP file (n1 in SEP header)
!! \param NZ Dimension in z for the SEP file (n1 in SEP header)
!! \param ni Dimension in x for a process.
!! \param nj Dimension in y for a process.
!! \param nk Dimension in z for a process -- usually equal to NZ.
!!
!! \note xyz SEP organized files are expected.

  subroutine read_sep_binary_mpiio(filename, NX, NY, NZ, ni, nj, nk, &
                                   imin, jmin, kmin, var)
  use my_mpi

  implicit none
  ! parameters
  character(len=*), intent(in) :: filename
  integer, intent(in) :: NX, NY, NZ, ni, nj, nk, imin, jmin, kmin
  real(kind=4), dimension(:, :, :), intent(inout) :: var

  ! variables
  integer :: fh
  integer, dimension(3) :: global_sizes, local_sizes, starting_idx
  integer(kind=MPI_OFFSET_KIND) :: displ
  integer :: subarray_type
  integer :: ier
  integer :: status(MPI_STATUS_SIZE)

  ! sizes for the subarrays
  global_sizes = (/ NX, NY, NZ /)
  local_sizes = (/ ni, nj, nk /)
  starting_idx = (/ imin-1, jmin-1, kmin-1 /)

  call synchronize_all()

  ! Create the 3D stencil for values consecutive in X but disjoint in Y and Z.
  call MPI_Type_create_subarray(3, global_sizes, local_sizes, starting_idx, &
                                MPI_ORDER_Fortran, MPI_REAL, subarray_type, ier)
  call MPI_Type_commit(subarray_type, ier)

  call MPI_File_open(my_local_mpi_comm_world, trim(adjustl(filename)), &
                     MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ier)

  ! View is set to match the stencil
  displ = 0
  call MPI_File_set_view(fh, displ, MPI_REAL, subarray_type, "native", &
                         MPI_INFO_NULL, ier)
  ! Number of elements is equal to the local size
  call MPI_File_read_all(fh, var, ni * nj * nk, MPI_REAL, status, ier)
  call MPI_File_close(fh, ier)

  end subroutine read_sep_binary_mpiio

!==============================================================================
!> \brief Assign values read from SEP files to GLL points
!!
!! \note from Yang Luo's external routines. Actually finds a nearby point.

  subroutine interpolate_sep_on_mesh(sep_var, xmin, ymin, ni, nj, NZ, &
                                   DX, DY, DZ, var)

  use generate_databases_par, only: NGLLX, NGLLY, NGLLZ, NSPEC => NSPEC_AB, &
                                    ibool, xstore, ystore, zstore, &
                                    CUSTOM_REAL
  implicit none
  !--- Parameters
  real(kind=4), dimension(:,:,:), intent(in) :: sep_var
  real, intent(in) :: xmin, ymin, DX, DY, DZ
  integer, intent(in) :: ni, nj, NZ
  real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: var
  !--- Variables
  integer :: iglob, i_target, j_target, k_target
  integer :: ispec, i, j, k

  do ispec=1,NSPEC
     do i=1,NGLLX
        do j=1,NGLLY
           do k=1,NGLLZ
              iglob=ibool(i,j,k,ispec);
              i_target=NINT((xstore(i, j, k, ispec)-xmin)/DX)+1
              j_target=NINT((ystore(i, j, k, ispec)-ymin)/DY)+1
              ! in SEP, z-axis is downward; in SPECFEM, z-axis is upward
              k_target=NINT(-zstore(i, j, k, ispec)/DZ)+1; !NOT using (z_temp-zmin)
              ! Stay in the computational domain
              if (i_target < 1)  i_target=  1
              if (j_target < 1)  j_target=  1
              if (k_target < 1)  k_target=  1
              if (i_target > ni)  i_target= ni
              if (j_target > nj)  j_target= nj
              if (k_target > NZ)  k_target= NZ

              var(i,j,k,ispec) = sep_var(i_target,j_target,k_target)
           enddo
        enddo
     enddo
  enddo

  end subroutine interpolate_sep_on_mesh

!==============================================================================
!> Find offsets and number of elements to read from the SEP file according
!! to the slice topology.

  subroutine find_slice_bounds_sep(NX, NY, NZ, OX, OY, OZ, DX, DY, DZ, &
                                 xmin, ymin, imin, jmin, kmin, ni, nj, nk)

  use generate_databases_par, only: xstore, ystore, zstore

  implicit none
  ! Parameters
  integer, intent(in) :: NX, NY, NZ
  real, intent(in) :: OX, OY, OZ, DX, DY, DZ
  real, intent(out) :: xmin, ymin
  integer, intent(out) :: imin, jmin, kmin, ni, nj, nk
  ! Variables
  real :: xmax, ymax, zmin, zmax
  integer :: imax, jmax, kmax

  ! Actual bounds for the current slice mesh
  xmin=minval(xstore); xmax=maxval(xstore);
  ymin=minval(ystore); ymax=maxval(ystore);
  ! in SEP, z-axis is downward; in SPECFEM, z-axis is upward
  zmin=minval(-zstore); zmax=maxval(-zstore);

  ! Bounds for the SEP model
  imin=floor((xmin-OX)/DX+1); imax=ceiling((xmax-OX)/DX+1);
  jmin=floor((ymin-OY)/DY+1); jmax=ceiling((ymax-OY)/DY+1);
  kmin=floor((zmin-OZ)/DZ+1); kmax=ceiling((zmax-OZ)/DZ+1);
  if (imin < 1)  imin= 1;  if (jmin < 1)  jmin= 1;  if (kmin < 1)  kmin= 1;
  if (imax > NX) imin=NX;  if (jmax > NY) jmax=NY;  if (kmax > NZ) kmax=NZ;
  ! Number of SEP indexes for the current slice
  ni=imax-imin+1
  nj=jmax-jmin+1
  nk=kmax-kmin+1

  ! Corrected bounds for the current slice when looking for SEP indexes
  xmin=DX*(imin-1)+OX; xmax=DX*(imax-1)+OX;
  ymin=DY*(jmin-1)+OY; ymax=DY*(jmax-1)+OY;
  zmin=DZ*(kmin-1)+OZ; zmax=DZ*(kmax-1)+OZ;

  end subroutine find_slice_bounds_sep


!==============================================================================
!> Make sure that each element is fully acoustic or fully elastic.

  subroutine correct_sep_interface()

  use generate_databases_par, only: NGLLX, NGLLY, NGLLZ, NSPEC => NSPEC_AB
  use create_regions_mesh_ext_par, only: rhostore, rho_vp, rho_vs, &
                                         ispec_is_acoustic, ispec_is_elastic

  implicit none
  integer :: ispec, i, j, k
  integer :: num_acoustic_pts, num_elastic_pts

  do ispec = 1, NSPEC
    num_acoustic_pts = 0
    num_elastic_pts = 0
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          if (rho_vs(i, j, k, ispec) == 0) then
            num_acoustic_pts = num_acoustic_pts + 1
          else
            num_elastic_pts = num_elastic_pts + 1
          endif
        enddo ! k
      enddo ! j
    enddo ! i
    if (num_acoustic_pts > num_elastic_pts) then
      ispec_is_acoustic(ispec) = .true.
      ! Z axis is up. Bottom points might include the elastic material
      if (any(rho_vs(:, :, 1, ispec) /= 0.0)) then
        rho_vs(:, :, :, ispec) = 0.0
        rhostore(:, :, :, ispec) = minval(rhostore(:, :, :, ispec))
        rho_vp(:, :, :, ispec) = minval(rho_vp(:, :, :, ispec))
      endif
    else
      ispec_is_elastic(ispec) = .true.
      ! Z axis is up. Top points might include the acoustic interface
      if (any(rho_vs(:, :, NGLLZ, ispec) == 0.0)) then
        rho_vs(:, :, NGLLZ, ispec) = rho_vs(:, :, NGLLZ-1, ispec)
        rhostore(:, :, NGLLZ, ispec) = rhostore(:, :, NGLLZ-1, ispec)
        rho_vp(:, :, NGLLZ, ispec) = rho_vp(:, :, NGLLZ-1, ispec)
      endif
    endif
  enddo ! ispec

  end subroutine correct_sep_interface

end module model_sep_mod

