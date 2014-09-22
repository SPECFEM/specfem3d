module model_sep_mode
  implicit none
contains

!==============================================================================
subroutine model_sep()
  use create_regions_mesh_ext_par, only rhostore, rho_vp, rho_vs

  real(kind=4), allocatable, dimension(:,:,:) :: vp_sep, vs_sep, rho_sep
  integer :: NX, NY, NZ
  real :: OX, OY, OZ, DX, DY, DZ
  character(len=512) :: sep_bin_vp
  real :: xmin, ymin
  integer :: imin, jmin, kmin, ni, nj, nk

  ! Read SEP header for each SEP file
  !    -> n1 (NX), n2 (NY), n3 (NZ)
  !    -> o1 (O1), o2 (O2), o3 (O3)
  !    -> d1 (D1), d2 (D2), d3 (D3)
  !    -> in (filename)
  !    -> Check for unit* to be meter
  ! Parse only one of the header, vp is the most likely to be present
  ! It might be useful to make sure that all SEP files have coherent dimensions.
  parse_sep_vs_header(sep_header_name_vp,                 &
                      NX, NY, NZ, OX, OY, OZ, D1, D2, D3, &
                      sep_bin_vp)

  ! Match slice topography and parts of the SEP file.
  call find_slice_bounds_sep(NX, NY, NZ, OX, OY, OZ, DX, DY, DZ, &
                             xmin, ymin, imin, jmin, kmin, ni, nj, nk)

  ! Read available SEP files, assign default values for unfound files.
  allocate(vp_sep(ni, nj, NZ))
  allocate(vs_sep(ni, nj, NZ))
  allocate(rho_sep(ni, nj, NZ))

  call read_sep_binary_mpiio(FILE_VP, NX, NY, NZ, ni, nj, NZ, vp_sep)
  call read_sep_binary_mpiio(FILE_VS, NX, NY, NZ, ni, nj, NZ, vs_sep)
  call read_sep_binary_mpiio(FILE_RHO, NX, NY, NZ, ni, nj, NZ, rho_sep)

  ! Interpolate SEP values on meshfem mesh.
  call interpolate_sep_on_mesh(rho_sep, xmin, ymin, ni, nj, NZ, &
                               DX, DY, DZ, rhostore)
  call interpolate_sep_on_mesh(vp_sep, xmin, ymin, ni, nj, NZ, &
                               DX, DY, DZ, rho_vp)
  call interpolate_sep_on_mesh(vs_sep, xmin, ymin, ni, nj, NZ, &
                               DX, DY, DZ, rho_vs)
  ! SPECFEM expects rho*velocity
  rho_vp = rho_vp * rhostore
  rho_vs = rho_vs * rhostore

  deallocate(vp_sep)
  deallocate(vs_sep)
  deallocate(rho_sep)
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
subroutine read_sep_binary_mpiio(filename, NX, NY, NZ, ni, nj, nk, var)
  use mpi
  ! Parameters
  character(len=*), intent(in) :: filename
  integer, intent(in) :: NX, NY, NZ, ni, nj, nk
  real(kind =4), dimension(:, :, :), intent(inout) :: var

  !Variables
  integer :: fh
  integer, dimension(3) :: global_sizes, local_sizes, starting_idx
  integer(kind=MPI_OFFSET_KIND) :: displ
  integer :: subarray_type
  integer :: mpi_status, ier

  ! Sizes for the subarrays
  global_sizes = (/ NX, NY, NZ /)
  local_sizes = (/ ni, nj, nk /)
  starting_idx = (/ imin-1, jmin-1, kmin-1 /)

  ! Create the 3D stencil for values consecutive in X but disjoint in Y and Z.
  call MPI_Type_create_subarray(3, global_sizes, local_sizes, starting_idx, &
                                MPI_ORDER_FORTRAN, MPI_REAL, subarray_type, ier)
  call MPI_Type_commit(subarray_type, ier)
  
  call MPI_File_open(MPI_COMM_WORLD, trim(adjustl(filename)), &
                     MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ier)
  ! View is set to match the stencil 
  displ = 0
  call MPI_File_set_view(fh, displ, MPI_REAL, subarray_type, "native", &
                         MPI_INFO_NULL, ier)
  ! Number of element is equal to the local size                       
  call MPI_File_read_all(fh, var, ni * nj * nk, &
                         MPI_REAL, mpi_status, ier)
  call MPI_File_close(fh, ier)
end subroutine read_sep_binary_mpiio

!==============================================================================
!> \brief Assign values read from SEP files to GLL points
!!
!! \note from Yang Luo's external routines. Actually finds a nearby point.
subroutine interpolate_sep_on_mesh(sep_var, xmin, ymin, ni, nj, NZ, &
                                   DX, DY, DZ, var)
  use generate_database_par, only: NGLLX, NGLLY, NGLLZ, NSPEC=>NSPEC_AB &
                                   ibool
  !--- Parameters
  real(kind=4), dimension(:,:,:), intent(in) :: sep_var 
  real, intent(in) :: xmin, ymin, DX, DY, DZ
  integer, intent(in) :: ni, nj, NZ
  real, dimension(:,:,:,:) :: var
  !--- Variables
  integer :: iglob, i_target, j_target, k_target
  integer :: ispec, i, j, k                                   

  do ispec=1,NSPEC
     do i=1,NGLLX
        do j=1,NGLLY
           do k=1,NGLLZ
              iglob=ibool(i,j,k,ispec);
              i_target=NINT((xstore(iglob)-xmin)/DX)+1
              j_target=NINT((ystore(iglob)-ymin)/DY)+1
              ! in SEP, z-axis is downward; in SPECFEM, z-axis is upward
              k_target=NINT(-zstore(iglob)/DZ)+1; !NOT using (z_temp-zmin)
              ! Stay in the computational domain
              if(i_target< 1)  i_target=  1
              if(j_target< 1)  j_target=  1
              if(k_target< 1)  k_target=  1
              if(i_target>ni)  i_target= ni
              if(j_target>nj)  j_target= nj
              if(k_target>NZ)  k_target= NZ
              
              var(i,j,k,ispec) = sep_var(i_target,j_target,k_target)
           enddo
        enddo
     enddo
  enddo
end subroutine interpolate_sep_on_mesh

!==============================================================================
!>
subroutine find_slice_bounds_sep(NX, NY, NZ, OX, OY, OZ, DX, DY, DZ, &
                                 xmin, ymin, imin, jmin, kmin, ni, nj, nk)
  use generate_database_par, only: xstor, ystore, zstore
  ! Parameters
  integer, intent(in) :: NX, NY, NZ
  real, intent(in) :: OX, OY, OZ, DX, DY, DZ
  real, intent(out) :: xmin, ymin
  integer, intent(out) :: imin, jmin, kmin, ni, nj, nk

  ! Actual bounds for the current slice mesh
  xmin=minval(xstore); xmax=maxval(xstore);
  ymin=minval(ystore); ymax=maxval(ystore);
  ! in SEP, z-axis is downward; in SPECFEM, z-axis is upward
  zmin=minval(-zstore); zmax=maxval(-zstore);

  ! Bounds for the SEP model
  imin=floor((xmin-OX)/DX+1); imax=ceiling((xmax-OX)/DX+1);
  jmin=floor((ymin-OY)/DY+1); jmax=ceiling((ymax-OY)/DY+1);
  kmin=floor((zmin-OZ)/DZ+1); kmax=ceiling((zmax-OZ)/DZ+1);
  if(imin<1)  imin= 1;  if(jmin<1)  jmin= 1;  if(kmin<1)  kmin= 1;
  if(imax>NX) imin=NX;  if(jmax>NY) jmax=NY;  if(kmax>NZ) kmax=NZ;
  ! Number of SEP indexes for the current slice 
  ni=imax-imin+1
  nj=jmax-jmin+1
  nk=kmax-kmin+1

  ! Corrected bounds for the current slice when looking for SEP indexes
  xmin=DX*(imin-1)+OX; xmax=DX*(imax-1)+OX;
  ymin=DY*(jmin-1)+OY; ymax=DY*(jmax-1)+OY;
  zmin=DZ*(kmin-1)+OZ; zmax=DZ*(kmax-1)+OZ;

end subroutine find_slice_bounds_sep

end module model_sep_mode

