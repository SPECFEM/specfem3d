!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
!
! United States and French Government Sponsorship Acknowledged.

! creation of arrays for movie and shakemap routines for external meshes

  subroutine setup_movie_meshes()

  use specfem_par
  use specfem_par_movie
  implicit none

  integer :: i,j,k,ispec,iglob,ier
  integer :: ipoin,nfaces_org
  character(len=256):: filename

! initializes mesh arrays for movies and shakemaps
  allocate(nfaces_perproc_surface_ext_mesh(NPROC), &
          faces_surface_offset_ext_mesh(NPROC),stat=ier)
  if( ier /= 0 ) stop 'error allocating array for movie faces'

  nfaces_org = nfaces_surface_ext_mesh
  if (nfaces_surface_ext_mesh == 0) then
    ! dummy arrays
    if (USE_HIGHRES_FOR_MOVIES) then
      allocate(faces_surface_ext_mesh(NGLLX*NGLLY,1), &
        store_val_x_external_mesh(NGLLX*NGLLY*1), &
        store_val_y_external_mesh(NGLLX*NGLLY*1), &
        store_val_z_external_mesh(NGLLX*NGLLY*1), &
        store_val_ux_external_mesh(NGLLX*NGLLY*1), &
        store_val_uy_external_mesh(NGLLX*NGLLY*1), &
        store_val_uz_external_mesh(NGLLX*NGLLY*1),stat=ier)
      if( ier /= 0 ) stop 'error allocating dummy arrays for highres movie'
    else
      allocate(faces_surface_ext_mesh(NGNOD2D_FOUR_CORNERS,1), &
        store_val_x_external_mesh(NGNOD2D_FOUR_CORNERS*1), &
        store_val_y_external_mesh(NGNOD2D_FOUR_CORNERS*1), &
        store_val_z_external_mesh(NGNOD2D_FOUR_CORNERS*1), &
        store_val_ux_external_mesh(NGNOD2D_FOUR_CORNERS*1), &
        store_val_uy_external_mesh(NGNOD2D_FOUR_CORNERS*1), &
        store_val_uz_external_mesh(NGNOD2D_FOUR_CORNERS*1),stat=ier)
      if( ier /= 0 ) stop 'error allocating dummy arrays for lowres movie'
    endif
  else
    if (USE_HIGHRES_FOR_MOVIES) then
      allocate(faces_surface_ext_mesh(NGLLX*NGLLY,nfaces_surface_ext_mesh), &
        store_val_x_external_mesh(NGLLX*NGLLY*nfaces_surface_ext_mesh), &
        store_val_y_external_mesh(NGLLX*NGLLY*nfaces_surface_ext_mesh), &
        store_val_z_external_mesh(NGLLX*NGLLY*nfaces_surface_ext_mesh), &
        store_val_ux_external_mesh(NGLLX*NGLLY*nfaces_surface_ext_mesh), &
        store_val_uy_external_mesh(NGLLX*NGLLY*nfaces_surface_ext_mesh), &
        store_val_uz_external_mesh(NGLLX*NGLLY*nfaces_surface_ext_mesh),stat=ier)
      if( ier /= 0 ) stop 'error allocating arrays for highres movie'
    else
      allocate(faces_surface_ext_mesh(NGNOD2D_FOUR_CORNERS,nfaces_surface_ext_mesh), &
        store_val_x_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_ext_mesh), &
        store_val_y_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_ext_mesh), &
        store_val_z_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_ext_mesh), &
        store_val_ux_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_ext_mesh), &
        store_val_uy_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_ext_mesh), &
        store_val_uz_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_ext_mesh),stat=ier)
      if( ier /= 0 ) stop 'error allocating arrays for lowres movie'
    endif
  endif
  store_val_ux_external_mesh(:) = 0._CUSTOM_REAL
  store_val_uy_external_mesh(:) = 0._CUSTOM_REAL
  store_val_uz_external_mesh(:) = 0._CUSTOM_REAL

  ! number of surface faces for all partitions together
  call sum_all_i(nfaces_surface_ext_mesh,nfaces_surface_glob_ext_mesh)

  ! arrays used for collected/gathered fields
  if (myrank == 0) then
    if (USE_HIGHRES_FOR_MOVIES) then
      allocate(store_val_x_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh), &
        store_val_y_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh), &
        store_val_z_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh), &
        store_val_ux_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh), &
        store_val_uy_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh), &
        store_val_uz_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh),stat=ier)
      if( ier /= 0 ) stop 'error allocating arrays for highres movie'
    else
      allocate(store_val_x_all_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_glob_ext_mesh), &
        store_val_y_all_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_glob_ext_mesh), &
        store_val_z_all_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_glob_ext_mesh), &
        store_val_ux_all_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_glob_ext_mesh), &
        store_val_uy_all_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_glob_ext_mesh), &
        store_val_uz_all_external_mesh(NGNOD2D_FOUR_CORNERS*nfaces_surface_glob_ext_mesh),stat=ier)
      if( ier /= 0 ) stop 'error allocating arrays for lowres movie'
    endif
  endif
  call gather_all_i(nfaces_surface_ext_mesh,1,nfaces_perproc_surface_ext_mesh,1,NPROC)

  ! array offsets
  faces_surface_offset_ext_mesh(1) = 0
  do i = 2, NPROC
    faces_surface_offset_ext_mesh(i) = sum(nfaces_perproc_surface_ext_mesh(1:i-1))
  enddo
  if (USE_HIGHRES_FOR_MOVIES) then
    faces_surface_offset_ext_mesh(:) = faces_surface_offset_ext_mesh(:)*NGLLX*NGLLY
  else
    faces_surface_offset_ext_mesh(:) = faces_surface_offset_ext_mesh(:)*NGNOD2D_FOUR_CORNERS
  endif

! stores global indices of GLL points on the surface to array faces_surface_ext_mesh
  if( MOVIE_TYPE == 2 ) then

    allocate( faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array faces_surface_ext_mesh_ispec'

    ! stores global indices
    nfaces_surface_ext_mesh = 0
    do ispec = 1, NSPEC_AB

      if (ispec_is_surface_external_mesh(ispec)) then

        ! zmin face
        iglob = ibool(2,2,1,ispec)
        if (iglob_is_surface_external_mesh(iglob)) then
          nfaces_surface_ext_mesh = nfaces_surface_ext_mesh + 1
          faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh) = ispec
          if (USE_HIGHRES_FOR_MOVIES) then
            ipoin =0
            do j = NGLLY, 1, -1
              do i = 1, NGLLX
                ipoin = ipoin+1
                faces_surface_ext_mesh(ipoin,nfaces_surface_ext_mesh) = ibool(i,j,1,ispec)
              enddo
            enddo
          else
            faces_surface_ext_mesh(1,nfaces_surface_ext_mesh) = ibool(1,1,1,ispec)
            faces_surface_ext_mesh(2,nfaces_surface_ext_mesh) = ibool(1,NGLLY,1,ispec)
            faces_surface_ext_mesh(3,nfaces_surface_ext_mesh) = ibool(NGLLX,NGLLY,1,ispec)
            faces_surface_ext_mesh(4,nfaces_surface_ext_mesh) = ibool(NGLLX,1,1,ispec)
          endif
        endif
        ! zmax face
        iglob = ibool(2,2,NGLLZ,ispec)
        if (iglob_is_surface_external_mesh(iglob)) then
          nfaces_surface_ext_mesh = nfaces_surface_ext_mesh + 1
          faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh) = ispec
          if (USE_HIGHRES_FOR_MOVIES) then
            ipoin =0
            do j = 1, NGLLY
              do i = 1, NGLLX
                ipoin = ipoin+1
                faces_surface_ext_mesh(ipoin,nfaces_surface_ext_mesh) = ibool(i,j,NGLLZ,ispec)
              enddo
            enddo
          else
            faces_surface_ext_mesh(1,nfaces_surface_ext_mesh) = ibool(1,1,NGLLZ,ispec)
            faces_surface_ext_mesh(2,nfaces_surface_ext_mesh) = ibool(NGLLX,1,NGLLZ,ispec)
            faces_surface_ext_mesh(3,nfaces_surface_ext_mesh) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
            faces_surface_ext_mesh(4,nfaces_surface_ext_mesh) = ibool(1,NGLLY,NGLLZ,ispec)
          endif
        endif
        ! ymin face
        iglob = ibool(2,1,2,ispec)
        if (iglob_is_surface_external_mesh(iglob)) then
          nfaces_surface_ext_mesh = nfaces_surface_ext_mesh + 1
          faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh) = ispec
          if (USE_HIGHRES_FOR_MOVIES) then
            ipoin =0
            do k = 1, NGLLZ
              do i = 1, NGLLX
                ipoin = ipoin+1
                faces_surface_ext_mesh(ipoin,nfaces_surface_ext_mesh) = ibool(i,1,k,ispec)
              enddo
            enddo
          else
            faces_surface_ext_mesh(1,nfaces_surface_ext_mesh) = ibool(1,1,1,ispec)
            faces_surface_ext_mesh(2,nfaces_surface_ext_mesh) = ibool(NGLLX,1,1,ispec)
            faces_surface_ext_mesh(3,nfaces_surface_ext_mesh) = ibool(NGLLX,1,NGLLZ,ispec)
            faces_surface_ext_mesh(4,nfaces_surface_ext_mesh) = ibool(1,1,NGLLZ,ispec)
          endif
        endif
        ! ymax face
        iglob = ibool(2,NGLLY,2,ispec)
        if (iglob_is_surface_external_mesh(iglob)) then
          nfaces_surface_ext_mesh = nfaces_surface_ext_mesh + 1
          faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh) = ispec
          if (USE_HIGHRES_FOR_MOVIES) then
            ipoin =0
            do k = 1, NGLLZ
              do i = NGLLX, 1, -1
                ipoin = ipoin+1
                faces_surface_ext_mesh(ipoin,nfaces_surface_ext_mesh) = ibool(i,NGLLY,k,ispec)
              enddo
            enddo
          else
            faces_surface_ext_mesh(1,nfaces_surface_ext_mesh) = ibool(NGLLX,NGLLY,1,ispec)
            faces_surface_ext_mesh(2,nfaces_surface_ext_mesh) = ibool(1,NGLLY,1,ispec)
            faces_surface_ext_mesh(3,nfaces_surface_ext_mesh) = ibool(1,NGLLY,NGLLZ,ispec)
            faces_surface_ext_mesh(4,nfaces_surface_ext_mesh) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
          endif
        endif
        ! xmin face
        iglob = ibool(1,2,2,ispec)
        if (iglob_is_surface_external_mesh(iglob)) then
          nfaces_surface_ext_mesh = nfaces_surface_ext_mesh + 1
          faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh) = ispec
          if (USE_HIGHRES_FOR_MOVIES) then
            ipoin =0
            do k = 1, NGLLZ
              do j = NGLLY, 1, -1
                ipoin = ipoin+1
                faces_surface_ext_mesh(ipoin,nfaces_surface_ext_mesh) = ibool(1,j,k,ispec)
              enddo
           enddo
          else
            faces_surface_ext_mesh(1,nfaces_surface_ext_mesh) = ibool(1,NGLLY,1,ispec)
            faces_surface_ext_mesh(2,nfaces_surface_ext_mesh) = ibool(1,1,1,ispec)
            faces_surface_ext_mesh(3,nfaces_surface_ext_mesh) = ibool(1,1,NGLLZ,ispec)
            faces_surface_ext_mesh(4,nfaces_surface_ext_mesh) = ibool(1,NGLLY,NGLLZ,ispec)
          endif
        endif
        ! xmax face
        iglob = ibool(NGLLX,2,2,ispec)
        if (iglob_is_surface_external_mesh(iglob)) then
          nfaces_surface_ext_mesh = nfaces_surface_ext_mesh + 1
          faces_surface_ext_mesh_ispec(nfaces_surface_ext_mesh) = ispec
          if (USE_HIGHRES_FOR_MOVIES) then
            ipoin =0
            do k = 1, NGLLZ
              do j = 1, NGLLY
                ipoin = ipoin+1
                faces_surface_ext_mesh(ipoin,nfaces_surface_ext_mesh) = ibool(NGLLX,j,k,ispec)
              enddo
           enddo
          else
            faces_surface_ext_mesh(1,nfaces_surface_ext_mesh) = ibool(NGLLX,1,1,ispec)
            faces_surface_ext_mesh(2,nfaces_surface_ext_mesh) = ibool(NGLLX,NGLLY,1,ispec)
            faces_surface_ext_mesh(3,nfaces_surface_ext_mesh) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
            faces_surface_ext_mesh(4,nfaces_surface_ext_mesh) = ibool(NGLLX,1,NGLLZ,ispec)
          endif
        endif
      endif
    enddo ! NSPEC_AB

    ! checks number of faces
    if( nfaces_surface_ext_mesh /= nfaces_org ) then
      print*,'error number of movie faces: ',nfaces_surface_ext_mesh,nfaces_org
      call exit_mpi(myrank,'error number of faces')
    endif
  endif ! MOVIE_TYPE == 2

  ! user output
  if (myrank == 0) then
    if( PLOT_CROSS_SECTIONS ) then
      write(IMAIN,*) 'movie cross-sections:'
    else
      write(IMAIN,*) 'movie surface:'
    endif
    write(IMAIN,*) '  nfaces_surface_ext_mesh:',nfaces_surface_ext_mesh
    write(IMAIN,*) '  nfaces_perproc_surface_ext_mesh:',nfaces_perproc_surface_ext_mesh
    write(IMAIN,*) '  nfaces_surface_glob_ext_mesh:',nfaces_surface_glob_ext_mesh

    ! updates number of surface elements in an include file for the movies
    if( nfaces_surface_glob_ext_mesh > 0 ) then
      filename = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)) // '/surface_from_mesher.h'
      open(unit=IOUT,file=trim(filename),status='unknown')
      write(IOUT,*) '!'
      write(IOUT,*) '! this is the parameter file for static compilation for movie creation'
      write(IOUT,*) '!'
      write(IOUT,*) '! number of elements containing surface faces '
      write(IOUT,*) '! ---------------'
      write(IOUT,*)
      write(IOUT,*) 'integer,parameter :: NSPEC_SURFACE_EXT_MESH = ',nfaces_surface_glob_ext_mesh
      write(IOUT,*)
      close(IOUT)
    endif

  endif

  ! for gathering movie data
  if (USE_HIGHRES_FOR_MOVIES) then
    ! hi-res movies output all gll points on surface
    nfaces_perproc_surface_ext_mesh(:) = nfaces_perproc_surface_ext_mesh(:)*NGLLX*NGLLY
    nfaces_surface_ext_mesh_points = nfaces_surface_ext_mesh*NGLLX*NGLLY
    nfaces_surface_glob_em_points = nfaces_surface_glob_ext_mesh*NGLLX*NGLLY
  else
    ! low-res movies only output at element corners
    nfaces_perproc_surface_ext_mesh(:) = nfaces_perproc_surface_ext_mesh(:)*NGNOD2D_FOUR_CORNERS
    nfaces_surface_ext_mesh_points = nfaces_surface_ext_mesh*NGNOD2D_FOUR_CORNERS
    nfaces_surface_glob_em_points = nfaces_surface_glob_ext_mesh*NGNOD2D_FOUR_CORNERS
  endif

  end subroutine setup_movie_meshes

