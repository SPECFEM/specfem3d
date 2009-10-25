!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

! initializes mesh arrays for movies and shakemaps
  allocate(nfaces_perproc_surface_ext_mesh(NPROC))
  allocate(faces_surface_offset_ext_mesh(NPROC))
  if (nfaces_surface_external_mesh == 0) then
    if (USE_HIGHRES_FOR_MOVIES) then
      allocate(faces_surface_external_mesh(NGLLX*NGLLY,1))
      allocate(store_val_x_external_mesh(NGLLX*NGLLY*1))
      allocate(store_val_y_external_mesh(NGLLX*NGLLY*1))
      allocate(store_val_z_external_mesh(NGLLX*NGLLY*1))
      allocate(store_val_ux_external_mesh(NGLLX*NGLLY*1))
      allocate(store_val_uy_external_mesh(NGLLX*NGLLY*1))
      allocate(store_val_uz_external_mesh(NGLLX*NGLLY*1))
    else
      allocate(faces_surface_external_mesh(NGNOD2D,1))
      allocate(store_val_x_external_mesh(NGNOD2D*1))
      allocate(store_val_y_external_mesh(NGNOD2D*1))
      allocate(store_val_z_external_mesh(NGNOD2D*1))
      allocate(store_val_ux_external_mesh(NGNOD2D*1))
      allocate(store_val_uy_external_mesh(NGNOD2D*1))
      allocate(store_val_uz_external_mesh(NGNOD2D*1))
    endif
  else
    if (USE_HIGHRES_FOR_MOVIES) then
      allocate(faces_surface_external_mesh(NGLLX*NGLLY,nfaces_surface_external_mesh))
      allocate(store_val_x_external_mesh(NGLLX*NGLLY*nfaces_surface_external_mesh))
      allocate(store_val_y_external_mesh(NGLLX*NGLLY*nfaces_surface_external_mesh))
      allocate(store_val_z_external_mesh(NGLLX*NGLLY*nfaces_surface_external_mesh))
      allocate(store_val_ux_external_mesh(NGLLX*NGLLY*nfaces_surface_external_mesh))
      allocate(store_val_uy_external_mesh(NGLLX*NGLLY*nfaces_surface_external_mesh))
      allocate(store_val_uz_external_mesh(NGLLX*NGLLY*nfaces_surface_external_mesh))
    else
      allocate(faces_surface_external_mesh(NGNOD2D,nfaces_surface_external_mesh))
      allocate(store_val_x_external_mesh(NGNOD2D*nfaces_surface_external_mesh))
      allocate(store_val_y_external_mesh(NGNOD2D*nfaces_surface_external_mesh))
      allocate(store_val_z_external_mesh(NGNOD2D*nfaces_surface_external_mesh))
      allocate(store_val_ux_external_mesh(NGNOD2D*nfaces_surface_external_mesh))
      allocate(store_val_uy_external_mesh(NGNOD2D*nfaces_surface_external_mesh))
      allocate(store_val_uz_external_mesh(NGNOD2D*nfaces_surface_external_mesh))
    endif
  endif

! number of surface faces for all partitions together
  call sum_all_i(nfaces_surface_external_mesh,nfaces_surface_glob_ext_mesh)
  
  if (myrank == 0) then
    if (USE_HIGHRES_FOR_MOVIES) then
      allocate(store_val_x_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh))
      allocate(store_val_y_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh))
      allocate(store_val_z_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh))
      allocate(store_val_ux_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh))
      allocate(store_val_uy_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh))
      allocate(store_val_uz_all_external_mesh(NGLLX*NGLLY*nfaces_surface_glob_ext_mesh))
    else
      allocate(store_val_x_all_external_mesh(NGNOD2D*nfaces_surface_glob_ext_mesh))
      allocate(store_val_y_all_external_mesh(NGNOD2D*nfaces_surface_glob_ext_mesh))
      allocate(store_val_z_all_external_mesh(NGNOD2D*nfaces_surface_glob_ext_mesh))
      allocate(store_val_ux_all_external_mesh(NGNOD2D*nfaces_surface_glob_ext_mesh))
      allocate(store_val_uy_all_external_mesh(NGNOD2D*nfaces_surface_glob_ext_mesh))
      allocate(store_val_uz_all_external_mesh(NGNOD2D*nfaces_surface_glob_ext_mesh))
    endif
  endif
  call gather_all_i(nfaces_surface_external_mesh,1,nfaces_perproc_surface_ext_mesh,1,NPROC)

  faces_surface_offset_ext_mesh(1) = 0
  do i = 2, NPROC
    faces_surface_offset_ext_mesh(i) = sum(nfaces_perproc_surface_ext_mesh(1:i-1))
  enddo
  if (USE_HIGHRES_FOR_MOVIES) then
    faces_surface_offset_ext_mesh(:) = faces_surface_offset_ext_mesh(:)*NGLLX*NGLLY
  else
    faces_surface_offset_ext_mesh(:) = faces_surface_offset_ext_mesh(:)*NGNOD2D
  endif

! stores global indices of GLL points on the surface to array faces_surface_external_mesh
  nfaces_surface_external_mesh = 0
  do ispec = 1, NSPEC_AB
    if (ispec_is_surface_external_mesh(ispec)) then
      iglob = ibool(2,2,1,ispec)
      if (iglob_is_surface_external_mesh(iglob)) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
        if (USE_HIGHRES_FOR_MOVIES) then
          ipoin =0
          do j = NGLLY, 1, -1
            do i = 1, NGLLX
              ipoin = ipoin+1
              faces_surface_external_mesh(ipoin,nfaces_surface_external_mesh) = ibool(i,j,1,ispec)
            enddo
          enddo
        else
          faces_surface_external_mesh(1,nfaces_surface_external_mesh) = ibool(1,1,1,ispec)
          faces_surface_external_mesh(2,nfaces_surface_external_mesh) = ibool(1,NGLLY,1,ispec)
          faces_surface_external_mesh(3,nfaces_surface_external_mesh) = ibool(NGLLX,NGLLY,1,ispec)
          faces_surface_external_mesh(4,nfaces_surface_external_mesh) = ibool(NGLLX,1,1,ispec)
        endif
      endif
      iglob = ibool(2,2,NGLLZ,ispec)
      if (iglob_is_surface_external_mesh(iglob)) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
        if (USE_HIGHRES_FOR_MOVIES) then
          ipoin =0
          do j = 1, NGLLY
            do i = 1, NGLLX
              ipoin = ipoin+1
              faces_surface_external_mesh(ipoin,nfaces_surface_external_mesh) = ibool(i,j,NGLLZ,ispec)
            enddo
          enddo
        else
          faces_surface_external_mesh(1,nfaces_surface_external_mesh) = ibool(1,1,NGLLZ,ispec)
          faces_surface_external_mesh(2,nfaces_surface_external_mesh) = ibool(NGLLX,1,NGLLZ,ispec)
          faces_surface_external_mesh(3,nfaces_surface_external_mesh) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
          faces_surface_external_mesh(4,nfaces_surface_external_mesh) = ibool(1,NGLLY,NGLLZ,ispec)
        endif
      endif
      iglob = ibool(2,1,2,ispec)
      if (iglob_is_surface_external_mesh(iglob)) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
        if (USE_HIGHRES_FOR_MOVIES) then
          ipoin =0
          do k = 1, NGLLZ
            do i = 1, NGLLX
              ipoin = ipoin+1
              faces_surface_external_mesh(ipoin,nfaces_surface_external_mesh) = ibool(i,1,k,ispec)
            enddo
          enddo
        else
          faces_surface_external_mesh(1,nfaces_surface_external_mesh) = ibool(1,1,1,ispec)
          faces_surface_external_mesh(2,nfaces_surface_external_mesh) = ibool(NGLLX,1,1,ispec)
          faces_surface_external_mesh(3,nfaces_surface_external_mesh) = ibool(NGLLX,1,NGLLZ,ispec)
          faces_surface_external_mesh(4,nfaces_surface_external_mesh) = ibool(1,1,NGLLZ,ispec)
        endif
      endif
      iglob = ibool(2,NGLLY,2,ispec)
      if (iglob_is_surface_external_mesh(iglob)) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
        if (USE_HIGHRES_FOR_MOVIES) then
          ipoin =0
          do k = 1, NGLLZ
            do i = NGLLX, 1, -1
              ipoin = ipoin+1
              faces_surface_external_mesh(ipoin,nfaces_surface_external_mesh) = ibool(i,NGLLY,k,ispec)
            enddo
          enddo
        else
          faces_surface_external_mesh(1,nfaces_surface_external_mesh) = ibool(NGLLX,NGLLY,1,ispec)
          faces_surface_external_mesh(2,nfaces_surface_external_mesh) = ibool(1,NGLLY,1,ispec)
          faces_surface_external_mesh(3,nfaces_surface_external_mesh) = ibool(1,NGLLY,NGLLZ,ispec)
          faces_surface_external_mesh(4,nfaces_surface_external_mesh) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
        endif
      endif
      iglob = ibool(1,2,2,ispec)
      if (iglob_is_surface_external_mesh(iglob)) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
        if (USE_HIGHRES_FOR_MOVIES) then
          ipoin =0
          do k = 1, NGLLZ
            do j = NGLLY, 1, -1
              ipoin = ipoin+1
              faces_surface_external_mesh(ipoin,nfaces_surface_external_mesh) = ibool(1,j,k,ispec)
            enddo
         enddo
        else
          faces_surface_external_mesh(1,nfaces_surface_external_mesh) = ibool(1,NGLLY,1,ispec)
          faces_surface_external_mesh(2,nfaces_surface_external_mesh) = ibool(1,1,1,ispec)
          faces_surface_external_mesh(3,nfaces_surface_external_mesh) = ibool(1,1,NGLLZ,ispec)
          faces_surface_external_mesh(4,nfaces_surface_external_mesh) = ibool(1,NGLLY,NGLLZ,ispec)
        endif
      endif
      iglob = ibool(NGLLX,2,2,ispec)
      if (iglob_is_surface_external_mesh(iglob)) then
        nfaces_surface_external_mesh = nfaces_surface_external_mesh + 1
        if (USE_HIGHRES_FOR_MOVIES) then
          ipoin =0
          do k = 1, NGLLZ
            do j = 1, NGLLY
              ipoin = ipoin+1
              faces_surface_external_mesh(ipoin,nfaces_surface_external_mesh) = ibool(NGLLX,j,k,ispec)
            enddo
         enddo
        else
          faces_surface_external_mesh(1,nfaces_surface_external_mesh) = ibool(NGLLX,1,1,ispec)
          faces_surface_external_mesh(2,nfaces_surface_external_mesh) = ibool(NGLLX,NGLLY,1,ispec)
          faces_surface_external_mesh(3,nfaces_surface_external_mesh) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
          faces_surface_external_mesh(4,nfaces_surface_external_mesh) = ibool(NGLLX,1,NGLLZ,ispec)
        endif
      endif

    endif
  enddo ! NSPEC_AB

  if (myrank == 0) then 
    write(IMAIN,*) 'movie: nfaces_surface_external_mesh   = ',nfaces_surface_external_mesh
    write(IMAIN,*) 'movie: nfaces_perproc_surface_ext_mesh = ',nfaces_perproc_surface_ext_mesh
    write(IMAIN,*) 'movie: nfaces_surface_glob_ext_mesh    = ',nfaces_surface_glob_ext_mesh
  endif

  
  end subroutine
  
  
  
  
  
