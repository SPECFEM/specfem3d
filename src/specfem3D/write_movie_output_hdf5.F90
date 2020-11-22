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


#include "config.fh"

  subroutine write_movie_output_h5()

  use specfem_par
  use specfem_par_movie
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  if (.not. MOVIE_SIMULATION) return

  ! gets resulting array values onto CPU
  if (GPU_MODE .and. &
    ( &
      CREATE_SHAKEMAP .or. &
      ( MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) .or. &
      ( MOVIE_VOLUME .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) .or. &
      ( PNM_IMAGE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) &
     )) then
    ! acoustic domains
    if (ACOUSTIC_SIMULATION) then
      ! transfers whole fields
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                potential_dot_acoustic,potential_dot_dot_acoustic,Mesh_pointer)
    endif
    ! elastic domains
    if (ELASTIC_SIMULATION) then
      ! transfers whole fields
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc, accel, Mesh_pointer)
    endif
  endif

  ! computes SHAKING INTENSITY MAP
  if (CREATE_SHAKEMAP) then
    call wmo_create_shakemap_h5()
  endif

  if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    ! wait
    if(NIONOD > 0) call wait_all_send()

    ! saves MOVIE on the SURFACE
    if (MOVIE_SURFACE) then
      call wmo_movie_surface_output_h5()
    endif

    ! saves MOVIE in full 3D MESH
    if (MOVIE_VOLUME) then
      call wmo_movie_volume_output_h5()
    endif

    ! creates cross-section PNM image
    if (PNM_IMAGE) then
      call write_PNM_create_image() ! this is excluded from this file. need to check
    endif

  endif

  end subroutine write_movie_output_h5


  subroutine wait_all_send()
  use specfem_par

  implicit none

  integer :: ireq

  ! wait till all mpi_isends are finished
  if (n_req_surf /= 0) then
    do ireq=1, n_req_surf
      call wait_req(req_dump_surf(ireq))
    enddo
  endif
  ! wait till all mpi_isends are finished
  if (n_req_vol /= 0) then
    do ireq=1, n_req_vol
      call wait_req(req_dump_vol(ireq))
    enddo
  endif

  n_req_surf = 0; n_req_vol = 0

  call synchronize_all()

  end subroutine wait_all_send


!================================================================

  subroutine wmo_movie_surface_output_h5()

! output of surface moviedata files

! option MOVIE_TYPE == 1: only at top, free surface
!        MOVIE_TYPE == 2: for all external, outer mesh surfaces

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie
  use phdf5_utils
  implicit none

  ! temporary array for single elements
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: val_element
  real(kind=CUSTOM_REAL),dimension(1):: dummy
  integer :: ispec2D,ispec,ipoin,iglob,ier,ia,ireq
  integer :: npoin_elem

  ! surface points for single face
  if (USE_HIGHRES_FOR_MOVIES) then
    npoin_elem = NGLLX*NGLLY
  else
    npoin_elem = NGNOD2D_FOUR_CORNERS
  endif

  ! saves surface velocities
  do ispec2D = 1,nfaces_surface
    ispec = faces_surface_ispec(ispec2D)

    if (ispec_is_acoustic(ispec)) then
      ! acoustic elements
      if (SAVE_DISPLACEMENT) then
        ! displacement vector
        call compute_gradient_in_acoustic(ispec,potential_acoustic,val_element)
      else
        ! velocity vector
        call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,val_element)
      endif

      ! all surface element points
      do ipoin = 1, npoin_elem
        ia = npoin_elem * (ispec2D - 1) + ipoin
        iglob = faces_surface_ibool(ipoin,ispec2D)
        ! puts displ/velocity values into storage array
        call wmo_get_vel_vector(ispec,iglob,ia,val_element)
      enddo

    else if (ispec_is_elastic(ispec)) then
      ! elastic elements
      ! all surface element points
      do ipoin = 1, npoin_elem
        ia = npoin_elem * (ispec2D - 1) + ipoin
        iglob = faces_surface_ibool(ipoin,ispec2D)

        if (SAVE_DISPLACEMENT) then
          ! velocity x,y,z-components
          store_val_ux(ia) = displ(1,iglob)
          store_val_uy(ia) = displ(2,iglob)
          store_val_uz(ia) = displ(3,iglob)
        else
          ! velocity x,y,z-components
          store_val_ux(ia) = veloc(1,iglob)
          store_val_uy(ia) = veloc(2,iglob)
          store_val_uz(ia) = veloc(3,iglob)
        endif
      enddo
    endif
  enddo

  if (NIONOD > 0) then
    ! send surface body to io node
    call isend_cr_inter(store_val_ux,nfaces_surface_points,0,io_tag_surface_ux,req_dump_surf(1))
    call isend_cr_inter(store_val_uy,nfaces_surface_points,0,io_tag_surface_uy,req_dump_surf(2))
    call isend_cr_inter(store_val_uz,nfaces_surface_points,0,io_tag_surface_uz,req_dump_surf(3))

    n_req_surf = 3
  else
    call write_movie_surface_noserv()
  endif

end subroutine wmo_movie_surface_output_h5


subroutine write_movie_surface_noserv()
  use io_server
  use phdf5_utils
  use specfem_par
  use specfem_par_movie
  use shared_parameters

  implicit none

  integer                                           :: ier, nfaces_actual, nfaces_aug=(NGLLX-1)*(NGLLY-1),nnodes_per_face_aug=4
  integer                                           :: len_array_aug, len_array_aug_proc
  character(len=64)                                 :: dset_name, group_name, tempstr
  real(kind=CUSTOM_REAL)                            :: aug_factor
  logical :: if_corrective = .true.

  ! initialize h5 and xdmf  if it == 0
  type(h5io) :: h5
  h5 = h5io()
  fname_h5_data_surf = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_surface.h5"


  ! get the offset info from main rank
  if (it == NTSTEP_BETWEEN_FRAMES) then
    call bcast_all_cr(faces_surface_offset,size(faces_surface_offset))
  endif

  if(myrank == 0) then
    size_surf_array = size(store_val_x_all)
    call bcast_all_singlei(size_surf_array)
  else
    call bcast_all_singlei(size_surf_array)
  endif

  nfaces_actual = size_surf_array/(NGLLX*NGLLY)
  len_array_aug = nfaces_actual*nfaces_aug*nnodes_per_face_aug

  ! initialization of h5 file
  call h5_init(h5, fname_h5_data_surf)

  if (it == NTSTEP_BETWEEN_FRAMES .and. myrank == 0) then
    ! create a hdf5 file
    call h5_create_file(h5)

    ! save xyz coordinate array
    group_name = "surf_coord"
    call h5_create_group(h5, group_name)
    call h5_open_group(h5, group_name)

    ! low resolution output
    if (.not. USE_HIGHRES_FOR_MOVIES) then
      dset_name = "x"
      call h5_write_dataset_1d_d(h5, dset_name, store_val_x_all)
      call h5_close_dataset(h5)
      dset_name = "y"
      call h5_write_dataset_1d_d(h5, dset_name, store_val_y_all)
      call h5_close_dataset(h5)
      dset_name = "z"
      call h5_write_dataset_1d_d(h5, dset_name, store_val_z_all)
      call h5_close_dataset(h5)

      ! write xdmf header
      call write_xdmf_surface_header(size_surf_array,pos_xdmf_surf)

      ! high resolution output
    else
      allocate(surf_x_aug(len_array_aug),stat=ier)
      allocate(surf_y_aug(len_array_aug),stat=ier)
      allocate(surf_z_aug(len_array_aug),stat=ier)

      ! nfaces*25nodes => n*16faces*4
      dset_name = "x"
      call recompose_for_hires(store_val_x_all,surf_x_aug)
      call h5_write_dataset_1d_d(h5, dset_name, surf_x_aug)
      call h5_close_dataset(h5)
      dset_name = "y"
      call recompose_for_hires(store_val_y_all,surf_y_aug)
      call h5_write_dataset_1d_d(h5, dset_name, surf_y_aug)
      call h5_close_dataset(h5)
      dset_name = "z"
      call recompose_for_hires(store_val_z_all,surf_z_aug)
      call h5_write_dataset_1d_d(h5, dset_name, surf_z_aug)
      call h5_close_dataset(h5)

      ! write xdmf header
      call write_xdmf_surface_header(len_array_aug,pos_xdmf_surf)

      deallocate(surf_x_aug, surf_y_aug, surf_z_aug, stat= ier)
    endif

    call h5_close_group(h5)
    call h5_close_file(h5)

  endif

  call synchronize_all()


  ! write dataset in one single h5 file in collective mode.
  ! create a group for each io step
  write(tempstr, "(i6.6)") it
  group_name = "it_"//tempstr

  if (myrank == 0) then
    call h5_open_file(h5)
    call h5_create_group(h5, group_name)
    call h5_open_group(h5, group_name)
    if(.not. USE_HIGHRES_FOR_MOVIES)then
      call h5_create_dataset_gen_in_group(h5,"ux",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"uy",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"uz",(/size_surf_array/),1,CUSTOM_REAL)
    else
      call h5_create_dataset_gen_in_group(h5,"ux",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"uy",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"uz",(/len_array_aug/),1,CUSTOM_REAL)
   endif
    call h5_close_group(h5)
    call h5_close_file(h5)
  endif

  call synchronize_all()
  call h5_open_file_p(h5)
  call h5_open_group(h5, group_name)

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "ux"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, store_val_ux, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "uy"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, store_val_uy, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "uz"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, store_val_uz, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)

    ! write xdmf body
    if (myrank==0) call write_xdmf_surface_body(it, size_surf_array,pos_xdmf_surf)

  else
    nfaces_actual      = size(store_val_ux)/(NGLLX*NGLLY)
    len_array_aug_proc = nfaces_actual*nfaces_aug*nnodes_per_face_aug ! augmented array length for each proc
    aug_factor         = real(len_array_aug)/real(size_surf_array)

    allocate(surf_ux_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uy_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uz_aug(len_array_aug_proc),stat=ier)

    dset_name = "ux"
    call recompose_for_hires(store_val_ux, surf_ux_aug)
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            surf_ux_aug, (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "uy"
    call recompose_for_hires(store_val_uy, surf_uy_aug)
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            surf_uy_aug, (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "uz"
    call recompose_for_hires(store_val_uz, surf_uz_aug)
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            surf_uz_aug, (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)

    ! write xdmf body
    if(myrank==0) call write_xdmf_surface_body(it, len_array_aug, pos_xdmf_surf)

    deallocate(surf_ux_aug, surf_uy_aug, surf_uz_aug, stat= ier)
  endif

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine write_movie_surface_noserv

subroutine wmo_create_shakemap_h5()

! creation of shakemap file
!
! option MOVIE_TYPE == 1: uses horizontal peak-ground values, only at top, free surface
!        MOVIE_TYPE == 2: uses norm of vector as peak-ground, for all external, outer mesh surfaces

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie
  implicit none

  ! temporary array for single elements
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: displ_element,veloc_element,accel_element
  integer :: ipoin,ispec,iglob,ispec2D,ia
  integer :: npoin_elem

  ! note: shakemap arrays are initialized to zero after allocation

  ! surface points for single face
  if (USE_HIGHRES_FOR_MOVIES) then
    npoin_elem = NGLLX*NGLLY
  else
    npoin_elem = NGNOD2D_FOUR_CORNERS
  endif

  ! determines displacement, velocity and acceleration maximum amplitudes
  do ispec2D = 1,nfaces_surface
    ispec = faces_surface_ispec(ispec2D)

    if (ispec_is_acoustic(ispec)) then
      ! acoustic elements

      ! computes displ/veloc/accel for local element
      ! displacement vector
      call compute_gradient_in_acoustic(ispec,potential_acoustic,displ_element)
      ! velocity vector
      call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,veloc_element)
      ! accel ?
      call compute_gradient_in_acoustic(ispec,potential_dot_dot_acoustic,accel_element)

      ! all surface element points
      do ipoin = 1, npoin_elem
        ia = npoin_elem * (ispec2D - 1) + ipoin
        iglob = faces_surface_ibool(ipoin,ispec2D)

        if (MOVIE_TYPE == 1) then
          ! only top surface
          ! horizontal peak-ground value
          call wmo_get_max_vector_top(ispec,iglob,ia,displ_element,veloc_element,accel_element)
        else
          ! all outer surfaces
          ! norm of particle displ/veloc/accel vector
          call wmo_get_max_vector_norm(ispec,iglob,ia,displ_element,veloc_element,accel_element)
        endif
      enddo

    else if (ispec_is_elastic(ispec)) then
      ! elastic elements
      ! all surface element points
      do ipoin = 1, npoin_elem
        ia = npoin_elem * (ispec2D - 1) + ipoin
        iglob = faces_surface_ibool(ipoin,ispec2D)

        if (MOVIE_TYPE == 1) then
          ! only top surface, using horizontal peak-ground value
          ! horizontal displacement
          shakemap_ux(ia) = max(shakemap_ux(ia),abs(displ(1,iglob)),abs(displ(2,iglob)))
          ! horizontal velocity
          shakemap_uy(ia) = max(shakemap_uy(ia),abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          ! horizontal acceleration
          shakemap_uz(ia) = max(shakemap_uz(ia),abs(accel(1,iglob)),abs(accel(2,iglob)))
        else
          ! all outer surfaces, using norm of particle displ/veloc/accel vector
          ! saves norm of displacement,velocity and acceleration vector
          ! norm of displacement
          shakemap_ux(ia) = max(shakemap_ux(ia),sqrt(displ(1,iglob)**2 + displ(2,iglob)**2 + displ(3,iglob)**2))
          ! norm of velocity
          shakemap_uy(ia) = max(shakemap_uy(ia),sqrt(veloc(1,iglob)**2 + veloc(2,iglob)**2 + veloc(3,iglob)**2))
          ! norm of acceleration
          shakemap_uz(ia) = max(shakemap_uz(ia),sqrt(accel(1,iglob)**2 + accel(2,iglob)**2 + accel(3,iglob)**2))
        endif
      enddo
    else
      ! other element types not supported yet
      call exit_MPI(myrank,'Invalid element for shakemap, only acoustic or elastic elements are supported for now')
    endif

  enddo

  ! saves shakemap only at the end of the simulation
  if (it == NSTEP) call wmo_save_shakemap_h5()

  end subroutine wmo_create_shakemap_h5


!================================================================

  subroutine wmo_save_shakemap_h5()

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: req

  if (NIONOD > 0) then
    ! send surface body to io node
    call isend_cr_inter(shakemap_ux,nfaces_surface_points,0,io_tag_shake_ux,req)
    call isend_cr_inter(shakemap_uy,nfaces_surface_points,0,io_tag_shake_uy,req)
    call isend_cr_inter(shakemap_uz,nfaces_surface_points,0,io_tag_shake_uz,req)
  else
    call write_shakemap_noserv()
  endif
end subroutine wmo_save_shakemap_h5


subroutine write_shakemap_noserv()
  use io_server
  use phdf5_utils
  use specfem_par
  use specfem_par_movie

  implicit none

  integer                                           :: ier, nfaces_actual, nfaces_aug=(NGLLX-1)*(NGLLY-1),nnodes_per_face_aug=4
  integer                                           :: len_array_aug, len_array_aug_proc
  character(len=64)                                 :: dset_name, group_name, tempstr
  real(kind=CUSTOM_REAL)                            :: aug_factor
  logical :: if_corrective = .true.

  ! initialize h5 and xdmf  if it == 0
  type(h5io) :: h5
  h5 = h5io()
  fname_h5_data_surf = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/shakemap.h5"

  ! get the offset info from main rank
  call bcast_all_cr(faces_surface_offset,size(faces_surface_offset))

  if(myrank == 0) then
    size_surf_array = size(store_val_x_all)
    call bcast_all_singlei(size_surf_array)
  else
    call bcast_all_singlei(size_surf_array)
  endif

  nfaces_actual = size_surf_array/(NGLLX*NGLLY)
  len_array_aug = nfaces_actual*nfaces_aug*nnodes_per_face_aug

  ! initialization of h5 file
  call h5_init(h5, fname_h5_data_surf)

  if (myrank == 0) then
    ! create a hdf5 file
    call h5_create_file(h5)

    ! save xyz coordinate array
    group_name = "surf_coord"
    call h5_create_group(h5, group_name)
    call h5_open_group(h5, group_name)

    ! low resolution output
    if (.not. USE_HIGHRES_FOR_MOVIES) then
      dset_name = "x"
      call h5_write_dataset_1d_d(h5, dset_name, store_val_x_all)
      call h5_close_dataset(h5)
      dset_name = "y"
      call h5_write_dataset_1d_d(h5, dset_name, store_val_y_all)
      call h5_close_dataset(h5)
      dset_name = "z"
      call h5_write_dataset_1d_d(h5, dset_name, store_val_z_all)
      call h5_close_dataset(h5)

      ! write xdmf header
      call write_xdmf_shakemap(size_surf_array)

      ! high resolution output
    else
      allocate(surf_x_aug(len_array_aug),stat=ier)
      allocate(surf_y_aug(len_array_aug),stat=ier)
      allocate(surf_z_aug(len_array_aug),stat=ier)

      ! nfaces*25nodes => n*16faces*4
      dset_name = "x"
      call recompose_for_hires(store_val_x_all,surf_x_aug)
      call h5_write_dataset_1d_d(h5, dset_name, surf_x_aug)
      call h5_close_dataset(h5)
      dset_name = "y"
      call recompose_for_hires(store_val_y_all,surf_y_aug)
      call h5_write_dataset_1d_d(h5, dset_name, surf_y_aug)
      call h5_close_dataset(h5)
      dset_name = "z"
      call recompose_for_hires(store_val_z_all,surf_z_aug)
      call h5_write_dataset_1d_d(h5, dset_name, surf_z_aug)
      call h5_close_dataset(h5)

      ! write xdmf header
      call write_xdmf_shakemap(len_array_aug)

      deallocate(surf_x_aug, surf_y_aug, surf_z_aug, stat= ier)
    endif

    call h5_close_group(h5)
    call h5_close_file(h5)

  endif

  call synchronize_all()


  ! write dataset in one single h5 file in collective mode.
  ! create a group for each io step
  write(tempstr, "(i6.6)") it
  group_name = "shakemap"

  if (myrank == 0) then
    call h5_open_file(h5)
    call h5_create_group(h5, group_name)
    call h5_open_group(h5, group_name)
    if(.not. USE_HIGHRES_FOR_MOVIES)then
      call h5_create_dataset_gen_in_group(h5,"shakemap_ux",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"shakemap_uy",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"shakemap_uz",(/size_surf_array/),1,CUSTOM_REAL)
    else
      call h5_create_dataset_gen_in_group(h5,"shakemap_ux",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"shakemap_uy",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group(h5,"shakemap_uz",(/len_array_aug/),1,CUSTOM_REAL)
   endif
    call h5_close_group(h5)
    call h5_close_file(h5)
  endif

  call synchronize_all()
  call h5_open_file_p(h5)
  call h5_open_group(h5, group_name)

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "shakemap_ux"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, shakemap_ux, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "shakemap_uy"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, shakemap_uy, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "shakemap_uz"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, shakemap_uz, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)

  else
    nfaces_actual      = size(shakemap_ux)/(NGLLX*NGLLY)
    len_array_aug_proc = nfaces_actual*nfaces_aug*nnodes_per_face_aug ! augmented array length for each proc
    aug_factor         = real(len_array_aug)/real(size_surf_array)

    allocate(surf_ux_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uy_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uz_aug(len_array_aug_proc),stat=ier)

    dset_name = "shakemap_ux"
    call recompose_for_hires(shakemap_ux, surf_ux_aug)
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            surf_ux_aug, (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "shakemap_uy"
    call recompose_for_hires(shakemap_uy, surf_uy_aug)
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            surf_uy_aug, (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "shakemap_uz"
    call recompose_for_hires(shakemap_uz, surf_uz_aug)
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            surf_uz_aug, (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)

    deallocate(surf_ux_aug, surf_uy_aug, surf_uz_aug, stat= ier)
  endif

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine write_shakemap_noserv

!=====================================================================

subroutine wmo_movie_volume_output_h5()

! outputs movie files for div, curl and velocity

  use specfem_par
  use constants, only: dest_ionod
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic
  use specfem_par_movie
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB)               :: d_p,div_glob
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: veloc_element
  ! divergence and curl only in the global nodes
!  real(kind=CUSTOM_REAL),dimension(:),allocatable          :: div_glob
  integer,dimension(:),allocatable                         :: valence
  integer                                                  :: ispec,ier,iglob,req
  character(len=3)                                         :: channel
  character(len=1)                                         :: compx,compy,compz
  character(len=MAX_STRING_LEN)                            :: outputname

  integer :: iproc
  logical :: pressure_io=.false., divglob_io=.false., div_io=.false., &
             veloc_io=.false., curl_io=.false.

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif
  integer :: req_count,ireq
  req_count=1

  if (NIONOD == 0 .and. it == NTSTEP_BETWEEN_FRAMES) then
    call prepare_vol_movie_noserv()
  endif

    ! gets component characters: X/Y/Z or E/N/Z
  call write_channel_name(1,channel)
  compx(1:1) = channel(3:3) ! either X or E
  call write_channel_name(2,channel)
  compy(1:1) = channel(3:3) ! either Y or N
  call write_channel_name(3,channel)
  compz(1:1) = channel(3:3) ! Z

  ! saves velocity here to avoid static offset on displacement for movies
  velocity_x_on_node(:) = 0._CUSTOM_REAL
  velocity_y_on_node(:) = 0._CUSTOM_REAL
  velocity_z_on_node(:) = 0._CUSTOM_REAL

  if (ACOUSTIC_SIMULATION) then
    ! uses velocity_x,.. as temporary arrays to store velocity on all GLL points
    do ispec = 1,NSPEC_AB
      ! only acoustic elements
      if (.not. ispec_is_acoustic(ispec)) cycle
      ! calculates velocity
      call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,veloc_element)
      !velocity_x(:,:,:,ispec) = veloc_element(1,:,:,:)
      !velocity_y(:,:,:,ispec) = veloc_element(2,:,:,:)
      !velocity_z(:,:,:,ispec) = veloc_element(3,:,:,:)
      call elm2node_base(ispec,veloc_element, velocity_x_on_node, velocity_y_on_node, velocity_z_on_node)
    enddo

    ! outputs pressure field for purely acoustic simulations
    if (.not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) then
      d_p(:) = 0._CUSTOM_REAL
      do ispec = 1,NSPEC_AB
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          d_p(iglob) = - potential_dot_dot_acoustic(iglob)
        ENDDO_LOOP_IJK
      enddo

      if (NIONOD > 0) then
        ! send pressure_loc
        call isend_cr_inter(d_p,NGLOB_AB,dest_ionod,io_tag_vol_pres,req_dump_vol(req_count))
        req_count = req_count+1
      else
        ! direct io
        pressure_io=.true.
        call write_vol_data_noserv(d_p,"pressure")
      endif
    endif
  endif ! acoustic

  if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
    ! allocate array for global points
    allocate(valence(NGLOB_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2005')
    if (ier /= 0) stop 'error allocating arrays for movie div and curl'

    ! saves full snapshot data to local disk
    if (ELASTIC_SIMULATION) then
      ! calculates divergence and curl of velocity field
      call wmo_movie_div_curl_h5(veloc, &
                              div_glob,valence, &
                              div_on_node,curl_x_on_node,curl_y_on_node,curl_z_on_node, &
                              velocity_x_on_node,velocity_y_on_node,velocity_z_on_node, &
                              ispec_is_elastic)

      if (NIONOD > 0) then
        ! send div_glob
        call isend_cr_inter(div_glob,NGLOB_AB,dest_ionod,io_tag_vol_divglob,req_dump_vol(req_count))
        req_count = req_count+1
      else
        ! direct io
        divglob_io=.true.
        call write_vol_data_noserv(div_glob,"div_glob")
      endif
    endif ! elastic

    ! saves full snapshot data to local disk
    if (POROELASTIC_SIMULATION) then
      ! calculates divergence and curl of velocity field
      call wmo_movie_div_curl_h5(velocs_poroelastic, &
                              div_glob,valence, &
                              div_on_node,curl_x_on_node,curl_y_on_node,curl_z_on_node, &
                              velocity_x_on_node,velocity_y_on_node,velocity_z_on_node, &
                              ispec_is_poroelastic)
    endif ! poroelastic

    deallocate(valence)

    ! div and curl on elemental level
    ! writes our divergence
    if (NIONOD > 0) then
      call isend_cr_inter(div_on_node,NGLOB_AB,dest_ionod,io_tag_vol_div,req_dump_vol(req_count))
      req_count = req_count+1

      ! writes out curl
      call isend_cr_inter(curl_x_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curlx,req_dump_vol(req_count))
      req_count = req_count+1
      call isend_cr_inter(curl_y_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curly,req_dump_vol(req_count))
      req_count = req_count+1
      call isend_cr_inter(curl_z_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curlz,req_dump_vol(req_count))
      req_count = req_count+1
    else
      ! direct io
      div_io=.true.
      curl_io=.true.
      call write_vol_data_noserv(div_on_node,"div")
      call write_vol_data_noserv(curl_x_on_node,"curl_x")
      call write_vol_data_noserv(curl_y_on_node,"curl_y")
      call write_vol_data_noserv(curl_z_on_node,"curl_z")
    endif
  endif ! elastic or poroelastic

  ! velocity
  if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
    if (NIONOD > 0) then
      call isend_cr_inter(velocity_x_on_node,NGLOB_AB,dest_ionod,io_tag_vol_velox,req_dump_vol(req_count))
      req_count = req_count+1
      call isend_cr_inter(velocity_y_on_node,NGLOB_AB,dest_ionod,io_tag_vol_veloy,req_dump_vol(req_count))
      req_count = req_count+1
      call isend_cr_inter(velocity_z_on_node,NGLOB_AB,dest_ionod,io_tag_vol_veloz,req_dump_vol(req_count))
      req_count = req_count+1
    else
      ! direct io
      veloc_io=.true.
      call write_vol_data_noserv(velocity_x_on_node,"velo_x")
      call write_vol_data_noserv(velocity_y_on_node,"velo_y")
      call write_vol_data_noserv(velocity_z_on_node,"velo_z")
    endif
  endif

  ! store the number of mpi_isend reqs
  n_req_vol = req_count-1

  if (it==NTSTEP_BETWEEN_FRAMES .and. myrank==0 .and. NIONOD == 0) then
    ! create xdmf header
    call write_xdmf_vol_noserv(pressure_io, divglob_io, div_io, veloc_io, curl_io)
  endif

end subroutine wmo_movie_volume_output_h5


subroutine prepare_vol_movie_noserv()
  use specfem_par
  use specfem_par_movie

  use io_server
  use phdf5_utils

  implicit none

  integer, dimension(9,NSPEC_AB*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)) :: elm_conn_loc

  !character(len=64) :: fname_h5_data_vol
  character(len=64) :: group_name
  character(len=64) :: dset_name
  type(h5io)        :: h5
  logical           :: if_collect = .true. ! in io_server.f90
  integer           :: iproc, ier, nglob_all, nspec_all

  allocate(nglob_par_proc_nio(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nglob_par_proc')
  if (ier /= 0) stop 'error allocating arrays for nglob_par_proc'
   allocate(nelm_par_proc_nio(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nelm_par_proc')
  if (ier /= 0) stop 'error allocating arrays for nelm_par_proc'
  allocate(nglob_offset(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nglob_offset')
  if (ier /= 0) stop 'error allocating arrays for nglob_offset'
  allocate(nelm_offset(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nelm_offset')
  if (ier /= 0) stop 'error allocating arrays for nelm_offset'



  ! initialize h5 object
  h5 = h5io()
  fname_h5_data_vol = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_volume.h5"
  call h5_init(h5, fname_h5_data_vol)

  ! group for storing node coordinates and mesh element connectivity
  group_name = "mesh"

  call gather_all_all_singlei((/NSPEC_AB/),nelm_par_proc_nio,NPROC)
  call gather_all_all_singlei((/NGLOB_AB/),nglob_par_proc_nio,NPROC)

  nglob_all = sum(nglob_par_proc_nio(0:NPROC-1))
  nspec_all = sum(nelm_par_proc_nio(0:NPROC-1))

  nglob_offset(0) = 0
  nelm_offset(0)  = 0
  do iproc = 1, NPROC-1
    nglob_offset(iproc) = sum(nglob_par_proc_nio(0:iproc-1))
    nelm_offset(iproc)  = sum(nelm_par_proc_nio(0:iproc-1))
  enddo

  do iproc=1,NPROC-1
    nelm_offset(iproc) = nelm_offset(iproc)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)
  enddo

  ! create connectivity dataset
  call get_conn_for_movie(NSPEC_AB, elm_conn_loc, nglob_offset(myrank))

  if(myrank==0) then

   ! create a hdf5 file
    call h5_create_file(h5)

    ! create coordinate dataset
    call h5_open_or_create_group(h5, group_name)

    dset_name = "elm_conn"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/9,nspec_all*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)/), 2, 1)
    dset_name = "x"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all/), 1, CUSTOM_REAL)
    dset_name = "y"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all/), 1, CUSTOM_REAL)
    dset_name = "z"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all/), 1, CUSTOM_REAL)

    call h5_close_group(h5)
    call h5_close_file(h5)
  endif

  call h5_open_file_p(h5)
  call h5_open_group(h5,group_name)


  dset_name = "elm_conn"
  call h5_write_dataset_2d_i_collect_hyperslab_in_group(h5,dset_name,elm_conn_loc,(/0,nelm_offset(myrank)/),if_collect)
  dset_name = "x"
  call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,xstore,(/nglob_offset(myrank)/),if_collect)
  dset_name = "y"
  call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,ystore,(/nglob_offset(myrank)/),if_collect)
  dset_name = "z"
  call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,zstore,(/nglob_offset(myrank)/),if_collect)

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine prepare_vol_movie_noserv


subroutine write_vol_data_noserv(darr, dset_name)
  use specfem_par
  use specfem_par_movie
  use phdf5_utils

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: darr
  character(len=*), intent(in)                     :: dset_name
  character(len=64)                                :: fname_h5_data_vol, group_name
  character(len=12)                                :: tempstr
  type(h5io)                                       :: h5
  integer                                          :: rank=1, dim
  logical                                          :: if_collective = .true.

  h5 = h5io()
  fname_h5_data_vol = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_volume.h5"
  call h5_init(h5, fname_h5_data_vol)

  write(tempstr, "(i6.6)") it
  group_name = "it_"//tempstr

  dim = sum(nglob_par_proc_nio(:))

  if (myrank==0) then
    ! create it group
    call h5_open_file(h5)
    ! check if group_name exists
    call h5_open_or_create_group(h5, group_name)

    ! create dataset
    call h5_create_dataset_gen_in_group(h5, dset_name, (/dim/), rank, CUSTOM_REAL)
    call h5_close_group(h5)
    call h5_close_file(h5)
  endif
  call synchronize_all()

  ! collective write
  call h5_open_file_p(h5)
  call h5_open_group(h5,group_name)
  call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
            darr, (/nglob_offset(myrank)/), if_collective)

  call h5_close_group(h5)
  call h5_close_file(h5)
end subroutine write_vol_data_noserv


! this h5 version returns the arrays in nodes base in stead of element base
! by the requirement of xdmf3
subroutine wmo_movie_div_curl_h5(veloc, &
                                div_glob,valence, &
                                div_on_node,curl_x_on_node,curl_y_on_node,curl_z_on_node, &
                                velocity_x_on_node,velocity_y_on_node,velocity_z_on_node, &
                                ispec_is)

! calculates div, curl and velocity

  use constants

  use specfem_par, only: NSPEC_AB,NGLOB_AB,ibool,hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,irregular_element_number, &
                          xix_regular

  implicit none


  ! velocity field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB),intent(in) :: veloc

  ! divergence and curl only in the global nodes
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(out) :: div_glob
  integer,dimension(NGLOB_AB),intent(out) :: valence

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(out) :: div_on_node, curl_x_on_node, curl_y_on_node, curl_z_on_node
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: velocity_x_on_node,velocity_y_on_node,velocity_z_on_node
  logical,dimension(NSPEC_AB),intent(in) :: ispec_is

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: veloc_element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dvxdxl,dvxdyl, &
                                dvxdzl,dvydxl,dvydyl,dvydzl,dvzdxl,dvzdyl,dvzdzl
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l
  integer :: ispec,ispec_irreg,i,j,k,l,iglob

  ! initializes
  div_glob(:) = 0.0_CUSTOM_REAL
  valence(:) = 0

  ! loops over elements
  do ispec = 1,NSPEC_AB
    ! only selected elements
    if (.not. ispec_is(ispec)) cycle

    ispec_irreg = irregular_element_number(ispec)

    ! loads veloc to local element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          veloc_element(1,i,j,k) = veloc(1,iglob)
          veloc_element(2,i,j,k) = veloc(2,iglob)
          veloc_element(3,i,j,k) = veloc(3,iglob)
        enddo
      enddo
    enddo

    ! calculates divergence and curl of velocity field
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here

!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called
!! DK DK Oct 2018: however this curl and div calculation routine for movies is almost never called

          do l = 1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1l = tempx1l + veloc_element(1,l,j,k)*hp1
            tempy1l = tempy1l + veloc_element(2,l,j,k)*hp1
            tempz1l = tempz1l + veloc_element(3,l,j,k)*hp1

            hp2 = hprime_yy(j,l)
            tempx2l = tempx2l + veloc_element(1,i,l,k)*hp2
            tempy2l = tempy2l + veloc_element(2,i,l,k)*hp2
            tempz2l = tempz2l + veloc_element(3,i,l,k)*hp2

            hp3 = hprime_zz(k,l)
            tempx3l = tempx3l + veloc_element(1,i,j,l)*hp3
            tempy3l = tempy3l + veloc_element(2,i,j,l)*hp3
            tempz3l = tempz3l + veloc_element(3,i,j,l)*hp3
          enddo

          if (ispec_irreg /= 0) then ! irregular element
            ! get derivatives of ux, uy and uz with respect to x, y and z
            xixl = xix(i,j,k,ispec_irreg)
            xiyl = xiy(i,j,k,ispec_irreg)
            xizl = xiz(i,j,k,ispec_irreg)
            etaxl = etax(i,j,k,ispec_irreg)
            etayl = etay(i,j,k,ispec_irreg)
            etazl = etaz(i,j,k,ispec_irreg)
            gammaxl = gammax(i,j,k,ispec_irreg)
            gammayl = gammay(i,j,k,ispec_irreg)
            gammazl = gammaz(i,j,k,ispec_irreg)

            dvxdxl(i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
            dvxdyl(i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
            dvxdzl(i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

            dvydxl(i,j,k) = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
            dvydyl(i,j,k) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
            dvydzl(i,j,k) = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

            dvzdxl(i,j,k) = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
            dvzdyl(i,j,k) = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
            dvzdzl(i,j,k) = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

          else ! regular element

            dvxdxl(i,j,k) = xix_regular*tempx1l
            dvxdyl(i,j,k) = xix_regular*tempx2l
            dvxdzl(i,j,k) = xix_regular*tempx3l

            dvydxl(i,j,k) = xix_regular*tempy1l
            dvydyl(i,j,k) = xix_regular*tempy2l
            dvydzl(i,j,k) = xix_regular*tempy3l

            dvzdxl(i,j,k) = xix_regular*tempz1l
            dvzdyl(i,j,k) = xix_regular*tempz2l
            dvzdzl(i,j,k) = xix_regular*tempz3l

          endif

        enddo
      enddo
    enddo

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool(i,j,k,ispec)

          ! divergence \nabla \cdot \bf{v}
          div_on_node(iglob) = dvxdxl(i,j,k) + dvydyl(i,j,k) + dvzdzl(i,j,k)

          ! curl
          curl_x_on_node(iglob) = dvzdyl(i,j,k) - dvydzl(i,j,k)
          curl_y_on_node(iglob) = dvxdzl(i,j,k) - dvzdxl(i,j,k)
          curl_z_on_node(iglob) = dvydxl(i,j,k) - dvxdyl(i,j,k)

          ! velocity field
          velocity_x_on_node(iglob) = veloc_element(1,i,j,k)
          velocity_y_on_node(iglob) = veloc_element(2,i,j,k)
          velocity_z_on_node(iglob) = veloc_element(3,i,j,k)

          valence(iglob) = valence(iglob)+1
          div_glob(iglob) = div_glob(iglob) + div_on_node(iglob)
        enddo
      enddo
    enddo
  enddo !NSPEC_AB

  do i = 1,NGLOB_AB
    ! checks if point has a contribution
    ! note: might not be the case for points in acoustic elements
    if (valence(i) /= 0) then
      ! averages by number of contributions
      div_glob(i) = div_glob(i)/valence(i)
    endif
  enddo

  end subroutine wmo_movie_div_curl_h5


  subroutine elm2node_base(ispec,elm_base, node_base_x, node_base_y, node_base_z)
    use specfem_par
    use specfem_par_movie

    implicit none

    integer, intent(in) :: ispec
    real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: elm_base
    real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(inout) :: node_base_x, node_base_y, node_base_z
    integer :: i,j,k,iglob

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          ! velocity field
          node_base_x(iglob) = elm_base(1,i,j,k)
          node_base_y(iglob) = elm_base(2,i,j,k)
          node_base_z(iglob) = elm_base(3,i,j,k)
        enddo
      enddo
    enddo
  end subroutine elm2node_base



  subroutine get_conn_for_movie(nspec,elm_conn,o)
    use specfem_par
    implicit none

    integer, intent(in)                                                    :: nspec
    integer, dimension(9,nspec*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)), intent(out) :: elm_conn
    integer, intent(in) :: o ! node id offset (starting global element id of each proc)

    integer :: ispec,ii,iglob,icub,jcub,kcub,cell_type=9, dp=2

    do ispec=1, nspec
      ! extract information from full GLL grid
      ! node order follows vtk format

      do icub=0,NGLLX-2
        do jcub=0,NGLLY-2
          do kcub=0,NGLLZ-2
            ii = 1+(ispec-1)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1) + (icub*(NGLLY-1)*(NGLLZ-1)+jcub*(NGLLZ-1)+kcub)
            elm_conn(1, ii)  = cell_type
            elm_conn(2, ii)  = ibool(icub+1,jcub+1,kcub+1,ispec)-1 +o! node id starts 0 in xdmf rule
            elm_conn(3, ii)  = ibool(icub+2,jcub+1,kcub+1,ispec)-1 +o
            elm_conn(4, ii)  = ibool(icub+2,jcub+2,kcub+1,ispec)-1 +o
            elm_conn(5, ii)  = ibool(icub+1,jcub+2,kcub+1,ispec)-1 +o
            elm_conn(6, ii)  = ibool(icub+1,jcub+1,kcub+2,ispec)-1 +o
            elm_conn(7, ii)  = ibool(icub+2,jcub+1,kcub+2,ispec)-1 +o
            elm_conn(8, ii)  = ibool(icub+2,jcub+2,kcub+2,ispec)-1 +o
            elm_conn(9, ii)  = ibool(icub+1,jcub+2,kcub+2,ispec)-1 +o
          enddo
        enddo
      enddo
    enddo

  end subroutine get_conn_for_movie


subroutine write_xdmf_vol_noserv(pressure_io, divglob_io, div_io, veloc_io, curl_io)
  use specfem_par
  use specfem_par_movie
  use io_server
  implicit none

  logical, intent(in) :: pressure_io, divglob_io, div_io, veloc_io, curl_io
  character(len=20)                         :: proc_str, it_str,nelm, nglo
  integer                                   :: iiout, nout, i, ii

  ! writeout xdmf file for volume movie
  fname_xdmf_vol = trim(OUTPUT_FILES)//"/movie_volume.xmf"

  open(unit=xdmf_vol, file=fname_xdmf_vol, recl=256)

  ! definition of topology and geometry
  ! refer only control nodes (8 or 27) as a coarse output
  ! data array need to be extracted from full data array on gll points
  write(xdmf_vol,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_vol,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_vol,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
  write(xdmf_vol,*) '<Domain name="mesh">'
  ! loop for writing information of mesh partitions
  nelm=i2c(sum(nelm_par_proc_nio(:))*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))
  nglo=i2c(sum(nglob_par_proc_nio(:)))

  write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm)//'">'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                         //trim(nelm)//' 9">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/elm_conn'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '</Topology>'
  write(xdmf_vol,*) '<Geometry GeometryType="X_Y_Z">'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/x'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/y'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/z'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '</Geometry>'
  write(xdmf_vol,*) '<Grid Name="time_col" GridType="Collection" CollectionType="Temporal">'

  do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ii = i*NTSTEP_BETWEEN_FRAMES
    write(it_str, "(i6.6)") ii


    write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
    write(xdmf_vol,*) '  <Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
    write(xdmf_vol,*) '  <Topology Reference="/Xdmf/Domain/Topology" />'
    write(xdmf_vol,*) '  <Geometry Reference="/Xdmf/Domain/Geometry" />'

    if (pressure_io) then
      write(xdmf_vol,*) '  <Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/pressure'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (divglob_io) then
      write(xdmf_vol,*) '  <Attribute Name="div_glob" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/div_glob'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (div_io) then
      write(xdmf_vol,*) '  <Attribute Name="div" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/div'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (veloc_io) then
      write(xdmf_vol,*) '  <Attribute Name="velo_x" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/velo_x'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="velo_y" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/velo_y'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="velo_z" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/velo_z'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (curl_io) then
      write(xdmf_vol,*) '  <Attribute Name="curl_x" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/curl_x'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="curl_y" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/curl_y'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="curl_z" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/curl_z'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    write(xdmf_vol,*) '</Grid>'
  enddo

  write(xdmf_vol,*) '</Grid>'
  write(xdmf_vol,*) '</Domain>'
  write(xdmf_vol,*) '</Xdmf>'

  close(xdmf_vol)

end subroutine write_xdmf_vol_noserv


