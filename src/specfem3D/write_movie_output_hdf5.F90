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

  ! saves MOVIE on the SURFACE
  if (MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    call wmo_movie_surface_output_h5()
  endif

  ! saves MOVIE in full 3D MESH
  if (MOVIE_VOLUME .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    call wmo_movie_volume_output_h5()
  endif

  ! creates cross-section PNM image
  if (PNM_IMAGE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    call write_PNM_create_image() ! this is excluded from this file. need to check
  endif

  end subroutine write_movie_output_h5



!================================================================

  subroutine wmo_movie_surface_output_h5()

! output of surface moviedata files
! IO is done by only a master process

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
  integer :: ispec2D,ispec,ipoin,iglob,ier,ia,req
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

 ! send surface body to io node
  call isend_cr_inter(store_val_ux,nfaces_surface_points,0,io_tag_surface_ux,req)
  call isend_cr_inter(store_val_uy,nfaces_surface_points,0,io_tag_surface_uy,req)
  call isend_cr_inter(store_val_uz,nfaces_surface_points,0,io_tag_surface_uz,req)

end subroutine wmo_movie_surface_output_h5


subroutine wmo_create_shakemap_h5()

! creation of shapemap file
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
  integer :: ier, req
  real(kind=CUSTOM_REAL),dimension(1):: dummy

 ! send surface body to io node
  call isend_cr_inter(shakemap_ux,nfaces_surface_points,0,io_tag_shake_ux,req)
  call isend_cr_inter(shakemap_uy,nfaces_surface_points,0,io_tag_shake_uy,req)
  call isend_cr_inter(shakemap_uz,nfaces_surface_points,0,io_tag_shake_uz,req)

end subroutine wmo_save_shakemap_h5


!=====================================================================

subroutine wmo_movie_volume_output_h5()

! outputs movie files for div, curl and velocity

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic
  use specfem_par_movie
  use phdf5_utils
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB)               :: send_array
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: veloc_element
  ! divergence and curl only in the global nodes
  real(kind=CUSTOM_REAL),dimension(:),allocatable          :: div_glob
  integer,dimension(:),allocatable                         :: valence
  integer                                                  :: ispec,ier,iglob,req,arrsize
  character(len=3)                                         :: channel
  character(len=1)                                         :: compx,compy,compz
  character(len=MAX_STRING_LEN)                            :: outputname
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  send_array(:) = 0._CUSTOM_REAL

  ! gets component characters: X/Y/Z or E/N/Z
  call write_channel_name(1,channel)
  compx(1:1) = channel(3:3) ! either X or E
  call write_channel_name(2,channel)
  compy(1:1) = channel(3:3) ! either Y or N
  call write_channel_name(3,channel)
  compz(1:1) = channel(3:3) ! Z

  ! saves velocity here to avoid static offset on displacement for movies
  velocity_x(:,:,:,:) = 0._CUSTOM_REAL
  velocity_y(:,:,:,:) = 0._CUSTOM_REAL
  velocity_z(:,:,:,:) = 0._CUSTOM_REAL

  if (ACOUSTIC_SIMULATION) then
    ! uses velocity_x,.. as temporary arrays to store velocity on all GLL points
    do ispec = 1,NSPEC_AB
      ! only acoustic elements
      if (.not. ispec_is_acoustic(ispec)) cycle
      ! calculates velocity
      call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,veloc_element)
      velocity_x(:,:,:,ispec) = veloc_element(1,:,:,:)
      velocity_y(:,:,:,ispec) = veloc_element(2,:,:,:)
      velocity_z(:,:,:,ispec) = veloc_element(3,:,:,:)
    enddo

    ! outputs pressure field for purely acoustic simulations
    if (.not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) then
      do ispec = 1,NSPEC_AB
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          send_array(iglob) = - potential_dot_dot_acoustic(iglob)
        ENDDO_LOOP_IJK
      enddo

      ! send pressure_loc
      arrsize = NGLOB_AB
      call isend_cr_inter(send_array,arrsize,0,io_tag_vol_pres,req)
      call wait_req(req)
    endif
  endif ! acoustic

  if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
    ! allocate array for global points
    allocate(div_glob(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2004')
    allocate(valence(NGLOB_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2005')
    if (ier /= 0) stop 'error allocating arrays for movie div and curl'

    ! saves full snapshot data to local disk
    if (ELASTIC_SIMULATION) then
      ! calculates divergence and curl of velocity field
      call wmo_movie_div_curl(veloc, &
                              div_glob,valence, &
                              div,curl_x,curl_y,curl_z, &
                              velocity_x,velocity_y,velocity_z, &
                              ispec_is_elastic)

      ! send div_glob
      arrsize = NGLOB_AB
      call isend_cr_inter(div_glob,arrsize,0,io_tag_vol_divglob,req)
      call wait_req(req)
    endif ! elastic

    ! saves full snapshot data to local disk
    if (POROELASTIC_SIMULATION) then
      ! calculates divergence and curl of velocity field
      call wmo_movie_div_curl(velocs_poroelastic, &
                              div_glob,valence, &
                              div,curl_x,curl_y,curl_z, &
                              velocity_x,velocity_y,velocity_z, &
                              ispec_is_poroelastic)
    endif ! poroelastic

    deallocate(div_glob,valence)

    ! div and curl on elemental level
    ! writes our divergence
    arrsize = NGLOB_AB
    call convert_elmbase_to_nodebase(div,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_div,req)
    call wait_req(req)
 

    ! writes out curl
    call convert_elmbase_to_nodebase(curl_x,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_curlx,req)
    call wait_req(req)
    call convert_elmbase_to_nodebase(curl_y,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_curly,req)
    call wait_req(req)
    call convert_elmbase_to_nodebase(curl_z,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_curlz,req)
    call wait_req(req)
 

  endif

  ! velocity
  if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then

    arrsize = NGLOB_AB
    call convert_elmbase_to_nodebase(velocity_x,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_velox,req)
    call wait_req(req)
    call convert_elmbase_to_nodebase(velocity_y,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_veloy,req)
    call wait_req(req)
    call convert_elmbase_to_nodebase(velocity_z,send_array)
    call isend_cr_inter(send_array,arrsize,0,io_tag_vol_veloz,req)
    call wait_req(req)
 
  endif

end subroutine wmo_movie_volume_output_h5

! # TODO: check if the vectorization loop definition is fine
subroutine convert_elmbase_to_nodebase(arr_in, arr_out)
  use specfem_par

  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB), intent(in) :: arr_in
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(out)                  :: arr_out
  integer                                                                   :: ispec,iglob
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  do ispec=1, NSPEC_AB
    DO_LOOP_IJK
      iglob = ibool(INDEX_IJK,ispec)
      arr_out(iglob) = arr_in(INDEX_IJK,ispec)
    ENDDO_LOOP_IJK
  enddo
end subroutine convert_elmbase_to_nodebase