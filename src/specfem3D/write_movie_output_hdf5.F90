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
    call wait_all_send()

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

 ! send surface body to io node
  call isend_cr_inter(store_val_ux,nfaces_surface_points,0,io_tag_surface_ux,req_dump_surf(1))
  call isend_cr_inter(store_val_uy,nfaces_surface_points,0,io_tag_surface_uy,req_dump_surf(2))
  call isend_cr_inter(store_val_uz,nfaces_surface_points,0,io_tag_surface_uz,req_dump_surf(3))

  n_req_surf = 3

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
  integer :: req

 ! send surface body to io node
  call isend_cr_inter(shakemap_ux,nfaces_surface_points,0,io_tag_shake_ux,req)
  call isend_cr_inter(shakemap_uy,nfaces_surface_points,0,io_tag_shake_uy,req)
  call isend_cr_inter(shakemap_uz,nfaces_surface_points,0,io_tag_shake_uz,req)

end subroutine wmo_save_shakemap_h5


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
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB)               :: d_p
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: veloc_element
  ! divergence and curl only in the global nodes
  real(kind=CUSTOM_REAL),dimension(:),allocatable          :: div_glob
  integer,dimension(:),allocatable                         :: valence
  integer                                                  :: ispec,ier,iglob,req
  character(len=3)                                         :: channel
  character(len=1)                                         :: compx,compy,compz
  character(len=MAX_STRING_LEN)                            :: outputname
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif
  integer :: req_count,ireq
  req_count=1

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
!      velocity_x(:,:,:,ispec) = veloc_element(1,:,:,:)
!      velocity_y(:,:,:,ispec) = veloc_element(2,:,:,:)
!      velocity_z(:,:,:,ispec) = veloc_element(3,:,:,:)
    enddo

    ! outputs pressure field for purely acoustic simulations
    if (.not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) then
      do ispec = 1,NSPEC_AB
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          d_p(iglob) = - potential_dot_dot_acoustic(iglob)
        ENDDO_LOOP_IJK
      enddo

      ! send pressure_loc
      call isend_cr_inter(d_p,NGLOB_AB,dest_ionod,io_tag_vol_pres,req_dump_vol(req_count))
      req_count = req_count+1
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
      call wmo_movie_div_curl_h5(veloc, &
                              div_glob,valence, &
                              div_on_node,curl_x_on_node,curl_y_on_node,curl_z_on_node, &
                              velocity_x_on_node,velocity_y_on_node,velocity_z_on_node, &
                              ispec_is_elastic)

      ! send div_glob
      call isend_cr_inter(div_glob,NGLOB_AB,dest_ionod,io_tag_vol_divglob,req_dump_vol(req_count))
      req_count = req_count+1
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

    deallocate(div_glob,valence)

    ! div and curl on elemental level
    ! writes our divergence
    call isend_cr_inter(div_on_node,NGLOB_AB,dest_ionod,io_tag_vol_div,req_dump_vol(req_count))
    req_count = req_count+1
 
    ! writes out curl
    call isend_cr_inter(curl_x_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curlx,req_dump_vol(req_count))
    req_count = req_count+1
    call isend_cr_inter(curl_y_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curly,req_dump_vol(req_count))
    req_count = req_count+1
    call isend_cr_inter(curl_z_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curlz,req_dump_vol(req_count))
    req_count = req_count+1
  endif

  ! velocity
  if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
    call isend_cr_inter(velocity_x_on_node,NGLOB_AB,dest_ionod,io_tag_vol_velox,req_dump_vol(req_count))
    req_count = req_count+1
    call isend_cr_inter(velocity_y_on_node,NGLOB_AB,dest_ionod,io_tag_vol_veloy,req_dump_vol(req_count))
    req_count = req_count+1
    call isend_cr_inter(velocity_z_on_node,NGLOB_AB,dest_ionod,io_tag_vol_veloz,req_dump_vol(req_count))
    req_count = req_count+1
  endif

  ! store the number of mpi_isend reqs
  n_req_vol = req_count-1

end subroutine wmo_movie_volume_output_h5


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

