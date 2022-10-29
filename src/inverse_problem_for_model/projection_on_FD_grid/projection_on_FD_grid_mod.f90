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

module projection_on_FD_grid

  !! IMPORT SPECFEM VARIABLES
  use specfem_par, only: NGLLX, NGLLY, NGLLZ, NDIM, NSPEC_AB, &
       NGLOB_AB, ibool, xstore, ystore, zstore,  NUM_ITER, NGNOD, xigll, yigll, zigll, NPROC, HUGEVAL, &
       MIDX, MIDY, MIDZ


  !! IMPORT INVERSE_PROBLEM VARIABLES
  use inverse_problem_par

  implicit none

  !! fd grid parameters
  integer,                public                                   :: nx_fd_proj, ny_fd_proj, nz_fd_proj
  real(kind=CUSTOM_REAL), public                                   :: hx_fd_proj, hy_fd_proj, hz_fd_proj
  real(kind=CUSTOM_REAL), public                                   :: ox_fd_proj, oy_fd_proj, oz_fd_proj

  !! locals
  integer,                private                                  :: ifd, jfd, kfd
  integer,                private                                  :: igg, jgg, kgg
  real(kind=CUSTOM_REAL), private                                  :: xfd, yfd, zfd

  real(kind=CUSTOM_REAL), private, dimension(:,:,:), allocatable   :: valence , valence_tmp, model_on_FD_grid_tmp

contains

!--------------------------------------------------------------------------------------------------------------------
!  Projection FD grid model on SEM mesh
!--------------------------------------------------------------------------------------------------------------------

  subroutine Project_model_FD_grid2SEM(model_on_SEM_mesh, model_on_FD_grid, myrank)

  integer,                                                 intent(in)    :: myrank
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(inout) :: model_on_SEM_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:),   allocatable, intent(inout) :: model_on_FD_grid

  integer                                                :: ispec, iglob, i_fd, j_fd, k_fd
  real(kind=CUSTOM_REAL)                                 :: vinterp, v1, v2, v3, v4, v5, v6, v7, v8
  real(kind=CUSTOM_REAL)                                 :: x_sem, y_sem, z_sem
  real(kind=CUSTOM_REAL)                                 :: x_loc, y_loc, z_loc


  if (myrank == 0) then
     write(INVERSE_LOG_FILE,*)
     write(INVERSE_LOG_FILE,*) '    - > projection FD grid / SEM mesh  '
     write(INVERSE_LOG_FILE,*)
  endif


  do ispec =1, NSPEC_AB

     do kgg = 1, NGLLZ
        do jgg = 1, NGLLY
           do igg = 1, NGLLX

              iglob = ibool(igg,jgg,kgg,ispec)
              x_sem = xstore(iglob)
              y_sem = ystore(iglob)
              z_sem = zstore(iglob)

              !! get value of model on FD grid by trilinear interpolation
              i_fd = floor((x_sem -  ox_fd_proj)/ hx_fd_proj) + 1
              j_fd = floor((y_sem -  oy_fd_proj)/ hy_fd_proj) + 1
              k_fd = floor((z_sem -  oz_fd_proj)/ hz_fd_proj) + 1

              if (i_fd <= nx_fd_proj .and. j_fd <= ny_fd_proj .and. k_fd <= nz_fd_proj .and. &
                   i_fd > 0 .and. j_fd > 0 .and. k_fd > 0          ) then

                 v1 = model_on_FD_grid(i_fd,     j_fd,     k_fd    )
                 v2 = v1; v3 = v1; v4 = v1; v5 = v1; v6 = v1; v7 = v1; v8 = v1

                 if (i_fd < nx_fd_proj) v2 = model_on_FD_grid(i_fd + 1, j_fd,     k_fd    )
                 if (i_fd < nx_fd_proj .and. j_fd < ny_fd_proj) v3 = model_on_FD_grid(i_fd + 1, j_fd + 1, k_fd    )
                 if (j_fd < ny_fd_proj) v4 = model_on_FD_grid(i_fd,     j_fd + 1, k_fd    )

                 if (k_fd < nz_fd_proj) then
                    v5 = model_on_FD_grid(i_fd,     j_fd,     k_fd + 1)
                    if (i_fd < nx_fd_proj) v6 = model_on_FD_grid(i_fd + 1, j_fd,     k_fd + 1)
                    if (i_fd < nx_fd_proj .and. j_fd < ny_fd_proj)   v7 = model_on_FD_grid(i_fd + 1, j_fd + 1, k_fd + 1)
                    if (j_fd < ny_fd_proj) v8 = model_on_FD_grid(i_fd,     j_fd + 1, k_fd + 1)
                 endif

                 x_loc = x_sem - (ox_fd_proj + real( i_fd - 1, CUSTOM_REAL) * hx_fd_proj)
                 y_loc = y_sem - (oy_fd_proj + real( j_fd - 1, CUSTOM_REAL) * hy_fd_proj)
                 z_loc = z_sem - (oz_fd_proj + real( k_fd - 1, CUSTOM_REAL) * hz_fd_proj)

                 call TrilinearInterp(Vinterp, x_loc, y_loc, z_loc, v1, v2, v3, v4, v5, v6, v7, v8, &
                                      hx_fd_proj, hy_fd_proj, hz_fd_proj)
              else

                 Vinterp = 0.

              endif

              model_on_SEM_mesh(igg,jgg,kgg,ispec)=Vinterp

           enddo
        enddo
     enddo
  enddo

  end subroutine Project_model_FD_grid2SEM


!--------------------------------------------------------------------------------------------------------------------
!  read fd grid parameters !!
!--------------------------------------------------------------------------------------------------------------------

  subroutine read_fd_grid_parameters_for_projection()

  integer ier

  open(676,file='fd_proj_grid.txt',iostat=ier)

  if (ier /= 0) then
    print *,'Error opening file "fd_proj_grid.txt", that defines properties of the regular grid'
    stop
  endif

  ! read fd grid parameters !!
  read(676, *)  ox_fd_proj, oy_fd_proj, oz_fd_proj
  read(676, *)  hx_fd_proj, hy_fd_proj, hz_fd_proj
  read(676, *)  nx_fd_proj, ny_fd_proj, nz_fd_proj
  close(676)

  end subroutine read_fd_grid_parameters_for_projection


!--------------------------------------------------------------------------------------------------------------------
!  pre-compute interpolation coeffs for projection on FD GRID
!--------------------------------------------------------------------------------------------------------------------

  subroutine compute_interpolation_coeff_FD_SEM(projection_fd, myrank)


    type(profd)  ,                          intent(inout)  :: projection_fd
    integer,                                intent(in)     :: myrank

    integer                                                :: ispec_selected, ispec, iglob
    double precision                                       :: xi_loc, eta_loc, gamma_loc
    double precision                                       :: x_found,  y_found,  z_found
    double precision                                       :: x_to_locate, y_to_locate, z_to_locate
    real(kind=CUSTOM_REAL)                                 :: distance_min_glob,distance_max_glob
    real(kind=CUSTOM_REAL)                                 :: elemsize_min_glob,elemsize_max_glob
    real(kind=CUSTOM_REAL)                                 :: x_min_glob, x_max_glob
    real(kind=CUSTOM_REAL)                                 :: y_min_glob, y_max_glob
    real(kind=CUSTOM_REAL)                                 :: z_min_glob, z_max_glob
    real(kind=CUSTOM_REAL)                                 :: xmin, xmax
    real(kind=CUSTOM_REAL)                                 :: ymin, ymax
    real(kind=CUSTOM_REAL)                                 :: zmin, zmax
    integer                                                :: imin, imax
    integer                                                :: jmin, jmax
    integer                                                :: kmin, kmax
    integer,                 dimension(NGNOD)              :: iaddx, iaddy, iaddz
    double precision,        dimension(NGLLX)              :: hxis, hpxis
    double precision,        dimension(NGLLY)              :: hetas, hpetas
    double precision,        dimension(NGLLZ)              :: hgammas, hpgammas
    integer                                                :: nb_fd_point_loc, ier
    integer, dimension(:,:,:), allocatable                 :: point_already_found
    double precision, dimension(:,:,:), allocatable        :: xi_in_fd, eta_in_fd, gamma_in_fd

    ! read fd grid parameters !!
    call read_fd_grid_parameters_for_projection()

    ! get mesh properties (mandatory before calling locate_point_in_mesh_simple)
    call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

    nb_fd_point_loc = 0

    allocate(point_already_found(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 144')
    allocate(xi_in_fd(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 145')
    allocate(eta_in_fd(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 146')
    allocate(gamma_in_fd(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 147')

    point_already_found(:,:,:) = 0

    !! loop over elements
    do ispec = 1, NSPEC_AB
!!$       if (DEBUG_MODE) write(IIDD,*) isepc, NSPEC_AB
       xmin =  HUGEVAL
       xmax = -HUGEVAL
       ymin =  HUGEVAL
       ymax = -HUGEVAL
       zmin =  HUGEVAL
       zmax = -HUGEVAL
       !! get boundary to
       do kgg = 1, NGLLZ
          do jgg = 1, NGLLY
             do igg = 1, NGLLX
                iglob = ibool(igg,jgg,kgg,ispec)

                xmin  = min( xmin, xstore(iglob))
                xmax  = max( xmax, xstore(iglob))

                ymin  = min( ymin, ystore(iglob))
                ymax  = max( ymax, ystore(iglob))

                zmin  = min( zmin, zstore(iglob))
                zmax  = max( zmax, zstore(iglob))

             enddo
          enddo
       enddo

       kmin =  1+ int((zmin  - oz_fd_proj) / hz_fd_proj)
       kmax =  1+ int((zmax  - oz_fd_proj) / hz_fd_proj)
       jmin =  1+ int((ymin  - oy_fd_proj) / hy_fd_proj)
       jmax =  1+ int((ymax  - oy_fd_proj) / hy_fd_proj)
       imin =  1+ int((xmin  - ox_fd_proj) / hx_fd_proj)
       imax =  1+ int((xmax  - ox_fd_proj) / hx_fd_proj)

        if (DEBUG_MODE) then
          write(IIDD,*) ' projection SEM2FD : boundary element'
          write(IIDD,*) ispec, NSPEC_AB
          write(IIDD,*) xmin, xmax
          write(IIDD,*) ymin, ymax
          write(IIDD,*) zmin, zmax
          write(IIDD,*)
          write(IIDD,*) imin, imax
          write(IIDD,*) jmin, jmax
          write(IIDD,*) kmin, kmax
          write(IIDD,*)
       endif

       !! loop on fd grid to count the numbers of point in fd grid that live in my mesh partition
       do kfd = kmin, kmax
          zfd = oz_fd_proj + (kfd-1)*hz_fd_proj
          do jfd = jmin, jmax
             yfd = oy_fd_proj + (jfd-1)*hy_fd_proj
             do ifd = imin, imax
                xfd = ox_fd_proj + (ifd-1)*hx_fd_proj


                x_to_locate = xfd
                y_to_locate = yfd
                z_to_locate = zfd

                !! compute element ispec where is the point (xfd, yfd, zfd)
                call locate_point_in_element(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
                     ispec_selected, xi_loc, eta_loc, gamma_loc, x_found, y_found, z_found, myrank, ispec)

!!$                !! compute islice MPI partition where is the point  (xfd, yfd, zfd)
!!$                call locate_MPI_slice_and_bcast_to_all_1(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
!!$                     xi_loc, eta_loc, gamma_loc, ispec_selected, islice_selected,  distance_from_target, myrank)

!!$                if (DEBUG_MODE) then
!!$                   write(IIDD,*)
!!$                   write(IIDD,*)  ifd, jfd, kfd
!!$                   write(IIDD,*)  xi_loc, eta_loc, gamma_loc
!!$                   write(IIDD,*)  point_already_found(ifd, jfd, kfd)
!!$                endif

                if (abs(xi_loc) < 1.05d0 .and. abs(eta_loc) < 1.05d0 .and. abs(gamma_loc) < 1.05d0) then
                   if (point_already_found(ifd, jfd, kfd) == 0) then
                      point_already_found(ifd, jfd, kfd)=ispec_selected
                      xi_in_fd(ifd, jfd, kfd)=xi_loc
                      eta_in_fd(ifd, jfd, kfd)=eta_loc
                      gamma_in_fd(ifd, jfd, kfd)=gamma_loc
                      nb_fd_point_loc = nb_fd_point_loc + 1
                   endif
                endif

             enddo
          enddo
       enddo
    enddo

    if (DEBUG_MODE) write(IIDD,*) ' END  compute_interpolation_coeff_FD_SEM step 1',  nb_fd_point_loc

    !! allocate projection structure
    projection_fd%nx=nx_fd_proj
    projection_fd%ny=ny_fd_proj
    projection_fd%nz=nz_fd_proj

    projection_fd%hx=hx_fd_proj
    projection_fd%hy=hy_fd_proj
    projection_fd%hz=hz_fd_proj

    projection_fd%ox=ox_fd_proj
    projection_fd%oy=oy_fd_proj
    projection_fd%oz=oz_fd_proj

    projection_fd%nb_fd_point=nb_fd_point_loc

    allocate(projection_fd%ispec_selected(nb_fd_point_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 148')
    allocate(projection_fd%index_on_fd_grid(3,nb_fd_point_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 149')
    allocate(projection_fd%hxi(NGLLX,nb_fd_point_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 150')
    allocate(projection_fd%heta(NGLLX,nb_fd_point_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 151')
    allocate(projection_fd%hgamma(NGLLX,nb_fd_point_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 152')

    nb_fd_point_loc=0
    do kfd = 1, nz_fd_proj
       do jfd = 1, ny_fd_proj
          do ifd = 1, nx_fd_proj

             if (point_already_found(ifd, jfd, kfd) > 0 ) then

                 nb_fd_point_loc =  nb_fd_point_loc + 1

                 xi_loc=xi_in_fd(ifd, jfd, kfd)
                 eta_loc=eta_in_fd(ifd, jfd, kfd)
                 gamma_loc=gamma_in_fd(ifd, jfd, kfd)
                 ispec_selected = point_already_found(ifd, jfd, kfd)

                 call lagrange_any(xi_loc, NGLLX, xigll, hxis, hpxis)
                 call lagrange_any(eta_loc, NGLLY, yigll, hetas, hpetas)
                 call lagrange_any(gamma_loc, NGLLZ, zigll, hgammas, hpgammas)

                 projection_fd%ispec_selected(nb_fd_point_loc) = ispec_selected
                 projection_fd%hxi(:,nb_fd_point_loc) = hxis(:)
                 projection_fd%heta(:,nb_fd_point_loc) = hetas(:)
                 projection_fd%hgamma(:,nb_fd_point_loc) = hgammas(:)
                 projection_fd%index_on_fd_grid(1,nb_fd_point_loc)=ifd
                 projection_fd%index_on_fd_grid(2,nb_fd_point_loc)=jfd
                 projection_fd%index_on_fd_grid(3,nb_fd_point_loc)=kfd


              endif

           enddo
        enddo
     enddo

!!$    point_already_found(:,:,:)=0
!!$    !! loop over elements
!!$    do ispec = 1, NSPEC_AB
!!$
!!$       xmin =  HUGEVAL
!!$       xmax = -HUGEVAL
!!$       ymin =  HUGEVAL
!!$       ymax = -HUGEVAL
!!$       zmin =  HUGEVAL
!!$       zmax = -HUGEVAL
!!$       !! get boundary to
!!$       do kgg = 1, NGLLZ
!!$          do jgg = 1, NGLLY
!!$             do igg = 1, NGLLX
!!$                iglob = ibool(igg,jgg,kgg,ispec)
!!$
!!$                xmin  = min( xmin, xstore(iglob))
!!$                xmax  = max( xmax, xstore(iglob))
!!$
!!$                ymin  = min( ymin, ystore(iglob))
!!$                ymax  = max( ymax, ystore(iglob))
!!$
!!$                zmin  = min( zmin, zstore(iglob))
!!$                zmax  = max( zmax, zstore(iglob))
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       kmin =  1+ (zmin  - oz_fd_proj) / hz_fd_proj
!!$       kmax =  1+ (zmax  - oz_fd_proj) / hz_fd_proj
!!$       jmin =  1+ (ymin  - oy_fd_proj) / hy_fd_proj
!!$       jmax =  1+ (ymax  - oy_fd_proj) / hy_fd_proj
!!$       imin =  1+ (xmin  - ox_fd_proj) / hx_fd_proj
!!$       imax =  1+ (xmax  - ox_fd_proj) / hx_fd_proj
!!$
!!$
!!$       !! loop on fd grid to fill proj_fd structure
!!$       do kfd = kmin, kmax
!!$          zfd = oz_fd_proj + (kfd-1)*hz_fd_proj
!!$          do jfd = jmin, jmax
!!$             yfd = oy_fd_proj + (jfd-1)*hy_fd_proj
!!$             do ifd = imin, imax
!!$                xfd = ox_fd_proj + (ifd-1)*hx_fd_proj
!!$
!!$                x_to_locate = xfd
!!$                y_to_locate = yfd
!!$                z_to_locate = zfd
!!$
!!$                !! compute element ispec where is the point (xfd, yfd, zfd)
!!$                call locate_point_in_element(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
!!$                     ispec_selected, xi_loc, eta_loc, gamma_loc, x_found, y_found, z_found, myrank, ispec)
!!$
!!$                !! compute islice MPI partition where is the point  (xfd, yfd, zfd)
!!$                call locate_MPI_slice_and_bcast_to_all_1(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
!!$                     xi_loc, eta_loc, gamma_loc, ispec_selected, islice_selected,  distance_from_target, myrank)
!!$
!!$                if (abs(xi_loc) < 1.05d0 .and. abs(eta_loc) < 1.05d0 .and. abs(gamma_loc) < 1.05d0) then
!!$                   if (islice_selected == myrank .and. point_already_found(ifd, jfd, kfd) == 0) then
!!$                      ! compute Lagrange polynomials at the source location
!!$                      call lagrange_any(xi_loc, NGLLX, xigll, hxis, hpxis)
!!$                      call lagrange_any(eta_loc, NGLLY, yigll, hetas, hpetas)
!!$                      call lagrange_any(gamma_loc, NGLLZ, zigll, hgammas, hpgammas)
!!$
!!$                      nb_fd_point_loc = nb_fd_point_loc + 1
!!$                      point_already_found(ifd, jfd, kfd) = 1
!!$
!!$                      projection_fd%ispec_selected(nb_fd_point_loc) = ispec_selected
!!$                      projection_fd%hxi(:,nb_fd_point_loc) = hxis(:)
!!$                      projection_fd%heta(:,nb_fd_point_loc) = hetas(:)
!!$                      projection_fd%hgamma(:,nb_fd_point_loc) = hgammas(:)
!!$                      projection_fd%index_on_fd_grid(1,nb_fd_point_loc)=ifd
!!$                      projection_fd%index_on_fd_grid(2,nb_fd_point_loc)=jfd
!!$                      projection_fd%index_on_fd_grid(3,nb_fd_point_loc)=kfd
!!$
!!$                if (DEBUG_MODE) then
!!$                   write(IIDD, *)
!!$                   write(IIDD, *) ' point :',  nb_fd_point_loc
!!$                   write(IIDD, *) ' index :', ifd, jfd, kfd
!!$                   write(IIDD, *) 'x :', x_to_locate,  x_found
!!$                   write(IIDD, *) 'y :', y_to_locate,  y_found
!!$                   write(IIDD, *) 'z :', z_to_locate,  z_found
!!$                   write(IIDD, *) 'dist ', distance_from_target, islice_selected
!!$                   write(IIDD, *)
!!$                endif
!!$
!!$                   endif
!!$                endif
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo

    if (DEBUG_MODE) write(IIDD,*) ' END  compute_interpolation_coeff_FD_SEM step 2', nb_fd_point_loc
    !! on pourrait faire un all_reduce de point_already_found pour voir si on les a tous trouves et s'il y a pas de doublons
!!$    if (DEBUG_MODE) then
!!$       write(IIDD, * )
!!$       do kfd = 1, nz_fd_proj
!!$          do jfd = 1, ny_fd_proj
!!$             do ifd = 1, nx_fd_proj
!!$                if (point_already_found(ifd, jfd, kfd) > 0) then
!!$                   write(IIDD, * ) ifd, jfd, kfd
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo
!!$       write(IIDD, * )
!!$    endif

    deallocate(point_already_found)
    deallocate( xi_in_fd, eta_in_fd, gamma_in_fd)

    if (DEBUG_MODE) write(IIDD,*) ' END  compute_interpolation_coeff_FD_SEM'

  end subroutine compute_interpolation_coeff_FD_SEM


!--------------------------------------------------------------------------------------------------------------------
!  perform projection on FD grid
!--------------------------------------------------------------------------------------------------------------------

  subroutine Project_model_SEM2FD_grid(model_on_SEM_mesh, model_on_FD_grid, projection_fd, myrank)

     type(profd),                                               intent(in)     :: projection_fd
     integer,                                                   intent(in)     :: myrank
     real(kind=CUSTOM_REAL),   dimension(:,:,:,:), allocatable, intent(in)     :: model_on_SEM_mesh
     real(kind=CUSTOM_REAL),   dimension(:,:,:),   allocatable, intent(inout)  :: model_on_FD_grid

     integer                                                                   :: igrid, ispec, iglob, ier
     double precision,         dimension(NGLLX)                                :: hxis
     double precision,         dimension(NGLLY)                                :: hetas
     double precision,         dimension(NGLLZ)                                :: hgammas
     double precision                                                          :: val

     allocate(valence( projection_fd%nx,projection_fd%ny, projection_fd%nz ),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 153')
     allocate(valence_tmp( projection_fd%nx,projection_fd%ny, projection_fd%nz ),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 154')
     allocate(model_on_FD_grid_tmp( projection_fd%nx, projection_fd%ny, projection_fd%nz ),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 155')

     valence_tmp(:,:,:) = 0._CUSTOM_REAL
     model_on_FD_grid_tmp(:,:,:) = 0._CUSTOM_REAL

     do igrid = 1, projection_fd%nb_fd_point

        ispec= projection_fd%ispec_selected(igrid)
        hxis(:) =  projection_fd%hxi(:,igrid)
        hetas(:) =  projection_fd%heta(:,igrid)
        hgammas(:) =  projection_fd%hgamma(:,igrid)
        ifd =  projection_fd%index_on_fd_grid(1,igrid)
        jfd =  projection_fd%index_on_fd_grid(2,igrid)
        kfd =  projection_fd%index_on_fd_grid(3,igrid)
        val = 0.d0
        do kgg = 1,NGLLZ
           do jgg = 1,NGLLY
              do igg = 1,NGLLX
                 iglob = ibool(igg,jgg,kgg,ispec)
                 val = val + hxis(igg)*hetas(jgg)*hgammas(kgg)*dble(model_on_SEM_mesh(igg,jgg,kgg,ispec))
              enddo
           enddo
        enddo
        !write(*,*) ifd, projection_fd%nx
        model_on_FD_grid_tmp(ifd, jfd, kfd) = val
        valence_tmp(ifd, jfd, kfd) =  valence_tmp(ifd, jfd, kfd) + 1

!!$        if (DEBUG_MODE) then
!!$           write(IIDD, *) igrid, ifd, jfd, kfd, valence_tmp(ifd, jfd, kfd)
!!$        endif

     enddo

     model_on_FD_grid(:,:,:)=0.
     call sum_all_all_cr_array(model_on_FD_grid_tmp(1,1,1),  model_on_FD_grid(1,1,1), &
                               projection_fd%nx*projection_fd%ny*projection_fd%nz)
     valence(:,:,:) = 0.
     call sum_all_all_cr_array(valence_tmp(1,1,1),  valence(1,1,1), projection_fd%nx*projection_fd%ny*projection_fd%nz)

     !! check if valence is not == 0
     if (myrank == 0) then
        write(IIDD, *)
        write(IIDD, *) 'Check interp'
        write(IIDD, *)
        do kfd = 1, projection_fd%nz
           do jfd = 1, projection_fd%ny
              do ifd = 1, projection_fd%nx
                 if (valence(ifd, jfd, kfd) < 1.d0 ) then
                    !write(*,*) 'Projection failled : Warning point ', ifd, jfd, kfd, ' not found '
                    !write(*,*)  ox_fd_proj+(ifd-1)*hx_fd_proj, oy_fd_proj+(jfd-1)*hy_fd_proj, oz_fd_proj+(kfd-1)*hz_fd_proj
!!$                    write(IIDD,*)  ifd, jfd, kfd, valence(ifd, jfd, kfd)
                 else
                    model_on_FD_grid(ifd, jfd, kfd) = model_on_FD_grid(ifd, jfd, kfd) / valence(ifd, jfd, kfd)
                 endif
              enddo
           enddo
        enddo
     endif

     deallocate(valence, valence_tmp, model_on_FD_grid_tmp)

   end subroutine Project_model_SEM2FD_grid


!--------------------------------------------------------------------------------------------------------------------
!  locate MPI slice which contains the point and bcast to all
!--------------------------------------------------------------------------------------------------------------------

  subroutine locate_MPI_slice_and_bcast_to_all_1(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
                                                 xi, eta, gamma, ispec_selected, islice_selected, distance_from_target, myrank)

    integer,                                        intent(in)        :: myrank
    integer,                                        intent(inout)     :: ispec_selected, islice_selected
    double precision,                               intent(inout)     :: x_found,  y_found,  z_found
    double precision,                               intent(in)        :: x_to_locate, y_to_locate, z_to_locate
    double precision,                               intent(inout)     :: xi, eta, gamma, distance_from_target

    double precision,   dimension(:,:), allocatable                   :: distance_from_target_all
    double precision,   dimension(:,:), allocatable                   :: xi_all, eta_all, gamma_all
    double precision,   dimension(:),   allocatable                   :: x_array_found
    double precision,   dimension(:),   allocatable                   :: y_array_found
    double precision,   dimension(:),   allocatable                   :: z_array_found
    integer,            dimension(:,:), allocatable                   :: ispec_selected_all
    integer                                                           :: iproc, ier

    !! to avoid compler error when calling gather_all*
    double precision,  dimension(1)                                   :: distance_from_target_dummy
    double precision,  dimension(1)                                   :: xi_dummy, eta_dummy, gamma_dummy
    double precision,  dimension(1)                                   :: x_dummy, y_dummy, z_dummy
    integer,           dimension(1)                                   :: ispec_selected_dummy, islice_selected_dummy

    allocate(distance_from_target_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 156')
    allocate(xi_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 157')
    allocate(eta_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 158')
    allocate(gamma_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 159')
    allocate(x_array_found(0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 160')
    allocate(y_array_found(0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 161')
    allocate(z_array_found(0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 162')

    allocate(ispec_selected_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 163')

    distance_from_target = sqrt( (x_to_locate - x_found)**2&
                                +(y_to_locate - y_found)**2&
                                +(z_to_locate - z_found)**2)
!!$    if (DEBUG_MODE) then
!!$       write(IIDD, *) 'x :', x_to_locate,  x_found
!!$       write(IIDD, *) 'y :', y_to_locate,  y_found
!!$       write(IIDD, *) 'z :', z_to_locate,  z_found
!!$    endif

    !! this is just to avoid a compiler error
    distance_from_target_dummy(1)=distance_from_target
    x_dummy(1)=x_found
    y_dummy(1)=y_found
    z_dummy(1)=z_found
    xi_dummy(1)=xi
    eta_dummy(1)=eta
    gamma_dummy(1)=gamma
    ispec_selected_dummy(1)=ispec_selected

    ! gather all on myrank=0
    call gather_all_dp(distance_from_target_dummy, 1, distance_from_target_all, 1, NPROC)
    call gather_all_dp(xi_dummy,    1,  xi_all,    1,  NPROC)
    call gather_all_dp(eta_dummy,   1,  eta_all,   1,  NPROC)
    call gather_all_dp(gamma_dummy, 1,  gamma_all, 1,  NPROC)
    call gather_all_dp(x_dummy, 1,  x_array_found, 1,  NPROC)
    call gather_all_dp(y_dummy, 1,  y_array_found, 1,  NPROC)
    call gather_all_dp(z_dummy, 1,  z_array_found, 1,  NPROC)
    call gather_all_i(ispec_selected_dummy, 1, ispec_selected_all, 1, NPROC)

    ! find the slice and element to put the source
    if (myrank == 0) then

       distance_from_target = HUGEVAL

       do iproc=0, NPROC-1
          if (distance_from_target >= distance_from_target_all(1,iproc)) then
             distance_from_target =  distance_from_target_all(1,iproc)
             islice_selected_dummy(1) = iproc
             ispec_selected_dummy(1) = ispec_selected_all(1,iproc)
             xi_dummy(1)    = xi_all(1,iproc)
             eta_dummy(1)   = eta_all(1,iproc)
             gamma_dummy(1) = gamma_all(1,iproc)
             distance_from_target_dummy(1)=distance_from_target
             x_dummy(1) = x_array_found(iproc)
             y_dummy(1) = y_array_found(iproc)
             z_dummy(1) = z_array_found(iproc)
          endif
       enddo

    endif

    ! bcast form myrank=0
    call bcast_all_i(islice_selected_dummy,1)
    call bcast_all_i(ispec_selected_dummy,1)
    call bcast_all_dp(xi_dummy,1)
    call bcast_all_dp(eta_dummy,1)
    call bcast_all_dp(gamma_dummy,1)
    call bcast_all_dp(distance_from_target_dummy,1)
    call bcast_all_dp(distance_from_target_all,NPROC)
    call bcast_all_dp(x_dummy,1)
    call bcast_all_dp(y_dummy,1)
    call bcast_all_dp(z_dummy,1)
    !! it was just to avoid compler error
    islice_selected=islice_selected_dummy(1)
    ispec_selected=ispec_selected_dummy(1)
    xi=xi_dummy(1)
    eta=eta_dummy(1)
    gamma=gamma_dummy(1)
    distance_from_target=distance_from_target_dummy(1)
    x_found=x_dummy(1)
    y_found=y_dummy(1)
    z_found=z_dummy(1)

!!$    if (DEBUG_MODE .and. distance_from_target > 1.) then
!!$       write(IIDD, *) ' warning point no correctly localized  :'
!!$       write(IIDD, *) 'x :', x_to_locate,  x_found
!!$       write(IIDD, *) 'y :', y_to_locate,  y_found
!!$       write(IIDD, *) 'z :', z_to_locate,  z_found
!!$       write(IIDD, *) 'dist ', distance_from_target_all(1,islice_selected), islice_selected
!!$       write(IIDD, *)
!!$    endif

    deallocate(distance_from_target_all, xi_all, eta_all, gamma_all, x_array_found, y_array_found, z_array_found)
    deallocate(ispec_selected_all)

    !if (myrank==0)  write(INVERSE_LOG_FILE,*) ' bcast parameters to all slices  : passed'

  end subroutine locate_MPI_slice_and_bcast_to_all_1


!--------------------------------------------------------------------------------------------------------------------
!  locate point in mesh.
!--------------------------------------------------------------------------------------------------------------------
  subroutine locate_point_in_mesh_simple(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
                                         ispec_selected, xi_found, eta_found, gamma_found, x_found, y_found, z_found, myrank)

! note: this routine differs slightly from the one in locate_point.f90
!       by "simply" finding the best element using its inner GLL points

    double precision,                   intent(in)     :: x_to_locate, y_to_locate, z_to_locate
    real(kind=CUSTOM_REAL),             intent(in)     ::elemsize_max_glob
    integer,                            intent(in)     :: myrank
    integer,          dimension(NGNOD), intent(in)     :: iaddx,iaddy,iaddz
    double precision,                   intent(inout)  :: x_found,  y_found,  z_found
    double precision,                   intent(inout)  :: xi_found, eta_found, gamma_found
    integer,                            intent(inout)  :: ispec_selected

    ! locals
    integer                                            :: iter_loop , ispec, iglob, i, j, k
    double precision                                   :: x_target, y_target, z_target
    ! location search
    logical                                            :: located_target
    double precision                                   :: typical_size_squared, dist_squared
    double precision                                   :: distmin_squared
    double precision                                   :: x,y,z
    double precision                                   :: xi,eta,gamma,dx,dy,dz,dxi,deta
    double precision                                   :: xixs,xiys,xizs
    double precision                                   :: etaxs,etays,etazs
    double precision                                   :: gammaxs,gammays,gammazs, dgamma
    ! coordinates of the control points of the surface element
    double precision, dimension(NGNOD)                 :: xelm,yelm,zelm
    integer                                            :: ia,iax,iay,iaz
    integer                                            :: ix_initial_guess, iy_initial_guess, iz_initial_guess

    ! sets typical element size for search
    typical_size_squared =  elemsize_max_glob
    ! use 10 times the distance as a criterion for source detection
    typical_size_squared = (10. * typical_size_squared)**2

    ! INITIALIZE LOCATION --------
    x_target=x_to_locate
    y_target=y_to_locate
    z_target=z_to_locate
    ! flag to check that we located at least one target element
    located_target = .false.
    ispec_selected   = 1    !! first element by default
    ix_initial_guess = 1
    iy_initial_guess = 1
    iz_initial_guess = 1
    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    !! find the element candidate that may contain the target point
    do ispec = 1, NSPEC_AB

       iglob = ibool(MIDX,MIDY,MIDZ,ispec)
       dist_squared = (x_target- dble(xstore(iglob)))**2 &
            + (y_target - dble(ystore(iglob)))**2 &
            + (z_target - dble(zstore(iglob)))**2
       if (dist_squared > typical_size_squared) cycle ! exclude elements that are too far from target

       ! find closest GLL point form target
       do k=2, NGLLZ-1
          do j=2, NGLLY-1
             do i=2, NGLLX-1

                iglob=ibool(i,j,k,ispec)
                dist_squared = (x_target - dble(xstore(iglob)))**2 &
                     + (y_target - dble(ystore(iglob)))**2 &
                     + (z_target - dble(zstore(iglob)))**2

                if (dist_squared < distmin_squared) then

                   distmin_squared = dist_squared
                   ispec_selected  = ispec
                   ix_initial_guess = i
                   iy_initial_guess = j
                   iz_initial_guess = k
                   located_target = .true.

                   x_found = xstore(iglob)
                   y_found = ystore(iglob)
                   z_found = zstore(iglob)

                endif

             enddo
          enddo
       enddo

    enddo

    ! general coordinate of initial guess
    xi    = xigll(ix_initial_guess)
    eta   = yigll(iy_initial_guess)
    gamma = zigll(iz_initial_guess)

    ! define coordinates of the control points of the element
    do ia=1,NGNOD
       iax = 0
       iay = 0
       iaz = 0
       if (iaddx(ia) == 0) then
          iax = 1
       else if (iaddx(ia) == 1) then
          iax = (NGLLX+1)/2
       else if (iaddx(ia) == 2) then
          iax = NGLLX
       else
          call exit_MPI(myrank,'incorrect value of iaddx')
       endif

       if (iaddy(ia) == 0) then
          iay = 1
       else if (iaddy(ia) == 1) then
          iay = (NGLLY+1)/2
       else if (iaddy(ia) == 2) then
          iay = NGLLY
       else
          call exit_MPI(myrank,'incorrect value of iaddy')
       endif

       if (iaddz(ia) == 0) then
          iaz = 1
       else if (iaddz(ia) == 1) then
          iaz = (NGLLZ+1)/2
       else if (iaddz(ia) == 2) then
          iaz = NGLLZ
       else
          call exit_MPI(myrank,'incorrect value of iaddz')
       endif

       iglob = ibool(iax,iay,iaz,ispec_selected)
       xelm(ia) = dble(xstore(iglob))
       yelm(ia) = dble(ystore(iglob))
       zelm(ia) = dble(zstore(iglob))

    enddo

!!$    if (DEBUG_MODE) then
!!$       write(IIDD,*)
!!$       write(IIDD,'(a24,3f20.5)') 'Start to locate point :', x_target, y_target, z_target
!!$
!!$    endif

    ! iterate to solve the non linear system
    do iter_loop = 1, NUM_ITER

     ! recompute jacobian for the new point
       call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                               xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

!!$       if (iter_loop==1 .and. DEBUG_MODE) write(IIDD,'(a24,3f20.5)') 'Initilal guess        :', x, y, z

       ! compute distance to target location
       dx = - (x - x_target)
       dy = - (y - y_target)
       dz = - (z - z_target)

       ! compute increments
       dxi  = xixs*dx + xiys*dy + xizs*dz
       deta = etaxs*dx + etays*dy + etazs*dz
       dgamma = gammaxs*dx + gammays*dy + gammazs*dz

       ! update values
       xi = xi + dxi
       eta = eta + deta
       gamma = gamma + dgamma

       ! impose that we stay in that element
       ! (useful if user gives a point outside the mesh for instance)
       if (xi > 1.d0) xi     =  1.d0
       if (xi < -1.d0) xi     = -1.d0
       if (eta > 1.d0) eta    =  1.d0
       if (eta < -1.d0) eta    = -1.d0
       if (gamma > 1.d0) gamma  =  1.d0
       if (gamma < -1.d0) gamma  = -1.d0

!!$       if (DEBUG_MODE) write(IIDD, '(3f20.5)') dx, dy, dz

    enddo

    ! compute final coordinates of point found
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

    ! store xi,eta,gamma and x,y,z of point found
    ! note: xi/eta/gamma will be in range [-1,1]
    xi_found = xi
    eta_found = eta
    gamma_found = gamma

    x_found = x
    y_found = y
    z_found = z

  end subroutine locate_point_in_mesh_simple


!--------------------------------------------------------------------------------------------------------------------
!  locate point in element.
!--------------------------------------------------------------------------------------------------------------------
  subroutine locate_point_in_element(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
                                     ispec_selected, xi_found, eta_found, gamma_found, x_found, y_found, z_found, myrank, ispec)

    double precision,                   intent(in)     :: x_to_locate, y_to_locate, z_to_locate
    real(kind=CUSTOM_REAL),             intent(in)     :: elemsize_max_glob
    integer,                            intent(in)     :: myrank, ispec
    integer,          dimension(NGNOD), intent(in)     :: iaddx,iaddy,iaddz
    double precision,                   intent(inout)  :: x_found,  y_found,  z_found
    double precision,                   intent(inout)  :: xi_found, eta_found, gamma_found
    integer,                            intent(inout)  :: ispec_selected

    ! locals
    integer                                            :: iter_loop,  iglob, i, j, k
    double precision                                   :: x_target, y_target, z_target
    ! location search
    logical                                            :: located_target
    double precision                                   :: typical_size_squared, dist_squared
    double precision                                   :: distmin_squared
    double precision                                   :: x,y,z
    double precision                                   :: xi,eta,gamma,dx,dy,dz,dxi,deta
    double precision                                   :: xixs,xiys,xizs
    double precision                                   :: etaxs,etays,etazs
    double precision                                   :: gammaxs,gammays,gammazs, dgamma
    ! coordinates of the control points of the surface element
    double precision, dimension(NGNOD)                 :: xelm,yelm,zelm
    integer                                            :: ia,iax,iay,iaz
    integer                                            :: ix_initial_guess, iy_initial_guess, iz_initial_guess

    ! sets typical element size for search
    typical_size_squared =  elemsize_max_glob
    ! use 10 times the distance as a criterion for source detection
    typical_size_squared = (10. * typical_size_squared)**2

    ! INITIALIZE LOCATION --------
    x_target=x_to_locate
    y_target=y_to_locate
    z_target=z_to_locate
    ! flag to check that we located at least one target element
    located_target = .false.
    ispec_selected   = 1    !! first element by default
    ix_initial_guess = 1
    iy_initial_guess = 1
    iz_initial_guess = 1
    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    !! find the element candidate that may contain the target point
    !do ispec = 1, NSPEC_AB

       iglob = ibool(MIDX,MIDY,MIDZ,ispec)
       dist_squared = (x_target- dble(xstore(iglob)))**2 &
            + (y_target - dble(ystore(iglob)))**2 &
            + (z_target - dble(zstore(iglob)))**2
       !if (dist_squared > typical_size_squared) cycle ! exclude elements that are too far from target

       ! find closest GLL point form target
       do k=2, NGLLZ-1
          do j=2, NGLLY-1
             do i=2, NGLLX-1

                iglob=ibool(i,j,k,ispec)
                dist_squared = (x_target - dble(xstore(iglob)))**2 &
                     + (y_target - dble(ystore(iglob)))**2 &
                     + (z_target - dble(zstore(iglob)))**2

                if (dist_squared < distmin_squared) then

                   distmin_squared = dist_squared
                   ispec_selected  = ispec
                   ix_initial_guess = i
                   iy_initial_guess = j
                   iz_initial_guess = k
                   located_target = .true.

                   x_found = xstore(iglob)
                   y_found = ystore(iglob)
                   z_found = zstore(iglob)

                endif

             enddo
          enddo
       enddo

    !enddo

    ! general coordinate of initial guess
    xi    = xigll(ix_initial_guess)
    eta   = yigll(iy_initial_guess)
    gamma = zigll(iz_initial_guess)

    ! define coordinates of the control points of the element
    do ia=1,NGNOD
       iax = 0
       iay = 0
       iaz = 0
       if (iaddx(ia) == 0) then
          iax = 1
       else if (iaddx(ia) == 1) then
          iax = (NGLLX+1)/2
       else if (iaddx(ia) == 2) then
          iax = NGLLX
       else
          call exit_MPI(myrank,'incorrect value of iaddx')
       endif

       if (iaddy(ia) == 0) then
          iay = 1
       else if (iaddy(ia) == 1) then
          iay = (NGLLY+1)/2
       else if (iaddy(ia) == 2) then
          iay = NGLLY
       else
          call exit_MPI(myrank,'incorrect value of iaddy')
       endif

       if (iaddz(ia) == 0) then
          iaz = 1
       else if (iaddz(ia) == 1) then
          iaz = (NGLLZ+1)/2
       else if (iaddz(ia) == 2) then
          iaz = NGLLZ
       else
          call exit_MPI(myrank,'incorrect value of iaddz')
       endif

       iglob = ibool(iax,iay,iaz,ispec_selected)
       xelm(ia) = dble(xstore(iglob))
       yelm(ia) = dble(ystore(iglob))
       zelm(ia) = dble(zstore(iglob))

    enddo

!!$    if (DEBUG_MODE) then
!!$       write(IIDD,*)
!!$       write(IIDD,'(a24,3f20.5)') 'Start to locate point :', x_target, y_target, z_target
!!$
!!$    endif

    ! iterate to solve the non linear system
    do iter_loop = 1, NUM_ITER

     ! recompute jacobian for the new point
       call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                               xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

!!$       if (iter_loop==1 .and. DEBUG_MODE) write(IIDD,'(a24,3f20.5)') 'Initilal guess        :', x, y, z

       ! compute distance to target location
       dx = - (x - x_target)
       dy = - (y - y_target)
       dz = - (z - z_target)

       ! compute increments
       dxi  = xixs*dx + xiys*dy + xizs*dz
       deta = etaxs*dx + etays*dy + etazs*dz
       dgamma = gammaxs*dx + gammays*dy + gammazs*dz

       ! update values
       xi = xi + dxi
       eta = eta + deta
       gamma = gamma + dgamma

       ! impose that we stay in that element
       ! (useful if user gives a point outside the mesh for instance)
       if (xi > 1.d0) xi     =  1.d0
       if (xi < -1.d0) xi     = -1.d0
       if (eta > 1.d0) eta    =  1.d0
       if (eta < -1.d0) eta    = -1.d0
       if (gamma > 1.d0) gamma  =  1.d0
       if (gamma < -1.d0) gamma  = -1.d0

!!$       if (DEBUG_MODE) write(IIDD, '(3f20.5)') dx, dy, dz

    enddo

    ! compute final coordinates of point found
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

    ! store xi,eta,gamma and x,y,z of point found
    ! note: xi/eta/gamma will be in range [-1,1]
    xi_found = xi
    eta_found = eta
    gamma_found = gamma

    x_found = x
    y_found = y
    z_found = z

  end subroutine locate_point_in_element


!--------------------------------------------------------------------------------------------------------------------
!  trilinear interpolation
!--------------------------------------------------------------------------------------------------------------------

  subroutine TrilinearInterp(Vinterp,  x_loc, y_loc, z_loc, v1, v2, v3, v4, v5, v6, v7, v8, lx, ly, lz)

  real(kind=CUSTOM_REAL), intent(inout) :: Vinterp
  real(kind=CUSTOM_REAL), intent(in)    :: x_loc, y_loc, z_loc
  real(kind=CUSTOM_REAL), intent(in)    :: v1, v2, v3, v4, v5, v6, v7, v8, lx, ly, lz
  real(kind=CUSTOM_REAL)                :: dx, dy, dz

  dx = x_loc / lx
  dy = y_loc / ly
  dz = z_loc / lz

  Vinterp = &
  v1 * (1._CUSTOM_REAL - dx) * (1._CUSTOM_REAL - dy) * (1._CUSTOM_REAL - dz) + &
  v2 *                   dx  * (1._CUSTOM_REAL - dy) * (1._CUSTOM_REAL - dz) + &
  v3 *                   dx  *                   dy  * (1._CUSTOM_REAL - dz) + &
  v4 * (1._CUSTOM_REAL - dx) *                   dy  * (1._CUSTOM_REAL - dz) + &
  v5 * (1._CUSTOM_REAL - dx) * (1._CUSTOM_REAL - dy) *                   dz  + &
  v6 *                   dx  * (1._CUSTOM_REAL - dy) *                   dz  + &
  v7 *                   dx  *                   dy  *                   dz  + &
  v8 * (1._CUSTOM_REAL - dx) *                   dy  *                   dz

  end subroutine TrilinearInterp

end module projection_on_FD_grid
