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

module mesh_tools

  use specfem_par, only: CUSTOM_REAL, HUGEVAL, TINYVAL, NGNOD, NUM_ITER, NPROC, MAX_STRING_LEN, &
                         NGLLX, NGLLY, NGLLZ, NDIM, NSPEC_AB, NGLOB_AB, MIDX, MIDY, MIDZ, &
                         LOCAL_PATH, xigll, yigll, zigll, &
                         ibool, xstore, ystore, zstore

  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  !!----------------------------------------------------------------------------------------

  use inverse_problem_par

  implicit none

  PUBLIC  ::  compute_source_coeff,create_mass_matrices_Stacey_duplication_routine,  compute_force_elastic_arrays_source

!!$  PRIVATE ::  locate_MPI_slice_and_bcast_to_all

contains


!--------------------------------------------------------------------------------------------------------------------
! locate sources
!--------------------------------------------------------------------------------------------------------------------
!!$  subroutine locate_source(acqui_simu, myrank)
!!$
!!$    type(acqui), allocatable, dimension(:), intent(inout)  :: acqui_simu
!!$    integer,                                intent(in)     :: myrank
!!$
!!$    integer                                                :: ievent, nsrc_local, NEVENT
!!$    integer                                                :: ispec_selected_source, islice_selected_source
!!$    double precision                                       :: xi_source, eta_source, gamma_source
!!$    double precision                                       :: x_found,  y_found,  z_found
!!$    double precision                                       :: x_to_locate, y_to_locate, z_to_locate
!!$    real(kind=CUSTOM_REAL)                                 :: distance_min_glob,distance_max_glob
!!$    real(kind=CUSTOM_REAL)                                 :: elemsize_min_glob,elemsize_max_glob
!!$    real(kind=CUSTOM_REAL)                                 :: x_min_glob,x_max_glob
!!$    real(kind=CUSTOM_REAL)                                 :: y_min_glob,y_max_glob
!!$    real(kind=CUSTOM_REAL)                                 :: z_min_glob,z_max_glob
!!$    integer,                 dimension(NGNOD)              :: iaddx,iaddy,iaddz
!!$    double precision                                       :: distance_from_target
!!$
!!$200 format('         SOURCE : ', i5, '  LOCATED IN SLICE : ', i5, '  ELEMENT :', i10, '  ERROR IN LOCATION :', e15.6)
!!$300 format('         REAL POSITION : ', 3f20.5, '  FOUND POSITION : ', 3f20.5)
!!$
!!$    if (myrank == 0) then
!!$       write(INVERSE_LOG_FILE,*)
!!$       write(INVERSE_LOG_FILE,*) ' ... locate sources in specfem mesh :'
!!$    endif
!!$
!!$    NEVENT=acqui_simu(1)%nevent_tot
!!$
!!$    ! get mesh properties (mandatory before calling locate_point_in_mesh)
!!$    call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
!!$    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
!!$         x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
!!$         elemsize_min_glob,elemsize_max_glob, &
!!$         distance_min_glob,distance_max_glob)
!!$
!!$    do ievent = 1, NEVENT
!!$
!!$       x_to_locate = acqui_simu(ievent)%Xs
!!$       y_to_locate = acqui_simu(ievent)%Ys
!!$       z_to_locate = acqui_simu(ievent)%Zs
!!$
!!$       call locate_point_in_mesh(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
!!$            ispec_selected_source, xi_source, eta_source, gamma_source, x_found, y_found, z_found, myrank)
!!$
!!$       call locate_MPI_slice_and_bcast_to_all(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
!!$            xi_source, eta_source, gamma_source, ispec_selected_source, islice_selected_source, distance_from_target, myrank)
!!$
!!$       if (myrank == 0) then
!!$          write(INVERSE_LOG_FILE, 200) ievent, islice_selected_source, ispec_selected_source, distance_from_target
!!$          if (DEBUG_MODE) write(INVERSE_LOG_FILE, 300) x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found
!!$       endif
!!$
!!$       ! store in structure acqui
!!$       acqui_simu(ievent)%islice_slected_source=islice_selected_source
!!$       acqui_simu(ievent)%ispec_selected_source=ispec_selected_source
!!$
!!$       allocate(acqui_simu(ievent)%sourcearray(NDIM,NGLLX,NGLLY,NGLLZ))
!!$       nsrc_local=0
!!$      ! compute source array
!!$       if (myrank == islice_selected_source) then
!!$          nsrc_local = nsrc_local + 1
!!$          call compute_source_coeff(xi_source, eta_source, gamma_source, &
!!$               acqui_simu(ievent)%ispec_selected_source, acqui_simu(ievent)%sourcearray, &
!!$               acqui_simu(ievent)%Mxx, acqui_simu(ievent)%Myy, acqui_simu(ievent)%Mzz, &
!!$               acqui_simu(ievent)%Mxy, acqui_simu(ievent)%Mxz, acqui_simu(ievent)%Myz, &
!!$               acqui_simu(ievent)%Fx, acqui_simu(ievent)%Fy, acqui_simu(ievent)%Fz, &
!!$               acqui_simu(ievent)%source_type)
!!$       endif
!!$       acqui_simu(ievent)%nsources_local=nsrc_local
!!$
!!$    enddo
!!$
!!$    if (myrank == 0) then
!!$       write(INVERSE_LOG_FILE,*) ' ... locate sources passed'
!!$       write(INVERSE_LOG_FILE,*)
!!$    endif
!!$
!!$  end subroutine locate_source


!--------------------------------------------------------------------------------------------------------------------
!  locate MPI slice which contains the point and bcast to all  (TODO PUT IT ON OTHER MODULE)
!--------------------------------------------------------------------------------------------------------------------

  subroutine get_MPI_slice_and_bcast_to_all(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
                                            xi, eta, gamma, ispec_selected, islice_selected, distance_from_target, myrank)

    integer,                                        intent(in)        :: myrank
    integer,                                        intent(inout)     :: ispec_selected, islice_selected
    double precision,                               intent(in)        :: x_to_locate, y_to_locate, z_to_locate
    double precision,                               intent(inout)     :: x_found,  y_found,  z_found
    double precision,                               intent(inout)     :: xi, eta, gamma, distance_from_target

    double precision,   dimension(:,:), allocatable                   :: distance_from_target_all
    double precision,   dimension(:,:), allocatable                   :: xi_all, eta_all, gamma_all
    double precision,   dimension(:,:), allocatable                   ::  x_found_all, y_found_all, z_found_all
    integer,            dimension(:,:), allocatable                   :: ispec_selected_all
    integer                                                           :: iproc, ier

    !! to avoid compler error when calling gather_all*
    double precision,  dimension(1)                                   :: distance_from_target_dummy
    double precision,  dimension(1)                                   :: xi_dummy, eta_dummy, gamma_dummy
    double precision,  dimension(1)                                   :: x_found_dummy, y_found_dummy, z_found_dummy
    integer,           dimension(1)                                   :: ispec_selected_dummy, islice_selected_dummy

    allocate(distance_from_target_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 322')
    allocate(xi_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 323')
    allocate(eta_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 324')
    allocate(gamma_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 325')
    allocate(x_found_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 326')
    allocate(y_found_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 327')
    allocate(z_found_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 328')

    allocate(ispec_selected_all(1,0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 329')

    distance_from_target = sqrt( (x_to_locate - x_found)**2 + (y_to_locate - y_found)**2 + (z_to_locate - z_found)**2)

    !! this is just to avoid a compiler error
    distance_from_target_dummy(1)=distance_from_target
    xi_dummy(1)=xi
    eta_dummy(1)=eta
    gamma_dummy(1)=gamma
    ispec_selected_dummy(1)=ispec_selected
    x_found_dummy(1)=x_found
    y_found_dummy(1)=y_found
    z_found_dummy(1)=z_found

    ! gather all on myrank=0
    call gather_all_dp(distance_from_target_dummy, 1, distance_from_target_all, 1, NPROC)
    call gather_all_dp(xi_dummy,    1,  xi_all,    1,  NPROC)
    call gather_all_dp(eta_dummy,   1,  eta_all,   1,  NPROC)
    call gather_all_dp(gamma_dummy, 1,  gamma_all, 1,  NPROC)
    call gather_all_dp(x_found_dummy, 1,  x_found_all, 1,  NPROC)
    call gather_all_dp(y_found_dummy, 1,  y_found_all, 1,  NPROC)
    call gather_all_dp(z_found_dummy, 1,  z_found_all, 1,  NPROC)
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
             x_found_dummy(1)=x_found_all(1,iproc)
             y_found_dummy(1)=y_found_all(1,iproc)
             z_found_dummy(1)=z_found_all(1,iproc)
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
    call bcast_all_dp(x_found_dummy,1)
    call bcast_all_dp(y_found_dummy,1)
    call bcast_all_dp(z_found_dummy,1)

    !! it was just to avoid compler error
    islice_selected=islice_selected_dummy(1)
    ispec_selected=ispec_selected_dummy(1)
    xi=xi_dummy(1)
    eta=eta_dummy(1)
    gamma=gamma_dummy(1)
    x_found=x_found_dummy(1)
    y_found=y_found_dummy(1)
    z_found=z_found_dummy(1)
    distance_from_target=distance_from_target_dummy(1)

    deallocate(distance_from_target_all, xi_all, eta_all, gamma_all, x_found_all, y_found_all, z_found_all)
    deallocate(ispec_selected_all)

  end subroutine get_MPI_slice_and_bcast_to_all


!--------------------------------------------------------------------------------------------------------------------
!  locate point in mesh.
!--------------------------------------------------------------------------------------------------------------------

  subroutine get_point_in_mesh(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
                               ispec_selected, xi_found, eta_found, gamma_found, x_found, y_found, z_found, myrank)

    double precision,                   intent(in)     :: x_to_locate, y_to_locate, z_to_locate
    real(kind=CUSTOM_REAL),             intent(in)     :: elemsize_max_glob
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

    ! iterate to solve the non linear system
    do iter_loop = 1, NUM_ITER

     ! recompute jacobian for the new point
       call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                               xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

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

  end subroutine get_point_in_mesh


!---------------------------------------------------------------
! compute lagrange interpolation for point
!---------------------------------------------------------------

  subroutine  compute_source_coeff(xi,eta,gamma,ispec,interparray,Mxx,Myy,Mzz,Mxy,Mxz,Myz,factor_source,Fx,Fy,Fz,type,nu_source)

    real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out)    :: interparray
    double precision,                                         intent(in)     :: xi, eta, gamma
    double precision,                                         intent(in)     :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
    double precision,                                         intent(in)     :: factor_source
    double precision,                                         intent(in)     :: Fx,Fy,Fz
    double precision,       dimension(NDIM,NDIM)                             :: nu_source
    integer,                                                  intent(in)     :: ispec
    character(len=10),                                        intent(in)     :: type

    ! local
    double precision,       dimension(NGLLX)                                 :: hxis,hpxis
    double precision,       dimension(NGLLY)                                 :: hetas,hpetas
    double precision,       dimension(NGLLZ)                                 :: hgammas,hpgammas
    double precision                                                         :: norm,comp_x,comp_y,comp_z
    real(kind=CUSTOM_REAL)                                                   :: factor
    ! compute Lagrange polynomials at the source location
    call lagrange_any(xi,NGLLX,xigll,hxis,hpxis)
    call lagrange_any(eta,NGLLY,yigll,hetas,hpetas)
    call lagrange_any(gamma,NGLLZ,zigll,hgammas,hpgammas)

    select case (trim(adjustl(type)))

    case ('moment', 'shot')
       if (ispec_is_elastic(ispec)) then

          factor = 1.0
          call compute_arrays_source_cmt(ispec,interparray, &
               hxis,hetas,hgammas,hpxis,hpetas,hpgammas, &
               Mxx,Myy,Mzz,Mxy,Mxz,Myz)

       else if (ispec_is_acoustic(ispec)) then
          ! scalar moment of moment tensor values read in from CMTSOLUTION
          ! note: M0 by Dahlen and Tromp, eq. 5.91
          factor = 1.0/sqrt(2.0) * sqrt( Mxx**2 + Myy**2 + Mzz**2 &
               + 2*( Myz**2 + Mxz**2 + Mxy**2) )

          ! scales source such that it would be equivalent to explosion source moment tensor,
          ! where Mxx=Myy=Mzz, others Mxy,.. = zero, in equivalent elastic media
          ! (and getting rid of 1/sqrt(2) factor from scalar moment tensor definition above)
          factor = factor * sqrt(2.0) / sqrt(3.0)
          call compute_arrays_source_forcesolution(interparray,hxis,hetas,hgammas,factor,1.0d0,1.0d0,1.0d0,nu_source)

       else
          write(*,*) ' ABORT INVERSION: POINT SOURCE IS NOT IN ELASTIC OR ACOUSTIC DOMAIN'
          stop
       endif

    case('force')
       if (ispec_is_elastic(ispec)) then

         ! length of component vector
         norm = dsqrt( Fx**2 + Fy**2 + Fz**2 )
         ! checks norm of component vector
         if (norm < TINYVAL) stop 'error force point source: component vector has (almost) zero norm'

         ! normalizes given component vector
         comp_x = Fx/norm
         comp_y = Fy/norm
         comp_z = Fz/norm
       else if (ispec_is_acoustic(ispec)) then
         comp_x = 1.0d0
         comp_y = 1.0d0
         comp_z = 1.0d0
       else
          write(*,*) ' ABORT INVERSION: POINT SOURCE IS NOT IN ELASTIC OR ACOUSTIC DOMAIN'
          stop
       endif
       factor = sngl (factor_source)
       call compute_arrays_source_forcesolution(interparray,hxis,hetas,hgammas,factor,comp_x,comp_y,comp_z,nu_source)

    case default
       write(*,*) " Your source is not implemented in compute_source_coeff subroutine "
       stop

    end select

  end subroutine compute_source_coeff


!---------------------------------------------------------------
! compute mass matrix
!---------------------------------------------------------------

  subroutine create_mass_matrices_Stacey_duplication_routine()

    use constants, only: NGLLSQUARE
    use shared_parameters, only: NPROC, DT, STACEY_ABSORBING_CONDITIONS, PML_CONDITIONS, ELASTIC_SIMULATION,ACOUSTIC_SIMULATION

    use specfem_par, only: kappastore, rhostore, ibool, wxgll, wygll, wzgll, jacobianstore, &
         NGLLX, NGLLY, NGLLZ, NSPEC_AB, NGLOB_AB, num_abs_boundary_faces, abs_boundary_ispec, abs_boundary_ijk, &
         abs_boundary_normal, abs_boundary_jacobian2Dw, num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, my_neighbors_ext_mesh
    use specfem_par_elastic, only: rmass, rmassx, rmassy, rmassz, rho_vp, rho_vs, ispec_is_elastic
    use specfem_par_acoustic, only: rmass_acoustic, rmassz_acoustic, ispec_is_acoustic

    integer                           :: i, j, k, ispec, iglob, iface, igll, ier
    double precision                  :: weight
    real(kind=CUSTOM_REAL)            :: jacobianl, jacobianw
    real(kind=CUSTOM_REAL)            :: deltat, deltatover2
    real(kind=CUSTOM_REAL)            :: tx, ty, tz
    real(kind=CUSTOM_REAL)            :: nx, ny, nz, vn, sn

    ! acoustic domains
    if (ACOUSTIC_SIMULATION) then
       ! allocates memory
       if (.not. allocated(rmass_acoustic)) then
          allocate(rmass_acoustic(NGLOB_AB),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 330')
          if (ier /= 0) stop 'error allocating array rmass_acoustic'
       endif
       rmass_acoustic(:) = 0._CUSTOM_REAL

     ! returns acoustic mass matrix
       if (PML_CONDITIONS) then
          write(*,*) ' PML  not implemented yet '
          stop
          !call create_mass_matrices_pml_acoustic(nspec,ibool)
       else
          do ispec=1,nspec_ab
             if (ispec_is_acoustic(ispec)) then
                do k=1,NGLLZ
                   do j=1,NGLLY
                      do i=1,NGLLX
                         iglob = ibool(i,j,k,ispec)

                         weight = wxgll(i)*wygll(j)*wzgll(k)
                         jacobianl = jacobianstore(i,j,k,ispec)

                         ! distinguish between single and double precision for reals
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              real( dble(jacobianl) * weight / dble(kappastore(i,j,k,ispec)) ,kind=CUSTOM_REAL)
                      enddo
                   enddo
                enddo
             endif
          enddo
       endif
    endif

    if (ELASTIC_SIMULATION) then
       ! returns elastic mass matrix
       if (.not. allocated(rmass)) then
         allocate(rmass(NGLOB_AB),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 331')
       endif
       rmass(:) = 0._CUSTOM_REAL
       if (PML_CONDITIONS) then
          write(*,*) ' PML  not implemented yet '
          stop
          !call create_mass_matrices_pml_elastic(nspec,ibool)
       else
          do ispec=1,nspec_ab
             if (ispec_is_elastic(ispec)) then
                do k=1,NGLLZ
                   do j=1,NGLLY
                      do i=1,NGLLX
                         iglob = ibool(i,j,k,ispec)

                         weight = wxgll(i)*wygll(j)*wzgll(k)
                         jacobianl = jacobianstore(i,j,k,ispec)

                         rmass(iglob) = rmass(iglob) + &
                              real( dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)),kind=CUSTOM_REAL)
                      enddo
                   enddo
                enddo
             endif
          enddo
       endif
    endif

    if (STACEY_ABSORBING_CONDITIONS) then

       ! use the non-dimensional time step to make the mass matrix correction
       deltat = real(DT,kind=CUSTOM_REAL)
       deltatover2 = real(0.5d0*deltat,kind=CUSTOM_REAL)

        if (ELASTIC_SIMULATION) then
           rmassx(:) = 0._CUSTOM_REAL
           rmassy(:) = 0._CUSTOM_REAL
           rmassz(:) = 0._CUSTOM_REAL
        endif
        ! acoustic domains
        if (ACOUSTIC_SIMULATION) then
           if (.not. allocated(rmassz_acoustic)) then
             allocate(rmassz_acoustic(nglob_ab), stat=ier)
             if (ier /= 0) call exit_MPI_without_rank('error allocating array 332')
           endif
           rmassz_acoustic(:) = 0._CUSTOM_REAL
        endif

        do iface = 1,num_abs_boundary_faces
           ispec = abs_boundary_ispec(iface)
           if (ispec_is_elastic(ispec)) then
              ! reference GLL points on boundary face
              do igll = 1,NGLLSQUARE
                 ! gets local indices for GLL point
                 i = abs_boundary_ijk(1,igll,iface)
                 j = abs_boundary_ijk(2,igll,iface)
                 k = abs_boundary_ijk(3,igll,iface)

                 ! gets velocity
                 iglob = ibool(i,j,k,ispec)

                 ! gets associated normal
                 nx = abs_boundary_normal(1,igll,iface)
                 ny = abs_boundary_normal(2,igll,iface)
                 nz = abs_boundary_normal(3,igll,iface)

                 vn = deltatover2*(nx+ny+nz)

                 ! C*deltat/2 contributions
                 tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
                 ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
                 tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

                 ! gets associated, weighted jacobian
                 jacobianw = abs_boundary_jacobian2Dw(igll,iface)

                 ! assembles mass matrix on global points
                 rmassx(iglob) = rmassx(iglob) + tx*jacobianw
                 rmassy(iglob) = rmassy(iglob) + ty*jacobianw
                 rmassz(iglob) = rmassz(iglob) + tz*jacobianw
              enddo
           endif ! elastic

           ! acoustic element
           if (ispec_is_acoustic(ispec)) then

              ! reference GLL points on boundary face
              do igll = 1,NGLLSQUARE
                 ! gets local indices for GLL point
                 i = abs_boundary_ijk(1,igll,iface)
                 j = abs_boundary_ijk(2,igll,iface)
                 k = abs_boundary_ijk(3,igll,iface)

                 ! gets global index
                 iglob=ibool(i,j,k,ispec)

                 ! gets associated, weighted jacobian
                 jacobianw = abs_boundary_jacobian2Dw(igll,iface)

                 ! C * DT/2 contribution
                 !! need to divide by rho_vp
                 sn = deltatover2/ sqrt(kappastore(i,j,k,ispec) * rhostore(i,j,k,ispec))

                 rmassz_acoustic(iglob) = rmassz_acoustic(iglob) + jacobianw*sn
              enddo
           endif ! acoustic

        enddo

     endif

     !! need to create new masss matrix because the model is changed compared to
     !! the one stored in databases_mpi
     if (ELASTIC_SIMULATION) then
        if (STACEY_ABSORBING_CONDITIONS) then
           rmassx(:) = rmass(:) + rmassx(:)
           rmassy(:) = rmass(:) + rmassy(:)
           rmassz(:) = rmass(:) + rmassz(:)
        else
           rmassx(:) = rmass(:)
           rmassy(:) = rmass(:)
           rmassz(:) = rmass(:)
        endif
     endif

     if (ACOUSTIC_SIMULATION) then
        ! adds contributions
        if (STACEY_ABSORBING_CONDITIONS) then
           !if (USE_LDDRK) then
           !   rmass_acoustic(:) = rmass_acoustic(:)
           !else
              ! adds boundary contributions for newmark scheme
              rmass_acoustic(:) = rmass_acoustic(:) + rmassz_acoustic(:)
           !endif

        endif
     endif

     ! synchronize all the processes before assembling the mass matrix
     ! to make sure all the nodes have finished to read their databases
     call synchronize_all()

     if (ELASTIC_SIMULATION) then
        ! assemble mass matrix
        call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassx, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh)
        call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassy, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh)
        call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassz, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh)

        ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
        where(rmassx <= 0._CUSTOM_REAL) rmassx = 1._CUSTOM_REAL
        where(rmassy <= 0._CUSTOM_REAL) rmassy = 1._CUSTOM_REAL
        where(rmassz <= 0._CUSTOM_REAL) rmassz = 1._CUSTOM_REAL
        rmassx(:) = 1._CUSTOM_REAL / rmassx(:)
        rmassy(:) = 1._CUSTOM_REAL / rmassy(:)
        rmassz(:) = 1._CUSTOM_REAL / rmassz(:)

        ! not needed anymore
        deallocate(rmass)
     endif

     if (ACOUSTIC_SIMULATION) then
        call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_acoustic, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh)

        ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
        where(rmass_acoustic <= 0._CUSTOM_REAL) rmass_acoustic = 1._CUSTOM_REAL
        rmass_acoustic(:) = 1._CUSTOM_REAL / rmass_acoustic(:)
        ! not needed anymore
        if (STACEY_ABSORBING_CONDITIONS) deallocate(rmassz_acoustic)
     endif

   end subroutine create_mass_matrices_Stacey_duplication_routine

!---------------------------------------------------------------
! force source delta function in ispec in solid
!---------------------------------------------------------------

   subroutine compute_force_elastic_arrays_source(interparray, Fx, Fy, Fz, hxi, heta, hgamma)

     real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(inout)  :: interparray
     double precision,                                        intent(in)     :: Fx,Fy,Fz
     double precision,       dimension(NGLLX),                intent(in)     :: hxi
     double precision,       dimension(NGLLY),                intent(in)     :: heta
     double precision,       dimension(NGLLZ),                intent(in)     :: hgamma
     integer                                                                 :: i,j,k

     ! interpolates adjoint source onto GLL points within this element
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              interparray(1,i,j,k) = hxi(i) * heta(j) * hgamma(k) * Fx
              interparray(2,i,j,k) = hxi(i) * heta(j) * hgamma(k) * Fy
              interparray(3,i,j,k) = hxi(i) * heta(j) * hgamma(k) * Fz
           enddo
        enddo
     enddo

   end subroutine compute_force_elastic_arrays_source


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------
! trilinear interpolation with cheking
!---------------------------------------------------------------

   subroutine Get_value_by_trilinear_interp(interpolated_value,  x, y, z, regular_grid_array, nx, ny, nz, ox, oy, oz, hx, hy, hz)

      real(kind=CUSTOM_REAL),                                intent(inout)  :: interpolated_value
      real(kind=CUSTOM_REAL),                                intent(in)     :: x, y, z
      real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable, intent(in)     :: regular_grid_array
      real(kind=CUSTOM_REAL),                                intent(in)     :: ox, oy, oz, hx, hy, hz
      integer,                                               intent(in)     :: nx, ny, nz

      real(kind=CUSTOM_REAL)                                                :: x_loc, y_loc, z_loc
      real(kind=CUSTOM_REAL)                                                :: v1, v2, v3, v4, v5, v6, v7, v8
      integer                                                               :: i,j,k

      i = floor((x -  ox)/ hx) + 1
      j = floor((y -  oy)/ hy) + 1
      k = floor((z -  oz)/ hz) + 1

      if (i <= nx .and. j <= ny .and. k <= nz .and. &
           i > 0 .and. j > 0 .and. k > 0          ) then
         !! the point is insde regular_grid_array

         !! defaults values
         v1 = regular_grid_array(i,     j,     k    )
         v2 = v1; v3 = v1; v4 = v1; v5 = v1; v6 = v1; v7 = v1; v8 = v1

         if (i < nx) v2 = regular_grid_array(i + 1, j,     k    )
         if (i < nx .and. j < ny) v3 = regular_grid_array(i + 1, j + 1, k    )
         if (j < ny) v4 = regular_grid_array(i,     j + 1, k    )

         if (k < nz) then
            v5 = regular_grid_array(i,     j,     k + 1)
            if (i < nx) v6 = regular_grid_array(i + 1, j,     k + 1)
            if (i < nx .and. j < ny)   v7 = regular_grid_array(i + 1, j + 1, k + 1)
            if (j < ny) v8 = regular_grid_array(i,     j + 1, k + 1)
         endif

         x_loc = x - (ox + real( i - 1, CUSTOM_REAL) * hx)
         y_loc = y - (oy + real( j - 1, CUSTOM_REAL) * hy)
         z_loc = z - (oz + real( k - 1, CUSTOM_REAL) * hz)

         call TrilinrInterp(interpolated_value, x_loc, y_loc, z_loc, v1, v2, v3, v4, v5, v6, v7, v8, &
              hx, hy, hz)
      else
         write(*,*) " ERROR , point : ",x,y,z," is outside fd grid"
         stop
      endif

    end subroutine Get_value_by_trilinear_interp
!!----------------------------------------------------------------------------------------------------------
   subroutine TrilinrInterp(Vinterp,  x_loc, y_loc, z_loc, v1, v2, v3, v4, v5, v6, v7, v8, lx, ly, lz)

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

   end subroutine TrilinrInterp

end module mesh_tools
