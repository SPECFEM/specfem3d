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

!----
!----  locate_source finds the correct position of the source
!----

  subroutine locate_source(ibool,NSOURCES,myrank,NSPEC_AB,NGLOB_AB,xstore,ystore,zstore, &
                 xigll,yigll,zigll,NPROC, &
                 t_cmt,yr,jda,ho,mi,utm_x_source,utm_y_source, &
                 DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                 islice_selected_source,ispec_selected_source, &
                 xi_source,eta_source,gamma_source, &
                 UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                 PRINT_SOURCE_TIME_FUNCTION, &
                 nu_source,iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                 ispec_is_acoustic,ispec_is_elastic, &
                 num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  implicit none

  include "constants.h"

  integer NPROC,UTM_PROJECTION_ZONE
  integer NSPEC_AB,NGLOB_AB,NSOURCES

  logical PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION

  double precision DT

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  integer myrank

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic,ispec_is_elastic

  integer yr,jda,ho,mi

  double precision t_cmt(NSOURCES)
  double precision sec

  integer iprocloop

  integer i,j,k,ispec,iglob,iglob_selected,inode,iface,isource,imin,imax,jmin,jmax,kmin,kmax,igll,jgll,kgll
  integer iselected,jselected,iface_selected,iadjust,jadjust
  integer iproc(1)

  double precision, dimension(NSOURCES) :: utm_x_source,utm_y_source
  double precision dist
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta

  ! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

  ! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

  ! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz
  double precision dgamma

  double precision final_distance_source(NSOURCES)

  double precision x_target_source,y_target_source,z_target_source

  double precision,dimension(1) :: altitude_source,distmin_ele
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all
  double precision,dimension(4) :: elevation_node,dist_node

  integer islice_selected_source(NSOURCES)

  ! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

  integer ispec_selected_source(NSOURCES)

  integer ngather, ns, ne, ig, is, ng

  integer, dimension(NGATHER_SOURCES,0:NPROC-1) :: ispec_selected_source_all
  double precision, dimension(NGATHER_SOURCES,0:NPROC-1) :: xi_source_all,eta_source_all,gamma_source_all, &
     final_distance_source_all,x_found_source_all,y_found_source_all,z_found_source_all
  double precision, dimension(3,3,NGATHER_SOURCES,0:NPROC-1) :: nu_source_all

  double precision, dimension(:), allocatable :: tmp_local
  double precision, dimension(:,:),allocatable :: tmp_all_local

  double precision hdur(NSOURCES) !, hdur_gaussian(NSOURCES) !, t0

  double precision, dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(3,3,NSOURCES) :: nu_source

  double precision, dimension(NSOURCES) :: lat,long,depth
  double precision moment_tensor(6,NSOURCES)

  character(len=256) OUTPUT_FILES

  double precision, dimension(NSOURCES) :: x_found_source,y_found_source,z_found_source
  double precision, dimension(NSOURCES) :: elevation
  double precision distmin

  integer, dimension(:), allocatable :: tmp_i_local
  integer, dimension(:,:),allocatable :: tmp_i_all_local

  ! for surface locating and normal computing with external mesh
  integer :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL), dimension(3) :: u_vector,v_vector,w_vector
  logical, dimension(NGLOB_AB) :: iglob_is_surface_external_mesh
  logical, dimension(NSPEC_AB) :: ispec_is_surface_external_mesh
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  integer ix_initial_guess_source,iy_initial_guess_source,iz_initial_guess_source

  ! for calculation of source time function
  !integer it
  !double precision time_source
  !double precision, external :: comp_source_time_function

  integer, dimension(NSOURCES) :: idomain
  integer, dimension(NGATHER_SOURCES,0:NPROC-1) :: idomain_all
  

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  ! read all the sources
  call get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,NSOURCES)

  ! checks half-durations
  do isource = 1, NSOURCES
    ! null half-duration indicates a Heaviside
    ! replace with very short error function
    if(hdur(isource) < 5. * DT) hdur(isource) = 5. * DT
  enddo
  
  ! convert the half duration for triangle STF to the one for gaussian STF
  !hdur_gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  !t0 = - 1.5d0 * minval(t_cmt-hdur)

  ! define topology of the control element
  call usual_hex_nodes(iaddx,iaddy,iaddz)

  ! get MPI starting time
  time_start = wtime()

  ! loop on all the sources
  do isource = 1,NSOURCES

    !
    ! r -> z, theta -> -y, phi -> x
    !
    !  Mrr =  Mzz
    !  Mtt =  Myy
    !  Mpp =  Mxx
    !  Mrt = -Myz
    !  Mrp =  Mxz
    !  Mtp = -Mxy

    ! get the moment tensor
    Mzz(isource) = + moment_tensor(1,isource)
    Mxx(isource) = + moment_tensor(3,isource)
    Myy(isource) = + moment_tensor(2,isource)
    Mxz(isource) = + moment_tensor(5,isource)
    Myz(isource) = - moment_tensor(4,isource)
    Mxy(isource) = - moment_tensor(6,isource)

    ! gets UTM x,y
    call utm_geo(long(isource),lat(isource),utm_x_source(isource),utm_y_source(isource), &
                   UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

    ! get approximate topography elevation at source long/lat coordinates
    ! set distance to huge initial value
    distmin = HUGEVAL
    if(num_free_surface_faces > 0) then
    iglob_selected = 1
    ! loop only on points inside the element
    ! exclude edges to ensure this point is not shared with other elements
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1
    do iface=1,num_free_surface_faces
          do j=jmin,jmax
             do i=imin,imax

                ispec = free_surface_ispec(iface)
                igll = free_surface_ijk(1,(j-1)*NGLLY+i,iface)
                jgll = free_surface_ijk(2,(j-1)*NGLLY+i,iface)
                kgll = free_surface_ijk(3,(j-1)*NGLLY+i,iface)
                iglob = ibool(igll,jgll,kgll,ispec)

                ! keep this point if it is closer to the receiver
                dist = dsqrt((utm_x_source(isource)-dble(xstore(iglob)))**2 + &
                     (utm_y_source(isource)-dble(ystore(iglob)))**2)
                if(dist < distmin) then
                   distmin = dist
                   iglob_selected = iglob
                   iface_selected = iface
                   iselected = i
                   jselected = j
                   altitude_source(1) = zstore(iglob_selected)
                endif
             enddo
          enddo
          ! end of loop on all the elements on the free surface
       end do
!  weighted mean at current point of topography elevation of the four closest nodes   
!  set distance to huge initial value
       distmin = HUGEVAL
       do j=jselected,jselected+1
          do i=iselected,iselected+1
             inode = 1
             do jadjust=0,1
                do iadjust= 0,1
                   ispec = free_surface_ispec(iface_selected)
                   igll = free_surface_ijk(1,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
                   jgll = free_surface_ijk(2,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
                   kgll = free_surface_ijk(3,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
                   iglob = ibool(igll,jgll,kgll,ispec)

                   elevation_node(inode) = zstore(iglob)
                   dist_node(inode) = dsqrt((utm_x_source(isource)-dble(xstore(iglob)))**2 + &
                        (utm_y_source(isource)-dble(ystore(iglob)))**2)
                   inode = inode + 1
                end do
             end do
             dist = sum(dist_node)
             if(dist < distmin) then
                distmin = dist
                altitude_source(1) = (dist_node(1)/dist)*elevation_node(1) + &
                     (dist_node(2)/dist)*elevation_node(2) + &
                     (dist_node(3)/dist)*elevation_node(3) + &
                     (dist_node(4)/dist)*elevation_node(4) 
             endif
          end do
       end do
    end if
    !  MPI communications to determine the best slice
    distmin_ele(1)= distmin
    call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
    call gather_all_dp(altitude_source,1,elevation_all,1,NPROC)
    if(myrank == 0) then
       iproc = minloc(distmin_ele_all)
       altitude_source(1) = elevation_all(iproc(1))         
    end if
    call bcast_all_dp(altitude_source,1)  
    elevation(isource) = altitude_source(1)

    ! orientation consistent with the UTM projection
    !     East
    nu_source(1,1,isource) = 1.d0
    nu_source(1,2,isource) = 0.d0
    nu_source(1,3,isource) = 0.d0
    !     North
    nu_source(2,1,isource) = 0.d0
    nu_source(2,2,isource) = 1.d0
    nu_source(2,3,isource) = 0.d0
    !     Vertical
    nu_source(3,1,isource) = 0.d0
    nu_source(3,2,isource) = 0.d0
    nu_source(3,3,isource) = 1.d0

    x_target_source = utm_x_source(isource)
    y_target_source = utm_y_source(isource)
    !z_target_source = depth(isource)
    z_target_source =  - depth(isource)*1000.0d0 + elevation(isource)

    ! set distance to huge initial value
    distmin = HUGEVAL

    ispec_selected_source(isource) = 0

    do ispec=1,NSPEC_AB


      ! define the interval in which we look for points
      if(USE_FORCE_POINT_SOURCE) then
        imin = 1
        imax = NGLLX

        jmin = 1
        jmax = NGLLY

        kmin = 1
        kmax = NGLLZ

      else
        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1

        kmin = 2
        kmax = NGLLZ - 1
      endif

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            iglob = ibool(i,j,k,ispec)

            if (.not. SOURCES_CAN_BE_BURIED_EXT_MESH) then
              if ((.not. iglob_is_surface_external_mesh(iglob)) .or. (.not. ispec_is_surface_external_mesh(ispec))) then
                cycle
              endif
            endif

            !       keep this point if it is closer to the source
            dist=dsqrt((x_target_source-dble(xstore(iglob)))**2 &
                  +(y_target_source-dble(ystore(iglob)))**2 &
                  +(z_target_source-dble(zstore(iglob)))**2)
            if(dist < distmin) then
              distmin=dist
              ispec_selected_source(isource)=ispec
              ix_initial_guess_source = i
              iy_initial_guess_source = j
              iz_initial_guess_source = k

              ! store xi,eta,gamma and x,y,z of point found
              xi_source(isource) = dble(ix_initial_guess_source)
              eta_source(isource) = dble(iy_initial_guess_source)
              gamma_source(isource) = dble(iz_initial_guess_source)
              x_found_source(isource) = xstore(iglob)
              y_found_source(isource) = ystore(iglob)
              z_found_source(isource) = zstore(iglob)

              ! compute final distance between asked and found (converted to km)
              final_distance_source(isource) = dsqrt((x_target_source-x_found_source(isource))**2 + &
                (y_target_source-y_found_source(isource))**2 + (z_target_source-z_found_source(isource))**2)

            endif

          enddo
        enddo
      enddo

    ! end of loop on all the elements in current slice
    enddo

    if (ispec_selected_source(isource) == 0) then
      final_distance_source(isource) = HUGEVAL
    endif

    ! sets whether acoustic (1) or elastic (2)
    if( ispec_is_acoustic( ispec_selected_source(isource) ) ) then
      idomain(isource) = 1
    else if( ispec_is_elastic( ispec_selected_source(isource) ) ) then
      idomain(isource) = 2
    else
      idomain(isource) = 0
    endif

    ! get normal to the face of the hexaedra if receiver is on the surface
    if ((.not. SOURCES_CAN_BE_BURIED_EXT_MESH) .and. &
       .not. (ispec_selected_source(isource) == 0)) then
      pt0_ix = -1
      pt0_iy = -1
      pt0_iz = -1
      pt1_ix = -1
      pt1_iy = -1
      pt1_iz = -1
      pt2_ix = -1
      pt2_iy = -1
      pt2_iz = -1
      ! we get two vectors of the face (three points) to compute the normal
      if (xi_source(isource) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(1,2,2,ispec_selected_source(isource)))) then
        pt0_ix = 1
        pt0_iy = NGLLY
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = 1
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif
      if (xi_source(isource) == NGLLX .and. &
         iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec_selected_source(isource)))) then
        pt0_ix = NGLLX
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = NGLLX
        pt1_iy = NGLLY
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = 1
        pt2_iz = NGLLZ
      endif
      if (eta_source(isource) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,1,2,ispec_selected_source(isource)))) then
        pt0_ix = 1
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = NGLLX
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = 1
        pt2_iy = 1
        pt2_iz = NGLLZ
      endif
      if (eta_source(isource) == NGLLY .and. &
         iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec_selected_source(isource)))) then
        pt0_ix = NGLLX
        pt0_iy = NGLLY
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = NGLLY
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif
      if (gamma_source(isource) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,2,1,ispec_selected_source(isource)))) then
        pt0_ix = NGLLX
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = NGLLY
        pt2_iz = 1
      endif
      if (gamma_source(isource) == NGLLZ .and. &
         iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec_selected_source(isource)))) then
        pt0_ix = 1
        pt0_iy = 1
        pt0_iz = NGLLZ
        pt1_ix = NGLLX
        pt1_iy = 1
        pt1_iz = NGLLZ
        pt2_ix = 1
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif

      if (pt0_ix<0 .or.pt0_iy<0 .or. pt0_iz<0 .or. &
         pt1_ix<0 .or. pt1_iy<0 .or. pt1_iz<0 .or. &
         pt2_ix<0 .or. pt2_iy<0 .or. pt2_iz<0) then
        stop 'error in computing normal for sources.'
      endif

      u_vector(1) = xstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_source(isource))) &
         - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      u_vector(2) = ystore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_source(isource))) &
         - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      u_vector(3) = zstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_source(isource))) &
         - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      v_vector(1) = xstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_source(isource))) &
         - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      v_vector(2) = ystore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_source(isource))) &
         - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))
      v_vector(3) = zstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_source(isource))) &
         - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_source(isource)))

      ! cross product
      w_vector(1) = u_vector(2)*v_vector(3) - u_vector(3)*v_vector(2)
      w_vector(2) = u_vector(3)*v_vector(1) - u_vector(1)*v_vector(3)
      w_vector(3) = u_vector(1)*v_vector(2) - u_vector(2)*v_vector(1)

      ! normalize vector w
      w_vector(:) = w_vector(:)/sqrt(w_vector(1)**2+w_vector(2)**2+w_vector(3)**2)

      ! build the two other vectors for a direct base: we normalize u, and v=w^u
      u_vector(:) = u_vector(:)/sqrt(u_vector(1)**2+u_vector(2)**2+u_vector(3)**2)
      v_vector(1) = w_vector(2)*u_vector(3) - w_vector(3)*u_vector(2)
      v_vector(2) = w_vector(3)*u_vector(1) - w_vector(1)*u_vector(3)
      v_vector(3) = w_vector(1)*u_vector(2) - w_vector(2)*u_vector(1)

      ! build rotation matrice nu for seismograms
      !     East (u)
      nu_source(1,1,isource) = u_vector(1)
      nu_source(1,2,isource) = v_vector(1)
      nu_source(1,3,isource) = w_vector(1)
      !     North (v)
      nu_source(2,1,isource) = u_vector(2)
      nu_source(2,2,isource) = v_vector(2)
      nu_source(2,3,isource) = w_vector(2)
      !     Vertical (w)
      nu_source(3,1,isource) = u_vector(3)
      nu_source(3,2,isource) = v_vector(3)
      nu_source(3,3,isource) = w_vector(3)

    endif ! of if (.not. RECEIVERS_CAN_BE_BURIED_EXT_MESH)

! *******************************************
! find the best (xi,eta,gamma) for the source
! *******************************************

    if(.not. USE_FORCE_POINT_SOURCE) then

      ! use initial guess in xi, eta and gamma
      xi = xigll(ix_initial_guess_source)
      eta = yigll(iy_initial_guess_source)
      gamma = zigll(iz_initial_guess_source)

      ! define coordinates of the control points of the element
      do ia=1,NGNOD

        if(iaddx(ia) == 0) then
          iax = 1
        else if(iaddx(ia) == 1) then
          iax = (NGLLX+1)/2
        else if(iaddx(ia) == 2) then
          iax = NGLLX
        else
          call exit_MPI(myrank,'incorrect value of iaddx')
        endif

        if(iaddy(ia) == 0) then
          iay = 1
        else if(iaddy(ia) == 1) then
          iay = (NGLLY+1)/2
        else if(iaddy(ia) == 2) then
          iay = NGLLY
        else
          call exit_MPI(myrank,'incorrect value of iaddy')
        endif

        if(iaddz(ia) == 0) then
          iaz = 1
        else if(iaddz(ia) == 1) then
          iaz = (NGLLZ+1)/2
        else if(iaddz(ia) == 2) then
          iaz = NGLLZ
        else
          call exit_MPI(myrank,'incorrect value of iaddz')
        endif

        iglob = ibool(iax,iay,iaz,ispec_selected_source(isource))
        xelm(ia) = dble(xstore(iglob))
        yelm(ia) = dble(ystore(iglob))
        zelm(ia) = dble(zstore(iglob))

      enddo

      ! iterate to solve the non linear system
      do iter_loop = 1,NUM_ITER

        ! recompute jacobian for the new point
        call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

        ! compute distance to target location
        dx = - (x - x_target_source)
        dy = - (y - y_target_source)
        dz = - (z - z_target_source)

        ! compute increments
        dxi  = xix*dx + xiy*dy + xiz*dz
        deta = etax*dx + etay*dy + etaz*dz
        dgamma = gammax*dx + gammay*dy + gammaz*dz

        ! update values
        xi = xi + dxi
        eta = eta + deta
        gamma = gamma + dgamma

        ! impose that we stay in that element
        ! (useful if user gives a source outside the mesh for instance)
        if (xi > 1.d0) xi = 1.d0
        if (xi < -1.d0) xi = -1.d0
        if (eta > 1.d0) eta = 1.d0
        if (eta < -1.d0) eta = -1.d0
        if (gamma > 1.d0) gamma = 1.d0
        if (gamma < -1.d0) gamma = -1.d0

      enddo

      ! compute final coordinates of point found
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

      ! store xi,eta,gamma and x,y,z of point found
      ! note: xi/eta/gamma will be in range [-1,1]
      xi_source(isource) = xi
      eta_source(isource) = eta
      gamma_source(isource) = gamma
      x_found_source(isource) = x
      y_found_source(isource) = y
      z_found_source(isource) = z

      ! compute final distance between asked and found (converted to km)
      final_distance_source(isource) = dsqrt((x_target_source-x_found_source(isource))**2 + &
        (y_target_source-y_found_source(isource))**2 + (z_target_source-z_found_source(isource))**2)

    endif ! of if (.not. USE_FORCE_POINT_SOURCE)

  ! end of loop on all the sources
  enddo

  ! now gather information from all the nodes
  ngather = NSOURCES/NGATHER_SOURCES
  if (mod(NSOURCES,NGATHER_SOURCES)/= 0) ngather = ngather+1
  do ig = 1, ngather
    ns = (ig-1) * NGATHER_SOURCES + 1
    ne = min(ig*NGATHER_SOURCES, NSOURCES)
    ng = ne - ns + 1

    ispec_selected_source_all(:,:) = -1

    ! avoids warnings about temporary creations of arrays for function call by compiler
    allocate(tmp_i_local(ng),tmp_i_all_local(ng,0:NPROC-1))
    tmp_i_local(:) = ispec_selected_source(ns:ne)    
    call gather_all_i(tmp_i_local,ng,tmp_i_all_local,ng,NPROC)
    ispec_selected_source_all(1:ng,:) = tmp_i_all_local(:,:)

    ! acoustic/elastic domain
    tmp_i_local(:) = idomain(ns:ne)    
    call gather_all_i(tmp_i_local,ng,tmp_i_all_local,ng,NPROC)
    idomain_all(1:ng,:) = tmp_i_all_local(:,:)

    deallocate(tmp_i_local,tmp_i_all_local)
    
    ! avoids warnings about temporary creations of arrays for function call by compiler
    allocate(tmp_local(ng),tmp_all_local(ng,0:NPROC-1))    
    tmp_local(:) = xi_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    xi_source_all(1:ng,:) = tmp_all_local(:,:)
        
    tmp_local(:) = eta_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    eta_source_all(1:ng,:) = tmp_all_local(:,:)
    
    tmp_local(:) = gamma_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    gamma_source_all(1:ng,:) = tmp_all_local(:,:)        
    
    tmp_local(:) = final_distance_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    final_distance_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = x_found_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    x_found_source_all(1:ng,:) = tmp_all_local(:,:)

    tmp_local(:) = y_found_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    y_found_source_all(1:ng,:) = tmp_all_local(:,:)
    
    tmp_local(:) = z_found_source(ns:ne)
    call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
    z_found_source_all(1:ng,:) = tmp_all_local(:,:)

    do i=1,3
      do j=1,3
        tmp_local(:) = nu_source(i,j,ns:ne)
        call gather_all_dp(tmp_local,ng,tmp_all_local,ng,NPROC)
        nu_source_all(i,j,1:ng,:) = tmp_all_local(:,:)
      enddo
    enddo
    deallocate(tmp_local,tmp_all_local)

    ! this is executed by main process only
    if(myrank == 0) then

      ! check that the gather operation went well
      if(any(ispec_selected_source_all(1:ng,:) == -1)) call exit_MPI(myrank,'gather operation failed for source')

      ! loop on all the sources
      do is = 1,ng
        isource = ns + is - 1

        ! loop on all the results to determine the best slice
        distmin = HUGEVAL
        do iprocloop = 0,NPROC-1
          if(final_distance_source_all(is,iprocloop) < distmin) then
            distmin = final_distance_source_all(is,iprocloop)
            islice_selected_source(isource) = iprocloop
            ispec_selected_source(isource) = ispec_selected_source_all(is,iprocloop)
            xi_source(isource) = xi_source_all(is,iprocloop)
            eta_source(isource) = eta_source_all(is,iprocloop)
            gamma_source(isource) = gamma_source_all(is,iprocloop)
            x_found_source(isource) = x_found_source_all(is,iprocloop)
            y_found_source(isource) = y_found_source_all(is,iprocloop)
            z_found_source(isource) = z_found_source_all(is,iprocloop)
            nu_source(:,:,isource) = nu_source_all(:,:,isource,iprocloop)
            idomain(isource) = idomain_all(is,iprocloop)
          endif
        enddo
        final_distance_source(isource) = distmin

      enddo
    endif !myrank
  enddo ! ngather

  if (myrank == 0) then

    do isource = 1,NSOURCES

      if(SHOW_DETAILS_LOCATE_SOURCE .or. NSOURCES == 1) then

        write(IMAIN,*)
        write(IMAIN,*) '*************************************'
        write(IMAIN,*) ' locating source ',isource
        write(IMAIN,*) '*************************************'
        write(IMAIN,*)
        write(IMAIN,*) 'source located in slice ',islice_selected_source(isource)
        write(IMAIN,*) '               in element ',ispec_selected_source(isource)
        if( idomain(isource) == 1 ) then
          write(IMAIN,*) '               in acoustic domain'
        else if( idomain(isource) == 2 ) then
          write(IMAIN,*) '               in elastic domain'
        else
          write(IMAIN,*) '               in unknown domain'        
        endif
        
        write(IMAIN,*)
        if(USE_FORCE_POINT_SOURCE) then
          write(IMAIN,*) '  xi coordinate of source in that element: ',nint(xi_source(isource))
          write(IMAIN,*) '  eta coordinate of source in that element: ',nint(eta_source(isource))
          write(IMAIN,*) '  gamma coordinate of source in that element: ',nint(gamma_source(isource))
          write(IMAIN,*) 'nu1 = ',nu_source(1,:,isource)
          write(IMAIN,*) 'nu2 = ',nu_source(2,:,isource)
          write(IMAIN,*) 'nu3 = ',nu_source(3,:,isource)
          write(IMAIN,*) 'at (x,y,z) coordinates = ',x_found_source(isource),y_found_source(isource),z_found_source(isource)
        else
          write(IMAIN,*) '   xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) 'gamma coordinate of source in that element: ',gamma_source(isource)
        endif

        ! add message if source is a Heaviside
        if(hdur(isource) < 5.*DT) then
          write(IMAIN,*)
          write(IMAIN,*) 'Source time function is a Heaviside, convolve later'
          write(IMAIN,*)
        endif

        write(IMAIN,*)
        write(IMAIN,*) ' half duration: ',hdur(isource),' seconds'
        write(IMAIN,*) '    time shift: ',t_cmt(isource),' seconds'

        write(IMAIN,*)
        write(IMAIN,*) 'original (requested) position of the source:'
        write(IMAIN,*)
        write(IMAIN,*) '          latitude: ',lat(isource)
        write(IMAIN,*) '         longitude: ',long(isource)
        write(IMAIN,*)
        if( SUPPRESS_UTM_PROJECTION ) then
          write(IMAIN,*) '             x: ',utm_x_source(isource)
          write(IMAIN,*) '             y: ',utm_y_source(isource)
        else
          write(IMAIN,*) '         UTM x: ',utm_x_source(isource)
          write(IMAIN,*) '         UTM y: ',utm_y_source(isource)        
        endif
        write(IMAIN,*) '         depth: ',depth(isource),' km'
        write(IMAIN,*) 'topo elevation: ',elevation(isource)

        write(IMAIN,*)
        write(IMAIN,*) 'position of the source that will be used:'
        write(IMAIN,*)
        if( SUPPRESS_UTM_PROJECTION ) then
          write(IMAIN,*) '             x: ',x_found_source(isource)
          write(IMAIN,*) '             y: ',y_found_source(isource)
        else
          write(IMAIN,*) '         UTM x: ',x_found_source(isource)
          write(IMAIN,*) '         UTM y: ',y_found_source(isource)        
        endif
        write(IMAIN,*) '         depth: ',dabs(z_found_source(isource) - elevation(isource))/1000.,' km'
        write(IMAIN,*) '             z: ',z_found_source(isource)
        write(IMAIN,*)

        ! display error in location estimate
        write(IMAIN,*) 'error in location of the source: ',sngl(final_distance_source(isource)),' m'

        ! add warning if estimate is poor
        ! (usually means source outside the mesh given by the user)
        if(final_distance_source(isource) > 3000.d0) then
          write(IMAIN,*)
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
        endif

      endif  ! end of detailed output to locate source

      if(PRINT_SOURCE_TIME_FUNCTION) then
        write(IMAIN,*)
        write(IMAIN,*) 'printing the source-time function'
      endif

      ! checks CMTSOLUTION format for acoustic case
      if( idomain(isource) == 1 ) then
        if( Mxx(isource) /= Myy(isource) .or. Myy(isource) /= Mzz(isource) .or. &
           Mxy(isource) > TINYVAL .or. Mxz(isource) > TINYVAL .or. Myz(isource) > TINYVAL ) then
            write(IMAIN,*)
            write(IMAIN,*) ' error CMTSOLUTION format for acoustic source:'
            write(IMAIN,*) '   acoustic source needs explosive moment tensor with'
            write(IMAIN,*) '      Mrr = Mtt = Mpp '
            write(IMAIN,*) '   and '
            write(IMAIN,*) '      Mrt = Mrp = Mtp = zero'
            write(IMAIN,*)
            call exit_mpi(myrank,'error acoustic source')
        endif
      endif

! end of loop on all the sources
    enddo

! display maximum error in location estimate
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(maxval(final_distance_source)),' m'
    write(IMAIN,*)

  endif     ! end of section executed by main process only

! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_source,NSOURCES)
  call bcast_all_i(ispec_selected_source,NSOURCES)
  call bcast_all_dp(xi_source,NSOURCES)
  call bcast_all_dp(eta_source,NSOURCES)
  call bcast_all_dp(gamma_source,NSOURCES)

! elapsed time since beginning of source detection
  if(myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for detection of sources in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
  endif

  end subroutine locate_source

