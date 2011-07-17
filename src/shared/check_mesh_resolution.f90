!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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


  subroutine check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    kappastore,mustore,rho_vp,rho_vs, &
                                    DT, model_speed_max,min_resolved_period )

! check the mesh, stability and resolved period
!
! returns: maximum velocity in model ( model_speed_max ), minimum_period_resolved

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore,mustore,rho_vp,rho_vs
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  double precision :: DT
  real(kind=CUSTOM_REAL) :: model_speed_max,min_resolved_period

  ! local parameters
  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax,vpmin_glob,vpmax_glob,vsmin_glob,vsmax_glob
  real(kind=CUSTOM_REAL) :: distance_min,distance_max,distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: cmax,cmax_glob,pmax,pmax_glob
  real(kind=CUSTOM_REAL) :: dt_suggested,dt_suggested_glob,avg_distance

  logical:: DT_PRESENT

  integer :: myrank
  integer :: NSPEC_AB_global_min,NSPEC_AB_global_max,NSPEC_AB_global_sum
  integer :: NGLOB_AB_global_min,NGLOB_AB_global_max,NGLOB_AB_global_sum
  integer :: ispec,sizeprocs

  !********************************************************************************

  ! empirical choice for distorted elements to estimate time step and period resolved:
  ! courant number for time step estimate
  real(kind=CUSTOM_REAL),parameter :: COURANT_SUGGESTED = 0.3
  ! number of points per minimum wavelength for minimum period estimate
  real(kind=CUSTOM_REAL),parameter :: NPTS_PER_WAVELENGTH = 5

  !********************************************************************************
  !real(kind=CUSTOM_REAL),parameter :: NELEM_PER_WAVELENGTH = 1.5


  logical :: has_vs_zero

  ! initializations
  if( DT <= 0.0d0) then
    DT_PRESENT = .false.
  else
    DT_PRESENT = .true.
  endif

  vpmin_glob = HUGEVAL
  vpmax_glob = -HUGEVAL

  vsmin_glob = HUGEVAL
  vsmax_glob = -HUGEVAL

  distance_min_glob = HUGEVAL
  distance_max_glob = -HUGEVAL

  elemsize_min_glob = HUGEVAL
  elemsize_max_glob = -HUGEVAL

  cmax_glob = -HUGEVAL
  pmax_glob = -HUGEVAL

  dt_suggested_glob = HUGEVAL

  has_vs_zero = .false.

  ! checks courant number & minimum resolved period for each grid cell
  do ispec=1,NSPEC_AB

    ! determines minimum/maximum velocities within this element
    call get_vpvs_minmax(vpmin,vpmax,vsmin,vsmax,ispec,has_vs_zero, &
                        NSPEC_AB,kappastore,mustore,rho_vp,rho_vs)

    ! min/max for whole cpu partition
    vpmin_glob = min ( vpmin_glob, vpmin)
    vpmax_glob = max ( vpmax_glob, vpmax)

    vsmin_glob = min ( vsmin_glob, vsmin)
    vsmax_glob = max ( vsmax_glob, vsmax)

    ! computes minimum and maximum distance of neighbor GLL points in this grid cell
    call get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
                          NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    distance_min_glob = min( distance_min_glob, distance_min)
    distance_max_glob = max( distance_max_glob, distance_max)

    ! computes minimum and maximum size of this grid cell
    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                          NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    elemsize_min_glob = min( elemsize_min_glob, elemsize_min)
    elemsize_max_glob = max( elemsize_max_glob, elemsize_max)

    ! courant number
    ! based on minimum GLL point distance and maximum velocity
    ! i.e. on the maximum ratio of ( velocity / gridsize )
    if( DT_PRESENT ) then
      cmax = max( vpmax,vsmax ) * DT / distance_min
      cmax_glob = max(cmax_glob,cmax)
    endif

    ! suggested timestep
    dt_suggested = COURANT_SUGGESTED * distance_min / max( vpmax,vsmax )
    dt_suggested_glob = min( dt_suggested_glob, dt_suggested)

    ! estimation of minimum period resolved
    ! based on average GLL distance within element and minimum velocity
    !
    ! rule of thumb (Komatitsch et al. 2005):
    ! "average number of points per minimum wavelength in an element should be around 5."

    ! average distance between GLL points within this element
    avg_distance = elemsize_max / ( NGLLX - 1 )  ! since NGLLX = NGLLY = NGLLZ

    ! biggest possible minimum period such that number of points per minimum wavelength
    ! npts = ( min(vpmin,vsmin)  * pmax ) / avg_distance  is about ~ NPTS_PER_WAVELENGTH
    !
    ! note: obviously, this estimation depends on the choice of points per wavelength
    !          which is empirical at the moment.
    !          also, keep in mind that the minimum period is just an estimation and
    !          there is no such sharp cut-off period for valid synthetics.
    !          seismograms become just more and more inaccurate for periods shorter than this estimate.
    pmax = avg_distance / min( vpmin,vsmin ) * NPTS_PER_WAVELENGTH
    pmax_glob = max(pmax_glob,pmax)


    ! old: based on GLL distance, i.e. on maximum ratio ( gridspacing / velocity )
    !pmax = distance_max / min( vpmin,vsmin ) * NELEM_PER_WAVELENGTH
    !pmax_glob = max(pmax_glob,pmax)

  enddo

! determines global min/max values from all cpu partitions
  if( DT_PRESENT ) then
    ! courant number
    cmax = cmax_glob
    call max_all_cr(cmax,cmax_glob)
  endif

  ! minimum period
  pmax = pmax_glob
  call max_all_cr(pmax,pmax_glob)

  ! time step
  dt_suggested = dt_suggested_glob
  call min_all_cr(dt_suggested,dt_suggested_glob)

  ! Vp velocity
  vpmin = vpmin_glob
  vpmax = vpmax_glob
  call min_all_cr(vpmin,vpmin_glob)
  call max_all_cr(vpmax,vpmax_glob)

  ! Vs velocity
  vsmin = vsmin_glob
  if( has_vs_zero ) vsmin = 0.0

  vsmax = vsmax_glob
  call min_all_cr(vsmin,vsmin_glob)
  call max_all_cr(vsmax,vsmax_glob)

  ! checks velocities
  if( vpmin_glob <= 0.0_CUSTOM_REAL ) then
    call exit_mpi(myrank,"error: vp minimum velocity")
  endif
  if( vpmax_glob >= HUGEVAL ) then
    call exit_mpi(myrank,"error: vp maximum velocity")
  endif
  if( vsmin_glob < 0.0_CUSTOM_REAL ) then
    call exit_mpi(myrank,"error: vs minimum velocity")
  endif
  if( vsmax_glob >= HUGEVAL ) then
    call exit_mpi(myrank,"error: vs maximum velocity")
  endif

  ! GLL point distance
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  call min_all_cr(distance_min,distance_min_glob)
  call max_all_cr(distance_max,distance_max_glob)

  ! element size
  elemsize_min = elemsize_min_glob
  elemsize_max = elemsize_max_glob
  call min_all_cr(elemsize_min,elemsize_min_glob)
  call max_all_cr(elemsize_max,elemsize_max_glob)

  ! checks mesh
  if( distance_min_glob <= 0.0_CUSTOM_REAL ) then
    call exit_mpi(myrank,"error: GLL points minimum distance")
  endif
  if( distance_max_glob >= HUGEVAL ) then
    call exit_mpi(myrank,"error: GLL points maximum distance")
  endif
  if( elemsize_min_glob <= 0.0_CUSTOM_REAL ) then
    call exit_mpi(myrank,"error: element minimum size")
  endif
  if( elemsize_max_glob >= HUGEVAL ) then
    call exit_mpi(myrank,"error: element maximum size")
  endif

!! DK DK May 2009: added this to print the minimum and maximum number of elements
!! DK DK May 2009: and points in the CUBIT + SCOTCH mesh
  call min_all_i(NSPEC_AB,NSPEC_AB_global_min)
  call max_all_i(NSPEC_AB,NSPEC_AB_global_max)
  call sum_all_i(NSPEC_AB,NSPEC_AB_global_sum)

  call min_all_i(NGLOB_AB,NGLOB_AB_global_min)
  call max_all_i(NGLOB_AB,NGLOB_AB_global_max)
  call sum_all_i(NGLOB_AB,NGLOB_AB_global_sum)

  call world_size(sizeprocs)

! outputs infos
  if ( myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) '********'
    write(IMAIN,*) 'minimum and maximum number of elements'
    write(IMAIN,*) 'and points in the CUBIT + SCOTCH mesh:'
    write(IMAIN,*)
    write(IMAIN,*) 'NSPEC_AB_global_min = ',NSPEC_AB_global_min
    write(IMAIN,*) 'NSPEC_AB_global_max = ',NSPEC_AB_global_max
    !write(IMAIN,*) 'NSPEC_AB_global_mean = ',NSPEC_AB_global_sum / float(sizeprocs)
    write(IMAIN,*) 'NSPEC_AB_global_sum = ',NSPEC_AB_global_sum
    write(IMAIN,*)
    write(IMAIN,*) 'NGLOB_AB_global_min = ',NGLOB_AB_global_min
    write(IMAIN,*) 'NGLOB_AB_global_max = ',NGLOB_AB_global_max
    write(IMAIN,*) 'NGLOB_AB_global_sum = ',NGLOB_AB_global_sum
    write(IMAIN,*)
    write(IMAIN,*) '********'
    write(IMAIN,*) 'Model: P velocity min,max = ',vpmin_glob,vpmax_glob
    write(IMAIN,*) 'Model: S velocity min,max = ',vsmin_glob,vsmax_glob
    write(IMAIN,*) '********'
    write(IMAIN,*)
    write(IMAIN,*) '*********************************************'
    write(IMAIN,*) '*** Verification of simulation parameters ***'
    write(IMAIN,*) '*********************************************'
    write(IMAIN,*)
    write(IMAIN,*) '*** Max GLL point distance = ',distance_max_glob
    write(IMAIN,*) '*** Min GLL point distance = ',distance_min_glob
    write(IMAIN,*) '*** Max/min ratio = ',distance_max_glob/distance_min_glob
    write(IMAIN,*) '*** Max element size = ',elemsize_max_glob
    write(IMAIN,*) '*** Min element size = ',elemsize_min_glob
    write(IMAIN,*) '*** Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
    write(IMAIN,*)
    write(IMAIN,*) '*** Minimum period resolved = ',pmax_glob
    write(IMAIN,*) '*** Maximum suggested time step = ',dt_suggested_glob
    write(IMAIN,*)
    if( DT_PRESENT ) then
      write(IMAIN,*) '*** for DT : ',DT
      write(IMAIN,*) '*** Max stability for wave velocities = ',cmax_glob
      write(IMAIN,*)
    endif
  endif

  ! returns the maximum velocity
  if( myrank == 0 ) then
    if( vpmax_glob > vsmax_glob ) then
      model_speed_max = vpmax_glob
    else
      model_speed_max = vsmax_glob
    endif
  endif
  call bcast_all_cr(model_speed_max,1)

  ! returns minimum period
  if( myrank == 0 ) min_resolved_period = pmax_glob
  call bcast_all_cr(min_resolved_period,1)


  end subroutine check_mesh_resolution

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_vpvs_minmax(vpmin,vpmax,vsmin,vsmax,ispec,has_vs_zero, &
                            NSPEC_AB,kappastore,mustore,rho_vp,rho_vs)

! calculates the min/max size of the specified element (ispec)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax

  integer :: ispec
  logical :: has_vs_zero

  integer :: NSPEC_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
    kappastore,mustore,rho_vp,rho_vs

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ):: vp_elem,vs_elem
  real(kind=CUSTOM_REAL), dimension(1) :: val_min,val_max

  ! initializes
  vpmin = HUGEVAL
  vpmax = -HUGEVAL
  vsmin = HUGEVAL
  vsmax = -HUGEVAL

  ! vp
  where( rho_vp(:,:,:,ispec) > TINYVAL )
    vp_elem(:,:,:) = (FOUR_THIRDS * mustore(:,:,:,ispec) &
                    + kappastore(:,:,:,ispec)) / rho_vp(:,:,:,ispec)
  elsewhere
    vp_elem(:,:,:) = 0.0
  endwhere

  val_min = minval(vp_elem(:,:,:))
  val_max = maxval(vp_elem(:,:,:))

  vpmin = min(vpmin,val_min(1))
  vpmax = max(vpmax,val_max(1))

  ! vs
  where( rho_vs(:,:,:,ispec) > TINYVAL )
    vs_elem(:,:,:) = mustore(:,:,:,ispec) / rho_vs(:,:,:,ispec)
  elsewhere
    vs_elem(:,:,:) = 0.0
  endwhere

  val_min = minval(vs_elem(:,:,:))
  val_max = maxval(vs_elem(:,:,:))

  ! ignore fluid regions with Vs = 0
  if( val_min(1) > 0.0001 ) then
    vsmin = min(vsmin,val_min(1))
  else
    has_vs_zero = .true.
  endif
  vsmax = max(vsmax,val_max(1))

  end subroutine get_vpvs_minmax

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
                          NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

! calculates the min/max distances between neighboring GLL points within
! the specified element (ispec)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) :: distance_min,distance_max

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: dx,x0,y0,z0
  integer :: i,j,k,ii,jj,kk,iglob_a,iglob_b

  ! initializes
  distance_min = HUGEVAL
  distance_max = -HUGEVAL

  ! loops over all GLL points
  do k=1,NGLLZ-1
    do j=1,NGLLY-1
      do i=1,NGLLX-1
        iglob_a = ibool(i,j,k,ispec)
        x0 = xstore(iglob_a)
        y0 = ystore(iglob_a)
        z0 = zstore(iglob_a)

        ! loops over nearest neighbor points
        ! maybe a faster method could be found...
        do kk=k-1,k+1
          do jj=j-1,j+1
            do ii=i-1,i+1
              if( ii < 1 .or. jj < 1 .or. kk < 1 ) cycle
              ! distances between points
              iglob_b = ibool(ii,jj,kk,ispec)
              if( iglob_a /= iglob_b) then
                dx = sqrt( ( x0 - xstore(iglob_b) )**2 &
                          + ( y0 - ystore(iglob_b) )**2 &
                          + ( z0 - zstore(iglob_b) )**2 )
                if( dx < distance_min) distance_min = dx
                if( dx > distance_max) distance_max = dx
              endif
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo

  end subroutine get_GLL_minmaxdistance


!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                          NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

! calculates the min/max size of the specified element (ispec)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: dx,x0,y0,z0
  integer :: i,j,k,icorner,jcorner,iglob_a,iglob_b

  ! corners indices of reference cube faces
  ! shapes of arrays below
  integer,dimension(2),parameter :: corner_shape = (/3,NGNOD/)
  integer,dimension(3,NGNOD),parameter :: corner_ijk = &
    reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ, &
      NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ /),corner_shape)

  ! initializes
  elemsize_min = HUGEVAL
  elemsize_max = -HUGEVAL

  ! loops over corners
  do icorner=1,NGNOD
    i = corner_ijk(1,icorner)
    j = corner_ijk(2,icorner)
    k = corner_ijk(3,icorner)

    iglob_a = ibool(i,j,k,ispec)
    x0 = xstore(iglob_a)
    y0 = ystore(iglob_a)
    z0 = zstore(iglob_a)

    ! loops over all other corners
    do jcorner = icorner+1,NGNOD
      i = corner_ijk(1,jcorner)
      j = corner_ijk(2,jcorner)
      k = corner_ijk(3,jcorner)

      ! coordinates
      iglob_b = ibool(i,j,k,ispec)

      ! distances between points
      if( iglob_a /= iglob_b) then
        dx = sqrt( ( x0 - xstore(iglob_b) )**2 &
                  + ( y0 - ystore(iglob_b) )**2 &
                  + ( z0 - zstore(iglob_b) )**2 )
        if( dx < elemsize_min) elemsize_min = dx
        if( dx > elemsize_max) elemsize_max = dx
      endif

    enddo
  enddo

  end subroutine get_elem_minmaxsize
