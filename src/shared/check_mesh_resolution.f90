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


  subroutine check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    kappastore,mustore,rho_vp,rho_vs, &
                                    DT, model_speed_max )

! check the mesh, stability and resolved period 
!
! returns: maximum velocity in model ( model_speed_max )
  
  implicit none
  
  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore,mustore,rho_vp,rho_vs
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  double precision :: DT
  real(kind=CUSTOM_REAL) :: model_speed_max
  
  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ):: vp_elem,vs_elem
  real(kind=CUSTOM_REAL), dimension(1) :: val_min,val_max  
  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax,vpmin_glob,vpmax_glob,vsmin_glob,vsmax_glob
  real(kind=CUSTOM_REAL) :: distance_min,distance_max,distance_min_glob,distance_max_glob,dx !,dy,dz
  real(kind=CUSTOM_REAL) :: cmax,cmax_glob,pmax,pmax_glob
  real(kind=CUSTOM_REAL) :: dt_suggested,dt_suggested_glob  
  
  logical:: DT_PRESENT
  
  integer :: myrank  
  integer :: NSPEC_AB_global_min,NSPEC_AB_global_max,NSPEC_AB_global_sum
  integer :: NGLOB_AB_global_min,NGLOB_AB_global_max,NGLOB_AB_global_sum    
  integer :: i,j,k,ii,jj,kk,ispec,iglob_a,iglob_b,sizeprocs

! estimation of time step and period resolved
  real(kind=CUSTOM_REAL),parameter :: COURANT_SUGGESTED = 0.3
  real(kind=CUSTOM_REAL),parameter :: NELEM_PER_WAVELENGTH = 1.5
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

  cmax_glob = -HUGEVAL  
  pmax_glob = -HUGEVAL

  dt_suggested_glob = HUGEVAL

  has_vs_zero = .false.

! checks courant number & minimum resolved period for each grid cell
  do ispec=1,NSPEC_AB
          
! determines minimum/maximum velocities within this element
    vpmin = HUGEVAL
    vpmax = -HUGEVAL
    vsmin = HUGEVAL
    vsmax = -HUGEVAL
    ! vp
    where( rho_vp(:,:,:,ispec) > TINYVAL )
      vp_elem(:,:,:) = (FOUR_THIRDS * mustore(:,:,:,ispec) + kappastore(:,:,:,ispec)) / rho_vp(:,:,:,ispec)
    elsewhere
      vp_elem(:,:,:) = 0.0
    endwhere
    ! vs    
    where( rho_vs(:,:,:,ispec) > TINYVAL )
      vs_elem(:,:,:) = mustore(:,:,:,ispec) / rho_vs(:,:,:,ispec)
    elsewhere
      vs_elem(:,:,:) = 0.0
    endwhere

    val_min = minval(vp_elem(:,:,:))
    val_max = maxval(vp_elem(:,:,:))

    vpmin = min(vpmin,val_min(1))
    vpmax = max(vpmax,val_max(1))

    val_min = minval(vs_elem(:,:,:))
    val_max = maxval(vs_elem(:,:,:))
    
    ! ignore fluid regions with Vs = 0
    if( val_min(1) > 0.0001 ) then
      vsmin = min(vsmin,val_min(1))
    else
      has_vs_zero = .true.
    endif
    vsmax = max(vsmax,val_max(1))

    ! min/max for whole cpu partition
    vpmin_glob = min ( vpmin_glob, vpmin)
    vpmax_glob = max ( vpmax_glob, vpmax)
    
    vsmin_glob = min ( vsmin_glob, vsmin)
    vsmax_glob = max ( vsmax_glob, vsmax)
    
! compute minimum and maximum distance of GLL points in this grid cell          
    distance_min = HUGEVAL
    distance_max = -HUGEVAL
    
    ! loops over all GLL points
    do k=1,NGLLZ-1
      do j=1,NGLLY-1
        do i=1,NGLLX-1
          iglob_a = ibool(i,j,k,ispec)

          ! loops over nearest neighbor points
          ! maybe a faster method could be found...
          do kk=k-1,k+1  
            do jj=j-1,j+1
              do ii=i-1,i+1
                if( ii < 1 .or. jj < 1 .or. kk < 1 ) cycle
                ! distances between points
                iglob_b = ibool(ii,jj,kk,ispec)
                if( iglob_a /= iglob_b) then
                  dx = sqrt( ( xstore(iglob_a) - xstore(iglob_b) )**2 &
                          + ( ystore(iglob_a) - ystore(iglob_b) )**2 &
                          + ( zstore(iglob_a) - zstore(iglob_b) )**2 )
                  if( dx < distance_min) distance_min = dx
                  if( dx > distance_max) distance_max = dx
                endif
              enddo
            enddo
          enddo
              
        enddo
      enddo
    enddo
    
    distance_min_glob = min( distance_min_glob, distance_min)
    distance_max_glob = max( distance_max_glob, distance_max)
    
    ! courant number
    if( DT_PRESENT ) then
      cmax = max( vpmax,vsmax ) * DT / distance_min    
      cmax_glob = max(cmax_glob,cmax)
    endif
    
    ! suggested timestep
    dt_suggested = COURANT_SUGGESTED * distance_min / max( vpmax,vsmax )
    dt_suggested_glob = min( dt_suggested_glob, dt_suggested)          
                            
    ! estimation of minimum period resolved
    pmax = distance_max / min( vpmin,vsmin ) * NELEM_PER_WAVELENGTH  
    pmax_glob = max(pmax_glob,pmax)
    
  enddo

! determines global min/max values from all cpu partitions
  if( DT_PRESENT ) then
    cmax = cmax_glob
    call max_all_cr(cmax,cmax_glob)
  endif
  
  pmax = pmax_glob
  call min_all_cr(pmax,pmax_glob)

  dt_suggested = dt_suggested_glob
  call min_all_cr(dt_suggested,dt_suggested_glob)

  vpmin = vpmin_glob
  vpmax = vpmax_glob
  call min_all_cr(vpmin,vpmin_glob)
  call max_all_cr(vpmax,vpmax_glob)

  vsmin = vsmin_glob
  if( has_vs_zero ) vsmin = 0.0
  
  vsmax = vsmax_glob
  call min_all_cr(vsmin,vsmin_glob)
  call max_all_cr(vsmax,vsmax_glob)
  
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  call min_all_cr(distance_min,distance_min_glob)
  call max_all_cr(distance_max,distance_max_glob)


! checks mesh
  if( distance_min_glob <= 0.0_CUSTOM_REAL ) then
    call exit_mpi(myrank,"error: GLL points minimum distance")
  endif
  if( distance_max_glob >= HUGEVAL ) then
    call exit_mpi(myrank,"error: GLL points maximum distance")
  endif
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
  
  
  end subroutine
  