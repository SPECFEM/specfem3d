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

! for acoustic solver

  subroutine compute_add_sources_acoustic(NSPEC_AB,NGLOB_AB,potential_dot_dot_acoustic, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  xi_source,eta_source,gamma_source, &
                                  hdur,hdur_gaussian,t_cmt,dt,t0, &
                                  sourcearrays,kappastore,ispec_is_acoustic,&
                                  SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                                  nrec,islice_selected_rec,ispec_selected_rec, &
                                  nadj_rec_local,adj_sourcearrays,b_potential_dot_dot_acoustic )

  use specfem_par,only: PRINT_SOURCE_TIME_FUNCTION,stf_used_total  
  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,t_cmt 
  double precision :: dt,t0
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays 

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr,&
   comp_source_time_function_gauss

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT
  integer:: nrec
  integer,dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer:: nadj_rec_local
  real(kind=CUSTOM_REAL),dimension(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ):: adj_sourcearrays
  real(kind=CUSTOM_REAL),dimension(NGLOB_ADJOINT):: b_potential_dot_dot_acoustic
  
! local parameters
  double precision :: f0
  double precision :: stf 
  real(kind=CUSTOM_REAL) stf_used,stf_used_total_all,time_source 
  integer :: isource,iglob,ispec,i,j,k
  integer :: irec_local,irec

! plotting source time function
  if(PRINT_SOURCE_TIME_FUNCTION .and. .not. phase_is_inner ) then
    ! initializes total
    stf_used_total = 0.0_CUSTOM_REAL
  endif

! forward simulations    
  if (SIMULATION_TYPE == 1) then

    ! adds acoustic sources
    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if(myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
    
          if( ispec_is_acoustic(ispec) ) then

            if(USE_FORCE_POINT_SOURCE) then

              ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
              iglob = ibool(nint(xi_source(isource)), &
                             nint(eta_source(isource)), &
                             nint(gamma_source(isource)), &
                             ispec)
               
              f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format
              t0 = 1.2d0/f0

              if (it == 1 .and. myrank == 0) then
                print *,'using a source of dominant frequency ',f0
                print *,'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
                print *,'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
              endif

              ! gaussian source time function
              !stf_used = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

              ! we use nu_source(:,3) here because we want a source normal to the surface.
              ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
              stf_used = 1.d10 * comp_source_time_function_rickr(dble(it-1)*DT-t0-t_cmt(isource),f0)              

              ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
              ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
              ! to add minus the source to Chi_dot_dot to get plus the source in pressure:
              
              ! acoustic source for pressure gets divided by kappa              
              ! source contribution
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                - stf_used / kappastore(nint(xi_source(isource)), &
                                                        nint(eta_source(isource)), &
                                                        nint(gamma_source(isource)),ispec)
               
            else   

              ! gaussian source time 
              stf = comp_source_time_function_gauss(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))
              
              ! quasi-heaviside  
              !stf = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))
              
              ! distinguishes between single and double precision for reals
              if(CUSTOM_REAL == SIZE_REAL) then
                stf_used = sngl(stf)
              else
                stf_used = stf
              endif

              ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
              ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
              ! to add minus the source to Chi_dot_dot to get plus the source in pressure

              !     add source array
              do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! adds source contribution
                      ! note: acoustic source for pressure gets divided by kappa
                      iglob = ibool(i,j,k,ispec)
                      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                              - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)                          
                   enddo
                enddo
              enddo

            endif ! USE_FORCE_POINT_SOURCE

            stf_used_total = stf_used_total + stf_used
            
          endif ! ispec_is_elastic
        endif ! ispec_is_inner     
      endif ! myrank
    
    enddo ! NSOURCES
  endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_potential.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields: 
!       time for b_potential..( it ) corresponds to (NSTEP - it - 1 )*DT - t0  ...
!       since we start with saved wavefields b_potential..( 0 ) = potential..( NSTEP ) which correspond
!       to a time (NSTEP - 1)*DT - t0 
!       (see sources for simulation_type 1 and seismograms)
!       now, at the beginning of the time loop, the numerical Newark time scheme updates
!       the wavefields, that is b_potential..( it=1) corresponds now to time (NSTEP -1 - 1)*DT - t0
!
! let's define the start time t  to (1-1)*DT - t0 = -t0, and the end time T to (NSTEP-1)*DT - t0
! these are the start and end times of all seismograms
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms, 
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for it=1: (NSTEP -1 - 1)*DT - t0 for backward wavefields corresponds to time T-1
!                    and time (T-1) corresponds now to index (NSTEP -1) in the adjoint source array

! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then

    if( it < NSTEP ) then    
      ! receivers act as sources    
      irec_local = 0
      do irec = 1,nrec
        ! add the source (only if this proc carries the source)
        if(myrank == islice_selected_rec(irec)) then
          irec_local = irec_local + 1
          ! adds source array
          ispec = ispec_selected_rec(irec)
          do k = 1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                
                ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid                
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                              - adj_sourcearrays(irec_local,NSTEP-it,1,i,j,k) / kappastore(i,j,k,ispec)
              enddo
            enddo
          enddo
        endif        
      enddo ! nrec    
    endif ! it
  endif

! adjoint simulations
  if (SIMULATION_TYPE == 3) then  
    ! adds acoustic sources
    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if(myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
    
          if( ispec_is_acoustic(ispec) ) then

            if(USE_FORCE_POINT_SOURCE) then

              ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
              iglob = ibool(nint(xi_source(isource)), &
                             nint(eta_source(isource)), &
                             nint(gamma_source(isource)), &
                             ispec)
               
              f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format
              t0 = 1.2d0/f0

              if (it == 1 .and. myrank == 0) then
                print *,'using a source of dominant frequency ',f0
                print *,'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
                print *,'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
              endif

              ! gaussian source time function
              !stf_used = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

              ! we use nu_source(:,3) here because we want a source normal to the surface.
              ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
              stf_used = 1.d10 * comp_source_time_function_rickr(dble(NSTEP-it-1)*DT-t0-t_cmt(isource),f0) 

              ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
              ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
              ! to add minus the source to Chi_dot_dot to get plus the source in pressure:
              
              ! acoustic source for pressure gets divided by kappa
              ! source contribution
              b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                          - stf_used / kappastore(nint(xi_source(isource)), &
                                               nint(eta_source(isource)), &
                                               nint(gamma_source(isource)),ispec) 
               
            else   

              ! gaussian source time 
              stf = comp_source_time_function_gauss(dble(NSTEP-it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

              ! distinguishes between single and double precision for reals
              if(CUSTOM_REAL == SIZE_REAL) then
                stf_used = sngl(stf)
              else
                stf_used = stf
              endif

              ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
              ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
              ! to add minus the source to Chi_dot_dot to get plus the source in pressure

              !     add source array
              do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! adds source contribution
                      ! note: acoustic source for pressure gets divided by kappa
                      iglob = ibool(i,j,k,ispec)
                      b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                              - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)                          
                   enddo
                enddo
              enddo

            endif ! USE_FORCE_POINT_SOURCE
            
            stf_used_total = stf_used_total + stf_used
            
          endif ! ispec_is_elastic
        endif ! ispec_is_inner     
      endif ! myrank
    
    enddo ! NSOURCES
  endif 

  ! master prints out source time function to file
  if(PRINT_SOURCE_TIME_FUNCTION .and. phase_is_inner) then
    time_source = (it-1)*DT - t0
    call sum_all_cr(stf_used_total,stf_used_total_all)
    if( myrank == 0 ) write(IOSTF,*) time_source,stf_used_total_all
  endif

  
end subroutine compute_add_sources_acoustic
