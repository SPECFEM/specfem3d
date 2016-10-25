!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

! for poroelastic solver

  subroutine compute_add_sources_poroelastic(NSPEC_AB,NGLOB_AB, &
                                             accels,accelw, &
                                             rhoarraystore,phistore,tortstore, &
                                             ibool, &
                                             NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                                             sourcearrays, &
                                             ispec_is_poroelastic,SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                                             nrec,islice_selected_rec,ispec_selected_rec, &
                                             nadj_rec_local,adj_sourcearrays,b_accels,b_accelw, &
                                             NTSTEP_BETWEEN_READ_ADJSRC)

  use constants
  use specfem_par, only: xigll,yigll,zigll,xi_receiver,eta_receiver,gamma_receiver, &
                        station_name,network_name,adj_source_file, &
                        USE_FORCE_POINT_SOURCE, &
                        tshift_src,dt,t0, &
                        USE_LDDRK,istage, &
                        EXTERNAL_STF,user_source_time_function

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accels,accelw

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        phistore,tortstore
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays

  logical, dimension(NSPEC_AB) :: ispec_is_poroelastic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT
  integer:: nrec
  integer,dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer:: nadj_rec_local
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_ADJOINT) :: b_accels,b_accelw
  logical :: ibool_read_adj_arrays
  integer :: it_sub_adj,itime,NTSTEP_BETWEEN_READ_ADJSRC
  real(kind=CUSTOM_REAL),dimension(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ):: adj_sourcearrays

! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: adj_sourcearray
  real(kind=CUSTOM_REAL) stf_used

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_poroelastic

  integer :: isource,iglob,i,j,k,ispec,ier
  integer :: irec_local,irec
  real(kind=CUSTOM_REAL) :: phil,tortl,rhol_s,rhol_f,rhol_bar
  real(kind=CUSTOM_REAL) :: fac_s,fac_w

! forward simulations
  if (SIMULATION_TYPE == 1) then

    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_poroelastic(ispec)) then
          ! current time
          if (USE_LDDRK) then
            time_source_dble = dble(it-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(it-1)*DT - t0 - tshift_src(isource)
          endif

          ! determines source time function value
          stf = get_stf_poroelastic(time_source_dble,isource)

          !! VM VM add external source time function
          if (EXTERNAL_STF) then
            stf = user_source_time_function(it, isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          ! adds source array
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                ! get poroelastic parameters of current local GLL
                phil = phistore(i,j,k,ispec)
                tortl = tortstore(i,j,k,ispec)
                rhol_s = rhoarraystore(1,i,j,k,ispec)
                rhol_f = rhoarraystore(2,i,j,k,ispec)
                rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

                ! we distinguish between a single force which can be applied both in fluid and solid
                ! and a moment-tensor source which only makes sense for a solid
                if (USE_FORCE_POINT_SOURCE) then
                  ! single point force
                  ! the source is applied to both solid and fluid phase: bulk source.
                  fac_s = 1._CUSTOM_REAL - phil/tortl
                  fac_w = 1._CUSTOM_REAL - rhol_f/rhol_bar
                else
                  ! moment-tensor source can only be in solid phase
                  ! source in the solid phase only
                  fac_s = 1._CUSTOM_REAL
                  fac_w = - rhol_f/rhol_bar
                endif

                ! solid phase
                accels(:,iglob) = accels(:,iglob) &
                             + fac_s * sourcearrays(isource,:,i,j,k)*stf_used
                ! fluid phase
                accelw(:,iglob) = accelw(:,iglob) &
                             + fac_w * sourcearrays(isource,:,i,j,k)*stf_used
              enddo
            enddo
          enddo

        endif ! ispec_is_poroelastic
      endif ! myrank

    enddo ! NSOURCES
  endif ! forward

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_displ,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_displ( it ) corresponds to (NSTEP - it - 1)*DT - t0  ...
!       since we start with saved wavefields b_displ( 0) = displ( NSTEP ) which correspond
!       to a time (NSTEP - 1)*DT - t0
!       (see sources for simulation_type 1 and seismograms)
!       now, at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_displ( it=1) corresponds now to time (NSTEP -1 - 1)*DT - t0
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

    ! adds adjoint source in this partitions
    if (nadj_rec_local > 0) then

      ! add adjoint source following elastic block by block consideration
      ! read in adjoint sources block by block (for memory consideration)
      ! e.g., in exploration experiments, both the number of receivers (nrec) and
      ! the number of time steps (NSTEP) are huge,
      ! which may cause problems since we have a large array:
      !   adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this
      ! timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) ) !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries
      ! and interior --- indicated by 'iphase'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added
      ! twice
      if (ibool_read_adj_arrays) then

        ! allocates temporary source array
        allocate(adj_sourcearray(NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
        if (ier /= 0) stop 'error allocating array adj_sourcearray'

        !!! read ascii adjoint sources
        irec_local = 0
        do irec = 1, nrec
          ! compute source arrays
          if (myrank == islice_selected_rec(irec)) then
            irec_local = irec_local + 1
            ! reads in **net**.**sta**.**BH**.adj files
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
            call compute_arrays_adjoint_source(myrank,adj_source_file, &
                                               xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                                               adj_sourcearray, xigll,yigll,zigll, &
                                               it_sub_adj,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC)

            do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
              adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
            enddo

          endif
        enddo

        deallocate(adj_sourcearray)

      endif ! if (ibool_read_adj_arrays)

      if (it < NSTEP) then

        ! receivers act as sources
        irec_local = 0
        do irec = 1,nrec

          ! add the source (only if this proc carries the source)
          if (myrank == islice_selected_rec(irec)) then
            irec_local = irec_local + 1

            ispec = ispec_selected_rec(irec)
            if (ispec_is_poroelastic(ispec)) then
              ! adds source array
              do k = 1,NGLLZ
                do j = 1,NGLLY
                  do i = 1,NGLLX
                    iglob = ibool(i,j,k,ispec_selected_rec(irec))
                    ! get poroelastic parameters of current local GLL
                    phil = phistore(i,j,k,ispec_selected_rec(irec))
                    rhol_s = rhoarraystore(1,i,j,k,ispec_selected_rec(irec))
                    rhol_f = rhoarraystore(2,i,j,k,ispec_selected_rec(irec))
                    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

                    ! adjoint source is in the solid phase only since this is the only measurement
                    ! available

                    ! solid phase
                    accels(:,iglob) = accels(:,iglob)  &
                         + adj_sourcearrays(irec_local, &
                            NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC), &
                            :,i,j,k)
                    !
                    ! fluid phase
                    accelw(:,iglob) = accelw(:,iglob)  &
                         - rhol_f/rhol_bar * adj_sourcearrays(irec_local, &
                            NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC), &
                            :,i,j,k)
                  enddo
                enddo
              enddo
            endif ! ispec_is_poroelastic
          endif
        enddo ! nrec

      endif ! it
    endif ! nadj_rec_local
  endif !adjoint : if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then

! note:  b_displ() is read in after Newmark time scheme, thus
!           b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!           thus indexing is NSTEP - it , instead of NSTEP - it - 1

! adjoint simulations
  if (SIMULATION_TYPE == 3) then

    ! backward source reconstruction
    do isource = 1,NSOURCES

      ! add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_poroelastic(ispec)) then
          ! note: time step is now at NSTEP-it
          ! current time
          if (USE_LDDRK) then
            time_source_dble = dble(NSTEP-it)*DT - dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
          endif

          ! determines source time function value
          stf = get_stf_poroelastic(time_source_dble,isource)

          !! VM VM add external source time function
          if (EXTERNAL_STF) then
            ! time-reversed
            stf = user_source_time_function(NSTEP-it+1, isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          !  add source array
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                ! get poroelastic parameters of current local GLL
                phil = phistore(i,j,k,ispec)
                tortl = tortstore(i,j,k,ispec)
                rhol_s = rhoarraystore(1,i,j,k,ispec)
                rhol_f = rhoarraystore(2,i,j,k,ispec)
                rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f


                ! we distinguish between a single force which can be applied both in fluid and solid
                ! and a moment-tensor source which only makes sense for a solid
                if (USE_FORCE_POINT_SOURCE) then
                  ! single point force
                  ! the source is applied to both solid and fluid phase: bulk source.
                  fac_s = 1._CUSTOM_REAL - phil/tortl
                  fac_w = 1._CUSTOM_REAL - rhol_f/rhol_bar
                else
                  ! moment-tensor source can only be in solid phase
                  ! source in the solid phase only
                  fac_s = 1._CUSTOM_REAL
                  fac_w = - rhol_f/rhol_bar
                endif

                ! solid phase
                b_accels(:,iglob) = b_accels(:,iglob) &
                             + fac_s * sourcearrays(isource,:,i,j,k)*stf_used
                ! fluid phase
                b_accelw(:,iglob) = b_accelw(:,iglob) &
                             + fac_w * sourcearrays(isource,:,i,j,k)*stf_used
              enddo
            enddo
          enddo

        endif ! poroelastic
      endif ! myrank

    enddo ! NSOURCES
  endif ! adjoint

  end subroutine compute_add_sources_poroelastic

!
!=====================================================================
!

  double precision function get_stf_poroelastic(time_source_dble,isource)

! returns source time function value for specified time

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_RICKER_TIME_FUNCTION,hdur,hdur_Gaussian,DT

  implicit none

  double precision,intent(in) :: time_source_dble
  integer,intent(in) :: isource

  ! local parameters
  double precision :: stf

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr, &
    comp_source_time_function_gauss

  ! determines source time function value
  if (USE_FORCE_POINT_SOURCE) then
    ! single point force
    if (USE_RICKER_TIME_FUNCTION) then
      ! Ricker
      ! f0 has been stored in the hdur() array in the case of FORCESOLUTION,
      ! to use the same array as for CMTSOLUTION
      stf = comp_source_time_function_rickr(time_source_dble,hdur(isource))
    else
      ! Gaussian
      ! use a very small duration of 5*DT to mimic a Dirac in time
      stf = comp_source_time_function_gauss(time_source_dble,5.d0*DT)
    endif
  else
    ! moment-tensor
    if (USE_RICKER_TIME_FUNCTION) then
      ! Ricker
      stf = comp_source_time_function_rickr(time_source_dble,hdur(isource))
    else
      ! Gaussian
      ! since the source is a bulk source (applied to both fluid and solid parts)
      stf = comp_source_time_function_gauss(time_source_dble,hdur_Gaussian(isource))
    endif
  endif ! USE_FORCE_POINT_SOURCE

  ! return value
  get_stf_poroelastic = stf

  end function get_stf_poroelastic


