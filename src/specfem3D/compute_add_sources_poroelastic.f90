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

! for poroelastic solver

  subroutine compute_add_sources_poroelastic()

  use constants
  use specfem_par, only: station_name,network_name,USE_FORCE_POINT_SOURCE, &
                         tshift_src,dt,t0,USE_LDDRK,istage,USE_EXTERNAL_SOURCE_FILE,user_source_time_function, &
                         USE_BINARY_FOR_SEISMOGRAMS,ibool, &
                         UNDO_ATTENUATION_AND_OR_PML, &
                         NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                         sourcearrays,SIMULATION_TYPE,NSTEP, &
                         ispec_selected_rec, &
                         nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC, &
                         hxir_adjstore,hetar_adjstore,hgammar_adjstore,source_adjoint,number_adjsources_global,nadj_rec_local

  use specfem_par_poroelastic, only: b_accels_poroelastic,b_accelw_poroelastic,accels_poroelastic,accelw_poroelastic, &
                                      rhoarraystore,phistore,tortstore,ispec_is_poroelastic

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: stf_used,hlagrange

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_poroelastic

  logical :: ibool_read_adj_arrays
  integer :: isource,iglob,i,j,k,ispec,it_sub_adj
  integer :: irec_local,irec
  real(kind=CUSTOM_REAL) :: phil,tortl,rhol_s,rhol_f,rhol_bar
  real(kind=CUSTOM_REAL) :: fac_s,fac_w

  character(len=MAX_STRING_LEN) :: adj_source_file

  ! forward simulations
  if (SIMULATION_TYPE == 1) then
    ! ignore CMT sources for fault rupture simulations
    if (FAULT_SIMULATION) return

! openmp solver
!$OMP PARALLEL if (NSOURCES > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(isource,time_source_dble,stf_used,stf,iglob,ispec,i,j,k, &
!$OMP         phil,tortl,rhol_s,rhol_f,rhol_bar,fac_s,fac_w)

    ! adds poroelastic sources
!$OMP DO
    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_poroelastic(ispec)) then
          ! current time
          if (USE_LDDRK) then
            ! LDDRK
            ! note: the LDDRK scheme updates displacement after the stiffness computations and
            !       after adding boundary/coupling/source terms.
            !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
            !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
            time_source_dble = dble(it-1-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(it-1)*DT - t0 - tshift_src(isource)
          endif

          ! determines source time function value
          stf = get_stf_poroelastic(time_source_dble,isource)

          !! VM VM add external source time function
          if (USE_EXTERNAL_SOURCE_FILE) then
            stf = user_source_time_function(it, isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          ! adds source array
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
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
!$OMP ATOMIC
                accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) &
                             + fac_s * sourcearrays(isource,1,i,j,k)*stf_used
!$OMP ATOMIC
                accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) &
                             + fac_s * sourcearrays(isource,2,i,j,k)*stf_used
!$OMP ATOMIC
                accels_poroelastic(3,iglob) = accels_poroelastic(3,iglob) &
                             + fac_s * sourcearrays(isource,3,i,j,k)*stf_used

                ! fluid phase
!$OMP ATOMIC
                accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) &
                             + fac_w * sourcearrays(isource,1,i,j,k)*stf_used
!$OMP ATOMIC
                accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) &
                             + fac_w * sourcearrays(isource,2,i,j,k)*stf_used
!$OMP ATOMIC
                accelw_poroelastic(3,iglob) = accelw_poroelastic(3,iglob) &
                             + fac_w * sourcearrays(isource,3,i,j,k)*stf_used

              enddo
            enddo
          enddo

        endif ! ispec_is_poroelastic
      endif ! myrank
    enddo ! NSOURCES
!$OMP ENDDO
!$OMP END PARALLEL

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
        if (USE_BINARY_FOR_SEISMOGRAMS) stop 'Adjoint simulations not supported with .bin format, please use ASCII instead'
        ! ASCII format
        !!! read ascii adjoint sources
        do irec_local = 1, nadj_rec_local
          irec = number_adjsources_global(irec_local)
          ! compute source arrays
          ! reads in **net**.**sta**.**BH**.adj files
          adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
          call compute_arrays_adjoint_source(adj_source_file,irec_local)
        enddo
      endif ! if (ibool_read_adj_arrays)

      if (it < NSTEP) then
        ! receivers act as sources
        do irec_local = 1, nadj_rec_local
          irec = number_adjsources_global(irec_local)
          ! element index
          ispec = ispec_selected_rec(irec)
          if (ispec_is_poroelastic(ispec)) then
            ! adds source array
            do k = 1,NGLLZ
              do j = 1,NGLLY
                do i = 1,NGLLX
                  iglob = ibool(i,j,k,ispec)
                  ! get poroelastic parameters of current local GLL
                  phil = phistore(i,j,k,ispec)
                  rhol_s = rhoarraystore(1,i,j,k,ispec)
                  rhol_f = rhoarraystore(2,i,j,k,ispec)
                  rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

                  hlagrange = hxir_adjstore(i,irec_local) * hetar_adjstore(j,irec_local) * hgammar_adjstore(k,irec_local)

                  ! adjoint source is in the solid phase only since this is the only measurement
                  ! available

                  ! solid phase
                  accels_poroelastic(:,iglob) = accels_poroelastic(:,iglob) &
                    + source_adjoint(:,irec_local,NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)) * hlagrange
                  !
                  ! fluid phase
                  accelw_poroelastic(:,iglob) = accelw_poroelastic(:,iglob) &
                    - rhol_f/rhol_bar &
                    * source_adjoint(:,irec_local,NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)) * hlagrange
                enddo
              enddo
            enddo
          endif ! ispec_is_poroelastic

        enddo ! nrec

      endif ! it
    endif ! nadj_rec_local
  endif !adjoint : if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then

! note:  b_displ() is read in after Newmark time scheme, thus
!           b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!           thus indexing is NSTEP - it , instead of NSTEP - it - 1

  ! adjoint/backward simulations
  if (SIMULATION_TYPE == 3) then
    ! ignore CMT sources for fault rupture simulations
    if (FAULT_SIMULATION) return

    ! backward source reconstruction
    do isource = 1,NSOURCES

      ! add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_poroelastic(ispec)) then
          ! note: time step is now at NSTEP-it
          ! current time
          if (USE_LDDRK) then
            ! LDDRK
            ! note: the LDDRK scheme updates displacement after the stiffness computations and
            !       after adding boundary/coupling/source terms.
            !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
            !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
            if (UNDO_ATTENUATION_AND_OR_PML) then
              ! stepping moves forward from snapshot position
              time_source_dble = dble(NSTEP-it-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
            else
              time_source_dble = dble(NSTEP-it-1)*DT - dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
            endif
          else
            time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
          endif

          ! determines source time function value
          stf = get_stf_poroelastic(time_source_dble,isource)

          !! VM VM add external source time function
          if (USE_EXTERNAL_SOURCE_FILE) then
            ! time-reversed
            stf = user_source_time_function(NSTEP-it+1, isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          !  add source array
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
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
                b_accels_poroelastic(:,iglob) = b_accels_poroelastic(:,iglob) &
                             + fac_s * sourcearrays(isource,:,i,j,k)*stf_used
                ! fluid phase
                b_accelw_poroelastic(:,iglob) = b_accelw_poroelastic(:,iglob) &
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


