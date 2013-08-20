!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
                                  hdur,hdur_gaussian,tshift_src,dt,t0, &
                                  sourcearrays,kappastore,ispec_is_acoustic,&
                                  SIMULATION_TYPE,NSTEP, &
                                  nrec,islice_selected_rec,ispec_selected_rec, &
                                  nadj_rec_local,adj_sourcearrays,NTSTEP_BETWEEN_READ_ADJSRC)

  use specfem_par,only: PRINT_SOURCE_TIME_FUNCTION,stf_used_total, &
                        xigll,yigll,zigll,xi_receiver,eta_receiver,gamma_receiver,&
                        station_name,network_name,adj_source_file,nrec_local,number_receiver_global, &
                        pm1_source_encoding,nsources_local,USE_FORCE_POINT_SOURCE, &
                        USE_RICKER_TIME_FUNCTION
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
  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,hdur_tiny,tshift_src
  double precision :: dt,t0
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr,&
   comp_source_time_function_gauss

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP
  integer:: nrec
  integer,dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer:: nadj_rec_local
  logical :: ibool_read_adj_arrays
  integer :: it_sub_adj,itime,NTSTEP_BETWEEN_READ_ADJSRC
  real(kind=CUSTOM_REAL),dimension(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ):: &
    adj_sourcearrays

! local parameters
  double precision :: f0
  double precision :: stf
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: adj_sourcearray
  real(kind=CUSTOM_REAL) stf_used,stf_used_total_all,time_source
  integer :: isource,iglob,ispec,i,j,k,ier
  integer :: irec_local,irec

! adjoint sources in SU format
  integer :: it_start,it_end
  real(kind=CUSTOM_REAL) :: adj_temp(NSTEP)
  real(kind=CUSTOM_REAL) :: adj_src(NTSTEP_BETWEEN_READ_ADJSRC,NDIM)
  character(len=256) :: procname
  integer,parameter :: nheader=240      ! 240 bytes
  !integer(kind=2) :: i2head(nheader/2)  ! 2-byte-integer
  !integer(kind=4) :: i4head(nheader/4)  ! 4-byte-integer
  real(kind=4)    :: r4head(nheader/4)  ! 4-byte-real
  !equivalence (i2head,i4head,r4head)    ! share the same 240-byte-memory
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY),hgammar(NGLLZ), hpgammar(NGLLZ)

! plotting source time function
  if(PRINT_SOURCE_TIME_FUNCTION .and. .not. phase_is_inner ) then
    ! initializes total
    stf_used_total = 0.0_CUSTOM_REAL
  endif

! forward simulations
  if (SIMULATION_TYPE == 1 .and. nsources_local > 0) then

!way 2
     ! adds acoustic sources
     do isource = 1,NSOURCES

        !   add the source (only if this proc carries the source)
        if(myrank == islice_selected_source(isource)) then

           ispec = ispec_selected_source(isource)

           if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

              if( ispec_is_acoustic(ispec) ) then

                 if(USE_FORCE_POINT_SOURCE) then

                    f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing FORCESOLUTION file format

                    !if (it == 1 .and. myrank == 0) then
                    !  write(IMAIN,*) 'using a source of dominant frequency ',f0
                    !  write(IMAIN,*) 'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
                    !  write(IMAIN,*) 'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
                    !endif

                    if( USE_RICKER_TIME_FUNCTION ) then
                       stf_used = comp_source_time_function_rickr(dble(it-1)*DT-t0-tshift_src(isource),f0)
                    else
                       stf_used = comp_source_time_function_gauss(dble(it-1)*DT-t0-tshift_src(isource),hdur_tiny(isource))
                    endif

                    ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
                    ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
                    ! to add minus the source to Chi_dot_dot to get plus the source in pressure:

                    ! acoustic source for pressure gets divided by kappa
                    ! source contribution
                    do k=1,NGLLZ
                      do j=1,NGLLY
                        do i=1,NGLLX
                          iglob = ibool(i,j,k,ispec)
                          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                  - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)
                        enddo
                      enddo
                    enddo

                 else

                    if( USE_RICKER_TIME_FUNCTION ) then
                       stf = comp_source_time_function_rickr(dble(it-1)*DT-t0-tshift_src(isource),hdur(isource))
                    else
                       ! gaussian source time
                       stf = comp_source_time_function_gauss(dble(it-1)*DT-t0-tshift_src(isource),hdur_gaussian(isource))
                    endif

                    ! quasi-Heaviside
                    !stf = comp_source_time_function(dble(it-1)*DT-t0-tshift_src(isource),hdur_gaussian(isource))

                    ! source encoding
                    stf = stf * pm1_source_encoding(isource)

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

              endif ! ispec_is_acoustic
           endif ! ispec_is_inner
        endif ! myrank
     enddo ! NSOURCES
  endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1 )*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then

    ! adds adjoint source in this partitions
    if( nadj_rec_local > 0 ) then

      ! read in adjoint sources block by block (for memory consideration)
      ! e.g., in exploration experiments, both the number of receivers (nrec) and
      ! the number of time steps (NSTEP) are huge,
      ! which may cause problems since we have a large array:
      !     adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries and interior --- indicated by 'phase_is_inner'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while we calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added twice
      if (ibool_read_adj_arrays .and. (.not. phase_is_inner)) then

        ! allocates temporary source array
        allocate(adj_sourcearray(NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
        if( ier /= 0 ) stop 'error allocating array adj_sourcearray'

        if (.not. SU_FORMAT) then
           !!! read ascii adjoint sources
           irec_local = 0
           do irec = 1, nrec
             ! compute source arrays
             if (myrank == islice_selected_rec(irec)) then
               irec_local = irec_local + 1

               ! reads in **sta**.**net**.**LH**.adj files
               adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
               call compute_arrays_adjoint_source(myrank,adj_source_file, &
                                                  xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                                                  adj_sourcearray, xigll,yigll,zigll, &
                                                  it_sub_adj,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC)

               do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
                 adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
               enddo

             endif
           enddo
        else
          !!! read SU adjoint sources
          ! range of the block we need to read
          it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC + 1
          it_end   = it_start + NTSTEP_BETWEEN_READ_ADJSRC - 1
          write(procname,"(i4)") myrank
          open(unit=IIN_SU1, file=trim(adjustl(OUTPUT_FILES_PATH))//'../SEM/'//trim(adjustl(procname))//'_dx_SU.adj', &
                            status='old',access='direct',recl=240+4*(NSTEP),iostat = ier)
          if( ier /= 0 ) call exit_MPI(myrank,'file '//trim(adjustl(OUTPUT_FILES_PATH)) &
                                    //'../SEM/'//trim(adjustl(procname))//'_dx_SU.adj does not exit')

          do irec_local = 1,nrec_local
            irec = number_receiver_global(irec_local)
            read(IIN_SU1,rec=irec_local) r4head, adj_temp
            adj_src(:,1)=adj_temp(it_start:it_end)
            adj_src(:,2)=0.0  !TRIVIAL
            adj_src(:,3)=0.0  !TRIVIAL
            ! lagrange interpolators for receiver location
            call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
            call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
            call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
            ! interpolates adjoint source onto GLL points within this element
            do k = 1, NGLLZ
              do j = 1, NGLLY
                do i = 1, NGLLX
                  adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
                enddo
              enddo
            enddo
            do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
              adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
            enddo
          enddo
          close(IIN_SU1)
        endif !if (.not. SU_FORMAT)

        deallocate(adj_sourcearray)
      endif ! if(ibool_read_adj_arrays)

      if( it < NSTEP ) then
        ! receivers act as sources
        irec_local = 0
        do irec = 1,nrec
          ! add the source (only if this proc carries the source)
          if (myrank == islice_selected_rec(irec)) then
            irec_local = irec_local + 1

            ! adds source array
            ispec = ispec_selected_rec(irec)
            if( ispec_is_acoustic(ispec) ) then

              ! checks if element is in phase_is_inner run
              if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
                do k = 1,NGLLZ
                  do j = 1,NGLLY
                    do i = 1,NGLLX
                      iglob = ibool(i,j,k,ispec)
                      ! beware, for acoustic medium, a pressure source would be taking the negative
                      ! and divide by Kappa of the fluid;
                      ! this would have to be done when constructing the adjoint source.
                      !
                      ! note: we take the first component of the adj_sourcearrays
                      !          the idea is to have e.g. a pressure source, where all 3 components would be the same
                      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                 + adj_sourcearrays(irec_local, &
                                                    NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC), &
                                                    1,i,j,k)
                    enddo
                  enddo
                enddo
              endif ! phase_is_inner
            endif
          endif
        enddo ! nrec
    endif ! it
  endif ! nadj_rec_local > 0
endif

  ! master prints out source time function to file
  if(PRINT_SOURCE_TIME_FUNCTION .and. phase_is_inner) then
    time_source = (it-1)*DT - t0
    call sum_all_cr(stf_used_total,stf_used_total_all)
    if( myrank == 0 ) write(IOSTF,*) time_source,stf_used_total_all
  endif

  end subroutine compute_add_sources_acoustic

!
!=====================================================================
! for acoustic solver for back propagation wave field

  subroutine compute_add_sources_acoustic_bpwf(NSPEC_AB, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  hdur,hdur_gaussian,tshift_src,dt,t0, &
                                  sourcearrays,kappastore,ispec_is_acoustic,&
                                  SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                                  b_potential_dot_dot_acoustic)

  use specfem_par,only: PRINT_SOURCE_TIME_FUNCTION,stf_used_total, &
                        pm1_source_encoding,nsources_local,USE_FORCE_POINT_SOURCE, &
                        USE_RICKER_TIME_FUNCTION
  implicit none

  include "constants.h"

  integer :: NSPEC_AB

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,hdur_tiny,tshift_src
  double precision :: dt,t0
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr,&
   comp_source_time_function_gauss

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT
  real(kind=CUSTOM_REAL),dimension(NGLOB_ADJOINT):: b_potential_dot_dot_acoustic

! local parameters
  double precision :: f0
  double precision :: stf
  real(kind=CUSTOM_REAL) stf_used,stf_used_total_all,time_source
  integer :: isource,iglob,ispec,i,j,k

! adjoint sources in SU format
  integer,parameter :: nheader=240      ! 240 bytes

  ! checks if anything to do
  if( SIMULATION_TYPE /= 3 ) return

! plotting source time function
  if(PRINT_SOURCE_TIME_FUNCTION .and. .not. phase_is_inner ) then
    ! initializes total
    stf_used_total = 0.0_CUSTOM_REAL
  endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1 )*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! note:  b_potential() is read in after Newmark time scheme, thus
!           b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!           thus indexing is NSTEP - it , instead of NSTEP - it - 1

! adjoint simulations
  if( nsources_local > 0 ) then

     ! adds acoustic sources
     do isource = 1,NSOURCES

        !   add the source (only if this proc carries the source)
        if(myrank == islice_selected_source(isource)) then

           ispec = ispec_selected_source(isource)

           if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

              if( ispec_is_acoustic(ispec) ) then

                 if(USE_FORCE_POINT_SOURCE) then

                    f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing FORCESOLUTION file format

                    !if (it == 1 .and. myrank == 0) then
                    !  write(IMAIN,*) 'using a source of dominant frequency ',f0
                    !  write(IMAIN,*) 'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
                    !  write(IMAIN,*) 'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
                    !endif

                    if( USE_RICKER_TIME_FUNCTION ) then
                       stf_used = comp_source_time_function_rickr( &
                                  dble(NSTEP-it)*DT-t0-tshift_src(isource),f0)
                    else
                       stf_used = comp_source_time_function_gauss( &
                                  dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur_tiny(isource))
                    endif

                    ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
                    ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
                    ! to add minus the source to Chi_dot_dot to get plus the source in pressure:

                    ! acoustic source for pressure gets divided by kappa
                    ! source contribution
                    do k=1,NGLLZ
                      do j=1,NGLLY
                        do i=1,NGLLX
                          iglob = ibool(i,j,k,ispec)
                          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                  - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)
                        enddo
                      enddo
                    enddo

                 else
                    if( USE_RICKER_TIME_FUNCTION ) then
                       stf = comp_source_time_function_rickr( &
                             dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur(isource))
                    else
                       ! gaussian source time
                       stf = comp_source_time_function_gauss( &
                             dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur_gaussian(isource))
                    endif

                    ! quasi-Heaviside
                    !stf = comp_source_time_function(dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur_gaussian(isource))

                    ! source encoding
                    stf = stf * pm1_source_encoding(isource)

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

  end subroutine compute_add_sources_acoustic_bpwf

!
!=====================================================================
! for acoustic solver on GPU
  subroutine compute_add_sources_acoustic_GPU(NSPEC_AB,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,&
                                  hdur,hdur_gaussian,tshift_src,dt,t0, &
                                  ispec_is_acoustic,SIMULATION_TYPE,NSTEP, &
                                  nrec,islice_selected_rec,ispec_selected_rec, &
                                  nadj_rec_local,adj_sourcearrays, &
                                  NTSTEP_BETWEEN_READ_ADJSRC,Mesh_pointer )

  use specfem_par,only: PRINT_SOURCE_TIME_FUNCTION,stf_used_total, &
                        xigll,yigll,zigll,xi_receiver,eta_receiver,gamma_receiver,&
                        station_name,network_name,adj_source_file,nrec_local,number_receiver_global, &
                        nsources_local,USE_FORCE_POINT_SOURCE, &
                        USE_RICKER_TIME_FUNCTION
  implicit none

  include "constants.h"

  integer :: NSPEC_AB

! displacement and acceleration

! arrays with mesh parameters per slice

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! source
  integer :: NSOURCES,myrank,it
  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,hdur_tiny,tshift_src
  double precision :: dt,t0

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr,&
   comp_source_time_function_gauss

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP
  integer(kind=8) :: Mesh_pointer
  integer:: nrec
  integer,dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer:: nadj_rec_local
  logical :: ibool_read_adj_arrays
  integer :: it_sub_adj,itime,NTSTEP_BETWEEN_READ_ADJSRC
  real(kind=CUSTOM_REAL),dimension(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ):: &
    adj_sourcearrays

! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: adj_sourcearray
  real(kind=CUSTOM_REAL) stf_used_total_all,time_source
  integer :: isource,i,j,k,ier
  integer :: irec_local,irec
  double precision, dimension(NSOURCES) :: stf_pre_compute

! adjoint sources in SU format
  integer :: it_start,it_end
  real(kind=CUSTOM_REAL) :: adj_temp(NSTEP)
  real(kind=CUSTOM_REAL) :: adj_src(NTSTEP_BETWEEN_READ_ADJSRC,NDIM)
  character(len=256) :: procname
  integer,parameter :: nheader=240      ! 240 bytes
  !integer(kind=2) :: i2head(nheader/2)  ! 2-byte-integer
  !integer(kind=4) :: i4head(nheader/4)  ! 4-byte-integer
  real(kind=4)    :: r4head(nheader/4)  ! 4-byte-real
  !equivalence (i2head,i4head,r4head)    ! share the same 240-byte-memory
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY),hgammar(NGLLZ), hpgammar(NGLLZ)

! plotting source time function
  if(PRINT_SOURCE_TIME_FUNCTION .and. .not. phase_is_inner ) then
    ! initializes total
    stf_used_total = 0.0_CUSTOM_REAL
  endif

! forward simulations
  if (SIMULATION_TYPE == 1 .and. nsources_local > 0) then

!way 2
      if( NSOURCES > 0 ) then
         do isource = 1,NSOURCES
            ! precomputes source time function factor
            if(USE_FORCE_POINT_SOURCE) then
               if( USE_RICKER_TIME_FUNCTION ) then
                  stf_pre_compute(isource) = &
                       comp_source_time_function_rickr(dble(it-1)*DT-t0-tshift_src(isource),hdur(isource))
               else
                  stf_pre_compute(isource) = &
                       comp_source_time_function_gauss(dble(it-1)*DT-t0-tshift_src(isource),hdur_tiny(isource))
               endif
            else
               if( USE_RICKER_TIME_FUNCTION ) then
                  stf_pre_compute(isource) = &
                       comp_source_time_function_rickr(dble(it-1)*DT-t0-tshift_src(isource),hdur(isource))
               else
                  stf_pre_compute(isource) = &
                       comp_source_time_function_gauss(dble(it-1)*DT-t0-tshift_src(isource),hdur_gaussian(isource))
               endif
            endif
         enddo
         stf_used_total = stf_used_total + sum(stf_pre_compute(:))
         ! only implements SIMTYPE=1 and NOISE_TOM=0
         ! write(*,*) "fortran dt = ", dt
         ! change dt -> DT
         call compute_add_sources_ac_cuda(Mesh_pointer,phase_is_inner,NSOURCES,stf_pre_compute)
      endif
  endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1 )*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then

    ! adds adjoint source in this partitions
    if( nadj_rec_local > 0 ) then

      ! read in adjoint sources block by block (for memory consideration)
      ! e.g., in exploration experiments, both the number of receivers (nrec) and
      ! the number of time steps (NSTEP) are huge,
      ! which may cause problems since we have a large array:
      !     adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries and interior --- indicated by 'phase_is_inner'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while we calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added twice
      if (ibool_read_adj_arrays .and. (.not. phase_is_inner)) then

        ! allocates temporary source array
        allocate(adj_sourcearray(NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
        if( ier /= 0 ) stop 'error allocating array adj_sourcearray'

        if (.not. SU_FORMAT) then
           !!! read ascii adjoint sources
           irec_local = 0
           do irec = 1, nrec
             ! compute source arrays
             if (myrank == islice_selected_rec(irec)) then
               irec_local = irec_local + 1

               ! reads in **sta**.**net**.**LH**.adj files
               adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
               call compute_arrays_adjoint_source(myrank,adj_source_file, &
                                                  xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                                                  adj_sourcearray, xigll,yigll,zigll, &
                                                  it_sub_adj,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC)

               do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
                 adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
               enddo

             endif
           enddo
        else
          !!! read SU adjoint sources
          ! range of the block we need to read
          it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC + 1
          it_end   = it_start + NTSTEP_BETWEEN_READ_ADJSRC - 1
          write(procname,"(i4)") myrank
          open(unit=IIN_SU1, file=trim(adjustl(OUTPUT_FILES_PATH))//'../SEM/'//trim(adjustl(procname))//'_dx_SU.adj', &
                            status='old',access='direct',recl=240+4*(NSTEP),iostat = ier)
          if( ier /= 0 ) call exit_MPI(myrank,'file '//trim(adjustl(OUTPUT_FILES_PATH)) &
                                    //'../SEM/'//trim(adjustl(procname))//'_dx_SU.adj does not exit')

          do irec_local = 1,nrec_local
            irec = number_receiver_global(irec_local)
            read(IIN_SU1,rec=irec_local) r4head, adj_temp
            adj_src(:,1)=adj_temp(it_start:it_end)
            adj_src(:,2)=0.0  !TRIVIAL
            adj_src(:,3)=0.0  !TRIVIAL
            ! lagrange interpolators for receiver location
            call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
            call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
            call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
            ! interpolates adjoint source onto GLL points within this element
            do k = 1, NGLLZ
              do j = 1, NGLLY
                do i = 1, NGLLX
                  adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
                enddo
              enddo
            enddo
            do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
              adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
            enddo
          enddo
          close(IIN_SU1)
        endif !if (.not. SU_FORMAT)

        deallocate(adj_sourcearray)
      endif ! if(ibool_read_adj_arrays)

      if( it < NSTEP ) then
        ! receivers act as sources
        ! on GPU
        call add_sources_ac_sim_2_or_3_cuda(Mesh_pointer,adj_sourcearrays,phase_is_inner, &
                                           ispec_is_inner,ispec_is_acoustic, &
                                           ispec_selected_rec, &
                                           nrec, &
                                           NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC), &
                                           islice_selected_rec,nadj_rec_local, &
                                           NTSTEP_BETWEEN_READ_ADJSRC)
      endif ! it
    endif ! nadj_rec_local > 0
  endif

! note:  b_potential() is read in after Newmark time scheme, thus
!           b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!           thus indexing is NSTEP - it , instead of NSTEP - it - 1

! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. nsources_local > 0) then

      if( NSOURCES > 0 ) then
         do isource = 1,NSOURCES
            ! precomputes source time function factors
            if(USE_FORCE_POINT_SOURCE) then
               if( USE_RICKER_TIME_FUNCTION ) then
                  stf_pre_compute(isource) = &
                       comp_source_time_function_rickr(dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur(isource))
               else
                  stf_pre_compute(isource) = &
                       comp_source_time_function_gauss(dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur_tiny(isource))
               endif
            else
               if( USE_RICKER_TIME_FUNCTION ) then
                  stf_pre_compute(isource) = &
                       comp_source_time_function_rickr(dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur(isource))
               else
                  stf_pre_compute(isource) = &
                       comp_source_time_function_gauss(dble(NSTEP-it)*DT-t0-tshift_src(isource),hdur_gaussian(isource))
               endif
            endif
         enddo
         stf_used_total = stf_used_total + sum(stf_pre_compute(:))
         ! only implements SIMTYPE=3
         call compute_add_sources_ac_s3_cuda(Mesh_pointer,phase_is_inner,NSOURCES,stf_pre_compute)
      endif
  endif

  ! master prints out source time function to file
  if(PRINT_SOURCE_TIME_FUNCTION .and. phase_is_inner) then
    time_source = (it-1)*DT - t0
    call sum_all_cr(stf_used_total,stf_used_total_all)
    if( myrank == 0 ) write(IOSTF,*) time_source,stf_used_total_all
  endif

  end subroutine compute_add_sources_acoustic_GPU
