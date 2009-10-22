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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine setup_sources_receivers()

  use specfem_par
  implicit none

! write source and receiver VTK files for Paraview
  if (myrank == 0) then
    open(IOVTK,file=trim(OUTPUT_FILES)//'/sr.vtk',status='unknown')
    write(IOVTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOVTK,'(a)') 'Source and Receiver VTK file'
    write(IOVTK,'(a)') 'ASCII'
    write(IOVTK,'(a)') 'DATASET POLYDATA'
    ! LQY -- cannot figure out NSOURCES+nrec at this point
    write(IOVTK, '(a,i6,a)') 'POINTS ', 2, ' float'
  endif

! allocate arrays for source
  allocate(islice_selected_source(NSOURCES))
  allocate(ispec_selected_source(NSOURCES))
  allocate(Mxx(NSOURCES))
  allocate(Myy(NSOURCES))
  allocate(Mzz(NSOURCES))
  allocate(Mxy(NSOURCES))
  allocate(Mxz(NSOURCES))
  allocate(Myz(NSOURCES))
  allocate(xi_source(NSOURCES))
  allocate(eta_source(NSOURCES))
  allocate(gamma_source(NSOURCES))
  allocate(t_cmt(NSOURCES))
  allocate(hdur(NSOURCES))
  allocate(hdur_gaussian(NSOURCES))
  allocate(utm_x_source(NSOURCES))
  allocate(utm_y_source(NSOURCES))
  allocate(nu_source(3,3,NSOURCES))

! locate sources in the mesh
  call locate_source(ibool,NSOURCES,myrank,NSPEC_AB,NGLOB_AB, &
          xstore,ystore,zstore,xigll,yigll,zigll,NPROC, &
          sec,t_cmt,yr,jda,ho,mi,utm_x_source,utm_y_source, &
          NSTEP,DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
          islice_selected_source,ispec_selected_source, &
          xi_source,eta_source,gamma_source, &
          TOPOGRAPHY,UTM_PROJECTION_ZONE, &
          PRINT_SOURCE_TIME_FUNCTION, &
          nu_source,iglob_is_surface_external_mesh,ispec_is_surface_external_mesh)

  if(minval(t_cmt) /= 0.) call exit_MPI(myrank,'one t_cmt must be zero, others must be positive')

! filter source time function by Gaussian with hdur = HDUR_MOVIE when outputing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME .or. CREATE_SHAKEMAP) then
     hdur = sqrt(hdur**2 + HDUR_MOVIE**2)
     if(myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Each source is being convolved with HDUR_MOVIE = ',HDUR_MOVIE
        write(IMAIN,*)
     endif
  endif
! convert the half duration for triangle STF to the one for gaussian STF
  hdur_gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE

! define t0 as the earliest start time
  t0 = - 1.5d0 * minval(t_cmt-hdur)

!$$$$$$$$$$$$$$$$$$ RECEIVERS $$$$$$$$$$$$$$$$$$$$$

  if (SIMULATION_TYPE == 1) then
    call get_value_string(rec_filename, 'solver.STATIONS', 'DATA/STATIONS')

! get total number of stations
    open(unit=IIN,file=rec_filename,iostat=ios,status='old',action='read')
    nrec = 0
    do while(ios == 0)
      read(IIN,"(a)",iostat=ios) dummystring
      if(ios == 0) nrec = nrec + 1
    enddo
    close(IIN)
    if(nrec < 1) call exit_MPI(myrank,'need at least one receiver')

  else
    call get_value_string(rec_filename, 'solver.STATIONS', 'DATA/STATIONS_ADJOINT')
    call get_value_string(filtered_rec_filename, 'solver.STATIONS_FILTERED', 'DATA/STATIONS_ADJOINT_FILTERED')
    call station_filter(myrank,rec_filename,filtered_rec_filename,nrec, &
           LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)
    if (nrec < 1) call exit_MPI(myrank, 'adjoint simulation needs at least one receiver')
    call sync_all()
  endif

  if(myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) 'Total number of receivers = ', nrec
    else
      write(IMAIN,*) 'Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
  endif

  if(nrec < 1) call exit_MPI(myrank,'need at least one receiver')

! allocate memory for receiver arrays
  allocate(islice_selected_rec(nrec))
  allocate(ispec_selected_rec(nrec))
  allocate(xi_receiver(nrec))
  allocate(eta_receiver(nrec))
  allocate(gamma_receiver(nrec))
  allocate(station_name(nrec))
  allocate(network_name(nrec))
  allocate(nu(NDIM,NDIM,nrec))

! locate receivers in the mesh
  call locate_receivers(ibool,myrank,NSPEC_AB,NGLOB_AB, &
            xstore,ystore,zstore,xigll,yigll,zigll,rec_filename, &
            nrec,islice_selected_rec,ispec_selected_rec, &
            xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
            NPROC,utm_x_source(1),utm_y_source(1), &
            TOPOGRAPHY,UTM_PROJECTION_ZONE, &
            iglob_is_surface_external_mesh,ispec_is_surface_external_mesh )


!###################### SOURCE ARRAYS ################

  if (SIMULATION_TYPE == 1  .or. SIMULATION_TYPE == 3) then
    allocate(sourcearray(NDIM,NGLLX,NGLLY,NGLLZ))
    allocate(sourcearrays(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ))

! compute source arrays
    do isource = 1,NSOURCES

!   check that the source slice number is okay
      if(islice_selected_source(isource) < 0 .or. islice_selected_source(isource) > NPROC-1) &
            call exit_MPI(myrank,'something is wrong with the source slice number')

!   compute source arrays in source slice
      if(myrank == islice_selected_source(isource)) then
        call compute_arrays_source(ispec_selected_source(isource), &
              xi_source(isource),eta_source(isource),gamma_source(isource),sourcearray, &
              Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource), &
              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
              xigll,yigll,zigll,NSPEC_AB)
        sourcearrays(isource,:,:,:,:) = sourcearray(:,:,:,:)
      endif
    enddo
  endif

  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    nadj_rec_local = 0
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec))then
!   check that the source slice number is okay
        if(islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROC-1) &
              call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo
    allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLY,NGLLZ))
    if (nadj_rec_local > 0) allocate(adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ))
    irec_local = 0
    do irec = 1, nrec
!   compute only adjoint source arrays in the local slice
      if(myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
        call compute_arrays_adjoint_source(myrank, adj_source_file, &
              xi_receiver(irec), eta_receiver(irec), gamma_receiver(irec), &
              adj_sourcearray, xigll,yigll,zigll,NSTEP)

        adj_sourcearrays(irec_local,:,:,:,:,:) = adj_sourcearray(:,:,:,:,:)

      endif
    enddo
  endif

!--- select local receivers

! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    nrec_simulation = nrec
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if(myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

  if (nrec_local > 0) then
  ! allocate Lagrange interpolators for receivers
    allocate(hxir_store(nrec_local,NGLLX))
    allocate(hetar_store(nrec_local,NGLLY))
    allocate(hgammar_store(nrec_local,NGLLZ))

  ! define local to global receiver numbering mapping
    allocate(number_receiver_global(nrec_local))
    irec_local = 0
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      do irec = 1,nrec
      if(myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        number_receiver_global(irec_local) = irec
      endif
      enddo
    else
      do isource = 1,NSOURCES
        if(myrank == islice_selected_source(isource)) then
          irec_local = irec_local + 1
          number_receiver_global(irec_local) = isource
        endif
      enddo
    endif

  ! define and store Lagrange interpolators at all the receivers
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      do irec_local = 1,nrec_local
        irec = number_receiver_global(irec_local)
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
        hxir_store(irec_local,:) = hxir(:)
        hetar_store(irec_local,:) = hetar(:)
        hgammar_store(irec_local,:) = hgammar(:)
      enddo
    else
      allocate(hpxir_store(nrec_local,NGLLX))
      allocate(hpetar_store(nrec_local,NGLLY))
      allocate(hpgammar_store(nrec_local,NGLLZ))
      do irec_local = 1,nrec_local
        irec = number_receiver_global(irec_local)
        call lagrange_any(xi_source(irec),NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta_source(irec),NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma_source(irec),NGLLZ,zigll,hgammar,hpgammar)
        hxir_store(irec_local,:) = hxir(:)
        hetar_store(irec_local,:) = hetar(:)
        hgammar_store(irec_local,:) = hgammar(:)
        hpxir_store(irec_local,:) = hpxir(:)
        hpetar_store(irec_local,:) = hpetar(:)
        hpgammar_store(irec_local,:) = hpgammar(:)
      enddo
    endif
  endif ! nrec_local > 0

! check that the sum of the number of receivers in each slice is nrec
  call sum_all_i(nrec_local,nrec_tot_found)
  if(myrank == 0) then

    close(IOVTK)

    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all the slices'
    if(nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    else
      write(IMAIN,*) 'this total is okay'
    endif
    
    if(NSOURCES > 1) write(IMAIN,*) 'Using ',NSOURCES,' point sources'
    
  endif

  end subroutine