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

! locates sources and determines simulation start time t0
  call setup_sources()
 
! reads in stations file and locates receivers
  call setup_receivers()

! pre-compute source arrays
  call setup_sources_precompute_arrays()  

! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

! write source and receiver VTK files for Paraview
  call setup_sources_receivers_VTKfile()

! user output  
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all the slices'
    if(NSOURCES > 1) write(IMAIN,*) 'Using ',NSOURCES,' point sources'    
  endif

end subroutine setup_sources_receivers
  
!
!-------------------------------------------------------------------------------------------------
!  
  
subroutine setup_sources()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic  
  use specfem_par_movie
  implicit none
  
  double precision :: t0_ac
  integer :: yr,jda,ho,mi
  integer :: isource,ispec
  
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
!
! returns:  islice_selected_source & ispec_selected_source,
!                xi_source, eta_source & gamma_source 
  call locate_source(ibool,NSOURCES,myrank,NSPEC_AB,NGLOB_AB, &
          xstore,ystore,zstore,xigll,yigll,zigll,NPROC, &
          t_cmt,yr,jda,ho,mi,utm_x_source,utm_y_source, &
          DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
          islice_selected_source,ispec_selected_source, &
          xi_source,eta_source,gamma_source, &
          UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
          PRINT_SOURCE_TIME_FUNCTION, &
          nu_source,iglob_is_surface_external_mesh,ispec_is_surface_external_mesh,&
          ispec_is_acoustic,ispec_is_elastic, &
          num_free_surface_faces,free_surface_ispec,free_surface_ijk)

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

  ! uses an earlier start time if source is acoustic with a gaussian source time function
  t0_ac = 0.0d0
  do isource = 1,NSOURCES  
    if( myrank == islice_selected_source(isource) ) then    
      ispec = ispec_selected_source(isource)      
      if( ispec_is_acoustic(ispec) ) then
        t0_ac = - 3.0d0 * ( t_cmt(isource) - hdur(isource) )
        if(  t0_ac > t0 ) t0 = t0_ac
      endif
    endif
  enddo

  ! passes maximum value to all processes
  ! note: t0 is defined positive and will be subtracted from simulation time (it-1)*DT
  t0_ac = t0
  call max_all_all_dp(t0_ac,t0)

  ! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there
  call setup_sources_check_acoustic()
  
end subroutine setup_sources

!
!-------------------------------------------------------------------------------------------------
!  

  
subroutine setup_sources_check_acoustic()

! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there

  use specfem_par
  use specfem_par_acoustic
  implicit none
  
  integer :: isource,ixmin,ixmax,iymin,iymax,izmin,izmax,iface,ispec
  logical :: is_on,is_on_all

! outputs a warning in case of an acoustic source lying on the free surface
  do isource = 1,NSOURCES
    ! checks if source is close to face 
    is_on = .false. 
  
    ! only receivers in this process  
    if( myrank == islice_selected_source(isource) ) then

      ispec = ispec_selected_source(isource)
      ! only if receiver is in an acoustic element
      if( ispec_is_acoustic(ispec) ) then
                  
        ! checks with free surface face
        do iface = 1,num_free_surface_faces
  
          if( ispec == free_surface_ispec(iface) ) then
          
            ! determine face 
            ixmin = minval( free_surface_ijk(1,:,iface) )
            ixmax = maxval( free_surface_ijk(1,:,iface) )
           
            iymin = minval( free_surface_ijk(2,:,iface) )
            iymax = maxval( free_surface_ijk(2,:,iface) )
           
            izmin = minval( free_surface_ijk(3,:,iface) )
            izmax = maxval( free_surface_ijk(3,:,iface) )
           
            if( .not. USE_FORCE_POINT_SOURCE ) then
              ! xmin face 
              if(ixmin==1 .and. ixmax==1) then
                if( xi_source(isource) < -0.99d0) is_on = .true.
              ! xmax face 
              else if(ixmin==NGLLX .and. ixmax==NGLLX) then
                if( xi_source(isource) > 0.99d0) is_on = .true.
              ! ymin face 
              else if(iymin==1 .and. iymax==1) then
                if( eta_source(isource) < -0.99d0) is_on = .true.
              ! ymax face 
              else if(iymin==NGLLY .and. iymax==NGLLY) then
                if( eta_source(isource) > 0.99d0) is_on = .true.
              ! zmin face 
              else if(izmin==1 .and. izmax==1 ) then
                if( gamma_source(isource) < -0.99d0) is_on = .true.
              ! zmax face 
              else if(izmin==NGLLZ .and. izmax==NGLLZ ) then
                if( gamma_source(isource) > 0.99d0) is_on = .true.
              endif
            else
              ! note: for use_force_point_source xi/eta/gamma_source values are in the range [1,NGLL*]            
              ! xmin face 
              if(ixmin==1 .and. ixmax==1) then
                if( nint(xi_source(isource)) == 1) is_on = .true.
              ! xmax face 
              else if(ixmin==NGLLX .and. ixmax==NGLLX) then
                if( nint(xi_source(isource)) == NGLLX) is_on = .true.
              ! ymin face 
              else if(iymin==1 .and. iymax==1) then
                if( nint(eta_source(isource)) == 1) is_on = .true.
              ! ymax face 
              else if(iymin==NGLLY .and. iymax==NGLLY) then
                if( nint(eta_source(isource)) == NGLLY) is_on = .true.
              ! zmin face 
              else if(izmin==1 .and. izmax==1 ) then
                if( nint(gamma_source(isource)) == 1) is_on = .true.
              ! zmax face 
              else if(izmin==NGLLZ .and. izmax==NGLLZ ) then
                if( nint(gamma_source(isource)) ==NGLLZ) is_on = .true.
              endif              
            endif
            
          endif ! free_surface_ispec
        enddo ! iface
      endif ! ispec_is_acoustic
    endif ! islice_selected_rec
    
    ! user output    
    call any_all_l( is_on, is_on_all )
    if( myrank == 0 .and. is_on_all ) then       
      write(IMAIN,*) '**********************************************************************'
      write(IMAIN,*) '*** source: ',isource,'                                          ***'
      write(IMAIN,*) '*** Warning: acoustic source located exactly on the free surface ***'
      write(IMAIN,*) '*** will be zeroed                                                                           ***'
      write(IMAIN,*) '**********************************************************************'
      write(IMAIN,*)
    endif    
    
  enddo ! num_free_surface_faces


end subroutine setup_sources_check_acoustic

!
!-------------------------------------------------------------------------------------------------
!  

subroutine setup_receivers()

  use specfem_par
  use specfem_par_acoustic
  implicit none
  
  integer :: irec,isource !,ios
  
! reads in station file  
  if (SIMULATION_TYPE == 1) then
    call get_value_string(rec_filename, 'solver.STATIONS', 'DATA/STATIONS')
    call get_value_string(filtered_rec_filename, 'solver.STATIONS_FILTERED', 'DATA/STATIONS_FILTERED')
    call station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE,myrank,rec_filename,filtered_rec_filename,nrec, &
           LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)

    ! get total number of stations
    !open(unit=IIN,file=rec_filename,iostat=ios,status='old',action='read')
    !nrec = 0
    !do while(ios == 0)
    !  read(IIN,"(a)",iostat=ios) dummystring
    !  if(ios == 0) nrec = nrec + 1
    !enddo
    !close(IIN)
    if(nrec < 1) call exit_MPI(myrank,'need at least one receiver')
    call sync_all()

  else
    call get_value_string(rec_filename, 'solver.STATIONS', 'DATA/STATIONS_ADJOINT')
    call get_value_string(filtered_rec_filename, 'solver.STATIONS_FILTERED', 'DATA/STATIONS_ADJOINT_FILTERED')
    call station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE,myrank,rec_filename,filtered_rec_filename,nrec, &
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
            xstore,ystore,zstore,xigll,yigll,zigll,filtered_rec_filename, &
            nrec,islice_selected_rec,ispec_selected_rec, &
            xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
            NPROC,utm_x_source(1),utm_y_source(1), &
            UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
            iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
            num_free_surface_faces,free_surface_ispec,free_surface_ijk)

! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    nrec_simulation = nrec
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    ! adjoint simulation: receivers become adjoint sources
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if(myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

! checks if acoustic receiver is exactly on the free surface because pressure is zero there
  call setup_receivers_check_acoustic()
  
end subroutine setup_receivers


!
!-------------------------------------------------------------------------------------------------
!  

subroutine setup_receivers_check_acoustic()

! checks if acoustic receiver is exactly on the free surface because pressure is zero there

  use specfem_par
  use specfem_par_acoustic
  implicit none
  
  integer :: irec,ixmin,ixmax,iymin,iymax,izmin,izmax,iface,ispec
  logical :: is_on,is_on_all

! outputs a warning in case the receiver is lying on the free surface
  do irec = 1,nrec
  
    ! checks if receiver is close to face 
    is_on = .false. 
  
    ! only receivers in this process  
    if( myrank == islice_selected_rec(irec) ) then

      ispec = ispec_selected_rec(irec)
      ! only if receiver is in an acoustic element
      if( ispec_is_acoustic(ispec) ) then
        
        ! checks with free surface face
        do iface = 1,num_free_surface_faces
  
          if( ispec == free_surface_ispec(iface) ) then
          
            ! determine face 
            ixmin = minval( free_surface_ijk(1,:,iface) )
            ixmax = maxval( free_surface_ijk(1,:,iface) )
           
            iymin = minval( free_surface_ijk(2,:,iface) )
            iymax = maxval( free_surface_ijk(2,:,iface) )
           
            izmin = minval( free_surface_ijk(3,:,iface) )
            izmax = maxval( free_surface_ijk(3,:,iface) )
           
            ! xmin face 
            if(ixmin==1 .and. ixmax==1) then
              if( xi_receiver(irec) < -0.99d0) is_on = .true.
            ! xmax face 
            else if(ixmin==NGLLX .and. ixmax==NGLLX) then
              if( xi_receiver(irec) > 0.99d0) is_on = .true.
            ! ymin face 
            else if(iymin==1 .and. iymax==1) then
              if( eta_receiver(irec) < -0.99d0) is_on = .true.
            ! ymax face 
            else if(iymin==NGLLY .and. iymax==NGLLY) then
              if( eta_receiver(irec) > 0.99d0) is_on = .true.
            ! zmin face 
            else if(izmin==1 .and. izmax==1 ) then
              if( gamma_receiver(irec) < -0.99d0) is_on = .true.
            ! zmax face 
            else if(izmin==NGLLZ .and. izmax==NGLLZ ) then
              if( gamma_receiver(irec) > 0.99d0) is_on = .true.
            endif
                
          endif ! free_surface_ispec
        enddo ! iface
      endif ! ispec_is_acoustic
    endif ! islice_selected_rec
    
    ! user output    
    call any_all_l( is_on, is_on_all )
    if( myrank == 0 .and. is_on_all ) then       
      write(IMAIN,*) '**********************************************************************'
      write(IMAIN,*) '*** station:',irec,'                                          ***'
      write(IMAIN,*) '*** Warning: acoustic receiver located exactly on the free surface ***'
      write(IMAIN,*) '*** Warning: tangential component will be zero there               ***'
      write(IMAIN,*) '**********************************************************************'
      write(IMAIN,*)
    endif
        
  enddo ! num_free_surface_faces

end subroutine setup_receivers_check_acoustic


!
!-------------------------------------------------------------------------------------------------
!  
  
subroutine setup_sources_precompute_arrays()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  implicit none
  
  real(kind=CUSTOM_REAL) :: factor_source
  real(kind=CUSTOM_REAL) :: junk
  integer :: isource,ispec
  integer :: irec,irec_local
  integer :: icomp,itime,nadj_files_found,nadj_files_found_tot,ier
  character(len=3),dimension(NDIM) :: comp = (/ "BHN", "BHE", "BHZ" /)
  character(len=150) :: filename

  
! forward simulations  
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
      
        ispec = ispec_selected_source(isource)
        
        ! elastic moment tensor source
        if( ispec_is_elastic(ispec) ) then
          call compute_arrays_source(ispec, &
                        xi_source(isource),eta_source(isource),gamma_source(isource),sourcearray, &
                        Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource), &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        xigll,yigll,zigll,NSPEC_AB)
        endif
        
        ! acoustic case 
        if( ispec_is_acoustic(ispec) ) then
          ! scalar moment of moment tensor values read in from CMTSOLUTION 
          ! note: M0 by Dahlen and Tromp, eq. 5.91
          factor_source = 1.0/sqrt(2.0) * sqrt( Mxx(isource)**2 + Myy(isource)**2 + Mzz(isource)**2 &
                                    + 2*( Myz(isource)**2 + Mxz(isource)**2 + Mxy(isource)**2 ) )

          ! scales source such that it would be equivalent to explosion source moment tensor,
          ! where Mxx=Myy=Mzz, others Mxy,.. = zero, in equivalent elastic media
          ! (and getting rid of 1/sqrt(2) factor from scalar moment tensor definition above)
          factor_source = factor_source * sqrt(2.0) / sqrt(3.0)

          ! source array interpolated on all element gll points
          call compute_arrays_source_acoustic(xi_source(isource),eta_source(isource),gamma_source(isource),&
                        sourcearray,xigll,yigll,zigll,factor_source)
        endif
        
        ! stores source excitations
        sourcearrays(isource,:,:,:,:) = sourcearray(:,:,:,:)
        
      endif
    enddo
  endif

! ADJOINT simulations  
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
  
    ! counts local receivers which become adjoint sources
    nadj_rec_local = 0
    ! temporary counter to check if any files are found at all
    nadj_files_found = 0    
    do irec = 1,nrec
      if( myrank == islice_selected_rec(irec) ) then
        ! checks that the source slice number is okay
        if(islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROC-1) &
              call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')
              
        ! updates counter
        nadj_rec_local = nadj_rec_local + 1

        ! checks **sta**.**net**.**BH**.adj files for correct number of time steps
        adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
        do icomp = 1,NDIM
          filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
          open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
          if( ier == 0 ) then
            ! checks length of file
            itime = 0
            do while(ier == 0) 
              read(IIN,*,iostat=ier) junk,junk
              if( ier == 0 ) itime = itime + 1
            enddo
            if( itime /= NSTEP) &
              call exit_MPI(myrank,&
                'file '//trim(filename)//' has wrong length, please check with your simulation duration')
            nadj_files_found = nadj_files_found + 1
          endif
          close(IIN)
        enddo        
      endif
    enddo
    ! checks if any adjoint source files found at all
    call sum_all_i(nadj_files_found,nadj_files_found_tot)
    if( myrank == 0 ) then
      write(IMAIN,*)
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component traces found in all slices'
      if(nadj_files_found_tot == 0) &
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
    endif

    ! reads in adjoint source traces
    allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLY,NGLLZ))
    allocate(adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ))
    adj_sourcearrays = 0._CUSTOM_REAL
    adj_sourcearray = 0._CUSTOM_REAL

    ! pre-computes adjoint source arrays
    irec_local = 0
    do irec = 1, nrec
      ! computes only adjoint source arrays in the local slice
      if( myrank == islice_selected_rec(irec) ) then
        irec_local = irec_local + 1

        ! reads in **sta**.**net**.**BH**.adj files        
        adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
        
        call compute_arrays_adjoint_source(myrank, adj_source_file, &
                          xi_receiver(irec), eta_receiver(irec), gamma_receiver(irec), &
                          adj_sourcearray, xigll,yigll,zigll,NSTEP)

        adj_sourcearrays(irec_local,:,:,:,:,:) = adj_sourcearray(:,:,:,:,:)

      endif
    enddo
    ! frees temporary array
    deallocate(adj_sourcearray)
  endif

end subroutine setup_sources_precompute_arrays


!
!-------------------------------------------------------------------------------------------------
!  

subroutine setup_receivers_precompute_intp()

  use specfem_par
  implicit none
  
  integer :: irec,irec_local,isource
  
! stores local receivers interpolation factors
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
  if( myrank == 0 ) then
    if(nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    endif
  endif
  

end subroutine setup_receivers_precompute_intp
!
!-------------------------------------------------------------------------------------------------
!  

subroutine setup_sources_receivers_VTKfile()

  use specfem_par
  implicit none

  double precision :: shape3D(NGNOD)  
  double precision :: xil,etal,gammal
  double precision :: xmesh,ymesh,zmesh
  
  real(kind=CUSTOM_REAL),dimension(NGNOD) :: xelm,yelm,zelm  
  
  integer :: ia,ispec,isource,irec
  
  if (myrank == 0) then
    ! vtk file
    open(IOVTK,file=trim(OUTPUT_FILES)//'/sr.vtk',status='unknown')
    write(IOVTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOVTK,'(a)') 'Source and Receiver VTK file'
    write(IOVTK,'(a)') 'ASCII'
    write(IOVTK,'(a)') 'DATASET POLYDATA'
    write(IOVTK, '(a,i6,a)') 'POINTS ', NSOURCES+nrec, ' float'
  endif
  
  ! sources
  do isource=1,NSOURCES    
    ! spectral element id
    ispec = ispec_selected_source(isource)
    
    ! gets element ancor nodes
    if( myrank == islice_selected_source(isource) ) then
      ! find the coordinates of the eight corner nodes of the element
      call get_shape3D_element_corners(xelm,yelm,zelm,ispec,&
                      ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)

    endif
    ! master collects corner locations
    if( islice_selected_source(isource) /= 0 ) then
      if( myrank == 0 ) then
        call recvv_cr(xelm,NGNOD,islice_selected_source(isource),0)
        call recvv_cr(yelm,NGNOD,islice_selected_source(isource),0)
        call recvv_cr(zelm,NGNOD,islice_selected_source(isource),0)
      else if( myrank == islice_selected_source(isource) ) then
        call sendv_cr(xelm,NGNOD,0,0)
        call sendv_cr(yelm,NGNOD,0,0)
        call sendv_cr(zelm,NGNOD,0,0)
      endif
    endif
    
    if( myrank == 0 ) then
      ! get the 3-D shape functions
      xil = xi_source(isource)
      etal = eta_source(isource)
      gammal = gamma_source(isource)
      call get_shape3D_single(myrank,shape3D,xil,etal,gammal)            

      ! interpolates source locations
      xmesh = 0.0
      ymesh = 0.0
      zmesh = 0.0      
      do ia=1,NGNOD
        xmesh = xmesh + shape3D(ia)*xelm(ia)
        ymesh = ymesh + shape3D(ia)*yelm(ia)
        zmesh = zmesh + shape3D(ia)*zelm(ia)
      enddo

      ! writes out to VTK file
      write(IOVTK,*) xmesh,ymesh,zmesh
    endif
  enddo ! NSOURCES

  ! receivers
  do irec=1,nrec
    ispec = ispec_selected_rec(irec)
          
    ! find the coordinates of the eight corner nodes of the element
    if( myrank == islice_selected_rec(irec) ) then
      call get_shape3D_element_corners(xelm,yelm,zelm,ispec,&
                      ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)
    endif
    ! master collects corner locations
    if( islice_selected_rec(irec) /= 0 ) then
      if( myrank == 0 ) then
        call recvv_cr(xelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(yelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(zelm,NGNOD,islice_selected_rec(irec),0)
      else if( myrank == islice_selected_rec(irec) ) then
        call sendv_cr(xelm,NGNOD,0,0)
        call sendv_cr(yelm,NGNOD,0,0)
        call sendv_cr(zelm,NGNOD,0,0)
      endif
    endif

    if( myrank == 0 ) then
      ! get the 3-D shape functions
      if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then      
        xil = xi_receiver(irec)
        etal = eta_receiver(irec)
        gammal = gamma_receiver(irec)
      else
        xil = xi_source(irec)
        etal = eta_source(irec)
        gammal = gamma_source(irec)      
      endif
      call get_shape3D_single(myrank,shape3D,xil,etal,gammal)            
      
      ! interpolates receiver locations        
      xmesh = 0.0
      ymesh = 0.0
      zmesh = 0.0      
      do ia=1,NGNOD
        xmesh = xmesh + shape3D(ia)*xelm(ia)
        ymesh = ymesh + shape3D(ia)*yelm(ia)
        zmesh = zmesh + shape3D(ia)*zelm(ia)
      enddo

      ! writes out to VTK file
      write(IOVTK,*) xmesh,ymesh,zmesh      
    endif
  enddo
  
  ! closes vtk file
  if( myrank == 0 ) then
    write(IOVTK,*)
    close(IOVTK)
  endif

end subroutine setup_sources_receivers_VTKfile
