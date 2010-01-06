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

  subroutine iterate_time()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  
  implicit none

!
!   s t a r t   t i m e   i t e r a t i o n s
!

! synchronize all processes to make sure everybody is ready to start time loop
  call sync_all()
  if(myrank == 0) write(IMAIN,*) 'All processes are synchronized before time loop'

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
  endif

! create an empty file to monitor the start of the simulation
  if(myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown')
    write(IOUT,*) 'starting time loop'
    close(IOUT)
  endif

! get MPI starting time
  time_start = wtime()

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  do it = 1,NSTEP
  
! simulation status output and stability check
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5) then
      call it_check_stability()    
    endif
    
! update displacement using Newark time scheme
    call it_update_displacement_scheme()

! acoustic solver 
! (needs to be done first, before elastic one)
    if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic()
      
! elastic solver
    if( ELASTIC_SIMULATION ) call compute_forces_elastic()
    
! poroelastic solver
    if( POROELASTIC_SIMULATION ) stop 'poroelastic simulation not implemented yet'
    
! write the seismograms with time shift
    if (nrec_local > 0) then
      call it_write_seismograms()
    endif 

! resetting d/v/a/R/eps for the backward reconstruction with attenuation
    if (ATTENUATION ) then
      call it_store_attenuation_arrays()
    endif ! ATTENUATION

! shakemap creation
    if (EXTERNAL_MESH_CREATE_SHAKEMAP) then
      call it_create_shakemap_em()
    endif 

! movie file creation
    if(EXTERNAL_MESH_MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
      call it_create_movie_surface_em()
    endif

! save MOVIE on the SURFACE
    if(MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

      !stop 'DK DK MOVIE_SURFACE has been removed for now because we need a flag to detect the surface elements'

      call it_movie_surface_output_o()
    endif

! compute SHAKING INTENSITY MAP
    if(CREATE_SHAKEMAP) then

      !stop 'DK DK CREATE_SHAKEMAP has been removed for now because we need a flag to detect the surface elements'

      call it_create_shakemap_o()
    endif

! save MOVIE in full 3D MESH
    if(MOVIE_VOLUME .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
      call it_movie_volume_output()
    endif

! creates cross-section GIF image
    if(PNM_GIF_IMAGE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0 ) then
      call write_PNM_GIF_create_image()
    endif
!
!---- end of time iteration loop
!
  enddo   ! end of main time loop

  end subroutine iterate_time

  
!=====================================================================

  subroutine it_check_stability()

! computes the maximum of the norm of the displacement
! in all the slices using an MPI reduction
! and output timestamp file to check that simulation is running fine
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic  
  implicit none
  
  double precision :: tCPU,t_remain,t_total
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total
  
! compute maximum of norm of displacement in each slice
  if( ELASTIC_SIMULATION ) then
    Usolidnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2 + displ(3,:)**2))
  else 
    if( ACOUSTIC_SIMULATION ) then
      Usolidnorm = maxval(abs(potential_dot_dot_acoustic(:)))
    endif
  endif  
  
! compute the maximum of the maxima for all the slices using an MPI reduction
  call max_all_cr(Usolidnorm,Usolidnorm_all)

!! DK DK array not created yet for CUBIT
!   if (SIMULATION_TYPE == 3) then
!     b_Usolidnorm = maxval(sqrt(b_displ(1,:)**2 + b_displ(2,:)**2 + b_displ(3,:)**2))
!     call max_all_cr(b_Usolidnorm,b_Usolidnorm_all)
!   endif

  if(myrank == 0) then

    write(IMAIN,*) 'Time step # ',it
    write(IMAIN,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'

! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    if( ELASTIC_SIMULATION ) then
      write(IMAIN,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
    else 
      if( ACOUSTIC_SIMULATION ) then
        write(IMAIN,*) 'Max norm pressure P in all slices (m) = ',Usolidnorm_all    
      endif
    endif
!     if (SIMULATION_TYPE == 3) write(IMAIN,*) &
!           'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all

! compute estimated remaining simulation time
    t_remain = (NSTEP - it) * (tCPU/dble(it))
    int_t_remain = int(t_remain)
    ihours_remain = int_t_remain / 3600
    iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
    iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain
    write(IMAIN,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IMAIN,*) 'Time steps remaining = ',NSTEP - it
    write(IMAIN,*) 'Estimated remaining time in seconds = ',t_remain
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',t_total
    write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'

    if(it < 100) then
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '**** BEWARE: the above time estimates are not reliable'
      write(IMAIN,*) '**** because fewer than 100 iterations have been performed'
      write(IMAIN,*) '************************************************************'
    endif
    write(IMAIN,*)

! write time stamp file to give information about progression of simulation
    write(outputname,"('/timestamp',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'
    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    write(IOUT,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
!     if (SIMULATION_TYPE == 3) write(IOUT,*) &
!           'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
    close(IOUT)


! check stability of the code, exit if unstable
! negative values can occur with some compilers when the unstable value is greater
! than the greatest possible floating-point number of the machine
    if(Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0) &
        call exit_MPI(myrank,'forward simulation became unstable and blew up')
!     if(SIMULATION_TYPE == 3 .and. (b_Usolidnorm_all > STABILITY_THRESHOLD .or. b_Usolidnorm_all < 0)) &
!       call exit_MPI(myrank,'backward simulation became unstable and blew up')

  endif ! myrank
  
  end subroutine it_check_stability
  

!=====================================================================

  subroutine it_update_displacement_scheme()

! explicit Newark time scheme with acoustic & elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 delta_t chi_dot_dot(t+delta_t)
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where 
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!
! note that this stage calculates the predictor terms
!
!   for 
!   potential chi_dot(t+delta) requires + 1/2 delta_t chi_dot_dot(t+delta_t)
!                                   at a later stage (corrector) once where chi_dot_dot(t+delta) is calculated
!   and similar,
!   velocity v(t+delta_t) requires  + 1/2 delta_t a(t+delta_t)  
!                                   at a later stage once where a(t+delta) is calculated
! also:
!   boundary term B_elastic requires chi_dot_dot(t+delta)
!                                   thus chi_dot_dot has to be updated first before the elastic boundary term is considered
  
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  
  implicit none

! updates acoustic potentials
  if( ACOUSTIC_SIMULATION ) then
    potential_acoustic(:) = potential_acoustic(:) &
                            + deltat * potential_dot_acoustic(:) &
                            + deltatsqover2 * potential_dot_dot_acoustic(:)
    potential_dot_acoustic(:) = potential_dot_acoustic(:) &
                                + deltatover2 * potential_dot_dot_acoustic(:)
    potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
  endif

! updates elastic displacement and velocity
  if( ELASTIC_SIMULATION ) then
    displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
    accel(:,:) = 0._CUSTOM_REAL
    
    !! DK DK array not created yet for CUBIT
    ! if (SIMULATION_TYPE == 3) then
    !   b_displ(:,:) = b_displ(:,:) + b_deltat*b_veloc(:,:) + b_deltatsqover2*b_accel(:,:)
    !   b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)
    !   b_accel(:,:) = 0._CUSTOM_REAL
    ! endif
  endif


  end subroutine it_update_displacement_scheme
  
!=====================================================================

  subroutine it_write_seismograms()

! writes the seismograms with time shift
  
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic  
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element,veloc_element
  double precision :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,hlagrange
  integer :: irec_local,irec
  integer :: iglob,ispec,i,j,k

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

! perform the general interpolation using Lagrange polynomials
    if(FASTER_RECEIVERS_POINTS_ONLY) then
      ispec = ispec_selected_rec(irec)
      iglob = ibool(nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
         nint(gamma_receiver(irec)),ispec)

      ! elastic wave field   
      if( ispec_is_elastic(ispec) ) then
        dxd = dble(displ(1,iglob))
        dyd = dble(displ(2,iglob))
        dzd = dble(displ(3,iglob))
        vxd = dble(veloc(1,iglob))
        vyd = dble(veloc(2,iglob))
        vzd = dble(veloc(3,iglob))
        axd = dble(accel(1,iglob))
        ayd = dble(accel(2,iglob))
        azd = dble(accel(3,iglob))
      endif
      
      ! acoustic wave field
      if( ispec_is_acoustic(ispec) ) then
        ! displacement
        call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic, displ_element,&
                        hprime_xx,hprime_yy,hprime_zz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        ibool,rhostore)
        ! velocity
        call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                        potential_dot_acoustic, veloc_element,&
                        hprime_xx,hprime_yy,hprime_zz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        ibool,rhostore)
        ! displacement
        dxd = displ_element(1,nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
                                                    nint(gamma_receiver(irec)))
        dyd = displ_element(2,nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
                                                    nint(gamma_receiver(irec)))
        dzd = displ_element(3,nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
                                                    nint(gamma_receiver(irec)))
        ! velocity
        vxd = veloc_element(1,nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
                                                    nint(gamma_receiver(irec)))
        vyd = veloc_element(2,nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
                                                    nint(gamma_receiver(irec)))
        vzd = veloc_element(3,nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
                                                    nint(gamma_receiver(irec)))
        ! pressure
        axd = - potential_dot_dot_acoustic(iglob)
        ayd = - potential_dot_dot_acoustic(iglob)
        azd = - potential_dot_dot_acoustic(iglob)                                          
      endif ! acoustic
      
    else

      dxd = ZERO
      dyd = ZERO
      dzd = ZERO

      vxd = ZERO
      vyd = ZERO
      vzd = ZERO

      axd = ZERO
      ayd = ZERO
      azd = ZERO

      if (SIMULATION_TYPE == 1)  then

        ispec = ispec_selected_rec(irec)

        ! elastic wave field    
        if( ispec_is_elastic(ispec) ) then
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                
                ! receivers are always located at the surface of the mesh
                iglob = ibool(i,j,k,ispec)

                hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

                ! elastic wave field
                if( ispec_is_elastic(ispec) ) then
                  ! save displacement
                  dxd = dxd + dble(displ(1,iglob))*hlagrange
                  dyd = dyd + dble(displ(2,iglob))*hlagrange
                  dzd = dzd + dble(displ(3,iglob))*hlagrange
                  ! save velocity
                  vxd = vxd + dble(veloc(1,iglob))*hlagrange
                  vyd = vyd + dble(veloc(2,iglob))*hlagrange
                  vzd = vzd + dble(veloc(3,iglob))*hlagrange
                  ! save acceleration
                  axd = axd + dble(accel(1,iglob))*hlagrange
                  ayd = ayd + dble(accel(2,iglob))*hlagrange
                  azd = azd + dble(accel(3,iglob))*hlagrange
                endif
                
              enddo
            enddo
          enddo
        endif
        
        ! acoustic wave field
        if( ispec_is_acoustic(ispec) ) then
          ! displacement vector
          call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_acoustic, displ_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
          ! velocity vector
          call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, veloc_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
          ! interpolates vector field                 
          do k= 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)                          
                hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)              
                ! displacement
                dxd = dxd + hlagrange*displ_element(1,i,j,k)
                dyd = dyd + hlagrange*displ_element(2,i,j,k)
                dzd = dzd + hlagrange*displ_element(3,i,j,k)
                ! velocity
                vxd = vxd + hlagrange*veloc_element(1,i,j,k)
                vyd = vxd + hlagrange*veloc_element(2,i,j,k)
                vzd = vxd + hlagrange*veloc_element(3,i,j,k)
                ! pressure
                axd = axd - hlagrange*potential_dot_dot_acoustic(iglob)
                ayd = ayd - hlagrange*potential_dot_dot_acoustic(iglob)
                azd = azd - hlagrange*potential_dot_dot_acoustic(iglob)                  
              enddo
            enddo
          enddo
        endif ! acoustic
        
      else if (SIMULATION_TYPE == 2) then

        ! adjoint source is placed at receiver
        ispec = ispec_selected_source(irec)

        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX

              iglob = ibool(i,j,k,ispec)

              hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

              dxd = dxd + dble(displ(1,iglob))*hlagrange
              dyd = dyd + dble(displ(2,iglob))*hlagrange
              dzd = dzd + dble(displ(3,iglob))*hlagrange
              vxd = vxd + dble(veloc(1,iglob))*hlagrange
              vyd = vyd + dble(veloc(2,iglob))*hlagrange
              vzd = vzd + dble(veloc(3,iglob))*hlagrange
              axd = axd + dble(accel(1,iglob))*hlagrange
              ayd = ayd + dble(accel(2,iglob))*hlagrange
              azd = azd + dble(accel(3,iglob))*hlagrange

              displ_s(:,i,j,k) = displ(:,iglob)

            enddo
          enddo
        enddo

        !ispec = ispec_selected_source(irec)

        call compute_adj_source_frechet(displ_s,Mxx(irec),Myy(irec),Mzz(irec),Mxy(irec),Mxz(irec),Myz(irec),eps_s,eps_m_s, &
             hxir_store(irec_local,:),hetar_store(irec_local,:),hgammar_store(irec_local,:), &
             hpxir_store(irec_local,:),hpetar_store(irec_local,:),hpgammar_store(irec_local,:),hprime_xx,hprime_yy,hprime_zz, &
             xix(:,:,:,ispec),xiy(:,:,:,ispec),xiz(:,:,:,ispec),etax(:,:,:,ispec),etay(:,:,:,ispec),etaz(:,:,:,ispec), &
             gammax(:,:,:,ispec),gammay(:,:,:,ispec),gammaz(:,:,:,ispec))

        stf = comp_source_time_function(dble(NSTEP-it)*DT-t0-t_cmt(irec),hdur_gaussian(irec))
        stf_deltat = stf * deltat
        Mxx_der(irec_local) = Mxx_der(irec_local) + eps_s(1,1) * stf_deltat
        Myy_der(irec_local) = Myy_der(irec_local) + eps_s(2,2) * stf_deltat
        Mzz_der(irec_local) = Mzz_der(irec_local) + eps_s(3,3) * stf_deltat
        Mxy_der(irec_local) = Mxy_der(irec_local) + 2 * eps_s(1,2) * stf_deltat
        Mxz_der(irec_local) = Mxz_der(irec_local) + 2 * eps_s(1,3) * stf_deltat
        Myz_der(irec_local) = Myz_der(irec_local) + 2 * eps_s(2,3) * stf_deltat

        sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_s(:) * stf_deltat

      else if (SIMULATION_TYPE == 3) then
        
        ispec = ispec_selected_rec(irec)
        
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX

              iglob = ibool(i,j,k,ispec)

              hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

!! DK DK array not created yet for CUBIT
!         dxd = dxd + dble(b_displ(1,iglob))*hlagrange
!         dyd = dyd + dble(b_displ(2,iglob))*hlagrange
!         dzd = dzd + dble(b_displ(3,iglob))*hlagrange
!         vxd = vxd + dble(b_veloc(1,iglob))*hlagrange
!         vyd = vyd + dble(b_veloc(2,iglob))*hlagrange
!         vzd = vzd + dble(b_veloc(3,iglob))*hlagrange
!         axd = axd + dble(b_accel(1,iglob))*hlagrange
!         ayd = ayd + dble(b_accel(2,iglob))*hlagrange
!         azd = azd + dble(b_accel(3,iglob))*hlagrange
            enddo
          enddo
        enddo
      endif ! SIMULATION_TYPE

    endif ! FASTER_RECEIVERS_POINTS_ONLY

! store North, East and Vertical components

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      seismograms_d(:,irec_local,it) = sngl((nu(:,1,irec)*dxd + nu(:,2,irec)*dyd + nu(:,3,irec)*dzd))
      seismograms_v(:,irec_local,it) = sngl((nu(:,1,irec)*vxd + nu(:,2,irec)*vyd + nu(:,3,irec)*vzd))
      seismograms_a(:,irec_local,it) = sngl((nu(:,1,irec)*axd + nu(:,2,irec)*ayd + nu(:,3,irec)*azd))
    else
      seismograms_d(:,irec_local,it) = (nu(:,1,irec)*dxd + nu(:,2,irec)*dyd + nu(:,3,irec)*dzd)
      seismograms_v(:,irec_local,it) = (nu(:,1,irec)*vxd + nu(:,2,irec)*vyd + nu(:,3,irec)*vzd)
      seismograms_a(:,irec_local,it) = (nu(:,1,irec)*axd + nu(:,2,irec)*ayd + nu(:,3,irec)*azd)
    endif

    if (SIMULATION_TYPE == 2) seismograms_eps(:,:,irec_local,it) = eps_s(:,:)

  enddo ! nrec_local

! write the current or final seismograms
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      call write_seismograms(myrank,seismograms_d,number_receiver_global,station_name, &
            network_name,nrec,nrec_local,it,DT,NSTEP,t0,LOCAL_PATH,1)
      call write_seismograms(myrank,seismograms_v,number_receiver_global,station_name, &
            network_name,nrec,nrec_local,it,DT,NSTEP,t0,LOCAL_PATH,2)
      call write_seismograms(myrank,seismograms_a,number_receiver_global,station_name, &
            network_name,nrec,nrec_local,it,DT,NSTEP,t0,LOCAL_PATH,3)
    else
      call write_adj_seismograms(myrank,seismograms_d,number_receiver_global, &
            nrec_local,it,DT,NSTEP,t0,LOCAL_PATH,1)
    endif
  endif

  end subroutine it_write_seismograms


!================================================================
  
  subroutine it_store_attenuation_arrays()

! resetting d/v/a/R/eps for the backward reconstruction with attenuation
  
  use specfem_par
  use specfem_par_elastic
  
  implicit none

  if( it > 1 .and. it < NSTEP) then
    if (SIMULATION_TYPE == 3 .and. mod(NSTEP-it,NSTEP_Q_SAVE) == 0) then
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") NSTEP-it
      open(unit=27,file=trim(prname_Q)//trim(outputname),status='old',action='read',form='unformatted')
!! DK DK array not created yet for CUBIT
!   read(27) b_displ
!   read(27) b_veloc
!   read(27) b_accel
!   read(27) b_R_xx
!   read(27) b_R_yy
!   read(27) b_R_xy
!   read(27) b_R_xz
!   read(27) b_R_yz
!   read(27) b_epsilondev_xx
!   read(27) b_epsilondev_yy
!   read(27) b_epsilondev_xy
!   read(27) b_epsilondev_xz
!   read(27) b_epsilondev_yz
      close(27)
    else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. mod(it,NSTEP_Q_SAVE) == 0) then
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") it
      open(unit=27,file=trim(prname_Q)//trim(outputname),status='unknown',action='write',form='unformatted')
      write(27) displ
      write(27) veloc
      write(27) accel
      write(27) R_xx
      write(27) R_yy
      write(27) R_xy
      write(27) R_xz
      write(27) R_yz
      write(27) epsilondev_xx
      write(27) epsilondev_yy
      write(27) epsilondev_xy
      write(27) epsilondev_xz
      write(27) epsilondev_yz
      close(27)
    endif ! SIMULATION_TYPE
  endif ! it

  end subroutine it_store_attenuation_arrays
  
!================================================================
  
  subroutine it_create_shakemap_em()

! creation of shapemap file
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_movie
  implicit none
  
  integer :: ipoin,ispec,iglob,ispec2D

! initializes arrays for point coordinates
  if (it == 1) then
    store_val_ux_external_mesh(:) = -HUGEVAL
    store_val_uy_external_mesh(:) = -HUGEVAL
    store_val_uz_external_mesh(:) = -HUGEVAL
    do ispec2D = 1,nfaces_surface_ext_mesh
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, NGLLX*NGLLY
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = zstore(iglob)
        enddo
      else
        do ipoin = 1, 4
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = zstore(iglob)        
        enddo
      endif
    enddo
  endif

! stores displacement, velocity and acceleration amplitudes
  do ispec2D = 1,nfaces_surface_ext_mesh
    ispec = faces_surface_ext_mesh_ispec(ispec2D)    
    ! high-resolution
    if (USE_HIGHRES_FOR_MOVIES) then
      do ipoin = 1, NGLLX*NGLLY
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! saves norm of displacement,velocity and acceleration vector
        if( ispec_is_elastic(ispec) ) then            
          ! norm of displacement
          store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
               max(store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
               sqrt(displ(1,iglob)**2 + displ(2,iglob)**2 + displ(3,iglob)**2))
          ! norm of velocity     
          store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
               max(store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
               sqrt(veloc(1,iglob)**2 + veloc(2,iglob)**2 + veloc(3,iglob)**2))
          ! norm of acceleration     
          store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
               max(store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
               sqrt(accel(1,iglob)**2 + accel(2,iglob)**2 + accel(3,iglob)**2))
        endif
      enddo
    else
      ! low-resolution: only corner points outputted
      do ipoin = 1, 4
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! saves norm of displacement,velocity and acceleration vector
        if( ispec_is_elastic(ispec) ) then                    
          ! norm of displacement
          store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                max(store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                sqrt(displ(1,iglob)**2 + displ(2,iglob)**2 + displ(3,iglob)**2))
          ! norm of velocity      
          store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                max(store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                sqrt(veloc(1,iglob)**2 + veloc(2,iglob)**2 + veloc(3,iglob)**2))
          ! norm of acceleration
          store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                max(store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                sqrt(accel(1,iglob)**2 + accel(2,iglob)**2 + accel(3,iglob)**2))
        endif
      enddo
    endif
  enddo

! finalizes shakemap: master process collects all info   
  if (it == NSTEP) then
    if (USE_HIGHRES_FOR_MOVIES) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    else
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif

! creates shakemap file
    if(myrank == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata',status='unknown',form='unformatted')
      write(IOUT) store_val_x_all_external_mesh   ! x coordinates
      write(IOUT) store_val_y_all_external_mesh   ! y coordinates
      write(IOUT) store_val_z_all_external_mesh   ! z coordinates
      write(IOUT) store_val_ux_all_external_mesh  ! norm of displacement vector
      write(IOUT) store_val_uy_all_external_mesh  ! norm of velocity vector
      write(IOUT) store_val_uz_all_external_mesh  ! norm of acceleration vector
      close(IOUT)
    endif
  endif
  
  end subroutine it_create_shakemap_em
  
  
!================================================================

  subroutine it_create_movie_surface_em()

! creation of moviedata files  

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie  
  implicit none
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: veloc_element
  integer :: ispec2D,ispec,ipoin,iglob,i,j,k
  logical :: is_done
  
! initializes arrays for point coordinates
  if (it == NTSTEP_BETWEEN_FRAMES ) then
    do ispec2D = 1,nfaces_surface_ext_mesh
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, NGLLX*NGLLY
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = zstore(iglob)
        enddo
      else
        do ipoin = 1, 4
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = zstore(iglob)                  
        enddo
      endif
    enddo
  endif
  
! saves surface velocities
  do ispec2D = 1,nfaces_surface_ext_mesh
    ispec = faces_surface_ext_mesh_ispec(ispec2D)      

    if( ispec_is_acoustic(ispec) ) then
      ! velocity vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, veloc_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
    endif
    
    if (USE_HIGHRES_FOR_MOVIES) then
      do ipoin = 1, NGLLX*NGLLY
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! x,y,z coordinates
        !store_val_x_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = xstore(iglob)
        !store_val_y_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = ystore(iglob)
        !store_val_z_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = zstore(iglob)
        ! saves velocity vector        
        if( ispec_is_elastic(ispec) ) then
          ! velocity x,y,z-components
          store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc(1,iglob)
          store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc(2,iglob)
          store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc(3,iglob)
        endif
        ! acoustic pressure potential
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc_element(1,i,j,k)
                  store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc_element(2,i,j,k)
                  store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc_element(3,i,j,k)
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
          ! only pressure
          !store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = -potential_dot_dot_acoustic(iglob)
          !store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = -potential_dot_dot_acoustic(iglob)
          !store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = -potential_dot_dot_acoustic(iglob)        
        endif
      enddo
    else
      do ipoin = 1, 4
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! x,y,z coordinates
        !store_val_x_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = xstore(iglob)
        !store_val_y_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = ystore(iglob)
        !store_val_z_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = zstore(iglob)
        ! saves velocity vector        
        if( ispec_is_elastic(ispec) ) then
          ! velocity x,y,z-components
          store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc(1,iglob)
          store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc(2,iglob)
          store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc(3,iglob)      
        endif
        ! acoustic pressure potential
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc_element(1,i,j,k)
                  store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc_element(2,i,j,k)
                  store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc_element(3,i,j,k)
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
          ! only pressure
          !store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = -potential_dot_dot_acoustic(iglob)
          !store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = -potential_dot_dot_acoustic(iglob)
          !store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = -potential_dot_dot_acoustic(iglob)                
        endif
      enddo
    endif
  enddo

! master process collects all info
  if (USE_HIGHRES_FOR_MOVIES) then
    ! collects locations only once
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    endif
    ! updates/gathers velocity field (high-res)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
  else
    ! collects locations only once
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif
    ! updates/gathers velocity field (low-res)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
  endif

! file output
  if(myrank == 0) then
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
    write(IOUT) store_val_x_all_external_mesh   ! x coordinate
    write(IOUT) store_val_y_all_external_mesh   ! y coordinate
    write(IOUT) store_val_z_all_external_mesh   ! z coordinate
    write(IOUT) store_val_ux_all_external_mesh  ! velocity x-component
    write(IOUT) store_val_uy_all_external_mesh  ! velocity y-component
    write(IOUT) store_val_uz_all_external_mesh  ! velocity z-component
    close(IOUT)
  endif
  
  end subroutine it_create_movie_surface_em

    
!=====================================================================

  subroutine it_movie_surface_output_o()

! outputs moviedata files  
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_movie
  
  implicit none
  integer :: imin,imax,jmin,jmax,kmin,kmax,iface,igll
  integer :: ipoin,iloc
  integer :: ispec,i,j,k,iglob

! initializes arrays for point coordinates
  if (it == NTSTEP_BETWEEN_FRAMES ) then
    ipoin = 0
    do iface=1,num_free_surface_faces
      ispec = free_surface_ispec(iface)
      ! high_resolution
      if (USE_HIGHRES_FOR_MOVIES) then      
        do igll = 1, NGLLSQUARE
          ipoin = ipoin + 1
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)      
          iglob = ibool(i,j,k,ispec)
          ! coordinates
          store_val_x_external_mesh(ipoin) = xstore(iglob)
          store_val_y_external_mesh(ipoin) = ystore(iglob)
          store_val_z_external_mesh(ipoin) = zstore(iglob)
        enddo
      else
        imin = minval( free_surface_ijk(1,:,iface) )
        imax = maxval( free_surface_ijk(1,:,iface) )
        jmin = minval( free_surface_ijk(2,:,iface) )
        jmax = maxval( free_surface_ijk(2,:,iface) )
        kmin = minval( free_surface_ijk(3,:,iface) )
        kmax = maxval( free_surface_ijk(3,:,iface) )      
        do iloc = 1, NGNOD2D    
          ipoin = ipoin + 1
          ! corner points
          if( imin == imax ) then
            iglob = ibool(imin,iorderi(iloc),iorderj(iloc),ispec)
          else if( jmin == jmax ) then
            iglob = ibool(iorderi(iloc),jmin,iorderj(iloc),ispec)
          else
            iglob = ibool(iorderi(iloc),iorderj(iloc),kmin,ispec)
          endif
          ! coordinates
          store_val_x_external_mesh(ipoin) = xstore(iglob)
          store_val_y_external_mesh(ipoin) = ystore(iglob)
          store_val_z_external_mesh(ipoin) = zstore(iglob)
        enddo
      endif
    enddo
  endif

  
! outputs values at free surface
  ipoin = 0
  do iface=1,num_free_surface_faces
    ispec = free_surface_ispec(iface)
    ! high_resolution
    if (USE_HIGHRES_FOR_MOVIES) then      
      do igll = 1, NGLLSQUARE
        ipoin = ipoin + 1
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)      
        iglob = ibool(i,j,k,ispec)
        ! coordinates
        !store_val_x_external_mesh(ipoin) = xstore(iglob)
        !store_val_y_external_mesh(ipoin) = ystore(iglob)
        !store_val_z_external_mesh(ipoin) = zstore(iglob)
        ! elastic displacement/velocity
        if( ispec_is_elastic(ispec) ) then
          if(SAVE_DISPLACEMENT) then
             store_val_ux_external_mesh(ipoin) = displ(1,iglob)
             store_val_uy_external_mesh(ipoin) = displ(2,iglob)
             store_val_uz_external_mesh(ipoin) = displ(3,iglob)
          else
             store_val_ux_external_mesh(ipoin) = veloc(1,iglob)
             store_val_uy_external_mesh(ipoin) = veloc(2,iglob)
             store_val_uz_external_mesh(ipoin) = veloc(3,iglob)
          endif
        endif
      enddo
    else    
      imin = minval( free_surface_ijk(1,:,iface) )
      imax = maxval( free_surface_ijk(1,:,iface) )
      jmin = minval( free_surface_ijk(2,:,iface) )
      jmax = maxval( free_surface_ijk(2,:,iface) )
      kmin = minval( free_surface_ijk(3,:,iface) )
      kmax = maxval( free_surface_ijk(3,:,iface) )      
      do iloc = 1, NGNOD2D    
        ipoin = ipoin + 1
        ! corner points
        if( imin == imax ) then
          iglob = ibool(imin,iorderi(iloc),iorderj(iloc),ispec)
        else if( jmin == jmax ) then
          iglob = ibool(iorderi(iloc),jmin,iorderj(iloc),ispec)
        else
          iglob = ibool(iorderi(iloc),iorderj(iloc),kmin,ispec)
        endif
        ! coordinates
        !store_val_x_external_mesh(ipoin) = xstore(iglob)
        !store_val_y_external_mesh(ipoin) = ystore(iglob)
        !store_val_z_external_mesh(ipoin) = zstore(iglob)
        ! elastic displacement/velocity
        if( ispec_is_elastic(ispec) ) then
          if(SAVE_DISPLACEMENT) then
             store_val_ux_external_mesh(ipoin) = displ(1,iglob)
             store_val_uy_external_mesh(ipoin) = displ(2,iglob)
             store_val_uz_external_mesh(ipoin) = displ(3,iglob)
          else
             store_val_ux_external_mesh(ipoin) = veloc(1,iglob)
             store_val_uy_external_mesh(ipoin) = veloc(2,iglob)
             store_val_uz_external_mesh(ipoin) = veloc(3,iglob)
          endif
        endif
      enddo ! iloc
    endif
  enddo ! iface

! master process collects all info
  if (USE_HIGHRES_FOR_MOVIES) then
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    endif
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
  else
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
  endif

! file output
  if(myrank == 0) then
    write(outputname,"('/moviedata_free_surface',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
    write(IOUT) store_val_x_all_external_mesh   ! x coordinate
    write(IOUT) store_val_y_all_external_mesh   ! y coordinate
    write(IOUT) store_val_z_all_external_mesh   ! z coordinate
    write(IOUT) store_val_ux_all_external_mesh  ! velocity x-component
    write(IOUT) store_val_uy_all_external_mesh  ! velocity y-component
    write(IOUT) store_val_uz_all_external_mesh  ! velocity z-component
    close(IOUT)
  endif

! obsolete...
!  ispec = nmovie_points
!
!  call gather_all_cr(store_val_x,ispec,store_val_x_all,ispec,NPROC)
!  call gather_all_cr(store_val_y,ispec,store_val_y_all,ispec,NPROC)
!  call gather_all_cr(store_val_z,ispec,store_val_z_all,ispec,NPROC)
!  call gather_all_cr(store_val_ux,ispec,store_val_ux_all,ispec,NPROC)
!  call gather_all_cr(store_val_uy,ispec,store_val_uy_all,ispec,NPROC)
!  call gather_all_cr(store_val_uz,ispec,store_val_uz_all,ispec,NPROC)
!
!! save movie data to disk in home directory
!  if(myrank == 0) then
!    write(outputname,"('/moviedata',i6.6)") it
!    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
!    write(IOUT) store_val_x_all
!    write(IOUT) store_val_y_all
!    write(IOUT) store_val_z_all
!    write(IOUT) store_val_ux_all
!    write(IOUT) store_val_uy_all
!    write(IOUT) store_val_uz_all
!    close(IOUT)
!  endif

  end subroutine it_movie_surface_output_o
  
  
!=====================================================================

  subroutine it_create_shakemap_o()

! outputs shakemap file 
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_movie
  
  implicit none
  integer :: imin,imax,jmin,jmax,kmin,kmax,iface,igll,iloc,ipoin
  integer :: ispec,i,j,k,iglob

! outputs values on free surface  
  ipoin = 0
  do iface=1,num_free_surface_faces
    ispec = free_surface_ispec(iface)
    ! save all points for high resolution, or only four corners for low resolution
    if(USE_HIGHRES_FOR_MOVIES) then
      do igll = 1, NGLLSQUARE
        ipoin = ipoin + 1
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        store_val_x_external_mesh(ipoin) = xstore(iglob)
        store_val_y_external_mesh(ipoin) = ystore(iglob)
        store_val_z_external_mesh(ipoin) = zstore(iglob)
        ! todo: are we only interested in the absolute maximum of horizontal (E,N) components?
        if( ispec_is_elastic( ispec) ) then
          ! horizontal displacement
          store_val_ux_external_mesh(ipoin) = max(store_val_ux_external_mesh(ipoin),abs(displ(1,iglob)),abs(displ(2,iglob)))
          ! horizontal velocity
          store_val_uy_external_mesh(ipoin) = max(store_val_uy_external_mesh(ipoin),abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          ! horizontal acceleration
          store_val_uz_external_mesh(ipoin) = max(store_val_uz_external_mesh(ipoin),abs(accel(1,iglob)),abs(accel(2,iglob)))
        endif
      enddo    
    else
      imin = minval( free_surface_ijk(1,:,iface) )
      imax = maxval( free_surface_ijk(1,:,iface) )
      jmin = minval( free_surface_ijk(2,:,iface) )
      jmax = maxval( free_surface_ijk(2,:,iface) )
      kmin = minval( free_surface_ijk(3,:,iface) )
      kmax = maxval( free_surface_ijk(3,:,iface) )
      do iloc = 1, NGNOD2D
        ipoin = ipoin + 1
        ! corner points
        if( imin == imax ) then
          iglob = ibool(imin,iorderi(iloc),iorderj(iloc),ispec)
        else if( jmin == jmax ) then
          iglob = ibool(iorderi(iloc),jmin,iorderj(iloc),ispec)
        else
          iglob = ibool(iorderi(iloc),iorderj(iloc),kmin,ispec)
        endif        
        ! coordinates
        store_val_x_external_mesh(ipoin) = xstore(iglob)
        store_val_y_external_mesh(ipoin) = ystore(iglob)
        store_val_z_external_mesh(ipoin) = zstore(iglob)
        ! todo: are we only interested in the absolute maximum of horizontal (E,N) components?
        if( ispec_is_elastic( ispec) ) then        
          store_val_ux_external_mesh(ipoin) = max(store_val_ux_external_mesh(ipoin),abs(displ(1,iglob)),abs(displ(2,iglob)))
          store_val_uy_external_mesh(ipoin) = max(store_val_uy_external_mesh(ipoin),abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          store_val_uz_external_mesh(ipoin) = max(store_val_uz_external_mesh(ipoin),abs(accel(1,iglob)),abs(accel(2,iglob)))
        endif
      enddo
    endif ! USE_HIGHRES_FOR_MOVIES
  enddo

! save shakemap only at the end of the simulation
  if(it == NSTEP) then
    if (USE_HIGHRES_FOR_MOVIES) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    else
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif

! creates shakemap file
    if(myrank == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata_freesurface',status='unknown',form='unformatted')
      write(IOUT) store_val_x_all_external_mesh   ! x coordinates
      write(IOUT) store_val_y_all_external_mesh   ! y coordinates
      write(IOUT) store_val_z_all_external_mesh   ! z coordinates
      write(IOUT) store_val_ux_all_external_mesh  ! norm of displacement vector
      write(IOUT) store_val_uy_all_external_mesh  ! norm of velocity vector
      write(IOUT) store_val_uz_all_external_mesh  ! norm of acceleration vector
      close(IOUT)
    endif

! obsolete...  
!    ispec = nmovie_points
!    call gather_all_cr(store_val_x,ispec,store_val_x_all,ispec,NPROC)
!    call gather_all_cr(store_val_y,ispec,store_val_y_all,ispec,NPROC)
!    call gather_all_cr(store_val_z,ispec,store_val_z_all,ispec,NPROC)
!    call gather_all_cr(store_val_norm_displ,ispec,store_val_ux_all,ispec,NPROC)
!    call gather_all_cr(store_val_norm_veloc,ispec,store_val_uy_all,ispec,NPROC)
!    call gather_all_cr(store_val_norm_accel,ispec,store_val_uz_all,ispec,NPROC)
!
!! save movie data to disk in home directory
!    if(myrank == 0) then
!      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata',status='unknown',form='unformatted')
!      write(IOUT) store_val_x_all
!      write(IOUT) store_val_y_all
!      write(IOUT) store_val_z_all
!! this saves norm of displacement, velocity and acceleration
!! but we use the same ux, uy, uz arrays as for the movies to save memory
!      write(IOUT) store_val_ux_all
!      write(IOUT) store_val_uy_all
!      write(IOUT) store_val_uz_all
!      close(IOUT)
!    endif
!
  endif ! NTSTEP

  end subroutine it_create_shakemap_o

    
!=====================================================================

  subroutine it_movie_volume_output()

! outputs movie files for div, curl and velocity  
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie
  implicit none
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: veloc_element
  integer :: ispec,i,j,k,l,iglob
  
! save velocity here to avoid static offset on displacement for movies
  velocity_movie(:,:,:,:,:) = 0._CUSTOM_REAL
  
  if( ACOUSTIC_SIMULATION ) then
    ! uses div as temporary array to store velocity on all gll points
    do ispec=1,NSPEC_AB
      if( .not. ispec_is_acoustic(ispec) ) cycle

      ! calculates velocity
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                        potential_dot_acoustic, veloc_element,&
                        hprime_xx,hprime_yy,hprime_zz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        ibool,rhostore)
      velocity_movie(:,:,:,:,ispec) = veloc_element(:,:,:,:)
    enddo
  endif ! acoustic

! save full snapshot data to local disk
  if( ELASTIC_SIMULATION ) then

  ! calculate strain div and curl
    do ispec=1,NSPEC_AB
      if( .not. ispec_is_elastic(ispec) ) cycle
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempx3l = 0._CUSTOM_REAL

            tempy1l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL

            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL

            do l=1,NGLLX
              hp1 = hprime_xx(i,l)
              iglob = ibool(l,j,k,ispec)
              tempx1l = tempx1l + veloc(1,iglob)*hp1
              tempy1l = tempy1l + veloc(2,iglob)*hp1
              tempz1l = tempz1l + veloc(3,iglob)*hp1
  !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

  !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
              hp2 = hprime_yy(j,l)
              iglob = ibool(i,l,k,ispec)
              tempx2l = tempx2l + veloc(1,iglob)*hp2
              tempy2l = tempy2l + veloc(2,iglob)*hp2
              tempz2l = tempz2l + veloc(3,iglob)*hp2
  !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

  !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
              hp3 = hprime_zz(k,l)
              iglob = ibool(i,j,l,ispec)
              tempx3l = tempx3l + veloc(1,iglob)*hp3
              tempy3l = tempy3l + veloc(2,iglob)*hp3
              tempz3l = tempz3l + veloc(3,iglob)*hp3
            enddo

  !         get derivatives of ux, uy and uz with respect to x, y and z

            xixl = xix(i,j,k,ispec)
            xiyl = xiy(i,j,k,ispec)
            xizl = xiz(i,j,k,ispec)
            etaxl = etax(i,j,k,ispec)
            etayl = etay(i,j,k,ispec)
            etazl = etaz(i,j,k,ispec)
            gammaxl = gammax(i,j,k,ispec)
            gammayl = gammay(i,j,k,ispec)
            gammazl = gammaz(i,j,k,ispec)

            dvxdxl(i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
            dvxdyl(i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
            dvxdzl(i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

            dvydxl(i,j,k) = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
            dvydyl(i,j,k) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
            dvydzl(i,j,k) = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

            dvzdxl(i,j,k) = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
            dvzdyl(i,j,k) = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
            dvzdzl(i,j,k) = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

          enddo
        enddo
      enddo

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            div(i,j,k,ispec) = dvxdxl(i,j,k) + dvydyl(i,j,k) + dvzdzl(i,j,k)
            curl_x(i,j,k,ispec) = dvzdyl(i,j,k) - dvydzl(i,j,k)
            curl_y(i,j,k,ispec) = dvxdzl(i,j,k) - dvzdxl(i,j,k)
            curl_z(i,j,k,ispec) = dvydxl(i,j,k) - dvxdyl(i,j,k)
            
            iglob = ibool(i,j,k,ispec)
            velocity_movie(:,i,j,k,ispec) = veloc(:,iglob)
          enddo
        enddo
      enddo
    enddo !NSPEC_AB

    write(outputname,"('/proc',i6.6,'_div_it',i6.6,'.bin')") myrank,it
    open(unit=27,file='OUTPUT_FILES'//trim(outputname),status='unknown',form='unformatted')
    write(27) div
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_x_it',i6.6,'.bin')") myrank,it
    open(unit=27,file='OUTPUT_FILES'//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_x
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_y_it',i6.6,'.bin')") myrank,it
    open(unit=27,file='OUTPUT_FILES'//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_y
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_z_it',i6.6,'.bin')") myrank,it
    open(unit=27,file='OUTPUT_FILES'//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_z
    close(27)
    
    !write(outputname,"('veloc_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    !open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    !write(27) veloc
    !close(27)
  
  endif ! elastic
 
  if( ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION ) then
    write(outputname,"('/proc',i6.6,'_veloc_it',i6.6,'.bin')") myrank,it
    open(unit=27,file='OUTPUT_FILES'//trim(outputname),status='unknown',form='unformatted')
    write(27) velocity_movie
    close(27)  
  endif 
  
  end subroutine it_movie_volume_output
  