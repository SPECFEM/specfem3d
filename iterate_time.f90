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


!check stability
  do i=1,3
    Usolidnorm = maxval(abs(displ(i,:)))
    Usolidnorm_index = maxloc(abs(displ(i,:)))
    if(Usolidnorm > 1.e+15 ) then        
      print*,' stability issue:',myrank
      print*,'  norm: ',Usolidnorm,displ(i,Usolidnorm_index(1)),i
      print*,'  index: ',Usolidnorm_index(1)
      print*,'  x/y/z: ',xstore(Usolidnorm_index(1)),ystore(Usolidnorm_index(1)),zstore(Usolidnorm_index(1))
      print*,'  time step: ',it
      call exit_MPI(myrank,'forward simulation became unstable and blew up')
    endif
  enddo
! compute the maximum of the norm of the displacement
! in all the slices using an MPI reduction
! and output timestamp file to check that simulation is running fine
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5) then

! compute maximum of norm of displacement in each slice
    Usolidnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2 + displ(3,:)**2))

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
      write(IMAIN,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all
!     if (SIMULATION_TYPE == 3) write(IMAIN,*) &
!           'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
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

    endif
  endif





! update displacement using finite difference time scheme
  displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
  accel(:,:) = 0._CUSTOM_REAL

!! DK DK array not created yet for CUBIT
! if (SIMULATION_TYPE == 3) then
!   b_displ(:,:) = b_displ(:,:) + b_deltat*b_veloc(:,:) + b_deltatsqover2*b_accel(:,:)
!   b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)
!   b_accel(:,:) = 0._CUSTOM_REAL
! endif

! if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
!   ispec2D_moho_top = 0
!   ispec2D_moho_bot = 0
! endif

! assemble all the contributions between slices using MPI


    if(USE_DEVILLE_PRODUCTS) then
      call compute_forces_with_Deville(NSPEC_AB,NGLOB_AB,ATTENUATION,displ,accel,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
         hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
         kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh,.false., &
         NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,xi_source,eta_source,gamma_source,nu_source, &
         hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, & 
         one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,R_xx,R_yy,R_xy,R_xz,R_yz, &
         epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,iflag_attenuation_store,ABSORBING_CONDITIONS, &
         nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2DMAX_XMIN_XMAX_ext,NSPEC2DMAX_YMIN_YMAX_ext, &
         ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom, &
         nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
         veloc,rho_vp,rho_vs,jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom, &
         normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom) 
    else
      call compute_forces_no_Deville(NSPEC_AB,NGLOB_AB,displ,accel,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
         hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
         kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh,.false., &
         NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,xi_source,eta_source,gamma_source,nu_source,hdur,dt)
    endif

    call assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_AB,accel, &
         buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
         ninterfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
         request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    if(USE_DEVILLE_PRODUCTS) then
      call compute_forces_with_Deville(NSPEC_AB,NGLOB_AB,ATTENUATION,displ,accel,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
         hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
         kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh,.true., &
         NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,xi_source,eta_source,gamma_source,nu_source, &
         hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, & 
         one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,R_xx,R_yy,R_xy,R_xz,R_yz, &
         epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,iflag_attenuation_store,ABSORBING_CONDITIONS, &
         nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2DMAX_XMIN_XMAX_ext,NSPEC2DMAX_YMIN_YMAX_ext, &
         ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom, &
         nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
         veloc,rho_vp,rho_vs,jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom, &
         normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom)
    else
      call compute_forces_no_Deville(NSPEC_AB,NGLOB_AB,displ,accel,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
         hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
         kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh,.true., &
         NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,xi_source,eta_source,gamma_source,nu_source,hdur,dt)
    endif

    call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_AB,accel, &
         buffer_recv_vector_ext_mesh,ninterfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
         request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

!! DK DK May 2009: removed this because now each slice of a CUBIT + SCOTCH mesh
!! DK DK May 2009: has a different number of spectral elements and therefore
!! DK DK May 2009: only the general non-blocking MPI routines assemble_MPI_vector_ext_mesh_s
!! DK DK May 2009: and assemble_MPI_vector_ext_mesh_w above can be used.
!! DK DK May 2009: For adjoint runs below (SIMULATION_TYPE == 3) they should be used as well.
! if (SIMULATION_TYPE == 3) call assemble_MPI_vector(b_accel,iproc_xi,iproc_eta,addressing, &
!         iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
!         buffer_send_faces_vector,buffer_received_faces_vector,npoin2D_xi,npoin2D_eta, &
!         NPROC_XI,NPROC_ETA,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY)

! multiply by the inverse of the mass matrix
  accel(1,:) = accel(1,:)*rmass(:)
  accel(2,:) = accel(2,:)*rmass(:)
  accel(3,:) = accel(3,:)*rmass(:)

!! DK DK array not created yet for CUBIT
! if (SIMULATION_TYPE == 3) then
!   b_accel(1,:) = b_accel(1,:)*rmass(:)
!   b_accel(2,:) = b_accel(2,:)*rmass(:)
!   b_accel(3,:) = b_accel(3,:)*rmass(:)
! endif

  if(OCEANS) then

    stop 'DK DK oceans have been removed for now because we need a flag to detect the surface elements'

!   initialize the updates
    updated_dof_ocean_load(:) = .false.

! for surface elements exactly at the top of the model (ocean bottom)
    do ispec2D = 1,NSPEC2D_TOP

!! DK DK array not created yet for CUBIT      ispec = ibelm_top(ispec2D)

! only for DOFs exactly at the top of the model (ocean bottom)
      k = NGLLZ

      do j = 1,NGLLY
        do i = 1,NGLLX

! get global point number
          iglob = ibool(i,j,k,ispec)

! only update once
          if(.not. updated_dof_ocean_load(iglob)) then

! get normal
!! DK DK array not created yet for CUBIT            nx = normal_top(1,i,j,ispec2D)
!! DK DK array not created yet for CUBIT            ny = normal_top(2,i,j,ispec2D)
!! DK DK array not created yet for CUBIT            nz = normal_top(3,i,j,ispec2D)

! make updated component of right-hand side
! we divide by rmass() which is 1 / M
! we use the total force which includes the Coriolis term above
            force_normal_comp = (accel(1,iglob)*nx + &
                 accel(2,iglob)*ny + accel(3,iglob)*nz) / rmass(iglob)

            additional_term = (rmass_ocean_load(iglob) - rmass(iglob)) * force_normal_comp

            accel(1,iglob) = accel(1,iglob) + additional_term * nx
            accel(2,iglob) = accel(2,iglob) + additional_term * ny
            accel(3,iglob) = accel(3,iglob) + additional_term * nz

            if (SIMULATION_TYPE == 3) then
!! DK DK array not created yet for CUBIT
!             b_force_normal_comp = (b_accel(1,iglob)*nx + &
!                   b_accel(2,iglob)*ny + b_accel(3,iglob)*nz) / rmass(iglob)

              b_additional_term = (rmass_ocean_load(iglob) - rmass(iglob)) * b_force_normal_comp

!! DK DK array not created yet for CUBIT
!             b_accel(1,iglob) = b_accel(1,iglob) + b_additional_term * nx
!             b_accel(2,iglob) = b_accel(2,iglob) + b_additional_term * ny
!             b_accel(3,iglob) = b_accel(3,iglob) + b_additional_term * nz
            endif

!           done with this point
            updated_dof_ocean_load(iglob) = .true.

          endif

        enddo
      enddo
    enddo
  endif

  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)

!! DK DK array not created yet for CUBIT
! if (SIMULATION_TYPE == 3) b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)

! write the seismograms with time shift
  if (nrec_local > 0) then
  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

! perform the general interpolation using Lagrange polynomials
    if(FASTER_RECEIVERS_POINTS_ONLY) then

      iglob = ibool(nint(xi_receiver(irec)),nint(eta_receiver(irec)), &
           nint(gamma_receiver(irec)),ispec_selected_rec(irec))
      dxd = dble(displ(1,iglob))
      dyd = dble(displ(2,iglob))
      dzd = dble(displ(3,iglob))
      vxd = dble(veloc(1,iglob))
      vyd = dble(veloc(2,iglob))
      vzd = dble(veloc(3,iglob))
      axd = dble(accel(1,iglob))
      ayd = dble(accel(2,iglob))
      azd = dble(accel(3,iglob))

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

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

! receivers are always located at the surface of the mesh
            iglob = ibool(i,j,k,ispec_selected_rec(irec))

            hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)


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

          enddo
        enddo
      enddo

    else if (SIMULATION_TYPE == 2) then

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            iglob = ibool(i,j,k,ispec_selected_source(irec))

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

      ispec = ispec_selected_source(irec)

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

      do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool(i,j,k,ispec_selected_rec(irec))

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
    endif

    endif ! end of if(FASTER_RECEIVERS_POINTS_ONLY)

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

  enddo

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

  endif ! nrec_local

! resetting d/v/a/R/eps for the backward reconstruction with attenuation
  if (ATTENUATION .and. it > 1 .and. it < NSTEP) then
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
  endif
  endif

  if (EXTERNAL_MESH_CREATE_SHAKEMAP) then
    if (it == 1) then

      store_val_ux_external_mesh(:) = -HUGEVAL
      store_val_uy_external_mesh(:) = -HUGEVAL
      store_val_uz_external_mesh(:) = -HUGEVAL
      do ispec = 1,nfaces_surface_external_mesh
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, NGLLX*NGLLY
          store_val_x_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = xstore(faces_surface_external_mesh(ipoin,ispec))
          store_val_y_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = ystore(faces_surface_external_mesh(ipoin,ispec))
          store_val_z_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = zstore(faces_surface_external_mesh(ipoin,ispec))
        enddo
      else
        store_val_x_external_mesh(NGNOD2D*(ispec-1)+1) = xstore(faces_surface_external_mesh(1,ispec))
        store_val_x_external_mesh(NGNOD2D*(ispec-1)+2) = xstore(faces_surface_external_mesh(2,ispec))
        store_val_x_external_mesh(NGNOD2D*(ispec-1)+3) = xstore(faces_surface_external_mesh(3,ispec))
        store_val_x_external_mesh(NGNOD2D*(ispec-1)+4) = xstore(faces_surface_external_mesh(4,ispec))
        store_val_y_external_mesh(NGNOD2D*(ispec-1)+1) = ystore(faces_surface_external_mesh(1,ispec))
        store_val_y_external_mesh(NGNOD2D*(ispec-1)+2) = ystore(faces_surface_external_mesh(2,ispec))
        store_val_y_external_mesh(NGNOD2D*(ispec-1)+3) = ystore(faces_surface_external_mesh(3,ispec))
        store_val_y_external_mesh(NGNOD2D*(ispec-1)+4) = ystore(faces_surface_external_mesh(4,ispec))
        store_val_z_external_mesh(NGNOD2D*(ispec-1)+1) = zstore(faces_surface_external_mesh(1,ispec))
        store_val_z_external_mesh(NGNOD2D*(ispec-1)+2) = zstore(faces_surface_external_mesh(2,ispec))
        store_val_z_external_mesh(NGNOD2D*(ispec-1)+3) = zstore(faces_surface_external_mesh(3,ispec))
        store_val_z_external_mesh(NGNOD2D*(ispec-1)+4) = zstore(faces_surface_external_mesh(4,ispec))
      endif
      enddo
    endif

    do ispec = 1,nfaces_surface_external_mesh
    if (USE_HIGHRES_FOR_MOVIES) then
      do ipoin = 1, NGLLX*NGLLY
        store_val_ux_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = &
             max(store_val_ux_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin), &
             sqrt(displ(1,faces_surface_external_mesh(ipoin,ispec))**2 + &
             displ(2,faces_surface_external_mesh(ipoin,ispec))**2 + &
             displ(3,faces_surface_external_mesh(ipoin,ispec))**2))
        store_val_uy_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = &
             max(store_val_uy_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin), &
             sqrt(veloc(1,faces_surface_external_mesh(ipoin,ispec))**2 + &
             veloc(2,faces_surface_external_mesh(ipoin,ispec))**2 + &
             veloc(3,faces_surface_external_mesh(ipoin,ispec))**2))
        store_val_uz_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = &
             max(store_val_uz_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin), &
             sqrt(accel(1,faces_surface_external_mesh(ipoin,ispec))**2 + &
             accel(2,faces_surface_external_mesh(ipoin,ispec))**2 + &
             accel(3,faces_surface_external_mesh(ipoin,ispec))**2))

      enddo
    else
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+1) = &
           max(store_val_ux_external_mesh(NGNOD2D*(ispec-1)+1), &
           sqrt(displ(1,faces_surface_external_mesh(1,ispec))**2 + &
           displ(2,faces_surface_external_mesh(1,ispec))**2 + &
           displ(3,faces_surface_external_mesh(1,ispec))**2))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+2) = &
           max(store_val_ux_external_mesh(NGNOD2D*(ispec-1)+2), &
           sqrt(displ(1,faces_surface_external_mesh(2,ispec))**2 + &
           displ(2,faces_surface_external_mesh(2,ispec))**2 + &
           displ(3,faces_surface_external_mesh(2,ispec))**2))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+3) = &
           max(store_val_ux_external_mesh(NGNOD2D*(ispec-1)+3), &
           sqrt(displ(1,faces_surface_external_mesh(3,ispec))**2 + &
           displ(2,faces_surface_external_mesh(3,ispec))**2 + &
           displ(3,faces_surface_external_mesh(3,ispec))**2))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+4) = &
           max(store_val_ux_external_mesh(NGNOD2D*(ispec-1)+4), &
           sqrt(displ(1,faces_surface_external_mesh(4,ispec))**2 + &
           displ(2,faces_surface_external_mesh(4,ispec))**2 + &
           displ(3,faces_surface_external_mesh(4,ispec))**2))
     store_val_uy_external_mesh(NGNOD2D*(ispec-1)+1) = &
           max(store_val_uy_external_mesh(NGNOD2D*(ispec-1)+1), &
           sqrt(veloc(1,faces_surface_external_mesh(1,ispec))**2 + &
           veloc(2,faces_surface_external_mesh(1,ispec))**2 + &
           veloc(3,faces_surface_external_mesh(1,ispec))**2))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+2) = &
           max(store_val_uy_external_mesh(NGNOD2D*(ispec-1)+2), &
           sqrt(veloc(1,faces_surface_external_mesh(2,ispec))**2 + &
           veloc(2,faces_surface_external_mesh(2,ispec))**2 + &
           veloc(3,faces_surface_external_mesh(2,ispec))**2))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+3) = &
           max(store_val_uy_external_mesh(NGNOD2D*(ispec-1)+3), &
           sqrt(veloc(1,faces_surface_external_mesh(3,ispec))**2 + &
           veloc(2,faces_surface_external_mesh(3,ispec))**2 + &
           veloc(3,faces_surface_external_mesh(3,ispec))**2))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+4) = &
           max(store_val_uy_external_mesh(NGNOD2D*(ispec-1)+4), &
           sqrt(veloc(1,faces_surface_external_mesh(4,ispec))**2 + &
           veloc(2,faces_surface_external_mesh(4,ispec))**2 + &
           veloc(3,faces_surface_external_mesh(4,ispec))**2))
     store_val_uz_external_mesh(NGNOD2D*(ispec-1)+1) = &
           max(store_val_uz_external_mesh(NGNOD2D*(ispec-1)+1), &
           sqrt(accel(1,faces_surface_external_mesh(1,ispec))**2 + &
           accel(2,faces_surface_external_mesh(1,ispec))**2 + &
           accel(3,faces_surface_external_mesh(1,ispec))**2))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+2) = &
           max(store_val_uz_external_mesh(NGNOD2D*(ispec-1)+2), &
           sqrt(accel(1,faces_surface_external_mesh(2,ispec))**2 + &
           accel(2,faces_surface_external_mesh(2,ispec))**2 + &
           accel(3,faces_surface_external_mesh(2,ispec))**2))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+3) = &
           max(store_val_uz_external_mesh(NGNOD2D*(ispec-1)+3), &
           sqrt(accel(1,faces_surface_external_mesh(3,ispec))**2 + &
           accel(2,faces_surface_external_mesh(3,ispec))**2 + &
           accel(3,faces_surface_external_mesh(3,ispec))**2))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+4) = &
           max(store_val_uz_external_mesh(NGNOD2D*(ispec-1)+4), &
           sqrt(accel(1,faces_surface_external_mesh(4,ispec))**2 + &
           accel(2,faces_surface_external_mesh(4,ispec))**2 + &
           accel(3,faces_surface_external_mesh(4,ispec))**2))
    endif
    enddo

    if (it == NSTEP) then
    if (USE_HIGHRES_FOR_MOVIES) then
    call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    else
    call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif

    if(myrank == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata',status='unknown',form='unformatted')
      write(IOUT) store_val_x_all_external_mesh
      write(IOUT) store_val_y_all_external_mesh
      write(IOUT) store_val_z_all_external_mesh
      write(IOUT) store_val_ux_all_external_mesh
      write(IOUT) store_val_uy_all_external_mesh
      write(IOUT) store_val_uz_all_external_mesh
      close(IOUT)
    endif
    endif

 endif

  if(EXTERNAL_MESH_MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
! get coordinates of surface mesh and surface displacement
    do ispec = 1,nfaces_surface_external_mesh
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, NGLLX*NGLLY
          store_val_x_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = xstore(faces_surface_external_mesh(ipoin,ispec))
          store_val_y_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = ystore(faces_surface_external_mesh(ipoin,ispec))
          store_val_z_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = zstore(faces_surface_external_mesh(ipoin,ispec))
          store_val_ux_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = veloc(1,faces_surface_external_mesh(ipoin,ispec))
          store_val_uy_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = veloc(2,faces_surface_external_mesh(ipoin,ispec))
          store_val_uz_external_mesh(NGLLX*NGLLY*(ispec-1)+ipoin) = veloc(3,faces_surface_external_mesh(ipoin,ispec))
        enddo
      else
      store_val_x_external_mesh(NGNOD2D*(ispec-1)+1) = xstore(faces_surface_external_mesh(1,ispec))
      store_val_x_external_mesh(NGNOD2D*(ispec-1)+2) = xstore(faces_surface_external_mesh(2,ispec))
      store_val_x_external_mesh(NGNOD2D*(ispec-1)+3) = xstore(faces_surface_external_mesh(3,ispec))
      store_val_x_external_mesh(NGNOD2D*(ispec-1)+4) = xstore(faces_surface_external_mesh(4,ispec))
      store_val_y_external_mesh(NGNOD2D*(ispec-1)+1) = ystore(faces_surface_external_mesh(1,ispec))
      store_val_y_external_mesh(NGNOD2D*(ispec-1)+2) = ystore(faces_surface_external_mesh(2,ispec))
      store_val_y_external_mesh(NGNOD2D*(ispec-1)+3) = ystore(faces_surface_external_mesh(3,ispec))
      store_val_y_external_mesh(NGNOD2D*(ispec-1)+4) = ystore(faces_surface_external_mesh(4,ispec))
      store_val_z_external_mesh(NGNOD2D*(ispec-1)+1) = zstore(faces_surface_external_mesh(1,ispec))
      store_val_z_external_mesh(NGNOD2D*(ispec-1)+2) = zstore(faces_surface_external_mesh(2,ispec))
      store_val_z_external_mesh(NGNOD2D*(ispec-1)+3) = zstore(faces_surface_external_mesh(3,ispec))
      store_val_z_external_mesh(NGNOD2D*(ispec-1)+4) = zstore(faces_surface_external_mesh(4,ispec))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+1) = veloc(1,faces_surface_external_mesh(1,ispec))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+2) = veloc(1,faces_surface_external_mesh(2,ispec))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+3) = veloc(1,faces_surface_external_mesh(3,ispec))
      store_val_ux_external_mesh(NGNOD2D*(ispec-1)+4) = veloc(1,faces_surface_external_mesh(4,ispec))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+1) = veloc(2,faces_surface_external_mesh(1,ispec))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+2) = veloc(2,faces_surface_external_mesh(2,ispec))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+3) = veloc(2,faces_surface_external_mesh(3,ispec))
      store_val_uy_external_mesh(NGNOD2D*(ispec-1)+4) = veloc(2,faces_surface_external_mesh(4,ispec))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+1) = veloc(3,faces_surface_external_mesh(1,ispec))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+2) = veloc(3,faces_surface_external_mesh(2,ispec))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+3) = veloc(3,faces_surface_external_mesh(3,ispec))
      store_val_uz_external_mesh(NGNOD2D*(ispec-1)+4) = veloc(3,faces_surface_external_mesh(4,ispec))
      endif
    enddo

    if (USE_HIGHRES_FOR_MOVIES) then
    call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_external_mesh*NGLLX*NGLLY,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    else
    call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_external_mesh*NGNOD2D,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif

    if(myrank == 0) then
      write(outputname,"('/moviedata',i6.6)") it
      open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
      write(IOUT) store_val_x_all_external_mesh
      write(IOUT) store_val_y_all_external_mesh
      write(IOUT) store_val_z_all_external_mesh
      write(IOUT) store_val_ux_all_external_mesh
      write(IOUT) store_val_uy_all_external_mesh
      write(IOUT) store_val_uz_all_external_mesh
      close(IOUT)
    endif
  endif

! save MOVIE on the SURFACE
  if(MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

    stop 'DK DK MOVIE_SURFACE has been removed for now because we need a flag to detect the surface elements'

! get coordinates of surface mesh and surface displacement
    ipoin = 0

   k = NGLLZ
   if (USE_HIGHRES_FOR_MOVIES) then
     do ispec2D = 1,NSPEC2D_TOP
!! DK DK array not created yet for CUBIT       ispec = ibelm_top(ispec2D)
       do j = 1,NGLLY
         do i = 1,NGLLX
           ipoin = ipoin + 1
           iglob = ibool(i,j,k,ispec)
           store_val_x(ipoin) = xstore(iglob)
           store_val_y(ipoin) = ystore(iglob)
           store_val_z(ipoin) = zstore(iglob)
           if(SAVE_DISPLACEMENT) then
             store_val_ux(ipoin) = displ(1,iglob)
             store_val_uy(ipoin) = displ(2,iglob)
             store_val_uz(ipoin) = displ(3,iglob)
           else
             store_val_ux(ipoin) = veloc(1,iglob)
             store_val_uy(ipoin) = veloc(2,iglob)
             store_val_uz(ipoin) = veloc(3,iglob)
           endif
         enddo
       enddo
     enddo ! ispec_top
   else
     do ispec2D = 1,NSPEC2D_TOP
!! DK DK array not created yet for CUBIT       ispec = ibelm_top(ispec2D)
       do iloc = 1, NGNOD2D
         ipoin = ipoin + 1
         iglob = ibool(iorderi(iloc),iorderj(iloc),k,ispec)
         store_val_x(ipoin) = xstore(iglob)
         store_val_y(ipoin) = ystore(iglob)
         store_val_z(ipoin) = zstore(iglob)
         if(SAVE_DISPLACEMENT) then
           store_val_ux(ipoin) = displ(1,iglob)
           store_val_uy(ipoin) = displ(2,iglob)
           store_val_uz(ipoin) = displ(3,iglob)
         else
           store_val_ux(ipoin) = veloc(1,iglob)
           store_val_uy(ipoin) = veloc(2,iglob)
           store_val_uz(ipoin) = veloc(3,iglob)
         endif
       enddo
     enddo ! ispec_top
   endif

    ispec = nmovie_points

    call gather_all_cr(store_val_x,ispec,store_val_x_all,ispec,NPROC)
    call gather_all_cr(store_val_y,ispec,store_val_y_all,ispec,NPROC)
    call gather_all_cr(store_val_z,ispec,store_val_z_all,ispec,NPROC)
    call gather_all_cr(store_val_ux,ispec,store_val_ux_all,ispec,NPROC)
    call gather_all_cr(store_val_uy,ispec,store_val_uy_all,ispec,NPROC)
    call gather_all_cr(store_val_uz,ispec,store_val_uz_all,ispec,NPROC)

! save movie data to disk in home directory
    if(myrank == 0) then
      write(outputname,"('/moviedata',i6.6)") it
      open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
      write(IOUT) store_val_x_all
      write(IOUT) store_val_y_all
      write(IOUT) store_val_z_all
      write(IOUT) store_val_ux_all
      write(IOUT) store_val_uy_all
      write(IOUT) store_val_uz_all
      close(IOUT)
    endif

  endif

! compute SHAKING INTENSITY MAP
 if(CREATE_SHAKEMAP) then

    stop 'DK DK CREATE_SHAKEMAP has been removed for now because we need a flag to detect the surface elements'

    ipoin = 0
    k = NGLLZ

! save all points for high resolution, or only four corners for low resolution
    if(USE_HIGHRES_FOR_MOVIES) then

    do ispec2D = 1,NSPEC2D_TOP
!! DK DK array not created yet for CUBIT      ispec = ibelm_top(ispec2D)

! loop on all the points inside the element
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1
          iglob = ibool(i,j,k,ispec)
          store_val_x(ipoin) = xstore(iglob)
          store_val_y(ipoin) = ystore(iglob)
          store_val_z(ipoin) = zstore(iglob)
          store_val_norm_displ(ipoin) = max(store_val_norm_displ(ipoin),abs(displ(1,iglob)),abs(displ(2,iglob)))
          store_val_norm_veloc(ipoin) = max(store_val_norm_veloc(ipoin),abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          store_val_norm_accel(ipoin) = max(store_val_norm_accel(ipoin),abs(accel(1,iglob)),abs(accel(2,iglob)))
        enddo
      enddo
    enddo

    else
      do ispec2D = 1,NSPEC2D_TOP
!! DK DK array not created yet for CUBIT        ispec = ibelm_top(ispec2D)
        do iloc = 1, NGNOD2D
          ipoin = ipoin + 1
          iglob = ibool(iorderi(iloc),iorderj(iloc),k,ispec)
          store_val_x(ipoin) = xstore(iglob)
          store_val_y(ipoin) = ystore(iglob)
          store_val_z(ipoin) = zstore(iglob)
          store_val_norm_displ(ipoin) = max(store_val_norm_displ(ipoin),abs(displ(1,iglob)),abs(displ(2,iglob)))
          store_val_norm_veloc(ipoin) = max(store_val_norm_veloc(ipoin),abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          store_val_norm_accel(ipoin) = max(store_val_norm_accel(ipoin),abs(accel(1,iglob)),abs(accel(2,iglob)))
        enddo
      enddo
    endif

! save shakemap only at the end of the simulation
    if(it == NSTEP) then
    ispec = nmovie_points
    call gather_all_cr(store_val_x,ispec,store_val_x_all,ispec,NPROC)
    call gather_all_cr(store_val_y,ispec,store_val_y_all,ispec,NPROC)
    call gather_all_cr(store_val_z,ispec,store_val_z_all,ispec,NPROC)
    call gather_all_cr(store_val_norm_displ,ispec,store_val_ux_all,ispec,NPROC)
    call gather_all_cr(store_val_norm_veloc,ispec,store_val_uy_all,ispec,NPROC)
    call gather_all_cr(store_val_norm_accel,ispec,store_val_uz_all,ispec,NPROC)

! save movie data to disk in home directory
    if(myrank == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata',status='unknown',form='unformatted')
      write(IOUT) store_val_x_all
      write(IOUT) store_val_y_all
      write(IOUT) store_val_z_all
! this saves norm of displacement, velocity and acceleration
! but we use the same ux, uy, uz arrays as for the movies to save memory
      write(IOUT) store_val_ux_all
      write(IOUT) store_val_uy_all
      write(IOUT) store_val_uz_all
      close(IOUT)
    endif

    endif
  endif

! save MOVIE in full 3D MESH
  if(MOVIE_VOLUME .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

! save velocity here to avoid static offset on displacement for movies

! save full snapshot data to local disk

! calculate strain div and curl
    do ispec=1,NSPEC_AB

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
          enddo
        enddo
      enddo
    enddo

    write(outputname,"('div_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) div
    close(27)
    write(outputname,"('curl_x_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_x
    close(27)
    write(outputname,"('curl_y_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_y
    close(27)
    write(outputname,"('curl_z_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_z
    close(27)
    write(outputname,"('veloc_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) veloc
    close(27)

  endif

!
!---- end of time iteration loop
!
  enddo   ! end of main time loop



  end subroutine