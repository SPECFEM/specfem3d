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

  subroutine write_seismograms()

! writes the seismograms with time shift

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element,veloc_element,accel_element
  ! interpolated wavefield values
  double precision :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd
  ! receiver position
  double precision :: xi_r,eta_r,gamma_r

  integer :: irec_local,irec
  integer :: iglob,ispec,i,j,k

  ! adjoint locals
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM):: eps_s
  real(kind=CUSTOM_REAL),dimension(NDIM):: eps_m_s
  real(kind=CUSTOM_REAL):: stf_deltat
  double precision :: stf

  ! TODO: Test and Fix CUDA seismograms code.
  logical, parameter :: USE_CUDA_SEISMOGRAMS = .false.

  ! checks if anything to do
  if (.not. (nrec_local > 0 .or. (WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0))) return

  ! note: there might be some confusion about adjoint simulations, i.e. adjoint receivers and adjoint sources:
  !       for pure adjoint simulations (SIMULATION_TYPE == 2), CMT source locations become adjoint receivers
  !       for recording seismograms; station locations become (possible) adjoint sources.
  !       thus,
  !       1. adjoint sources are located at the receiver positions given in STATIONS file,
  !          'nadj_rec_local' is the number of local adjoint sources, i.e. station positions acting as adjoint source
  !       2. adjoint "receivers" are located at the CMT source positions given in CMTSOLUTION file,
  !          'nrec_local' is the number of local adjoint receivers, i.e. source positions acting as receivers for
  !          recording adjoint seismograms.

  ! remember for pure adjoint runs (see setup_receivers routine):
  ! - nrec_local              -> between 1 to NSOURCES
  ! - number_receiver_global  -> local to global mapping for "adjoint" receivers, i.e. at source positions
  ! - ispec_selected_source   -> element containing source position, which becomes an adjoint "receiver"

  ! gets resulting array values onto CPU
  if (GPU_MODE) then
    if (nrec_local > 0) then
      ! this transfers fields only in elements with stations for efficiency
      if (ACOUSTIC_SIMULATION) then
        ! only copy corresponding elements to CPU host
        ! timing: Elapsed time: 5.230904e-04
        call transfer_station_ac_from_device(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        Mesh_pointer,number_receiver_global, &
                        ispec_selected_rec,ispec_selected_source,ibool)

        ! alternative: transfers whole fields
        ! timing: Elapsed time: 4.138947e-03
        !call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
        !          potential_dot_acoustic,potential_dot_dot_acoustic,Mesh_pointer)
      endif

      ! this transfers fields only in elements with stations for efficiency
      if (ELASTIC_SIMULATION) then
        if (USE_CUDA_SEISMOGRAMS) then
          call transfer_seismograms_el_from_d(nrec_local,Mesh_pointer, &
                                             seismograms_d,seismograms_v,seismograms_a, &
                                             it)
        else
          call transfer_station_el_from_device(displ,veloc,accel, &
                                              b_displ,b_veloc,b_accel, &
                                              Mesh_pointer,number_receiver_global, &
                                              ispec_selected_rec,ispec_selected_source, &
                                              ibool)
        endif
        ! alternative: transfers whole fields
        !  call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc, accel, Mesh_pointer)
      endif
    endif
  endif

  if (.not. GPU_MODE .or. (GPU_MODE .and. (.not. USE_CUDA_SEISMOGRAMS))) then

    do irec_local = 1,nrec_local

      ! initializes wavefield values
      dxd = ZERO
      dyd = ZERO
      dzd = ZERO

      vxd = ZERO
      vyd = ZERO
      vzd = ZERO

      axd = ZERO
      ayd = ZERO
      azd = ZERO

      pd  = ZERO

      ! gets local receiver interpolators
      ! (1-D Lagrange interpolators)
      hxir(:) = hxir_store(irec_local,:)
      hetar(:) = hetar_store(irec_local,:)
      hgammar(:) = hgammar_store(irec_local,:)

      ! gets global number of that receiver
      irec = number_receiver_global(irec_local)

      ! spectral element in which the receiver is located
      if (SIMULATION_TYPE == 2) then
        ! adjoint "receivers" are located at CMT source positions
        ! note: we take here xi_source,.. when FASTER_RECEIVERS_POINTS_ONLY is set
        ispec = ispec_selected_source(irec)
        xi_r = xi_source(irec)
        eta_r = eta_source(irec)
        gamma_r = gamma_source(irec)
      else
        ! receiver located at station positions
        ispec = ispec_selected_rec(irec)
        xi_r = xi_receiver(irec)
        eta_r = eta_receiver(irec)
        gamma_r = gamma_receiver(irec)
      endif

      ! calculates interpolated wavefield values at receiver positions
      select case (SIMULATION_TYPE)
      case (1,2)
        ! forward simulations & pure adjoint simulations
        ! wavefields stored in displ,veloc,accel

        ! elastic wave field
        if (ispec_is_elastic(ispec)) then
          ! interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva(displ,veloc,accel,NGLOB_AB, &
                                        ispec,NSPEC_AB,ibool, &
                                        xi_r,eta_r,gamma_r, &
                                        hxir,hetar,hgammar, &
                                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)
        endif ! elastic

        ! acoustic wave field
        if (ispec_is_acoustic(ispec)) then
          ! displacement vector
          call compute_gradient_in_acoustic(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_acoustic,displ_element, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)

          ! velocity vector
          call compute_gradient_in_acoustic(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic,veloc_element, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)

          ! acceleration vector
          call compute_gradient_in_acoustic(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_dot_acoustic,accel_element, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)

          ! interpolates displ/veloc/accel/pressure at receiver locations
          call compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                                               potential_dot_dot_acoustic,potential_acoustic,NGLOB_AB, &
                                               ispec,NSPEC_AB,ibool, &
                                               xi_r,eta_r,gamma_r, &
                                               hxir,hetar,hgammar, &
                                               dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd,USE_TRICK_FOR_BETTER_PRESSURE)
        endif ! acoustic

        ! poroelastic wave field
        if (ispec_is_poroelastic(ispec)) then
          ! interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva(displs_poroelastic,velocs_poroelastic,accels_poroelastic,NGLOB_AB, &
                                        ispec,NSPEC_AB,ibool, &
                                        xi_r,eta_r,gamma_r, &
                                        hxir,hetar,hgammar, &
                                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)
        endif ! poroelastic

      case (3)
        ! adjoint/kernel simulations
        ! reconstructed forward wavefield stored in b_displ, b_veloc, b_accel

        ! elastic wave field
        if (ispec_is_elastic(ispec)) then
          ! backward field: interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva(b_displ,b_veloc,b_accel,NGLOB_ADJOINT, &
                                        ispec,NSPEC_AB,ibool, &
                                        xi_r,eta_r,gamma_r, &
                                        hxir,hetar,hgammar, &
                                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)
        endif ! elastic

        ! acoustic wave field
        if (ispec_is_acoustic(ispec)) then
          ! backward field: displacement vector
          call compute_gradient_in_acoustic(ispec,NSPEC_AB,NGLOB_ADJOINT, &
                          b_potential_acoustic,displ_element, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)

          ! backward field: velocity vector
          call compute_gradient_in_acoustic(ispec,NSPEC_AB,NGLOB_ADJOINT, &
                          b_potential_dot_acoustic,veloc_element, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)

          ! backward field: acceleration vector
          call compute_gradient_in_acoustic(ispec,NSPEC_AB,NGLOB_AB, &
                          b_potential_dot_dot_acoustic,accel_element, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)

          ! backward field: interpolates displ/veloc/accel/pressure at receiver locations
          call compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                                               b_potential_dot_dot_acoustic,b_potential_acoustic,NGLOB_ADJOINT, &
                                               ispec,NSPEC_AB,ibool, &
                                               xi_r,eta_r,gamma_r, &
                                               hxir,hetar,hgammar, &
                                               dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd,USE_TRICK_FOR_BETTER_PRESSURE)
        endif ! acoustic

      end select ! SIMULATION_TYPE

      ! additional calculations for pure adjoint simulations
      ! computes derivatives of source parameters
      if (SIMULATION_TYPE == 2) then

        ! elastic wave field
        if (ispec_is_elastic(ispec)) then
          ! stores elements displacement field
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)
                displ_element(:,i,j,k) = displ(:,iglob)
              enddo
            enddo
          enddo

          ! gets derivatives of local receiver interpolators
          hpxir(:) = hpxir_store(irec_local,:)
          hpetar(:) = hpetar_store(irec_local,:)
          hpgammar(:) = hpgammar_store(irec_local,:)

          ! computes the integrated derivatives of source parameters (M_jk and X_s)
          call compute_adj_source_frechet(displ_element,Mxx(irec),Myy(irec),Mzz(irec), &
                                          Mxy(irec),Mxz(irec),Myz(irec),eps_s,eps_m_s, &
                                          hxir,hetar,hgammar,hpxir,hpetar,hpgammar, &
                                          hprime_xx,hprime_yy,hprime_zz, &
                                          xix(:,:,:,ispec),xiy(:,:,:,ispec),xiz(:,:,:,ispec), &
                                          etax(:,:,:,ispec),etay(:,:,:,ispec),etaz(:,:,:,ispec), &
                                          gammax(:,:,:,ispec),gammay(:,:,:,ispec),gammaz(:,:,:,ispec))

          stf = comp_source_time_function(dble(NSTEP-it)*DT-t0-tshift_src(irec),hdur_Gaussian(irec))
          stf_deltat = stf * deltat

          Mxx_der(irec_local) = Mxx_der(irec_local) + eps_s(1,1) * stf_deltat
          Myy_der(irec_local) = Myy_der(irec_local) + eps_s(2,2) * stf_deltat
          Mzz_der(irec_local) = Mzz_der(irec_local) + eps_s(3,3) * stf_deltat
          Mxy_der(irec_local) = Mxy_der(irec_local) + 2 * eps_s(1,2) * stf_deltat
          Mxz_der(irec_local) = Mxz_der(irec_local) + 2 * eps_s(1,3) * stf_deltat
          Myz_der(irec_local) = Myz_der(irec_local) + 2 * eps_s(2,3) * stf_deltat

          sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_s(:) * stf_deltat
        endif ! elastic
      endif

      ! store North, East and Vertical components
      if (SIMULATION_TYPE == 2) then
        ! adjoint simulations
        ! adjoint "receiver" N/E/Z orientations given by nu_source array
        if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
          seismograms_d(:,irec_local,it) = real(nu_source(:,1,irec)*dxd &
                                              + nu_source(:,2,irec)*dyd &
                                              + nu_source(:,3,irec)*dzd,kind=CUSTOM_REAL)

        if (SAVE_SEISMOGRAMS_VELOCITY) &
          seismograms_v(:,irec_local,it) = real(nu_source(:,1,irec)*vxd &
                                              + nu_source(:,2,irec)*vyd &
                                              + nu_source(:,3,irec)*vzd,kind=CUSTOM_REAL)

        if (SAVE_SEISMOGRAMS_ACCELERATION) &
          seismograms_a(:,irec_local,it) = real(nu_source(:,1,irec)*axd &
                                              + nu_source(:,2,irec)*ayd &
                                              + nu_source(:,3,irec)*azd,kind=CUSTOM_REAL)
      else
        ! forward & kernel simulations
        if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
          seismograms_d(:,irec_local,it) = real(nu(:,1,irec)*dxd + nu(:,2,irec)*dyd + nu(:,3,irec)*dzd,kind=CUSTOM_REAL)

        if (SAVE_SEISMOGRAMS_VELOCITY) &
          seismograms_v(:,irec_local,it) = real(nu(:,1,irec)*vxd + nu(:,2,irec)*vyd + nu(:,3,irec)*vzd,kind=CUSTOM_REAL)

        if (SAVE_SEISMOGRAMS_ACCELERATION) &
          seismograms_a(:,irec_local,it) = real(nu(:,1,irec)*axd + nu(:,2,irec)*ayd + nu(:,3,irec)*azd,kind=CUSTOM_REAL)
      endif

      ! only one scalar in the case of pressure
      if (SAVE_SEISMOGRAMS_PRESSURE) &
        seismograms_p(1,irec_local,it) = real(pd,kind=CUSTOM_REAL)

      ! adjoint simulations
      if (SIMULATION_TYPE == 2) seismograms_eps(:,:,irec_local,it) = eps_s(:,:)

    enddo ! nrec_local

  endif

  ! write the current or final seismograms
  if ((mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) .and. .not. SU_FORMAT) then
    if (SIMULATION_TYPE == 2) then
      ! adjoint simulations
      if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
        call write_adj_seismograms_to_file(myrank,seismograms_d,number_receiver_global,nrec_local,it,DT,NSTEP,t0,1)
      if (SAVE_SEISMOGRAMS_VELOCITY) &
        call write_adj_seismograms_to_file(myrank,seismograms_v,number_receiver_global,nrec_local,it,DT,NSTEP,t0,2)
      if (SAVE_SEISMOGRAMS_ACCELERATION) &
        call write_adj_seismograms_to_file(myrank,seismograms_a,number_receiver_global,nrec_local,it,DT,NSTEP,t0,3)
      if (SAVE_SEISMOGRAMS_PRESSURE) &
        call write_adj_seismograms_to_file(myrank,seismograms_p,number_receiver_global,nrec_local,it,DT,NSTEP,t0,4)
    else
      ! forward & kernel simulations
      if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
        call write_seismograms_to_file(seismograms_d,1)
      if (SAVE_SEISMOGRAMS_VELOCITY) &
        call write_seismograms_to_file(seismograms_v,2)
      if (SAVE_SEISMOGRAMS_ACCELERATION) &
        call write_seismograms_to_file(seismograms_a,3)
      if (SAVE_SEISMOGRAMS_PRESSURE) &
        call write_seismograms_to_file(seismograms_p,4)
    endif
  endif

  ! write ONE binary file for all receivers (nrec_local) within one proc
  ! SU format, with 240-byte-header for each trace
  if ((mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) .and. SU_FORMAT) &
    call write_output_SU()

  end subroutine write_seismograms


!================================================================

! write seismograms to text files

  subroutine write_seismograms_to_file(seismograms,istore)

  use constants

  use specfem_par, only: myrank,number_receiver_global,station_name,network_name, &
          nrec,nrec_local,islice_selected_rec, &
          it,DT,NSTEP,t0,SIMULATION_TYPE,WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_SEISMOGRAMS

  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms

  ! local parameters
  integer irec,irec_local

  character(len=1) component

  ! parameters for master collects seismograms
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram
  integer :: nrec_local_received,NPROCTOT,total_seismos,receiver,sender
  integer :: iproc,ier
  integer,dimension(1) :: tmp_nrec_local_received,tmp_irec,tmp_nrec_local
  integer,dimension(:),allocatable:: islice_num_rec_local

  character(len=MAX_STRING_LEN) :: sisname

  ! saves displacement, velocity, acceleration, or pressure
  if (istore == 1) then
    component = 'd'
  else if (istore == 2) then
    component = 'v'
  else if (istore == 3) then
    component = 'a'
  else if (istore == 4) then
    component = 'p'
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

  allocate(one_seismogram(NDIM,NSTEP),stat=ier)
  if (ier /= 0) stop 'error while allocating one temporary seismogram'

  ! write out seismograms: all processes write their local seismograms themselves
  if (.not. WRITE_SEISMOGRAMS_BY_MASTER) then

    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      write(sisname,'(A,I5.5)') '/all_seismograms_node_',myrank
      if (USE_BINARY_FOR_SEISMOGRAMS) then
        open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown',form='unformatted',action='write')
      else
        open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
      endif
    endif

    ! loop on all the local receivers
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      ! writes out this seismogram
      one_seismogram = seismograms(:,irec_local,:)

      call write_one_seismogram(one_seismogram,irec, &
                                station_name,network_name,nrec, &
                                DT,t0,it,NSTEP,SIMULATION_TYPE, &
                                myrank,component,istore)

    enddo ! nrec_local

    ! create one large file instead of one small file per station to avoid file system overload
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

! only the master process does the writing of seismograms and
! collects the data from all other processes
  else ! if WRITE_SEISMOGRAMS_BY_MASTER

    if (myrank == 0) then
      ! on the master, gather all the seismograms

      ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
        write(sisname,'(A)') '/all_seismograms'
        if (USE_BINARY_FOR_SEISMOGRAMS) then
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown',form='unformatted',action='write')
        else
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
        endif
      endif

      total_seismos = 0

      ! loop on all the slices
      call world_size(NPROCTOT)

      ! counts number of local receivers for each slice
      allocate(islice_num_rec_local(0:NPROCTOT-1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'error allocating islice_num_rec_local')
      islice_num_rec_local(:) = 0
      do irec = 1,nrec
        iproc = islice_selected_rec(irec)
        islice_num_rec_local(iproc) = islice_num_rec_local(iproc) + 1
      enddo

      ! loop on all the slices
      do iproc = 0,NPROCTOT-1

        ! communicate only with processes which contain local receivers
        if (islice_num_rec_local(iproc) == 0) cycle

        ! receive except from proc 0, which is me and therefore I already have this value
        sender = iproc
        if (iproc /= 0) then
          call recv_i(tmp_nrec_local_received,1,sender,itag)
          nrec_local_received = tmp_nrec_local_received(1)
          if (nrec_local_received < 0) call exit_MPI(myrank,'error while receiving local number of receivers')
        else
          nrec_local_received = nrec_local
        endif

        if (nrec_local_received > 0) then
          do irec_local = 1,nrec_local_received
            ! receive except from proc 0, which is myself and therefore I already have these values
            if (iproc == 0) then
              ! get global number of that receiver
              irec = number_receiver_global(irec_local)
              one_seismogram(:,:) = seismograms(:,irec_local,:)
            else
              call recv_i(tmp_irec,1,sender,itag)
              irec = tmp_irec(1)
              if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')

              call recvv_cr(one_seismogram,NDIM*NSTEP,sender,itag)
            endif

            total_seismos = total_seismos + 1

            ! writes out this seismogram
            call write_one_seismogram(one_seismogram,irec, &
                                      station_name,network_name,nrec, &
                                      DT,t0,it,NSTEP,SIMULATION_TYPE, &
                                      myrank,component,istore)

          enddo ! nrec_local_received
        endif ! if (nrec_local_received > 0)
      enddo ! NPROCTOT-1
      deallocate(islice_num_rec_local)

      write(IMAIN,*) 'Component: .sem'//component
      write(IMAIN,*) '  total number of receivers saved is ',total_seismos,' out of ',nrec
      write(IMAIN,*)

      if (total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

     ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    else
      ! on the nodes, send the seismograms to the master
      receiver = 0
      tmp_nrec_local(1) = nrec_local
      call send_i(tmp_nrec_local,1,receiver,itag)
      if (nrec_local > 0) then
        do irec_local = 1,nrec_local
          ! get global number of that receiver
          irec = number_receiver_global(irec_local)
          tmp_irec(1) = irec
          call send_i(tmp_irec,1,receiver,itag)

          ! sends seismogram of that receiver
          one_seismogram(:,:) = seismograms(:,irec_local,:)
          call sendv_cr(one_seismogram,NDIM*NSTEP,receiver,itag)
        enddo
      endif
    endif ! myrank

  endif ! of if (WRITE_SEISMOGRAMS_BY_MASTER)

  deallocate(one_seismogram)

  end subroutine write_seismograms_to_file

!=====================================================================

  subroutine write_one_seismogram(one_seismogram,irec, &
              station_name,network_name,nrec, &
              DT,t0,it,NSTEP,SIMULATION_TYPE, &
              myrank,component,istore)

  use constants

  implicit none

  integer, intent(in) :: NSTEP,it,SIMULATION_TYPE,istore
  real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP) :: one_seismogram

  integer myrank
  double precision t0,DT

  integer :: nrec,irec
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name
  character(len=1) component

  ! local parameters
  integer :: iorientation,number_of_components
  integer :: length_station_name,length_network_name
  character(len=MAX_STRING_LEN) :: sisname,final_LOCAL_PATH
  character(len=3) :: channel

  ! see how many components we need to store: 1 for pressure, NDIM for a vector
  if (istore == 4) then ! this is for pressure
    number_of_components = 1
  else
    number_of_components = NDIM
  endif

  ! loop over each seismogram component
  do iorientation = 1,number_of_components

    ! gets channel name
    if (istore == 4) then ! this is for pressure
      call write_channel_name(istore,channel)
    else
      call write_channel_name(iorientation,channel)
    endif

    ! create the name of the seismogram file for each slice
    ! file name includes the name of the station, the network and the component
    length_station_name = len_trim(station_name(irec))
    length_network_name = len_trim(network_name(irec))

    ! check that length conforms to standard
    if (length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
      call exit_MPI(myrank,'wrong length of station name')

    if (length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
      call exit_MPI(myrank,'wrong length of network name')

    ! writes out **net**.**sta**.**BH**.sem* files
    write(sisname,"(a,'.',a,'.',a3,'.sem',a1)") network_name(irec)(1:length_network_name), &
       station_name(irec)(1:length_station_name),channel,component

    ! directory to store seismograms
    final_LOCAL_PATH = OUTPUT_FILES(1:len_trim(OUTPUT_FILES)) // '/'

    ! ASCII output format
    call write_output_ASCII_or_binary(one_seismogram, &
                                      NSTEP,it,SIMULATION_TYPE,DT,t0, &
                                      iorientation,sisname,final_LOCAL_PATH)

  enddo ! do iorientation

  end subroutine write_one_seismogram

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms_to_file(myrank,seismograms,number_receiver_global, &
                                           nrec_local,it,DT,NSTEP,t0,istore)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  implicit none

  integer :: myrank
  integer :: nrec_local,NSTEP,it,istore
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
  double precision :: t0,DT

  ! local parameters
  integer :: irec,irec_local
  integer :: iorientation,isample,number_of_components

  character(len=3) :: channel
  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

  ! saves displacement, velocity, acceleration, or pressure
  if (istore == 1) then
    component = 'd'
  else if (istore == 2) then
    component = 'v'
  else if (istore == 3) then
    component = 'a'
  else if (istore == 4) then
    component = 'p'
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! see how many components we need to store: 1 for pressure, NDIM for a vector
    if (istore == 4) then ! this is for pressure
      number_of_components = 1
    else
      number_of_components = NDIM
    endif

    ! loop over each seismogram component
    do iorientation = 1,number_of_components

      ! gets channel name
      if (istore == 4) then ! this is for pressure
        call write_channel_name(istore,channel)
      else
        call write_channel_name(iorientation,channel)
      endif

      ! create the name of the seismogram file for each slice
      ! file name includes the name of the station, the network and the component
      write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,channel,component

      ! save seismograms in text format with no subsampling.
      ! Because we do not subsample the output, this can result in large files
      ! if the simulation uses many time steps. However, subsampling the output
      ! here would result in a loss of accuracy when one later convolves
      ! the results with the source time function
      open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)),status='unknown')

      ! make sure we never write more than the maximum number of time steps
      ! subtract half duration of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
        ! distinguish between single and double precision for reals
        write(IOUT,*) real(dble(isample-1)*DT - t0,kind=CUSTOM_REAL),' ',seismograms(iorientation,irec_local,isample)
      enddo

      close(IOUT)

    enddo

  enddo

  end subroutine write_adj_seismograms_to_file

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2_to_file(myrank,seismograms,number_receiver_global,nrec_local,it,DT,NSTEP,t0)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  implicit none
  integer :: myrank
  integer :: nrec_local,NSTEP,it
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local,NSTEP) :: seismograms
  double precision :: t0,DT

  ! local parameters
  integer :: irec,irec_local
  integer :: idimval,jdimval,isample

  character(len=4) :: chn
  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

  component = 'd'

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do idimval = 1,NDIM
      do jdimval = idimval,NDIM

        ! strain channel name
        if (idimval == 1 .and. jdimval == 1) then
          chn = 'SNN'
        else if (idimval == 1 .and. jdimval == 2) then
          chn = 'SEN'
        else if (idimval == 1 .and. jdimval == 3) then
          chn = 'SEZ'
        else if (idimval == 2 .and. jdimval == 2) then
          chn = 'SEE'
        else if (idimval == 2 .and. jdimval == 3) then
          chn = 'SNZ'
        else if (idimval == 3 .and. jdimval == 3) then
          chn = 'SZZ'
        else
          call exit_MPI(myrank,'incorrect channel value')
        endif

        ! create the name of the seismogram file for each slice
        ! file name includes the name of the station, the network and the component
        write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,chn,component

        ! save seismograms in text format with no subsampling.
        ! Because we do not subsample the output, this can result in large files
        ! if the simulation uses many time steps. However, subsampling the output
        ! here would result in a loss of accuracy when one later convolves
        ! the results with the source time function
        open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)),status='unknown')

        ! make sure we never write more than the maximum number of time steps
        ! subtract half duration of the source to make sure travel time is correct
        do isample = 1,min(it,NSTEP)
          ! distinguish between single and double precision for reals
          write(IOUT,*) real(dble(isample-1)*DT - t0,kind=CUSTOM_REAL),' ',seismograms(jdimval,idimval,irec_local,isample)
        enddo

        close(IOUT)

      enddo ! jdimval
    enddo ! idimval
  enddo ! irec_local

  end subroutine write_adj_seismograms2_to_file

!=====================================================================

  subroutine write_channel_name(iorientation,channel)

  use specfem_par, only: DT,SUPPRESS_UTM_PROJECTION

  implicit none

  integer :: iorientation
  character(len=3) :: channel

  ! local parameters
  character(len=2) :: bic
  double precision:: sampling_rate

  ! gets band and instrument code
  sampling_rate = DT
  call band_instrument_code(sampling_rate,bic)

  ! sets channel name
  if (SUPPRESS_UTM_PROJECTION) then

    ! no UTM, pure Cartesian reference
    ! uses Cartesian X/Y/Z direction to denote channel
    select case (iorientation)
    case (1)
      channel = bic(1:2)//'X'
    case (2)
      channel = bic(1:2)//'Y'
    case (3)
      channel = bic(1:2)//'Z'
    case (4)
      channel = bic(1:2)//'P'  ! for pressure seismograms
    case default
      call exit_mpi(0,'error channel orientation value')
    end select

  else

    ! UTM conversion
    ! uses convention for N/E/Z to denote channel
    select case (iorientation)
    case (1)
      channel = bic(1:2)//'E'
    case (2)
      channel = bic(1:2)//'N'
    case (3)
      channel = bic(1:2)//'Z'
    case (4)
      channel = bic(1:2)//'P'  ! for pressure seismograms
    case default
      call exit_mpi(0,'error channel orientation value')
    end select

  endif

  end subroutine write_channel_name

!=====================================================================

  subroutine band_instrument_code(DT,bic)
  ! This subroutine is to choose the appropriate band and instrument codes for channel names of seismograms
  ! based on the IRIS convention (first two letters of channel codes, respectively,
  ! which were LH(Z/E/N) previously).
  ! For consistency with observed data, we now use the IRIS convention for band codes (first letter in channel codes) of
  ! SEM seismograms governed by their sampling rate.
  ! Instrument code (second letter in channel codes) is fixed to "X" which is assigned by IRIS for synthetic seismograms.
  ! See the manual for further explanations!
  ! Ebru Bozdag, November 2010
  implicit none
  double precision :: DT
  character(len=2) :: bic
  ! local parameter
  logical,parameter :: SUPPRESS_IRIS_CONVENTION = .false.

  ! see manual for ranges
  if (DT >= 1.0d0)  bic = 'LX'
  if (DT < 1.0d0 .and. DT > 0.1d0) bic = 'MX'
  if (DT <= 0.1d0 .and. DT > 0.0125d0) bic = 'BX'
  if (DT <= 0.0125d0 .and. DT > 0.004d0) bic = 'HX'
  if (DT <= 0.004d0 .and. DT > 0.001d0) bic = 'CX'
  if (DT <= 0.001d0) bic = 'FX'

  ! ignores IRIS convention, uses previous, constant band and instrument code
  if (SUPPRESS_IRIS_CONVENTION) then
    bic = 'BH'
  endif

 end subroutine band_instrument_code

