module adjoint_source
   use specfem_par
   use dsm_coupling
   use project_tomo_grid_mod, only : nx,ny,nz
   implicit none
   !include 'mpif.h'
   !include 'precision.h'
   integer NT_decimed,under_sampling,d_shift_tapper,nbadj
   real, allocatable :: all_sismograms(:,:,:,:)
   real(kind=CUSTOM_REAL), allocatable :: trav_time_ccZ(:,:,:),trav_ttZ(:,:,:),trav_time_ccX(:,:,:),trav_ttX(:,:,:)
   real(kind=CUSTOM_REAL), allocatable :: adj_src_TT(:,:,:),adj_src_all(:,:,:),seism_data(:,:,:,:),pick_data(:,:),sismo_buffer(:,:,:)
   real(kind=CUSTOM_REAL),allocatable  :: tr_data(:),tr_synth(:),tr_adj_src(:),t_w(:)
   real(kind=CUSTOM_REAL) sigma_data
   character(len=250), allocatable :: data_dir_name(:)
   character(len=250), allocatable :: data_pick_file(:)
   real(kind=CUSTOM_REAL),allocatable :: lwind(:),rwind(:),fr_low(:),fr_high(:)
   integer, allocatable :: icomp_to_invert(:),tapper_index(:,:,:)
   real(kind=CUSTOM_REAL) reg_par_1,reg_par_2,fitting_value,data_adjustment
   real(kind=CUSTOM_REAL) length_cor,length_cor0,sigma_model0
   character(len=2) reg_type
   logical READ_MODEL_IN_TOMO_GRID
   character(len=250) file_model_vp,file_model_vs,file_model_rho
   character(len=250) file_model_vp_prior,file_model_vs_prior,file_model_rho_prior
   integer, allocatable :: iflag_seism_data(:,:)
   !! inversion parameter
   integer niter
   character(len=100) my_debbug

  contains

    subroutine initialize_inversion()

      implicit none
      call read_inversion_parameters()
      call initialize_all_seismograms(NSTEP,nrec_local,nb_DSM_sources)
      call read_data_seism()
      call set_up_adjoint_sources_wksp()

    end subroutine initialize_inversion


    subroutine read_inversion_parameters()

      implicit none
      integer ier,i
      integer source_type
      character(len=300) sem_dir_save,sem_dir_read,tmp_char

      allocate(data_dir_name(nb_DSM_sources))
      allocate(data_pick_file(nb_DSM_sources))
      allocate(lwind(nb_DSM_sources))
      allocate(rwind(nb_DSM_sources))
      allocate(icomp_to_invert(nb_DSM_sources))
      allocate(fr_low(nb_DSM_sources),fr_high(nb_DSM_sources))

      if (myrank==0) then
        open(27,file='inversion.par')

        read(27,*) nb_DSM_sources
        write(*,*) nb_DSM_sources

        read(27,*) source_type
        write(*,*) source_type

        read(27,*) niter ! nb iterations
        write(*,*) niter

        do i=1,nb_DSM_sources

          read(27,'(a)') data_dir_name(i)  !! repertoire des data
          write(*,*) trim(data_dir_name(i))

          read(27,'(a)') data_pick_file(i) !! fichier pick
          write(*,*) trim(data_pick_file(i))

          read(27,*) lwind(i),rwind(i)     !! fenetre / pick
          write(*,*) lwind(i),rwind(i)

          read(27,*) fr_low(i),fr_high(i)  !! filtre
          write(*,*) fr_low(i),fr_high(i)

          read(27,*) icomp_to_invert(i)    !! composante a inverser
          write(*,*) icomp_to_invert(i)

          !write(*,'(a)') 'data :',trim(data_dir_name(i))

        enddo

        read(27,'(a2)') reg_type ! type de regularisation
        write(*,*) trim(reg_type)
        read(27,*) reg_par_1     ! coeff 1
        write(*,*) reg_par_1

        read(27,*) reg_par_2     ! coeff 2
        write(*,*) reg_par_2

        read(27,*) READ_MODEL_IN_TOMO_GRID
        write(*,*) READ_MODEL_IN_TOMO_GRID

        if (READ_MODEL_IN_TOMO_GRID) then
          read(27,'(a)') file_model_vp
          read(27,'(a)') file_model_vs
          read(27,'(a)') file_model_rho
        endif

        read(27,*) use_precond
        write(*,*) use_precond

        read(27,*) iparam_inv  ! =1 => (vp,vs), =2 =>(vp)
        write(*,*) iparam_inv

        read(27,*) use_norma
        write(*,*) use_norma
        if (use_norma) read(27,*) sigma_data

        scale_cost = 1._CUSTOM_REAL
        scale_pena = 1._CUSTOM_REAL
        if (use_norma) then
           scale_cost=1._CUSTOM_REAL/(deltat*nstep*(3+nb_DSM_sources+nrec)*sigma_data**2)
           reg_par_1  = reg_par_1  / (nx*ny*nz)
           reg_par_2  = reg_par_2  / (nx*ny*nz)
        endif

        if(reg_type=='CO') then
           read(27,'(a)') file_model_vp_prior
           read(27,'(a)') file_model_vs_prior
           read(27,'(a)') file_model_rho_prior
           read(27,*) sigma_data
           read(27,*) length_cor,length_cor0
           read(27,*) sigma_model0
           scale_cost = 1._CUSTOM_REAL/(deltat*nstep*(3+nb_DSM_sources+nrec)*sigma_data**2)
           scale_pena =   1._CUSTOM_REAL/(nx*ny*nz)
           reg_par_1  =  sigma_model0 * length_cor0 / length_cor
           reg_par_2  =  length_cor
           write(*,*) scale_cost,scale_pena
           write(*,*) reg_par_1,reg_par_2
        endif
        read(27,'(a)') sem_dir_save
        read(27,'(a)') sem_dir_read
        close(27)

        if (source_type==0) then
          COUPLING_WITH_DSM=.false.
          COUPLING_WITH_PLAN=.false.
        endif
        if (source_type==1) then
          COUPLING_WITH_DSM=.true.
          COUPLING_WITH_PLAN=.false.
        endif
        if (source_type==2) then
           COUPLING_WITH_DSM=.false.
           COUPLING_WITH_PLAN=.true.
        endif
      endif

     call mpi_bcast( iparam_inv,1 ,MPI_INTEGER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( niter,1 ,MPI_INTEGER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( icomp_to_invert,nb_DSM_sources ,MPI_INTEGER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( lwind, nb_DSM_sources,CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( rwind, nb_DSM_sources,CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( fr_low, nb_DSM_sources,CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( fr_high, nb_DSM_sources,CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( reg_par_1, 1, CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( reg_par_2, 1, CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( reg_type,2,MPI_CHARACTER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( data_dir_name, 250*nb_DSM_sources,MPI_CHARACTER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( data_pick_file, 250*nb_DSM_sources,MPI_CHARACTER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast(use_precond, 1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
     call mpi_bcast(COUPLING_WITH_DSM, 1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
     call mpi_bcast(COUPLING_WITH_PLAN, 1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
     call mpi_bcast(READ_MODEL_IN_TOMO_GRID, 1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( scale_cost, 1, CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( scale_pena, 1, CUSTOM_MPI_TYPE ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( sem_dir_save ,300 ,MPI_CHARACTER ,0,MPI_COMM_WORLD,ier)
     call mpi_bcast( sem_dir_read ,300 ,MPI_CHARACTER ,0,MPI_COMM_WORLD,ier)

     write(tmp_char,'(a13,i6.6)') 'model_saved_p',myrank
     SEM_FILE_OUTPUT_MODEL=trim(sem_dir_save)//trim(tmp_char)
     SEM_FILE_INTPUT_MODEL=trim(sem_dir_read)//trim(tmp_char)

     if (icomp_to_invert(1) == 6 ) then
       nbadj=2
     else
       nbadj=1
     endif
    end subroutine read_inversion_parameters



    subroutine read_data_seism()

      implicit none
      integer isrc,ista,irec_local,nbsrc,nbsta,icomp
      integer iti,ier
      real(kind=CUSTOM_REAL) xx,tt
      character(len=250) data_file
      character(len=3) channel
      character(len=1) component
      character(len=256) sismane,record_name
      real(kind=CUSTOM_REAL), allocatable :: pick_read(:)
      integer indx00,indx11
      real w0,w1

      allocate(seism_data(nstep,3,nrec_local,nb_DSM_sources))
      allocate(pick_data(2,nrec_local),pick_read(nb_DSM_sources))
      allocate(iflag_seism_data(nrec,nb_DSM_sources))
      allocate(tapper_index(4,nrec,nb_DSM_sources))
      allocate(trav_time_ccZ(2,nrec,nb_DSM_sources),trav_time_ccX(2,nrec,nb_DSM_sources))
      allocate(trav_ttZ(2,nrec,nb_DSM_sources),trav_ttX(2,nrec,nb_DSM_sources))

      seism_data=0._CUSTOM_REAL
      iflag_seism_data=0
      !d_shift_tapper=40
      !component='d'
      !do isrc=1,nb_DSM_sources
      !   do irec_local=1,nrec_local
      !     ista = number_receiver_global(irec_local)
      !     do icomp = 1, 3

      !       call write_channel_name1(icomp,channel)
             !data_file=trim(data_dir_name(isrc))//'/'//trim(station_name(ista))//'.'//trim(network_name(ista))
             !data_file=trim(station_name(ista))//'.'//trim(network_name(ista))

      !       write(sismane,"(a,'.',a,'.',a3,'.sem',a1)") trim(station_name(ista)),&
      !                     trim(network_name(ista)),channel,component

     !        data_file=trim(data_dir_name(isrc))//'/'//trim(sismane)
             !write(*,'(i4,a)') myrank,trim(data_file)


             !!  reading seismogram data
     !        open(27,file=trim(data_file))
     !        do iti=1,nstep
     !            read(27,*) tt,xx
     !            seism_data(iti,icomp,irec_local,isrc) = xx
     !        enddo
     !        close(27)


      !     enddo
      !  enddo
      !enddo
     d_shift_tapper=40
      if (myrank==0) then
        !! lecture des pick
        do isrc=1,nb_DSM_sources
           open(27,file=trim(data_pick_file(isrc)))
           do
             read(27,*,end=99) record_name,w0,w1,indx00,indx11
             call find_ista(record_name,ista)
             iflag_seism_data(ista,isrc)=1
             !write(*,*) trim(record_name),w0,w1,indx00,indx11
             indx11=int(rwind(isrc)/dt) + indx00 ! nouvelle def de la fin de fenetre (Tw apres le pointe)
             !write(*,*) int(rwind(isrc)/dt), indx00 , rwind(isrc), dt
             !write(*,*) int(lwind(isrc)/dt), indx00, lwind(isrc),dt
             indx00 = indx00 - int(lwind(isrc)/dt)

              tapper_index(1,ista,isrc)=max(indx00 - d_shift_tapper,0)
              tapper_index(2,ista,isrc)=indx00
              tapper_index(3,ista,isrc)=indx11
              tapper_index(4,ista,isrc)=min(indx11 + d_shift_tapper,nstep)
              !write(*,'(i5,3x,2(f10.4,1x),f10.5,4i5)')indx00, lwind(isrc),rwind(isrc),dt,tapper_index(:,ista,isrc)
           enddo
           99 close(27)
        enddo

      endif
      call mpi_bcast(tapper_index,4*nrec*nb_DSM_sources,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call mpi_bcast(iflag_seism_data,nrec*nb_DSM_sources,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      !d_shift_tapper=40
      component='d'
      do isrc=1,nb_DSM_sources
         do irec_local=1,nrec_local
           ista = number_receiver_global(irec_local)
           if (iflag_seism_data(ista,isrc)==1) then
            do icomp = 1, 3

             call write_channel_name1(icomp,channel)
             !data_file=trim(data_dir_name(isrc))//'/'//trim(station_name(ista))//'.'//trim(network_name(ista))
             !data_file=trim(station_name(ista))//'.'//trim(network_name(ista))

             write(sismane,"(a,'.',a,'.',a3,'.sem',a1)") trim(station_name(ista)),&
                           trim(network_name(ista)),channel,component

             data_file=trim(data_dir_name(isrc))//'/'//trim(sismane)
             !write(*,'(i4,a)') myrank,trim(data_file)


             !!  reading seismogram data
             open(27,file=trim(data_file))
             do iti=1,nstep
                 read(27,*) tt,xx
                 seism_data(iti,icomp,irec_local,isrc) = xx
             enddo
             close(27)


            enddo
          endif
        enddo
      enddo




  end subroutine read_data_seism



  subroutine find_ista(record_name,ista)
    implicit none
    character(len=256) record_name,sta_code
    character(len=1) ch1
    integer i,ista
    i=1
    ch1=record_name(i:i)
    do while (ch1 /= '.')
       i=i+1
       ch1=record_name(i:i)
    enddo
    sta_code=record_name(1:i-1)


    do ista=1,nrec
       !write(*,*) ista,trim(sta_code),'  ', trim(station_name(ista))
       if (trim(sta_code) == trim(station_name(ista))) return
    enddo
    ista=0
    write(*,*)  'STOP , station ',trim(sta_code), ' Not Found'
    stop
  end subroutine find_ista



  subroutine initialize_all_seismograms(Nt,nbsta,nbsrc)
      implicit none
      integer Nt,nbsta,nbsrc
      under_sampling=100
      NT_decimed = Nt / under_sampling

      if (nbsta > 0) then
        allocate(all_sismograms(nt,3,nbsta,nbsrc))
        allocate(adj_src_all(Nt,3,nbsta))
        allocate(adj_src_TT(Nt,3,nbsta))
        allocate(sismo_buffer(Nt,3,nbsrc))
      else
        allocate(all_sismograms(1,1,1,1))
        allocate(adj_src_all(1,1,1))
        allocate(adj_src_TT(1,1,1))
        allocate(sismo_buffer(1,1,1))
      endif

   end subroutine initialize_all_seismograms


   subroutine store_all_seismograms(isrc,iter,nbsta,nbsrc,sism,nt)
      implicit none
      integer ista, itt,nt
      integer nbsta,nbsrc,isrc,iter
      real(kind=CUSTOM_REAL) sism(3,nbsta,nt)

      if (nbsta > 0) then

         do ista=1,nbsta
            do itt=1,nt
               all_sismograms(itt,1,ista,isrc) = sism(1,ista,itt)
               all_sismograms(itt,2,ista,isrc) = sism(2,ista,itt)
               all_sismograms(itt,3,ista,isrc) = sism(3,ista,itt)
            enddo
         enddo

      endif

    end subroutine store_all_seismograms

!------------------------------------------------------------------

!------------------------------------------------------------------
   subroutine write_sismo_in_disk(iter)
     implicit none
     integer iter,irec,ista,icomp,irec_local,isrc
     character(len=10) station_file
     irec=iter+1
     !write(*,*) 'myrank ', irec
     do irec_local=1,nrec_local
         ista = number_receiver_global(irec_local)
         do isrc=1,nb_DSM_sources
            do icomp=1,3
               sismo_buffer(:,icomp,isrc)=all_sismograms(:,icomp,irec_local,isrc)
            enddo
         enddo
         write(station_file,"(a,'.',a3)") trim(station_name(ista)),'bin'
         open(27,file=trim(station_file),access='direct',recl=CUSTOM_REAL*3*nstep*nb_DSM_sources)
         write(27,rec=irec) sismo_buffer
         close(27)
     enddo



   end subroutine write_sismo_in_disk

!------------------------------------------------------------------

    subroutine set_up_adjoint_sources_wksp()

      implicit none

      allocate(tr_data(nstep),tr_synth(nstep),tr_adj_src(nstep),t_w(nstep))


    end subroutine set_up_adjoint_sources_wksp

!-----------------------------------------------------------------
    subroutine store_adj_src(isrc,iadj)
      implicit none
      integer isrc,iadj


      if (icomp_to_invert(isrc) < 5) return
      if (icomp_to_invert(isrc) == 5) then

         adj_src_all(:,:,:) = adj_src_TT(:,:,:)

      endif

     if (icomp_to_invert(isrc) == 6 .and. iadj==1) then
        adj_src_all(:,:,:) = 0.
        adj_src_all(:,1,:) = adj_src_TT(:,1,:)
     endif

     if (icomp_to_invert(isrc) == 6 .and. iadj==2) then
        adj_src_all(:,:,:) = 0.
        adj_src_all(:,3,:) = adj_src_TT(:,3,:)
     endif

   end subroutine store_adj_src
!-----------------------------------------------------------------
    subroutine define_adjoint_sources(isrc)

     implicit none
     integer ierr
     integer i,i0,i1,i2,i3
     integer ishift_time
     integer norder,irek,npts
     integer isrc,irec_local,ista,icomp,nbadj
     real(kind=CUSTOM_REAL) Weigth,ddt,adjust,adjust_time,tshift,ccm
     !character(len=100) my_debbug

     ddt=dt
     Weigth=1.
     norder=4
     irek=1
     nbadj=1
     npts=nstep
     adjust=0.
     adjust_time=0.
     trav_time_ccZ(:,:,isrc)=0.
     trav_time_ccX(:,:,isrc)=0.
     adj_src_all=0._CUSTOM_REAL
     adj_src_TT=0._CUSTOM_REAL

     do irec_local=1,nrec_local

          ista = number_receiver_global(irec_local)
          !write(*,*) 'myrank ', myrank, ista,isrc
          i0=tapper_index(1,ista,isrc)
          i1=tapper_index(2,ista,isrc)
          i2=tapper_index(3,ista,isrc)
          i3=tapper_index(4,ista,isrc)

          !write(*,*) myrank,i0,i1,i2,i3

          !calcul du tapper
          call tapper_window_W(t_w,i0,i1,i2,i3,nstep,Weigth)

   if (iflag_seism_data(ista,isrc)==1)  then

       do icomp = 1, 3


           tr_data(:)    = seism_data(:,icomp,irec_local,isrc)
           tr_synth(:)   = seismograms_d(icomp,irec_local,:)
           tr_adj_src(:) = t_w(:)*(tr_synth(:) - tr_data(:))
           !write(my_debbug,'(a4,i4.4,i4.4)')'nof_',ista,isrc
           !open(666,file=trim(my_debbug))
           !write(666,*) tr_adj_src(:) !adj_src_all(:,icomp,irec_local)
           !close(666)
           tr_synth(:) = tr_adj_src(:)
           call  bwfilt (tr_synth, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))

           adj_src_all(:,icomp,irec_local)=tr_adj_src(:)

           !write(my_debbug,'(a4,i4.4,i4.4)')'adj_',ista,isrc
           !open(666,file=trim(my_debbug))
           !write(666,*) ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc)
           !write(666,*) adj_src_all(:,icomp,irec_local)
           !close(666)
       enddo

      if (icomp_to_invert(isrc) == 1) then
          adj_src_all(:,2,irec_local) = 0.
          adj_src_all(:,3,irec_local) = 0.
      endif

      if (icomp_to_invert(isrc) == 2) then
          adj_src_all(:,1,irec_local) = 0.
          adj_src_all(:,3,irec_local) = 0.
      endif

      if (icomp_to_invert(isrc) == 3) then
          adj_src_all(:,1,irec_local) = 0.
          adj_src_all(:,2,irec_local) = 0.
      endif


     adjust = adjust + sum(adj_src_all(:,1,irec_local)**2) + &
                       sum(adj_src_all(:,2,irec_local)**2) + &
                       sum(adj_src_all(:,3,irec_local)**2)



     if (icomp_to_invert(isrc) == 5) then !! adjoint tomo for Z
         icomp=3
         tr_data(:)    = t_w(:)*seism_data(:,icomp,irec_local,isrc)
         tr_synth(:)   = t_w(:)*seismograms_d(icomp,irec_local,:)

         call  bwfilt (tr_synth, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))
         tr_synth(:) = tr_adj_src(:)

         call  bwfilt (tr_data, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))
         tr_data(:) = tr_adj_src(:)

         call  calc_max_corr(tr_data, tr_synth, npts, i0, i3, ishift_time,ccm)

         adj_src_TT(:,1,irec_local) = 0.
         adj_src_TT(:,2,irec_local) = 0.
         tshift = -ishift_time*ddt
         call calc_adj_tr(tshift,tr_synth, tr_adj_src, tr_data, npts, ddt)
         adj_src_TT(:,icomp,irec_local)=tr_adj_src(:)

         adjust_time = adjust_time + 0.5 * tshift**2
         trav_time_ccZ(1,ista,isrc)=tshift
         trav_time_ccZ(2,ista,isrc)=ccm

     endif

     if (icomp_to_invert(isrc) == 6) then !! adjoint tomo for comp X and Z
         nbadj=2

         adj_src_TT(:,1,irec_local) = 0.

         icomp=1
         tr_data(:)    = t_w(:)*seism_data(:,icomp,irec_local,isrc)
         tr_synth(:)   = t_w(:)*seismograms_d(icomp,irec_local,:)

         call  bwfilt (tr_synth, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))
         tr_synth(:) = tr_adj_src(:)

         call  bwfilt (tr_data, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))
         tr_data(:) = tr_adj_src(:)

         call  calc_max_corr(tr_data, tr_synth, npts, i0, i3, ishift_time,ccm)

         tshift = -ishift_time*ddt
         call calc_adj_tr(tshift,tr_synth, tr_adj_src, tr_data, npts, ddt)
         adj_src_TT(:,icomp,irec_local)=tr_adj_src(:)

         adjust_time = adjust_time + 0.5 * tshift**2
         trav_time_ccX(1,ista,isrc)=tshift
         trav_time_ccX(2,ista,isrc)=ccm

         icomp=3
         tr_data(:)    = t_w(:)*seism_data(:,icomp,irec_local,isrc)
         tr_synth(:)   = t_w(:)*seismograms_d(icomp,irec_local,:)

         call  bwfilt (tr_synth, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))
         tr_synth(:) = tr_adj_src(:)

         call  bwfilt (tr_data, tr_adj_src, ddt, npts, irek, norder, fr_low(isrc),fr_high(isrc))
         tr_data(:) = tr_adj_src(:)

         call  calc_max_corr(tr_data, tr_synth, npts, i0, i3, ishift_time,ccm)

         tshift = -ishift_time*ddt
         call calc_adj_tr(tshift,tr_synth, tr_adj_src, tr_data, npts, ddt)
         adj_src_TT(:,icomp,irec_local)=tr_adj_src(:)



         adjust_time = adjust_time + 0.5 * tshift**2
         trav_time_ccZ(1,ista,isrc)=tshift
         trav_time_ccZ(2,ista,isrc)=ccm

     endif


     endif
     enddo

     if (icomp_to_invert(isrc) > 4) then
       adjust = adjust_time * scale_cost
       call mpi_reduce(trav_time_ccZ,trav_ttZ,2*nb_DSM_sources*nrec,CUSTOM_MPI_TYPE ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
       call mpi_reduce(trav_time_ccX,trav_ttX,2*nb_DSM_sources*nrec,CUSTOM_MPI_TYPE ,MPI_SUM,0, MPI_COMM_WORLD, ierr)

       if (myrank==0) then
          open(666,file='trav_time.out',access='append')
          write(666,*) '----- '
          do i=1,nrec
             if (icomp_to_invert(isrc) == 5) write(666,*) i,isrc,trav_ttZ(1,i,isrc),trav_ttZ(2,i,isrc)
             if (icomp_to_invert(isrc) == 6) write(666,'(2i10,4f20.10)') i,isrc,trav_ttZ(1,i,isrc),trav_ttZ(2,i,isrc),&
              trav_ttX(1,i,isrc),trav_ttX(2,i,isrc)
          enddo
          close(666)
       endif
     else
       adjust = 0.5 * adjust * deltat * scale_cost
     endif

     call mpi_reduce(adjust,data_adjustment,1,CUSTOM_MPI_TYPE ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(data_adjustment,1,CUSTOM_MPI_TYPE,0, MPI_COMM_WORLD, ierr)
     fitting_value = fitting_value + data_adjustment
     !write(*,*) 'fitting value : ',myrank,isrc,fittiing_value,adjust,scale_cost
!!$     write(*,*) 'size adj' , Nstep,
!!$   open(27,file='debug_adj_src',access='direct',recl=4*Nstep*3*nrec_local)
!!$   write(27,rec=1)  adj_src_all
!!$   close(27)

  end subroutine define_adjoint_sources

  subroutine calc_adj_tr(tshift,tr_synth,adj_src,vl,npts,dt)
    implicit none
    integer npts,i
    real(kind=CUSTOM_REAL) tshift,tr_synth(npts),adj_src(npts),vl(npts)
    real(kind=CUSTOM_REAL) dt,norm

    do i=2,npts-1
      vl(i) =  (tr_synth(i+1) - tr_synth(i-1) ) / (2*dt)
    enddo
    vl(1)=(tr_synth(2) - tr_synth(1) ) / (dt)
    vl(npts)=(tr_synth(npts) - tr_synth(npts-1) ) / (dt)

    norm = dt * sum(vl(:)*vl(:))
    adj_src(:) =  tshift * vl(:) / norm
    !write(*,*) "tshift ",tshift
  end subroutine calc_adj_tr


  subroutine calc_max_corr(s,d,n,i_inf,i_sup,ishift,cc)
    implicit none
    integer n,i_inf,i_sup,ishift
    integer j,ileft,iright,ilag,k
    real(kind=CUSTOM_REAL) s(n),d(n)
    real(kind=CUSTOM_REAL) norm0,norm1,norm,ccmax,cc

    ishift=0;k=0;

    ileft=-int((i_sup-i_inf)/2)
    iright=int((i_sup-i_inf)/2)

    norm0=sum(s(i_inf:i_sup)**2)
    norm1=sum(d(i_inf:i_sup)**2)
    norm=sqrt(norm1)*sqrt(norm0)

    ccmax=-1e30;

    do ilag=ileft,0
       k=k+1
       cc=0.
       do j=1-ilag,n
          cc=cc+s(j)*d(j+ilag)
       enddo
       if (cc > ccmax) then
          ccmax=cc
          ishift=ilag
       endif
    enddo

    do ilag=1,iright
       k=k+1
       cc=0.
       do j=1,n-ilag
          cc=cc+s(j)*d(j+ilag)
       enddo
       if (cc > ccmax) then
          ccmax=cc
          ishift=ilag
       endif
    enddo
    !if (myrank==0) then
     !open(667,file='debug_tt',access='append')
     write(*,*) myrank,ccmax,norm0,norm1,norm,ccmax/norm
     !close(667)
    !endif
    cc=ccmax/norm

  end subroutine calc_max_corr



  subroutine tapper_window_W(t_w,i0,i1,i2,i3,nstep,W)
    implicit none
    integer i,i0,i1,i2,i3,nstep
    real(kind=CUSTOM_REAL) t_w(nstep),omega,phi,pi,W

    PI = 3.1415926d0
    t_w(:)=0.


    omega = pi / (i1 - i0)
    phi = pi / 2 - omega * i1
   !write(*,*) i0,i1,W,omega,phi
   do i = i0, i1
       t_w(i) = W*(0.5 + 0.5 *sin(omega * i + phi))
   enddo
   do i = i1+1,i2-1
     t_w(i)=W
   enddo
   omega = pi / (i3 - i2)
   phi = pi/2 - omega * i2
   do i= i2,i3
     t_w(i) = W*(0.5 + 0.5 * sin(omega * i + phi))
   enddo
end subroutine




!#############################################################
subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

  ! recursive filtering of data with butterworth filter
  ! x: input array
  ! y: output array
  ! dt: time increment
  ! n: number of data points

  ! irek=0: forward filtering only
  ! irek=1: forward and backward filtering

  ! norder: order of butterworth filter
  ! norder=0: only filtering, no determination of coefficients
  ! norder<0: no starplots of transfer function and impulse response

  ! f1: low cutoff frequency (Hz)
  ! f1=0: low pass filter

  ! f2: high cutoff frequency (Hz)
  ! f2>0.5/dt: high pass filter

  implicit none
  integer n
  real(kind=CUSTOM_REAL) x(n),y(n)
  real(kind=CUSTOM_REAL), dimension (10) ::  a, b1, b2
  real(kind=CUSTOM_REAL) :: dt,f1,f2
  integer :: iunit, npoles,norder,irek,lx
  !real(kind(0d0)) :: x(n),y(n)

   iunit = 3

   if(norder/=0) then
      npoles=iabs(norder)
      !determination of filter coefficients
      call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
      if(norder>=0) then
         !plot of transfer function and impuulse response
         lx = 100
         !filtering
      endif
   endif


   if(n/=0) then
      call rekurs(x,y,n,a,b1,b2,npoles,irek)
   endif
   return
 end subroutine bwfilt



!---------------------------------------------------------------


subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag/=0: forward and backward filtering

  implicit none

  real(kind=CUSTOM_REAL), dimension(10) :: z,z1,z2 ,a,b1,b2
  real(kind=CUSTOM_REAL) ::  x1,x2
  integer :: ndat, npoles, iflag, n,i
  real(kind=CUSTOM_REAL) :: x(ndat), y(ndat)

  !forward

  x1 = 0.d0
  x2 = 0.d0

  do i = 1, npoles
     z1(i) = 0.d0
     z2(i) = 0.d0
  enddo

  do n = 1, ndat
     z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
     do i = 2, npoles
        z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
     enddo
     x2=x1
     x1=x(n)
     do i = 1, npoles
        z2(i) =z1(i)
        z1(i) =z(i)
     enddo
     y(n) = z(npoles)
  enddo

  if(iflag==0) then
     return
  endif

  !backward

  x1 =0.d0
  x2 =0.d0

  do i = 1, npoles
     z1(i) = 0.d0
     z2(i) = 0.d0
  enddo

  do n = ndat, 1, -1
     z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
     do i =2, npoles
        z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
     enddo
     x2=x1
     x1=y(n)
     do i = 1,npoles
        z2(i)=z1(i)
        z1(i)=z(i)
     enddo
     y(n) = z(npoles)
  enddo
  return
end subroutine rekurs






!---------------------------------------------------------------


subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter

  real(kind=CUSTOM_REAL),dimension(10) :: a,b1,b2
  complex(kind=CUSTOM_REAL) :: s(20), t1,t2,p
  real(kind=CUSTOM_REAL), parameter :: pi = 3.141592653589793d0
  real(kind=CUSTOM_REAL) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles


  if(npoles>10) then
     stop ' npoles greater than 10: STOP '
  endif

  d2= 2.d0/dt
  w1=d2*tan(2.d0*pi*f1/d2)
  w2=d2*tan(2.d0*pi*f2/d2)
  w0=0.5*(w2-w1)

  i=1
  npol2=npoles/2+1
  do n =1,npoles
     p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
     t1 = p*cmplx(w0,0.d0)
     t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
     s(i)=t1+t2
     s(i+1)=t1-t2
     i=i+2
  enddo

  do n=1,npoles
     ssum=2*real(s(n))
     sprod=dble(s(n)*conjg(s(n)))
     fact1=d2*d2-d2*ssum+sprod
     fact2=2.d0*(sprod-d2*d2)
     fact3=d2*d2+d2*ssum+sprod
     a(n)=2.d0*d2*w0/fact1
     b1(n)=fact2/fact1
     b2(n)=fact3/fact1
  enddo
  return
end subroutine bpcoeff


!=====================================================================

  subroutine write_channel_name1(iorientation,channel)

  use specfem_par,only: DT,SUPPRESS_UTM_PROJECTION
  implicit none

  integer :: iorientation
  character(len=3) :: channel

  ! local parameters
  character(len=2) :: bic
  double precision:: sampling_rate

  ! gets band and instrument code
  sampling_rate = DT
  call band_instrument_code1(sampling_rate,bic)

  ! sets channel name
  if( SUPPRESS_UTM_PROJECTION ) then

    ! no UTM, pure Cartesian reference
    ! uses Cartesian X/Y/Z direction to denote channel
    select case(iorientation)
    case(1)
      channel = bic(1:2)//'X'
    case(2)
      channel = bic(1:2)//'Y'
    case(3)
      channel = bic(1:2)//'Z'
    case default
      call exit_mpi(0,'error channel orientation value')
    end select

  else

    ! UTM conversion
    ! uses convention for N/E/Z to denote channel
    select case(iorientation)
    case(1)
      channel = bic(1:2)//'E'
    case(2)
      channel = bic(1:2)//'N'
    case(3)
      channel = bic(1:2)//'Z'
    case default
      call exit_mpi(0,'error channel orientation value')
    end select

  endif

  end subroutine write_channel_name1

!=====================================================================



!=====================================================================

  subroutine band_instrument_code1(DT,bic)
  ! This subroutine is to choose the appropriate band and instrument codes for
  ! channel names of seismograms
  ! based on the IRIS convention (first two letters of channel codes,
  ! respectively,
  ! which were LH(Z/E/N) previously).
  ! For consistency with observed data, we now use the IRIS convention for band
  ! codes (first letter in channel codes) of
  ! SEM seismograms governed by their sampling rate.
  ! Instrument code (second letter in channel codes) is fixed to "X" which is
  ! assigned by IRIS for synthetic seismograms.
  ! See the manual for further explanations!
  ! Ebru, November 2010
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
  if( SUPPRESS_IRIS_CONVENTION ) then
    bic = 'BH'
  endif

 end subroutine band_instrument_code1


end module adjoint_source
