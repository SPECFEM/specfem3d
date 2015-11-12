  module reading_field

    contains

      subroutine read_mesh(isim)
        use global_parameters
        character(len=250) file_to_read
        character(len=4)  appmynum
        real(kind=CUSTOM_REAL), allocatable :: ssol(:,:,:),zsol(:,:,:)
        integer isim,iproc,iel,indx_stored


        !isim=1 !! read mesh just for the first simulation
        ! number of elements for each procs
        allocate(nb_stored(0: nbproc-1))

        nel=0

        ! count total elements
        do iproc=0, nbproc-1

           write(file_to_read,'(a10,i5.5,a4)')'parameters',iproc,'.par'
           open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'/'//trim(file_to_read))
           read(10,*) nb_stored(iproc),ibeg,iend
           close(10)
           nel = nel + nb_stored(iproc)
           write(*,*) iproc,nb_stored(iproc),ibeg,iend
           write(*,*) nel
           write(*,*)

        enddo
        allocate(scoor(ibeg:iend,ibeg:iend,nel),zcoor(ibeg:iend,ibeg:iend,nel))
        allocate(data_read(ibeg:iend,ibeg:iend,nel))
        allocate(depth_ele(nel))
        ! reading mesh  !! LE PATH DE MESH DEPEND DU TYPE DE SOURCE ie adapter
        ! working_axisem_dir:
        indx_stored=1
        do iproc=0,nbproc-1

           if (nb_stored(iproc) > 0) then
              call define_io_appendix(appmynum,iproc)
              allocate(ssol(ibeg:iend,ibeg:iend, nb_stored(iproc)),zsol(ibeg:iend,ibeg:iend, nb_stored(iproc)))
              file_to_read=trim(working_axisem_dir)//trim(simdir(isim))//'Data/mesh_sol_'//appmynum//'.dat'
              open(10,file=trim(file_to_read),FORM="UNFORMATTED")
              read(10) ssol(ibeg:iend,ibeg:iend,:),zsol(ibeg:iend,ibeg:iend,:)
!!$              write(*,*) '----------'
!!$              write(*,*) ssol(ibeg:iend,ibeg:iend,:)
!!$              write(*,*) '----------'
              close(10)
              scoor(ibeg:iend,ibeg:iend,indx_stored:indx_stored+nb_stored(iproc)-1)= ssol(ibeg:iend,ibeg:iend,:)
              zcoor(ibeg:iend,ibeg:iend,indx_stored:indx_stored+nb_stored(iproc)-1)= zsol(ibeg:iend,ibeg:iend,:)
              indx_stored=indx_stored+nb_stored(iproc)
              deallocate(ssol,zsol)

           endif
        enddo

        do iel=1,nel
           depth_ele(iel) = sqrt(scoor(2,2,iel)**2+zcoor(2,2,iel)**2)
        enddo



      end subroutine read_mesh

      subroutine read_veloc_field_and_interpol(isim)
        use mpi_mod
        use global_parameters
        use post_processing
        use writing_mod
        real(kind=SINGLE_REAL), allocatable :: data_to_read(:,:,:),stack_v(:),stalta(:)
        integer itime,iproc,indx_stored(3),isim,ifield,i
        integer ivx,ivy,ivz
        integer, allocatable :: iunit(:,:)
        character(len=4) appmynum
        character(len=256) fichier
        integer nlta,nsta
        real(kind=SINGLE_REAL) Energy_1,Energy_0,thres

        nlta=100
        nsta=20
        thres=0.1

        !real(kind=SINGLE_REAL) f1,f2,phi
        !isim=1 !! hardcoded

        allocate(iunit(0:nbproc-1,3))


        ! unit file
        i=150
        do ifield=1,3
           do iproc=0, nbproc-1
              i=i+1
              iunit(iproc,ifield)=i
           enddo
        enddo

        i=i+1
        ivx=i
        i=i+1
        ivy=i
        i=i+1
        ivz=i



        !write(*,*)  iunit
        write(*,*) ' src_type(:,1) ', myrank, src_type(:,1)
        write(*,*) ' src_type(:,2) ', myrank, src_type(:,2)
        call compute_prefactor(src_type(isim,1),src_type(isim,2),isim)

        if (myrank == 0) then
          write(fichier,'(a6,a15)') '/Data/',output_veloc_name(1)
          open(ivx,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_veloc_name(2)
          open(ivy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_veloc_name(3)
          open(ivz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")

          write(*,*) 'nbrec to write', nbrec
          write(*,*) 'nbrec to read ',sum(nb_stored(:))
          write(*,*) 'nt to write ', ntime


          write(ivx) nbrec,ntime
          write(ivy) nbrec,ntime
          write(ivz) nbrec,ntime


         data_rec=0.

         do ifield=1,3
           write(*,*) ifield
           if (trim(src_type(isim,1))=='monopole' .and. ifield==2) then
              write(*,*) 'monopole => not up'
              cycle
           endif

           ! open files
           do iproc=0, nbproc-1
               call define_io_appendix(appmynum,iproc)
               write(fichier,'(a6,a15,a1)') '/Data/',input_veloc_name(ifield),'_'
               open(unit=iunit(iproc,ifield), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                   //appmynum//'.bindat', &
                   FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
           enddo
         enddo
        endif

        allocate(stack_v(ntime),stalta(ntime))

         !data_rec=0.
         Energy_1=0.
         Energy_0=0.
         ! read files
         do itime=1,ntime
            !write(*,*) ' reading field :', itime
            data_rec=0.
            Energy_1=0.

            indx_stored=1
            !! us
            ifield=1
            if (myrank ==0) then
               do iproc=0, nbproc-1
                  if (nb_stored(iproc) > 0) then

                     allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                     read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                     data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                          = data_to_read(ibeg:iend,ibeg:iend,:)
                     indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                     Energy_1 = Energy_1 + sum(data_to_read(ibeg:iend,ibeg:iend,:)**2)
                     deallocate(data_to_read)

                  endif
               enddo
            endif
            call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            call interpol_field(ifield)

            !write(*,*)  data_read(1,1,1),data_rec(1,1)
            if (.not.(trim(src_type(isim,1))=='monopole')) then
               ifield=2
               if (myrank == 0) then
                  do iproc=0, nbproc-1
                     if (nb_stored(iproc) > 0) then

                        allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                        read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                        data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                             = data_to_read(ibeg:iend,ibeg:iend,:)
                        indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                        Energy_1 = Energy_1 + sum(data_to_read(ibeg:iend,ibeg:iend,:)**2)
                        deallocate(data_to_read)

                     endif
                  enddo
               endif
               call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
               call interpol_field(ifield)
            endif
            !! uz
            ifield=3
            if (myrank ==0) then
               do iproc=0, nbproc-1
                  if (nb_stored(iproc) > 0) then

                     allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                     read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                     data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                          = data_to_read(ibeg:iend,ibeg:iend,:)
                     indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                     Energy_1 = Energy_1 + sum(data_to_read(ibeg:iend,ibeg:iend,:)**2)
                     deallocate(data_to_read)

                  endif
               enddo
               !if (itime==1) then
               !   Energy_0=Energy_1
               !else
               !   write(*,*) itime,Energy_1,Energy_0,dtt*pick_tresh
               !   if ((Energy_1-Energy_0) > dtt*pick_tresh) then
               !      ipick=itime
               !   else
               !      Energy_0=Energy_1
               !   endif
               !endif
               stack_v(itime)=Energy_1
            endif
            call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            call interpol_field(ifield)



           call  compute_3D_cyl()
           !if (irecmin < 10 .and. irecmax > 10) write(*,*) '1',data_rec(10,2)

           call rotate2cartesian_with_source_in_pole()

           !if (irecmin < 10 .and. irecmax > 10) write(*,*) '2',data_rec(10,2)

           call rotate_back_source() ! cartesien dans le repere terrestre global

           !if (irecmin < 10 .and. irecmax > 10) write(*,*) '3',data_rec(10,2)

           call rotate_back_to_local_cart() ! cartesien local

           call rotate_from_chunk_azimuth() !! VM VM add rotation azimuth chunk
           !if (irecmin < 10 .and. irecmax > 10) write(*,*) '4',data_rec(10,2)

           call reduce_mpi_veloc()
           !if (myrank == 0) write(*,*) '5',data_rec(10,2)
           if (myrank == 0) call write_veloc_or_displ_3D(ivx,ivy,ivz)



         enddo ! pas de temps

         if (myrank ==0) then
         ! close files
         do ifield=1,3

            do iproc=0, nbproc-1
               close(iunit(iproc,ifield))
            enddo
         enddo
         close(ivx)
         close(ivy)
         close(ivz)

         call substalta(stack_v, nsta, nlta, stalta,ntime)
         call pick(stalta,ntime,ipick,thres)
         open(ivx,file='pick_guess.txt')
         write(ivx,*) ipick,ipick*dtt
         close(ivx)
         open(ivx,file='StaLta.txt')
         do ipick=1,ntime
             write(ivx,*) ipick*dtt,stack_v(ipick),stalta(ipick)
         enddo
         close(ivx)
         deallocate(stack_v,stalta)

        endif




      end subroutine read_veloc_field_and_interpol

      subroutine pick(stalta,n,i,thres)
        use global_parameters
        integer i,n
        real(kind=SINGLE_REAL) stalta(n),thres

        do i=1,n
          if (stalta(i) >= thres) exit
        enddo

      end subroutine pick

      subroutine read_stress_field_and_interpol(isim)
        use mpi_mod
        use global_parameters
        use post_processing
        use writing_mod
        real(kind=SINGLE_REAL), allocatable ::  data_to_read(:,:,:)
        integer itime,iproc,indx_stored(6),isim,ifield,i
        integer isxx,isyy,iszz,isxy,isxz,isyz
        integer, allocatable :: iunit(:,:)
        character(len=4) appmynum
        character(len=256) fichier
        !real(kind=SINGLE_REAL) f1,f2,phi
        !isim=1 !! hardcoded

        allocate(iunit(0:nbproc-1,6))


        ! unit file
        i=150
        do ifield=1,6
           do iproc=0, nbproc-1
              i=i+1
              iunit(iproc,ifield)=i
           enddo
        enddo

        i=i+1
        isxx=i
        i=i+1
        isyy=i
        i=i+1
        iszz=i
        i=i+1
        isxy=i
        i=i+1
        isxz=i
        i=i+1
        isyz=i

        !write(*,*)  iunit
        !write(*,*) src_type
        call compute_prefactor(src_type(isim,1),src_type(isim,2),isim)
        if (myrank == 0) then
          write(fichier,'(a6,a15)') '/Data/',output_stress_name(1)
          open(isxx,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_stress_name(2)
          open(isyy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_stress_name(3)
          open(iszz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_stress_name(4)
          open(isxy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_stress_name(5)
          open(isxz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
          write(fichier,'(a6,a15)') '/Data/',output_stress_name(6)
          open(isyz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")

          write(isxx) nbrec,ntime
          write(isyy) nbrec,ntime
          write(iszz) nbrec,ntime
          write(isxy) nbrec,ntime
          write(isxz) nbrec,ntime
          write(isyz) nbrec,ntime

          data_rec=0.
          do ifield=1,6
             write(*,*) ifield
             if (trim(src_type(isim,1))=='monopole' .and. (ifield==4 .or. ifield==6)) then
                write(*,*) 'monopole => not up'
                cycle
             endif
             ! open files
             do iproc=0, nbproc-1
                call define_io_appendix(appmynum,iproc)
                write(fichier,'(a6,a15,a1)') '/Data/',input_stress_name(ifield),'_'
                open(unit=iunit(iproc,ifield), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                   //appmynum//'.bindat', &
                   FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
             enddo
          enddo
        endif

        ! read files
        do itime=1,ntime
           stress_rec=0.
           indx_stored=1
           !! s11
           ifield=1
           if (myrank == 0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then

                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)

                  ! write(*,*) data_to_read(ibeg:iend,ibeg:iend,:)

                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                      = data_to_read(ibeg:iend,ibeg:iend,:)

                   !write(*,*)  data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1)

                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)

                endif
             enddo
           endif
           call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
           call interpol_stress(ifield)


           !! s22
           ifield=2
           if (myrank==0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then

                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)

                endif
             enddo
           endif
           call  mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
           call interpol_stress(ifield)


           !! s33
           ifield=3
           if (myrank==0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then

                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                      = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)

                endif

             enddo
           endif
           call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
           call interpol_stress(ifield)

           !! s12
           if (.not.(trim(src_type(isim,1))=='monopole')) then
              ifield=4
              if (myrank==0) then
                do iproc=0, nbproc-1
                   if (nb_stored(iproc) > 0) then

                    allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                    read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                    data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                         = data_to_read(ibeg:iend,ibeg:iend,:)
                    indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                    deallocate(data_to_read)

                  endif
                enddo
              endif
              call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
              call interpol_stress(ifield)
           endif


           !! s13
           ifield=5
           if (myrank == 0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then

                 allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                 read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                 data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                      = data_to_read(ibeg:iend,ibeg:iend,:)
                 indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                 deallocate(data_to_read)

                endif

             enddo
           endif
           call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
           call interpol_stress(ifield)


           !! s23
           if (.not.(trim(src_type(isim,1))=='monopole')) then
              ifield=6
              if (myrank ==0) then
                do iproc=0, nbproc-1
                   if (nb_stored(iproc) > 0) then

                    allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                    read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                    data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                         = data_to_read(ibeg:iend,ibeg:iend,:)
                    indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                    deallocate(data_to_read)

                   endif
                enddo
              endif
              call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
              call interpol_stress(ifield)
           endif


           call  compute_stress_3D_cyl()

           call rotate2cartesian_with_source_in_pole_stress()

           call rotate_back_source_stress() ! cartesien dans le repere tereste global

           call rotate_back_to_local_cart_stress() ! cartesien local

           call rotate_from_chunk_azimuth_stress() !! VM VM add rotation azimuth chunk

           call reduce_mpi_stress()

           !write(*,*) 'appel a write_stress3D'
           if (myrank == 0) call  write_stress3D(isxx,isyy,iszz,isxy,isxz,isyz)

     enddo

     ! close files
     do ifield=1,6
        do iproc=0, nbproc-1
           close(iunit(iproc,ifield))
        enddo
     enddo

     close(isxx)
     close(isyy)
     close(iszz)
     close(isxy)
     close(isxz)
     close(isyz)

   end subroutine read_stress_field_and_interpol

!! CD CD ========================================================================
!! CD CD add this ===============================================================

  subroutine read_displ_field_and_interpol(isim)

  use mpi_mod
  use global_parameters
  use post_processing
  use writing_mod

  real(kind=SINGLE_REAL), allocatable :: data_to_read(:,:,:),stack_v(:),stalta(:)
  real(kind=SINGLE_REAL) Energy_1,Energy_0,thres

  integer itime,iproc,indx_stored(3),isim,ifield,i
  integer iux,iuy,iuz
  integer nlta,nsta
  integer, allocatable :: iunit(:,:)

  character(len = 4) appmynum
  character(len = 256) fichier

!
!---
!

  nlta  = 100
  nsta  = 20
  thres = 0.1

  allocate(iunit(0:nbproc-1,3))


  ! unit file
  i = 150
  do ifield=1,3
    do iproc=0, nbproc-1
      i = i+1
      iunit(iproc,ifield) = i
    enddo
  enddo

  i   = i+1
  iux = i
  i   = i+1
  iuy = i
  i   = i+1
  iuz = i

  write(*,*) ' src_type(:,1) ', myrank, src_type(:,1)
  write(*,*) ' src_type(:,2) ', myrank, src_type(:,2)

  call compute_prefactor(src_type(isim,1),src_type(isim,2),isim)

  if (myrank == 0) then
    write(fichier,'(a6,a15)') '/Data/',output_displ_name(1)
    open(iux,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")

    write(fichier,'(a6,a15)') '/Data/',output_displ_name(2)
    open(iuy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")

    write(fichier,'(a6,a15)') '/Data/',output_displ_name(3)
    open(iuz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")

    write(*,*) 'nbrec to write', nbrec
    write(*,*) 'nbrec to read ', sum(nb_stored(:))
    write(*,*) 'nt to write ', ntime

    write(iux) nbrec,ntime
    write(iuy) nbrec,ntime
    write(iuz) nbrec,ntime

    data_rec = 0.

    do ifield=1,3
      write(*,*) ifield
      if (trim(src_type(isim,1))=='monopole' .and. ifield==2) then
        write(*,*) 'monopole => not up'
        cycle
      endif

      ! open files
      do iproc=0, nbproc-1

        call define_io_appendix(appmynum,iproc)
        write(fichier,'(a6,a15,a1)') '/Data/',input_displ_name(ifield),'_'
        open(unit=iunit(iproc,ifield), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                   //appmynum//'.bindat', FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
      enddo
    enddo
  endif

  allocate(stack_v(ntime),stalta(ntime))

  Energy_1 = 0.
  Energy_0 = 0.

  ! read files
  do itime=1,ntime

    data_rec    = 0.
    Energy_1    = 0.

    indx_stored = 1

    ifield      = 1 !! ====> us

    if (myrank == 0) then

      do iproc=0, nbproc-1
        if (nb_stored(iproc) > 0) then

          allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
          read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
          data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                  = data_to_read(ibeg:iend,ibeg:iend,:)
          indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
          Energy_1 = Energy_1 + sum(data_to_read(ibeg:iend,ibeg:iend,:)**2)
          deallocate(data_to_read)

        endif
      enddo

    endif

    call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call interpol_field(ifield)

    if (.not.(trim(src_type(isim,1)) == 'monopole')) then

      ifield = 2 !! ====> up

      if (myrank == 0) then
        do iproc = 0,nbproc-1
          if (nb_stored(iproc) > 0) then

            allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
            read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
            data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                    = data_to_read(ibeg:iend,ibeg:iend,:)
            indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
            Energy_1 = Energy_1 + sum(data_to_read(ibeg:iend,ibeg:iend,:)**2)
            deallocate(data_to_read)

          endif
        enddo

      endif

      call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call interpol_field(ifield)

    endif

    ifield = 3 !! ====> uz

    if (myrank == 0) then
      do iproc=0, nbproc-1
        if (nb_stored(iproc) > 0) then

          allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
          read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
          data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                  = data_to_read(ibeg:iend,ibeg:iend,:)
          indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
          Energy_1 = Energy_1 + sum(data_to_read(ibeg:iend,ibeg:iend,:)**2)
          deallocate(data_to_read)

        endif
      enddo

      stack_v(itime) = Energy_1

    endif

    call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call interpol_field(ifield)

    call  compute_3D_cyl()

    call rotate2cartesian_with_source_in_pole()
    call rotate_back_source() ! cartesien dans le repere terrestre global
    call rotate_back_to_local_cart() ! cartesien local
    call rotate_from_chunk_azimuth() !! VM VM add rotation azimuth chunk

    call reduce_mpi_veloc()

    if (myrank == 0) call write_veloc_or_displ_3D(iux,iuy,iuz)

  enddo ! pas de temps

  if (myrank ==0) then

    ! close files
    do ifield=1,3
      do iproc=0, nbproc-1
        close(iunit(iproc,ifield))
      enddo
    enddo
    close(iux)
    close(iuy)
    close(iuz)

!!         call substalta(stack_v, nsta, nlta, stalta,ntime)
!!         call pick(stalta,ntime,ipick,thres)
!!         open(iux,file='pick_guess.txt')
!!         write(iux,*) ipick,ipick*dtt
!!         close(iux)
!!         open(iux,file='StaLta.txt')
!!         do ipick=1,ntime
!!             write(iux,*) ipick*dtt,stack_v(ipick),stalta(ipick)
!!         enddo
!!         close(iux)
    deallocate(stack_v,stalta)

  endif


  end subroutine read_displ_field_and_interpol

!================================================================================
!================================================================================


      !-----------------------------------------------------------------------------
      subroutine define_io_appendix(app,iproc)
        !
        ! Defines the 4 digit character string appended to any
        ! data or io file related to process myid.
        !
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        implicit none
        integer, intent(in)           :: iproc
        character(len=4), intent(out) :: app

        write(app,"(I4.4)") iproc

      end subroutine define_io_appendix
      !=============================================================================
!================================================================================
!! Compute STA/LTA for picking
  subroutine substalta(sig, nsta, nlta, stalta,m)

    use global_parameters
    integer m
    real(kind=SINGLE_REAL) sig(m)
    real(kind=SINGLE_REAL) stalta(m)

    integer, intent(in) :: nsta, nlta

    integer :: nsta_1, nlta_1, i,j

    real(kind=SINGLE_REAL), dimension(:), allocatable :: sta, lta, pad_sta, pad_lta, tmp1,tmp2

    !m = size(sig)

    nsta_1 = nsta - 1
    nlta_1 = nlta - 1
    !write(*,*) m,nsta_1,nlta_1
    allocate(sta(m))
    allocate(lta(m))
    allocate(tmp1(m))
    allocate(tmp2(m))
    allocate(pad_sta(nsta_1))
    allocate(pad_lta(nlta_1))
    sta = 0.
    lta = 0.
    pad_sta = 0.
    pad_lta = 1.

    !*** compute the short time average (STA)
    !do i=1,nsta
    !   tmp1(1:nsta_1) = pad_sta(:)
    !   tmp1(nsta_1+1:m) = sig(i:m - nsta_1 + i-1)**2
    !   sta = sta + tmp1
    !enddo


    do i=nsta,m
       do j=i-nsta_1,i
         sta(i) = sta(i) + sig(j)
       enddo
    enddo
    sta = sta / nsta

    !*** compute the long time average (LTA)
    !do i =1,nlta
    !   tmp2(1:nlta_1) = pad_lta(:)
    !   tmp2(nlta_1+1:m) = sig(i:m - nlta_1 + i-1)**2
    !   lta = lta + tmp2
    !enddo
    do i=nlta,m
      do j=i-nlta_1,i
          lta(i)=lta(i)+sig(j)
      enddo
    enddo
    lta = lta / nlta

    sta(1:nsta_1) = 0.
    lta(1:nlta_1) = 1.

    do i=1,m
       if (lta(i) < 1e-10) then
          lta(i) = 1.
          sta(i) = 0.
       endif
    enddo
    stalta = sta / lta

  end subroutine substalta

  end module reading_field
