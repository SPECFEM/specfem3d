!###################################################################################################################################
!                                                                                                                                  !
!                                                                                                                                  !
!    MODULE input_output,  PURPOSE : READ ACQUISITON FOR FWI                                                                       !
!                                                                                                                                  !
!                                                                                                                                  !
!                   -- read all sources   (tractions, moment, forces, ...)                                                         !
!                   -- read all stations  (list of stations and position for each event )                                        !
!                   -- read all data      (components to be inverted)                                                              !
!                                                                                                                                  !
!    all information is stored in type : acqui                                                                                     !
!                                                                                                                                  !
!                                                                                                                                  !
!                                                                                                                                  !
! NOTE : MPI is based on WORLD define in specfem : my_local_mpi_comm_world and my_local_mpi_comm_for_bcast                         !
!                                                                                                                                  !
!  August 2016                                                                                                                     !
!  Vadim Monteiller, CNRS Marseille, France                                                                                        !
!                                                                                                                                  !
!##################################################################################################################################!

module input_output

  !! IMPORT VARIABLES FROM SPECFEM ---------------------------------------------------------
  use constants, only: mygroup

  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS, BROADCAST_SAME_MESH_AND_MODEL, ANISOTROPY

  use specfem_par, only: CUSTOM_REAL, HUGEVAL, NGNOD, NUM_ITER, NPROC, MAX_STRING_LEN, &
                                  NGLLX, NGLLY, NGLLZ, NDIM, NSPEC_AB, NGLOB_AB, MIDX, MIDY, MIDZ, &
                                  LOCAL_PATH, xigll, yigll, zigll, &
                                  ibool, xstore, ystore, zstore, &
                                  xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
                                  myrank, USE_SOURCES_RECEIVERS_Z,INVERSE_FWI_FULL_PROBLEM



  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  !!----------------------------------------------------------------------------------------

  use inverse_problem_par
  use mesh_tools
  use IO_model

  implicit none

  PUBLIC  :: init_input_output_mod, SetUpInversion, get_mode_running, dump_adjoint_sources, write_bin_sismo_on_disk, &
             WirteOutputs

  PRIVATE :: read_acqui_file, read_inver_file, store_default_acqui_values, is_blank_line, bcast_all_acqui, &
             get_stations, read_data_gather, create_name_database_inversion, read_and_distribute_events_for_simultaneous_runs

  !! DEFINITION OF PRIVATE VARIABLES
  integer,                PRIVATE       :: NEVENT
  real(kind=CUSTOM_REAL), PRIVATE       :: fl, fh
contains

!################################## ROUTINES THAT CAN BE CALLED FROM MAIN  #########################################################

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------------
!> \brief read input data, sources, receiver and store all in acqui_simu and send to the concerned MPI slices
!-------------------------------------------------------------------------------------------------------------------
  subroutine init_input_output_mod(inversion_param, acqui_simu, myrank)


    integer,                                         intent(in)    ::  myrank
    type(acqui),  dimension(:), allocatable,         intent(inout) ::  acqui_simu
    type(inver),                                     intent(inout) ::  inversion_param

    ! locals
    integer                                                        :: ievent
    character(len=MAX_LEN_STRING)                                  :: name_file
    character(len=MAX_LEN_STRING)                                  :: acqui_file, inver_file
    real(kind=CUSTOM_REAL)                                         :: elemsize_min_glob,elemsize_max_glob
    real(kind=CUSTOM_REAL)                                         :: distance_min_glob,distance_max_glob

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '          *********************************************'
       write(INVERSE_LOG_FILE,*) '          ***       READING INPUT PARAMETERS        ***'
       write(INVERSE_LOG_FILE,*) '          *********************************************'
       write(INVERSE_LOG_FILE,*)
    endif


    acqui_file=inversion_param%input_acqui_file
    inver_file=inversion_param%input_inver_file

    if (DEBUG_MODE) then
       write(name_file, '( "debug_file_input_ouput_module_rank",i10.10) ') myrank
       open(IIDD,file=trim(prefix_to_path)//trim(name_file))
    endif

    !! read acqui file
    if (myrank == 0) then
       call read_acqui_file(acqui_file, acqui_simu, myrank)
       call read_inver_file(inver_file, acqui_simu, inversion_param, myrank)
       call get_point_source(acqui_simu)
       call get_stations(acqui_simu)
    endif

    if (myrank == 0) call flush_iunit(INVERSE_LOG_FILE)

    call bcast_all_acqui(acqui_simu,  inversion_param, myrank)
    call locate_source(acqui_simu, myrank)
  !  call locate_receiver(acqui_simu, myrank)

    if (myrank == 0) call flush_iunit(INVERSE_LOG_FILE)

    !! not need to read data for only forward simulation
    if (.not. inversion_param%only_forward) call read_data_gather(acqui_simu, myrank)


    !! create name for outputs
    do ievent=1,acqui_simu(1)%nevent_tot
       call create_name_database_inversion(acqui_simu(ievent)%prname_inversion, myrank, ievent, LOCAL_PATH)
    enddo

    if (myrank == 0) call flush_iunit(INVERSE_LOG_FILE)

    !!! We can read input model from external files instead of get model stored in specfem databases -----------------------------
    !!! otherwise the model read from mesher cubit or meshfem3D is kept.
    if (inversion_param%input_SEM_model .and. inversion_param%input_FD_model) then

       if (myrank == 0) then
          write(*,*) ' ERROR : You want to read input model both in SEM and FD grid '
          write(*,*) ' Please choose one of the solution  '
       endif

       stop
    else

       if (inversion_param%input_FD_model) then

          if (myrank == 0) then
             write(INVERSE_LOG_FILE,*)
             write(INVERSE_LOG_FILE,*) '          *********************************************'
             write(INVERSE_LOG_FILE,*) '          ***      READING EXTERNAL FD MODEL        ***'
             write(INVERSE_LOG_FILE,*) '          *********************************************'
             write(INVERSE_LOG_FILE,*)
          endif

          call ReadInputFDmodel(inversion_param)

       endif

       if (inversion_param%input_SEM_model) then

          if (myrank == 0) then
             write(INVERSE_LOG_FILE,*)
             write(INVERSE_LOG_FILE,*) '          *********************************************'
             write(INVERSE_LOG_FILE,*) '          ***      READING EXTERNAL SEM MODEL       ***'
             write(INVERSE_LOG_FILE,*) '          *********************************************'
             write(INVERSE_LOG_FILE,*)
          endif

          call ReadInputSEMmodel(inversion_param)

       endif

       
       if (inversion_param%input_SEM_prior) then

          if (myrank == 0) then
             write(INVERSE_LOG_FILE,*)
             write(INVERSE_LOG_FILE,*) '          *********************************************'
             write(INVERSE_LOG_FILE,*) '          ***   READING EXTERNAL SEM PRIOR MODEL    ***'
             write(INVERSE_LOG_FILE,*) '          *********************************************'
             write(INVERSE_LOG_FILE,*)
          endif

          call ReadInputSEMpriormodel(inversion_param)

       endif


    endif

    !! get dimension of mesh
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB, &
         ibool,xstore,ystore,zstore, &
         inversion_param%xmin,inversion_param%xmax, &
         inversion_param%ymin,inversion_param%ymax, &
         inversion_param%zmin,inversion_param%zmax, &
         elemsize_min_glob,elemsize_max_glob, &
         distance_min_glob,distance_max_glob)

    !!-------------------------------------------------------------------------------------------------------------------------

    !!! write what we read to debug : check what is stored in type acqui_simu
    if (DEBUG_MODE) then

       ! check what I read
       write(IIDD,*)
       write(IIDD,*) '######################################################'
       write(IIDD,*)
       write(IIDD,*) 'MYRANK :' ,myrank, ' READ : '
       do ievent=1,acqui_simu(1)%nevent_tot
          write(IIDD,*)
          write(IIDD,*) 'event ',ievent
          write(IIDD,'(a)') trim(acqui_simu(ievent)%event_name)
          write(IIDD,'(a)') trim(acqui_simu(ievent)%source_file)
          write(IIDD,'(a)') trim(acqui_simu(ievent)%station_file)
          write(IIDD,'(a)') trim(acqui_simu(ievent)%data_file_gather)
          write(IIDD,*) ' EVENT POSITION : '
          write(IIDD,*) acqui_simu(ievent)%Xs,acqui_simu(ievent)%Ys, acqui_simu(ievent)%Zs

          select case (trim(adjustl(acqui_simu(ievent)%source_type)))
          case ('moment')
             write(IIDD,*) ' MOMENT TENSOR : '
             write(IIDD,*) acqui_simu(ievent)%Mxx,acqui_simu(ievent)%Myy, acqui_simu(ievent)%Mzz
             write(IIDD,*) acqui_simu(ievent)%Mxy,acqui_simu(ievent)%Myz, acqui_simu(ievent)%Myz
          case ('force')
             write(IIDD,*) ' FORCE : '
             write(IIDD,*) acqui_simu(ievent)%Fx,acqui_simu(ievent)%Fy, acqui_simu(ievent)%Fz
          end select

          write(IIDD,*) 'total station     :', acqui_simu(ievent)%nsta_tot
          write(IIDD,*) 'stations in slice :',acqui_simu(ievent)%nsta_slice
          if (acqui_simu(ievent)%nsta_slice > 0 .and. .not. inversion_param%only_forward) then
             write(IIDD,*) 'Check data stored : ', acqui_simu(ievent)%nsta_slice, acqui_simu(ievent)%Nt_data
             write(name_file,'(a16,i8.8,a1,i8.8)') 'Check_read_data_',ievent,'_',myrank
             write(IIDD,'(a,a)') 'Data that was read are dumped in file for checking: ',trim(name_file)
             open(IINN,file=trim(name_file),access='direct', &
                  recl=CUSTOM_REAL*acqui_simu(ievent)%nsta_slice*acqui_simu(ievent)%Nt_data)
             write(IINN,rec=1)  acqui_simu(ievent)%data_traces(:,:,1)
             write(IINN,rec=2)  acqui_simu(ievent)%data_traces(:,:,2)
             write(IINN,rec=3)  acqui_simu(ievent)%data_traces(:,:,3)
             close(IINN)
          endif
       enddo

    endif

  end subroutine init_input_output_mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!> read input inversion parameter file
!--------------------------------------------------------------------------------------------------------------------
  subroutine SetUpInversion(inversion_param, myrank)
    type(inver),                         intent(inout)    :: inversion_param
    integer,                             intent(in)       :: myrank
    character(len=MAX_LEN_STRING)                         :: acqui_file_ref

    write(6,*)
    write(6,*) '      SETUP INVERSION : ', NUMBER_OF_SIMULTANEOUS_RUNS, 'rank : ', myrank, ' group : ', mygroup
    write(6,*)
    call flush_iunit(6)

    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
       write(prefix_to_path,"('run',i4.4,'/')") mygroup + 1
       inversion_param%input_acqui_file=trim(prefix_to_path)//'/DATA/inverse_problem/acqui_file.txt'
       acqui_file_ref='./DATA/inverse_problem/acqui_file.txt'
       if (myrank == 0) then
          write(6,*) ' DISTRIBUTION OF EVENTS '
          call flush_iunit(6)
          !! only one process must do I/O (myrank=0, mygroup=0)
          if (mygroup == 0) call read_and_distribute_events_for_simultaneous_runs(NUMBER_OF_SIMULTANEOUS_RUNS, acqui_file_ref)
       endif
    else
       prefix_to_path='./'
       inversion_param%input_acqui_file='./DATA/inverse_problem/acqui_file.txt'
    endif

    inversion_param%input_inver_file='./DATA/inverse_problem/inversion_fwi_par'
    if (myrank == 0) then
       open(INVERSE_LOG_FILE,file=trim(prefix_to_path)//'output_inverse.txt')
       open(OUTPUT_ITERATION_FILE,file=trim(prefix_to_path)//'output_iterations.txt')
       open(OUTPUT_FWI_LOG,file=trim(prefix_to_path)//'output_fwi_log.txt')
       write(OUTPUT_ITERATION_FILE,'(a5,"|",a15,"|",a15,"|",a12,"|",a12,"|",a15)') &
            "iter "," cost function "," gradient norm ", " relat cost "," relat grad", " nb line search"
    endif

    !! to be sure that all processes in the MPI_COMM_WORLD wait for the end of deining new inputs files
    call synchronize_all_world()

  end subroutine SetUpInversion
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! read acqui_file, count events then distribute over groups and rewrite suitable acqui_file for each group
!--------------------------------------------------------------------------------------------------------------------
  subroutine read_and_distribute_events_for_simultaneous_runs(NUMBER_OF_SIMULTANEOUS_RUNS, acqui_file_ref)
    character(len=MAX_LEN_STRING),       intent(in)       :: acqui_file_ref
    integer,                             intent(in)       :: NUMBER_OF_SIMULTANEOUS_RUNS
    integer                                               :: number_of_events_in_acqui_file_ref
    integer                                               :: nevent_per_group, nevent_remained
    integer                                               :: igroup, ievent_in_group, ievent, ievent_global
    integer,  dimension(:),     allocatable               :: nevent_in_group
    character(len=MAX_LEN_STRING)                         :: line, prefix_to_path_tmp


    write(6,*)
    write(6,*)  ' NUMBER OF SIMULTANEOUS RUN > 0 '
    write(6,*)
    call flush_iunit(6)
    number_of_events_in_acqui_file_ref=0
    open(666,file=trim(acqui_file_ref))
    do
       read(666,'(a)',end=99) line
       !! no significant line
       if (is_blank_line(line)) cycle
       !! new event
       if (INDEX(line,'event_name') > 0) number_of_events_in_acqui_file_ref=number_of_events_in_acqui_file_ref+1
    enddo
99  close(666)


    nevent_per_group = number_of_events_in_acqui_file_ref / NUMBER_OF_SIMULTANEOUS_RUNS
    nevent_remained =  mod( number_of_events_in_acqui_file_ref, NUMBER_OF_SIMULTANEOUS_RUNS)

    allocate(nevent_in_group(NUMBER_OF_SIMULTANEOUS_RUNS))
    do ievent=1,NUMBER_OF_SIMULTANEOUS_RUNS
       if (ievent <= nevent_remained) then
          nevent_in_group(ievent)= nevent_per_group+1
       else
          nevent_in_group(ievent)= nevent_per_group
       endif
    enddo

    ievent_global = 0
    ievent_in_group=0
    igroup=1
    open(666,file=trim(acqui_file_ref))
    write(prefix_to_path_tmp,"('run',i4.4,'/')") igroup
    open(777, file=trim(prefix_to_path_tmp)//'DATA/inverse_problem/acqui_file.txt')
    do
       read(666,'(a)',end=999) line
       !! no significant line
       if (is_blank_line(line)) cycle

       !! new event
       if (INDEX(line,'event_name') > 0) then
           ievent_global = ievent_global + 1
           write(6,*)
           write(6,*) '   next event ', ievent_global
           ievent_in_group=ievent_in_group+1

       endif

       !! write lines related to the current event
       if (ievent_in_group > nevent_in_group(igroup)) then
          igroup=igroup+1
          write(6,*) ' group ', igroup
          write(prefix_to_path_tmp,"('run',i4.4,'/')") igroup
          close(777)
          open(777, file=trim(prefix_to_path_tmp)//'DATA/inverse_problem/acqui_file.txt')
          ievent_in_group=1
       endif
       write(777, '(a)') trim(line)
       write(6,*) trim(line)
       call flush_iunit(6)
    enddo
999  close(666)
    close(777)
    deallocate(nevent_in_group)

  end subroutine read_and_distribute_events_for_simultaneous_runs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! create the name of the database for inverse problem
!--------------------------------------------------------------------------------------------------------------------
  subroutine create_name_database_inversion(prname,iproc,ievent,LOCAL_PATH)

    implicit none

    integer,                       intent(in)   :: iproc, ievent
    ! name of the database file
    character(len=MAX_LEN_STRING), intent(inout) :: prname
    character(len=MAX_STRING_LEN), intent(in   ) :: LOCAL_PATH
    ! local
    character(len=MAX_STRING_LEN)                :: procname

    ! create the name for the database of the current slide and region
    write(procname,"('/proc',i6.6,'_',i6.6,'_')") iproc,ievent
    ! create full name with path
    prname = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // procname

  end subroutine create_name_database_inversion

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! write final outputs model in disk
!--------------------------------------------------------------------------------------------------------------------
  subroutine WirteOutputs(inversion_param)

    type(inver),                                                intent(in)    :: inversion_param


    !! FOR NOT ONLY ONE OUPTPUT BUT FOR FURTHER DEV WE WILL ADD OTHER

    !! writing the final solution in SEM mesh
    call WriteOuptutSEMmodel(inversion_param)


  end subroutine WirteOutputs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! get which option we want to run FWI or just direct problem
!--------------------------------------------------------------------------------------------------------------------
  subroutine get_mode_running(mode_running, inversion_param)

    use my_mpi
    use specfem_par, only: myrank
    include "precision.h"

    type(inver),                                    intent(inout) :: inversion_param
    character(len=MAX_LEN_STRING),                  intent(inout) :: mode_running
    integer                                                       :: ier

    if (myrank == 0) then
       call get_command_argument(1, mode_running)
    endif

    call MPI_BCAST(mode_running,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

    select case(trim(adjustl(mode_running)))

    case('forward')
       inversion_param%only_forward=.true.

    case('l-bfgs')
       inversion_param%only_forward=.false.

    case default

    end select

  end subroutine get_mode_running

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! master write adjoint sources gather to check
!----------------------------------------------------------------

  subroutine dump_adjoint_sources(iter, acqui_simu, myrank)

    integer,                                     intent(in)    :: iter, myrank
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu

    integer                                                    :: ievent
    character(len=MAX_LEN_STRING)                              :: name_file_tmp, ch_to_add

    do ievent=1, acqui_simu(1)%nevent_tot
       write(ch_to_add,'(i4.4,a4)') iter,'_adj'
       name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))
       call  write_bin_sismo_on_disk(ievent, acqui_simu, acqui_simu(ievent)%adjoint_sources,  name_file_tmp, myrank)
    enddo

  end subroutine dump_adjoint_sources
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! master write synthetics gather to check
!----------------------------------------------------------------
  subroutine dump_seismograms(iter, array_to_write,  acqui_simu, myrank)

    integer,                                                intent(in)    :: iter, myrank
    real(kind=CUSTOM_REAL),  dimension(:,:,:), allocatable, intent(in)    :: array_to_write
    type(acqui),  dimension(:), allocatable,                intent(inout) :: acqui_simu

    integer                                                               :: ievent
    character(len=MAX_LEN_STRING)                                         :: name_file_tmp, ch_to_add

    do ievent=1, acqui_simu(1)%nevent_tot
       write(ch_to_add,'(i4.4,a4)') iter,'_dir'
       name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))
       call  write_bin_sismo_on_disk(ievent, acqui_simu, array_to_write, name_file_tmp, myrank)
    enddo

  end subroutine dump_seismograms
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! master write synthetics gather to check
!----------------------------------------------------------------
  subroutine dump_filtered_data(iter, array_to_write,  acqui_simu, myrank)

    integer,                                                intent(in)    :: iter, myrank
    real(kind=CUSTOM_REAL),  dimension(:,:,:), allocatable, intent(in)    :: array_to_write
    type(acqui),  dimension(:), allocatable,                intent(inout) :: acqui_simu

    integer                                                               :: ievent
    character(len=MAX_LEN_STRING)                                         :: name_file_tmp, ch_to_add

    do ievent=1, acqui_simu(1)%nevent_tot
       write(ch_to_add,'(i4.4,a4)') iter,'_fil'
       name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))
       call  write_bin_sismo_on_disk(ievent, acqui_simu, array_to_write, name_file_tmp, myrank)
    enddo

  end subroutine dump_filtered_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! master write waveform synthetic data gather (need to choose between write_bin_sismo_on_disk or write_gather_on_disk)
!----------------------------------------------------------------
  subroutine write_bin_sismo_on_disk(ievent, acqui_simu, array_to_write, name_file_to_write, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                                intent(in)    :: myrank, ievent
    character(len=MAX_LEN_STRING),                          intent(in)    :: name_file_to_write
    real(kind=CUSTOM_REAL),  dimension(:,:,:), allocatable, intent(in)    :: array_to_write
    type(acqui),             dimension(:),     allocatable, intent(inout) :: acqui_simu

    integer                                                               :: idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                               :: tag, ier, nsta_irank, irank, icomp
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable                 :: Gather, Gather_loc
    integer,                dimension(:),     allocatable                 :: irec_global
    integer                                                               :: status(MPI_STATUS_SIZE)

     if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) '  ... Writing  synthetic data gather for event :  ', ievent
     endif

     if (myrank == 0) then
        NSTA=acqui_simu(ievent)%nsta_tot
        Nt=acqui_simu(ievent)%Nt_data
        allocate(Gather(NSTA,Nt,NDIM))
     endif


     do irank = 1, NPROC-1

        if (myrank == 0) then

           ! count the receiver in slice irank
           nsta_irank=0
           do irec = 1,  NSTA
              if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
           enddo
           if (nsta_irank > 0) then
              allocate(Gather_loc(nsta_irank,Nt,NDIM))  !! data to receive
              allocate(irec_global(nsta_irank))
              irec_local=0
              tag   = MPI_ANY_TAG

              call MPI_RECV(Gather_loc, Nt*nsta_irank*NDIM, CUSTOM_MPI_TYPE, irank, tag, my_local_mpi_comm_world, status,  ier)
              call MPI_RECV(irec_global, nsta_irank, MPI_INTEGER, irank, tag, my_local_mpi_comm_world, status,  ier)

              do icomp=1,NDIM
                 do irec_local = 1, nsta_irank
                    Gather(irec_global(irec_local), :, icomp) = Gather_loc(irec_local, :, icomp)
                 enddo
              enddo


              deallocate(Gather_loc)
              deallocate(irec_global)
           endif
        else
           if (myrank == irank .and. acqui_simu(ievent)%nsta_slice > 0) then
              NSTA_LOC=acqui_simu(ievent)%nsta_slice
              Nt=acqui_simu(ievent)%Nt_data
              allocate(Gather_loc(NSTA_LOC,Nt,NDIM))
              allocate(irec_global(NSTA_LOC))
              do irec_local = 1, NSTA_LOC
                 irec_global(irec_local) = acqui_simu(ievent)%number_receiver_global(irec_local)

                 !! choose the rigth seismograms_*
                 do icomp=1,NDIM
                    select case (trim(acqui_simu(ievent)%component(icomp)))
                    case ('PR')
                       Gather_loc(irec_local,:,icomp)=array_to_write(1,irec_local,:)
                    case ('UX')
                       Gather_loc(irec_local,:,icomp)=array_to_write(1,irec_local,:)
                    case ('UY')
                       Gather_loc(irec_local,:,icomp)=array_to_write(2,irec_local,:)
                    case ('UZ')
                       Gather_loc(irec_local,:,icomp)=array_to_write(3,irec_local,:)
                    end select
                 enddo

              enddo

              tag    = 2001
              call MPI_SEND(Gather_loc,  Nt*NSTA_LOC*NDIM, CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, ier)
              call MPI_SEND(irec_global, NSTA_LOC, CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, ier)

              deallocate(irec_global)
           endif

        endif

        !! not sure if need this sync
        call synchronize_all()

     enddo


     !!  write gather file
     if (myrank == 0) then
        do icomp=1,NDIM
           do irec_local = 1, acqui_simu(ievent)%nsta_slice
              !! choose the rigth seismograms_*
              select case (trim(acqui_simu(ievent)%component(icomp)))
              case ('PR')
                 Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(1,irec_local,:)
              case ('UX')
                 Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(1,irec_local,:)
              case ('UY')
                 Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(2,irec_local,:)
              case ('UZ')
                 Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(3,irec_local,:)
              end select
           enddo
        enddo
        open(IINN,file=trim(adjustl(name_file_to_write)), access='direct', recl=CUSTOM_REAL*Nt*NSTA, status='replace')

        !! write only the asked component or pressure
        irec=0
        do idim=1,NDIM

           select case (trim(acqui_simu(ievent)%component(idim)))
           case('UX', 'UY', 'UZ', 'PR')
              irec=irec+1
              write(IINN,rec=irec) Gather(:,:,idim)
           end select

        enddo
        close(IINN)
        deallocate(Gather)
     endif

   end subroutine write_bin_sismo_on_disk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! master read waveform data gather and bcast to MPI slice concerned
!----------------------------------------------------------------
  subroutine read_data_gather(acqui_simu, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                     intent(in)    :: myrank
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu

    integer                                                    :: ievent, idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                    :: tag, ier, nsta_irank, irank
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable      :: Gather, Gather_loc
    integer                                                    :: status(MPI_STATUS_SIZE)

    if (myrank == 0) write(INVERSE_LOG_FILE,'(/a17)') '... reading data '

    do ievent = 1, acqui_simu(1)%nevent_tot


       if (myrank == 0) then

          NSTA=acqui_simu(ievent)%nsta_tot
          Nt=acqui_simu(ievent)%Nt_data

          allocate(Gather(NSTA,Nt,NDIM))
          Gather(:,:,:) = 0._CUSTOM_REAL
          ! read gather file
          open(IINN,file=trim(adjustl(acqui_simu(ievent)%data_file_gather)), access='direct', &
               recl=CUSTOM_REAL*Nt*NSTA,status='old')
          !! read only the asked component or pressure
          irec=0
          do idim=1,NDIM

             select case (trim(acqui_simu(ievent)%component(idim)))
             case('UX', 'UY', 'UZ', 'PR')
                irec=irec+1
                read(IINN,rec=irec) Gather(:,:,idim)

             end select
          enddo
          close(IINN)

          !! store data gather in my slice if needed
          NSTA_LOC=acqui_simu(ievent)%nsta_slice
          allocate(acqui_simu(ievent)%data_traces(NSTA_LOC,Nt,NDIM))
          allocate(acqui_simu(ievent)%adjoint_sources(NDIM, NSTA_LOC, Nt))
          allocate(acqui_simu(ievent)%weight_trace(NDIM, NSTA_LOC))
          acqui_simu(ievent)%weight_trace(:,:)=1._CUSTOM_REAL

          if (VERBOSE_MODE .or. DEBUG_MODE)  allocate(acqui_simu(ievent)%synt_traces(NDIM, NSTA_LOC, Nt))

          irec_local=0
          do irec = 1, NSTA
             if (acqui_simu(ievent)%islice_selected_rec(irec) == myrank) then
                irec_local=irec_local+1
                acqui_simu(ievent)%data_traces(irec_local,:,:)=Gather(irec, :, :)
             endif
          enddo
       endif

       ! send gather to other MPI slices
       do irank = 1, NPROC-1

          if (myrank == 0) then !! then send

             ! count the receiver in slice irank
             nsta_irank=0
             do irec = 1,  NSTA
                if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
             enddo

             ! if there is receiver in slice irank then MPI send data
             if (nsta_irank > 0) then
                allocate(Gather_loc(nsta_irank,Nt,NDIM))  !! data to send
                irec_local=0
                do irec = 1, NSTA
                   if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) then
                      irec_local = irec_local + 1
                      Gather_loc(irec_local, :, :) = Gather(irec, :, :) !! store data to send
                   endif
                enddo
                  if (DEBUG_MODE) write(IIDD,*) 'myrank ', myrank , 'send to ', irank, ' :' , nsta_irank, Nt
                tag    = 2001
                call MPI_SEND(Gather_loc, Nt*nsta_irank*NDIM, CUSTOM_MPI_TYPE, irank, tag, my_local_mpi_comm_world, ier)

                deallocate(Gather_loc)

             endif

          else !! then receive gather

             if (myrank == irank .and. acqui_simu(ievent)%nsta_slice > 0) then
                NSTA_LOC=acqui_simu(ievent)%nsta_slice
                Nt=acqui_simu(ievent)%Nt_data
                allocate(Gather_loc(NSTA_LOC,Nt,NDIM),acqui_simu(ievent)%data_traces(NSTA_LOC,Nt,NDIM), &
                     acqui_simu(ievent)%adjoint_sources(NDIM, NSTA_LOC, Nt), acqui_simu(ievent)%weight_trace(NDIM, NSTA_LOC))

                if (VERBOSE_MODE .or. DEBUG_MODE) allocate(acqui_simu(ievent)%synt_traces(NDIM, NSTA_LOC, Nt))

                if (DEBUG_MODE) write(IIDD,*) 'myrank ',myrank,' wait for 0 :', NSTA_LOC,Nt
                tag   = MPI_ANY_TAG
                call MPI_RECV(Gather_loc,Nt*NSTA_LOC*NDIM,CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, status,  ier)
                !! store in acqui_simu
                acqui_simu(ievent)%data_traces(:,:,:)=Gather_loc(:,:,:)
                deallocate(Gather_loc)
             endif

          endif


       enddo

       if (myrank == 0) deallocate(Gather)

       call synchronize_all()

       !! set other parameters (in futrue work need to read any additional files)

       !! set frequency to invert
       acqui_simu(ievent)%freqcy_to_invert(:,1,:)=fl
       acqui_simu(ievent)%freqcy_to_invert(:,2,:)=fh

       acqui_simu(ievent)%fl_event=fl
       acqui_simu(ievent)%fh_event=fh

    enddo

    if (myrank == 0) write(INVERSE_LOG_FILE,'(a25//)') '... reading data : passed'

  end subroutine read_data_gather
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! read acquisiton input file and store acqui_simu type
!----------------------------------------------------------------
  subroutine read_acqui_file(acqui_file, acqui_simu, myrank)


    character(len=MAX_LEN_STRING),                  intent(in) ::  acqui_file
    integer,                                        intent(in) ::  myrank
    type(acqui),  dimension(:), allocatable,     intent(inout) ::  acqui_simu

    ! locals
    character(len=MAX_LEN_STRING)                              :: line, keyw, line_to_read
    integer                                                    :: ipos0, ipos1, ievent
    integer                                                    :: ier
    if (DEBUG_MODE)  write(IIDD,*) '       MYGROUP  ', mygroup, '    MYRANK ', myrank
    NEVENT=0
    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '     READING acquisition'
       write(INVERSE_LOG_FILE,*)
    endif

    !! 1/ read to count the number of events
    open(666,file=trim(acqui_file),iostat=ier)
    if (ier /= 0) then
       write(*,*) ' error opening  ', trim(acqui_file), ' mygoup ', mygroup
    else
       if (DEBUG_MODE) write(IIDD,*) ' opening  ', trim(acqui_file), ' mygoup ', mygroup
    endif
    do

       read(666,'(a)',end=99) line
       if (DEBUG_MODE) write(IIDD,'(a)') trim(line)
       if (is_blank_line(line)) cycle                 !! no significant line
       if (INDEX(line,'event_name') > 0) NEVENT=NEVENT+1  !! new event

    enddo
99  close(666)

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*) '       ALLOCATE  acquisition structure for ', NEVENT, ' events '
       write(INVERSE_LOG_FILE,*)
    endif

    !! 2/ allocate and store type(acqui) acqui_simu
    if (NEVENT > 0) then
       allocate(acqui_simu(NEVENT))
    else
        allocate(acqui_simu(1))
       write(*,*) 'ERROR NO EVENTS FOUND IN ACQUISITION FILE ',myrank, mygroup, trim(acqui_file)
       stop
    endif

    ! open event file
    open(666,file=trim(acqui_file))
    ievent=0
    do    !! loop on all lines

        !! READ AND STORE ALL ITEM RELATED TO EVENT -------------------------------
        do
           read(666,'(a)',end=999) line
           if (is_blank_line(line)) cycle

           !! INDICES TO READ line -----------------------------------------------
           ipos0=index(line,':')+1
           ipos1=index(line,'#')-1
           if (ipos1 < 0 ) ipos1=len_trim(line)

           !! STORE KEYWORD ITEM -------------------------------------------------
           keyw=trim(adjustl(line(1:ipos0-2)))

           !! DIFFERENT ITEM TO READ ---------------------------------------------
           select case (trim(keyw))

              case('event_name')
                 ievent=ievent+1
                 acqui_simu(ievent)%event_name=trim(adjustl(line(ipos0:ipos1)))
                 acqui_simu(ievent)%nevent_tot=NEVENT
                 call store_default_acqui_values(acqui_simu, ievent)

              case('moment', 'force')
                 acqui_simu(ievent)%source_file=trim(adjustl(line(ipos0:ipos1)))
                 acqui_simu(ievent)%source_type=trim(adjustl(keyw))
                 acqui_simu(ievent)%adjoint_source_type='L2_OIL_INDUSTRY'

              case('axisem', 'dsm', 'plane', 'fk')
                 acqui_simu(ievent)%source_file=trim(adjustl(line(ipos0:ipos1)))
                 acqui_simu(ievent)%source_type=trim(adjustl(keyw))
                 acqui_simu(ievent)%adjoint_source_type='L2_FWI_TELESEISMIC'

              case('traction_dir')
                 acqui_simu(ievent)%traction_dir=trim(adjustl(line(ipos0:ipos1)))

              case('station_file')
                 acqui_simu(ievent)%station_file=trim(adjustl(line(ipos0:ipos1)))

              case('data_file')
                 acqui_simu(ievent)%data_file_gather=trim(adjustl(line(ipos0:ipos1)))

              case ('component')
                 !! initiallize components
                 acqui_simu(ievent)%component(1)='  '
                 acqui_simu(ievent)%component(2)='  '
                 acqui_simu(ievent)%component(3)='  '
                 !need to add '00' in case of mising components (because fortran cannot let the missing component to '  ')
                 line_to_read=line(ipos0:ipos1)//' 00 00 00'

                 read(line_to_read,*)  acqui_simu(ievent)%component(1), &
                                            acqui_simu(ievent)%component(2), &
                                            acqui_simu(ievent)%component(3)

                 if (myrank == 0) then
                    write(INVERSE_LOG_FILE,*) 'event', ievent,  ' components : ', &
                         trim(acqui_simu(ievent)%component(1)),' ', &
                         trim(acqui_simu(ievent)%component(2)),' ', &
                         trim(acqui_simu(ievent)%component(3))
                 endif

              case('source_wavelet')
                 acqui_simu(ievent)%source_wavelet_file=trim(adjustl(line(ipos0:ipos1)))
                 acqui_simu(ievent)%external_source_wavelet=.true.

              case ('NSTEP')
                 read(line(ipos0:ipos1),*) acqui_simu(ievent)%Nt_data

              case ('DT')
                 read(line(ipos0:ipos1),*) acqui_simu(ievent)%dt_data

              case default
                 write(*,*) 'ERROR KEY WORD NOT MATCH : ', trim(keyw), ' in file ', trim(acqui_file)
                 exit

           end select

       enddo

    enddo

999 close(666)

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '     READING acquisition passed '
       write(INVERSE_LOG_FILE,*)
    endif

  end subroutine read_acqui_file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! read inversion input file and store acqui_simu type
!----------------------------------------------------------------
  subroutine read_inver_file(inver_file, acqui_simu, inversion_param, myrank)


    character(len=MAX_LEN_STRING),                  intent(in) ::  inver_file
    integer,                                        intent(in) ::  myrank
    type(acqui),  dimension(:), allocatable,     intent(inout) ::  acqui_simu
    type(inver),                                 intent(inout) ::  inversion_param

    ! locals
    character(len=MAX_LEN_STRING)                              :: line, keyw
    integer                                                    :: ipos0, ipos1


    open(666, file=trim(inver_file))
    do
       read(666,'(a)',end=99) line
       if (is_blank_line(line)) cycle

       !! INDICES TO READ line -----------------------------------------------
       ipos0=index(line,':')+1
       ipos1=index(line,'#')-1
       if (ipos1 < 0 ) ipos1=len_trim(line)

       !! STORE KEYWORD ITEM -------------------------------------------------
       keyw=trim(adjustl(line(1:ipos0-2)))

       !! DIFFERENT ITEM TO READ ---------------------------------------------
       select case (trim(keyw))

       case ('Niter')
          read(line(ipos0:ipos1),*)  inversion_param%Niter

       case('Niter_wolfe')
          read(line(ipos0:ipos1),*)  inversion_param%Niter_wolfe

       case('max_history_bfgs')
          read(line(ipos0:ipos1),*)  inversion_param%max_history_bfgs

       case('max_relative_pert')
           read(line(ipos0:ipos1),*)  inversion_param%max_relative_pert

       case('param_family')
          read(line(ipos0:ipos1),*) inversion_param%param_family

       case('fl')
          read(line(ipos0:ipos1),*) fl

       case('fh')
          read(line(ipos0:ipos1),*) fh

       case('input_sem_model')
          read(line(ipos0:ipos1),*)  inversion_param%input_sem_model

       case('input_sem_prior')
          read(line(ipos0:ipos1),*)  inversion_param%input_sem_prior

       case('output_model')
          read(line(ipos0:ipos1),*)  inversion_param%output_model

       case('input_fd_model')
          read(line(ipos0:ipos1),*)  inversion_param%input_fd_model

       case ('taper')
          inversion_param%use_taper=.true.
          read(line(ipos0:ipos1),*) inversion_param%xmin_taper, inversion_param%xmax_taper, &
               inversion_param%ymin_taper, inversion_param%ymax_taper, &
               inversion_param%zmin_taper, inversion_param%zmax_taper

       case('shin_precond')
          read(line(ipos0:ipos1),*) inversion_param%shin_precond

       case('energy_precond')
          read(line(ipos0:ipos1),*) inversion_param%energy_precond

       case('z2_precond')
          read(line(ipos0:ipos1),*) inversion_param%z2_precond

       case('z_precond')
          read(line(ipos0:ipos1),*) inversion_param%aPrc, &
                                    inversion_param%zPrc1, &
                                    inversion_param%zPrc2

          inversion_param%z_precond=.true.

       
       case('relat_grad')
          read(line(ipos0:ipos1),*) inversion_param%relat_grad

       case('relat_cost')
          read(line(ipos0:ipos1),*) inversion_param%relat_cost

       case('dump_model_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_model_at_each_iteration

       case('dump_gradient_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_gradient_at_each_iteration

       case('dump_descent_direction_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_descent_direction_at_each_iteration

       case('use_tk_fd_regularisation')
          read(line(ipos0:ipos1),*) inversion_param%weight_Tikonov
          inversion_param%use_regularisation_FD_Tikonov=.true.
   
       case('use_tk_sem_regularisation')
          allocate(inversion_param%smooth_weight(inversion_param%NinvPar))
          read(line(ipos0:ipos1),*) inversion_param%smooth_weight(1), &
                                    inversion_param%smooth_weight(2), &
                                    inversion_param%smooth_weight(3)

          inversion_param%use_regularisation_SEM_Tikonov=.true.

       case('use_tk_sem_damping')
          allocate(inversion_param%damp_weight(inversion_param%NinvPar))
          read(line(ipos0:ipos1),*) inversion_param%damp_weight(1), &
                                    inversion_param%damp_weight(2), &
                                    inversion_param%damp_weight(3)
          inversion_param%use_damping_SEM_Tikonov=.true.

       case default
          write(*,*) 'ERROR KEY WORD NOT MATCH : ', trim(keyw), ' in file ', trim(inver_file)
          exit

       end select

    enddo

99  close(666)

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '     READ  ', trim(inver_file)
       write(INVERSE_LOG_FILE,*) '     Nb tot events ', acqui_simu(1)%nevent_tot
       write(INVERSE_LOG_FILE,*)
    endif

    if (VERBOSE_MODE .or. DEBUG_MODE) then
       inversion_param%dump_model_at_each_iteration=.true.
       inversion_param%dump_gradient_at_each_iteration=.true.
       inversion_param%dump_descent_direction_at_each_iteration=.true.
    endif

  end subroutine read_inver_file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------
! define default values
!----------------------------------------------------------------
  subroutine store_default_acqui_values(acqui_simu, ievent)
    type(acqui), allocatable, dimension(:), intent(inout)  :: acqui_simu
    integer,                                intent(in)     :: ievent

    acqui_simu(ievent)%source_file  ='none'
    acqui_simu(ievent)%traction_dir ='none'
    acqui_simu(ievent)%data_file_gather ='none'
    acqui_simu(ievent)%source_type  ='none'
    acqui_simu(ievent)%station_file ='none'
    acqui_simu(ievent)%adjoint_source_type='none'

    acqui_simu(ievent)%Xs=0.d0
    acqui_simu(ievent)%Ys=0.d0
    acqui_simu(ievent)%Zs=0.d0

    acqui_simu(ievent)%Mxx=0.d0
    acqui_simu(ievent)%Myy=0.d0
    acqui_simu(ievent)%Mzz=0.d0
    acqui_simu(ievent)%Mxy=0.d0
    acqui_simu(ievent)%Mxz=0.d0
    acqui_simu(ievent)%Myz=0.d0

    acqui_simu(ievent)%Fx=0.d0
    acqui_simu(ievent)%Fy=0.d0
    acqui_simu(ievent)%Fz=0.d0

    acqui_simu(ievent)%t_shift=0.d0
    acqui_simu(ievent)%hdur=0.d0

    acqui_simu(ievent)%nsources_local=0

  end subroutine store_default_acqui_values

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! true if blank line or commented line (ie that begins with #)
!----------------------------------------------------------------
  logical function is_blank_line(line)
    character(len=MAX_LEN_STRING), intent(in) :: line
    is_blank_line=.false.
    if (len(trim(adjustl(line))) == 0) is_blank_line = .true.
    if (INDEX(trim(adjustl(line)),'#') == 1) is_blank_line = .true.
  end function is_blank_line
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! true if blank line or commented line (ie that begins with #)
!----------------------------------------------------------------
  subroutine remove_blank_in_line(line, line_without_blank)
    character(len=MAX_LEN_STRING), intent(in)    :: line
    character(len=MAX_LEN_STRING), intent(inout) :: line_without_blank
    integer                                      :: n, i, k

    n=len_trim(adjustl(line))
    k=0
    do i=1, n
       if (line(i:i+1) == ' ') then
          ! it's blanc character then do nothing
       else
          k=k+1
          line_without_blank(k:k+1)=line(i:i+1)
       endif
    enddo
  end subroutine remove_blank_in_line
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------
! read station file in specfem format
!------------------------------------------------------------
  subroutine get_stations(acqui_simu)

    use constants, only: NDIM

    type(acqui), allocatable, dimension(:), intent(inout)  :: acqui_simu

    integer                                                :: ievent, irec, nsta, nrec_loc
    character(len=MAX_LEN_STRING)                          :: rec_filename,filtered_rec_filename
    write(INVERSE_LOG_FILE,*)
    write(INVERSE_LOG_FILE,*) '     READING stations '

    INVERSE_FWI_FULL_PROBLEM = .true.

    ! loop on all events
    do ievent = 1, NEVENT
       ! open station file
       rec_filename=trim(adjustl(acqui_simu(ievent)%station_file))
       filtered_rec_filename=rec_filename(1:len_trim(rec_filename))//'_FILTERED'

       call station_filter(rec_filename,filtered_rec_filename,nsta)

       acqui_simu(ievent)%nsta_tot=nsta
       allocate(acqui_simu(ievent)%station_name(nsta),acqui_simu(ievent)%network_name(nsta))
       allocate(acqui_simu(ievent)%position_station(3,nsta))
       allocate(acqui_simu(ievent)%xi_rec(nsta), &
                acqui_simu(ievent)%eta_rec(nsta), &
                acqui_simu(ievent)%gamma_rec(nsta))
       allocate(acqui_simu(ievent)%islice_selected_rec(nsta), &
                acqui_simu(ievent)%ispec_selected_rec(nsta), &
                acqui_simu(ievent)%number_receiver_global(nsta), &
                acqui_simu(ievent)%nu(NDIM,NDIM,nsta))

       ! reads STATIONS_FILTERED file, locates receivers in the mesh and compute Lagrange interpolators
       call locate_receivers(filtered_rec_filename,nsta,acqui_simu(ievent)%islice_selected_rec, &
                             acqui_simu(ievent)%ispec_selected_rec, &
                             acqui_simu(ievent)%xi_rec,acqui_simu(ievent)%eta_rec,acqui_simu(ievent)%gamma_rec, &
                             acqui_simu(ievent)%station_name,acqui_simu(ievent)%network_name,acqui_simu(ievent)%nu,1.0d0,1.0d0)

       nrec_loc = 0
       do irec = 1, nsta
         if (myrank == acqui_simu(ievent)%islice_selected_rec(irec)) then
             nrec_loc=nrec_loc+1
             acqui_simu(ievent)%number_receiver_global(nrec_loc)=irec
         endif
       enddo

       acqui_simu(ievent)%nsta_slice=nrec_loc

       if (acqui_simu(ievent)%nsta_slice > 0) then
          allocate(acqui_simu(ievent)%hxi    (NGLLX,nrec_loc))
          allocate(acqui_simu(ievent)%heta   (NGLLY,nrec_loc))
          allocate(acqui_simu(ievent)%hgamma (NGLLZ,nrec_loc))
          allocate(acqui_simu(ievent)%hpxi   (NGLLX,nrec_loc))
          allocate(acqui_simu(ievent)%hpeta  (NGLLY,nrec_loc))
          allocate(acqui_simu(ievent)%hpgamma(NGLLZ,nrec_loc))
          allocate(acqui_simu(ievent)%freqcy_to_invert(NDIM,2,nrec_loc))
       else
          allocate(acqui_simu(ievent)%hxi    (1,1))
          allocate(acqui_simu(ievent)%heta   (1,1))
          allocate(acqui_simu(ievent)%hgamma (1,1))
          allocate(acqui_simu(ievent)%hpxi   (1,1))
          allocate(acqui_simu(ievent)%hpeta  (1,1))
          allocate(acqui_simu(ievent)%hpgamma(1,1))
          allocate(acqui_simu(ievent)%freqcy_to_invert(NDIM,2,1))
       endif

       nrec_loc = 0
       do irec=1, acqui_simu(ievent)%nsta_tot
          if (myrank == acqui_simu(ievent)%islice_selected_rec(irec)) then

             nrec_loc = nrec_loc + 1

             ! compute Lagrange polynomials at the receiver location
             call lagrange_any(acqui_simu(ievent)%xi_rec(irec),NGLLX,xigll, &
                               acqui_simu(ievent)%hxi(1,nrec_loc),acqui_simu(ievent)%hpxi(1,nrec_loc))
             call lagrange_any(acqui_simu(ievent)%eta_rec(irec),NGLLY,yigll, &
                               acqui_simu(ievent)%heta(1,nrec_loc),acqui_simu(ievent)%hpeta(1,nrec_loc))
             call lagrange_any(acqui_simu(ievent)%gamma_rec(irec),NGLLZ,zigll, &
                               acqui_simu(ievent)%hgamma(1,nrec_loc),acqui_simu(ievent)%hpgamma(1,nrec_loc))

          endif
       enddo

    enddo

    write(INVERSE_LOG_FILE,*) '     READING stations passed '
    write(INVERSE_LOG_FILE,*)

 end subroutine get_stations

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------
! read source parameter file
!------------------------------------------------------------
  subroutine get_point_source(acqui_simu)
    type(acqui), allocatable, dimension(:), intent(inout)  :: acqui_simu
    ! locals
    character(len=MAX_LEN_STRING)                          :: string
    integer                                                :: ievent, ier
    integer                                                :: i
    real(kind=CUSTOM_REAL)                                 :: dt_dummy

    write(INVERSE_LOG_FILE,*)
    write(INVERSE_LOG_FILE,*) '     READING sources parameters '

    do ievent=1,acqui_simu(1)%nevent_tot

       select case (acqui_simu(ievent)%source_type)

       case('moment')
          open(IINN,file=trim(acqui_simu(ievent)%source_file),status='old',action='read',iostat=ier)
          if (ier /= 0) then
             print *,'Error opening file: ',trim(acqui_simu(ievent)%source_file)
             stop 'Error opening CMTSOLUTION file'
          endif


          ! gets header line
          read(IINN,"(a256)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading header line in event ',ievent
             stop 'Error reading header line in station in CMTSOLUTION file'
          endif

          ! skips empty lines
          do while (len_trim(string) == 0)
             read(IINN,"(a256)",iostat=ier) string
             if (ier /= 0) then
                print *, 'Error reading header line in event ',ievent
                stop 'Error reading header line in station in CMTSOLUTION file'
             endif
          enddo

          ! ignore line with event name
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading event name in event ',ievent
             stop 'Error reading event name in station in CMTSOLUTION file'
          endif

          ! read time shift
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading time shift in event ',ievent
             stop 'Error reading time shift in station in CMTSOLUTION file'
          endif
          read(string(12:len_trim(string)),*) acqui_simu(ievent)%t_shift
          !write(*,*) 'read t_shift ' , acqui_simu(ievent)%t_shift
          ! read half duration
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading half duration in event ',ievent
             stop 'Error reading half duration in station in CMTSOLUTION file'
          endif
          read(string(15:len_trim(string)),*) acqui_simu(ievent)%hdur
          !write(*,*) 'read hdur ' , acqui_simu(ievent)%hdur
          ! read latitude
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading latitude in event ',ievent
             stop 'Error reading latitude in station in CMTSOLUTION file'
          endif
          read(string(10:len_trim(string)),*) acqui_simu(ievent)%Ys
          !write(*,*) 'read Ys ' , acqui_simu(ievent)%Ys
          ! read longitude
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading longitude in event ',ievent
             stop 'Error reading longitude in station in CMTSOLUTION file'
          endif
          read(string(11:len_trim(string)),*) acqui_simu(ievent)%Xs
          !write(*,*) 'read Xs ' , acqui_simu(ievent)%Xs
          ! read depth
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading depth in event ',ievent
             stop 'Error reading depth in station in CMTSOLUTION file'
          endif
          read(string(7:len_trim(string)),*) acqui_simu(ievent)%Zs
          !write(*,*) 'read Zs ' , acqui_simu(ievent)%Zs

          ! seismic moment tensor
          ! CMTSOLUTION: components given in dyne-cm
          ! read Mrr
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading Mrr in event ',ievent
             stop 'Error reading Mrr in station in CMTSOLUTION file'
          endif
          read(string(5:len_trim(string)),*)  acqui_simu(ievent)%Mzz
          !write(*,*)  acqui_simu(ievent)%Mzz
          ! read Mtt
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading Mtt in event ',ievent
             stop 'Error reading Mtt in station in CMTSOLUTION file'
          endif
          read(string(5:len_trim(string)),*) acqui_simu(ievent)%Myy
          !write(*,*)  acqui_simu(ievent)%Myy
          ! read Mpp
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading Mpp in event ',ievent
             stop 'Error reading Mpp in station in CMTSOLUTION file'
          endif
          read(string(5:len_trim(string)),*) acqui_simu(ievent)%Mxx
          !write(*,*)  acqui_simu(ievent)%Mxx
          ! read Mrt
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading Mrt in event ',ievent
             stop 'Error reading Mrt in station in CMTSOLUTION file'
          endif
          read(string(5:len_trim(string)),*) acqui_simu(ievent)%Myz
          acqui_simu(ievent)%Myz = - acqui_simu(ievent)%Myz
          !write(*,*)  acqui_simu(ievent)%Myz
          ! read Mrp
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading Mrp in event ',ievent
             stop 'Error reading Mrp in station in CMTSOLUTION file'
          endif
          read(string(5:len_trim(string)),*) acqui_simu(ievent)%Mxz
          !write(*,*)  acqui_simu(ievent)%Mxz
          ! read Mtp
          read(IINN,"(a)",iostat=ier) string
          if (ier /= 0) then
             print *, 'Error reading Mtp in event ',ievent
             stop 'Error reading Mtp in station in CMTSOLUTION file'
          endif
          read(string(5:len_trim(string)),*) acqui_simu(ievent)%Mxy
          acqui_simu(ievent)%Mxy = - acqui_simu(ievent)%Mxy
          !write(*,*)  acqui_simu(ievent)%Mxy
          close(IINN)

          ! to be consistent with specfem
          acqui_simu(ievent)%Mxx = acqui_simu(ievent)%Mxx * 1.d-7
          acqui_simu(ievent)%Myy = acqui_simu(ievent)%Myy * 1.d-7
          acqui_simu(ievent)%Mzz = acqui_simu(ievent)%Mzz * 1.d-7
          acqui_simu(ievent)%Mxy = acqui_simu(ievent)%Mxy * 1.d-7
          acqui_simu(ievent)%Mxz = acqui_simu(ievent)%Mxz * 1.d-7
          acqui_simu(ievent)%Myz = acqui_simu(ievent)%Myz * 1.d-7

          if (acqui_simu(ievent)%external_source_wavelet) then
             allocate(acqui_simu(ievent)%source_wavelet(acqui_simu(ievent)%Nt_data,1))
             open(IINN, file=trim(acqui_simu(ievent)%source_wavelet_file))
             do i=1,acqui_simu(ievent)%Nt_data
                read(IINN, *) dt_dummy, acqui_simu(ievent)%source_wavelet(i,1)
             enddo
             close(IINN)
          endif

       case('force')
          print *, 'Abort not implemented yet : FORCESOLUTION in event ',ievent
          stop

       case default
          !! nothing to do

       end select


    enddo

    write(INVERSE_LOG_FILE,*) '     READING sources parameters passed '

  end subroutine get_point_source

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------
! MPI BCAST of acqui_simu
!---------------------------------------------------------------

  subroutine bcast_all_acqui(acqui_simu, inversion_param, myrank)
    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem
    type(acqui), allocatable, dimension(:), intent(inout)  :: acqui_simu
    type(inver),                            intent(inout)  :: inversion_param
    integer,                                intent(in)     :: myrank
    integer                                                :: nsta_tot, i,ier

    if (myrank == 0) NEVENT=acqui_simu(1)%nevent_tot
    call MPI_BCAST(NEVENT,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

    if (myrank > 0) allocate(acqui_simu(NEVENT))

    do i=1,NEVENT

       ! events
       acqui_simu(i)%nevent_tot=NEVENT
       ! bcast file paths
       call MPI_BCAST(acqui_simu(i)%source_file,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%traction_dir,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%data_file_gather,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%source_type,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%adjoint_source_type,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%station_file,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%event_name,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%component, 6,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

       call MPI_BCAST(acqui_simu(i)%Nt_data,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%dt_data,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

       ! bcast source parameters
       call MPI_BCAST(acqui_simu(i)%Xs,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Ys,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Zs,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

       call MPI_BCAST(acqui_simu(i)%Mxx,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Myy,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Mzz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Mxy,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Mxz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Myz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

       call MPI_BCAST(acqui_simu(i)%Fx,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Fy,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%Fz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

       call MPI_BCAST(acqui_simu(i)%t_shift,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%hdur,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%source_wavelet_file, MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

       call MPI_BCAST(acqui_simu(i)%external_source_wavelet, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
       if (acqui_simu(i)%external_source_wavelet) then
          if (myrank > 0) allocate(acqui_simu(i)%source_wavelet(acqui_simu(i)%Nt_data,1))
           call MPI_BCAST(acqui_simu(i)%source_wavelet,acqui_simu(i)%Nt_data, CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       endif

       ! stations
       if (myrank == 0) nsta_tot=acqui_simu(i)%nsta_tot
       call  MPI_BCAST(nsta_tot,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
       acqui_simu(i)%nsta_tot=nsta_tot
       if (myrank > 0) then
          allocate(acqui_simu(i)%station_name(nsta_tot),acqui_simu(i)%network_name(nsta_tot))
          allocate(acqui_simu(i)%position_station(3,nsta_tot))
       endif

       call MPI_BCAST(acqui_simu(i)%station_name, nsta_tot* MAX_LENGTH_STATION_NAME,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%network_name, nsta_tot* MAX_LENGTH_NETWORK_NAME,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(i)%position_station, 3*nsta_tot,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    enddo

    ! inverse params
    call MPI_BCAST(inversion_param%Niter,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%Niter_wolfe,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%max_history_bfgs,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(fl,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(fh,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%max_relative_pert,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%relat_grad,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%aPrc,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zPrc1,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zPrc2,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%relat_cost,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%weight_Tikonov,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_damping_SEM_Tikonov,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_regularisation_FD_Tikonov,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_regularisation_SEM_Tikonov,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

    if (myrank >0 .and. inversion_param%use_regularisation_SEM_Tikonov ) then 
       allocate(inversion_param%smooth_weight(inversion_param%NinvPar))
    end if
    if (inversion_param%use_regularisation_SEM_Tikonov) &
    call MPI_BCAST(inversion_param%smooth_weight(1),inversion_param%NinvPar,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

    if (myrank >0 .and. inversion_param%use_damping_SEM_Tikonov ) then 
       allocate(inversion_param%damp_weight(inversion_param%NinvPar))
    end if
    if (inversion_param%use_damping_SEM_Tikonov) &
    call MPI_BCAST(inversion_param%damp_weight(1),inversion_param%NinvPar,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

    call MPI_BCAST(inversion_param%output_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%input_fd_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%input_sem_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%input_sem_prior,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_taper, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%shin_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%energy_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%z2_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%z_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier) 
    call MPI_BCAST(inversion_param%dump_model_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%dump_gradient_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%dump_descent_direction_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%param_family, MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
    
    ! user taper on gradient (MASK)
    call MPI_BCAST(inversion_param%xmin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%xmax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%ymin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%ymax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zmin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zmax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    
   

  end subroutine bcast_all_acqui

end module input_output
