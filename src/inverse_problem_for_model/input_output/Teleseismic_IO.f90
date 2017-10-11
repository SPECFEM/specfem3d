module Teleseismic_IO_mod

  ! from specfem
  use constants, only: mygroup
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NGNOD

  ! from inversion
  use inverse_problem_par
  use mesh_tools
  use passive_imaging_format_mod, only: gather, read_pif_header_file, read_binary_data,   &
                                        read_binary_source_signature, get_data_component, &
                                        calc_delta_dist_baz, calc_dist_baz_cart, lowcase
  integer, private :: NEVENT

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------------
!> reading event file information for teleseismic inversion
!-------------------------------------------------------------------------------------------------------------------
  subroutine read_acqui_teleseismic_file(acqui_file, acqui_simu, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    character(len=MAX_LEN_STRING),                  intent(in) ::  acqui_file
    integer,                                        intent(in) ::  myrank
    type(acqui),  dimension(:), allocatable,     intent(inout) ::  acqui_simu


    ! locals
    character(len=MAX_LEN_STRING)                              :: line, keyw, filename !, line_to_read
    integer                                                    :: ipos0, ipos1, ievent
    integer                                                    :: ier, nsta, nt, ncomp, ista

    double precision                                           :: baz, dist, gcarc
    
    type(gather), dimension(:), allocatable :: mygather
    
    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '     READING teleseismic acquisition file '
       write(INVERSE_LOG_FILE,*)
    endif

    NEVENT=0
    ! only master reads acqui file
    if (myrank == 0) then
       do
          !! 1/ read to count the number of events
          open(666,file=trim(acqui_file),iostat=ier)
          if (ier /= 0) then
             write(*,*) ' error opening  ', trim(acqui_file), ' mygoup ', mygroup
          else
             if (DEBUG_MODE) write(IIDD,*) ' opening  ', trim(acqui_file), ' mygoup ', mygroup
          endif

          read(666,'(a)',end=99) line
          if (DEBUG_MODE) write(IIDD,'(a)') trim(line)
          !if (is_blank_line(line)) cycle                     !! no significant line
          if (INDEX(line,'event_name') > 0) NEVENT=NEVENT+1  !! new event
       enddo
99     close(666)

       !! 2/ allocate and store type(acqui) acqui_simu
       if (NEVENT > 0) then
          allocate(acqui_simu(NEVENT))
          allocate(mygather(NEVENT))
       else
          allocate(acqui_simu(1))
          write(*,*) 'ERROR NO EVENTS FOUND IN ACQUISITION FILE ',myrank, mygroup, trim(acqui_file)
          stop
       endif

       ! open event file
       open(666,file=trim(acqui_file))
       ievent=0
       do
          read(666,'(a)',end=999) line
          !if (is_blank_line(line)) cycle
          ! INDICES TO READ line -----------------------------------------------
          ipos0=index(line,':')+1
          ipos1=index(line,'#')-1
          if (ipos1 < 0 ) ipos1=len_trim(line)

          !! STORE KEYWORD ITEM -------------------------------------------------
          keyw     = trim(adjustl(line(1:ipos0-2)))
          filename = trim(adjustl(line(ipos0:ipos1)))

          select case (trim(keyw))
          case('event_name')
             
             ievent=ievent+1

             !*** Read pif header file
          
             call read_pif_header_file(filename,mygather(ievent))
             
             !*** Fill acquisition structure
             ! Get primary informations
             nsta = mygather(ievent)%hdr%nsta
             nt   = mygather(ievent)%hdr%nt
             
             ! Fill primary informations from header
             acqui_simu(ievent)%nevent_tot           = NEVENT
             acqui_simu(ievent)%nsta_tot             = nsta
             acqui_simu(ievent)%Nt_data              = nt
             acqui_simu(ievent)%dt_data              = mygather(ievent)%hdr%dt
             acqui_simu(ievent)%event_name           = mygather(ievent)%hdr%event_name
             acqui_simu(ievent)%source_type_physical = mygather(ievent)%hdr%source_type
             acqui_simu(ievent)%source_type_modeling = mygather(ievent)%hdr%modeling_tool
             acqui_simu(ievent)%source_file          = mygather(ievent)%hdr%modeling_path
             acqui_simu(ievent)%Origin_chunk_lat     = mygather(ievent)%hdr%mesh_origin(1)
             acqui_simu(ievent)%Origin_chunk_lon     = mygather(ievent)%hdr%mesh_origin(2)
             acqui_simu(ievent)%Origin_chunk_azi     = mygather(ievent)%hdr%mesh_origin(3)
             acqui_simu(ievent)%is_time_pick         = mygather(ievent)%hdr%is_pick
             acqui_simu(ievent)%time_window          = mygather(ievent)%hdr%is_window
             acqui_simu(ievent)%time_before_pick     = mygather(ievent)%hdr%tbef
             acqui_simu(ievent)%time_after_pick      = mygather(ievent)%hdr%taft
             acqui_simu(ievent)%station_coord_system = mygather(ievent)%hdr%coord_sys
             acqui_simu(ievent)%source_wavelet_file  = mygather(ievent)%hdr%estimated_src

             ! Some checks about source wavelet and traction
             if (acqui_simu(ievent)%source_wavelet_file /= 'undef') then
                acqui_simu(ievent)%external_source_wavelet=.true.
                ! note sure if i should use this one.. 
                allocate(acqui_simu(ievent)%user_source_time_function(1,nt)) 
                call read_binary_source_signature(acqui_simu(ievent)%source_wavelet_file, &
                                                                                      nt, &
                                       acqui_simu(ievent)%user_source_time_function(1,:))
             end if
             
             select case (adjustl(trim(acqui_simu(ievent)%source_type_modeling)))
             case('axisem','dsm')
                acqui_simu(ievent)%traction_dir = mygather(ievent)%hdr%modeling_path
             end select

             ! Determine which components
             call get_data_component(mygather(ievent)%hdr%data_type, &
                                     mygather(ievent)%hdr%data_comp, &
                                                              ncomp, &
                                       acqui_simu(ievent)%component)
             
             ! Fill source informations
             acqui_simu(ievent)%event_lat   = mygather(ievent)%source%lat
             acqui_simu(ievent)%event_lon   = mygather(ievent)%source%lon
             acqui_simu(ievent)%event_depth = mygather(ievent)%source%ele
             acqui_simu(ievent)%xshot       = mygather(ievent)%source%x
             acqui_simu(ievent)%yshot       = mygather(ievent)%source%y
             !acqui_simu(ievent)%zshot       = mygather(ievent)%source%z  !use ele instead of z because it
             acqui_simu(ievent)%zshot       = mygather(ievent)%source%ele !is given wrt to earth surface
             
             ! Allocate stations array and fill arrays
             allocate(acqui_simu(ievent)%station_name(nsta))
             allocate(acqui_simu(ievent)%network_name(nsta))
             allocate(acqui_simu(ievent)%position_station(3,nsta))
             allocate(acqui_simu(ievent)%read_station_position(3,nsta)) !! actually geographic coord
             acqui_simu(ievent)%station_name(:)       = mygather(ievent)%stations(:)%name
             acqui_simu(ievent)%network_name(:)       = mygather(ievent)%stations(:)%ntwk
             acqui_simu(ievent)%position_station(1,:) = mygather(ievent)%stations(:)%x
             acqui_simu(ievent)%position_station(2,:) = mygather(ievent)%stations(:)%y
             ! again, i'm not sure here between z and ele... use ele for now
             ! acqui_simu(ievent)%position_station(3,:) = mygather(ievent)%stations(:)%z
             acqui_simu(ievent)%position_station(3,:) = mygather(ievent)%stations(:)%ele
             acqui_simu(ievent)%read_station_position(1,:) = mygather(ievent)%stations(:)%lat
             acqui_simu(ievent)%read_station_position(2,:) = mygather(ievent)%stations(:)%lon
             acqui_simu(ievent)%read_station_position(2,:) = mygather(ievent)%stations(:)%ele


             ! Use time picks if needed
             if (acqui_simu(ievent)%is_time_pick) then
                allocate(acqui_simu(ievent)%time_pick(nsta))
                acqui_simu(ievent)%time_pick(:) = mygather(ievent)%stations(:)%tpick
             end if

             ! Compute baz etc.
             allocate(acqui_simu(ievent)%baz(nsta))
             allocate(acqui_simu(ievent)%dist(nsta))
             allocate(acqui_simu(ievent)%gcarc(nsta)) !! not used great circle arc not used
             ! allocate((acqui_simu(ievent)%inc(nsta))   !! not used incidence angle
             select case(trim(adjustl(acqui_simu(ievent)%source_type_modeling)))
             case('pointsource') ! then local source
                do ista = 1, nsta
                   call calc_dist_baz_cart(mygather(ievent)%source%x,           &
                                           mygather(ievent)%source%y,           &
                                           mygather(ievent)%stations(ista)%x,   &
                                           mygather(ievent)%stations(ista)%y,   &       
                                           dist,                               & 
                                           baz)
                   acqui_simu(ievent)%baz(ista)  = real(baz,kind=custom_real)
                   acqui_simu(ievent)%dist(ista) = real(dist,kind=custom_real)
                end do
             case('fk','axisem','dsm') ! then teleseismic source
                do ista = 1, nsta
                   call calc_delta_dist_baz(mygather(ievent)%source%lat,           &
                                            mygather(ievent)%source%lon,           &
                                            mygather(ievent)%stations(ista)%lat,   &
                                            mygather(ievent)%stations(ista)%lon,   &
                                            gcarc,                                 &
                                            dist,                                  & 
                                            baz)
                   acqui_simu(ievent)%gcarc(ista) = gcarc
                   acqui_simu(ievent)%dist(ista)  = dist
                   acqui_simu(ievent)%baz(ista)   = baz
                end do
             end select

             ! Correct back-azimuth from mesh orientation (to check, depends on how we considere azi)
             !!acqui_simu(ievent)%baz(:) = acqui_simu(ievent)%baz(:) &
             !!                          - acqui_simu(ievent)%Origin_chunk_azi
             
          end select
       enddo

999    close(666)

    endif


    ! master broadcasts read values
    call mpi_bcast(nevent, 1, mpi_integer, 0, my_local_mpi_comm_world, ier)
    if (myrank > 0) allocate(acqui_simu(NEVENT))
    do ievent = 1, NEVENT

       ! broadcast integers
       call mpi_bcast(acqui_simu(ievent)%nevent_tot, 1,     mpi_integer, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%nsta_tot,   1,     mpi_integer, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%nt_data,    1,     mpi_integer, 0, &
            my_local_mpi_comm_world, ier)

       ! broadcast strings
       call mpi_bcast(acqui_simu(ievent)%event_name,           max_len_string, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%source_type_physical,            256, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%source_type_modeling,            256, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%source_file,          max_len_string, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%station_coord_system, max_len_string, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%source_wavelet_file,  max_len_string, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%traction_dir,         max_len_string, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%component,                         6, mpi_character, 0, &
            my_local_mpi_comm_world, ier) 

       ! broadcast reals
       call mpi_bcast(acqui_simu(ievent)%dt_data,             1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%origin_chunk_lat,    1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%origin_chunk_lon,    1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%origin_chunk_azi,    1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%time_before_pick,    1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%time_after_pick,     1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%event_lat,           1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%event_lon,           1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%event_depth,         1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%xshot,               1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%yshot,               1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%zshot,               1, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)

       
       ! broadcast logicals
       call mpi_bcast(acqui_simu(ievent)%is_time_pick,        1, mpi_logical, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%time_window,         1, mpi_logical, 0, &
            my_local_mpi_comm_world, ier)

       ! conditional broadcasts for wavelet
       if (acqui_simu(ievent)%source_wavelet_file /= 'undef') then
          ! Source time functions
          call mpi_bcast(acqui_simu(ievent)%external_source_wavelet, 1, mpi_logical, 0, &
               my_local_mpi_comm_world, ier)
          allocate(acqui_simu(ievent)%user_source_time_function(1,acqui_simu(ievent)%nt_data)) 
          call mpi_bcast(acqui_simu(ievent)%user_source_time_function, &
                         acqui_simu(ievent)%nt_data,                   &
                         custom_mpi_type, 0, my_local_mpi_comm_world, ier)
       end if

       ! conditional broadcast for allocatable arrays
       allocate(acqui_simu(ievent)%station_name(acqui_simu(ievent)%nsta_tot))
       allocate(acqui_simu(ievent)%network_name(acqui_simu(ievent)%nsta_tot))
       allocate(acqui_simu(ievent)%position_station(3,acqui_simu(ievent)%nsta_tot))
       allocate(acqui_simu(ievent)%read_station_position(3,acqui_simu(ievent)%nsta_tot))
       !! actually geographic coord
       nsta = acqui_simu(ievent)%nsta_tot

       call mpi_bcast(acqui_simu(ievent)%station_name, MAX_LENGTH_STATION_NAME*nsta, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%network_name, MAX_LENGTH_STATION_NAME*nsta, mpi_character, 0, &
            my_local_mpi_comm_world, ier)
       
       call mpi_bcast(acqui_simu(ievent)%position_station,      3*nsta, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%read_station_position, 3*nsta, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)

       if (acqui_simu(ievent)%is_time_pick) then
          allocate(acqui_simu(ievent)%time_pick(nsta))
          call mpi_bcast(acqui_simu(ievent)%time_pick,      nsta, custom_mpi_type, 0, &
               my_local_mpi_comm_world, ier)
       end if
       allocate(acqui_simu(ievent)%baz(nsta))
       allocate(acqui_simu(ievent)%dist(nsta))
       allocate(acqui_simu(ievent)%gcarc(nsta))
       call mpi_bcast(acqui_simu(ievent)%baz,               nsta, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%gcarc,             nsta, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
       call mpi_bcast(acqui_simu(ievent)%dist,              nsta, custom_mpi_type, 0, &
            my_local_mpi_comm_world, ier)
    enddo

    !! to do do not forget to bcast all acqui_simu structure to other MPI slices
    !! ......



  end subroutine read_acqui_teleseismic_file

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! master read waveform data gather and bcast to MPI slice concerned
!----------------------------------------------------------------
  subroutine read_pif_data_gather(acqui_simu, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                     intent(in)    :: myrank
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu

    integer                                                    :: ievent, idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                    :: tag, ier, nsta_irank, irank
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable      :: Gather, Gather_loc
    integer                                                    :: status(MPI_STATUS_SIZE)
    real(kind=CUSTOM_REAL)                                     :: dummy_real

    nb_traces_tot=0.

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

             ! First check if displacement, velocity, acceleration or pressure
             select case (lowcase(trim(acqui_simu(ievent)%component(idim))))
             case('d')

             case('v')

             case('a')
                
             case('p')

             end select

             ! Check available data
             select case (lowcase(trim(acqui_simu(ievent)%component(idim))))
             case('z')

             case('x','y')

             case('r','t')

             case('l','q')

             end select
             
             nb_traces_tot=nb_traces_tot+NSTA
             irec=irec+1
             read(IINN,rec=irec) Gather(:,:,idim)
             
             ! Rotate from original data coordinate system ((x,y,z),(zen),(rtz)) to mesh one (xyz)
             select case (lowcase(trim(acqui_simu(ievent)%component(idim))))
             case('xyz')  ! data are already in the mesh coordinate system
                
                ! Data rotation not required
                
             case('zen')  ! data are in standard coordinate system

                ! Data rotation required to pass in mesh system (zen -> xyz)
                call define_mesh_rotation_matrix(lat0,lon0,azi0)
                call rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)

             case('rtz')  !  dataare in the souce receiver coordinate system

                ! Data rotation required (baz-azi) (rtz -> zne)
                call rotate_ZRT_to_ZNE(vz2,vr,vt,vz,vn,ve,nrec,nt,bazi)
                
                ! Data rotation required to pass in mesh system (zen -> xyz)
                call define_mesh_rotation_matrix(lat0,lon0,azi0)
                call rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)
                
             case('lqt')  ! data are in the ray coordinate system

                ! Data rotation required (baz-azi and incidence angle) (rtz -> zen)
                call rotate_LQT_to_ZNE(vl,vq,vt,vz,vn,ve,nrec,nt,bazi,inci)
                
                ! Data rotation required to pass in mesh system (zen -> xyz)
                call define_mesh_rotation_matrix(lat0,lon0,azi0)
                call rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)
                
             end select


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
       !!acqui_simu(ievent)%freqcy_to_invert(:,1,:)=fl
       !!acqui_simu(ievent)%freqcy_to_invert(:,2,:)=fh

       !! get band pass filter values if needed
       if ( use_band_pass_filter) then
          acqui_simu(ievent)%Nfrq=NIFRQ
          acqui_simu(ievent)%band_pass_filter=use_band_pass_filter
          allocate(acqui_simu(ievent)%fl_event(acqui_simu(ievent)%Nfrq))
          allocate(acqui_simu(ievent)%fh_event(acqui_simu(ievent)%Nfrq))
          acqui_simu(ievent)%fl_event(:)=fl(:)
          acqui_simu(ievent)%fh_event(:)=fh(:)
          !! WARNING WARNING
          !! this is for telesismic case for now only one
          !! frequency is allowed (todo fix it)
          acqui_simu(ievent)%freqcy_to_invert(:,1,:)=fl(1)
          acqui_simu(ievent)%freqcy_to_invert(:,2,:)=fh(1)
       endif

    enddo

    call MPI_BCAST(nb_traces_tot,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    dummy_real=nb_traces_tot
    call sum_all_all_cr_for_simulatenous_runs(dummy_real,nb_traces_tot,1)

    if (myrank == 0) write(INVERSE_LOG_FILE,'(a25//)') '... reading data : passed'

  end subroutine read_data_gather
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------------
!> store arrays that needed for specfem to be able to use stations in mesh
!-------------------------------------------------------------------------------------------------------------------
  subroutine setup_teleseismic_stations(acqui_simu, myrank)

    !! from  acqui_simu%position_station assumed to be in Cartesian coordinate
    !!
    !! compute all arrays specific to specfem :
    !!
    !! xi_rec,eta_rec,gamma_rec
    !! islice_selected_rec
    !! ispec_selected_rec
    !! nsta_slice
    !! number_receiver_global
    !! nu
    !! hxi, hpxi
    !! heta, hpeta
    !! hgamma, hpgamma
    !!
    !! ------------------------------------------------------------------------------

    type(acqui), allocatable, dimension(:), intent(inout)  :: acqui_simu
    integer,                                intent(in)     :: myrank

    integer                                                :: ievent, ireceiver, nsta_slice, irec_local, NSTA, NEVENT
    integer                                                :: ispec_selected, islice_selected
    double precision                                       :: xi_receiver, eta_receiver, gamma_receiver
    double precision                                       :: x_found,  y_found,  z_found
    double precision                                       :: x_to_locate, y_to_locate, z_to_locate
    real(kind=CUSTOM_REAL)                                 :: distance_min_glob,distance_max_glob
    real(kind=CUSTOM_REAL)                                 :: elemsize_min_glob,elemsize_max_glob
    real(kind=CUSTOM_REAL)                                 :: x_min_glob,x_max_glob
    real(kind=CUSTOM_REAL)                                 :: y_min_glob,y_max_glob
    real(kind=CUSTOM_REAL)                                 :: z_min_glob,z_max_glob
    integer,                 dimension(NGNOD)              :: iaddx,iaddy,iaddz
    double precision,        dimension(NGLLX)              :: hxis,hpxis
    double precision,        dimension(NGLLY)              :: hetas,hpetas
    double precision,        dimension(NGLLZ)              :: hgammas,hpgammas
    double precision                                       :: distance_from_target


    ! get mesh properties (mandatory before calling locate_point_in_mesh)
    call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
         x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
         elemsize_min_glob,elemsize_max_glob, &
         distance_min_glob,distance_max_glob)
    if (DEBUG_MODE) then
       write (IIDD, *)
       write (IIDD, *) ' LOCATE TELESEISMIC STATIONS ---------------------'
       write (IIDD, *)
    endif

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) ' ... locate stations in specfem mesh :'
       write(INVERSE_LOG_FILE,*)
       call flush_iunit(INVERSE_LOG_FILE)
    endif

    NEVENT=acqui_simu(1)%nevent_tot
    do ievent = 1, NEVENT
       if (myrank == 0) then
            write(INVERSE_LOG_FILE,*) ' ... locate stations for source ', ievent
            call flush_iunit(INVERSE_LOG_FILE)
       endif
       NSTA = acqui_simu(ievent)%nsta_tot
       allocate(acqui_simu(ievent)%xi_rec(NSTA), &
                acqui_simu(ievent)%eta_rec(NSTA), &
                acqui_simu(ievent)%gamma_rec(NSTA))

       allocate(acqui_simu(ievent)%islice_selected_rec(NSTA), &
                acqui_simu(ievent)%ispec_selected_rec(NSTA), &
                acqui_simu(ievent)%number_receiver_global(NSTA))

       acqui_simu(ievent)%number_receiver_global(:)=-1
       acqui_simu(ievent)%ispec_selected_rec(:)=-1
       acqui_simu(ievent)%islice_selected_rec(:)=-1

       nsta_slice=0
       do ireceiver = 1, NSTA

          x_to_locate = acqui_simu(ievent)%position_station(1,ireceiver)
          y_to_locate = acqui_simu(ievent)%position_station(2,ireceiver)
          z_to_locate = acqui_simu(ievent)%position_station(3,ireceiver)

          call get_point_in_mesh(x_to_locate, y_to_locate, z_to_locate, iaddx, iaddy, iaddz, elemsize_max_glob, &
               ispec_selected, xi_receiver, eta_receiver, gamma_receiver, x_found, y_found, z_found, myrank)

          call get_MPI_slice_and_bcast_to_all(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
               xi_receiver, eta_receiver, gamma_receiver, ispec_selected, islice_selected,  distance_from_target, myrank)

          acqui_simu(ievent)%xi_rec(ireceiver)=xi_receiver
          acqui_simu(ievent)%eta_rec(ireceiver)=eta_receiver
          acqui_simu(ievent)%gamma_rec(ireceiver)=gamma_receiver
          acqui_simu(ievent)%islice_selected_rec(ireceiver)=islice_selected
          acqui_simu(ievent)%ispec_selected_rec(ireceiver)=ispec_selected

          if (myrank == islice_selected) then
             nsta_slice=nsta_slice+1
             acqui_simu(ievent)%number_receiver_global(nsta_slice)=ireceiver
          endif

          if (DEBUG_MODE) then
             call flush_iunit(IIDD)
             if (distance_from_target > 1000. ) then
                write (IIDD, '(a29,f20.5,5x,i5)') 'WARNING : error location sta ', distance_from_target, ireceiver
                write (IIDD, '(a17,3f20.5)') 'desired position :', x_to_locate, y_to_locate, z_to_locate
                write (IIDD, '(a17,3f20.5)') 'position found   :', x_found, y_found, z_found
                write (IIDD, '(a15,3f13.5)') ' xi eta gamma ', xi_receiver, eta_receiver, gamma_receiver
                write (IIDD, *)
             else
                write (IIDD, '(a22,f20.5,5x,i5)') 'LOCATED STA : error location sta ', distance_from_target, ireceiver
                write (IIDD, '(a17,3f20.5)') 'desired position :', x_to_locate, y_to_locate, z_to_locate
                write (IIDD, '(a17,3f20.5)') 'position found   :', x_found, y_found, z_found
                write (IIDD, '(a15,3f13.5)') ' xi eta gamma ', xi_receiver, eta_receiver, gamma_receiver
                write (IIDD, *)
             endif

          endif

       enddo

       acqui_simu(ievent)%nsta_slice=nsta_slice

    enddo

   ! store local information
    do ievent = 1, NEVENT

       if (acqui_simu(ievent)%nsta_slice > 0) then
          allocate(acqui_simu(ievent)%hxi(NGLLX,acqui_simu(ievent)%nsta_slice))
          allocate(acqui_simu(ievent)%heta(NGLLY,acqui_simu(ievent)%nsta_slice))
          allocate(acqui_simu(ievent)%hgamma(NGLLZ,acqui_simu(ievent)%nsta_slice))
          allocate(acqui_simu(ievent)%hpxi(NGLLX,acqui_simu(ievent)%nsta_slice))
          allocate(acqui_simu(ievent)%hpeta(NGLLY,acqui_simu(ievent)%nsta_slice))
          allocate(acqui_simu(ievent)%hpgamma(NGLLZ,acqui_simu(ievent)%nsta_slice))
          allocate(acqui_simu(ievent)%freqcy_to_invert(NDIM,2,acqui_simu(ievent)%nsta_slice))
       else
          allocate(acqui_simu(ievent)%hxi(1,1))
          allocate(acqui_simu(ievent)%heta(1,1))
          allocate(acqui_simu(ievent)%hgamma(1,1))
          allocate(acqui_simu(ievent)%hpxi(1,1))
          allocate(acqui_simu(ievent)%hpeta(1,1))
          allocate(acqui_simu(ievent)%hpgamma(1,1))
          allocate(acqui_simu(ievent)%freqcy_to_invert(NDIM,2,1))
       endif

       irec_local=0
       do ireceiver=1, acqui_simu(ievent)%nsta_tot
          if (myrank == acqui_simu(ievent)%islice_selected_rec(ireceiver)) then

             irec_local = irec_local + 1

             xi_receiver=acqui_simu(ievent)%xi_rec(ireceiver)
             eta_receiver=acqui_simu(ievent)%eta_rec(ireceiver)
             gamma_receiver=acqui_simu(ievent)%gamma_rec(ireceiver)

             ! compute Lagrange polynomials at the source location
             call lagrange_any(xi_receiver,NGLLX,xigll,hxis,hpxis)
             call lagrange_any(eta_receiver,NGLLY,yigll,hetas,hpetas)
             call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammas,hpgammas)

             acqui_simu(ievent)%hxi(:,irec_local)=hxis(:)
             acqui_simu(ievent)%heta(:,irec_local)=hetas(:)
             acqui_simu(ievent)%hgamma(:,irec_local)=hgammas(:)

             acqui_simu(ievent)%hpxi(:,irec_local)=hpxis(:)
             acqui_simu(ievent)%hpeta(:,irec_local)=hpetas(:)
             acqui_simu(ievent)%hpgamma(:,irec_local)=hpgammas(:)


          endif
       enddo

    enddo

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) ' ... locate stations passed '
       write(INVERSE_LOG_FILE,*)
       call flush_iunit(INVERSE_LOG_FILE)
    endif

  end subroutine setup_teleseismic_stations

  
end module Teleseismic_IO_mod
