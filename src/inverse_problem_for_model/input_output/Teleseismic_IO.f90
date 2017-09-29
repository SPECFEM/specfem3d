module Teleseismic_IO_mod

  ! from specfem 
  use constants, only: mygroup
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NGNOD

  ! from inversion 
  use inverse_problem_par
  use mesh_tools

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
    character(len=MAX_LEN_STRING)                              :: line, keyw !, line_to_read
    integer                                                    :: ipos0, ipos1, ievent
    integer                                                    :: ier
    
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
       end do
99     close(666)
       
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
       do 
          read(666,'(a)',end=999) line
          !if (is_blank_line(line)) cycle
          ! INDICES TO READ line -----------------------------------------------
          ipos0=index(line,':')+1
          ipos1=index(line,'#')-1
          if (ipos1 < 0 ) ipos1=len_trim(line)
          
          !! STORE KEYWORD ITEM -------------------------------------------------
          keyw=trim(adjustl(line(1:ipos0-2)))
          select case (trim(keyw))
             
          case('event_name')
             ievent=ievent+1
             !! this name is the name of file contains all information on teleseismic event
             acqui_simu(ievent)%event_name=trim(adjustl(line(ipos0:ipos1)))
             acqui_simu(ievent)%nevent_tot=NEVENT

             !! to do define this routine  .... for reading all info for event 
             call read_one_teleseismic_event()
            
          end select
       end do

999    close(666)

    end if

    
    ! master broadcasts read values
    call MPI_BCAST(NEVENT,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    if (myrank > 0) allocate(acqui_simu(NEVENT))
    do ievent = 1, NEVENT
       call MPI_BCAST(acqui_simu(ievent)%event_name,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
    end do

    !! to do do not forget to bcast all acqui_simu structure to other mpi slices  
    !! ......
    
    

  end subroutine read_acqui_teleseismic_file


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------------
!> subroutine to read one event file 
!-------------------------------------------------------------------------------------------------------------------
  subroutine read_one_teleseismic_event()

    !! reading input event teleseismic file
    
    !!
    !!
    write(*,*) ' teleseismic event is under construction then we cannnot provide it for now ....'
    !! need to compute cartesian coordinate of stations : acqui_simu(ievent)%position_station(3,nsta_tot)
    !!
    !!

  end subroutine read_one_teleseismic_event

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------------
!> store arrays that needed for specfem to be able to use stations in mesh 
!-------------------------------------------------------------------------------------------------------------------
  subroutine setup_teleseismic_stations(acqui_simu, myrank)
    
    !! from  acqui_simu%position_station assumed to be in cartesian coordinate
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
