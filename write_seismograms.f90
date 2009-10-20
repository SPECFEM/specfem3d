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

! write seismograms to text files

  subroutine write_seismograms(myrank,seismograms,number_receiver_global, &
               station_name,network_name,nrec,nrec_local, &
               it,DT,NSTEP,hdur,LOCAL_PATH,istore)

  implicit none

  include "constants.h"

  integer nrec,nrec_local,NSTEP,it,myrank,istore
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,irec_local,length_station_name,length_network_name
  integer iorientation,irecord,isample

  character(len=4) chn
  character(len=1) component
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

! parameters for master collects seismograms  
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram
  integer :: nrec_local_received,NPROCTOT,total_seismos,receiver,sender
  integer :: iproc,ier
   
! save displacement, velocity or acceleration
  if(istore == 1) then
    component = 'd'
  else if(istore == 2) then
    component = 'v'
  else if(istore == 3) then
    component = 'a'
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

! all the processes write their local seismograms themselves
  if( .not. WRITE_SEISMOGRAMS_BY_MASTER ) then

    do irec_local = 1,nrec_local

! get global number of that receiver
      irec = number_receiver_global(irec_local)

! save three components of displacement vector
      irecord = 1

      do iorientation = 1,NDIM

        if(iorientation == 1) then
          chn = 'BHE'
        else if(iorientation == 2) then
          chn = 'BHN'
        else if(iorientation == 3) then
          chn = 'BHZ'
        else
          call exit_MPI(myrank,'incorrect channel value')
        endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component
        length_station_name = len_trim(station_name(irec))
        length_network_name = len_trim(network_name(irec))

! check that length conforms to standard
        if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
           call exit_MPI(myrank,'wrong length of station name')

        if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
           call exit_MPI(myrank,'wrong length of network name')

        write(sisname,"(a,'.',a,'.',a3,'.sem',a1)") station_name(irec)(1:length_station_name),&
           network_name(irec)(1:length_network_name),chn,component

! directory to store seismograms
        if( USE_OUTPUT_FILES_PATH ) then      
          final_LOCAL_PATH = 'OUTPUT_FILES'//'/'        
        else      
! suppress white spaces if any
          clean_LOCAL_PATH = adjustl(LOCAL_PATH)
! create full final local path
          final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'        
        endif
      
            
! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')

! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
        do isample = 1,min(it,NSTEP)
          if(irecord == 1) then
! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(iorientation,irec_local,isample)
            else
              write(IOUT,*) dble(isample-1)*DT - hdur,' ',seismograms(iorientation,irec_local,isample)
            endif
          else
            call exit_MPI(myrank,'incorrect record label')
          endif
        enddo

        close(IOUT)

      enddo ! NDIM

    enddo ! nrec_local

! now only the master process does the writing of seismograms and
! collects the data from all other processes
  else ! WRITE_SEISMOGRAMS_BY_MASTER

    allocate(one_seismogram(NDIM,NSTEP),stat=ier)
    if(ier /= 0) stop 'error while allocating one temporary seismogram'

  
    if(myrank == 0) then ! on the master, gather all the seismograms

      total_seismos = 0

      ! loop on all the slices
      call world_size(NPROCTOT)      
      do iproc = 0,NPROCTOT-1

        ! receive except from proc 0, which is me and therefore I already have this value
        sender = iproc
        if(iproc /= 0) then
          call recv_i(nrec_local_received,1,sender,itag)
          if(nrec_local_received < 0) call exit_MPI(myrank,'error while receiving local number of receivers')
        else
          nrec_local_received = nrec_local
        endif
         
        if (nrec_local_received > 0) then
          do irec_local = 1,nrec_local_received
            ! receive except from proc 0, which is myself and therefore I already have these values
            if(iproc == 0) then
              ! get global number of that receiver
              irec = number_receiver_global(irec_local)
              one_seismogram(:,:) = seismograms(:,irec_local,:)
            else
              call recv_i(irec,1,sender,itag)
              if(irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')
              
              call recvv_cr(one_seismogram,NDIM*NSTEP,sender,itag)
            endif

            total_seismos = total_seismos + 1

! save three components of displacement vector
            irecord = 1

            do iorientation = 1,NDIM

              if(iorientation == 1) then
                chn = 'BHE'
              else if(iorientation == 2) then
                chn = 'BHN'
              else if(iorientation == 3) then
                chn = 'BHZ'
              else
                call exit_MPI(myrank,'incorrect channel value')
              endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component
              length_station_name = len_trim(station_name(irec))
              length_network_name = len_trim(network_name(irec))

! check that length conforms to standard
              if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
                call exit_MPI(myrank,'wrong length of station name')

              if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
                call exit_MPI(myrank,'wrong length of network name')

              write(sisname,"(a,'.',a,'.',a3,'.sem',a1)") station_name(irec)(1:length_station_name),&
                network_name(irec)(1:length_network_name),chn,component

! directory to store seismograms
              if( USE_OUTPUT_FILES_PATH ) then      
                final_LOCAL_PATH = 'OUTPUT_FILES'//'/'        
              else      
! suppress white spaces if any
                clean_LOCAL_PATH = adjustl(LOCAL_PATH)
! create full final local path
                final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'        
              endif

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
              open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')

! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
              do isample = 1,min(it,NSTEP)
                if(irecord == 1) then
! distinguish between single and double precision for reals
                  if(CUSTOM_REAL == SIZE_REAL) then
                    write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',one_seismogram(iorientation,isample)
                  else
                    write(IOUT,*) dble(isample-1)*DT - hdur,' ',one_seismogram(iorientation,isample)
                  endif
                else
                  call exit_MPI(myrank,'incorrect record label')
                endif
              enddo

              close(IOUT)

            enddo ! NDIM
          enddo ! nrec_local_received
        endif ! if(nrec_local_received > 0 )
      enddo ! NPROCTOT-1

      write(IMAIN,*) 'Component: .sem'//component
      write(IMAIN,*) '  total number of receivers saved is ',total_seismos,' out of ',nrec
      write(IMAIN,*)

      if(total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

    else  ! on the nodes, send the seismograms to the master
       receiver = 0
       call send_i(nrec_local,1,receiver,itag)
       if (nrec_local > 0) then
         do irec_local = 1,nrec_local
           ! get global number of that receiver
           irec = number_receiver_global(irec_local)
           call send_i(irec,1,receiver,itag)
           
           ! sends seismogram of that receiver
           one_seismogram(:,:) = seismograms(:,irec_local,:)
           call sendv_cr(one_seismogram,NDIM*NSTEP,receiver,itag)
         enddo
       endif
    endif ! myrank
  
    deallocate(one_seismogram)
    
  endif ! WRITE_SEISMOGRAMS_BY_MASTER

  end subroutine write_seismograms

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms(myrank,seismograms,number_receiver_global, &
               nrec_local,it,DT,NSTEP,hdur,LOCAL_PATH,istore)

  implicit none

  include "constants.h"

  integer nrec_local,NSTEP,it,myrank,istore
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH


  integer irec,irec_local
  integer iorientation,irecord,isample

  character(len=4) chn
  character(len=1) component
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

! save displacement, velocity or acceleration
  if(istore == 1) then
    component = 'd'
  else if(istore == 2) then
    component = 'v'
  else if(istore == 3) then
    component = 'a'
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

! save three components of displacement vector
    irecord = 1

    do iorientation = 1,NDIM

      if(iorientation == 1) then
        chn = 'BHE'
      else if(iorientation == 2) then
        chn = 'BHN'
      else if(iorientation == 3) then
        chn = 'BHZ'
      else
        call exit_MPI(myrank,'incorrect channel value')
      endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component

      write(sisname,"(a,i5.5,'.',a,'.',a3,'.sem',a1)") 'S',irec_local,&
           'NT',chn,component

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')

! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
        if(irecord == 1) then
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(iorientation,irec_local,isample)
          else
            write(IOUT,*) dble(isample-1)*DT - hdur,' ',seismograms(iorientation,irec_local,isample)
          endif
        else
          call exit_MPI(myrank,'incorrect record label')
        endif
      enddo

      close(IOUT)

      enddo

  enddo

  end subroutine write_adj_seismograms

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2(myrank,seismograms,number_receiver_global, &
               nrec_local,it,DT,NSTEP,hdur,LOCAL_PATH)

  implicit none

  include "constants.h"

  integer nrec_local,NSTEP,it,myrank
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH


  integer irec,irec_local
  integer idim,jdim,irecord,isample

  character(len=4) chn
  character(len=1) component
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

! save three components of displacement vector
    irecord = 1

    do idim = 1, 3
      do jdim = idim, 3

      if(idim == 1 .and. jdim == 1) then
        chn = 'SNN'
      else if(idim == 1 .and. jdim == 2) then
        chn = 'SEN'
      else if(idim == 1 .and. jdim == 3) then
        chn = 'SEZ'
      else if(idim == 2 .and. jdim == 2) then
        chn = 'SEE'
      else if(idim == 2 .and. jdim == 3) then
        chn = 'SNZ'
      else if(idim == 3 .and. jdim == 3) then
        chn = 'SZZ'
      else
        call exit_MPI(myrank,'incorrect channel value')
      endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component

      write(sisname,"(a,i5.5,'.',a,'.',a3,'.sem',a1)") 'S',irec_local,&
           'NT',chn,component

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')

! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
        if(irecord == 1) then
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(jdim,idim,irec_local,isample)
          else
            write(IOUT,*) dble(isample-1)*DT - hdur,' ',seismograms(jdim,idim,irec_local,isample)
          endif
        else
          call exit_MPI(myrank,'incorrect record label')
        endif
      enddo

      close(IOUT)

    enddo

  enddo

end do

end subroutine write_adj_seismograms2
