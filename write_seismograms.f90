!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! write seismograms to text files

  subroutine write_seismograms(myrank,seismograms,number_receiver_global, &
               station_name,network_name,nrec,nrec_local, &
               it,DT,NSTEP,hdur,LOCAL_PATH)

  implicit none

  include "constants.h"

  integer nrec,nrec_local,NSTEP,it,myrank
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(3,nrec_local,NSTEP) :: seismograms
  character(len=8), dimension(nrec) :: station_name,network_name
  double precision hdur,DT
  character(len=150) LOCAL_PATH

  integer irec,irec_local,length_station_name,length_network_name
  integer iorientation,irecord,isample

  character(len=4) chn
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

! save three components of displacement vector
    irecord = 1

    do iorientation = 1,3

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
      if(length_station_name < 1 .or. length_station_name > 8 .or. &
         length_network_name < 1 .or. length_network_name > 2) &
           call exit_MPI(myrank,'wrong length of station or network name')

    if(length_network_name == 1) then

      if(length_station_name == 1) then
        write(sisname,"(a1,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 2) then
        write(sisname,"(a2,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 3) then
        write(sisname,"(a3,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 4) then
        write(sisname,"(a4,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 5) then
        write(sisname,"(a5,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 6) then
        write(sisname,"(a6,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 7) then
        write(sisname,"(a7,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else
        write(sisname,"(a8,'.',a1,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      endif

    else

      if(length_station_name == 1) then
        write(sisname,"(a1,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 2) then
        write(sisname,"(a2,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 3) then
        write(sisname,"(a3,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 4) then
        write(sisname,"(a4,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 5) then
        write(sisname,"(a5,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 6) then
        write(sisname,"(a6,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else if(length_station_name == 7) then
        write(sisname,"(a7,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      else
        write(sisname,"(a8,'.',a2,'.',a3,'.semd')") station_name(irec),network_name(irec),chn
      endif

    endif

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname,status='unknown')

! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
        if(irecord == 1) then
          write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(iorientation,irec_local,isample)
        else
          call exit_MPI(myrank,'incorrect record label')
        endif
      enddo

      close(IOUT)

      enddo

  enddo

  end subroutine write_seismograms

