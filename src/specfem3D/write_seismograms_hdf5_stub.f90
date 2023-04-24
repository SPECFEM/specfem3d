!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine write_seismograms_h5()

  end subroutine write_seismograms_h5


!================================================================

! write seismograms to text files

  subroutine write_seismograms_to_file_h5(seismograms,istore)

  use constants
  use specfem_par
  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms


  end subroutine write_seismograms_to_file_h5

!=====================================================================

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms_to_file_h5(seismograms,istore)

  use constants
  use specfem_par

  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismograms

  end subroutine write_adj_seismograms_to_file_h5

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2_to_file_h5(myrank,seismograms_eps,number_receiver_global,nrec_local,it,DT,NSTEP,t0)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  implicit none
  integer :: myrank
  integer :: nrec_local,NSTEP,it
  integer, dimension(nrec_local) :: number_receiver_global
  ! note: seismograms here is still an array of size *,*,*,NSTEP
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local,NSTEP) :: seismograms_eps
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
          write(IOUT,*) real(dble(isample-1)*DT - t0,kind=CUSTOM_REAL),' ',seismograms_eps(jdimval,idimval,irec_local,isample)
        enddo

        close(IOUT)

      enddo ! jdimval
    enddo ! idimval
  enddo ! irec_local

  end subroutine write_adj_seismograms2_to_file_h5

