!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  subroutine create_serial_name_database(prname,iproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! create name of the database for serial codes (AVS_DX and codes to check buffers)

  implicit none

  include "constants.h"

  integer iproc,iregion_code,NPROCTOT

! name of the database file
  character(len=150) prname,procname,base_path,serial_prefix, &
      LOCAL_PATH,OUTPUT_FILES

  integer iprocloop
  integer, dimension(:), allocatable :: num_active_proc

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_reg',i1,'_')") iproc,iregion_code

! on a machine with local disks, path on frontend can be different from local paths
  if(.not. LOCAL_PATH_IS_ALSO_GLOBAL) then

! allocate array for active processors
    allocate(num_active_proc(0:NPROCTOT-1))

! read filtered file with name of active machines
    open(unit=48,file=trim(OUTPUT_FILES)//'/filtered_machines.txt',status='old',action='read')
    do iprocloop = 0,NPROCTOT-1
      read(48,*) num_active_proc(iprocloop)
    enddo
    close(48)

! create the serial prefix pointing to the correct machine
    write(serial_prefix,"('/auto/scratch_n',i6.6,'/')") num_active_proc(iproc)

! suppress everything until the last "/" to define the base name of local path
! this is system dependent since it assumes the disks are mounted remotely
    base_path = LOCAL_PATH(index(LOCAL_PATH,'/',.true.)+1:len_trim(LOCAL_PATH))

! create full name with path
    prname = trim(serial_prefix) // trim(base_path) // procname

! deallocate array
    deallocate(num_active_proc)

! on shared-memory machines, global path is the same as local path
  else

! create full name with path
    prname = trim(LOCAL_PATH) // procname

  endif

  end subroutine create_serial_name_database

