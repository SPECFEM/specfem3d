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

  subroutine create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC)

! create name of the database for serial codes (AVS_DX and codes to check buffers)

  implicit none

  include "constants.h"

  integer iproc,NPROC

! name of the database file
  character(len=150) prname,procname,LOCAL_PATH,clean_LOCAL_PATH,serial_prefix

  integer iprocloop,nproc_max_loop
  integer, dimension(:), allocatable :: num_active_proc

  nproc_max_loop = NPROC-1

! create the name for the database of the current slide and region
  write(procname,10) iproc
 10 format('/proc',i4.4,'_')

! on a Beowulf-type machine, path on frontend can be different from local paths
  if(USE_LOCAL_PATH_BEOWULF) then

! allocate array for active processors
    allocate(num_active_proc(0:nproc_max_loop))

! read filtered file with name of active machines
    open(unit=48,file='OUTPUT_FILES/filtered_machines.txt',status='old')
    do iprocloop = 0,nproc_max_loop
      read(48,*) num_active_proc(iprocloop)
    enddo
    close(48)

! create the serial prefix pointing to the correct machine
    write(serial_prefix,20) num_active_proc(iproc)
 20 format('/auto/scratch_n',i3.3,'/')

! suppress everything until the last "/" to define the base name of local path
! this is system dependent since it assumes the disks are mounted
! as on our Beowulf (Unix and NFS)
    clean_LOCAL_PATH = LOCAL_PATH(index(LOCAL_PATH,'/',.true.)+1:len_trim(LOCAL_PATH))

! create full name with path
    prname = serial_prefix(1:len_trim(serial_prefix)) // clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // procname

! deallocate array
    deallocate(num_active_proc)

! on shared-memory machines, global path is the same as local path
  else

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full name with path
    prname = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // procname

  endif

  end subroutine create_serial_name_database

