
! DK DK determine if we compile in serial or in MPI mode
#ifndef USE_SERIAL
#define USE_MPI
#endif

! DK DK if we compile in serial mode, we do not need any routine from this file
#ifdef USE_MPI

  subroutine mgetd(dv,ml,in,mbfd,dpos,dbuff)

! perform function 'read real(kind=8) from message buffer'.

  use mpi

  implicit none

  integer :: mbfd
  integer :: dpos
  real(kind=8), dimension(mbfd) :: dbuff

  integer ml,in
  real(kind=8) dv(*)

  if ( mbfd < dpos/8+ml ) stop 'mbfd is too small (mgetd).'

  call mpi_unpack( dbuff,mbfd*8,dpos,dv,ml,mpi_double_precision,mpi_comm_world,in )

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mgeti(iv,ml,in,mbfi,ipos,ibuff)

! perform function 'read integer from message buffer'.

  use mpi

  implicit none

  integer :: mbfi
  integer :: ipos
  integer, dimension(mbfi) :: ibuff

  integer iv(*),ml,in

  if ( mbfi < ipos/4+ml ) stop 'mbfi is too small (mgeti).'

  call mpi_unpack( ibuff,mbfi*4,ipos,iv,ml,mpi_integer,mpi_comm_world,in )

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mputd(dv,ml,in,mbfd,dpos,dbuff)

! perform function 'write real(kind=8) to message buffer'.

  use mpi

  implicit none

  integer :: mbfd
  integer :: dpos
  real(kind=8), dimension(mbfd) :: dbuff

  integer ml,in
  real(kind=8) dv(*)

  if ( mbfd < dpos/8+ml ) stop 'mbfd is too small (mputd).'

  call mpi_pack( dv,ml,mpi_double_precision,dbuff,mbfd*8,dpos,mpi_comm_world,in )

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mputi(iv,ml,in,mbfi,ipos,ibuff)

! perform function 'write integer to message buffer'.

  use mpi

  implicit none

  integer :: mbfi
  integer :: ipos
  integer, dimension(mbfi) :: ibuff

  integer iv(*),ml,in

  if ( mbfi < ipos/4+ml ) stop 'mbfi is too small (mputi).'

  call mpi_pack( iv,ml,mpi_integer,ibuff,mbfi*4,ipos,mpi_comm_world,in )

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mrecv(mt,in,ierr,mbfi,mbfd,ipos,dpos,param,ibuff,dbuff)

! perform function 'receive message'.

  use mpi

  implicit none

  integer :: mbfi,mbfd
  integer :: ipos,dpos
  integer, dimension(4) :: param
  integer, dimension(mbfi) :: ibuff
  real(kind=8), dimension(mbfd) :: dbuff

  integer mt,in
  integer ierr,numi,numd,np

  integer, dimension(mpi_status_size) :: status

  if ( in >= 0 ) then
    np=in
!! DK DK +100000 here is just an optional message tag
    call mpi_recv( param,4,mpi_integer,np,mt+100000,mpi_comm_world,status,ierr )
    ipos = 0
    dpos = 0
    numi = param(1) / 4
    numd = param(2) / 8
    if( numi /= 0 ) call mpi_recv( ibuff,numi,mpi_integer,np,mt+110000,mpi_comm_world,status,ierr )
    if( numd /= 0 ) call mpi_recv( dbuff,numd,mpi_double_precision,np,mt+120000,mpi_comm_world,status,ierr )
  else
    stop 'Invalid arguments (mrecv).'
  endif

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine msend(mn,mt,ierr,mbfi,mbfd,size,ipos,dpos,param,ibuff,dbuff)

! perform function 'send message'.

  use mpi

  implicit none

  integer :: mbfi,mbfd
  integer :: size,ipos,dpos
  integer, dimension(4) :: param
  integer, dimension(mbfi) :: ibuff
  real(kind=8), dimension(mbfd) :: dbuff

  integer mn,mt,ierr
  integer np

  if( mn < 0 ) then
    param(1) = ipos
    param(2) = dpos

    do np = 1, size-1
!! DK DK +100000 here is just an optional message tag
      call mpi_send( param,4,mpi_integer,np,mt+100000, mpi_comm_world,ierr )
      if( ipos /= 0 ) call mpi_send( ibuff,ipos/4,mpi_integer,np,mt+110000, mpi_comm_world,ierr )
      if( dpos /= 0 ) call mpi_send( dbuff,dpos/8,mpi_double_precision,np, mt+120000,mpi_comm_world,ierr )
    enddo
  endif

  if( mn >= 0 ) then
    param(1) = ipos
    param(2) = dpos
!! DK DK +100000 here is just an optional message tag
    call mpi_send( param,4,mpi_integer,mn,mt+100000, mpi_comm_world,ierr )
    if( ipos /= 0 ) call mpi_send( ibuff,ipos/4,mpi_integer,mn,mt+110000, mpi_comm_world,ierr )
    if( dpos /= 0 ) call mpi_send( dbuff,dpos/8,mpi_double_precision,mn, mt+120000,mpi_comm_world,ierr )
  endif

  end

#endif

