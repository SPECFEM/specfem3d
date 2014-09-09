cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mgetd(dv,ml,in)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'read real*8 from message buffer'.
      include 'mpicom.h'
      integer ml,in
      real*8 dv(*)
c
      if ( mbfd.lt.dpos/8+ml ) stop 'mbfd is too small (mgetd).'
c
      call mpi_unpack( dbuff,mbfd*8,dpos,dv,ml,mpi_double_precision,
     1                 mpi_comm_world,in )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mgeti(iv,ml,in)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'read integer from message buffer'.
      include 'mpicom.h'
      integer iv(*),ml,in
c
      if ( mbfi.lt.ipos/4+ml ) stop 'mbfi is too small (mgeti).'
c
      call mpi_unpack( ibuff,mbfi*4,ipos,iv,ml,mpi_integer,
     1                 mpi_comm_world,in )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine minit
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'initiate processes'.
c
c ----------------------------------------------------------------------
c
      include 'mpicom.h'
c
c ----------------------------------------------------------------------
c
      integer in
c
      call mpi_init( in )
      call mpi_comm_rank( mpi_comm_world,rank,in )
      call mpi_comm_size( mpi_comm_world,size,in )
c
      ipos=0
      dpos=0
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mputd(dv,ml,in)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'write real*8 to message buffer'.
      include 'mpicom.h'
      integer ml,in
      real*8 dv(*)
c
      if ( mbfd.lt.dpos/8+ml ) stop 'mbfd is too small (mputd).'
c
      call mpi_pack( dv,ml,mpi_double_precision,dbuff,mbfd*8,dpos,
     1               mpi_comm_world,in )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mputi(iv,ml,in)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'write integer to message buffer'.
      include 'mpicom.h'
      integer iv(*),ml,in
c
      if ( mbfi.lt.ipos/4+ml ) stop 'mbfi is too small (mputi).'
c
      call mpi_pack( iv,ml,mpi_integer,ibuff,mbfi*4,ipos,
     1               mpi_comm_world,in )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mrecv(mt,in,ierr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'receive message'.
c
c ----------------------------------------------------------------------
c
      include 'mpicom.h'
      integer mt,in
      integer ierr,numi,numd,np
c
c ----------------------------------------------------------------------
c
      if ( in.ge.0 ) then
        np=in
        call mpi_recv( param,4,mpi_integer,np,mt+100000,
     1                 mpi_comm_world,status,ierr )
        ipos=0
        dpos=0
        numi=param(1)/4
        numd=param(2)/8
        if( numi.ne.0 ) then
          call mpi_recv( ibuff,numi,mpi_integer,np,mt+110000,
     1                   mpi_comm_world,status,ierr )
        endif
        if( numd.ne.0 ) then
          call mpi_recv( dbuff,numd,mpi_double_precision,np,mt+120000,
     1                   mpi_comm_world,status,ierr )
        endif
      else
        stop 'Invalid arguments (mrecv).'
      endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine msend(mn,mt,ierr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'send message'.
c
c ----------------------------------------------------------------------
c
      include 'mpicom.h'
      integer mn,mt,ierr
      integer np
c
c ----------------------------------------------------------------------
c
      if( mn.lt.0 ) then
        param(1)=ipos
        param(2)=dpos
        do 100 np = 1, size-1
          call mpi_send( param,4,mpi_integer,np,mt+100000,
     1                   mpi_comm_world,ierr )
          if( ipos.ne.0 ) then
            call mpi_send( ibuff,ipos/4,mpi_integer,np,mt+110000,
     1                     mpi_comm_world,ierr )
          endif
          if( dpos.ne.0 ) then
            call mpi_send( dbuff,dpos/8,mpi_double_precision,np,
     1                     mt+120000,mpi_comm_world,ierr )
          endif
  100   continue
      endif
      if( mn.ge.0 ) then
        param(1)=ipos
        param(2)=dpos
        call mpi_send( param,4,mpi_integer,mn,mt+100000,
     1                 mpi_comm_world,ierr )
        if( ipos.ne.0 ) then
          call mpi_send( ibuff,ipos/4,mpi_integer,mn,mt+110000,
     1                   mpi_comm_world,ierr )
        endif
        if( dpos.ne.0 ) then
          call mpi_send( dbuff,dpos/8,mpi_double_precision,mn,
     1                   mt+120000,mpi_comm_world,ierr )
        endif
      endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine msndi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'initialize send buffer'.
c
c ----------------------------------------------------------------------
c
      include 'mpicom.h'
c
c ----------------------------------------------------------------------
c
      ipos=0
      dpos=0
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mterm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        perform function 'exit communication software'.
      include 'mpicom.h'
      integer in
c
      call mpi_finalize(in)
c
      return
      end
