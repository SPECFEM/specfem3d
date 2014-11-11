      integer mbfi,mbfd
      parameter( mbfi=20000,mbfd=800000 )
      include  'mpif.h'
      integer   rank,size,status(mpi_status_size),ipos,dpos
      integer   param(4),ibuff(mbfi)
      real*8    dbuff(mbfd)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common / mpi_area / ibuff,dbuff,rank,size,ipos,dpos
