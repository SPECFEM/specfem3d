
  program jfdkfd

!! DK DK regrid Hauksson on regular grid in lat-long for So-Cal

  implicit none

  include "../../constants.h"

!! DK DK size of smoothing window
  integer, parameter :: NSIZE = 2
  double precision, dimension(16,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: value_new,value_old

  integer i,ival,j,ivalue,iloop,jloop,jval

  double precision meanval

! read raw grid in UTM space
  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
    read(*,*) (value_old(ivalue,i,j),ivalue=1,16)
  enddo
  enddo

! smooth with window
  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
    do ivalue = 1,16
      meanval = 0.d0
      do iloop = i-NSIZE,i+NSIZE
      do jloop = j-NSIZE,j+NSIZE
        ival = iloop
        jval = jloop
        if(ival<1) ival = 1
        if(jval<1) jval = 1
        if(ival>NGRID_NEW_HAUKSSON) ival = NGRID_NEW_HAUKSSON
        if(jval>NGRID_NEW_HAUKSSON) jval = NGRID_NEW_HAUKSSON
        meanval = meanval + value_old(ivalue,ival,jval)
      enddo
      enddo
      meanval = meanval / dble((2*NSIZE+1)**2)
      value_new(ivalue,i,j) = meanval
    enddo
    !print *,(value_new(ivalue,i,j),ivalue=1,16)
    write(*,'(3f7.4)') (value_new(ivalue,i,j),ivalue=1,16)
  enddo
  enddo

  end

