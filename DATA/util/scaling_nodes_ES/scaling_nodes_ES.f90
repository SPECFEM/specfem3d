
! study scaling and optimal parameters on the ES

  program scaling_nodes_ES

  implicit none

! ensure 100 % load balancing on the ES (i.e., even values of NPROC_XI only)
  logical, parameter :: ONLY_100PERCENT = .false.

! store arrays in single or double precision on the ES
  logical, parameter :: DOUBLE_PRECISION_ES = .false.

  integer NPROC_XI,procs,rounded_es_nodes,mem_avail
  integer imult,NER_BASE,NER
  integer ifirst,ilast,istep
  real es_nodes,percent_ES,ratio_nodes
  real period,memory,ratio_mem

  integer, parameter :: NUMBER_VALUES_DISPLAYED = 30

  if(DOUBLE_PRECISION_ES) then
    print *
    print *,'assuming double precision used on the ES'
    print *
  else
    print *
    print *,'assuming single precision used on the ES'
    print *
  endif

!!! DK DK UGLY
   stop 'need to modify this code for basins, only one chunk, and NPROC_XI can differ from NPROC_ETA (same for NER_XI)'

! only even numbers use 100 % of the processors
! max value of NPROC_XI on the ES is 29
  if(ONLY_100PERCENT) then
    ifirst = 6
    ilast = 29
    istep = 2
  else
    ifirst = 5
    ilast = 29
    istep = 1
  endif

  do NPROC_XI=ifirst,ilast,istep

    procs = 6*NPROC_XI**2
    es_nodes = procs / 8.
    rounded_es_nodes = int(procs/8. + 0.5)
    percent_ES = procs / 5120.
    mem_avail = rounded_es_nodes * 16
    ratio_nodes = es_nodes / rounded_es_nodes

    print *
    print *,'-----------------------------------------------------------'
    print *
    print *,'NPROC_XI = ',NPROC_XI
    print *,'processors = ',procs
    print *,'ES nodes = ',es_nodes
    print *,'rounded ES nodes = ',rounded_es_nodes
    print *,'percentage of the ES = ',100.*percent_ES,' %'
    print *,'memory available (Gb) = ',mem_avail
    print *,'ratio ES_nodes / rounded_ES_nodes = ',100.*ratio_nodes,' %'

    NER_BASE = 16 * NPROC_XI

    do imult = 1,NUMBER_VALUES_DISPLAYED
      NER = NER_BASE * imult
! on our cluster, reference is NER = 240, accurate to 18 s, uses 52.5 Gb (70 % of 75)
      period = 18. / (NER/240.)
      memory = 52.5 * (NER/240.)**3
      if(DOUBLE_PRECISION_ES) memory = memory * 2.
      ratio_mem = 100.*memory/mem_avail
!! DK DK display interesting periods only
  if(period < 17.) then
      if(ratio_mem >= 91. .and. ratio_mem < 100.) then
        write(*,200) imult,NER
        print *,'period (s) = ',period,'  memory (Gb) = ',memory
        print *,' usage (%) = ',ratio_mem,' #### TOO HIGH #### period ',period
      else if(ratio_mem >= 85. .and. ratio_mem < 100.) then
        write(*,200) imult,NER
        print *,'period (s) = ',period,'  memory (Gb) = ',memory
        print *,' usage (%) = ',ratio_mem,' **** excellent **** period ',period
      else if(ratio_mem >= 75. .and. ratio_mem < 100.) then
        write(*,200) imult,NER
        print *,'period (s) = ',period,'  memory (Gb) = ',memory
        print *,' usage (%) = ',ratio_mem,' **** very good **** period ',period
      else if(ratio_mem >= 65. .and. ratio_mem < 100.) then
        write(*,200) imult,NER
        print *,'period (s) = ',period,'  memory (Gb) = ',memory
        print *,' usage (%) = ',ratio_mem,' **** correct **** period ',period
      else if(ratio_mem >= 50. .and. ratio_mem < 100.) then
        write(*,200) imult,NER
        print *,'period (s) = ',period,'  memory (Gb) = ',memory
        print *,' usage (%) = ',ratio_mem,' **** fair **** period ',period
      else if(ratio_mem >= 100.) then
!! DK DK do not display if memory usage is above 100 %
!! DK DK suppressed        write(*,200) imult,NER
!! DK DK suppressed        print *,'period (s) = ',period,'  memory (Gb) = ',memory
!! DK DK suppressed        print *,' usage (%) = ',ratio_mem,' ######## impossible'
      else
!! DK DK do not display if memory usage is too low
!        write(*,200) imult,NER
!        print *,'period (s) = ',period,'  memory (Gb) = ',memory
!        print *,' usage (%) = ',ratio_mem
      endif
  endif
    enddo

 200 format('possible NER_',i2.2,' = ',i6)

  enddo

  print *
  print *

  end program scaling_nodes_ES

