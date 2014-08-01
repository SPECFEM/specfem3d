!
! $Id:$
!
!----------------------------------------------------------------------
! Code that uses the libtau library to output arrival times for phases
!----------------------------------------------------------------------

  subroutine ttimes(dist_deg,depth,nphases,names,times)

  implicit none
  integer, parameter :: MAX_PHASES=60

  real, intent(in) :: dist_deg, depth

  integer, intent(out) :: nphases 
  character*8, dimension(*), intent(out) :: names
  double precision, dimension(*), intent(out) :: times

  ! legacy variables needed to use the libtau routines
  logical prnt(3)
  character*8 phlst(10)
  real usrc(2)
  real, dimension(MAX_PHASES) :: dtdd,dtdh,dddp
  real, dimension(MAX_PHASES) :: times_sngl
  character*262 modnam
  character*256 iaspmod


  ! ask for all phases
  phlst(1)="all"
  prnt(1)=.false.
  prnt(2)=.false.
  prnt(3)=.false.
  call getenv('IASPMODEL', iaspmod)
  if (trim(iaspmod) == '') then
    modnam='iasp91'
  else
    modnam=iaspmod
  endif
  call tabin(1,modnam)
  call brnset(1,phlst,prnt)
  call depset(depth,usrc)

  call trtm(dist_deg,MAX_PHASES,nphases,times_sngl,dtdd,dtdh,dddp,names)
  times(1:nphases) = dble(times_sngl(1:nphases))
  
  !write(*,*) 'Found ', nphases, ' phases:'

  end subroutine ttimes
