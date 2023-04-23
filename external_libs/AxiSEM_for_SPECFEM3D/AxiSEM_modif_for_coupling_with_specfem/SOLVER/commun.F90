!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
!> This is the communication module which loads/sorts data
!! to exchange/examine over the processors.
!!
!! APRIL 2007:
!! At this level, we only call parallel routines, but do not invoke MPI here.
module commun

  use global_parameters
  use commpi
  use data_proc

  implicit none
  !public :: pdistsum_solid_4
  public :: pdistsum_solid, pdistsum_solid_1D, pdistsum_fluid

  public :: assembmass_sum_solid, assembmass_sum_fluid ! assemble and sum massmat
  public :: broadcast_int, broadcast_int_arr
  public :: broadcast_dble, broadcast_char, broadcast_log
  public :: pinit, pend
  public :: pmin, pmax, pmax_int, psum, psum_int, psum_dble
  public :: barrier, pcheck
  public :: comm_elem_number
  private

contains

!-----------------------------------------------------------------------------------------
!> Wrapper routine to use the standard solid communication routine also with
!! arrays with one rank less, i.e. skalar fields
subroutine pdistsum_solid_1D(vec)
  ! @TODO: maybe one could just use a reshape in the calling routines instead
  use data_mesh, only: npol, nel_solid
  real(kind=realkind), intent(inout) :: vec(0:npol,0:npol,1:nel_solid,1:1)

  call pdistsum_solid(vec)
end subroutine pdistsum_solid_1D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This is a driver routine to perform the assembly of field f of dimension nc
!! defined in the solid. The assembly/direct stiffness summation is composed of
!! the "gather" and "scatter" operations, i.e. to add up all element-edge
!! contributions at the global stage and place them back into local arrays.
!! Nissen-Meyer et al., GJI 2007, "A 2-D spectral-element method...", section 4.
!! If nproc>1, the asynchronous messaging scheme is invoked to additionally
!! sum & exchange values on processor boundary points.
subroutine pdistsum_solid(vec, phase)

  use data_mesh, only: gvec_solid, npol, nel_solid, igloc_solid
  use data_time, only: idmpi, iclockmpi
  use clocks_mod

  real(kind=realkind), intent(inout) :: vec(0:,0:,:,:)
  integer, optional, intent(in)      :: phase
  ! phase == 0 => do all (default)
  ! phase == 1 => feed buffer start communication, gather/scatter
  ! phase == 2 => wait for communication to finish, extract buffer
  integer                            :: phase_loc
  integer                            :: nc
  integer                            :: iel, jpol, ipol, idest, ipt

  nc = size(vec, dim=4)
  phase_loc = 0 ! Default value
  if (present(phase)) phase_loc = phase

  if (phase_loc == 0 .or. phase_loc == 1) then

#ifndef serial
     if (nproc > 1) then
        iclockmpi = tick()
        call feed_buffer_solid(vec, nc)
        call send_recv_buffers_solid(nc)
        iclockmpi = tick(id=idmpi, since=iclockmpi)
     endif
#endif

     gvec_solid = 0
     ! Gather element boundaries
     ipt = 1
     do iel = 1, nel_solid
        jpol = 0
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           gvec_solid(idest,1:nc) = gvec_solid(idest,1:nc) + vec(ipol,jpol,iel,1:nc)
           ipt = ipt + 1
        enddo

        do jpol = 1, npol-1
           ipol = 0
           idest = igloc_solid(ipt)
           gvec_solid(idest,1:nc) = gvec_solid(idest,1:nc) + vec(ipol,jpol,iel,1:nc)
           ipt = ipt + npol

           ipol = npol
           idest = igloc_solid(ipt)
           gvec_solid(idest,1:nc) = gvec_solid(idest,1:nc) + vec(ipol,jpol,iel,1:nc)
           ipt = ipt + 1
        enddo

        jpol = npol
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           gvec_solid(idest,1:nc) = gvec_solid(idest,1:nc) + vec(ipol,jpol,iel,1:nc)
           ipt = ipt + 1
        enddo
     enddo

     ! Scatter
     ipt = 1
     do iel = 1, nel_solid
        jpol = 0
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,1:nc) = gvec_solid(idest,1:nc)
           ipt = ipt + 1
        enddo

        do jpol = 1, npol-1
           ipol = 0
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,1:nc) = gvec_solid(idest,1:nc)
           ipt = ipt + npol

           ipol = npol
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,1:nc) = gvec_solid(idest,1:nc)
           ipt = ipt + 1
        enddo

        jpol = npol
        do ipol = 0, npol
           idest = igloc_solid(ipt)
           vec(ipol,jpol,iel,1:nc) = gvec_solid(idest,1:nc)
           ipt = ipt + 1
        enddo
     enddo
  endif

#ifndef serial
  if (nproc > 1) then
     if (phase_loc == 0 .or. phase_loc == 2) then
        iclockmpi = tick()
        call extract_from_buffer_solid(vec,nc)
        iclockmpi = tick(id=idmpi, since=iclockmpi)
     endif
  endif ! nproc > 1
#endif

end subroutine pdistsum_solid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This is a driver routine to perform the assembly of field f of dimension nc
!! defined in the fluid. The assembly/direct stiffness summation is composed of
!! the "gather" and "scatter" operations, i.e. to add up all element-edge
!! contributions at the global stage and place them back into local arrays.
!! Nissen-Meyer et al., GJI 2007, "A 2-D spectral-element method...", section 4.
subroutine pdistsum_fluid(vec, phase)

  use data_mesh, only: igloc_fluid
  use data_time, only: idmpi, iclockmpi
  use clocks_mod
  use data_mesh, only: gvec_fluid, npol, nel_fluid

  real(kind=realkind), intent(inout) :: vec(0:npol,0:npol,nel_fluid)
  integer, optional, intent(in)      :: phase
  ! phase == 0 => do all (default)
  ! phase == 1 => feed buffer start communication, gather/scatter
  ! phase == 2 => wait for communication to finish, extract buffer
  integer                            :: phase_loc
  integer                            :: iel, jpol, ipol, idest, ipt

  phase_loc = 0 ! Default value
  if (present(phase)) phase_loc = phase

  if (phase_loc == 0 .or. phase_loc == 1) then

#ifndef serial
     if (nproc > 1) then
        iclockmpi = tick()
        call feed_buffer_fluid(vec)
        call send_recv_buffers_fluid()
        iclockmpi = tick(id=idmpi, since=iclockmpi)
     endif
#endif

     gvec_fluid(:) = 0.d0

     ! Gather
     ipt = 1
     do iel = 1, nel_fluid
        jpol = 0
        do ipol = 0, npol
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
           ipt = ipt + 1
        enddo

        do jpol = 1, npol-1
           ipol = 0
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
           ipt = ipt + npol

           ipol = npol
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
           ipt = ipt + 1
        enddo

        jpol = npol
        do ipol = 0, npol
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = gvec_fluid(idest) + vec(ipol,jpol,iel)
           ipt = ipt + 1
        enddo
     enddo

     ! Scatter
     ipt = 1
     do iel = 1, nel_fluid
        jpol = 0
        do ipol = 0, npol
           idest = igloc_fluid(ipt)
           vec(ipol,jpol,iel) = gvec_fluid(idest)
           ipt = ipt + 1
        enddo

        do jpol = 1, npol-1
           ipol = 0
           idest = igloc_fluid(ipt)
           vec(ipol,jpol,iel) = gvec_fluid(idest)
           ipt = ipt + npol

           ipol = npol
           idest = igloc_fluid(ipt)
           vec(ipol,jpol,iel) = gvec_fluid(idest)
           ipt = ipt + 1
        enddo

        jpol = npol
        do ipol = 0, npol
           idest = igloc_fluid(ipt)
           vec(ipol,jpol,iel) = gvec_fluid(idest)
           ipt = ipt + 1
        enddo
     enddo
  endif

#ifndef serial
  if (nproc > 1) then
     if (phase_loc == 0 .or. phase_loc == 2) then
        iclockmpi = tick()
        call extract_from_buffer_fluid(vec)
        iclockmpi = tick(id=idmpi, since=iclockmpi)
     endif
  endif
#endif


end subroutine pdistsum_fluid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine assembmass_sum_solid(f1,res)

  use data_mesh, only: igloc_solid
  use data_mesh, only: gvec_solid, npol, nel_solid

  real(kind=realkind), intent(in)   :: f1(0:,0:,:)
  real(kind=dp)   , intent(out)     :: res
  integer                           :: ipt, idest, iel, ipol, jpol

  !!@TODO Optimise for npol = 4
  !! MvD: I guess that is not necessary, because it is not called int the time loop
  res = 0.d0
  gvec_solid(:,1) = 0.d0
  do iel = 1, nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
           idest = igloc_solid(ipt)
           gvec_solid(idest,1) = gvec_solid(idest,1) + f1(ipol,jpol,iel)
        enddo
     enddo
  enddo
  res = res + sum(gvec_solid(:,1))
#ifndef serial
  if (nproc > 1) res=ppsum_dble(res)
#endif

end subroutine assembmass_sum_solid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine assembmass_sum_fluid(f1,res)

  use data_mesh, only: igloc_fluid
  use data_mesh, only: gvec_fluid
  use data_mesh, only: gvec_solid, npol, nel_fluid

  real(kind=realkind), intent(in)   :: f1(0:,0:,:)
  real(kind=dp)   , intent(out)     :: res
  integer ipt, idest
  integer iel, ipol, jpol

  res = 0.d0

  gvec_fluid(:) = 0.d0
  do iel = 1, nel_fluid
     do ipol = 0, npol
        do jpol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
           idest = igloc_fluid(ipt)
           gvec_fluid(idest) = gvec_fluid(idest) + f1(ipol,jpol,iel)
        enddo
     enddo
  enddo
  res = res + sum(gvec_fluid)

#ifndef serial
  if (nproc > 1) res = ppsum_dble(res)
#endif

end subroutine assembmass_sum_fluid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pinit

  use data_io, only: define_io_appendix
  integer           :: ioerr, nproc_mesh
  character(len=20) :: dbname

  ! Get mesh number of processors

  dbname = 'Mesh/meshdb.dat0000'

  open(1000, file=trim(dbname), FORM="UNFORMATTED", &
             STATUS="OLD", POSITION="REWIND", IOSTAT=ioerr)
  if (ioerr /= 0) then
     write(*,*) 'Could not open mesh file ', trim(dbname)
     stop
  endif
  read(1000) nproc_mesh
  close(1000)

#ifndef serial
  ! Start message passing interface if parallel simulation
  if (nproc_mesh > 1) then
     call ppinit
     if (nproc_mesh /= nproc) then
        write(*,*) mynum, 'Problem with number of processors!'
        write(*,*) mynum, 'Mesh constructed for:', nproc_mesh
        write(*,*) mynum, 'Job submission for:', nproc
        stop
     endif
  else
     nproc = nproc_mesh
     mynum = 0
  endif
#endif

#ifdef serial
  if (nproc_mesh /= 1) then
        write(*,*) 'ERROR: Solver compiled with SERIAL flag, but mesh has nproc > 1: ', &
                    nproc_mesh
        stop
  endif

  nproc = 1
  mynum = 0
#endif

  lpr = .false.
  if (nproc > 1) then
     if (mynum == nproc/2-1) lpr = .true.
  else
     lpr = .true.
  endif

  call define_io_appendix(appmynum, mynum)

  procstrg = 'Proc '//appmynum(3:4)//' '

  if (lpr) write(*,'(a,i5)') '    Initialized run for nproc =', nproc

end subroutine pinit
!! End message passing interface if parallel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pend

#ifndef serial
  if (nproc > 1) call ppend
#endif

end subroutine pend
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine broadcast_char(input_char,input_proc)

  character(*), intent(inout)   :: input_char
  integer, intent(in)           :: input_proc

#ifndef serial
  if (nproc > 1) call pbroadcast_char(input_char,input_proc)
#endif

end subroutine broadcast_char
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine broadcast_log(input_log,input_proc)

  integer, intent(in)    :: input_proc
  logical, intent(inout) :: input_log

#ifndef serial
  if (nproc > 1) call pbroadcast_log(input_log,input_proc)
#endif

end subroutine broadcast_log
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine broadcast_int(input_int,input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int

#ifndef serial
  if (nproc > 1) call pbroadcast_int(input_int, input_proc)
#endif

end subroutine broadcast_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine broadcast_int_arr(input_int,input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int(:)

#ifndef serial
  if (nproc > 1) call pbroadcast_int_arr(input_int, input_proc)
#endif

end subroutine broadcast_int_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine broadcast_dble(input_dble,input_proc)

  integer, intent(in)             :: input_proc
  real(kind=dp)   , intent(inout) :: input_dble

#ifndef serial
  if (nproc > 1) call pbroadcast_dble(input_dble,input_proc)
#endif

end subroutine broadcast_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function pmin(scal)

  real(kind=dp)    :: scal

  pmin = scal
#ifndef serial
  if (nproc > 1) pmin = ppmin(scal)
#endif

end function pmin
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function pmax(scal)

  real(kind=dp)    :: scal

  pmax = scal
#ifndef serial
  if (nproc > 1) pmax = ppmax(scal)
#endif

end function pmax
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function pmax_int(scal)

  integer :: scal

  pmax_int=scal
#ifndef serial
  if (nproc > 1) pmax_int = ppmax_int(scal)
#endif

end function pmax_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=realkind) function psum(scal)

  real(kind=realkind) :: scal

  psum = scal
#ifndef serial
  if (nproc > 1) psum = ppsum(scal)
#endif

end function psum
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function psum_int(scal)

  integer :: scal

  psum_int = scal
#ifndef serial
  if (nproc > 1) psum_int = ppsum_int(scal)
#endif

end function psum_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function psum_dble(scal)

  real(kind=dp)    :: scal

  psum_dble = scal
#ifndef serial
  if (nproc > 1) psum_dble = ppsum_dble(scal)
#endif

end function psum_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine barrier

#ifndef serial
  if (nproc > 1) call pbarrier
#endif

end subroutine barrier
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pcheck(test, errmsg)

  logical, intent(in)            :: test
  character(len=*), intent(in)   :: errmsg

#ifndef serial
  if (nproc > 1) call ppcheck(test, errmsg)
#endif
  if (nproc == 1 .and. test) then
     print '(/,a,/,/,a,/)', 'ERROR:', trim(parse_nl(errmsg))
     stop
  endif

end subroutine pcheck
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine comm_elem_number(my_elems, glob_elems, my_first, my_last, var_name)
! < Communicates the number of elements this processor has to the others and
!! retrieves global number of local first and last element

  use data_io, only: verbose
  integer, intent(in)              :: my_elems
  integer, intent(out)             :: glob_elems, my_first, my_last
  integer                          :: all_elems(0:nproc-1), iproc
  character(len=*), optional       :: var_name

  if (nproc > 1) then
#ifndef serial
     if (verbose > 1 .and. mynum == 0) then
       if (present(var_name)) then
         write(*,"(A,A)") '   Communicating local element numbers of var ', trim(var_name)
       else
         write(*,"(A)") '   Communicating local element numbers'
       endif
       call flush(6)
     endif

     all_elems(mynum) = my_elems

     do iproc = 0, nproc-1
        call pbroadcast_int(all_elems(iproc), iproc)
     enddo

     if (my_elems == 0) then
         my_first = 1
         my_last  = 1
     else
         my_first = sum(all_elems(0:mynum-1)) + 1
         my_last  = sum(all_elems(0:mynum))
     endif

     if (verbose > 1) then
         do iproc = 0, nproc-1
             if (iproc == mynum) then
               if (my_elems == 0) then
                   write(*,"('   Proc:', I5, ' has no elements of this type')") mynum
               else
                   write(*,"('   Proc:', I5, ', first elem:', I10, ', last elem:', I10)" ) &
                         mynum, my_first, my_last
               endif
             endif
             call flush(6)
             call barrier
         enddo
     endif

     glob_elems = sum(all_elems(:))
     call flush(6)
     call barrier
#endif
  else
     my_first = 1
     my_last  = my_elems
     glob_elems = my_elems
  endif

end subroutine comm_elem_number
!-----------------------------------------------------------------------------------------

end module commun
!=========================================================================================
