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
module clocks_mod

!-----------------------------------------------------------------------
!                  CLOCKS module (for timing code sections)
!
!
! ORIGINAL AUTHOR: V. Balaji (vb@gfdl.gov)
!                  SGI/GFDL Princeton University
!
! ADAPTED BY:      S. Stahler and M. v. Driel
!                  to adhere Fortran 2003 standards
!
!-----------------------------------------------------------------------

    implicit none

    private

    integer, parameter            :: sp = selected_real_kind(6, 37)
    integer, parameter            :: dp = selected_real_kind(15, 307)
    integer, parameter            :: qp = selected_real_kind(33, 4931)

    integer                       :: ticks_per_sec, max_ticks, ref_tick, start_tick, end_tick
    real(kind=dp)                 :: tick_rate
    integer, parameter            :: max_clocks = 256
    integer                       :: clock_num = 0
    logical                       :: clocks_initialized = .false.

    !clocks are stored in this internal type
    type                          :: clock
       character(len=32)          :: name
       integer                    :: ticks, calls
    end type clock

    type(clock)                   :: clocks(0:max_clocks)

    public                        :: clocks_init, clocks_exit, get_clock, clock_id, tick

    character(len=256), private   :: &
         version = 'Clocks.F90,v 3.0 2014/10/19'

  contains


!-----------------------------------------------------------------------------------------
subroutine clocks_init(flag)
    !initialize clocks module
    !if flag is set, only print if flag=0
    !for instance, flag could be set to pe number by the calling program
    !to have only PE 0 in a parallel run print clocks
    integer, intent(in), optional :: flag
    integer :: i
    logical :: verbose

    verbose = .false.

    if ( PRESENT(flag) ) verbose = flag == 0

    if ( clocks_initialized ) return
    clocks_initialized = .true.

    !initialize clocks and reference tick
    call system_clock( ref_tick, ticks_per_sec, max_ticks )
    tick_rate = 1./ticks_per_sec
    start_tick = ref_tick
    if ( verbose ) then
        write(*,*) '    CLOCKS module '//trim(version)
        write(*,*) '    Realtime clock resolution=', tick_rate, '(', &
                   ticks_per_sec, ' ticks/sec)'
    endif

    !default clock name is Clock001, etc

    if ( verbose ) then
       do i = 1,max_clocks
          write( clocks(i)%name,"(a5,i3.3)") 'Clock', i
       enddo
    endif

    clocks%ticks = 0
    clocks%calls = 0

    return
end subroutine clocks_init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function clock_id(name)
    !return an ID for a new or existing clock
    integer                         :: clock_id
    character(len=*), intent(in)    :: name

    if (.not. clocks_initialized ) call clocks_init()
    clock_id = 0
    do while( trim(name) /= trim(clocks(clock_id)%name) )
       clock_id = clock_id + 1
       if ( clock_id > clock_num ) then
           if ( clock_num == max_clocks ) then
               print *, 'CLOCKS ERROR: you are requesting too many clocks, max clocks=', &
                         max_clocks
               return
           else
               clock_num = clock_id
               clocks(clock_id)%name = name
               return
           endif
       endif
    enddo
    return
end function clock_id
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function tick( string, id, name, since )
    integer                                 :: tick
    character(len=*), intent(in), optional  :: string
    character(len=*), intent(in), optional  :: name
    integer, intent(in), optional           :: id, since
    integer                                 :: current_tick, nid

    !take time first, so that this routine's overhead isn't included
    call system_clock(current_tick)
    if (.not. clocks_initialized ) call clocks_init()

    !ref_tick is the clock value at the last call to tick (or clocks_init)
    !unless superseded by the since argument.
    if ( PRESENT(since) )ref_tick = since

    !correct ref_tick in the unlikely event of clock rollover
    if ( current_tick < ref_tick )ref_tick = ref_tick - max_ticks

    if ( PRESENT(string) ) then
        !print time since reference tick
        print '(a,f14.6)', &
            'CLOCKS: '//trim(string), (current_tick-ref_tick)*tick_rate
    else if ( PRESENT(id) ) then
        !accumulate time on clock id
        if ( 0 < id .and. id <= max_clocks ) then
            clocks(id)%ticks = clocks(id)%ticks + current_tick - ref_tick
            clocks(id)%calls = clocks(id)%calls + 1
        else
            print *, 'CLOCKS ERROR: invalid id=', id
        endif
    else if ( PRESENT(name) ) then
        nid = clock_id(name)
    !accumulate time on clock id
        if ( 0 < nid .and. nid <= max_clocks ) then
            clocks(nid)%ticks = clocks(nid)%ticks + current_tick - ref_tick
            clocks(nid)%calls = clocks(nid)%calls + 1
        else
            print *, 'CLOCKS ERROR: invalid id=', nid
        endif
    endif
    !reset reference tick
    call system_clock(ref_tick)
    tick = ref_tick

    return
end function tick
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_clock( id, ticks, calls, total_time, time_per_call )
    integer, intent(in)                  :: id
    integer, intent(out), optional       :: ticks, calls
    real(kind=dp), intent(out), optional :: total_time, time_per_call

    if ( 0 < id .and. id <= max_clocks ) then
        if ( PRESENT(ticks) )ticks = clocks(id)%ticks
        if ( PRESENT(calls) ) calls = clocks(id)%calls
        if ( PRESENT(total_time) )total_time = clocks(id)%ticks*tick_rate
        if ( PRESENT(time_per_call) )time_per_call = &
             clocks(id)%ticks*tick_rate/clocks(id)%calls
    else
        print *, 'CLOCKS ERROR: invalid id=', id
    endif

    return
end subroutine get_clock
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine clocks_exit(flag)
    !print all cumulative clocks
    !if flag is set, only print if flag=0
    !for instance, flag could be set to pe number by the calling program
    !to have only PE 0 print clocks
    integer, intent(in), optional :: flag
    integer                       :: i
    !total_time is for one clock
    !cumul_time is total measured time between clocks_init and clocks_exit
    real(kind=dp)                 :: total_time, time_per_call, cumul_time

    if ( PRESENT(flag) ) then
        if ( flag /= 0 )return
    endif

    call system_clock(end_tick)
    cumul_time = (end_tick-start_tick)*tick_rate
    write(*,"(32x,a)") '           calls        t_call       t_total t_frac'
    do i = 1, max_clocks
       if ( clocks(i)%calls /= 0 ) then
           total_time = clocks(i)%ticks*tick_rate
           time_per_call = total_time/clocks(i)%calls
           write(*,"(a40,i8,2f14.6,f7.3)") &
                'CLOCKS: '//clocks(i)%name, &
                clocks(i)%calls, time_per_call, total_time, total_time/cumul_time
       endif
    enddo
    write(*,"(a,f14.6)") 'CLOCKS: Total measured time: ', cumul_time

    return

end subroutine clocks_exit
!-----------------------------------------------------------------------------------------

end module clocks_mod
!=========================================================================================
