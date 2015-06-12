!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
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
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!=======================
module clocks_mod
!=======================

!-----------------------------------------------------------------------
!                  CLOCKS module (for timing code sections)
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
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
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

    implicit none

    private

    integer, private              :: ticks_per_sec, max_ticks, ref_tick, start_tick, end_tick
    real(kind=4), private         :: tick_rate
    integer, private, parameter   :: max_clocks = 256
    integer, private              :: clock_num = 0
    logical, private              :: clocks_initialized = .FALSE.

    !clocks are stored in this internal type
    type, private                 :: clock
       character(len=32)          :: name
       integer                    :: ticks, calls
    end type clock

    type(clock)                   :: clocks(0:max_clocks)

    public                        :: clocks_init, clocks_exit, get_clock, clock_id, tick

    character(len=256), private   :: &
         version = '$Id: clocks.F90,v 2.2 2001/02/14 19:06:12 vb Exp $'

  contains

!-----------------------------------------------------------------------
subroutine clocks_init(flag)
    !initialize clocks module
    !if flag is set, only print if flag=0
    !for instance, flag could be set to pe number by the calling program
    !to have only PE 0 in a parallel run print clocks
    integer, intent(in), optional :: flag
    integer :: i
    logical :: verbose

    verbose = .FALSE.

    if( PRESENT(flag) ) verbose = flag.EQ.0

    if( clocks_initialized ) return
    clocks_initialized = .TRUE.

    !initialize clocks and reference tick
    call system_clock( ref_tick, ticks_per_sec, max_ticks )
    tick_rate = 1./ticks_per_sec
    start_tick = ref_tick
    if( verbose )then
        write(6,*) '    CLOCKS module '//trim(version)
        write(6,*) '    Realtime clock resolution=', tick_rate, '(', &
                   ticks_per_sec, ' ticks/sec)'
    end if

    !default clock name is Clock001, etc

    if( verbose ) then
       do i = 1,max_clocks
          write( clocks(i)%name,"(a5,i3.3)") 'Clock', i
       end do
    end if

    clocks%ticks = 0
    clocks%calls = 0

    return
end subroutine clocks_init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function clock_id(name)
    !return an ID for a new or existing clock
    integer                         :: clock_id
    character(len=*), intent(in)    :: name

    if( .NOT.clocks_initialized ) call clocks_init()
    clock_id = 0
    do while( trim(name).NE.trim(clocks(clock_id)%name) )
       clock_id = clock_id + 1
       if( clock_id.GT.clock_num )then
           if( clock_num.EQ.max_clocks )then
               print *, &
                    'CLOCKS ERROR: you are requesting too many clocks, max clocks=', max_clocks
               return
           else
               clock_num = clock_id
               clocks(clock_id)%name = name
               return
           end if
       end if
    end do
    return
end function clock_id
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function tick( string, id, name, since )
    integer                                 :: tick
    character(len=*), intent(in), optional  :: string
    character(len=*), intent(in), optional  :: name
    integer, intent(in), optional           :: id, since
    integer                                 :: current_tick, nid

    !take time first, so that this routine's overhead isn't included
    call system_clock(current_tick)
    if( .NOT.clocks_initialized )call clocks_init()

    !ref_tick is the clock value at the last call to tick (or clocks_init)
    !unless superseded by the since argument.
    if( PRESENT(since) )ref_tick = since

    !correct ref_tick in the unlikely event of clock rollover
    if( current_tick.LT.ref_tick )ref_tick = ref_tick - max_ticks

    if( PRESENT(string) )then
        !print time since reference tick
        print '(a,f14.6)', &
            'CLOCKS: '//trim(string), (current_tick-ref_tick)*tick_rate
    else if( PRESENT(id) )then
        !accumulate time on clock id
        if( 0.LT.id .AND. id.LE.max_clocks )then
            clocks(id)%ticks = clocks(id)%ticks + current_tick - ref_tick
            clocks(id)%calls = clocks(id)%calls + 1
        else
            print *, 'CLOCKS ERROR: invalid id=', id
        end if
    else if( PRESENT(name) )then
        nid = clock_id(name)
    !accumulate time on clock id
        if( 0.LT.nid .AND. nid.LE.max_clocks )then
            clocks(nid)%ticks = clocks(nid)%ticks + current_tick - ref_tick
            clocks(nid)%calls = clocks(nid)%calls + 1
        else
            print *, 'CLOCKS ERROR: invalid id=', nid
        end if
    end if
    !reset reference tick
    call system_clock(ref_tick)
    tick = ref_tick

    return
end function tick
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine get_clock( id, ticks, calls, total_time, time_per_call )
    integer, intent(in)                 :: id
    integer, intent(out), optional      :: ticks, calls
    real(kind=4), intent(out), optional :: total_time, time_per_call

    if( 0.LT.id .AND. id.LE.max_clocks )then
        if( PRESENT(ticks) )ticks = clocks(id)%ticks
        if( PRESENT(calls) )calls = clocks(id)%calls
        if( PRESENT(total_time) )total_time = clocks(id)%ticks*tick_rate
        if( PRESENT(time_per_call) )time_per_call = &
             clocks(id)%ticks*tick_rate/clocks(id)%calls
    else
        print *, 'CLOCKS ERROR: invalid id=', id
    end if

    return
end subroutine get_clock
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine clocks_exit(flag)
    !print all cumulative clocks
    !if flag is set, only print if flag=0
    !for instance, flag could be set to pe number by the calling program
    !to have only PE 0 print clocks
    integer, intent(in), optional :: flag
    integer :: i
    !total_time is for one clock
    !cumul_time is total measured time between clocks_init and clocks_exit
    real(kind=4) :: total_time, time_per_call, cumul_time

    if( PRESENT(flag) )then
        if( flag.NE.0 )return
    end if

    call system_clock(end_tick)
    cumul_time = (end_tick-start_tick)*tick_rate
    write(6,"(32x,a)") '           calls        t_call       t_total t_frac'
    do i = 1, max_clocks
       if( clocks(i)%calls.NE.0 )then
           total_time = clocks(i)%ticks*tick_rate
           time_per_call = total_time/clocks(i)%calls
           write(6,"(a40,i8,2f14.6,f7.3)") &
                'CLOCKS: '//clocks(i)%name, &
                clocks(i)%calls, time_per_call, total_time, total_time/cumul_time
       end if
    end do
    write(6,"(a,f14.6)") 'CLOCKS: Total measured time: ', cumul_time

    return

end subroutine clocks_exit
!--------------------------------------------------------------------

!=======================
end module clocks_mod
!=======================
