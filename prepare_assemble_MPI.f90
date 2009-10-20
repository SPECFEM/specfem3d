!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


  subroutine prepare_assemble_MPI (nelmnts,ibool, &
                                   knods, ngnode, &
                                   npoin, &
                                   ninterface, max_interface_size, &
                                   my_nelmnts_neighbours, my_interfaces, &
                                   ibool_interfaces_asteroid, &
                                   nibool_interfaces_asteroid &
                                   )

! returns ibool_interfaces_asteroid with the global indices (as defined in ibool) &
! returns nibool_interfaces_asteroid with the number of points in ibool_interfaces_asteroid
! for all points on the (surface) interface defined by knods, ninterface,my_nelmnts_neighbours and my_interfaces

  implicit none

  include 'constants.h'

  integer, intent(in)  :: nelmnts, npoin, ngnode
  integer, dimension(ngnode,nelmnts), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLY,NGLLZ,nelmnts), intent(in)  :: ibool

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_nelmnts_neighbours
  integer, dimension(6,max_interface_size,ninterface)  :: my_interfaces
  integer, dimension(NGLLX*NGLLX*max_interface_size,ninterface)  :: &
       ibool_interfaces_asteroid
  integer, dimension(ninterface)  :: &
       nibool_interfaces_asteroid

  integer  :: num_interface
  integer  :: ispec_interface

  logical, dimension(:),allocatable  :: mask_ibool_asteroid

  integer  :: ixmin, ixmax
  integer  :: iymin, iymax
  integer  :: izmin, izmax
  integer, dimension(ngnode)  :: n
  integer  :: e1, e2, e3, e4
  integer  :: type
  integer  :: ispec

  integer  :: k
  integer  :: npoin_interface_asteroid

  integer  :: ix,iy,iz,ier

  allocate( mask_ibool_asteroid(npoin), stat=ier); if( ier /= 0) stop 'error allocating array'

  ibool_interfaces_asteroid(:,:) = 0
  nibool_interfaces_asteroid(:) = 0

  do num_interface = 1, ninterface
     npoin_interface_asteroid = 0
     mask_ibool_asteroid(:) = .false.

     do ispec_interface = 1, my_nelmnts_neighbours(num_interface)
        ! element with an interface
        ispec = my_interfaces(1,ispec_interface,num_interface)
        ! type of interface
        type = my_interfaces(2,ispec_interface,num_interface)
        ! nodes of face/edge
        do k = 1, ngnode
           n(k) = knods(k,ispec)
        end do
        e1 = my_interfaces(3,ispec_interface,num_interface)
        e2 = my_interfaces(4,ispec_interface,num_interface)
        e3 = my_interfaces(5,ispec_interface,num_interface)
        e4 = my_interfaces(6,ispec_interface,num_interface)
        call get_edge(ngnode, n, type, e1, e2, e3, e4, ixmin, ixmax, iymin, iymax, izmin, izmax)

        do iz = min(izmin,izmax), max(izmin,izmax)
           do iy = min(iymin,iymax), max(iymin,iymax)
              do ix = min(ixmin,ixmax), max(ixmin,ixmax)

                 if(.not. mask_ibool_asteroid(ibool(ix,iy,iz,ispec))) then
                    mask_ibool_asteroid(ibool(ix,iy,iz,ispec)) = .true.
                    npoin_interface_asteroid = npoin_interface_asteroid + 1
                    ibool_interfaces_asteroid(npoin_interface_asteroid,num_interface)=&
                         ibool(ix,iy,iz,ispec)
                 end if
              end do
           end do
        end do

     end do
     nibool_interfaces_asteroid(num_interface) = npoin_interface_asteroid


  end do

  deallocate( mask_ibool_asteroid )
  
end subroutine prepare_assemble_MPI

!
!----
!

subroutine get_edge ( ngnode, n, type, e1, e2, e3, e4, ixmin, ixmax, iymin, iymax, izmin, izmax )

  implicit none

  include "constants.h"

  integer, intent(in)  :: ngnode
  integer, dimension(ngnode), intent(in)  :: n
  integer, intent(in)  :: type, e1, e2, e3, e4
  integer, intent(out)  :: ixmin, ixmax, iymin, iymax, izmin, izmax

  integer, dimension(4) :: en
  integer :: valence, i

   if ( type == 1 ) then
     if ( e1 == n(1) ) then
        ixmin = 1
        ixmax = 1
        iymin = 1
        iymax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = 1
        iymax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = NGLLY
        iymax = NGLLY
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(4) ) then
        ixmin = 1
        ixmax = 1
        iymin = NGLLY
        iymax = NGLLY
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(5) ) then
        ixmin = 1
        ixmax = 1
        iymin = 1
        iymax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(6) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = 1
        iymax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(7) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = NGLLY
        iymax = NGLLY
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(8) ) then
        ixmin = 1
        ixmax = 1
        iymin = NGLLY
        iymax = NGLLY
        izmin = NGLLZ
        izmax = NGLLZ
     end if
  else
     if ( type == 2 ) then
        if ( e1 ==  n(1) ) then
           ixmin = 1
           iymin = 1
           izmin = 1
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(2) ) then
           ixmin = NGLLX
           iymin = 1
           izmin = 1
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if

        end if
        if ( e1 == n(3) ) then
           ixmin = NGLLX
           iymin = NGLLY
           izmin = 1
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(4) ) then
           ixmin = 1
           iymin = NGLLY
           izmin = 1
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(5) ) then
           ixmin = 1
           iymin = 1
           izmin = NGLLZ
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(6) ) then
           ixmin = NGLLX
           iymin = 1
           izmin = NGLLZ
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(7) ) then
           ixmin = NGLLX
           iymin = NGLLY
           izmin = NGLLZ
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(8) ) then
           ixmin = 1
           iymin = NGLLY
           izmin = NGLLZ
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if

     else
        if (type == 4) then
           en(1) = e1
           en(2) = e2
           en(3) = e3
           en(4) = e4

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = NGLLX
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLZ
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = NGLLY
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = NGLLZ
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           endif

        else
           stop 'ERROR get_edge'
        endif

     end if
  end if

end subroutine get_edge

