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

!=================
module linked_list
!=================

  use global_parameters
  implicit none

  private
  public :: list, link

  type :: link
     private
     type(link), pointer        :: next => null()   ! next link in list
     type(link), pointer        :: prev => null()   ! next link in list
     real(kind=realkind), public, allocatable :: ldata(:,:)
     contains
     procedure, pass :: getNextLink  ! return next pointer
     procedure, pass :: setNextLink  ! set next pointer
     procedure, pass :: getPrevLink  ! return next pointer
     procedure, pass :: setPrevLink  ! set next pointer
     procedure, pass :: init
  end type link

  type :: list
     private
     type(link), pointer :: firstLink => null()    ! first link in list
     type(link), pointer :: lastLink => null()     ! last link in list
     type(link), pointer :: currentLink => null()  ! iterator
     integer              :: length = 0
     contains
     procedure, pass :: free            ! empty the list and free the memory
     procedure, pass :: append          ! append an element to the end of the list
     procedure, pass :: insert          ! insert an element to beginning of the list
     procedure, pass :: getFirst        ! return first element and set iterator
                                        ! to first element
     procedure, pass :: getLast         ! return last element and set iterator
                                        ! to last element
     procedure, pass :: getCurrent      ! iterator, can be moved with getNext
                                        ! if not set, returns first element
     procedure, pass :: resetCurrent    ! reset the iterator
     procedure, pass :: getNext         ! get the next element
                                        ! if current not set, returns first element
     procedure, pass :: getPrev         ! get the previous element
                                        ! if current not set, returns last element
     procedure, pass :: getLength       ! return length of the list
     procedure, pass :: eol             ! return .true. if current iterator element is
                                        ! the last element
     procedure, pass :: bol             ! return .true. if current iterator element is
                                        ! the first element
  end type list

contains

!-----------------------------------------------------------------------------------------
function getNextLink(this)
  class(link)           :: this
  type(link), pointer   :: getNextLink
  getNextLink => this%next
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setNextLink(this,next)
  class(link)           :: this
  type(link), pointer   :: next
  this%next => next
end subroutine setNextLink
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getPrevLink(this)
  class(link)           :: this
  type(link), pointer   :: getPrevLink
  getPrevLink => this%prev
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setPrevLink(this, prev)
  class(link)           :: this
  type(link), pointer   :: prev
  this%prev => prev
end subroutine setPrevLink
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init(this, ldata, next, prev)
  class(link)             :: this
  real(kind=realkind)     :: ldata(:,:)
  type(link), pointer     :: next, prev

  this%next => next
  this%prev => prev
  ! Make a copy of the data !

  allocate(this%ldata(size(ldata, dim=1),size(ldata, dim=2)))
  this%ldata = ldata
end subroutine init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine append(this, ldata)
  class(list)             :: this
  real(kind=realkind)     :: ldata(:,:)
  type(link), pointer     :: newLink

  allocate(newLink)

  if (.not. associated(this%firstLink)) then
     call newLink%init(ldata, null(), null())
     this%firstLink => newLink
     this%lastLink => this%firstLink
  else
     call newLink%init(ldata, null(), this%lastLink)
     call this%lastLink%setNextLink(newLink)
     this%lastLink => newLink
  end if
   
  this%length = this%length + 1
end subroutine append
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine insert(this, ldata)
  class(list)             :: this
  real(kind=realkind)     :: ldata(:,:)
  type(link), pointer     :: newLink

  allocate(newLink)
  
  if (.not. associated(this%firstLink)) then
     call newLink%init(ldata, null(), null())
     this%firstLink => newLink
     this%lastLink => this%firstLink
  else
     call newLink%init(ldata, this%firstLink, null())
     call this%firstLink%setPrevLink(newLink)
     this%firstLink => newLink
  end if

  this%length = this%length + 1
end subroutine insert
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getFirst(this)
  class(list)           :: this
  class(link),pointer   :: getFirst

  if(associated(this%firstLink)) then
     getFirst => this%firstLink
     this%currentLink => this%firstLink
  else
     stop 'trying to go access data, but list is empty'
  endif
end function getFirst
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getLast(this)
  class(list)           :: this
  class(link),pointer   :: getLast

  if(associated(this%lastLink)) then
     getLast => this%lastLink
     this%currentLink => this%lastLink
  else
     stop 'trying to go access data, but list is empty'
  endif
end function getLast
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getCurrent(this)
  class(list)           :: this
  class(link),pointer   :: getCurrent

  if (.not. associated(this%currentLink)) then
     if(associated(this%firstLink)) then
        this%currentLink => this%firstLink
     else
        stop 'trying to go access data, but list is empty'
     endif
  endif
  getCurrent => this%currentLink
end function getCurrent
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine resetCurrent(this)
  class(list)           :: this
  this%currentLink => null()
end subroutine resetCurrent
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getNext(this)
  class(list)           :: this
  class(link),pointer   :: getNext

  if (.not. associated(this%currentLink)) then
     this%currentLink => this%firstLink
     getNext => this%currentLink
  elseif (associated(this%currentLink%getNextLink())) then
     this%currentLink => this%currentLink%getNextLink()
     getNext => this%currentLink
  else
     stop 'trying to go beyond last element in list'
  end if 
end function getNext
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function getPrev(this)
  class(list)           :: this
  class(link),pointer   :: getPrev

  if (.not. associated(this%currentLink)) then
     this%currentLink => this%lastLink
     getPrev => this%currentLink
  elseif (associated(this%currentLink%getPrevLink())) then
     this%currentLink => this%currentLink%getPrevLink()
     getPrev => this%currentLink
  else
     stop 'trying to go beyond first element in list'
  end if 
end function getPrev
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function getLength(this)
  class(list)           :: this
  getLength = this%length
end function getLength
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
logical function eol(this)
  class(list)           :: this
  ! MvD: not sure about logics here: currentLink => null beeing the default,
  !      same for first and last link, so what about initialization of current link?
  eol = associated(this%currentLink, target=this%lastLink) .or. (this%length == 0)
end function eol
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
logical function bol(this)
  class(list)           :: this
  ! MvD: same as in eol
  bol = associated(this%currentLink, target=this%firstLink) .or. (this%length == 0)
end function bol
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine free(this)
  class(list)               :: this
  class(link), pointer      :: current, next

  if(associated(this%firstLink)) then
     next => this%firstLink
     do while ( associated(next) )
        current => next
        next => current%getNextLink()
        deallocate( current%ldata )
        deallocate( current )
     enddo
     this%firstLink => null()
  endif

  this%length = 0
end subroutine free
!-----------------------------------------------------------------------------------------

!=====================
end module linked_list
!=====================


!program test_list
!  use linked_list
!  use global_parameters
!  implicit none
!
!  type(list)            :: l
!  class(link), pointer  :: ll
!  real(kind=realkind)   :: a(2,2), b(1,2)
!
!  a = 0
!  b = 1
!  write(6,*) 'eol', l%eol()
!  write(6,*) 'bol', l%bol()
!  call l%append(a)
!  call l%append(b)
!
!  write(6,*) 'first'
!  ll => l%getFirst()
!  write(6,*) ll%ldata
!  a = 2
!  call l%append(a(1:1,1:1))
!  a = -1
!  call l%insert(a)
!  ll => l%getFirst()
!
!  write(6,*) ll%ldata
!  ll => l%getNext()
!  write(6,*) ll%ldata
!
!  call l%resetCurrent()
!  do while(.not. l%eol())
!    ll => l%getNext()
!    write(6,*) 'fw loop', ll%ldata
!    ll%ldata = 10
!  enddo
!
!  call l%resetCurrent()
!  do while(.not. l%eol())
!    ll => l%getNext()
!    write(6,*) 'fw loop', ll%ldata
!  enddo
!end program
