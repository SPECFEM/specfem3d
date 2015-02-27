!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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


subroutine parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

  use postprocess_par,only: MAX_STRING_LEN, MAX_KERNEL_NAMES

  implicit none

  character,parameter :: delimiter = ','
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited, kernel_names(MAX_KERNEL_NAMES)
  integer :: iker, nker

  iker = 1
  call strtok(kernel_names_comma_delimited, delimiter, kernel_names(iker))
  do while (kernel_names(iker) /= char(0))
     iker = iker + 1
     call strtok(char(0), delimiter, kernel_names(iker))
  enddo
  nker = iker-1

end subroutine


!
!-------------------------------------------------------------------------------------------------
!

! The following utility function was modified from http://fortranwiki.org/fortran/show/strtok
!
subroutine strtok (source_string, delimiter, token)

!     @(#) Tokenize a string in a similar manner to C routine strtok(3c).
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in delimiter.
!
!             then, if the returned value is not equal to char(0), keep calling until it is
!             with SOURCE_STRING set to char(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning char(0).
!
!     Input:  source_string =   Source string to tokenize.
!             delimiter    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!     Output: strtok()
!
!     LIMITATIONS:
!     can not be called with a different string until current string is totally processed, even from different procedures

      use postprocess_par,only: MAX_STRING_LEN

!     PARAMETERS:
      character(len=MAX_STRING_LEN), intent(in)  :: source_string
      character(len=1), intent(in)  :: delimiter
      character(len=MAX_STRING_LEN), intent(out) :: token

!     SAVED VALUES:
      character(len=MAX_STRING_LEN),save :: saved_string
      integer,save :: isaved_start  ! points to beginning of unprocessed data
      integer,save :: isource_len   ! length of original input string

!     LOCAL VALUES:
      integer :: ibegin        ! beginning of token to return
      integer :: ifinish       ! end of token to return

      ! initialize stored copy of input string and pointer into input string on first call
      if (source_string(1:1) /= char(0)) then
          isaved_start = 1                 ! beginning of unprocessed data
          saved_string = source_string     ! save input string from first call in series
          isource_len = LEN(saved_string)  ! length of input string from first call
      endif

      ibegin = isaved_start

      do
         if ( (ibegin <= isource_len) .AND. (index(delimiter,saved_string(ibegin:ibegin)) /= 0)) then
             ibegin = ibegin + 1
         else
             exit
         endif
      enddo

      if (ibegin > isource_len) then
          token = char(0)
          RETURN
      endif

      ifinish = ibegin

      do
         if ((ifinish <= isource_len) .AND.  (index(delimiter,saved_string(ifinish:ifinish)) == 0)) then
             ifinish = ifinish + 1
         else
             exit
         endif
      enddo

      !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
      token = saved_string(ibegin:ifinish-1)
      isaved_start = ifinish

end subroutine strtok

