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

!> General variables pertaining to process identification
!===================
module data_proc
!===================
!

  implicit none
  public 

  integer          :: nproc              !< Number of total processors
  integer          :: mynum              !< Local processor label, from 0 to nproc-1
  character(len=4) :: appnproc, appmynum !< processor-identifying file extension
  logical          :: lpr                !< last processor logical flag, for write stdout 
  character(len=8) :: procstrg           !< String containing mynum to include in writes

!=======================
end module data_proc
!=======================
