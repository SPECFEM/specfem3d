!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_value_integer(value_to_get, name, default_value)

  implicit none

  integer value_to_get, default_value
  character(len=*) name

  call unused_string(name)

  value_to_get = default_value

  end subroutine get_value_integer

!--------------------

  subroutine get_value_double_precision(value_to_get, name, default_value)

  implicit none

  double precision value_to_get, default_value
  character(len=*) name

  call unused_string(name)

  value_to_get = default_value

  end subroutine get_value_double_precision

!--------------------

  subroutine get_value_logical(value_to_get, name, default_value)

  implicit none

  logical value_to_get, default_value
  character(len=*) name

  call unused_string(name)

  value_to_get = default_value

  end subroutine get_value_logical

!--------------------

  subroutine get_value_string(value_to_get, name, default_value)

  implicit none

  character(len=*) value_to_get, default_value
  character(len=*) name

  call unused_string(name)

  value_to_get = default_value

  end subroutine get_value_string
