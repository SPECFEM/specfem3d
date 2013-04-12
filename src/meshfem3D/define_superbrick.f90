!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

! define the superbrick that implements the symmetric four-to-two mesh doubling.
! Generated automatically by a script: UTILS/doubling_brick/define_superbrick.pl

  subroutine define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb

  x_superbrick(1) = 3.d0 / 2.d0
  y_superbrick(1) = 1.d0
  z_superbrick(1) = 2.d0

  x_superbrick(2) = 3.d0 / 2.d0
  y_superbrick(2) = 1.d0
  z_superbrick(2) = 3.d0 / 2.d0

  x_superbrick(3) = 3.d0 / 2.d0
  y_superbrick(3) = 3.d0 / 2.d0
  z_superbrick(3) = 3.d0 / 2.d0

  x_superbrick(4) = 3.d0 / 2.d0
  y_superbrick(4) = 3.d0 / 2.d0
  z_superbrick(4) = 2.d0

  x_superbrick(5) = 2.d0
  y_superbrick(5) = 1.d0
  z_superbrick(5) = 2.d0

  x_superbrick(6) = 2.d0
  y_superbrick(6) = 1.d0
  z_superbrick(6) = 1.d0

  x_superbrick(7) = 2.d0
  y_superbrick(7) = 3.d0 / 2.d0
  z_superbrick(7) = 1.d0

  x_superbrick(8) = 2.d0
  y_superbrick(8) = 3.d0 / 2.d0
  z_superbrick(8) = 2.d0

  x_superbrick(9) = 3.d0 / 2.d0
  y_superbrick(9) = 2.d0
  z_superbrick(9) = 1.d0

  x_superbrick(10) = 3.d0 / 2.d0
  y_superbrick(10) = 2.d0
  z_superbrick(10) = 2.d0

  x_superbrick(11) = 2.d0
  y_superbrick(11) = 2.d0
  z_superbrick(11) = 1.d0 / 2.d0

  x_superbrick(12) = 2.d0
  y_superbrick(12) = 2.d0
  z_superbrick(12) = 2.d0

  x_superbrick(13) = 1.d0
  y_superbrick(13) = 1.d0
  z_superbrick(13) = 1.d0

  x_superbrick(14) = 1.d0
  y_superbrick(14) = 1.d0
  z_superbrick(14) = 1.d0 / 2.d0

  x_superbrick(15) = 1.d0
  y_superbrick(15) = 2.d0
  z_superbrick(15) = 1.d0 / 2.d0

  x_superbrick(16) = 1.d0
  y_superbrick(16) = 2.d0
  z_superbrick(16) = 1.d0

  x_superbrick(17) = 3.d0 / 2.d0
  y_superbrick(17) = 1.d0
  z_superbrick(17) = 1.d0

  x_superbrick(18) = 2.d0
  y_superbrick(18) = 1.d0
  z_superbrick(18) = 1.d0 / 2.d0

  x_superbrick(19) = 1.d0
  y_superbrick(19) = 1.d0
  z_superbrick(19) = 3.d0 / 2.d0

  x_superbrick(20) = 1.d0
  y_superbrick(20) = 1.d0
  z_superbrick(20) = 2.d0

  x_superbrick(21) = 1.d0
  y_superbrick(21) = 3.d0 / 2.d0
  z_superbrick(21) = 3.d0 / 2.d0

  x_superbrick(22) = 1.d0
  y_superbrick(22) = 3.d0 / 2.d0
  z_superbrick(22) = 2.d0

  x_superbrick(23) = 1.d0
  y_superbrick(23) = 2.d0
  z_superbrick(23) = 2.d0

  x_superbrick(24) = 1.d0
  y_superbrick(24) = 1.d0
  z_superbrick(24) = 0.d0

  x_superbrick(25) = 2.d0
  y_superbrick(25) = 1.d0
  z_superbrick(25) = 0.d0

  x_superbrick(26) = 2.d0
  y_superbrick(26) = 2.d0
  z_superbrick(26) = 0.d0

  x_superbrick(27) = 1.d0
  y_superbrick(27) = 2.d0
  z_superbrick(27) = 0.d0

  x_superbrick(28) = 3.d0 / 2.d0
  y_superbrick(28) = 1.d0 / 2.d0
  z_superbrick(28) = 3.d0 / 2.d0

  x_superbrick(29) = 3.d0 / 2.d0
  y_superbrick(29) = 1.d0 / 2.d0
  z_superbrick(29) = 2.d0

  x_superbrick(30) = 2.d0
  y_superbrick(30) = 1.d0 / 2.d0
  z_superbrick(30) = 1.d0

  x_superbrick(31) = 2.d0
  y_superbrick(31) = 1.d0 / 2.d0
  z_superbrick(31) = 2.d0

  x_superbrick(32) = 3.d0 / 2.d0
  y_superbrick(32) = 0.d0
  z_superbrick(32) = 1.d0

  x_superbrick(33) = 3.d0 / 2.d0
  y_superbrick(33) = 0.d0
  z_superbrick(33) = 2.d0

  x_superbrick(34) = 2.d0
  y_superbrick(34) = 0.d0
  z_superbrick(34) = 1.d0 / 2.d0

  x_superbrick(35) = 2.d0
  y_superbrick(35) = 0.d0
  z_superbrick(35) = 2.d0

  x_superbrick(36) = 1.d0
  y_superbrick(36) = 0.d0
  z_superbrick(36) = 1.d0 / 2.d0

  x_superbrick(37) = 1.d0
  y_superbrick(37) = 0.d0
  z_superbrick(37) = 1.d0

  x_superbrick(38) = 1.d0
  y_superbrick(38) = 1.d0 / 2.d0
  z_superbrick(38) = 3.d0 / 2.d0

  x_superbrick(39) = 1.d0
  y_superbrick(39) = 1.d0 / 2.d0
  z_superbrick(39) = 2.d0

  x_superbrick(40) = 1.d0
  y_superbrick(40) = 0.d0
  z_superbrick(40) = 2.d0

  x_superbrick(41) = 2.d0
  y_superbrick(41) = 0.d0
  z_superbrick(41) = 0.d0

  x_superbrick(42) = 1.d0
  y_superbrick(42) = 0.d0
  z_superbrick(42) = 0.d0

  x_superbrick(43) = 1.d0 / 2.d0
  y_superbrick(43) = 1.d0
  z_superbrick(43) = 2.d0

  x_superbrick(44) = 1.d0 / 2.d0
  y_superbrick(44) = 1.d0
  z_superbrick(44) = 3.d0 / 2.d0

  x_superbrick(45) = 1.d0 / 2.d0
  y_superbrick(45) = 3.d0 / 2.d0
  z_superbrick(45) = 3.d0 / 2.d0

  x_superbrick(46) = 1.d0 / 2.d0
  y_superbrick(46) = 3.d0 / 2.d0
  z_superbrick(46) = 2.d0

  x_superbrick(47) = 0.d0
  y_superbrick(47) = 1.d0
  z_superbrick(47) = 2.d0

  x_superbrick(48) = 0.d0
  y_superbrick(48) = 1.d0
  z_superbrick(48) = 1.d0

  x_superbrick(49) = 0.d0
  y_superbrick(49) = 3.d0 / 2.d0
  z_superbrick(49) = 1.d0

  x_superbrick(50) = 0.d0
  y_superbrick(50) = 3.d0 / 2.d0
  z_superbrick(50) = 2.d0

  x_superbrick(51) = 1.d0 / 2.d0
  y_superbrick(51) = 2.d0
  z_superbrick(51) = 1.d0

  x_superbrick(52) = 1.d0 / 2.d0
  y_superbrick(52) = 2.d0
  z_superbrick(52) = 2.d0

  x_superbrick(53) = 0.d0
  y_superbrick(53) = 2.d0
  z_superbrick(53) = 1.d0 / 2.d0

  x_superbrick(54) = 0.d0
  y_superbrick(54) = 2.d0
  z_superbrick(54) = 2.d0

  x_superbrick(55) = 1.d0 / 2.d0
  y_superbrick(55) = 1.d0
  z_superbrick(55) = 1.d0

  x_superbrick(56) = 0.d0
  y_superbrick(56) = 1.d0
  z_superbrick(56) = 1.d0 / 2.d0

  x_superbrick(57) = 0.d0
  y_superbrick(57) = 1.d0
  z_superbrick(57) = 0.d0

  x_superbrick(58) = 0.d0
  y_superbrick(58) = 2.d0
  z_superbrick(58) = 0.d0

  x_superbrick(59) = 1.d0 / 2.d0
  y_superbrick(59) = 1.d0 / 2.d0
  z_superbrick(59) = 3.d0 / 2.d0

  x_superbrick(60) = 1.d0 / 2.d0
  y_superbrick(60) = 1.d0 / 2.d0
  z_superbrick(60) = 2.d0

  x_superbrick(61) = 0.d0
  y_superbrick(61) = 1.d0 / 2.d0
  z_superbrick(61) = 1.d0

  x_superbrick(62) = 0.d0
  y_superbrick(62) = 1.d0 / 2.d0
  z_superbrick(62) = 2.d0

  x_superbrick(63) = 1.d0 / 2.d0
  y_superbrick(63) = 0.d0
  z_superbrick(63) = 1.d0

  x_superbrick(64) = 1.d0 / 2.d0
  y_superbrick(64) = 0.d0
  z_superbrick(64) = 2.d0

  x_superbrick(65) = 0.d0
  y_superbrick(65) = 0.d0
  z_superbrick(65) = 1.d0 / 2.d0

  x_superbrick(66) = 0.d0
  y_superbrick(66) = 0.d0
  z_superbrick(66) = 2.d0

  x_superbrick(67) = 0.d0
  y_superbrick(67) = 0.d0
  z_superbrick(67) = 0.d0

  ibool_superbrick(1, 1) = 2
  ibool_superbrick(2, 1) = 6
  ibool_superbrick(3, 1) = 7
  ibool_superbrick(4, 1) = 3
  ibool_superbrick(5, 1) = 1
  ibool_superbrick(6, 1) = 5
  ibool_superbrick(7, 1) = 8
  ibool_superbrick(8, 1) = 4

  ibool_superbrick(1, 2) = 3
  ibool_superbrick(2, 2) = 7
  ibool_superbrick(3, 2) = 11
  ibool_superbrick(4, 2) = 9
  ibool_superbrick(5, 2) = 4
  ibool_superbrick(6, 2) = 8
  ibool_superbrick(7, 2) = 12
  ibool_superbrick(8, 2) = 10

  ibool_superbrick(1, 3) = 14
  ibool_superbrick(2, 3) = 18
  ibool_superbrick(3, 3) = 11
  ibool_superbrick(4, 3) = 15
  ibool_superbrick(5, 3) = 13
  ibool_superbrick(6, 3) = 17
  ibool_superbrick(7, 3) = 9
  ibool_superbrick(8, 3) = 16

  ibool_superbrick(1, 4) = 19
  ibool_superbrick(2, 4) = 2
  ibool_superbrick(3, 4) = 3
  ibool_superbrick(4, 4) = 21
  ibool_superbrick(5, 4) = 20
  ibool_superbrick(6, 4) = 1
  ibool_superbrick(7, 4) = 4
  ibool_superbrick(8, 4) = 22

  ibool_superbrick(1, 5) = 17
  ibool_superbrick(2, 5) = 18
  ibool_superbrick(3, 5) = 11
  ibool_superbrick(4, 5) = 9
  ibool_superbrick(5, 5) = 2
  ibool_superbrick(6, 5) = 6
  ibool_superbrick(7, 5) = 7
  ibool_superbrick(8, 5) = 3

  ibool_superbrick(1, 6) = 21
  ibool_superbrick(2, 6) = 3
  ibool_superbrick(3, 6) = 9
  ibool_superbrick(4, 6) = 16
  ibool_superbrick(5, 6) = 22
  ibool_superbrick(6, 6) = 4
  ibool_superbrick(7, 6) = 10
  ibool_superbrick(8, 6) = 23

  ibool_superbrick(1, 7) = 13
  ibool_superbrick(2, 7) = 17
  ibool_superbrick(3, 7) = 9
  ibool_superbrick(4, 7) = 16
  ibool_superbrick(5, 7) = 19
  ibool_superbrick(6, 7) = 2
  ibool_superbrick(7, 7) = 3
  ibool_superbrick(8, 7) = 21

  ibool_superbrick(1, 8) = 24
  ibool_superbrick(2, 8) = 25
  ibool_superbrick(3, 8) = 26
  ibool_superbrick(4, 8) = 27
  ibool_superbrick(5, 8) = 14
  ibool_superbrick(6, 8) = 18
  ibool_superbrick(7, 8) = 11
  ibool_superbrick(8, 8) = 15

  ibool_superbrick(1, 9) = 28
  ibool_superbrick(2, 9) = 30
  ibool_superbrick(3, 9) = 6
  ibool_superbrick(4, 9) = 2
  ibool_superbrick(5, 9) = 29
  ibool_superbrick(6, 9) = 31
  ibool_superbrick(7, 9) = 5
  ibool_superbrick(8, 9) = 1

  ibool_superbrick(1, 10) = 32
  ibool_superbrick(2, 10) = 34
  ibool_superbrick(3, 10) = 30
  ibool_superbrick(4, 10) = 28
  ibool_superbrick(5, 10) = 33
  ibool_superbrick(6, 10) = 35
  ibool_superbrick(7, 10) = 31
  ibool_superbrick(8, 10) = 29

  ibool_superbrick(1, 11) = 36
  ibool_superbrick(2, 11) = 34
  ibool_superbrick(3, 11) = 18
  ibool_superbrick(4, 11) = 14
  ibool_superbrick(5, 11) = 37
  ibool_superbrick(6, 11) = 32
  ibool_superbrick(7, 11) = 17
  ibool_superbrick(8, 11) = 13

  ibool_superbrick(1, 12) = 38
  ibool_superbrick(2, 12) = 28
  ibool_superbrick(3, 12) = 2
  ibool_superbrick(4, 12) = 19
  ibool_superbrick(5, 12) = 39
  ibool_superbrick(6, 12) = 29
  ibool_superbrick(7, 12) = 1
  ibool_superbrick(8, 12) = 20

  ibool_superbrick(1, 13) = 32
  ibool_superbrick(2, 13) = 34
  ibool_superbrick(3, 13) = 18
  ibool_superbrick(4, 13) = 17
  ibool_superbrick(5, 13) = 28
  ibool_superbrick(6, 13) = 30
  ibool_superbrick(7, 13) = 6
  ibool_superbrick(8, 13) = 2

  ibool_superbrick(1, 14) = 37
  ibool_superbrick(2, 14) = 32
  ibool_superbrick(3, 14) = 28
  ibool_superbrick(4, 14) = 38
  ibool_superbrick(5, 14) = 40
  ibool_superbrick(6, 14) = 33
  ibool_superbrick(7, 14) = 29
  ibool_superbrick(8, 14) = 39

  ibool_superbrick(1, 15) = 37
  ibool_superbrick(2, 15) = 32
  ibool_superbrick(3, 15) = 17
  ibool_superbrick(4, 15) = 13
  ibool_superbrick(5, 15) = 38
  ibool_superbrick(6, 15) = 28
  ibool_superbrick(7, 15) = 2
  ibool_superbrick(8, 15) = 19

  ibool_superbrick(1, 16) = 42
  ibool_superbrick(2, 16) = 41
  ibool_superbrick(3, 16) = 25
  ibool_superbrick(4, 16) = 24
  ibool_superbrick(5, 16) = 36
  ibool_superbrick(6, 16) = 34
  ibool_superbrick(7, 16) = 18
  ibool_superbrick(8, 16) = 14

  ibool_superbrick(1, 17) = 48
  ibool_superbrick(2, 17) = 44
  ibool_superbrick(3, 17) = 45
  ibool_superbrick(4, 17) = 49
  ibool_superbrick(5, 17) = 47
  ibool_superbrick(6, 17) = 43
  ibool_superbrick(7, 17) = 46
  ibool_superbrick(8, 17) = 50

  ibool_superbrick(1, 18) = 49
  ibool_superbrick(2, 18) = 45
  ibool_superbrick(3, 18) = 51
  ibool_superbrick(4, 18) = 53
  ibool_superbrick(5, 18) = 50
  ibool_superbrick(6, 18) = 46
  ibool_superbrick(7, 18) = 52
  ibool_superbrick(8, 18) = 54

  ibool_superbrick(1, 19) = 56
  ibool_superbrick(2, 19) = 14
  ibool_superbrick(3, 19) = 15
  ibool_superbrick(4, 19) = 53
  ibool_superbrick(5, 19) = 55
  ibool_superbrick(6, 19) = 13
  ibool_superbrick(7, 19) = 16
  ibool_superbrick(8, 19) = 51

  ibool_superbrick(1, 20) = 44
  ibool_superbrick(2, 20) = 19
  ibool_superbrick(3, 20) = 21
  ibool_superbrick(4, 20) = 45
  ibool_superbrick(5, 20) = 43
  ibool_superbrick(6, 20) = 20
  ibool_superbrick(7, 20) = 22
  ibool_superbrick(8, 20) = 46

  ibool_superbrick(1, 21) = 56
  ibool_superbrick(2, 21) = 55
  ibool_superbrick(3, 21) = 51
  ibool_superbrick(4, 21) = 53
  ibool_superbrick(5, 21) = 48
  ibool_superbrick(6, 21) = 44
  ibool_superbrick(7, 21) = 45
  ibool_superbrick(8, 21) = 49

  ibool_superbrick(1, 22) = 45
  ibool_superbrick(2, 22) = 21
  ibool_superbrick(3, 22) = 16
  ibool_superbrick(4, 22) = 51
  ibool_superbrick(5, 22) = 46
  ibool_superbrick(6, 22) = 22
  ibool_superbrick(7, 22) = 23
  ibool_superbrick(8, 22) = 52

  ibool_superbrick(1, 23) = 55
  ibool_superbrick(2, 23) = 13
  ibool_superbrick(3, 23) = 16
  ibool_superbrick(4, 23) = 51
  ibool_superbrick(5, 23) = 44
  ibool_superbrick(6, 23) = 19
  ibool_superbrick(7, 23) = 21
  ibool_superbrick(8, 23) = 45

  ibool_superbrick(1, 24) = 57
  ibool_superbrick(2, 24) = 24
  ibool_superbrick(3, 24) = 27
  ibool_superbrick(4, 24) = 58
  ibool_superbrick(5, 24) = 56
  ibool_superbrick(6, 24) = 14
  ibool_superbrick(7, 24) = 15
  ibool_superbrick(8, 24) = 53

  ibool_superbrick(1, 25) = 61
  ibool_superbrick(2, 25) = 59
  ibool_superbrick(3, 25) = 44
  ibool_superbrick(4, 25) = 48
  ibool_superbrick(5, 25) = 62
  ibool_superbrick(6, 25) = 60
  ibool_superbrick(7, 25) = 43
  ibool_superbrick(8, 25) = 47

  ibool_superbrick(1, 26) = 65
  ibool_superbrick(2, 26) = 63
  ibool_superbrick(3, 26) = 59
  ibool_superbrick(4, 26) = 61
  ibool_superbrick(5, 26) = 66
  ibool_superbrick(6, 26) = 64
  ibool_superbrick(7, 26) = 60
  ibool_superbrick(8, 26) = 62

  ibool_superbrick(1, 27) = 65
  ibool_superbrick(2, 27) = 36
  ibool_superbrick(3, 27) = 14
  ibool_superbrick(4, 27) = 56
  ibool_superbrick(5, 27) = 63
  ibool_superbrick(6, 27) = 37
  ibool_superbrick(7, 27) = 13
  ibool_superbrick(8, 27) = 55

  ibool_superbrick(1, 28) = 59
  ibool_superbrick(2, 28) = 38
  ibool_superbrick(3, 28) = 19
  ibool_superbrick(4, 28) = 44
  ibool_superbrick(5, 28) = 60
  ibool_superbrick(6, 28) = 39
  ibool_superbrick(7, 28) = 20
  ibool_superbrick(8, 28) = 43

  ibool_superbrick(1, 29) = 65
  ibool_superbrick(2, 29) = 63
  ibool_superbrick(3, 29) = 55
  ibool_superbrick(4, 29) = 56
  ibool_superbrick(5, 29) = 61
  ibool_superbrick(6, 29) = 59
  ibool_superbrick(7, 29) = 44
  ibool_superbrick(8, 29) = 48

  ibool_superbrick(1, 30) = 63
  ibool_superbrick(2, 30) = 37
  ibool_superbrick(3, 30) = 38
  ibool_superbrick(4, 30) = 59
  ibool_superbrick(5, 30) = 64
  ibool_superbrick(6, 30) = 40
  ibool_superbrick(7, 30) = 39
  ibool_superbrick(8, 30) = 60

  ibool_superbrick(1, 31) = 63
  ibool_superbrick(2, 31) = 37
  ibool_superbrick(3, 31) = 13
  ibool_superbrick(4, 31) = 55
  ibool_superbrick(5, 31) = 59
  ibool_superbrick(6, 31) = 38
  ibool_superbrick(7, 31) = 19
  ibool_superbrick(8, 31) = 44

  ibool_superbrick(1, 32) = 67
  ibool_superbrick(2, 32) = 42
  ibool_superbrick(3, 32) = 24
  ibool_superbrick(4, 32) = 57
  ibool_superbrick(5, 32) = 65
  ibool_superbrick(6, 32) = 36
  ibool_superbrick(7, 32) = 14
  ibool_superbrick(8, 32) = 56


  iboun_sb(:,:) = .false.

  iboun_sb(1,2) = .true.
  iboun_sb(1,6) = .true.
  iboun_sb(2,2) = .true.
  iboun_sb(2,4) = .true.
  iboun_sb(2,6) = .true.
  iboun_sb(3,4) = .true.
  iboun_sb(4,6) = .true.
  iboun_sb(5,2) = .true.
  iboun_sb(6,4) = .true.
  iboun_sb(6,6) = .true.
  iboun_sb(8,2) = .true.
  iboun_sb(8,4) = .true.
  iboun_sb(8,5) = .true.
  iboun_sb(9,2) = .true.
  iboun_sb(9,6) = .true.
  iboun_sb(10,2) = .true.
  iboun_sb(10,3) = .true.
  iboun_sb(10,6) = .true.
  iboun_sb(11,3) = .true.
  iboun_sb(12,6) = .true.
  iboun_sb(13,2) = .true.
  iboun_sb(14,3) = .true.
  iboun_sb(14,6) = .true.
  iboun_sb(16,2) = .true.
  iboun_sb(16,3) = .true.
  iboun_sb(16,5) = .true.
  iboun_sb(17,1) = .true.
  iboun_sb(17,6) = .true.
  iboun_sb(18,1) = .true.
  iboun_sb(18,4) = .true.
  iboun_sb(18,6) = .true.
  iboun_sb(19,4) = .true.
  iboun_sb(20,6) = .true.
  iboun_sb(21,1) = .true.
  iboun_sb(22,4) = .true.
  iboun_sb(22,6) = .true.
  iboun_sb(24,1) = .true.
  iboun_sb(24,4) = .true.
  iboun_sb(24,5) = .true.
  iboun_sb(25,1) = .true.
  iboun_sb(25,6) = .true.
  iboun_sb(26,1) = .true.
  iboun_sb(26,3) = .true.
  iboun_sb(26,6) = .true.
  iboun_sb(27,3) = .true.
  iboun_sb(28,6) = .true.
  iboun_sb(29,1) = .true.
  iboun_sb(30,3) = .true.
  iboun_sb(30,6) = .true.
  iboun_sb(32,1) = .true.
  iboun_sb(32,3) = .true.
  iboun_sb(32,5) = .true.

  end subroutine define_superbrick

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_superbrick_one_layer(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb

  x_superbrick(1) = 3.d0 / 2.d0
  y_superbrick(1) = 1.d0
  z_superbrick(1) = 1.d0

  x_superbrick(2) = 3.d0 / 2.d0
  y_superbrick(2) = 1.d0
  z_superbrick(2) = 2.d0 / 3.d0

  x_superbrick(3) = 3.d0 / 2.d0
  y_superbrick(3) = 3.d0 / 2.d0
  z_superbrick(3) = 2.d0 / 3.d0

  x_superbrick(4) = 3.d0 / 2.d0
  y_superbrick(4) = 3.d0 / 2.d0
  z_superbrick(4) = 1.d0

  x_superbrick(5) = 2.d0
  y_superbrick(5) = 1.d0
  z_superbrick(5) = 1.d0

  x_superbrick(6) = 2.d0
  y_superbrick(6) = 1.d0
  z_superbrick(6) = 1.d0 / 3.d0

  x_superbrick(7) = 2.d0
  y_superbrick(7) = 3.d0 / 2.d0
  z_superbrick(7) = 1.d0 / 3.d0

  x_superbrick(8) = 2.d0
  y_superbrick(8) = 3.d0 / 2.d0
  z_superbrick(8) = 1.d0

  x_superbrick(9) = 3.d0 / 2.d0
  y_superbrick(9) = 2.d0
  z_superbrick(9) = 1.d0 / 3.d0

  x_superbrick(10) = 3.d0 / 2.d0
  y_superbrick(10) = 2.d0
  z_superbrick(10) = 1.d0

  x_superbrick(11) = 2.d0
  y_superbrick(11) = 2.d0
  z_superbrick(11) = 0.d0

  x_superbrick(12) = 2.d0
  y_superbrick(12) = 2.d0
  z_superbrick(12) = 1.d0

  x_superbrick(13) = 1.d0
  y_superbrick(13) = 1.d0
  z_superbrick(13) = 1.d0 / 3.d0

  x_superbrick(14) = 1.d0
  y_superbrick(14) = 1.d0
  z_superbrick(14) = 0.d0

  x_superbrick(15) = 1.d0
  y_superbrick(15) = 2.d0
  z_superbrick(15) = 0.d0

  x_superbrick(16) = 1.d0
  y_superbrick(16) = 2.d0
  z_superbrick(16) = 1.d0 / 3.d0

  x_superbrick(17) = 3.d0 / 2.d0
  y_superbrick(17) = 1.d0
  z_superbrick(17) = 1.d0 / 3.d0

  x_superbrick(18) = 2.d0
  y_superbrick(18) = 1.d0
  z_superbrick(18) = 0.d0

  x_superbrick(19) = 1.d0
  y_superbrick(19) = 1.d0
  z_superbrick(19) = 2.d0 / 3.d0

  x_superbrick(20) = 1.d0
  y_superbrick(20) = 1.d0
  z_superbrick(20) = 1.d0

  x_superbrick(21) = 1.d0
  y_superbrick(21) = 3.d0 / 2.d0
  z_superbrick(21) = 2.d0 / 3.d0

  x_superbrick(22) = 1.d0
  y_superbrick(22) = 3.d0 / 2.d0
  z_superbrick(22) = 1.d0

  x_superbrick(23) = 1.d0
  y_superbrick(23) = 2.d0
  z_superbrick(23) = 1.d0

  x_superbrick(24) = 3.d0 / 2.d0
  y_superbrick(24) = 1.d0 / 2.d0
  z_superbrick(24) = 2.d0 / 3.d0

  x_superbrick(25) = 3.d0 / 2.d0
  y_superbrick(25) = 1.d0 / 2.d0
  z_superbrick(25) = 1.d0

  x_superbrick(26) = 2.d0
  y_superbrick(26) = 1.d0 / 2.d0
  z_superbrick(26) = 1.d0 / 3.d0

  x_superbrick(27) = 2.d0
  y_superbrick(27) = 1.d0 / 2.d0
  z_superbrick(27) = 1.d0

  x_superbrick(28) = 3.d0 / 2.d0
  y_superbrick(28) = 0.d0
  z_superbrick(28) = 1.d0 / 3.d0

  x_superbrick(29) = 3.d0 / 2.d0
  y_superbrick(29) = 0.d0
  z_superbrick(29) = 1.d0

  x_superbrick(30) = 2.d0
  y_superbrick(30) = 0.d0
  z_superbrick(30) = 0.d0

  x_superbrick(31) = 2.d0
  y_superbrick(31) = 0.d0
  z_superbrick(31) = 1.d0

  x_superbrick(32) = 1.d0
  y_superbrick(32) = 0.d0
  z_superbrick(32) = 0.d0

  x_superbrick(33) = 1.d0
  y_superbrick(33) = 0.d0
  z_superbrick(33) = 1.d0 / 3.d0

  x_superbrick(34) = 1.d0
  y_superbrick(34) = 1.d0 / 2.d0
  z_superbrick(34) = 2.d0 / 3.d0

  x_superbrick(35) = 1.d0
  y_superbrick(35) = 1.d0 / 2.d0
  z_superbrick(35) = 1.d0

  x_superbrick(36) = 1.d0
  y_superbrick(36) = 0.d0
  z_superbrick(36) = 1.d0

  x_superbrick(37) = 1.d0 / 2.d0
  y_superbrick(37) = 1.d0
  z_superbrick(37) = 1.d0

  x_superbrick(38) = 1.d0 / 2.d0
  y_superbrick(38) = 1.d0
  z_superbrick(38) = 2.d0 / 3.d0

  x_superbrick(39) = 1.d0 / 2.d0
  y_superbrick(39) = 3.d0 / 2.d0
  z_superbrick(39) = 2.d0 / 3.d0

  x_superbrick(40) = 1.d0 / 2.d0
  y_superbrick(40) = 3.d0 / 2.d0
  z_superbrick(40) = 1.d0

  x_superbrick(41) = 0.d0
  y_superbrick(41) = 1.d0
  z_superbrick(41) = 1.d0

  x_superbrick(42) = 0.d0
  y_superbrick(42) = 1.d0
  z_superbrick(42) = 1.d0 / 3.d0

  x_superbrick(43) = 0.d0
  y_superbrick(43) = 3.d0 / 2.d0
  z_superbrick(43) = 1.d0 / 3.d0

  x_superbrick(44) = 0.d0
  y_superbrick(44) = 3.d0 / 2.d0
  z_superbrick(44) = 1.d0

  x_superbrick(45) = 1.d0 / 2.d0
  y_superbrick(45) = 2.d0
  z_superbrick(45) = 1.d0 / 3.d0

  x_superbrick(46) = 1.d0 / 2.d0
  y_superbrick(46) = 2.d0
  z_superbrick(46) = 1.d0

  x_superbrick(47) = 0.d0
  y_superbrick(47) = 2.d0
  z_superbrick(47) = 0.d0

  x_superbrick(48) = 0.d0
  y_superbrick(48) = 2.d0
  z_superbrick(48) = 1.d0

  x_superbrick(49) = 1.d0 / 2.d0
  y_superbrick(49) = 1.d0
  z_superbrick(49) = 1.d0 / 3.d0

  x_superbrick(50) = 0.d0
  y_superbrick(50) = 1.d0
  z_superbrick(50) = 0.d0

  x_superbrick(51) = 1.d0 / 2.d0
  y_superbrick(51) = 1.d0 / 2.d0
  z_superbrick(51) = 2.d0 / 3.d0

  x_superbrick(52) = 1.d0 / 2.d0
  y_superbrick(52) = 1.d0 / 2.d0
  z_superbrick(52) = 1.d0

  x_superbrick(53) = 0.d0
  y_superbrick(53) = 1.d0 / 2.d0
  z_superbrick(53) = 1.d0 / 3.d0

  x_superbrick(54) = 0.d0
  y_superbrick(54) = 1.d0 / 2.d0
  z_superbrick(54) = 1.d0

  x_superbrick(55) = 1.d0 / 2.d0
  y_superbrick(55) = 0.d0
  z_superbrick(55) = 1.d0 / 3.d0

  x_superbrick(56) = 1.d0 / 2.d0
  y_superbrick(56) = 0.d0
  z_superbrick(56) = 1.d0

  x_superbrick(57) = 0.d0
  y_superbrick(57) = 0.d0
  z_superbrick(57) = 0.d0

  x_superbrick(58) = 0.d0
  y_superbrick(58) = 0.d0
  z_superbrick(58) = 1.d0

  ibool_superbrick(1, 1) = 2
  ibool_superbrick(2, 1) = 6
  ibool_superbrick(3, 1) = 7
  ibool_superbrick(4, 1) = 3
  ibool_superbrick(5, 1) = 1
  ibool_superbrick(6, 1) = 5
  ibool_superbrick(7, 1) = 8
  ibool_superbrick(8, 1) = 4

  ibool_superbrick(1, 2) = 3
  ibool_superbrick(2, 2) = 7
  ibool_superbrick(3, 2) = 11
  ibool_superbrick(4, 2) = 9
  ibool_superbrick(5, 2) = 4
  ibool_superbrick(6, 2) = 8
  ibool_superbrick(7, 2) = 12
  ibool_superbrick(8, 2) = 10

  ibool_superbrick(1, 3) = 14
  ibool_superbrick(2, 3) = 18
  ibool_superbrick(3, 3) = 11
  ibool_superbrick(4, 3) = 15
  ibool_superbrick(5, 3) = 13
  ibool_superbrick(6, 3) = 17
  ibool_superbrick(7, 3) = 9
  ibool_superbrick(8, 3) = 16

  ibool_superbrick(1, 4) = 19
  ibool_superbrick(2, 4) = 2
  ibool_superbrick(3, 4) = 3
  ibool_superbrick(4, 4) = 21
  ibool_superbrick(5, 4) = 20
  ibool_superbrick(6, 4) = 1
  ibool_superbrick(7, 4) = 4
  ibool_superbrick(8, 4) = 22

  ibool_superbrick(1, 5) = 17
  ibool_superbrick(2, 5) = 18
  ibool_superbrick(3, 5) = 11
  ibool_superbrick(4, 5) = 9
  ibool_superbrick(5, 5) = 2
  ibool_superbrick(6, 5) = 6
  ibool_superbrick(7, 5) = 7
  ibool_superbrick(8, 5) = 3

  ibool_superbrick(1, 6) = 21
  ibool_superbrick(2, 6) = 3
  ibool_superbrick(3, 6) = 9
  ibool_superbrick(4, 6) = 16
  ibool_superbrick(5, 6) = 22
  ibool_superbrick(6, 6) = 4
  ibool_superbrick(7, 6) = 10
  ibool_superbrick(8, 6) = 23

  ibool_superbrick(1, 7) = 13
  ibool_superbrick(2, 7) = 17
  ibool_superbrick(3, 7) = 9
  ibool_superbrick(4, 7) = 16
  ibool_superbrick(5, 7) = 19
  ibool_superbrick(6, 7) = 2
  ibool_superbrick(7, 7) = 3
  ibool_superbrick(8, 7) = 21

  ibool_superbrick(1, 8) = 24
  ibool_superbrick(2, 8) = 26
  ibool_superbrick(3, 8) = 6
  ibool_superbrick(4, 8) = 2
  ibool_superbrick(5, 8) = 25
  ibool_superbrick(6, 8) = 27
  ibool_superbrick(7, 8) = 5
  ibool_superbrick(8, 8) = 1

  ibool_superbrick(1, 9) = 28
  ibool_superbrick(2, 9) = 30
  ibool_superbrick(3, 9) = 26
  ibool_superbrick(4, 9) = 24
  ibool_superbrick(5, 9) = 29
  ibool_superbrick(6, 9) = 31
  ibool_superbrick(7, 9) = 27
  ibool_superbrick(8, 9) = 25

  ibool_superbrick(1, 10) = 32
  ibool_superbrick(2, 10) = 30
  ibool_superbrick(3, 10) = 18
  ibool_superbrick(4, 10) = 14
  ibool_superbrick(5, 10) = 33
  ibool_superbrick(6, 10) = 28
  ibool_superbrick(7, 10) = 17
  ibool_superbrick(8, 10) = 13

  ibool_superbrick(1, 11) = 34
  ibool_superbrick(2, 11) = 24
  ibool_superbrick(3, 11) = 2
  ibool_superbrick(4, 11) = 19
  ibool_superbrick(5, 11) = 35
  ibool_superbrick(6, 11) = 25
  ibool_superbrick(7, 11) = 1
  ibool_superbrick(8, 11) = 20

  ibool_superbrick(1, 12) = 28
  ibool_superbrick(2, 12) = 30
  ibool_superbrick(3, 12) = 18
  ibool_superbrick(4, 12) = 17
  ibool_superbrick(5, 12) = 24
  ibool_superbrick(6, 12) = 26
  ibool_superbrick(7, 12) = 6
  ibool_superbrick(8, 12) = 2

  ibool_superbrick(1, 13) = 33
  ibool_superbrick(2, 13) = 28
  ibool_superbrick(3, 13) = 24
  ibool_superbrick(4, 13) = 34
  ibool_superbrick(5, 13) = 36
  ibool_superbrick(6, 13) = 29
  ibool_superbrick(7, 13) = 25
  ibool_superbrick(8, 13) = 35

  ibool_superbrick(1, 14) = 33
  ibool_superbrick(2, 14) = 28
  ibool_superbrick(3, 14) = 17
  ibool_superbrick(4, 14) = 13
  ibool_superbrick(5, 14) = 34
  ibool_superbrick(6, 14) = 24
  ibool_superbrick(7, 14) = 2
  ibool_superbrick(8, 14) = 19

  ibool_superbrick(1, 15) = 42
  ibool_superbrick(2, 15) = 38
  ibool_superbrick(3, 15) = 39
  ibool_superbrick(4, 15) = 43
  ibool_superbrick(5, 15) = 41
  ibool_superbrick(6, 15) = 37
  ibool_superbrick(7, 15) = 40
  ibool_superbrick(8, 15) = 44

  ibool_superbrick(1, 16) = 43
  ibool_superbrick(2, 16) = 39
  ibool_superbrick(3, 16) = 45
  ibool_superbrick(4, 16) = 47
  ibool_superbrick(5, 16) = 44
  ibool_superbrick(6, 16) = 40
  ibool_superbrick(7, 16) = 46
  ibool_superbrick(8, 16) = 48

  ibool_superbrick(1, 17) = 50
  ibool_superbrick(2, 17) = 14
  ibool_superbrick(3, 17) = 15
  ibool_superbrick(4, 17) = 47
  ibool_superbrick(5, 17) = 49
  ibool_superbrick(6, 17) = 13
  ibool_superbrick(7, 17) = 16
  ibool_superbrick(8, 17) = 45

  ibool_superbrick(1, 18) = 38
  ibool_superbrick(2, 18) = 19
  ibool_superbrick(3, 18) = 21
  ibool_superbrick(4, 18) = 39
  ibool_superbrick(5, 18) = 37
  ibool_superbrick(6, 18) = 20
  ibool_superbrick(7, 18) = 22
  ibool_superbrick(8, 18) = 40

  ibool_superbrick(1, 19) = 50
  ibool_superbrick(2, 19) = 49
  ibool_superbrick(3, 19) = 45
  ibool_superbrick(4, 19) = 47
  ibool_superbrick(5, 19) = 42
  ibool_superbrick(6, 19) = 38
  ibool_superbrick(7, 19) = 39
  ibool_superbrick(8, 19) = 43

  ibool_superbrick(1, 20) = 39
  ibool_superbrick(2, 20) = 21
  ibool_superbrick(3, 20) = 16
  ibool_superbrick(4, 20) = 45
  ibool_superbrick(5, 20) = 40
  ibool_superbrick(6, 20) = 22
  ibool_superbrick(7, 20) = 23
  ibool_superbrick(8, 20) = 46

  ibool_superbrick(1, 21) = 49
  ibool_superbrick(2, 21) = 13
  ibool_superbrick(3, 21) = 16
  ibool_superbrick(4, 21) = 45
  ibool_superbrick(5, 21) = 38
  ibool_superbrick(6, 21) = 19
  ibool_superbrick(7, 21) = 21
  ibool_superbrick(8, 21) = 39

  ibool_superbrick(1, 22) = 53
  ibool_superbrick(2, 22) = 51
  ibool_superbrick(3, 22) = 38
  ibool_superbrick(4, 22) = 42
  ibool_superbrick(5, 22) = 54
  ibool_superbrick(6, 22) = 52
  ibool_superbrick(7, 22) = 37
  ibool_superbrick(8, 22) = 41

  ibool_superbrick(1, 23) = 57
  ibool_superbrick(2, 23) = 55
  ibool_superbrick(3, 23) = 51
  ibool_superbrick(4, 23) = 53
  ibool_superbrick(5, 23) = 58
  ibool_superbrick(6, 23) = 56
  ibool_superbrick(7, 23) = 52
  ibool_superbrick(8, 23) = 54

  ibool_superbrick(1, 24) = 57
  ibool_superbrick(2, 24) = 32
  ibool_superbrick(3, 24) = 14
  ibool_superbrick(4, 24) = 50
  ibool_superbrick(5, 24) = 55
  ibool_superbrick(6, 24) = 33
  ibool_superbrick(7, 24) = 13
  ibool_superbrick(8, 24) = 49

  ibool_superbrick(1, 25) = 51
  ibool_superbrick(2, 25) = 34
  ibool_superbrick(3, 25) = 19
  ibool_superbrick(4, 25) = 38
  ibool_superbrick(5, 25) = 52
  ibool_superbrick(6, 25) = 35
  ibool_superbrick(7, 25) = 20
  ibool_superbrick(8, 25) = 37

  ibool_superbrick(1, 26) = 57
  ibool_superbrick(2, 26) = 55
  ibool_superbrick(3, 26) = 49
  ibool_superbrick(4, 26) = 50
  ibool_superbrick(5, 26) = 53
  ibool_superbrick(6, 26) = 51
  ibool_superbrick(7, 26) = 38
  ibool_superbrick(8, 26) = 42

  ibool_superbrick(1, 27) = 55
  ibool_superbrick(2, 27) = 33
  ibool_superbrick(3, 27) = 34
  ibool_superbrick(4, 27) = 51
  ibool_superbrick(5, 27) = 56
  ibool_superbrick(6, 27) = 36
  ibool_superbrick(7, 27) = 35
  ibool_superbrick(8, 27) = 52

  ibool_superbrick(1, 28) = 55
  ibool_superbrick(2, 28) = 33
  ibool_superbrick(3, 28) = 13
  ibool_superbrick(4, 28) = 49
  ibool_superbrick(5, 28) = 51
  ibool_superbrick(6, 28) = 34
  ibool_superbrick(7, 28) = 19
  ibool_superbrick(8, 28) = 38

  iboun_sb(:,:) = .false.
  iboun_sb(1,2) = .true.
  iboun_sb(1,6) = .true.
  iboun_sb(2,2) = .true.
  iboun_sb(2,4) = .true.
  iboun_sb(2,6) = .true.
  iboun_sb(3,4) = .true.
  iboun_sb(3,5) = .true.
  iboun_sb(4,6) = .true.
  iboun_sb(5,2) = .true.
  iboun_sb(6,4) = .true.
  iboun_sb(6,6) = .true.
  iboun_sb(8,2) = .true.
  iboun_sb(8,6) = .true.
  iboun_sb(9,2) = .true.
  iboun_sb(9,3) = .true.
  iboun_sb(9,6) = .true.
  iboun_sb(10,3) = .true.
  iboun_sb(10,5) = .true.
  iboun_sb(11,6) = .true.
  iboun_sb(12,2) = .true.
  iboun_sb(13,3) = .true.
  iboun_sb(13,6) = .true.
  iboun_sb(15,1) = .true.
  iboun_sb(15,6) = .true.
  iboun_sb(16,1) = .true.
  iboun_sb(16,4) = .true.
  iboun_sb(16,6) = .true.
  iboun_sb(17,4) = .true.
  iboun_sb(17,5) = .true.
  iboun_sb(18,6) = .true.
  iboun_sb(19,1) = .true.
  iboun_sb(20,4) = .true.
  iboun_sb(20,6) = .true.
  iboun_sb(22,1) = .true.
  iboun_sb(22,6) = .true.
  iboun_sb(23,1) = .true.
  iboun_sb(23,3) = .true.
  iboun_sb(23,6) = .true.
  iboun_sb(24,3) = .true.
  iboun_sb(24,5) = .true.
  iboun_sb(25,6) = .true.
  iboun_sb(26,1) = .true.
  iboun_sb(27,3) = .true.
  iboun_sb(27,6) = .true.

  end subroutine define_superbrick_one_layer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_basic_doubling_brick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb,case_num)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
  integer :: case_num

  SELECT CASE (case_num)
      CASE (1)
          x_superbrick(1) = 1.d0 / 2.d0
          y_superbrick(1) = 1.d0
          z_superbrick(1) = 2.d0

          x_superbrick(2) = 1.d0 / 2.d0
          y_superbrick(2) = 1.d0
          z_superbrick(2) = 3.d0 / 2.d0

          x_superbrick(3) = 1.d0 / 2.d0
          y_superbrick(3) = 1.d0 / 2.d0
          z_superbrick(3) = 3.d0 / 2.d0

          x_superbrick(4) = 1.d0 / 2.d0
          y_superbrick(4) = 1.d0 / 2.d0
          z_superbrick(4) = 2.d0

          x_superbrick(5) = 0.d0
          y_superbrick(5) = 1.d0
          z_superbrick(5) = 2.d0

          x_superbrick(6) = 0.d0
          y_superbrick(6) = 1.d0
          z_superbrick(6) = 1.d0

          x_superbrick(7) = 0.d0
          y_superbrick(7) = 1.d0 / 2.d0
          z_superbrick(7) = 1.d0

          x_superbrick(8) = 0.d0
          y_superbrick(8) = 1.d0 / 2.d0
          z_superbrick(8) = 2.d0

          x_superbrick(9) = 1.d0 / 2.d0
          y_superbrick(9) = 0.d0
          z_superbrick(9) = 1.d0

          x_superbrick(10) = 1.d0 / 2.d0
          y_superbrick(10) = 0.d0
          z_superbrick(10) = 2.d0

          x_superbrick(11) = 0.d0
          y_superbrick(11) = 0.d0
          z_superbrick(11) = 1.d0 / 2.d0

          x_superbrick(12) = 0.d0
          y_superbrick(12) = 0.d0
          z_superbrick(12) = 2.d0

          x_superbrick(13) = 1.d0
          y_superbrick(13) = 1.d0
          z_superbrick(13) = 1.d0

          x_superbrick(14) = 1.d0
          y_superbrick(14) = 1.d0
          z_superbrick(14) = 1.d0 / 2.d0

          x_superbrick(15) = 1.d0
          y_superbrick(15) = 0.d0
          z_superbrick(15) = 1.d0 / 2.d0

          x_superbrick(16) = 1.d0
          y_superbrick(16) = 0.d0
          z_superbrick(16) = 1.d0

          x_superbrick(17) = 1.d0 / 2.d0
          y_superbrick(17) = 1.d0
          z_superbrick(17) = 1.d0

          x_superbrick(18) = 0.d0
          y_superbrick(18) = 1.d0
          z_superbrick(18) = 1.d0 / 2.d0

          x_superbrick(19) = 1.d0
          y_superbrick(19) = 1.d0
          z_superbrick(19) = 3.d0 / 2.d0

          x_superbrick(20) = 1.d0
          y_superbrick(20) = 1.d0
          z_superbrick(20) = 2.d0

          x_superbrick(21) = 1.d0
          y_superbrick(21) = 1.d0 / 2.d0
          z_superbrick(21) = 3.d0 / 2.d0

          x_superbrick(22) = 1.d0
          y_superbrick(22) = 1.d0 / 2.d0
          z_superbrick(22) = 2.d0

          x_superbrick(23) = 1.d0
          y_superbrick(23) = 0.d0
          z_superbrick(23) = 2.d0

          x_superbrick(24) = 1.d0
          y_superbrick(24) = 1.d0
          z_superbrick(24) = 0.d0

          x_superbrick(25) = 0.d0
          y_superbrick(25) = 1.d0
          z_superbrick(25) = 0.d0

          x_superbrick(26) = 0.d0
          y_superbrick(26) = 0.d0
          z_superbrick(26) = 0.d0

          x_superbrick(27) = 1.d0
          y_superbrick(27) = 0.d0
          z_superbrick(27) = 0.d0

          ibool_superbrick(1, 1) = 7
          ibool_superbrick(2, 1) = 3
          ibool_superbrick(3, 1) = 2
          ibool_superbrick(4, 1) = 6
          ibool_superbrick(5, 1) = 8
          ibool_superbrick(6, 1) = 4
          ibool_superbrick(7, 1) = 1
          ibool_superbrick(8, 1) = 5

          ibool_superbrick(1, 2) = 11
          ibool_superbrick(2, 2) = 9
          ibool_superbrick(3, 2) = 3
          ibool_superbrick(4, 2) = 7
          ibool_superbrick(5, 2) = 12
          ibool_superbrick(6, 2) = 10
          ibool_superbrick(7, 2) = 4
          ibool_superbrick(8, 2) = 8

          ibool_superbrick(1, 3) = 11
          ibool_superbrick(2, 3) = 15
          ibool_superbrick(3, 3) = 14
          ibool_superbrick(4, 3) = 18
          ibool_superbrick(5, 3) = 9
          ibool_superbrick(6, 3) = 16
          ibool_superbrick(7, 3) = 13
          ibool_superbrick(8, 3) = 17

          ibool_superbrick(1, 4) = 3
          ibool_superbrick(2, 4) = 21
          ibool_superbrick(3, 4) = 19
          ibool_superbrick(4, 4) = 2
          ibool_superbrick(5, 4) = 4
          ibool_superbrick(6, 4) = 22
          ibool_superbrick(7, 4) = 20
          ibool_superbrick(8, 4) = 1

          ibool_superbrick(1, 5) = 11
          ibool_superbrick(2, 5) = 9
          ibool_superbrick(3, 5) = 17
          ibool_superbrick(4, 5) = 18
          ibool_superbrick(5, 5) = 7
          ibool_superbrick(6, 5) = 3
          ibool_superbrick(7, 5) = 2
          ibool_superbrick(8, 5) = 6

          ibool_superbrick(1, 6) = 9
          ibool_superbrick(2, 6) = 16
          ibool_superbrick(3, 6) = 21
          ibool_superbrick(4, 6) = 3
          ibool_superbrick(5, 6) = 10
          ibool_superbrick(6, 6) = 23
          ibool_superbrick(7, 6) = 22
          ibool_superbrick(8, 6) = 4

          ibool_superbrick(1, 7) = 9
          ibool_superbrick(2, 7) = 16
          ibool_superbrick(3, 7) = 13
          ibool_superbrick(4, 7) = 17
          ibool_superbrick(5, 7) = 3
          ibool_superbrick(6, 7) = 21
          ibool_superbrick(7, 7) = 19
          ibool_superbrick(8, 7) = 2

          ibool_superbrick(1, 8) = 26
          ibool_superbrick(2, 8) = 27
          ibool_superbrick(3, 8) = 24
          ibool_superbrick(4, 8) = 25
          ibool_superbrick(5, 8) = 11
          ibool_superbrick(6, 8) = 15
          ibool_superbrick(7, 8) = 14
          ibool_superbrick(8, 8) = 18

          iboun_sb(:,:) = .false.
          iboun_sb(1,1) = .true.
          iboun_sb(1,4) = .true.
          iboun_sb(1,6) = .true.
          iboun_sb(2,1) = .true.
          iboun_sb(2,3) = .true.
          iboun_sb(2,6) = .true.
          iboun_sb(3,2) = .true.
          iboun_sb(3,3) = .true.
          iboun_sb(3,4) = .true.
          iboun_sb(4,2) = .true.
          iboun_sb(4,4) = .true.
          iboun_sb(4,6) = .true.
          iboun_sb(5,1) = .true.
          iboun_sb(5,4) = .true.
          iboun_sb(6,2) = .true.
          iboun_sb(6,3) = .true.
          iboun_sb(6,6) = .true.
          iboun_sb(7,2) = .true.
          iboun_sb(7,4) = .true.
          iboun_sb(8,1) = .true.
          iboun_sb(8,2) = .true.
          iboun_sb(8,3) = .true.
          iboun_sb(8,4) = .true.
          iboun_sb(8,5) = .true.
      CASE (2)
          x_superbrick(1) = 1.d0 / 2.d0
          y_superbrick(1) = 0.d0
          z_superbrick(1) = 2.d0

          x_superbrick(2) = 1.d0 / 2.d0
          y_superbrick(2) = 0.d0
          z_superbrick(2) = 3.d0 / 2.d0

          x_superbrick(3) = 1.d0 / 2.d0
          y_superbrick(3) = 1.d0 / 2.d0
          z_superbrick(3) = 3.d0 / 2.d0

          x_superbrick(4) = 1.d0 / 2.d0
          y_superbrick(4) = 1.d0 / 2.d0
          z_superbrick(4) = 2.d0

          x_superbrick(5) = 0.d0
          y_superbrick(5) = 0.d0
          z_superbrick(5) = 2.d0

          x_superbrick(6) = 0.d0
          y_superbrick(6) = 0.d0
          z_superbrick(6) = 1.d0

          x_superbrick(7) = 0.d0
          y_superbrick(7) = 1.d0 / 2.d0
          z_superbrick(7) = 1.d0

          x_superbrick(8) = 0.d0
          y_superbrick(8) = 1.d0 / 2.d0
          z_superbrick(8) = 2.d0

          x_superbrick(9) = 1.d0 / 2.d0
          y_superbrick(9) = 1.d0
          z_superbrick(9) = 1.d0

          x_superbrick(10) = 1.d0 / 2.d0
          y_superbrick(10) = 1.d0
          z_superbrick(10) = 2.d0

          x_superbrick(11) = 0.d0
          y_superbrick(11) = 1.d0
          z_superbrick(11) = 1.d0 / 2.d0

          x_superbrick(12) = 0.d0
          y_superbrick(12) = 1.d0
          z_superbrick(12) = 2.d0

          x_superbrick(13) = 1.d0
          y_superbrick(13) = 0.d0
          z_superbrick(13) = 1.d0

          x_superbrick(14) = 1.d0
          y_superbrick(14) = 0.d0
          z_superbrick(14) = 1.d0 / 2.d0

          x_superbrick(15) = 1.d0
          y_superbrick(15) = 1.d0
          z_superbrick(15) = 1.d0 / 2.d0

          x_superbrick(16) = 1.d0
          y_superbrick(16) = 1.d0
          z_superbrick(16) = 1.d0

          x_superbrick(17) = 1.d0 / 2.d0
          y_superbrick(17) = 0.d0
          z_superbrick(17) = 1.d0

          x_superbrick(18) = 0.d0
          y_superbrick(18) = 0.d0
          z_superbrick(18) = 1.d0 / 2.d0

          x_superbrick(19) = 1.d0
          y_superbrick(19) = 0.d0
          z_superbrick(19) = 3.d0 / 2.d0

          x_superbrick(20) = 1.d0
          y_superbrick(20) = 0.d0
          z_superbrick(20) = 2.d0

          x_superbrick(21) = 1.d0
          y_superbrick(21) = 1.d0 / 2.d0
          z_superbrick(21) = 3.d0 / 2.d0

          x_superbrick(22) = 1.d0
          y_superbrick(22) = 1.d0 / 2.d0
          z_superbrick(22) = 2.d0

          x_superbrick(23) = 1.d0
          y_superbrick(23) = 1.d0
          z_superbrick(23) = 2.d0

          x_superbrick(24) = 1.d0
          y_superbrick(24) = 0.d0
          z_superbrick(24) = 0.d0

          x_superbrick(25) = 0.d0
          y_superbrick(25) = 0.d0
          z_superbrick(25) = 0.d0

          x_superbrick(26) = 0.d0
          y_superbrick(26) = 1.d0
          z_superbrick(26) = 0.d0

          x_superbrick(27) = 1.d0
          y_superbrick(27) = 1.d0
          z_superbrick(27) = 0.d0

          ibool_superbrick(1, 1) = 6
          ibool_superbrick(2, 1) = 2
          ibool_superbrick(3, 1) = 3
          ibool_superbrick(4, 1) = 7
          ibool_superbrick(5, 1) = 5
          ibool_superbrick(6, 1) = 1
          ibool_superbrick(7, 1) = 4
          ibool_superbrick(8, 1) = 8

          ibool_superbrick(1, 2) = 7
          ibool_superbrick(2, 2) = 3
          ibool_superbrick(3, 2) = 9
          ibool_superbrick(4, 2) = 11
          ibool_superbrick(5, 2) = 8
          ibool_superbrick(6, 2) = 4
          ibool_superbrick(7, 2) = 10
          ibool_superbrick(8, 2) = 12

          ibool_superbrick(1, 3) = 18
          ibool_superbrick(2, 3) = 14
          ibool_superbrick(3, 3) = 15
          ibool_superbrick(4, 3) = 11
          ibool_superbrick(5, 3) = 17
          ibool_superbrick(6, 3) = 13
          ibool_superbrick(7, 3) = 16
          ibool_superbrick(8, 3) = 9

          ibool_superbrick(1, 4) = 2
          ibool_superbrick(2, 4) = 19
          ibool_superbrick(3, 4) = 21
          ibool_superbrick(4, 4) = 3
          ibool_superbrick(5, 4) = 1
          ibool_superbrick(6, 4) = 20
          ibool_superbrick(7, 4) = 22
          ibool_superbrick(8, 4) = 4

          ibool_superbrick(1, 5) = 18
          ibool_superbrick(2, 5) = 17
          ibool_superbrick(3, 5) = 9
          ibool_superbrick(4, 5) = 11
          ibool_superbrick(5, 5) = 6
          ibool_superbrick(6, 5) = 2
          ibool_superbrick(7, 5) = 3
          ibool_superbrick(8, 5) = 7

          ibool_superbrick(1, 6) = 3
          ibool_superbrick(2, 6) = 21
          ibool_superbrick(3, 6) = 16
          ibool_superbrick(4, 6) = 9
          ibool_superbrick(5, 6) = 4
          ibool_superbrick(6, 6) = 22
          ibool_superbrick(7, 6) = 23
          ibool_superbrick(8, 6) = 10

          ibool_superbrick(1, 7) = 17
          ibool_superbrick(2, 7) = 13
          ibool_superbrick(3, 7) = 16
          ibool_superbrick(4, 7) = 9
          ibool_superbrick(5, 7) = 2
          ibool_superbrick(6, 7) = 19
          ibool_superbrick(7, 7) = 21
          ibool_superbrick(8, 7) = 3

          ibool_superbrick(1, 8) = 25
          ibool_superbrick(2, 8) = 24
          ibool_superbrick(3, 8) = 27
          ibool_superbrick(4, 8) = 26
          ibool_superbrick(5, 8) = 18
          ibool_superbrick(6, 8) = 14
          ibool_superbrick(7, 8) = 15
          ibool_superbrick(8, 8) = 11

          iboun_sb(:,:) = .false.
          iboun_sb(1,1) = .true.
          iboun_sb(1,3) = .true.
          iboun_sb(1,6) = .true.
          iboun_sb(2,1) = .true.
          iboun_sb(2,4) = .true.
          iboun_sb(2,6) = .true.
          iboun_sb(3,2) = .true.
          iboun_sb(3,3) = .true.
          iboun_sb(3,4) = .true.
          iboun_sb(4,2) = .true.
          iboun_sb(4,3) = .true.
          iboun_sb(4,6) = .true.
          iboun_sb(5,1) = .true.
          iboun_sb(5,3) = .true.
          iboun_sb(6,2) = .true.
          iboun_sb(6,4) = .true.
          iboun_sb(6,6) = .true.
          iboun_sb(7,2) = .true.
          iboun_sb(7,3) = .true.
          iboun_sb(8,1) = .true.
          iboun_sb(8,2) = .true.
          iboun_sb(8,3) = .true.
          iboun_sb(8,4) = .true.
          iboun_sb(8,5) = .true.
      CASE (3)
          x_superbrick(1) = 1.d0 / 2.d0
          y_superbrick(1) = 1.d0
          z_superbrick(1) = 2.d0

          x_superbrick(2) = 1.d0 / 2.d0
          y_superbrick(2) = 1.d0
          z_superbrick(2) = 3.d0 / 2.d0

          x_superbrick(3) = 1.d0 / 2.d0
          y_superbrick(3) = 1.d0 / 2.d0
          z_superbrick(3) = 3.d0 / 2.d0

          x_superbrick(4) = 1.d0 / 2.d0
          y_superbrick(4) = 1.d0 / 2.d0
          z_superbrick(4) = 2.d0

          x_superbrick(5) = 1.d0
          y_superbrick(5) = 1.d0
          z_superbrick(5) = 2.d0

          x_superbrick(6) = 1.d0
          y_superbrick(6) = 1.d0
          z_superbrick(6) = 1.d0

          x_superbrick(7) = 1.d0
          y_superbrick(7) = 1.d0 / 2.d0
          z_superbrick(7) = 1.d0

          x_superbrick(8) = 1.d0
          y_superbrick(8) = 1.d0 / 2.d0
          z_superbrick(8) = 2.d0

          x_superbrick(9) = 1.d0 / 2.d0
          y_superbrick(9) = 0.d0
          z_superbrick(9) = 1.d0

          x_superbrick(10) = 1.d0 / 2.d0
          y_superbrick(10) = 0.d0
          z_superbrick(10) = 2.d0

          x_superbrick(11) = 1.d0
          y_superbrick(11) = 0.d0
          z_superbrick(11) = 1.d0 / 2.d0

          x_superbrick(12) = 1.d0
          y_superbrick(12) = 0.d0
          z_superbrick(12) = 2.d0

          x_superbrick(13) = 0.d0
          y_superbrick(13) = 1.d0
          z_superbrick(13) = 1.d0

          x_superbrick(14) = 0.d0
          y_superbrick(14) = 1.d0
          z_superbrick(14) = 1.d0 / 2.d0

          x_superbrick(15) = 0.d0
          y_superbrick(15) = 0.d0
          z_superbrick(15) = 1.d0 / 2.d0

          x_superbrick(16) = 0.d0
          y_superbrick(16) = 0.d0
          z_superbrick(16) = 1.d0

          x_superbrick(17) = 1.d0 / 2.d0
          y_superbrick(17) = 1.d0
          z_superbrick(17) = 1.d0

          x_superbrick(18) = 1.d0
          y_superbrick(18) = 1.d0
          z_superbrick(18) = 1.d0 / 2.d0

          x_superbrick(19) = 0.d0
          y_superbrick(19) = 1.d0
          z_superbrick(19) = 3.d0 / 2.d0

          x_superbrick(20) = 0.d0
          y_superbrick(20) = 1.d0
          z_superbrick(20) = 2.d0

          x_superbrick(21) = 0.d0
          y_superbrick(21) = 1.d0 / 2.d0
          z_superbrick(21) = 3.d0 / 2.d0

          x_superbrick(22) = 0.d0
          y_superbrick(22) = 1.d0 / 2.d0
          z_superbrick(22) = 2.d0

          x_superbrick(23) = 0.d0
          y_superbrick(23) = 0.d0
          z_superbrick(23) = 2.d0

          x_superbrick(24) = 0.d0
          y_superbrick(24) = 1.d0
          z_superbrick(24) = 0.d0

          x_superbrick(25) = 1.d0
          y_superbrick(25) = 1.d0
          z_superbrick(25) = 0.d0

          x_superbrick(26) = 1.d0
          y_superbrick(26) = 0.d0
          z_superbrick(26) = 0.d0

          x_superbrick(27) = 0.d0
          y_superbrick(27) = 0.d0
          z_superbrick(27) = 0.d0

          ibool_superbrick(1, 1) = 3
          ibool_superbrick(2, 1) = 7
          ibool_superbrick(3, 1) = 6
          ibool_superbrick(4, 1) = 2
          ibool_superbrick(5, 1) = 4
          ibool_superbrick(6, 1) = 8
          ibool_superbrick(7, 1) = 5
          ibool_superbrick(8, 1) = 1

          ibool_superbrick(1, 2) = 9
          ibool_superbrick(2, 2) = 11
          ibool_superbrick(3, 2) = 7
          ibool_superbrick(4, 2) = 3
          ibool_superbrick(5, 2) = 10
          ibool_superbrick(6, 2) = 12
          ibool_superbrick(7, 2) = 8
          ibool_superbrick(8, 2) = 4

          ibool_superbrick(1, 3) = 15
          ibool_superbrick(2, 3) = 11
          ibool_superbrick(3, 3) = 18
          ibool_superbrick(4, 3) = 14
          ibool_superbrick(5, 3) = 16
          ibool_superbrick(6, 3) = 9
          ibool_superbrick(7, 3) = 17
          ibool_superbrick(8, 3) = 13

          ibool_superbrick(1, 4) = 21
          ibool_superbrick(2, 4) = 3
          ibool_superbrick(3, 4) = 2
          ibool_superbrick(4, 4) = 19
          ibool_superbrick(5, 4) = 22
          ibool_superbrick(6, 4) = 4
          ibool_superbrick(7, 4) = 1
          ibool_superbrick(8, 4) = 20

          ibool_superbrick(1, 5) = 9
          ibool_superbrick(2, 5) = 11
          ibool_superbrick(3, 5) = 18
          ibool_superbrick(4, 5) = 17
          ibool_superbrick(5, 5) = 3
          ibool_superbrick(6, 5) = 7
          ibool_superbrick(7, 5) = 6
          ibool_superbrick(8, 5) = 2

          ibool_superbrick(1, 6) = 16
          ibool_superbrick(2, 6) = 9
          ibool_superbrick(3, 6) = 3
          ibool_superbrick(4, 6) = 21
          ibool_superbrick(5, 6) = 23
          ibool_superbrick(6, 6) = 10
          ibool_superbrick(7, 6) = 4
          ibool_superbrick(8, 6) = 22

          ibool_superbrick(1, 7) = 16
          ibool_superbrick(2, 7) = 9
          ibool_superbrick(3, 7) = 17
          ibool_superbrick(4, 7) = 13
          ibool_superbrick(5, 7) = 21
          ibool_superbrick(6, 7) = 3
          ibool_superbrick(7, 7) = 2
          ibool_superbrick(8, 7) = 19

          ibool_superbrick(1, 8) = 27
          ibool_superbrick(2, 8) = 26
          ibool_superbrick(3, 8) = 25
          ibool_superbrick(4, 8) = 24
          ibool_superbrick(5, 8) = 15
          ibool_superbrick(6, 8) = 11
          ibool_superbrick(7, 8) = 18
          ibool_superbrick(8, 8) = 14

          iboun_sb(:,:) = .false.
          iboun_sb(1,2) = .true.
          iboun_sb(1,4) = .true.
          iboun_sb(1,6) = .true.
          iboun_sb(2,2) = .true.
          iboun_sb(2,3) = .true.
          iboun_sb(2,6) = .true.
          iboun_sb(3,1) = .true.
          iboun_sb(3,3) = .true.
          iboun_sb(3,4) = .true.
          iboun_sb(4,1) = .true.
          iboun_sb(4,4) = .true.
          iboun_sb(4,6) = .true.
          iboun_sb(5,2) = .true.
          iboun_sb(5,4) = .true.
          iboun_sb(6,1) = .true.
          iboun_sb(6,3) = .true.
          iboun_sb(6,6) = .true.
          iboun_sb(7,1) = .true.
          iboun_sb(7,4) = .true.
          iboun_sb(8,1) = .true.
          iboun_sb(8,2) = .true.
          iboun_sb(8,3) = .true.
          iboun_sb(8,4) = .true.
          iboun_sb(8,5) = .true.
      CASE (4)
          x_superbrick(1) = 1.d0 / 2.d0
          y_superbrick(1) = 0.d0
          z_superbrick(1) = 2.d0

          x_superbrick(2) = 1.d0 / 2.d0
          y_superbrick(2) = 0.d0
          z_superbrick(2) = 3.d0 / 2.d0

          x_superbrick(3) = 1.d0 / 2.d0
          y_superbrick(3) = 1.d0 / 2.d0
          z_superbrick(3) = 3.d0 / 2.d0

          x_superbrick(4) = 1.d0 / 2.d0
          y_superbrick(4) = 1.d0 / 2.d0
          z_superbrick(4) = 2.d0

          x_superbrick(5) = 1.d0
          y_superbrick(5) = 0.d0
          z_superbrick(5) = 2.d0

          x_superbrick(6) = 1.d0
          y_superbrick(6) = 0.d0
          z_superbrick(6) = 1.d0

          x_superbrick(7) = 1.d0
          y_superbrick(7) = 1.d0 / 2.d0
          z_superbrick(7) = 1.d0

          x_superbrick(8) = 1.d0
          y_superbrick(8) = 1.d0 / 2.d0
          z_superbrick(8) = 2.d0

          x_superbrick(9) = 1.d0 / 2.d0
          y_superbrick(9) = 1.d0
          z_superbrick(9) = 1.d0

          x_superbrick(10) = 1.d0 / 2.d0
          y_superbrick(10) = 1.d0
          z_superbrick(10) = 2.d0

          x_superbrick(11) = 1.d0
          y_superbrick(11) = 1.d0
          z_superbrick(11) = 1.d0 / 2.d0

          x_superbrick(12) = 1.d0
          y_superbrick(12) = 1.d0
          z_superbrick(12) = 2.d0

          x_superbrick(13) = 0.d0
          y_superbrick(13) = 0.d0
          z_superbrick(13) = 1.d0

          x_superbrick(14) = 0.d0
          y_superbrick(14) = 0.d0
          z_superbrick(14) = 1.d0 / 2.d0

          x_superbrick(15) = 0.d0
          y_superbrick(15) = 1.d0
          z_superbrick(15) = 1.d0 / 2.d0

          x_superbrick(16) = 0.d0
          y_superbrick(16) = 1.d0
          z_superbrick(16) = 1.d0

          x_superbrick(17) = 1.d0 / 2.d0
          y_superbrick(17) = 0.d0
          z_superbrick(17) = 1.d0

          x_superbrick(18) = 1.d0
          y_superbrick(18) = 0.d0
          z_superbrick(18) = 1.d0 / 2.d0

          x_superbrick(19) = 0.d0
          y_superbrick(19) = 0.d0
          z_superbrick(19) = 3.d0 / 2.d0

          x_superbrick(20) = 0.d0
          y_superbrick(20) = 0.d0
          z_superbrick(20) = 2.d0

          x_superbrick(21) = 0.d0
          y_superbrick(21) = 1.d0 / 2.d0
          z_superbrick(21) = 3.d0 / 2.d0

          x_superbrick(22) = 0.d0
          y_superbrick(22) = 1.d0 / 2.d0
          z_superbrick(22) = 2.d0

          x_superbrick(23) = 0.d0
          y_superbrick(23) = 1.d0
          z_superbrick(23) = 2.d0

          x_superbrick(24) = 0.d0
          y_superbrick(24) = 0.d0
          z_superbrick(24) = 0.d0

          x_superbrick(25) = 1.d0
          y_superbrick(25) = 0.d0
          z_superbrick(25) = 0.d0

          x_superbrick(26) = 1.d0
          y_superbrick(26) = 1.d0
          z_superbrick(26) = 0.d0

          x_superbrick(27) = 0.d0
          y_superbrick(27) = 1.d0
          z_superbrick(27) = 0.d0

          ibool_superbrick(1, 1) = 2
          ibool_superbrick(2, 1) = 6
          ibool_superbrick(3, 1) = 7
          ibool_superbrick(4, 1) = 3
          ibool_superbrick(5, 1) = 1
          ibool_superbrick(6, 1) = 5
          ibool_superbrick(7, 1) = 8
          ibool_superbrick(8, 1) = 4

          ibool_superbrick(1, 2) = 3
          ibool_superbrick(2, 2) = 7
          ibool_superbrick(3, 2) = 11
          ibool_superbrick(4, 2) = 9
          ibool_superbrick(5, 2) = 4
          ibool_superbrick(6, 2) = 8
          ibool_superbrick(7, 2) = 12
          ibool_superbrick(8, 2) = 10

          ibool_superbrick(1, 3) = 14
          ibool_superbrick(2, 3) = 18
          ibool_superbrick(3, 3) = 11
          ibool_superbrick(4, 3) = 15
          ibool_superbrick(5, 3) = 13
          ibool_superbrick(6, 3) = 17
          ibool_superbrick(7, 3) = 9
          ibool_superbrick(8, 3) = 16

          ibool_superbrick(1, 4) = 19
          ibool_superbrick(2, 4) = 2
          ibool_superbrick(3, 4) = 3
          ibool_superbrick(4, 4) = 21
          ibool_superbrick(5, 4) = 20
          ibool_superbrick(6, 4) = 1
          ibool_superbrick(7, 4) = 4
          ibool_superbrick(8, 4) = 22

          ibool_superbrick(1, 5) = 17
          ibool_superbrick(2, 5) = 18
          ibool_superbrick(3, 5) = 11
          ibool_superbrick(4, 5) = 9
          ibool_superbrick(5, 5) = 2
          ibool_superbrick(6, 5) = 6
          ibool_superbrick(7, 5) = 7
          ibool_superbrick(8, 5) = 3

          ibool_superbrick(1, 6) = 21
          ibool_superbrick(2, 6) = 3
          ibool_superbrick(3, 6) = 9
          ibool_superbrick(4, 6) = 16
          ibool_superbrick(5, 6) = 22
          ibool_superbrick(6, 6) = 4
          ibool_superbrick(7, 6) = 10
          ibool_superbrick(8, 6) = 23

          ibool_superbrick(1, 7) = 13
          ibool_superbrick(2, 7) = 17
          ibool_superbrick(3, 7) = 9
          ibool_superbrick(4, 7) = 16
          ibool_superbrick(5, 7) = 19
          ibool_superbrick(6, 7) = 2
          ibool_superbrick(7, 7) = 3
          ibool_superbrick(8, 7) = 21

          ibool_superbrick(1, 8) = 24
          ibool_superbrick(2, 8) = 25
          ibool_superbrick(3, 8) = 26
          ibool_superbrick(4, 8) = 27
          ibool_superbrick(5, 8) = 14
          ibool_superbrick(6, 8) = 18
          ibool_superbrick(7, 8) = 11
          ibool_superbrick(8, 8) = 15

          iboun_sb(:,:) = .false.
          iboun_sb(1,2) = .true.
          iboun_sb(1,3) = .true.
          iboun_sb(1,6) = .true.
          iboun_sb(2,2) = .true.
          iboun_sb(2,4) = .true.
          iboun_sb(2,6) = .true.
          iboun_sb(3,1) = .true.
          iboun_sb(3,3) = .true.
          iboun_sb(3,4) = .true.
          iboun_sb(4,1) = .true.
          iboun_sb(4,3) = .true.
          iboun_sb(4,6) = .true.
          iboun_sb(5,2) = .true.
          iboun_sb(5,3) = .true.
          iboun_sb(6,1) = .true.
          iboun_sb(6,4) = .true.
          iboun_sb(6,6) = .true.
          iboun_sb(7,1) = .true.
          iboun_sb(7,3) = .true.
          iboun_sb(8,1) = .true.
          iboun_sb(8,2) = .true.
          iboun_sb(8,3) = .true.
          iboun_sb(8,4) = .true.
          iboun_sb(8,5) = .true.
  END SELECT

  end subroutine define_basic_doubling_brick

