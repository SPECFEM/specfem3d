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
!> Routines for general matrix-matrix and matrix-vector multiplication. Called a
!! bazillion times, presumably fast.
module unrolled_loops

  use global_parameters, only: realkind

  implicit none
  public

  contains

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! Size is fixed to npol x npol
pure subroutine mxm(a,b,c)

  use data_mesh, only: npol
  use global_parameters, only: realkind

  real(kind=realkind), intent(in)  :: a(0: ,0: ),b(0: ,0: ) ! < Input matrices
  real(kind=realkind), intent(out) :: c(0: ,0: )            ! < Result
  integer                          :: i, j

  do j = 0, npol
     do i = 0, npol
        c(i,j) = sum(a(i,:) * b(:,j))
     enddo
  enddo

end subroutine mxm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies vector a leftwise to matrix b to have vector c.
!! Size is fixed to npol x npol
pure subroutine vxm(a,b,c)

  use data_mesh, only: npol
  use global_parameters, only: realkind

  real(kind=realkind), intent(in)  :: a(0: )         ! < Vector a
  real(kind=realkind), intent(in)  :: b(0: ,0: )     ! < Matrix b
  real(kind=realkind), intent(out) :: c(0: )         ! < Resulting vector c
  integer                          :: j

  do j = 0, npol
     c(j) = sum(a * b(:,j))
  enddo

end subroutine vxm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine mxm_cg4_sparse_a(a,b,c)
   ! mxm for sparse a as found for coarse grained memory variables cg4

   real(kind=realkind), intent(in)  :: a(1:4), b(0:4,0:4)
   real(kind=realkind), intent(out) :: c(0:4,0:4)
   integer                          :: j

   ! c ist sparse, so initialization does matter
   c = 0

   do j = 0, 4
     c(1,j) = &
        + a(1) * b(1,j) &
        + a(2) * b(3,j)

     c(3,j) = &
        + a(3) * b(1,j) &
        + a(4) * b(3,j)
   enddo

end subroutine mxm_cg4_sparse_a
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine mxm_cg4_sparse_b(a,b,c)
   ! mxm for sparse b as found for coarse grained memory variables cg4

   real(kind=realkind), intent(in)  :: a(0:4,0:4), b(1:4)
   real(kind=realkind), intent(out) :: c(0:4,0:4)
   integer i

   ! c ist sparse, so initialization does matter
   c = 0

   do i = 0, 4
     c(i,1) = &
        + a(i,1) * b(1) &
        + a(i,3) * b(3)

     c(i,3) = &
        + a(i,1) * b(2) &
        + a(i,3) * b(4)
   enddo

end subroutine mxm_cg4_sparse_b
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine mxm_cg4_sparse_c(a,b,c)


  real(kind=realkind), intent(in)  :: a(0:,0:), b(0:,0:)
  real(kind=realkind), intent(out) :: c(1:4)

  c(1) = &
     + a(1,0) * b(0,1) &
     + a(1,1) * b(1,1) &
     + a(1,2) * b(2,1) &
     + a(1,3) * b(3,1) &
     + a(1,4) * b(4,1)

  c(2) = &
     + a(1,0) * b(0,3) &
     + a(1,1) * b(1,3) &
     + a(1,2) * b(2,3) &
     + a(1,3) * b(3,3) &
     + a(1,4) * b(4,3)

  c(3) = &
     + a(3,0) * b(0,1) &
     + a(3,1) * b(1,1) &
     + a(3,2) * b(2,1) &
     + a(3,3) * b(3,1) &
     + a(3,4) * b(4,1)

  c(4) = &
     + a(3,0) * b(0,3) &
     + a(3,1) * b(1,3) &
     + a(3,2) * b(2,3) &
     + a(3,3) * b(3,3) &
     + a(3,4) * b(4,3)

end subroutine mxm_cg4_sparse_c
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! Size is fixed to 4x4
pure subroutine mxm_4(a,b,c)

  use global_parameters, only: realkind

  real(kind=realkind), intent(in)  :: a(0:4,0:4),b(0:4,0:4) ! < Input matrices
  real(kind=realkind), intent(out) :: c(0:4,0:4)            ! < Result
  integer                          :: i

  do i = 0, 4
     c(i,0) = sum(a(i,:) * b(:,0))
  enddo
  do i = 0, 4
     c(i,1) = sum(a(i,:) * b(:,1))
  enddo
  do i = 0, 4
     c(i,2) = sum(a(i,:) * b(:,2))
  enddo
  do i = 0, 4
     c(i,3) = sum(a(i,:) * b(:,3))
  enddo
  do i = 0, 4
     c(i,4) = sum(a(i,:) * b(:,4))
  enddo

end subroutine mxm_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies vector a leftwise to matrix b to have vector c.
!! Size is fixed to npol x npol
pure subroutine vxm_4(a,b,c)

  use global_parameters, only: realkind

  real(kind=realkind), intent(in)  :: a(0:4)         ! < Vector a
  real(kind=realkind), intent(in)  :: b(0:4,0:4)     ! < Matrix b
  real(kind=realkind), intent(out) :: c(0:4)         ! < Resulting vector c
  integer                          :: j

  do j = 0, 4
     c(j) = sum(a * b(:,j))
  enddo

end subroutine vxm_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function outerprod(a,b)
  ! outer product (dyadic) from numerical recipes

  real(kind=realkind), dimension(:), intent(in)     :: a, b
  real(kind=realkind), dimension(size(a),size(b))   :: outerprod

  outerprod = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))
end function outerprod
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function outerprod_4(a,b)
  ! outer product (dyadic) from numerical recipes

  real(kind=realkind), dimension(0:4), intent(in)   :: a, b
  real(kind=realkind), dimension(0:4,0:4)           :: outerprod_4

  outerprod_4 = spread(a, dim=2, ncopies=5) * spread(b, dim=1, ncopies=5)
end function outerprod_4
!-----------------------------------------------------------------------------------------


end module unrolled_loops
!=========================================================================================
