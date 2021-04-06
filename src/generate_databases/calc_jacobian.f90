!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


  subroutine calc_jacobian(myrank,xix_elem,xiy_elem,xiz_elem, &
                           etax_elem,etay_elem,etaz_elem, &
                           gammax_elem,gammay_elem,gammaz_elem,jacobian_elem, &
                           xelm,yelm,zelm,dershape3D)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ZERO,MAX_STRING_LEN
  use generate_databases_par, only: NGNOD,OUTPUT_FILES

  implicit none

  integer myrank

  double precision,intent(in) :: dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD),intent(in) :: xelm,yelm,zelm

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    xix_elem,xiy_elem,xiz_elem,etax_elem,etay_elem,etaz_elem, &
    gammax_elem,gammay_elem,gammaz_elem,jacobian_elem

  ! local parameters
  integer :: i,j,k,ia
  double precision :: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision :: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision :: jacobian
  character(len=MAX_STRING_LEN) :: filename

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        xxi = ZERO
        xeta = ZERO
        xgamma = ZERO
        yxi = ZERO
        yeta = ZERO
        ygamma = ZERO
        zxi = ZERO
        zeta = ZERO
        zgamma = ZERO

        do ia = 1,NGNOD
          xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
          xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
          xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)

          yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
          yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
          ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)

          zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
          zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
          zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
        enddo

        jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
                   xeta*(yxi*zgamma-ygamma*zxi) + &
                   xgamma*(yxi*zeta-yeta*zxi)

        ! check that the Jacobian transform is invertible, i.e. that the Jacobian never becomes negative or null
        if (jacobian <= ZERO) then
          print *,'Error: rank ',myrank,'found invalid element with negative Jacobian ',jacobian
          print *,'  NGNOD = ',NGNOD
          print *,'  point (i,j,k) = ',i,j,k
          print *,'  element points x/y/z: '
          do ia = 1,NGNOD
            print *,'  ',xelm(ia),yelm(ia),zelm(ia)
          enddo
          write(filename,'(a,i6.6,a)') trim(OUTPUT_FILES)//'/error_proc',myrank,'_element_with_invalid_jacobian'
          call write_VTK_data_points_elem(NGNOD,xelm,yelm,zelm,jacobian,filename)
          print *,'  written out:',trim(filename)
          print *,'Please check your mesh...'
          call exit_MPI(myrank,'Error negative or null 3D Jacobian found')
        endif

        ! invert the relation (Fletcher p. 50 vol. 2)
        xix = (yeta*zgamma-ygamma*zeta) / jacobian
        xiy = (xgamma*zeta-xeta*zgamma) / jacobian
        xiz = (xeta*ygamma-xgamma*yeta) / jacobian
        etax = (ygamma*zxi-yxi*zgamma) / jacobian
        etay = (xxi*zgamma-xgamma*zxi) / jacobian
        etaz = (xgamma*yxi-xxi*ygamma) / jacobian
        gammax = (yxi*zeta-yeta*zxi) / jacobian
        gammay = (xeta*zxi-xxi*zeta) / jacobian
        gammaz = (xxi*yeta-xeta*yxi) / jacobian

        ! compute and store the jacobian for the solver
        jacobian = 1.d0 / (xix*(etay*gammaz-etaz*gammay) &
                        -xiy*(etax*gammaz-etaz*gammax) &
                        +xiz*(etax*gammay-etay*gammax))

        ! save the derivatives and the jacobian
        ! distinguish between single and double precision for reals
        xix_elem(i,j,k) = real(xix,kind=CUSTOM_REAL)
        xiy_elem(i,j,k) = real(xiy,kind=CUSTOM_REAL)
        xiz_elem(i,j,k) = real(xiz,kind=CUSTOM_REAL)
        etax_elem(i,j,k) = real(etax,kind=CUSTOM_REAL)
        etay_elem(i,j,k) = real(etay,kind=CUSTOM_REAL)
        etaz_elem(i,j,k) = real(etaz,kind=CUSTOM_REAL)
        gammax_elem(i,j,k) = real(gammax,kind=CUSTOM_REAL)
        gammay_elem(i,j,k) = real(gammay,kind=CUSTOM_REAL)
        gammaz_elem(i,j,k) = real(gammaz,kind=CUSTOM_REAL)
        jacobian_elem(i,j,k) = real(jacobian,kind=CUSTOM_REAL)

      enddo
    enddo
  enddo

  end subroutine calc_jacobian

!
!-------------------------------------------------------------------------------------------------
!

  subroutine calc_coords(x_elem,y_elem,z_elem,xelm,yelm,zelm,shape3D)

  use generate_databases_par, only: NGNOD,NGLLX,NGLLY,NGLLZ,ZERO

  implicit none

  double precision :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: x_elem,y_elem,z_elem

  !local
  integer :: i,j,k,ia
  double precision :: xmesh,ymesh,zmesh

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        xmesh = ZERO
        ymesh = ZERO
        zmesh = ZERO

        do ia = 1,NGNOD
          xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
          ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
          zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
        enddo

        x_elem(i,j,k) = xmesh
        y_elem(i,j,k) = ymesh
        z_elem(i,j,k) = zmesh

      enddo
    enddo
  enddo

  end subroutine calc_coords

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_element_regularity(xelm,yelm,zelm,any_regular_elem,cube_edge_size_squared, &
                                      nspec_irregular,ispec,nspec,irregular_element_number,ANY_FAULT_IN_THIS_PROC)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,myrank
  use generate_databases_par, only: NGNOD,USE_MESH_COLORING_GPU
  use create_regions_mesh_ext_par, only: dershape3D

  implicit none
  real, dimension(NGNOD),intent(in) :: xelm,yelm,zelm

  logical, intent(inout) :: any_regular_elem
  double precision, intent(inout) :: cube_edge_size_squared
  integer, intent(inout) :: nspec_irregular

  integer, intent(in) :: ispec,nspec
  integer, dimension(nspec),intent(inout) :: irregular_element_number
  logical, intent(in) :: ANY_FAULT_IN_THIS_PROC

  ! local parameters
  double precision :: dist1_sq,dist2_sq,dist3_sq
  double precision :: threshold
  double precision,parameter :: threshold_percentage = 1.e-5
  double precision,parameter :: threshold_zero = 1.e-25

  logical :: eqx1,eqx2,eqx3,eqx4,eqx5,eqx6
  logical :: eqy1,eqy2,eqy3,eqy4,eqy5,eqy6
  logical :: eqz1,eqz2,eqz3,eqz4,eqz5,eqz6
  logical :: is_regular_element
  logical, external :: is_equal_number

  ! jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix_reg,xiy_reg,xiz_reg,etax_reg,etay_reg,etaz_reg, &
                                                          gammax_reg,gammay_reg,gammaz_reg,jacobian_reg
  double precision, dimension(NGNOD) :: xelm_dble,yelm_dble,zelm_dble

  ! by default, we assume to have a perfect regular shape (cube)

  ! checks if the potential cube has the same size as the previous ones
  dist1_sq = (xelm(2)-xelm(1))**2 + (yelm(2)-yelm(1))**2 +(zelm(2)-zelm(1))**2

  ! sets irregular element (default for NGNOD 27 and others)
  threshold = threshold_percentage * cube_edge_size_squared
  if (NGNOD == 27 &
      .or. ANY_FAULT_IN_THIS_PROC &
      .or. USE_MESH_COLORING_GPU &
      .or. (any_regular_elem .and. ( abs(dist1_sq - cube_edge_size_squared) > threshold ))) then
    ! not a regular shape
    irregular_element_number(ispec) = ispec - (nspec - nspec_irregular)
    return
  endif

  ! number convention:
  !
  !              8                 7
  !              o ------------- o
  !             /|              /|
  !            / |             / |
  !          5o---------------o6 |
  !           |  |            |  |
  !           | 4o----------- |--o3
  !           | /             | /
  !           |/              |/
  !     z     o---------------o
  !     | y   1                 2
  !     |/
  !     o---> x
  !

  ! checks x-plane (using epsilon() function to determine almost negligible differences)
  eqx1 = is_equal_number(xelm(1),xelm(4))
  eqx2 = is_equal_number(xelm(1),xelm(5))
  eqx3 = is_equal_number(xelm(1),xelm(8))
  eqx4 = is_equal_number(xelm(2),xelm(3))
  eqx5 = is_equal_number(xelm(2),xelm(6))
  eqx6 = is_equal_number(xelm(2),xelm(7))
  ! checks y-plane
  eqy1 = is_equal_number(yelm(1),yelm(2))
  eqy2 = is_equal_number(yelm(1),yelm(5))
  eqy3 = is_equal_number(yelm(1),yelm(6))
  eqy4 = is_equal_number(yelm(3),yelm(4))
  eqy5 = is_equal_number(yelm(3),yelm(7))
  eqy6 = is_equal_number(yelm(3),yelm(8))
  ! checks z-plane
  eqz1 = is_equal_number(zelm(1),zelm(2))
  eqz2 = is_equal_number(zelm(1),zelm(3))
  eqz3 = is_equal_number(zelm(1),zelm(4))
  eqz4 = is_equal_number(zelm(5),zelm(6))
  eqz5 = is_equal_number(zelm(5),zelm(7))
  eqz6 = is_equal_number(zelm(5),zelm(8))

  ! checks if the element is a cube (following numbering convention in a 8 nodes element)
  ! with numerically negligible differences (to avoid representation issues of float values)
  if (eqx1 .and. eqx2 .and. eqx3 .and. eqx4 .and. eqx5 .and. eqx6 .and. &
      eqy1 .and. eqy2 .and. eqy3 .and. eqy4 .and. eqy5 .and. eqy6 .and. &
      eqz1 .and. eqz2 .and. eqz3 .and. eqz4 .and. eqz5 .and. eqz6 ) then
  ! or direct float comparisons
  !if (xelm(1) == xelm(4) .and. xelm(1) == xelm(5) .and. xelm(1) == xelm(8) .and. &
  !    xelm(2) == xelm(3) .and. xelm(2) == xelm(6) .and. xelm(2) == xelm(7) .and. &
  !    yelm(1) == yelm(2) .and. yelm(1) == yelm(5) .and. yelm(1) == yelm(6) .and. &
  !    yelm(3) == yelm(4) .and. yelm(3) == yelm(7) .and. yelm(3) == yelm(8) .and. &
  !    zelm(1) == zelm(2) .and. zelm(1) == zelm(3) .and. zelm(1) == zelm(4) .and. &
  !    zelm(5) == zelm(6) .and. zelm(5) == zelm(7) .and. zelm(5) == zelm(8) ) then

    dist2_sq = (xelm(5)-xelm(1))**2 + (yelm(5)-yelm(1))**2 + (zelm(5)-zelm(1))**2
    dist3_sq = (xelm(4)-xelm(1))**2 + (yelm(4)-yelm(1))**2 + (zelm(4)-zelm(1))**2

    threshold = threshold_percentage * dist1_sq

    is_regular_element = .false.

    ! checks derivatives of shape function
    if (abs(dist2_sq - dist1_sq) < threshold .and. abs(dist3_sq - dist1_sq) < threshold) then
      ! regular shape
      ! checks shape functions and orientation
      xelm_dble(:) = dble(xelm(:))
      yelm_dble(:) = dble(yelm(:))
      zelm_dble(:) = dble(zelm(:))

      ! jacobian and derivatives of mapping
      call calc_jacobian(myrank,xix_reg,xiy_reg,xiz_reg, &
                         etax_reg,etay_reg,etaz_reg, &
                         gammax_reg,gammay_reg,gammaz_reg, &
                         jacobian_reg,xelm_dble,yelm_dble,zelm_dble,dershape3D)

      ! only xix == etay == gammaz are non-zero for regular elements
      ! check
      if ((abs(xix_reg(1,1,1) - etay_reg(1,1,1)) < threshold_zero) .and. &
          (abs(xix_reg(1,1,1) - gammaz_reg(1,1,1)) < threshold_zero) .and. &
           abs(xiy_reg(1,1,1)) < threshold_zero .and. abs(xiz_reg(1,1,1)) < threshold_zero .and. &
           abs(etax_reg(1,1,1)) < threshold_zero .and. abs(etaz_reg(1,1,1)) < threshold_zero .and. &
           abs(gammax_reg(1,1,1)) < threshold_zero .and. abs(gammay_reg(1,1,1)) < threshold_zero) then
        ! regular shape
        is_regular_element = .true.
      else
        ! debug
        !print *,'debug: regular element should have xix == etay == gammaz ',xix_reg(1,1,1),etay_reg(1,1,1),gammaz_reg(1,1,1)
        !print *,'  xix    ',xix_reg(1,1,1),xiy_reg(1,1,1),xiz_reg(1,1,1)
        !print *,'  etax   ',etax_reg(1,1,1),etay_reg(1,1,1),etaz_reg(1,1,1)
        !print *,'  gammax ',gammax_reg(1,1,1),gammay_reg(1,1,1),gammaz_reg(1,1,1)
        ! non-regular (might be due to difference in orientation)
        is_regular_element = .false.
      endif
    endif

    if (is_regular_element) then
      ! regular shape
      irregular_element_number(ispec) = 0
      ! test if first cube found in mesh
      if (.not. any_regular_elem ) then
        cube_edge_size_squared = dist1_sq
        any_regular_elem = .true.
      endif
      ! regular shape (perfect cube), decreases number of irregular elements
      nspec_irregular = nspec_irregular - 1
    else
      ! not a regular shape
      irregular_element_number(ispec) = ispec - (nspec - nspec_irregular)
    endif

  else
    ! not a regular shape
    irregular_element_number(ispec) = ispec - (nspec - nspec_irregular)
  endif

  end subroutine check_element_regularity

!
!-------------------------------------------------------------------------------------------------
!

  function is_equal_number(a,b)

! compares two real values and determines if they are equal (within numerical precision)
! to avoid direct comparisons a == b which can fail for rounding issues

  implicit none

  real, intent(in) :: a,b
  logical :: is_equal_number

  ! note: epsilon() intrinsic function
  !       returns a positive number that is almost negligible
  !       example: If x is of type REAL(4), EPSILON (X) has the value 2**-23
  if (abs(a - b) <= epsilon(a)) then
    is_equal_number = .true.
  else
    is_equal_number = .false.
  endif

  end function is_equal_number
