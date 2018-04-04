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

  use generate_databases_par, only: NGNOD,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,SIZE_REAL,ZERO

  implicit none

  integer myrank

  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    xix_elem,xiy_elem,xiz_elem,etax_elem,etay_elem,etaz_elem, &
    gammax_elem,gammay_elem,gammaz_elem,jacobian_elem


  integer i,j,k,ia
  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision jacobian

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO

      do ia=1,NGNOD
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
      if (jacobian <= ZERO) call exit_MPI(myrank,'Error negative or null 3D Jacobian found')

!     invert the relation (Fletcher p. 50 vol. 2)
      xix = (yeta*zgamma-ygamma*zeta) / jacobian
      xiy = (xgamma*zeta-xeta*zgamma) / jacobian
      xiz = (xeta*ygamma-xgamma*yeta) / jacobian
      etax = (ygamma*zxi-yxi*zgamma) / jacobian
      etay = (xxi*zgamma-xgamma*zxi) / jacobian
      etaz = (xgamma*yxi-xxi*ygamma) / jacobian
      gammax = (yxi*zeta-yeta*zxi) / jacobian
      gammay = (xeta*zxi-xxi*zeta) / jacobian
      gammaz = (xxi*yeta-xeta*yxi) / jacobian

!     compute and store the jacobian for the solver
      jacobian = 1. / (xix*(etay*gammaz-etaz*gammay) &
                      -xiy*(etax*gammaz-etaz*gammax) &
                      +xiz*(etax*gammay-etay*gammax))

!     save the derivatives and the jacobian

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



  subroutine calc_coords(x_elem,y_elem,z_elem,xelm,yelm,zelm,shape3D)

  use generate_databases_par, only: NGNOD,NGLLX,NGLLY,NGLLZ,ZERO

  implicit none

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: x_elem,y_elem,z_elem

  !local
  integer i,j,k,ia
  double precision xmesh,ymesh,zmesh

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      xmesh = ZERO
      ymesh = ZERO
      zmesh = ZERO

      do ia=1,NGNOD
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


  subroutine check_element_regularity(xelm,yelm,zelm,any_regular_elem,cube_edge_size_squared, &
                                      nspec_irregular,ispec,nspec,irregular_element_number,ANY_FAULT_IN_THIS_PROC)

  use generate_databases_par, only: NGNOD,CUSTOM_REAL,USE_MESH_COLORING_GPU

  real, dimension(NGNOD) :: xelm,yelm,zelm
  logical                            :: any_regular_elem,ANY_FAULT_IN_THIS_PROC
  real                               :: cube_edge_size_squared
  integer                            :: nspec_irregular,ispec,nspec
  integer, dimension(nspec)          :: irregular_element_number

  !local
  real                               :: dist1_sq,dist2_sq,dist3_sq

  !checks if the potential cube has the same size as the previous ones
  dist1_sq = (xelm(2)-xelm(1))**2 + (yelm(2)-yelm(1))**2 +(zelm(2)-zelm(1))**2
  if (NGNOD == 27 .or. ANY_FAULT_IN_THIS_PROC .or. USE_MESH_COLORING_GPU .or. &
     (any_regular_elem .and. ( abs(dist1_sq - cube_edge_size_squared) > (1e-5)*cube_edge_size_squared ))) then
    irregular_element_number(ispec) = ispec - (nspec - nspec_irregular)
    return
  endif

  ! checks if the element is a cube (following numbering convention in a 8 nodes element)
  if (xelm(1) == xelm(4) .and. xelm(1) == xelm(5) .and. xelm(1) == xelm(8) .and. &
      xelm(2) == xelm(3) .and. xelm(2) == xelm(6) .and. xelm(2) == xelm(7) .and. &
      yelm(1) == yelm(2) .and. yelm(1) == yelm(5) .and. yelm(1) == yelm(6) .and. &
      yelm(3) == yelm(4) .and. yelm(3) == yelm(7) .and. yelm(3) == yelm(8) .and. &
      zelm(1) == zelm(2) .and. zelm(1) == zelm(3) .and. zelm(1) == zelm(4) .and. &
      zelm(5) == zelm(6) .and. zelm(5) == zelm(7) .and. zelm(5) == zelm(8) ) then

    dist2_sq = (xelm(5)-xelm(1))**2 + (yelm(5)-yelm(1))**2 +(zelm(5)-zelm(1))**2
    dist3_sq = (xelm(4)-xelm(1))**2 + (yelm(4)-yelm(1))**2 +(zelm(4)-zelm(1))**2

    if (abs(dist2_sq - dist1_sq) < 1e-5*dist1_sq .and. abs(dist3_sq - dist1_sq) < 1e-5*dist1_sq) then

      ! test if first cube found in mesh
      if (.not. any_regular_elem ) then
        cube_edge_size_squared = dist1_sq
        any_regular_elem = .true.
      endif
      nspec_irregular = nspec_irregular - 1
    else
      irregular_element_number(ispec) = ispec - (nspec - nspec_irregular)
    endif

  else
    irregular_element_number(ispec) = ispec - (nspec - nspec_irregular)
  endif

  end subroutine check_element_regularity
