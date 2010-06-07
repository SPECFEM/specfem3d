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

subroutine compute_boundary_kernel()


! isotropic topography kernel computation
! compare with Tromp et al. (2005), eq. (25), or see Liu & Tromp (2008), eq. (65)        

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL):: kernel_moho_top,kernel_moho_bot
  integer :: i,j,k
  integer :: ispec2D,igll,jgll
  integer :: ispec_top,ispec_bot,iglob_top,iglob_bot
  logical :: is_done
      
  ! loops over top/bottom elements of moho surface
  do ispec2D = 1, NSPEC2D_MOHO
    ispec_top = ibelm_moho_top(ispec2D)
    ispec_bot = ibelm_moho_bot(ispec2D)

    ! elements on both sides available  
    if( ispec_top > 0 .and. ispec_bot > 0 ) then
      ! loops over surface
      do igll=1,NGLLSQUARE
        i = ijk_moho_top(1,igll,ispec2D)
        j = ijk_moho_top(2,igll,ispec2D)
        k = ijk_moho_top(3,igll,ispec2D)            
        iglob_top = ibool(i,j,k,ispec_top)

        ! computes contribution from top element
        call compute_boundary_kernel_elem( kernel_moho_top, &
                    mustore(i,j,k,ispec_top), &
                    kappastore(i,j,k,ispec_top),rho_vs(i,j,k,ispec_top), &
                    accel(:,iglob_top),b_displ(:,iglob_top), &
                    dsdx_top(:,:,i,j,k,ispec2D),b_dsdx_top(:,:,i,j,k,ispec2D), &
                    normal_moho_top(:,igll,ispec2D) )

        ! finds corresponding global node in bottom element
        is_done = .false.
        do jgll = 1,NGLLSQUARE
          i = ijk_moho_bot(1,jgll,ispec2D)
          j = ijk_moho_bot(2,jgll,ispec2D)
          k = ijk_moho_bot(3,jgll,ispec2D)
          iglob_bot = ibool(i,j,k,ispec_bot)
        
          if( iglob_bot /= iglob_top ) cycle
          ! iglob_top == iglob_bot!

          ! computes contribution from bottom element
          call compute_boundary_kernel_elem( kernel_moho_bot, &
                      mustore(i,j,k,ispec_bot), &
                      kappastore(i,j,k,ispec_bot),rho_vs(i,j,k,ispec_bot), &
                      accel(:,iglob_bot),b_displ(:,iglob_bot), &
                      dsdx_bot(:,:,i,j,k,ispec2D),b_dsdx_bot(:,:,i,j,k,ispec2D), &
                      normal_moho_bot(:,jgll,ispec2D) )

          ! note: kernel point position: indices given by ijk_moho_top(:,igll,ispec2D)            
          moho_kl(igll,ispec2D) = moho_kl(igll,ispec2D) &
                                + (kernel_moho_top - kernel_moho_bot) * deltat
          
          ! kernel done for this point
          is_done = .true.
        enddo
        
        ! checks
        if( .not. is_done ) then
          print*,'error : moho kernel not computed'
          print*,'ispec:',ispec_top,ispec_bot,iglob_top,i,j,k
          call exit_mpi(myrank,'error moho kernel computation')
        endif
        
      enddo

    ! only one element available
    ! e.g. free-surface: see Tromp et al. (2005), eq. (28)
    else if( ispec_bot > 0 .or. ispec_top > 0 ) then

      ! loops over surface
      do igll=1,NGLLSQUARE

        if( ispec_top > 0 ) then
          i = ijk_moho_top(1,igll,ispec2D)
          j = ijk_moho_top(2,igll,ispec2D)
          k = ijk_moho_top(3,igll,ispec2D)            
          iglob_top = ibool(i,j,k,ispec_top)
          
          ! computes contribution from top element          
          call compute_boundary_kernel_elem( kernel_moho_top, &
                    mustore(i,j,k,ispec_top), &
                    kappastore(i,j,k,ispec_top),rho_vs(i,j,k,ispec_top), &
                    accel(:,iglob_top),b_displ(:,iglob_top), &
                    dsdx_top(:,:,i,j,k,ispec2D),b_dsdx_top(:,:,i,j,k,ispec2D), &
                    normal_moho_top(:,igll,ispec2D) )

          ! note: kernel point position igll: indices given by ijk_moho_top(:,igll,ispec2D)            
          moho_kl(igll,ispec2D) = moho_kl(igll,ispec2D) + kernel_moho_top * deltat

        else
          i = ijk_moho_bot(1,igll,ispec2D)
          j = ijk_moho_bot(2,igll,ispec2D)
          k = ijk_moho_bot(3,igll,ispec2D)            
          iglob_bot = ibool(i,j,k,ispec_bot)
          
          ! computes contribution from bottom element          
          call compute_boundary_kernel_elem( kernel_moho_bot, &
                    mustore(i,j,k,ispec_bot), &
                    kappastore(i,j,k,ispec_bot),rho_vs(i,j,k,ispec_bot), &
                    accel(:,iglob_bot),b_displ(:,iglob_bot), &
                    dsdx_bot(:,:,i,j,k,ispec2D),b_dsdx_bot(:,:,i,j,k,ispec2D), &
                    normal_moho_bot(:,igll,ispec2D) )
                    
          ! note: kernel point position igll: indices given by ijk_moho_bot(:,igll,ispec2D)            
          moho_kl(igll,ispec2D) = moho_kl(igll,ispec2D) - kernel_moho_bot * deltat
          
        endif
      enddo
    endif          
  enddo ! ispec2D      


end subroutine compute_boundary_kernel

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_boundary_kernel_elem(kernel, mul, kappal, rho_vsl, &
                                        accel, b_displ, ds, b_ds, norm)

! compute the boundary kernel contribution from one side of the boundary
! see e.g.: Tromp et al. (2005), eq. (25), or Liu & Tromp (2008), eq. (65)

  implicit none
  include 'constants.h'

  real(kind=CUSTOM_REAL)  kernel, mul, kappal, rho_vsl
  real(kind=CUSTOM_REAL) :: accel(NDIM), b_displ(NDIM), ds(NDIM,NDIM), b_ds(NDIM,NDIM), norm(NDIM)

  real(kind=CUSTOM_REAL) :: eps3, eps(NDIM,NDIM), epsdev(NDIM,NDIM), normal(NDIM,1)
  real(kind=CUSTOM_REAL) :: b_eps3, b_eps(NDIM,NDIM), b_epsdev(NDIM,NDIM)
  real(kind=CUSTOM_REAL) :: temp1(NDIM,NDIM), rhol, kl(1,1), one_matrix(1,1)


  normal(:,1) = norm
  one_matrix(1,1) = ONE

  ! adjoint strain (epsilon) trace
  eps3 = ds(1,1) + ds(2,2) + ds(3,3)

  ! adjoint strain tensor
  eps(1,1) = ds(1,1)
  eps(2,2) = ds(2,2)
  eps(3,3) = ds(3,3)
  eps(1,2) = (ds(1,2) + ds(2,1))/2
  eps(1,3) = (ds(1,3) + ds(3,1))/2
  eps(2,3) = (ds(2,3) + ds(3,2))/2
  eps(2,1) = eps(1,2)
  eps(3,1) = eps(1,3)
  eps(3,2) = eps(2,3)

  ! adjoint deviatoric strain component
  epsdev = eps
  epsdev(1,1) = eps(1,1) - eps3 / 3
  epsdev(2,2) = eps(2,2) - eps3 / 3
  epsdev(3,3) = eps(3,3) - eps3 / 3


  ! backward/reconstructed-forward strain (epsilon) trace
  b_eps3 = b_ds(1,1) + b_ds(2,2) + b_ds(3,3)

  ! backward/reconstructed-forward strain tensor
  b_eps(1,1) = b_ds(1,1)
  b_eps(2,2) = b_ds(2,2)
  b_eps(3,3) = b_ds(3,3)
  b_eps(1,2) = (b_ds(1,2) + b_ds(2,1))/2
  b_eps(1,3) = (b_ds(1,3) + b_ds(3,1))/2
  b_eps(2,3) = (b_ds(2,3) + b_ds(3,2))/2
  b_eps(2,1) = b_eps(1,2)
  b_eps(3,1) = b_eps(1,3)
  b_eps(3,2) = b_eps(2,3)

  ! backward/reconstructed-forward deviatoric strain
  b_epsdev = b_eps
  b_epsdev(1,1) = b_eps(1,1) - b_eps3 / 3
  b_epsdev(2,2) = b_eps(2,2) - b_eps3 / 3
  b_epsdev(3,3) = b_eps(3,3) - b_eps3 / 3

  ! matrix multiplication
  temp1 = matmul(epsdev,b_epsdev)

  ! density value
  rhol = rho_vsl ** 2 / mul

  ! isotropic kernel value 
  ! see e.g.: Tromp et al. (2005), eq. (25), or Liu & Tromp 2008, eq. (65)
  kl = ( rhol * dot_product(accel(:), b_displ(:)) + kappal * eps3 * b_eps3 &
       + 2 * mul * (temp1(1,1) + temp1(2,2) + temp1(3,3)) ) * one_matrix &
       - kappal *  matmul(transpose(normal),matmul(eps,normal)) * b_eps3 &
       - kappal *  matmul(transpose(normal),matmul(b_eps,normal)) * eps3 &
       - 2 * mul * matmul(transpose(normal), matmul(matmul(b_epsdev,ds), normal)) &
       - 2 * mul * matmul(transpose(normal), matmul(matmul(epsdev,b_ds), normal))

  kernel = kl(1,1)

end subroutine compute_boundary_kernel_elem
