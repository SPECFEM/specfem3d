!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

subroutine compute_boundary_kernel(kernel, mul, kappal, rho_vsl, accel, b_displ, ds, b_ds, norm)

! compute the boundary kernel contribution from one side of the boundary

  implicit none
  include 'constants.h'
 
  real(kind=CUSTOM_REAL)  kernel, mul, kappal, rho_vsl
  real(kind=CUSTOM_REAL) :: accel(NDIM), b_displ(NDIM), ds(NDIM,NDIM), b_ds(NDIM,NDIM), norm(NDIM)

  real(kind=CUSTOM_REAL) :: eps3, eps(NDIM,NDIM), epsdev(NDIM,NDIM), normal(NDIM,1)
  real(kind=CUSTOM_REAL) :: b_eps3, b_eps(NDIM,NDIM), b_epsdev(NDIM,NDIM)
  real(kind=CUSTOM_REAL) :: temp1(NDIM,NDIM), rhol, kl(1,1), one_matrix(1,1)
   

  normal(:,1) = norm
  one_matrix(1,1) = ONE

  eps3 = ds(1,1) + ds(2,2) + ds(3,3)
  
  eps(1,1) = ds(1,1)
  eps(2,2) = ds(2,2)
  eps(3,3) = ds(3,3)
  eps(1,2) = (ds(1,2) + ds(2,1))/2
  eps(1,3) = (ds(1,3) + ds(3,1))/2
  eps(2,3) = (ds(2,3) + ds(3,2))/2
  eps(2,1) = eps(1,2)
  eps(3,1) = eps(1,3)
  eps(3,2) = eps(2,3)

  epsdev = eps
  epsdev(1,1) = eps(1,1) - eps3 / 3
  epsdev(2,2) = eps(2,2) - eps3 / 3
  epsdev(3,3) = eps(3,3) - eps3 / 3
 
 
  b_eps3 = b_ds(1,1) + b_ds(2,2) + b_ds(3,3)
  
  b_eps(1,1) = b_ds(1,1)
  b_eps(2,2) = b_ds(2,2)
  b_eps(3,3) = b_ds(3,3)
  b_eps(1,2) = (b_ds(1,2) + b_ds(2,1))/2
  b_eps(1,3) = (b_ds(1,3) + b_ds(3,1))/2
  b_eps(2,3) = (b_ds(2,3) + b_ds(3,2))/2
  b_eps(2,1) = b_eps(1,2)
  b_eps(3,1) = b_eps(1,3)
  b_eps(3,2) = b_eps(2,3)

  b_epsdev = b_eps
  b_epsdev(1,1) = b_eps(1,1) - b_eps3 / 3
  b_epsdev(2,2) = b_eps(2,2) - b_eps3 / 3
  b_epsdev(3,3) = b_eps(3,3) - b_eps3 / 3

  temp1 = matmul(epsdev,b_epsdev)

  rhol = rho_vsl ** 2 / mul

  kl = (rhol * dot_product(accel(:), b_displ(:)) &
             + kappal * eps3 * b_eps3 &
             + 2 * mul * (temp1(1,1) + temp1(2,2) + temp1(3,3))) * one_matrix &
             - kappal *  matmul(transpose(normal),matmul(eps,normal)) * b_eps3 &
             - kappal *  matmul(transpose(normal),matmul(b_eps,normal)) * eps3 &
             - 2 * mul * matmul(transpose(normal), matmul(matmul(b_epsdev,ds), normal)) &
             - 2 * mul * matmul(transpose(normal), matmul(matmul(epsdev,b_ds), normal))
             
  kernel = kl(1,1)

end subroutine compute_boundary_kernel
