module elastic_tensor_tools_mod
  
  use interpolation_mod, only: si, sp, di, dp, cp, hp, deg2rad, rad2deg

  implicit none

  integer(kind=si), dimension(:,:), allocatable :: ind_vec2tens, ind_vec2tens_voigt
  
contains

  !================================================================================
  ! Define rotation matrix with euler angles :
  !    a = angle from   x-axis (rotation around   z-axis), (xyz       ->    x'y'z') 
  !    b = angle from  z'-axis (rotation around  y'-axis), (x'y'z'    -> x''y''z'')
  !    c = angle from y''-axis (rotation around x''-axis), (x''y''z'' ->       XYZ)
  ! Rotmat = matrix from ref xyz (e.g. cartesian or vti) to rotated new one
  subroutine define_rotation_matrix(a0,b0,c0,rotmat,rotmat_t)

    real(kind=dp),    intent(in) :: a0, b0, c0
    real(kind=dp)                :: a, b, c
    real(kind=dp), dimension(3,3), intent(out) :: rotmat, rotmat_t

    !* Pass in radian
    a = deg2rad * a0
    b = deg2rad * b0
    c = deg2rad * c0

    !* First column
    rotmat(1,1) =  cos(a) * sin(b)
    rotmat(2,1) = -sin(a) * cos(c) + cos(a) * cos(b) * sin(c)
    rotmat(3,1) =  sin(a) * sin(c) + cos(a) * cos(b) * cos(c)

    !* Second column
    rotmat(1,2) =  sin(a) * sin(b)
    rotmat(2,2) =  cos(a) * cos(c) + sin(a) * cos(b) * sin(c)
    rotmat(3,2) = -cos(a) * sin(c) + sin(a) * cos(b) * cos(c)

    !* Third column
    rotmat(1,3) =  cos(b)
    rotmat(2,3) = -sin(b) * sin(c)   
    rotmat(3,3) = -sin(b) * cos(c)

    !* Define reverse rotation
    rotmat_t    = transpose(rotmat)

  end subroutine define_rotation_matrix
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define bond rotation matrices to rotate a voigt tensor
  !    (see e.g. Auld 1973, for bond matrix definition secion D pages 73-76)
  ! HERE THE STRESS ONE => can be used to rotate the stress vector or stiffness tensor
  subroutine define_bond_stress_matrix(rotmat,bond)

    real(kind=dp), dimension(3,3), intent(in)  :: rotmat
    real(kind=dp), dimension(6,6), intent(out) :: bond  

    ! First column
    bond(1,1) = rotmat(1,1)*rotmat(1,1)
    bond(2,1) = rotmat(2,1)*rotmat(2,1)
    bond(3,1) = rotmat(3,1)*rotmat(3,1)
    bond(4,1) = rotmat(2,1)*rotmat(3,1)
    bond(5,1) = rotmat(3,1)*rotmat(1,1)
    bond(6,1) = rotmat(1,1)*rotmat(2,1)
    
    ! Second column
    bond(1,2) = rotmat(1,2)*rotmat(1,2)
    bond(2,2) = rotmat(2,2)*rotmat(2,2)
    bond(3,2) = rotmat(3,2)*rotmat(3,2)
    bond(4,2) = rotmat(2,2)*rotmat(3,2)
    bond(5,2) = rotmat(3,2)*rotmat(1,2)
    bond(6,2) = rotmat(1,2)*rotmat(2,2)

    ! Third column
    bond(1,3) = rotmat(1,3)*rotmat(1,3)
    bond(2,3) = rotmat(2,3)*rotmat(2,3)
    bond(3,3) = rotmat(3,3)*rotmat(3,3)
    bond(4,3) = rotmat(2,3)*rotmat(3,3)
    bond(5,3) = rotmat(3,3)*rotmat(1,3)
    bond(6,3) = rotmat(1,3)*rotmat(2,3)

    ! Fourth column
    bond(1,4) = 2._dp*rotmat(1,2)*rotmat(1,3)
    bond(2,4) = 2._dp*rotmat(2,2)*rotmat(2,3)
    bond(3,4) = 2._dp*rotmat(3,2)*rotmat(3,3)
    bond(4,4) = rotmat(2,2)*rotmat(3,3) + rotmat(2,3)*rotmat(3,2)
    bond(5,4) = rotmat(1,2)*rotmat(3,3) + rotmat(1,3)*rotmat(3,2)
    bond(6,4) = rotmat(1,2)*rotmat(2,3) + rotmat(1,3)*rotmat(2,2)

    ! Fifth column
    bond(1,5) = 2._dp*rotmat(1,3)*rotmat(1,1)
    bond(2,5) = 2._dp*rotmat(2,3)*rotmat(2,1)
    bond(3,5) = 2._dp*rotmat(3,3)*rotmat(3,1)
    bond(4,5) = rotmat(2,1)*rotmat(3,3) + rotmat(2,3)*rotmat(3,1)
    bond(5,5) = rotmat(1,3)*rotmat(3,1) + rotmat(1,1)*rotmat(3,3)
    bond(6,5) = rotmat(1,3)*rotmat(2,1) + rotmat(1,1)*rotmat(2,3)

    ! Sixth column
    bond(1,6) = 2._dp*rotmat(1,1)*rotmat(1,2)
    bond(2,6) = 2._dp*rotmat(2,1)*rotmat(2,2)
    bond(3,6) = 2._dp*rotmat(3,1)*rotmat(3,2)
    bond(4,6) = rotmat(2,2)*rotmat(3,1) + rotmat(2,1)*rotmat(3,2)
    bond(5,6) = rotmat(1,1)*rotmat(3,2) + rotmat(1,2)*rotmat(3,1)
    bond(6,6) = rotmat(1,1)*rotmat(2,2) + rotmat(1,2)*rotmat(2,1)

  end subroutine define_bond_stress_matrix
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Define bond rotation matrices to rotate a voigt tensor
  !    (see e.g. Auld 1973, for bond matrix definition secion D pages 73-76)
  ! HERE THE STRAIN ONE => can be used to rotate the strain vector or compliance tensor
  subroutine define_bond_strain_matrix(rotmat,bond)

    real(kind=dp), dimension(3,3), intent(in)  :: rotmat
    real(kind=dp), dimension(6,6), intent(out) :: bond  
    
    ! First column
    bond(1,1) = rotmat(1,1)*rotmat(1,1)
    bond(2,1) = rotmat(2,1)*rotmat(2,1)
    bond(3,1) = rotmat(3,1)*rotmat(3,1)
    bond(4,1) = 2._dp*rotmat(2,1)*rotmat(3,1)
    bond(5,1) = 2._dp*rotmat(3,1)*rotmat(1,1)
    bond(6,1) = 2._dp*rotmat(1,1)*rotmat(2,1)
    
    ! Second column
    bond(1,2) = rotmat(1,2)*rotmat(1,2)
    bond(2,2) = rotmat(2,2)*rotmat(2,2)
    bond(3,2) = rotmat(3,2)*rotmat(3,2)
    bond(4,2) = 2._dp*rotmat(2,2)*rotmat(3,2)
    bond(5,2) = 2._dp*rotmat(3,2)*rotmat(1,2)
    bond(6,2) = 2._dp*rotmat(1,2)*rotmat(2,2)

    ! Third column
    bond(1,3) = rotmat(1,3)*rotmat(1,3)
    bond(2,3) = rotmat(2,3)*rotmat(2,3)
    bond(3,3) = rotmat(3,3)*rotmat(3,3)
    bond(4,3) = 2._dp*rotmat(2,3)*rotmat(3,3)
    bond(5,3) = 2._dp*rotmat(3,3)*rotmat(1,3)
    bond(6,3) = 2._dp*rotmat(1,3)*rotmat(2,3)

    ! Fourth column
    bond(1,4) = rotmat(1,2)*rotmat(1,3)
    bond(2,4) = rotmat(2,2)*rotmat(2,3)
    bond(3,4) = rotmat(3,2)*rotmat(3,3)
    bond(4,4) = rotmat(2,2)*rotmat(3,3) + rotmat(2,3)*rotmat(3,2)
    bond(5,4) = rotmat(1,2)*rotmat(3,3) + rotmat(1,3)*rotmat(3,2)
    bond(6,4) = rotmat(1,2)*rotmat(2,3) + rotmat(1,3)*rotmat(2,2)

    ! Fifth column
    bond(1,5) = rotmat(1,3)*rotmat(1,1)
    bond(2,5) = rotmat(2,3)*rotmat(2,1)
    bond(3,5) = rotmat(3,3)*rotmat(3,1)
    bond(4,5) = rotmat(2,1)*rotmat(3,3) + rotmat(2,3)*rotmat(3,1)
    bond(5,5) = rotmat(1,3)*rotmat(3,1) + rotmat(1,1)*rotmat(3,3)
    bond(6,5) = rotmat(1,3)*rotmat(2,1) + rotmat(1,1)*rotmat(2,3)

    ! Sixth column
    bond(1,6) = rotmat(1,1)*rotmat(1,2)
    bond(2,6) = rotmat(2,1)*rotmat(2,2)
    bond(3,6) = rotmat(3,1)*rotmat(3,2)
    bond(4,6) = rotmat(2,2)*rotmat(3,1) + rotmat(2,1)*rotmat(3,2)
    bond(5,6) = rotmat(1,1)*rotmat(3,2) + rotmat(1,2)*rotmat(3,1)
    bond(6,6) = rotmat(1,1)*rotmat(2,2) + rotmat(1,2)*rotmat(2,1)

  end subroutine define_bond_strain_matrix
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rotation of first order tensor (vector actually)
  subroutine rotate_vector(rotmat,vi,vi_r)
    
    real(kind=dp), dimension(3),   intent(in)  :: vi
    real(kind=dp), dimension(3,3), intent(in)  :: rotmat
    real(kind=dp), dimension(3),   intent(out) :: vi_r

    integer(kind=si) :: i, ip

    vi_r = 0._dp
    
    do ip = 1, 3
       do i = 1, 3
          vi_r(ip) = vi_r(ip) + rotmat(ip,i)*vi(i)
       end do
    end do
        
  end subroutine rotate_vector
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Rotation of second order tensor (not efficient but corresponds to definition)
  subroutine rotate_second_order_tensor(rotmat,cij,cij_r)
    
    real(kind=dp), dimension(3,3), intent(in)  :: cij
    real(kind=dp), dimension(3,3), intent(in)  :: rotmat
    real(kind=dp), dimension(3,3), intent(out) :: cij_r

    integer(kind=si) :: i, j, ip, jp

    cij_r = 0._dp
    
    do jp = 1, 3
       do ip = 1, 3
          do j = 1, 3
             do i = 1, 3
                cij_r(ip,jp) = cij_r(ip,jp) + rotmat(ip,i)*rotmat(jp,j)*cij(i,j)
             end do
          end do
       end do
    end do
    
  end subroutine rotate_second_order_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rotation of second order tensor (very not efficient but corresponds to definition)
  subroutine rotate_fourth_order_tensor(rotmat,cijkl,cijkl_r)

    real(kind=dp), dimension(3,3,3,3), intent(in)  :: cijkl
    real(kind=dp),     dimension(3,3), intent(in)  :: rotmat
    real(kind=dp), dimension(3,3,3,3), intent(out) :: cijkl_r

    integer(kind=si) :: i, j, k, l, ip, jp, kp, lp

    cijkl_r = 0._dp
    
    do kp = 1, 3
       do lp = 1, 3
          do jp = 1, 3
             do ip = 1, 3
                do l = 1, 3
                   do k = 1, 3
                      do j = 1, 3
                         do i = 1, 3
                            cijkl_r(ip,jp,kp,lp) = cijkl_r(ip,jp,kp,lp)        &
                                                 + rotmat(ip,i) * rotmat(jp,j) &
                                                 * rotmat(kp,k) * rotmat(lp,l) * cijkl(i,j,k,l)
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    

  end subroutine rotate_fourth_order_tensor
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Rotation of fourth order tensor in voigt matrix with bond matrix
  !    (see e.g. Auld 1973, for bond matrix definition)
  subroutine rotate_tensor_with_bond_matrix(tensor,bond,tensor_r)

    real(kind=dp), dimension(6,6), intent(in)  :: tensor
    real(kind=dp), dimension(6,6), intent(in)  :: bond
    real(kind=dp), dimension(6,6), intent(out) :: tensor_r
    
    real(kind=dp), dimension(6,6)  :: tensor_tmp, bond_t
    integer(kind=si)               :: i, j, k

    ! Get transpose of bond
    bond_t = transpose(bond)

    ! First comute CM^t
    tensor_tmp = 0._dp
    do j = 1, 6
       do k = 1, 6
          do i = 1,6
             tensor_tmp(i,j) = tensor_tmp(i,j) + tensor(i,k) * bond_t(k,j)
          end do
       end do
    end do

    ! Then Compute M*(CM^t)
    tensor_r = 0._dp
    do j = 1, 6
       do k = 1, 6
          do i = 1,6
             tensor_r(i,j) = tensor_r(i,j) + bond_t(i,j) * tensor_tmp(k,j) 
          end do
       end do
    end do
    
  end subroutine rotate_tensor_with_bond_matrix
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Get dilatational stiffness tensor according to Browaeys and Chevrot (2004)
  !   (four-rank stiffness tensor has two two-rank tensors contractions)  
  subroutine get_dilatational_stiffness_tensor(cij,dilatational)

    real(kind=dp), dimension(6,6), intent(in)  :: cij
    real(kind=dp), dimension(6,6), intent(out) :: dilatational

    ! First column
    dilatational(1,1) = cij(1,1) + cij(1,2) +cij(1,3) 
    dilatational(2,1) = cij(1,6) + cij(2,6) +cij(3,6) 
    dilatational(3,1) = cij(1,5) + cij(2,5) +cij(3,5) 

    ! Second column
    dilatational(1,2) = dilatational(2,1)
    dilatational(2,2) = cij(1,2) + cij(2,2) + cij(3,2) 
    dilatational(3,2) = cij(1,4) + cij(2,4) + cij(3,4)

    ! Third column
    dilatational(1,3) = dilatational(3,1)
    dilatational(2,3) = dilatational(3,2)
    dilatational(3,3) = cij(1,3) + cij(2,3) + cij(3,3) 
    
  end subroutine get_dilatational_stiffness_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Get voigt stiffness tensor according to Browaeys and Chevrot (2004)
  subroutine get_voigt_stiffness_tensor(cij,voigt)
    
    real(kind=dp), dimension(6,6), intent(in)  :: cij
    real(kind=dp), dimension(6,6), intent(out) :: voigt
    
    ! First column
    voigt(1,1) = cij(1,1) + cij(6,6) +cij(5,5) 
    voigt(2,1) = cij(1,6) + cij(2,6) +cij(4,5) 
    voigt(3,1) = cij(1,5) + cij(3,5) +cij(4,6) 

    ! Second column
    voigt(1,2) = voigt(2,1)
    voigt(2,2) = cij(6,6) + cij(2,2) + cij(4,4) 
    voigt(3,2) = cij(2,4) + cij(3,4) + cij(5,6)

    ! Third column
    voigt(1,3) = voigt(3,1)
    voigt(2,3) = voigt(3,2)
    voigt(3,3) = cij(5,5) + cij(4,4) + cij(3,3) 
    
  end subroutine get_voigt_stiffness_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Pass from elastic tensor to the elastic vector defined by Browaeys and Chevrot (2004)
  subroutine elastic_tensor_to_elastic_vector


  end subroutine elastic_tensor_to_elastic_vector
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Pass from elastic vector to elastic tensor (see Browaeys and Chevrot (2004))
  subroutine elastic_vector_to_elastic_tensor


  end subroutine elastic_vector_to_elastic_tensor
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Pass triclinic elastic vector to isotropic one according to fedorov (1968)
  subroutine get_isotropic_part_fedorov(triclinic,isotropic)

    real(kind=dp), dimension(21), intent(in)  :: triclinic
    real(kind=dp), dimension(21), intent(out) :: isotropic

    real(kind=dp) :: kappa, mu, lambda, c11, c22, c33, c23, c13, c12, c44, c55, c66, lp2m

    c11 = triclinic(1);   c22=triclinic(2);   c33=triclinic(3);
    c23 = triclinic(4);   c13=triclinic(5);   c12=triclinic(6);
    c44 = triclinic(7);   c55=triclinic(8);   c66=triclinic(9);

    kappa  = (c11 + c22 + c33 + 2._dp*(c12 + c13 + c23)) / 9._dp
    mu     = (2._dp*(c11 + c22 + c33 - c12 - c23 - c13) + 6._dp * (c44 + c55 + c66)) / 30._dp
    lambda = kappa - (2._dp*mu /3._dp)

    lp2m = lambda + 2._dp * mu
    
    isotropic(:)   = 0._dp
    isotropic(1:3) = lp2m
    isotropic(4:6) = lambda
    isotropic(7:9) = mu
    
  end subroutine get_isotropic_part_fedorov
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define isotropic tensor with analytical formula (not efficient but still usefull)
  subroutine define_isotropic_tensor_1(lambda,mu,tensor)

    real(kind=dp), intent(in)                 :: lambda, mu
    real(kind=dp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl

    tensor(:) = 0._dp
    
    !* Loop over cij components
    do ipar = 1, 9  ! not necessary to go until 21 here

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)
       
       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)
       
       !*** Fill elastic vector
       tensor(ipar) = lambda*dij*dkl + mu*(dik*djl + dil*djk) 
       
    end do

  end subroutine define_isotropic_tensor_1
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define hexagonal tensor with analytical fromula from cij values
  subroutine define_hexagonal_tensor_1(c11,c33,c44,c66,c13,s,tensor)

    real(kind=dp), intent(in)               :: c11, c33, c44, c66, c13
    real(kind=dp), dimension(3), intent(in) :: s
    
    real(kind=dp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl
    real(kind=dp)    :: si, sj, sk, sl
    real(kind=dp)    :: sisj, sksl, sjsk, sisl, sisk, sjsl
    real(kind=dp)    :: sisjsk, sisjsl, sisksl, sjsksl, sisjsksl
    
    !*** Init
    tensor(:) = 0._dp

    !*** Precompute direction cosines
    si = s(i)
    sj = s(j)
    sk = s(k)
    sl = s(l)
    sisj = si * sj
    sksl = sk * sl
    sjsk = sj * sk
    sisl = si * sl
    sisk = si * sk
    sjsl = sj * sl
    sisjsk = sisj * sk
    sisjsl = sisj * sl
    sisksl = sisk * sl
    sjsksl = sjsk * sl
    sisjsksl = sisjsk * sl

    !* Loop over cij components
    do ipar = 1, 21  ! not necessary to go until 21 here

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)
       
       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)
       
       !*** Fill elastic vector
       tensor(ipar) = (c11 - 2._dp*c66) * dij*dkl + c66 * (dik*djl + dil*djk) &
            + (c13 - c11 + 2._dp*c66) * (dij*sksl + dkl*sisj)                  &
            + (c44 - c66) * (dik*sjsl + dil*sjsk + djk*sisl + djl*sisk)        &
            + (c11 + c33 - 2._dp*c13 - 4._dp*c44)*sisjsksl
       
    end do
    
  end subroutine define_hexagonal_tensor_1
  !--------------------------------------------------------------------------------


  !================================================================================
  ! Define hexagonal tensor with analytical fromula from thomsen parameters
  subroutine define_hexagonal_tensor_2(c33,c44,eps,del,gam,s,tensor)

    real(kind=dp), intent(in)               :: c33, c44, eps, del, gam
    real(kind=dp), dimension(3), intent(in) :: s
    
    real(kind=dp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl
    real(kind=dp)    :: si, sj, sk, sl
    real(kind=dp)    :: sisj, sksl, sjsk, sisl, sisk, sjsl
    real(kind=dp)    :: sisjsk, sisjsl, sisksl, sjsksl, sisjsksl
    
    !*** Init
    tensor(:) = 0._dp

    !*** Precompute direction cosines
    si = s(i)
    sj = s(j)
    sk = s(k)
    sl = s(l)
    sisj = si * sj
    sksl = sk * sl
    sjsk = sj * sk
    sisl = si * sl
    sisk = si * sk
    sjsl = sj * sl
    sisjsk = sisj * sk
    sisjsl = sisj * sl
    sisksl = sisk * sl
    sjsksl = sjsk * sl
    sisjsksl = sisjsk * sl

    !* Loop over cij components
    do ipar = 1, 21  ! not necessary to go until 21 here

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)
       
       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)
       
       !*** Fill elastic vector
       tensor(ipar) = c33 * dij*dkl + c44 * (dik*djl + dil*djk - 2._dp*dij*dkl) &
            + 2._dp * eps * c33 * (dij*dkl - dij*sksl - dkl*sisj + sisjsksl)    &
            +         del * c33 * (dij*sksl + dkl*sisj - 2._dp*sisjsksl)        &
            + 2._dp * gam * c44 * (-2._dp*dij*dkl + dik*djl + dil*djk           &
                                   +2._dp*dij*sksl + 2._dp*dkl*sisj             &
                                   - dik*sjsl - dil*sjsk - djk*sisl - djl*sisk)
                              
    end do
    
  end subroutine define_hexagonal_tensor_2
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define orthorhombic tensor with analytical fromula from cij values
  !    a, b, and c are normal vector to the three mirror symmetry planes
  subroutine define_orthorhombic_tensor1(c11,c22,c33,c23,c13,c12,c44,c55,c66,a,b,c,tensor)

    real(kind=dp), intent(in)               :: c11, c22, c33, c44, c55, c66, c13, c23, c12
    real(kind=dp), dimension(3), intent(in) :: a, b, c
    
    real(kind=dp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl
    real(kind=dp)    :: aiajakal, bibjbkbl, cicjckcl
    real(kind=dp)    :: aiaj, akal, bibj, bkbl
    real(kind=dp)    :: aiak, aial, ajal, ajak
    real(kind=dp)    :: bibk, bibl, bjbl, bjbk
    real(kind=dp)    :: ai, aj, ak, al, bi, bj, bk, bl, ci, cj, ck, cl
    
    !*** Init
    tensor(:) = 0._dp

    !*** Precompute direction cosines
    ai = a(i);   aj = a(j);   ak = a(k);   al = a(l);
    bi = b(i);   bj = b(j);   bk = b(k);   bl = b(l);
    ci = c(i);   cj = c(j);   ck = c(k);   cl = c(l);
    
    aiaj = ai*aj;   akal = ak*al;   aiak = ai*ak;   aial = ai*al;
    bibj = bi*bj;   bkbl = bk*bl;   ajal = aj*al;   ajak = aj*ak;
    bibk = bi*bk;   bibl = bi*bl;   bjbl = bj*bl;   bjbk = bj*bk;
    
    aiajakal = ai*aj*ak*al
    bibjbkbl = bi*bj*bk*bl
    cicjckcl = ci*cj*ck*cl

    !* Loop over cij components
    do ipar = 1, 21  ! not necessary to go until 21 here

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)
       
       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)
       
       !*** Fill elastic vector
       tensor(ipar) = (c13 + c23 - c12) * dij*dkl + (c55 - c66 + c44) * (dik*djl + dil*djk) &
            + (c12 - c23) * (aiaj*dkl + dij*akal) + (c12-c13)*(bibj*dkl + dij*bkbl)         &
            + (c66 - c44) * (aiak*djl + aial*djk + ajal*dik + ajak*dil)                     &
            + (c66 - c55) * (bibk*djl + bibl*djk + bjbl*dik + bjbk*dil)                     &
            + (c11 - c13 + c23 - c12 + 2._dp*( c44 - c55 - c66)) * aiajakal                 &
            + (c22 + c13 - c23 - c12 + 2._dp*(-c44 + c55 - c66)) * bibjbkbl                 &
            + (c33 - c13 - c23 + c12 + 2._dp*(-c44 - c55 + c66)) * cicjckcl
    end do
    
  end subroutine define_orthorhombic_tensor1
  !--------------------------------------------------------------------------------

  
  !================================================================================
  ! Routine for initialisation of tensor, voigt matrix and elastic vector indexing
  subroutine define_indexing_vec_to_tens
    
    !*** Define indexing to pass from tensor to vector (see browaeys and chevrot 2004)
    !    IN TENSOR INDEXING
    if (.not.allocated(ind_vec2tens)) then

       allocate(ind_vec2tens(4,21))
       
       ! Group 1
       ind_vec2tens(1:4,1)  = (/ 1, 1, 1, 1 /)
       ind_vec2tens(1:4,2)  = (/ 2, 2, 2, 2 /)
       ind_vec2tens(1:4,3)  = (/ 3, 3, 3, 3 /)
       
       ! Group 2
       ind_vec2tens(1:4,4)  = (/ 2, 2, 3, 3 /)
       ind_vec2tens(1:4,5)  = (/ 1, 1, 3, 3 /)
       ind_vec2tens(1:4,6)  = (/ 1, 1, 2, 2 /)
       
       ! Group 3
       ind_vec2tens(1:4,7)  = (/ 2, 3, 2, 3 /)
       ind_vec2tens(1:4,8)  = (/ 1, 3, 1, 3 /)
       ind_vec2tens(1:4,9)  = (/ 1, 2, 1, 2 /)
       
       ! Group 4
       ind_vec2tens(1:4,10) = (/ 1, 1, 2, 3 /)
       ind_vec2tens(1:4,11) = (/ 2, 2, 1, 3 /)
       ind_vec2tens(1:4,12) = (/ 3, 3, 1, 2 /)
       
       ! Group 5
       ind_vec2tens(1:4,13) = (/ 3, 3, 2, 3 /)
       ind_vec2tens(1:4,14) = (/ 1, 1, 1, 3 /)
       ind_vec2tens(1:4,15) = (/ 2, 2, 1, 2 /)
       ind_vec2tens(1:4,16) = (/ 2, 2, 2, 3 /)
       ind_vec2tens(1:4,17) = (/ 3, 3, 1, 3 /)
       ind_vec2tens(1:4,18) = (/ 1, 1, 1, 2 /)
       
       ! Group 6
       ind_vec2tens(1:4,19) = (/ 1, 3, 1, 2 /)
       ind_vec2tens(1:4,20) = (/ 2, 3, 1, 2 /)
       ind_vec2tens(1:4,21) = (/ 2, 3, 1, 3 /)

    end if
       
    !*** Define indexing to pass from tensor to vector (see browaeys and chevrot 2004)
    !    IN VOIGT INDEXING
    if (.not.allocated(ind_vec2tens_voigt)) then

       allocate(ind_vec2tens_voigt(2,21))
       
       ! Group 1
       ind_vec2tens_voigt(1:2,1)  = (/ 1, 1 /)
       ind_vec2tens_voigt(1:2,2)  = (/ 2, 2 /)
       ind_vec2tens_voigt(1:2,3)  = (/ 3, 3 /)
       
       ! Group 2
       ind_vec2tens_voigt(1:2,4)  = (/ 2, 3 /)
       ind_vec2tens_voigt(1:2,5)  = (/ 1, 3 /)
       ind_vec2tens_voigt(1:2,6)  = (/ 1, 2 /)
       
       ! Group 3
       ind_vec2tens_voigt(1:2,7)  = (/ 4, 4 /)
       ind_vec2tens_voigt(1:2,8)  = (/ 5, 5 /)
       ind_vec2tens_voigt(1:2,9)  = (/ 6, 6 /)
       
       ! Group 4
       ind_vec2tens_voigt(1:2,10) = (/ 1, 4 /)
       ind_vec2tens_voigt(1:2,11) = (/ 2, 5 /)
       ind_vec2tens_voigt(1:2,12) = (/ 3, 6 /)
       
       ! Group 5
       ind_vec2tens_voigt(1:2,13) = (/ 3, 4 /)
       ind_vec2tens_voigt(1:2,14) = (/ 1, 5 /)
       ind_vec2tens_voigt(1:2,15) = (/ 2, 6 /)
       ind_vec2tens_voigt(1:2,16) = (/ 2, 4 /)
       ind_vec2tens_voigt(1:2,17) = (/ 3, 5 /)
       ind_vec2tens_voigt(1:2,18) = (/ 1, 6 /)
       
       ! Group 6
       ind_vec2tens_voigt(1:2,19) = (/ 5, 6 /)
       ind_vec2tens_voigt(1:2,20) = (/ 4, 6 /)
       ind_vec2tens_voigt(1:2,21) = (/ 4, 5 /)

    end if
       
  end subroutine define_indexing_vec_to_tens
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Get direction cosines (th = dip angle from vertical, ph = azimuth from north (y))
  subroutine get_direction_cosines(th,ph,s)

    real(kind=dp),               intent(in)  :: th, ph
    real(kind=dp), dimension(3), intent(out) :: s
    
    real(kind=dp) :: thrad, phrad

    thrad = th*deg2rad
    phrad = ph*deg2rad

    s(1) = sin(phrad) * sin(thrad)
    s(2) = cos(phrad) * sin(thrad)
    s(3) = cos(thrad)
    
  end subroutine get_direction_cosines
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Get angles of symetry axis from direction cosines
  !     (th = dip angle from vertical, ph = azimuth from north (y))
  subroutine get_symmetry_angles(s,th,ph)

    real(kind=dp), dimension(3), intent(in)  :: s
    real(kind=dp),               intent(out) :: th, ph
    
    real(kind=dp) :: thrad, phrad, r

    r = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
    ph = atan2(s(1),s(2))
    th = acos(s(3)/r)
    
    th = thrad*rad2deg
    ph = phrad*rad2deg
    
  end subroutine get_symmetry_angles
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Small function for kronecker delta function
  function delta(i,j) result(d)

    integer(kind=si), intent(in) :: i, j
    integer(kind=si)             :: d
    
    if (i==j) then
       d = 1
    else
       d = 0
    end if

  end function delta
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Compute square norm of vector
  function square_norm_vector(vector) result(norm)

    real(kind=dp), dimension(3), intent(in) :: vector
    real(kind=dp)                           :: norm

    integer(kind=si) :: i

    norm = 0._dp

    do i = 1, 3
       norm = norm + vector(i)*vector(i)
    end do
            
  end function square_norm_vector
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Compute square norm of a second order tensor
  function square_norm_tensor_2(tensor) result(norm)

    real(kind=dp), dimension(3,3), intent(in) :: tensor
    real(kind=dp)                             :: norm

    integer(kind=si) :: i, j

    norm = 0._dp

    do j = 1, 3
       do i = 1, 3
          norm = norm + tensor(i,j)*tensor(i,j)
       end do
    end do
        
  end function square_norm_tensor_2
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Compute square norm of a second order tensor
  function square_norm_tensor_4(tensor) result(norm)

    real(kind=dp), dimension(3,3,3,3), intent(in) :: tensor
    real(kind=dp)                                 :: norm

    integer(kind=si) :: i, j, k, l

    norm = 0._dp

    do l = 1,3
       do k = 1, 3
          do j = 1, 3
             do i = 1, 3
                norm = norm + tensor(i,j,k,l)*tensor(i,j,k,l)
             end do
          end do
       end do
    end do
    
  end function square_norm_tensor_4
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Pass a voigt matrix to elastic vector
  function transform_voigt_matrix_to_vector(cij) result(vi)

    real(kind=dp), dimension(6,6), intent(in) :: cij
    real(kind=dp), dimension(21)              :: vi
    
    real(kind=dp) :: sqrt_two, two_sqrt_two
    
    ! Constants
    sqrt_two     = sqrt(2._dp)
    two_sqrt_two = 2._dp * sqrt_two

    ! Group 1
    vi( 1) = cij(1,1)
    vi( 2) = cij(2,2)
    vi( 3) = cij(3,3)
    
    ! Group 2
    vi( 4) = cij(2,3) * sqrt_two
    vi( 5) = cij(1,3) * sqrt_two
    vi( 6) = cij(1,2) * sqrt_two
    
    ! Group 3
    vi( 7) = cij(4,4)
    vi( 8) = cij(5,5)
    vi( 9) = cij(6,6)
    
    ! Group 4
    vi(10) = cij(1,4)
    vi(11) = cij(2,5)
    vi(12) = cij(3,6)
    
    ! Group 5
    vi(13) = cij(3,4)
    vi(14) = cij(1,5)
    vi(15) = cij(2,6)
    vi(16) = cij(2,4)
    vi(17) = cij(3,5)
    vi(18) = cij(1,6)
    
    ! Group 6
    vi(19) = cij(5,6) * two_sqrt_two
    vi(20) = cij(4,6) * two_sqrt_two
    vi(21) = cij(4,5) * two_sqrt_two
    
  end function transform_voigt_matrix_to_vector
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Pass a voigt tensor to elastic vector
  function transform_kelvin_tensor_to_vector(cij) result(vi)

    real(kind=dp), dimension(6,6), intent(in) :: cij
    real(kind=dp), dimension(21)              :: vi
    
    real(kind=dp) :: sqrt_two
    
    ! Constants
    sqrt_two     = sqrt(2._dp)

    ! Group 1
    vi( 1) = cij(1,1)
    vi( 2) = cij(2,2)
    vi( 3) = cij(3,3)
    
    ! Group 2
    vi( 4) = cij(2,3) * sqrt_two
    vi( 5) = cij(1,3) * sqrt_two
    vi( 6) = cij(1,2) * sqrt_two
    
    ! Group 3
    vi( 7) = cij(4,4)
    vi( 8) = cij(5,5)
    vi( 9) = cij(6,6)
    
    ! Group 4
    vi(10) = cij(1,4) * sqrt_two
    vi(11) = cij(2,5) * sqrt_two
    vi(12) = cij(3,6) * sqrt_two
    
    ! Group 5
    vi(13) = cij(3,4) * sqrt_two
    vi(14) = cij(1,5) * sqrt_two
    vi(15) = cij(2,6) * sqrt_two
    vi(16) = cij(2,4) * sqrt_two
    vi(17) = cij(3,5) * sqrt_two
    vi(18) = cij(1,6) * sqrt_two
    
    ! Group 6
    vi(19) = cij(5,6) * sqrt_two
    vi(20) = cij(4,6) * sqrt_two
    vi(21) = cij(4,5) * sqrt_two
    
  end function transform_kelvin_tensor_to_vector
  !--------------------------------------------------------------------------------


  !================================================================================
  ! Pass an elastic vector to a voigt tensor
  !   (see Browaeys and Chevrot (2004)
  function transform_vector_to_voigt_tensor(vi) result(cij)

    real(kind=dp), dimension(21),  intent(in)  :: vi
    real(kind=dp), dimension(6,6)              :: cij
    
    real(kind=dp) :: inv_sqrt_two, inv_two_sqrt_two

    ! Constants
    inv_sqrt_two     = 1._dp / sqrt(2._dp)
    inv_two_sqrt_two = 0.5_dp * inv_sqrt_two
     
    ! First column
    cij(1,1) = vi( 1)
    cij(2,1) = vi( 6) * inv_sqrt_two
    cij(3,1) = vi( 5) * inv_sqrt_two
    cij(4,1) = vi(10) * 0.5_dp
    cij(5,1) = vi(14) * 0.5_dp
    cij(6,1) = vi(18) * 0.5_dp

    ! Second column
    cij(1,2) = cij(2,1)
    cij(2,2) = vi( 2)
    cij(3,2) = vi( 4) * inv_sqrt_two
    cij(4,2) = vi(16) * 0.5_dp
    cij(5,2) = vi(11) * 0.5_dp
    cij(6,2) = vi(15) * 0.5_dp

    ! Third column
    cij(1,3) = cij(3,1)
    cij(2,3) = cij(3,2)
    cij(3,3) = vi( 3)
    cij(4,3) = vi(13) * 0.5_dp
    cij(5,3) = vi(17) * 0.5_dp
    cij(6,3) = vi(12) * 0.5_dp

    ! Fourth column
    cij(1,4) = cij(4,1)
    cij(2,4) = cij(4,2)
    cij(3,4) = cij(4,3)
    cij(4,4) = vi(7)  * 0.5_dp
    cij(5,4) = vi(21) * inv_two_sqrt_two
    cij(6,4) = vi(20) * inv_two_sqrt_two

    ! Fifth column
    cij(1,5) = cij(5,1)
    cij(2,5) = cij(5,2)
    cij(3,5) = cij(5,3)
    cij(4,5) = cij(5,4)
    cij(5,5) = vi( 8) * 0.5_dp
    cij(6,5) = vi(19) * inv_two_sqrt_two

    ! Sixth column
    cij(1,6) = cij(6,1)
    cij(2,6) = cij(6,2)
    cij(3,6) = cij(6,3)
    cij(4,6) = cij(6,4)
    cij(5,6) = cij(6,5)
    cij(6,6) = vi( 9) * 0.5_dp

  end function transform_vector_to_voigt_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Get Voigt m index from ij
  function voigt_index(i,j) result(m)

    integer(kind=si), intent(in) :: i, j
    integer(kind=si)             :: m, dij

    dij = delta(i,j) 

    m = dij*i + (1-dij)*(9-i-j) 

  end function voigt_index
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Transform fourth order tensor to second order voigt matrix
  function transform_tensor_fourth_to_voigt_matrix(cijkl) result(cij)

    real(kind=dp), dimension(3,3,3,3), intent(in) :: cijkl
    real(kind=dp), dimension(3,3)                 :: cij
    
    ! First column
    cij(1,1) =  cijkl(1,1,1,1)
    cij(2,1) = (cijkl(2,2,1,1) + cijkl(1,1,2,2)) * 0.5_dp
    cij(3,1) = (cijkl(3,3,1,1) + cijkl(1,1,3,3)) * 0.5_dp
    cij(4,1) = (cijkl(2,3,1,1) + cijkl(3,2,1,1) + &
                cijkl(1,1,2,3) + cijkl(1,1,3,2)) * 0.25_dp
    cij(5,1) = (cijkl(1,3,1,1) + cijkl(3,1,1,1) + &
                cijkl(1,1,1,3) + cijkl(1,1,3,1)) * 0.25_dp
    cij(6,1) = (cijkl(1,2,1,1) + cijkl(2,1,1,1) + &
                cijkl(1,1,2,1) + cijkl(1,1,1,2)) * 0.25_dp

    ! Second column
    cij(1,2) =  cij(2,1)
    cij(2,2) =  cijkl(2,2,2,2)
    cij(3,2) = (cijkl(3,3,2,2) + cijkl(2,2,3,3)) * 0.5_dp
    cij(4,2) = (cijkl(2,3,2,2) + cijkl(3,2,2,2) + &
                cijkl(2,2,2,3) + cijkl(2,2,3,2)) * 0.25_dp
    cij(5,2) = (cijkl(1,3,2,2) + cijkl(3,1,2,2) + &
                cijkl(2,2,1,3) + cijkl(2,2,3,1)) * 0.25_dp
    cij(6,2) = (cijkl(1,2,2,2) + cijkl(2,1,2,2) + &
                cijkl(2,2,1,2) + cijkl(2,2,2,1)) * 0.25_dp

    ! Third column
    cij(1,3) =  cij(3,1)
    cij(2,3) =  cij(3,2)
    cij(3,3) =  cijkl(3,3,3,3)
    cij(4,3) = (cijkl(2,3,3,3) + cijkl(3,2,3,3) + &
                cijkl(3,3,2,3) + cijkl(3,3,3,2)) * 0.25_dp
    cij(5,3) = (cijkl(1,3,3,3) + cijkl(3,1,3,3) + &
                cijkl(3,3,1,3) + cijkl(3,3,3,1)) * 0.25_dp
    cij(6,3) = (cijkl(1,2,3,3) + cijkl(2,1,3,3) + &
                cijkl(3,3,1,2) + cijkl(3,3,2,1)) * 0.25_dp

    ! Fourth column
    cij(1,4) =  cij(4,1)
    cij(2,4) =  cij(4,2)
    cij(3,4) =  cij(3,4)
    cij(4,4) = (cijkl(2,3,2,3) + cijkl(2,3,3,2) + &
                cijkl(3,2,2,3) + cijkl(3,2,3,2)) * 0.25_dp
    cij(5,4) = (cijkl(1,3,2,3) + cijkl(3,1,2,3) + &
                cijkl(1,3,3,2) + cijkl(3,1,3,2) + &
                cijkl(2,3,1,3) + cijkl(2,3,3,1) + &
                cijkl(3,2,1,3) + cijkl(3,2,3,1)) * 0.125_dp
    cij(6,4) = (cijkl(1,2,2,3) + cijkl(2,1,2,3) + &
                cijkl(1,2,3,2) + cijkl(2,1,3,2) + &
                cijkl(2,3,1,2) + cijkl(2,3,2,1) + &
                cijkl(3,2,1,2) + cijkl(3,2,2,1)) * 0.125_dp
    
    ! Fifth column
    cij(1,5) =  cij(5,1)
    cij(2,5) =  cij(5,2)
    cij(3,5) =  cij(5,3)
    cij(4,5) =  cij(5,4)
    cij(5,5) = (cijkl(1,3,1,3) + cijkl(1,3,3,1) + &
                cijkl(3,1,1,3) + cijkl(3,1,3,1)) * 0.25_dp
    cij(6,5) = (cijkl(1,3,2,3) + cijkl(3,1,2,3) + &
                cijkl(1,3,3,2) + cijkl(3,1,3,2) + &
                cijkl(2,3,1,3) + cijkl(2,3,3,1) + &
                cijkl(3,2,1,3) + cijkl(3,2,3,1)) * 0.125_dp

    ! Sixth column
    cij(1,6) =  cij(6,1)
    cij(2,6) =  cij(6,2)
    cij(3,6) =  cij(6,3)
    cij(4,6) =  cij(6,4)
    cij(5,6) =  cij(6,5)
    cij(6,6) = (cijkl(2,1,2,1) + cijkl(2,1,1,2) +  &
                cijkl(1,2,2,1) + cijkl(1,2,2,1)) * 0.25_dp
    
  end function transform_tensor_fourth_to_voigt_matrix
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Transform second order voigt tensor to fourth order tensor
  function transform_voigt_matrix_to_fourth_order_tensor(cmn) result(cijkl)

    real(kind=dp), dimension(3,3),    intent(in) :: cmn
    real(kind=dp), dimension(3,3,3,3)            :: cijkl

    integer(kind=si) :: i, j, k, l, m, n

    ! Make tensor
    cijkl = 0._dp
    do l = 1, 3
       do k = 1, 3

          ! Get second index
          n = voigt_index(k,l)
          
          do j = 1, 3
             do i = 1, 3

                ! Get first index
                m = voigt_index(i,j)                

                ! Fill tensor
                cijkl(i,j,k,l) = cmn(m,n)
                
             end do
          end do
       end do
    end do

  end function transform_voigt_matrix_to_fourth_order_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Transfrom voigt matrix 6x6 to kelvin tensor 6x6
  function voigt_matrix_to_kelvin_tensor(voigt) result(kelvin)

    real(kind=dp), dimension(6,6), intent(in) :: voigt
    real(kind=dp), dimension(6,6)             :: kelvin
    
    real(kind=dp) :: sqrt_two
    
    ! Constants
    sqrt_two     = sqrt(2._dp)

    ! Transform
    kelvin(1:3,1:3) = voigt(1:3,1:3)
    kelvin(1:3,4:6) = voigt(1:3,4:6) * sqrt_two
    kelvin(4:6,1:3) = voigt(4:6,1:3) * sqrt_two
    kelvin(4:6,4:6) = voigt(4:6,4:6) * 2._dp

  end function voigt_matrix_to_kelvin_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Transfrom voigt matrix 6x6 to kelvin tensor 6x6
  function kelvin_tensor_to_voigt_matrix(kelvin) result(voigt)

    real(kind=dp), dimension(6,6), intent(in) :: kelvin
    real(kind=dp), dimension(6,6)             :: voigt
    
    real(kind=dp) :: inv_sqrt_two
    
    ! Constants
    inv_sqrt_two = 1._dp / sqrt(2._dp)

    ! Transform
    voigt(1:3,1:3) = kelvin(1:3,1:3)
    voigt(1:3,4:6) = kelvin(1:3,4:6) * inv_sqrt_two
    voigt(4:6,1:3) = kelvin(4:6,1:3) * inv_sqrt_two
    voigt(4:6,4:6) = kelvin(4:6,4:6) * 0.5_dp

  end function kelvin_tensor_to_voigt_matrix
  !--------------------------------------------------------------------------------
  
end module elastic_tensor_tools_mod
