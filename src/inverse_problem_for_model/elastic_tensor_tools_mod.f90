module elastic_tensor_tools_mod
  
  use interpolation_mod, only: si, sp, di, dp, cp, hp, deg2rad, rad2deg

  implicit none

  integer(kind=si), dimension(:,:), allocatable :: ind_vec2tens, ind_vec2tens_voigt
  
contains

  !================================================================================
  ! Define rotation matrix with euler angles 
  subroutine define_rotation_matrix


  end subroutine define_rotation_matrix
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define bond rotation matrix to rotate a voigt tensor
  !    (see e.g. Auld 1973, for bond matrix definition)
  subroutine define_bond_matrix


  end subroutine define_bond_matrix
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rotation of second order tensor (not efficient but corresponds to definition)
  subroutine rotate_second_order_tensor


  end subroutine rotate_second_order_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rotation of second order tensor (very not efficient but corresponds to definition)
  subroutine rotate_fourth_order_tensor


  end subroutine rotate_fourth_order_tensor
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Rotation of fourth order tensor in voigt matrix with bond matrix
  !    (see e.g. Auld 1973, for bond matrix definition)
  subroutine rotate_tensor_with_bond_matrix


  end subroutine rotate_tensor_with_bond_matrix
  !--------------------------------------------------------------------------------
  
  !================================================================================
  ! Get dilatational stiffness tensor according to Browaeys and Chevrot (2004)
  !   (four-rank stiffness tensor has two two-rank tensors contractions)  
  subroutine get_dilatational_stiffness_tensor(cij,dilatational)

    real(kind=cp), dimension(6,6), intent(in)  :: cij
    real(kind=cp), dimension(6,6), intent(out) :: dilatational

    ! First column
    dilatational(1,1) = cij(1,1) + cij(1,2) +cij(1,3) 
    dilatational(2,1) = cij(1,6) + cij(2,6) +cij(3,6) 
    dilatational(3,1) = cij(1,5) + cij(2,5) +cij(3,5) 

    ! Second column
    dilatational(1,2) = dilatational(2,1)
    dilatational(2,2) = cij(1,2) + cij(2,2) + cij(3,2) 
    dilatational(3,2) = cij(1,4) + cij(2,4) + cij(3,4)

    ! Thirs column
    dilatational(1,3) = dilatational(3,1)
    dilatational(2,3) = dilatational(3,2)
    dilatational(3,3) = cij(1,3) + cij(2,3) + cij(3,3) 
    
  end subroutine get_dilatational_stiffness_tensor
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Get voigt stiffness tensor according to Browaeys and Chevrot (2004)
  subroutine get_voigt_stiffness_tensor(cij,voigt)
    
    real(kind=cp), dimension(6,6), intent(in)  :: cij
    real(kind=cp), dimension(6,6), intent(out) :: voigt
    
    ! First column
    voigt(1,1) = cij(1,1) + cij(6,6) +cij(5,5) 
    voigt(2,1) = cij(1,6) + cij(2,6) +cij(4,5) 
    voigt(3,1) = cij(1,5) + cij(3,5) +cij(4,6) 

    ! Second column
    voigt(1,2) = voigt(2,1)
    voigt(2,2) = cij(6,6) + cij(2,2) + cij(4,4) 
    voigt(3,2) = cij(2,4) + cij(3,4) + cij(5,6)

    ! Thirs column
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

    real(kind=cp), dimension(21), intent(in)  :: triclinic
    real(kind=cp), dimension(21), intent(out) :: isotropic

    real(kind=cp) :: kappa, mu, lambda, c11, c22, c33, c23, c13, c12, c44, c55, c66, lp2m

    c11 = triclinic(1);   c22=triclinic(2);   c33=triclinic(3);
    c23 = triclinic(4);   c13=triclinic(5);   c12=triclinic(6);
    c44 = triclinic(7);   c55=triclinic(8);   c66=triclinic(9);

    kappa  = (c11 + c22 + c33 + 2._cp*(c12 + c13 + c23)) / 9._cp
    mu     = (2._cp*(c11 + c22 + c33 - c12 - c23 - c13) + 6._cp * (c44 + c55 + c66)) / 30._cp
    lambda = kappa - (2._cp*mu /3._cp)

    lp2m = lambda + 2._cp * mu
    
    isotropic(:)   = 0._cp
    isotropic(1:3) = lp2m
    isotropic(4:6) = lambda
    isotropic(7:9) = mu
    
  end subroutine get_isotropic_part_fedorov
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define isotropic tensor with analytical formula (not efficient but still usefull)
  subroutine define_isotropic_tensor_1(lambda,mu,tensor)

    real(kind=cp), intent(in)                 :: lambda, mu
    real(kind=cp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl

    tensor(:) = 0._cp
    
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

    real(kind=cp), intent(in)               :: c11, c33, c44, c66, c13
    real(kind=cp), dimension(3), intent(in) :: s
    
    real(kind=cp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl
    real(kind=cp)    :: si, sj, sk, sl
    real(kind=cp)    :: sisj, sksl, sjsk, sisl, sisk, sjsl
    real(kind=cp)    :: sisjsk, sisjsl, sisksl, sjsksl, sisjsksl
    
    !*** Init
    tensor(:) = 0._cp

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
       tensor(ipar) = (c11 - 2._cp*c66) * dij*dkl + c66 * (dik*djl + dil*djk) &
            + (c13 - c11 + 2._cp*c66) * (dij*sksl + dkl*sisj)                  &
            + (c44 - c66) * (dik*sjsl + dil*sjsk + djk*sisl + djl*sisk)        &
            + (c11 + c33 - 2._cp*c13 - 4._cp*c44)*sisjsksl
       
    end do
    
  end subroutine define_hexagonal_tensor_1
  !--------------------------------------------------------------------------------


  !================================================================================
  ! Define hexagonal tensor with analytical fromula from thomsen parameters
  subroutine define_hexagonal_tensor_2(c33,c44,eps,del,gam,s,tensor)

    real(kind=cp), intent(in)               :: c33, c44, eps, del, gam
    real(kind=cp), dimension(3), intent(in) :: s
    
    real(kind=cp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl
    real(kind=cp)    :: si, sj, sk, sl
    real(kind=cp)    :: sisj, sksl, sjsk, sisl, sisk, sjsl
    real(kind=cp)    :: sisjsk, sisjsl, sisksl, sjsksl, sisjsksl
    
    !*** Init
    tensor(:) = 0._cp

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
       tensor(ipar) = c33 * dij*dkl + c44 * (dik*djl + dil*djk - 2._cp*dij*dkl) &
            + 2._cp * eps * c33 * (dij*dkl - dij*sksl - dkl*sisj + sisjsksl)    &
            +         del * c33 * (dij*sksl + dkl*sisj - 2._cp*sisjsksl)        &
            + 2._cp * gam * c44 * (-2._cp*dij*dkl + dik*djl + dil*djk           &
                                   +2._cp*dij*sksl + 2._cp*dkl*sisj             &
                                   - dik*sjsl - dil*sjsk - djk*sisl - djl*sisk)
                              
    end do
    
  end subroutine define_hexagonal_tensor_2
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Define orthorhombic tensor with analytical fromula from cij values
  !    a, b, and c are normal vector to the three mirror symmetry planes
  subroutine define_orthorhombic_tensor1(c11,c22,c33,c23,c13,c12,c44,c55,c66,a,b,c,tensor)

    real(kind=cp), intent(in)               :: c11, c22, c33, c44, c55, c66, c13, c23, c12
    real(kind=cp), dimension(3), intent(in) :: a, b, c
    
    real(kind=cp), dimension(21), intent(out) :: tensor

    integer(kind=si) :: ipar, i, j, k, l
    integer(kind=si) :: dij, dik, dil, djk, djl, dkl
    real(kind=cp)    :: aiajakal, bibjbkbl, cicjckcl
    real(kind=cp)    :: aiaj, akal, bibj, bkbl
    real(kind=cp)    :: aiak, aial, ajal, ajak
    real(kind=cp)    :: bibk, bibl, bjbl, bjbk
    real(kind=cp)    :: ai, aj, ak, al, bi, bj, bk, bl, ci, cj, ck, cl
    
    !*** Init
    tensor(:) = 0._cp

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

    real(kind=cp),               intent(in)  :: th, ph
    real(kind=cp), dimension(3), intent(out) :: s
    
    real(kind=cp) :: thrad, phrad

    thrad = th*deg2rad
    phrad = ph*deg2rad

    s(1) = sin(phrad) * sin(thrad)
    s(2) = cos(phrad) * sin(thrad)
    s(3) = cos(thrad)
    
  end subroutine get_direction_cosines
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Get angles of symétry axis from direction cosines
  !     (th = dip angle from vertical, ph = azimuth from north (y))
  subroutine get_symmetry_angles(s,th,ph)

    real(kind=cp), dimension(3), intent(in)  :: s
    real(kind=cp),               intent(out) :: th, ph
    
    real(kind=cp) :: thrad, phrad, r

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
  
end module elastic_tensor_tools_mod
