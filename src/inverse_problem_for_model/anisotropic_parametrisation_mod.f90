module anisotropic_parametrisation_mod

  use interpolation_mod, only: si, sp, di, dp, cp, hp
  use elastic_tensor_tools_mod, only: delta, ind_vec2tens, ind_vec2tens_voigt, &
                                      define_indexing_vec_to_tens
  implicit none

  real(kind=cp), dimension(9,21)                   :: hexa_dcij_dpref
  real(kind=cp), dimension(:,:,:,:,:), allocatable :: param_cij,  grad_cij
  real(kind=cp), dimension(:,:,:,:,:), allocatable :: param_pref, grad_pref

contains

  !================================================================================
  ! Gradients from cij to symmetry class by chain rule
  subroutine gradient_triclinic_to_isotropic_ref


  end subroutine gradient_triclinic_to_isotropic_ref

  subroutine gradient_triclinic_to_hexagonal_ref




  end subroutine gradient_triclinic_to_hexagonal_ref
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Parameters from symmetry class to cij
  subroutine parameter_isotropic_to_triclinic


  end subroutine parameter_isotropic_to_triclinic

  subroutine parameter_hexagonal_to_triclinic

    ! 1st, define parameters

  end subroutine parameter_hexagonal_to_triclinic
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Routine computing full cij tensor components from reference tti parameters
  !     param 1 : (c11,c33,c44,c66,c13,s1,s2,s3,rho)
  subroutine partial_derivative_param_ref(partial_derivative,physical_param)

    real(kind=cp), dimension (9),     intent(in) :: physical_param
    real(kind=cp), dimension (9,22), intent(out) :: partial_derivative

    integer(kind=si)            :: i, j, k, l, m, ipar
    integer(kind=si)            :: dij, dik, dil, djk, djl, dkl, dim, djm, dkm, dlm
    real(kind=cp)               :: rho, c11, c33, c44, c66, c13
    real(kind=cp)               :: si_, sj, sk, sl
    real(kind=cp)               :: sisj, sksl, sjsk, sisl, sisk, sjsl
    real(kind=cp)               :: sisjsk, sisjsl, sisksl, sjsksl, sisjsksl
    real(kind=cp), dimension(3) :: s

    !* Get direction vector and parameters
    c11    = physical_param(1)
    c33    = physical_param(2)
    c44    = physical_param(3)
    c66    = physical_param(4)
    c13    = physical_param(5)
    s(1:3) = physical_param(6:8)
    rho    = physical_param(9)

    !* Loop over cij components
    do ipar = 1, 21

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)

       !*** Precompute direction cosines
       si_ = s(i)
       sj = s(j)
       sk = s(k)
       sl = s(l)

       sisj = si_ * sj
       sksl = sk * sl
       sjsk = sj * sk
       sisl = si_ * sl
       sisk = si_ * sk
       sjsl = sj * sl

       sisjsk = sisj * sk
       sisjsl = sisj * sl
       sisksl = sisk * sl
       sjsksl = sjsk * sl

       sisjsksl = sisjsk * sl

       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)

       !*** Compute wrt vti components
       ! dcij_dc11
       partial_derivative(1,ipar) = dij*dkl - (dij*sksl + dkl*sisj) + sisjsksl

       ! dcij_dc33
       partial_derivative(2,ipar) = sisjsksl

       ! dcij_dc44
       partial_derivative(3,ipar) = dik*sjsl + dil*sjsk + djk*sisl + djl*sisk - 4._cp*sisjsksl

       ! dcij_dc66
       partial_derivative(4,ipar) = - 2._cp*dij*dkl + dik*djl + dil*djk         &
                                    + 2._dp * (dij*sk*sl + dkl*sisj)            &
                                    - (dik*sjsl + dil*sjsk + djk*sisl + djl*sisk)

       ! dcij_dc13
       partial_derivative(5,ipar) = dij*sksl + dkl*sisj - 2._cp*sisjsksl

       !*** Compute wrt direction cosines components
       ! dcij_ds1
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(6,ipar) = &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si_*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                           djk*si_*dlm + djk*dim*sl + djl*si_*dkm + djl*dim*sk) +              &
                          (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)

       ! dcij_ds2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(7,ipar) = &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si_*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                           djk*si_*dlm + djk*dim*sl + djl*si_*dkm + djl*dim*sk) +              &
                          (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)

       ! dcij_ds3
       m   = 3
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(8,ipar) =  &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si_*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                           djk*si_*dlm + djk*dim*sl + djl*si_*dkm + djl*dim*sk) +              &
                          (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)

       ! dcij_drho
       partial_derivative(9,ipar) = 0._cp

    enddo

    !*** Finally determine density partial derivatives
    partial_derivative(9,1:21) = 0._cp
    partial_derivative(9,  22) = 1._cp


  end subroutine partial_derivative_param_ref
  !--------------------------------------------------------------------------------



  !================================================================================
  !  Routine computing partial derivatives for thomsen tti parametrisation
  !     param 1 : (c33,c66,epsilon,delta,gamma,s1,s2,s3,rho)
  subroutine partial_derivative_param_edg(partial_derivative,physical_param)

    real(kind=cp), dimension (9),     intent(in) :: physical_param
    real(kind=cp), dimension (9,22), intent(out) :: partial_derivative

    integer(kind=si)            :: i, j, k, l, m, ipar
    integer(kind=si)            :: dij, dik, dil, djk, djl, dkl, dim, djm, dkm, dlm
    real(kind=cp)               :: rho, c33, c44, eps, del, gam
    real(kind=cp)               :: si_, sj, sk, sl
    real(kind=cp)               :: sisj, sksl, sjsk, sisl, sisk, sjsl
    real(kind=cp)               :: sisjsk, sisjsl, sisksl, sjsksl, sisjsksl
    real(kind=cp), dimension(3) :: s

    !* Get direction vector and parameters
    c33    = physical_param(1)
    c44    = physical_param(2)
    eps    = physical_param(3)
    del    = physical_param(4)
    gam    = physical_param(5)
    s(1:3) = physical_param(6:8)
    rho    = physical_param(9)

    !* Loop over cij components
    do ipar = 1, 21

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)

       !*** Precompute direction cosines
       si_ = s(i)
       sj = s(j)
       sk = s(k)
       sl = s(l)

       sisj = si_ * sj
       sksl = sk * sl
       sjsk = sj * sk
       sisl = si_ * sl
       sisk = si_ * sk
       sjsl = sj * sl

       sisjsk = sisj * sk
       sisjsl = sisj * sl
       sisksl = sisk * sl
       sjsksl = sjsk * sl

       sisjsksl = sisjsk * sl

       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)

       !*** Compute wrt vti components
       ! dcij_dc33
       partial_derivative(1,ipar) = dij*dkl &
                                  + 2._cp * eps * (dij*dkl - dij*sksl - dkl*sisj) + sisjsksl &
                                  +         del * (dij*sksl + dkl*sisj - 2._cp*sisjsksl)

       ! dcij_dc44
       partial_derivative(2,ipar) = (dik*djl + dil*djk - 2._cp*dij*dkl)                          &
                                  +  2._cp * gam * (-2._cp*dij*dkl + dik*djl + dil*djk           &
                                                    +2._cp*dij*sksl + 2._cp*dkl*sisj             &
                                                    - dik*sjsl - dil*sjsk - djk*sisl - djl*sisk)

       ! dcij_deps
       partial_derivative(3,ipar) = 2._cp * c33 * (dij*dkl - dij*sksl - dkl*sisj + sisjsksl)

       ! dcij_ddelta
       partial_derivative(4,ipar) = c33 * (dij*sksl + dkl*sisj - 2._cp*sisjsksl)

       ! dcij_dgamma
       partial_derivative(5,ipar) = 2._cp * c44 * (-2._cp*dij*dkl + dik*djl + dil*djk &
                                                   +2._cp*dij*sksl + 2._cp*dkl*sisj   &
                                                   - dik*sjsl - dil*sjsk - djk*sisl - djl*sisk)

       !*** Compute wrt direction cosines components
       ! dcij_ds1
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(6,ipar) = &
            + 2._cp * eps * c33 * (- dij*(sk*dlm + sl*dkm) - dkl*(si_*djm + sj*dim)              &
                                   + (sisjsk*dlm + sisjsl*dkm + sisksl*djm + sjsksl*dim))       &
            +         del * c33 * (dij*(sk*dlm + sl*dkm) + dkl*(si_*djm + sj*dim)                &
                                   - 2._cp*(sisjsk*dlm + sisjsl*dkm + sisksl*djm + sjsksl*dim)) &
            + 2._cp * gam * c44 * (+2._cp*dij*(sk*dlm + sl*dkm) + 2._cp*dkl*(si_*djm + sj*dim)   &
                                   - dik*(sj*dlm + sl*djm) - dil*(sj*dkm + sk*djm)              &
                                   - djk*(si_*dlm + sl*dim) - djl*(si_*dkm + sk*dim))

       ! dcij_ds2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(7,ipar) = &
            + 2._cp * eps * c33 * (- dij*(sk*dlm + sl*dkm) - dkl*(si_*djm + sj*dim)              &
                                   + (sisjsk*dlm + sisjsl*dkm + sisksl*djm + sjsksl*dim))       &
            +         del * c33 * (dij*(sk*dlm + sl*dkm) + dkl*(si_*djm + sj*dim)                &
                                   - 2._cp*(sisjsk*dlm + sisjsl*dkm + sisksl*djm + sjsksl*dim)) &
            + 2._cp * gam * c44 * (+2._cp*dij*(sk*dlm + sl*dkm) + 2._cp*dkl*(si_*djm + sj*dim)   &
                                   - dik*(sj*dlm + sl*djm) - dil*(sj*dkm + sk*djm)              &
                                   - djk*(si_*dlm + sl*dim) - djl*(si_*dkm + sk*dim))

       ! dcij_ds3
       m   = 3
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(8,ipar) = &
            + 2._cp * eps * c33 * (- dij*(sk*dlm + sl*dkm) - dkl*(si_*djm + sj*dim)              &
                                   + (sisjsk*dlm + sisjsl*dkm + sisksl*djm + sjsksl*dim))       &
            +         del * c33 * (dij*(sk*dlm + sl*dkm) + dkl*(si_*djm + sj*dim)                &
                                   - 2._cp*(sisjsk*dlm + sisjsl*dkm + sisksl*djm + sjsksl*dim)) &
            + 2._cp * gam * c44 * (+2._cp*dij*(sk*dlm + sl*dkm) + 2._cp*dkl*(si_*djm + sj*dim)   &
                                   - dik*(sj*dlm + sl*djm) - dil*(sj*dkm + sk*djm)              &
                                   - djk*(si_*dlm + sl*dim) - djl*(si_*dkm + sk*dim))

       ! dcij_drho
       partial_derivative(9,ipar) = 0._cp

    enddo

    !*** Finally determine density partial derivatives
    partial_derivative(9,1:21) = 0._cp
    partial_derivative(9,  22) = 1._cp


  end subroutine partial_derivative_param_edg
  !--------------------------------------------------------------------------------


 !================================================================================
  ! Routine computing full cij tensor components from reference orthorhombic parameters
  !     param 1 : (c11,c22,c33,c44,c55,c66,c23,c13,c12,a(1:3),b(1:3),c(1:3),rho)
  subroutine partial_derivative_param_ref_ortho(partial_derivative,physical_param)

    real(kind=cp), dimension (19),     intent(in) :: physical_param
    real(kind=cp), dimension (19,22), intent(out) :: partial_derivative

    integer(kind=si)            :: i, j, k, l, m, ipar
    integer(kind=si)            :: dij, dik, dil, djk, djl, dkl, dim, djm, dkm, dlm
    real(kind=cp)               :: rho, c11, c22, c33, c44, c55, c66, c23, c13, c12
    real(kind=cp)               :: ai, aj, ak, al
    real(kind=cp)               :: bi, bj, bk, bl
    real(kind=cp)               :: ci, cj, ck, cl
    real(kind=cp)               :: ajakal, bjbkbl, cjckcl
    real(kind=cp)               :: aiakal, bibkbl, cickcl
    real(kind=cp)               :: aiajal, bibjbl, cicjcl
    real(kind=cp)               :: aiajak, bibjbk, cicjck
    real(kind=cp)               :: aiajakal, bibjbkbl, cicjckcl
    real(kind=cp)               :: aiaj, akal, bibj, bkbl
    real(kind=cp)               :: aiak, aial, ajal, ajak
    real(kind=cp)               :: bibk, bibl, bjbl, bjbk
    real(kind=cp), dimension(3) :: a, b, c

    !* Get direction vector and parameters
    c11    = physical_param(1)
    c22    = physical_param(2)
    c33    = physical_param(3)
    c44    = physical_param(4)
    c55    = physical_param(5)
    c66    = physical_param(6)
    c23    = physical_param(7)
    c12    = physical_param(8)
    c13    = physical_param(9)
    a(1:3) = physical_param(10:12)
    b(1:3) = physical_param(13:15)
    c(1:3) = physical_param(16:18)
    rho    = physical_param(19)

    !* Loop over cij components
    do ipar = 1, 21

       !*** Get stiffness indexes
       i = ind_vec2tens(1,ipar)
       j = ind_vec2tens(2,ipar)
       k = ind_vec2tens(3,ipar)
       l = ind_vec2tens(4,ipar)

       !*** Precompute direction cosines
       ai = a(i);   aj = a(j);   ak = a(k);   al = a(l);
       bi = b(i);   bj = b(j);   bk = b(k);   bl = b(l);
       ci = c(i);   cj = c(j);   ck = c(k);   cl = c(l);

       aiaj = ai*aj;   akal = ak*al;   aiak = ai*ak;   aial = ai*al;
       bibj = bi*bj;   bkbl = bk*bl;   ajal = aj*al;   ajak = aj*ak;
       bibk = bi*bk;   bibl = bi*bl;   bjbl = bj*bl;   bjbk = bj*bk;

       aiajak = ai*aj*ak
       aiajal = ai*aj*al
       aiakal = ai*ak*al
       ajakal = aj*ak*al

       bibjbk = bi*bj*bk
       bibjbl = bi*bj*bl
       bibkbl = bi*bk*bl
       bjbkbl = bj*bk*bl

       cicjck = ci*cj*ck
       cicjcl = ci*cj*cl
       cickcl = ci*ck*cl
       cjckcl = cj*ck*cl

       aiajakal = ai*aj*ak*al
       bibjbkbl = bi*bj*bk*bl
       cicjckcl = ci*cj*ck*cl

       !*** Get kroneckers delta symbols
       dij = delta(i,j)
       dik = delta(i,k)
       dil = delta(i,l)
       djk = delta(j,k)
       djl = delta(j,l)
       dkl = delta(k,l)

       !*** Compute wrt vti components
       ! dcij_dc11
       partial_derivative(1,ipar) = aiajakal

       ! dcij_dc22
       partial_derivative(2,ipar) = bibjbkbl

       ! dcij_dc33
       partial_derivative(3,ipar) = cicjckcl

       ! dcij_dc44
       partial_derivative(4,ipar) = (dik*djl + dil*djk)                           &
                                  - (aiak*djl + aial*djk + ajal*dik + ajak*dil)   &
                                  + 2._dp*(aiajakal - bibjbkbl - cicjckcl)

       ! dcij_dc55
       partial_derivative(5,ipar) = (dik*djl + dil*djk)                           &
                                  - (bibk*djl + bibl*djk + bjbl*dik + bjbk*dil)   &
                                  + 2._dp*(-aiajakal + bibjbkbl - cicjckcl)

       ! dcij_dc66
       partial_derivative(6,ipar) = -(dik*djl + dil*djk)                           &
                                  + (aiak*djl + aial*djk + ajal*dik + ajak*dil)    &
                                  + (bibk*djl + bibl*djk + bjbl*dik + bjbk*dil)    &
                                  + 2._dp*(-aiajakal - bibjbkbl + cicjckcl)

       ! dcij_dc23
       partial_derivative(7,ipar) =  dij*dkl - (aiaj*dkl + dij*akal) &
                                             + (aiajakal - bibjbkbl - cicjckcl)

       ! dcij_dc13
       partial_derivative(8,ipar) =  dij*dkl - (bibj*dkl + dij*bkbl) &
                                             + (-aiajakal + bibjbkbl - cicjckcl)

       ! dcij_dc12
       partial_derivative(9,ipar) = -dij*dkl + (bibj*dkl + dij*bkbl) &
                                             + (-aiajakal - bibjbkbl + cicjckcl)


       !*** Compute wrt direction cosines components
       ! dcij_da1
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(10,ipar) = &
            + (c12 - c23) * (ai*djm*dkl + aj*dim*dkl + dij*ak*dlm + dij*al*dkm)  &
            + (c66 - c44) * (ai*djl*dkm + ak*djl*dim + al*djk*dim + ai*djk*dlm + &
            al*dik*djm + aj*dik*dim + ak*dil*djm + aj*dil*dkm) &
            + (c11 - c13 + c23 - c12 + 2._dp*( c44 - c55 - c66)) * (aiajak*dlm +  &
            aiajal*dkm + aiakal*djm + ajakal*dim)

       ! dcij_da2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(11,ipar) = &
            + (c12 - c23) * (ai*djm*dkl + aj*dim*dkl + dij*ak*dlm + dij*al*dkm) &
            + (c66 - c44) * (ai*djl*dkm + ak*djl*dim + al*djk*dim + ai*djk*dlm + &
            al*dik*djm + aj*dik*dim + ak*dil*djm + aj*dil*dkm) &
            + (c11 - c13 + c23 - c12 + 2._dp*( c44 - c55 - c66)) * &
            (aiajak*dlm +  aiajal*dkm + aiakal*djm + ajakal*dim)

       ! dcij_da3
       m   = 3
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(12,ipar) = &
            + (c12 - c23) * (ai*djm*dkl + aj*dim*dkl + dij*ak*dlm + dij*al*dkm) &
            + (c66 - c44) * (ai*djl*dkm + ak*djl*dim + al*djk*dim + ai*djk*dlm + &
            al*dik*djm + aj*dik*dim + ak*dil*djm + aj*dil*dkm) &
            + (c11 - c13 + c23 - c12 + 2._dp*( c44 - c55 - c66)) * (aiajak*dlm + &
            aiajal*dkm + aiakal*djm + ajakal*dim)

       ! dcij_db1
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(13,ipar) = &
            + (c12 - c23) * (bi*djm*dkl + bj*dim*dkl + dij*bk*dlm + dij*bl*dkm) &
            + (c66 - c55) * (bi*djl*dkm + bk*djl*dim + bl*djk*dim + bi*djk*dlm + &
            bl*dik*djm + bj*dik*dim + bk*dil*djm + bj*dil*dkm) &
            + (c22 + c13 - c23 - c12 + 2._dp*(-c44 + c55 - c66)) * (bibjbk*dlm +  &
            bibjbl*dkm + bibkbl*djm + bjbkbl*dim)

       ! dcij_db2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(14,ipar) = &
            + (c12 - c23) * (bi*djm*dkl + bj*dim*dkl + dij*bk*dlm + dij*bl*dkm) &
            + (c66 - c55) * (bi*djl*dkm + bk*djl*dim + bl*djk*dim + bi*djk*dlm + &
            bl*dik*djm + bj*dik*dim + bk*dil*djm + bj*dil*dkm) &
            + (c22 + c13 - c23 - c12 + 2._dp*(-c44 + c55 - c66)) * (bibjbk*dlm +  &
            bibjbl*dkm + bibkbl*djm + bjbkbl*dim)

       ! dcij_db3
       m   = 3
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(15,ipar) = &
            + (c12 - c23) * (bi*djm*dkl + bj*dim*dkl + dij*bk*dlm + dij*bl*dkm) &
            + (c66 - c55) * (bi*djl*dkm + bk*djl*dim + bl*djk*dim + bi*djk*dlm + &
            bl*dik*djm + bj*dik*dim + bk*dil*djm + bj*dil*dkm) &
            + (c22 + c13 - c23 - c12 + 2._dp*(-c44 + c55 - c66)) * (bibjbk*dlm +  &
            bibjbl*dkm + bibkbl*djm + bjbkbl*dim)

       ! dcij_dc1
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(16,ipar) = &
            + (c33 - c13 - c23 + c12 + 2._dp*(-c44 - c55 + c66)) * &
            (cicjck*dlm + cicjcl*dkm + cickcl*djm + cjckcl*dim)

       ! dcij_dc2
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(17,ipar) = &
            + (c33 - c13 - c23 + c12 + 2._dp*(-c44 - c55 + c66)) * &
            (cicjck*dlm + cicjcl*dkm + cickcl*djm + cjckcl*dim)

       ! dcij_dc3
       m   = 1
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)

       partial_derivative(18,ipar) = &
            + (c33 - c13 - c23 + c12 + 2._dp*(-c44 - c55 + c66)) * &
            (cicjck*dlm + cicjcl*dkm + cickcl*djm + cjckcl*dim)

       ! dcij_drho
       partial_derivative(19,ipar) = 0._cp

    enddo

    !*** Finally determine density partial derivatives
    partial_derivative(19,1:21) = 0._cp
    partial_derivative(19,  22) = 1._cp


  end subroutine partial_derivative_param_ref_ortho
  !--------------------------------------------------------------------------------


end module anisotropic_parametrisation_mod
