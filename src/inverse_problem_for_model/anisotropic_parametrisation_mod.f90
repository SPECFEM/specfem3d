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

    integer(kind=si)            :: i, j, k, l, m, ipar
    integer(kind=si)            :: dij, dik, dil, djk, djl, dkl, dim, djm, dkm, dlm
    real(kind=cp)               :: rho, c11, c33, c44, c66, c13
    real(kind=cp)               :: si, sj, sk, sl
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
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)
       
       ! dcij_ds2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(7,ipar) = &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)
       
       ! dcij_ds3
       m   = 3 
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(8,ipar) =  &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)

       ! dcij_drho
       partial_derivative(9,ipar) = 0._cp

    end do

    !*** Finally determine density partial derivatives
    partial_derivative(9,1:21) = 0._cp
    partial_derivative(9,  22) = 1._cp
    

  end subroutine partial_derivative_param_ref
  !--------------------------------------------------------------------------------


  
  !================================================================================
  ! Routine computing partial derivatives for reference tti parametrisation
  !     param 1 : (c11,c33,c44,c66,c13,s1,s2,s3,rho)
  subroutine partial_derivative_param_ref(partial_derivative,physical_param)

    integer(kind=si)            :: i, j, k, l, m, ipar
    integer(kind=si)            :: dij, dik, dil, djk, djl, dkl, dim, djm, dkm, dlm
    real(kind=cp)               :: rho, c11, c33, c44, c66, c13
    real(kind=cp)               :: si, sj, sk, sl
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
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)
       
       ! dcij_ds2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(7,ipar) = &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)
       
       ! dcij_ds3
       m   = 3 
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(8,ipar) =  &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)

       ! dcij_drho
       partial_derivative(9,ipar) = 0._cp

    end do

    !*** Finally determine density partial derivatives
    partial_derivative(9,1:21) = 0._cp
    partial_derivative(9,  22) = 1._cp
    

  end subroutine partial_derivative_param_ref
  !--------------------------------------------------------------------------------


  !================================================================================
  !  Routine computing partial derivatives for thomsen tti parametrisation
  !     param 1 : (c33,c66,eps,delta,gamma,s1,s2,s3,rho)
  subroutine partial_derivative_param_edg(partial_derivative,physical_param)

    integer(kind=si)            :: i, j, k, l, m, ipar
    integer(kind=si)            :: dij, dik, dil, djk, djl, dkl, dim, djm, dkm, dlm
    real(kind=cp)               :: rho, c33, c66, eps, delta, gamma
    real(kind=cp)               :: si, sj, sk, sl
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
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)
       
       ! dcij_ds2
       m   = 2
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(7,ipar) = &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)
       
       ! dcij_ds3
       m   = 3 
       dim = delta(i,m)
       djm = delta(j,m)
       dkm = delta(k,m)
       dlm = delta(l,m)
       partial_derivative(8,ipar) =  &
            (c13 - c11 + 2._cp*c66) *  (dij*sk*dlm + dij*dkm*sl + dkl*dim*sj + dkl*si*djm) + &
            (c44 - c66) * (dik*sj*dlm + dik*djm*sl + dil*sj*dkm + dil*djm*sk +               &
                        +  djk*si*dlm + djk*dim*sl + djl*si*dkm + djl*dim*sk) +              &
                        + (c11 + c33 - 2._cp*c13 - 4._cp*c44) * (sisjsk*dlm + sisjsl*dkm   + &
                                                                 sisksl*djm + sjsksl*dim)

       ! dcij_drho
       partial_derivative(9,ipar) = 0._cp

    end do

    !*** Finally determine density partial derivatives
    partial_derivative(9,1:21) = 0._cp
    partial_derivative(9,  22) = 1._cp
    

  end subroutine partial_derivative_param_edg
  !--------------------------------------------------------------------------------

  
end module anisotropic_parametrisation_mod
