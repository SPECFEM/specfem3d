!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

module iso_parameter_mod

  !! IMPORT VARIABLES FROM SPECFEM -----------------------------------------------------------------------------------------------


  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, NSPEC_AB, myrank, &
                                   rhostore, mustore, kappastore, FOUR_THIRDS

  use specfem_par_elastic, only: rho_kl,  mu_kl, kappa_kl, rho_vs, rho_vp
  use specfem_par_acoustic, only: rho_ac_kl, kappa_ac_kl

  use inverse_problem_par

  implicit none

  private
  public :: selector_iso_family, translate_from_iso_2_lame, translate_from_lame_2_iso, translate_lame_gradient_2_iso, &
            translate_from_iso_2_lame_ac, translate_from_lame_2_iso_ac, translate_lame_gradient_2_iso_ac

contains

!!================================================================================================================================

  subroutine selector_iso_family(inversion_param)

    type(inver),                                                  intent(inout)      :: inversion_param
    ! local
    integer :: ipar
    integer :: ipar_inv, ier
    logical, dimension(3) :: is_selected
    character(len=MAX_LEN_STRING), dimension(3) :: vti_family_name

    vti_family_name(1) = "rho"
    vti_family_name(2) = "vp"
    vti_family_name(3) = "vs"

    inversion_param%param_ref_name(1) = "density--(rho)"
    inversion_param%param_ref_name(2) = "Pwave-velocity--(vp)"
    inversion_param%param_ref_name(3) = "Swave-velocity--(vs)"

    is_selected(:) = .false.
    ipar_inv = 0
    inversion_param%NfamilyPar = 3

    !! look for wanted parameters
    do ipar = 1, inversion_param%NfamilyPar !! loop on all parameters : rho, vp, vs

       select case(trim(inversion_param%param_inv_name(ipar)))

       case('rho')
          if (.not. is_selected(1)) then
             ipar_inv = ipar_inv+1
             is_selected(1) = .true.
          endif

       case('vp')
          if (.not. is_selected(2)) then
             ipar_inv = ipar_inv+1
             is_selected(2) = .true.
          endif

       case('vs')
          if (.not. is_selected(3)) then
             ipar_inv = ipar_inv+1
             is_selected(3) = .true.
          endif

       end select

    enddo

    !! set wanted parameters in inversion structure
    inversion_param%NinvPar = ipar_inv
    allocate(inversion_param%Index_Invert(inversion_param%NinvPar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 569')
    inversion_param%Index_Invert(:) = 0

    ipar_inv = 0
    do ipar = 1, inversion_param%NfamilyPar !! loop on all parameters : rho, vp, vs
       if (is_selected(ipar)) then
          ipar_inv = ipar_inv+1
          inversion_param%Index_Invert(ipar_inv) = ipar
          inversion_param%param_inv_name(ipar_inv) = vti_family_name(ipar)
       endif
    enddo

  end subroutine selector_iso_family

!!================================================================================================================================

  subroutine translate_lame_gradient_2_iso(inversion_param, ispec, gradient)

    type(inver),                                  intent(in)        :: inversion_param
    integer,                                      intent(in)        :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:), intent(inout)     :: gradient

    integer                                                                       :: ipar, index_in_iso
    !!
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 3)                     :: model_iso, gradi_iso

    !! first we translate from cijkl -> iso
    call lame_2_iso(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2), model_iso(:,:,:,3))

    !! gradient in iso
    call grad_lame_2_iso(ispec,model_iso(:,:,:,1), model_iso(:,:,:,2), model_iso(:,:,:,3), &
                         gradi_iso(:,:,:,1), gradi_iso(:,:,:,2), gradi_iso(:,:,:,3))

    !! store just the gradient for inversible parameters
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       gradient(:,:,:,ipar) = gradi_iso(:,:,:, index_in_iso)
    enddo

  end subroutine translate_lame_gradient_2_iso


!!================================================================================================================================

  subroutine translate_from_iso_2_lame(inversion_param, ispec, model)

    type(inver),                                    intent(in)      :: inversion_param
    integer,                                        intent(in)      :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   intent(in)      :: model

    integer                                                         :: ipar, index_in_iso
    !! full inv parametrization
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 3)       :: model_iso

    !! (modeling -> inv)
    call lame_2_iso(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2), model_iso(:,:,:,3))

    !! We need to get just the inveter parameter and put them in right place
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       model_iso(:,:,:, index_in_iso) = model(:,:,:,ipar)
    enddo

    call iso_2_lame(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2), model_iso(:,:,:,3))

  end subroutine translate_from_iso_2_lame

!!================================================================================================================================

  subroutine translate_from_lame_2_iso(inversion_param, ispec, model)

    type(inver),                                    intent(in)      :: inversion_param
    integer,                                        intent(in)      :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   intent(inout)   :: model

    integer                                                         :: ipar, index_in_iso
    !! now in only one array
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 3)       :: model_iso

    !! modeling -> inv
    call lame_2_iso(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2), model_iso(:,:,:,3))

    !! We need to get just the inverted parameters and put them in right place
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       model(:,:,:,ipar) = model_iso(:,:,:, index_in_iso)
    enddo

  end subroutine translate_from_lame_2_iso

!!================================================================================================================================

  subroutine iso_2_lame(ispec, rho, vp, vs)

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp, vs
    integer                                                       :: ispec

    !! put new model directly in specfem database
    rhostore(:,:,:,ispec)   = rho(:,:,:)
    rho_vp(:,:,:,ispec)     = rho(:,:,:) * vp(:,:,:)
    rho_vs(:,:,:,ispec)     = rho(:,:,:) * vs(:,:,:)
    kappastore(:,:,:,ispec) = rho(:,:,:) * ( vp(:,:,:)**2  - FOUR_THIRDS * vs(:,:,:)**2 )
    mustore(:,:,:,ispec)    = rho(:,:,:) * vs(:,:,:)**2

  end subroutine iso_2_lame

!!================================================================================================================================

  subroutine lame_2_iso(ispec, rho, vp, vs)

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp, vs
    integer                                                       :: ispec

    rho(:,:,:) = rhostore(:,:,:,ispec)
    vp(:,:,:)  = (kappastore(:,:,:,ispec) + (4./3.) * mustore(:,:,:,ispec) ) / rho_vp(:,:,:,ispec)
    vs(:,:,:)  = mustore(:,:,:,ispec) /  rho_vs(:,:,:,ispec)

  end subroutine lame_2_iso

!!================================================================================================================================

  subroutine grad_lame_2_iso(ispec, rho, vp, vs, Grho, Gvp, Gvs)

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp, vs
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: Grho, Gvp, Gvs
    integer                                                       :: ispec

    !! finalize specfem kernels
    rho_kl(:,:,:,ispec)   = - rho_kl(:,:,:, ispec)
    mu_kl(:,:,:,ispec)    = - 2. *  mu_kl(:,:,:,ispec)
    kappa_kl(:,:,:,ispec) = - kappa_kl(:,:,:,ispec)

    Grho(:,:,:) =  rho_kl(:,:,:,ispec) + &
                   vs(:,:,:)**2 * mu_kl(:,:,:,ispec) + &
                  (vp(:,:,:)**2 - (4./3.) * vs(:,:,:)**2) * kappa_kl(:,:,:, ispec)

    Gvp(:,:,:)  = 2. * rho(:,:,:) * vp(:,:,:) * kappa_kl(:,:,:, ispec)

    Gvs(:,:,:)  = (-8./3.) * rho(:,:,:) * vs(:,:,:) * kappa_kl(:,:,:,ispec) + &
                  2. * rho(:,:,:) * vs(:,:,:) * mu_kl(:,:,:,ispec)

  end subroutine grad_lame_2_iso
!!================================================================================================================================

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! ===================================== ACOUSTIC CASE ========================================

  subroutine translate_lame_gradient_2_iso_ac(inversion_param, ispec, gradient)

    type(inver),                                    intent(in)      :: inversion_param
    integer,                                        intent(in)      :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   intent(inout)   :: gradient

    integer                                                         :: ipar, index_in_iso
    !! now in only one array
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 3)       :: model_iso, gradi_iso

    !! first we translate from cijkl -> vti
    call lame_2_iso_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2))

    !! gradient in Thomsen
    gradi_iso(:,:,:,:) = 0._CUSTOM_REAL
    ! 1 == Grho, 2 == Gvp
    call grad_lame_2_iso_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2), &
                            gradi_iso(:,:,:,1), gradi_iso(:,:,:,2))

    !! We need to get just the inverted parameters and put them in right place
    gradient(:,:,:, :) = 0._CUSTOM_REAL
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       gradient(:,:,:, ipar) = gradi_iso(:,:,:, index_in_iso)
    enddo

  end subroutine translate_lame_gradient_2_iso_ac

!!================================================================================================================================

  subroutine translate_from_iso_2_lame_ac(inversion_param, ispec, model)

    type(inver),                                    intent(in)      :: inversion_param
    integer,                                        intent(in)      :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   intent(in)      :: model

    integer                                                         :: ipar, index_in_iso
    !! full inv parametrization
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 3)       :: model_iso

    model_iso(:,:,:,:) = 0._CUSTOM_REAL
    !! (modeling -> inv)
    call lame_2_iso_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2))

    !! We need to get just the inverted parameters and put them in right place
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       model_iso(:,:,:, index_in_iso) = model(:,:,:,ipar)
    enddo

    call iso_2_lame_ac(ispec, model_iso(:,:,:,1) ,model_iso(:,:,:,2))

  end subroutine translate_from_iso_2_lame_ac

!!================================================================================================================================

  subroutine translate_from_lame_2_iso_ac(inversion_param, ispec, model)

    type(inver),                                    intent(in)      :: inversion_param
    integer,                                        intent(in)      :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   intent(inout)   :: model

    integer                                                         :: ipar, index_in_iso
    !! now in only one array
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 3)       :: model_iso

    model_iso(:,:,:,:) = 0._CUSTOM_REAL
    !! modeling -> inv
    call lame_2_iso_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2))

    !! We need to get just the inverted parameters and put them in right place
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       model(:,:,:, ipar) = model_iso(:,:,:, index_in_iso)
    enddo

  end subroutine translate_from_lame_2_iso_ac
!!================================================================================================================================

  subroutine iso_2_lame_ac(ispec, rho, vp)

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp
    integer                                                       :: ispec

    rhostore(:,:,:,ispec)   = rho(:,:,:)
    kappastore(:,:,:,ispec) = rho(:,:,:) * vp(:,:,:)**2

  end subroutine iso_2_lame_ac

!!================================================================================================================================

  subroutine lame_2_iso_ac(ispec, rho, vp)

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp
    integer                                                       :: ispec

    rho(:,:,:) = rhostore(:,:,:,ispec)
    vp(:,:,:)  = sqrt(kappastore(:,:,:,ispec) / rhostore(:,:,:,ispec))

  end subroutine lame_2_iso_ac

!!================================================================================================================================

  subroutine grad_lame_2_iso_ac(ispec, rho, vp, Grho, Gvp)

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: Grho, Gvp
    integer                                                       :: ispec

    !! put specfem kernel in **not** log (absolute kernels, rather than relative ones)
    !
    ! note: rho_ac_kl(:,:,:,ispec) contributions are still positive and without material factors up to this point.
    !       only in save_adjoint_kernels.f90 they would be added.
    !       thus, no need to divide by rhostore(:,:,:,ispec) or kappastore(:,:,:,ispec), but need to add minus sign

    !! finalize specfem kernel need to multiply by -1
    rho_ac_kl(:,:,:,ispec)   = - rho_ac_kl(:,:,:,ispec)
    kappa_ac_kl(:,:,:,ispec) = - kappa_ac_kl(:,:,:,ispec)

    ! gradients
    Grho(:,:,:) = rho_ac_kl(:,:,:,ispec) + kappa_ac_kl(:,:,:,ispec) * vp(:,:,:)**2
    Gvp(:,:,:)  = 2._CUSTOM_REAL * kappa_ac_kl(:,:,:,ispec) * rho(:,:,:) * vp(:,:,:)

  end subroutine grad_lame_2_iso_ac

end module iso_parameter_mod
