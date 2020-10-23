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

module vti_parameters_mod

  !! IMPORT VARIABLES FROM SPECFEM -----------------------------------------------------------------------------------------------


  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, NSPEC_AB, myrank, &
                                   rhostore, mustore, kappastore, FOUR_THIRDS

  use specfem_par_elastic, only: cijkl_kl, rho_kl, rho_vs, rho_vp, &
                                   c11store,c12store,c13store,c14store,c15store,c16store, &
                                   c22store,c23store,c24store,c25store,c26store,c33store, &
                                   c34store,c35store,c36store,c44store,c45store,c46store, &
                                   c55store,c56store,c66store

  use specfem_par_acoustic, only: rho_ac_kl, kappa_ac_kl


  use inverse_problem_par

  implicit none


  private
  public :: selector_vti_family, translate_from_vti_2_cijkl, translate_from_cijkl_2_vti, translate_cijkl_gradient_2_vti, &
            translate_from_vti_2_cijkl_ac, translate_from_cijkl_2_vti_ac, translate_cijkl_gradient_2_vti_ac

  contains

!!==========================================================================================================================

  subroutine selector_vti_family(inversion_param)
    type(inver),                                                  intent(inout)      :: inversion_param
    integer :: ipar
    integer :: ipar_inv, ier
    logical, dimension(6) :: is_selected
    character(len=MAX_LEN_STRING), dimension(6) :: vti_family_name

    vti_family_name(1)="rho"
    vti_family_name(2)="vp"
    vti_family_name(3)="vs"
    vti_family_name(4)="ep"
    vti_family_name(5)="de"
    vti_family_name(6)="gm"

    inversion_param%param_ref_name(1)="density--(rho)"
    inversion_param%param_ref_name(2)="Pwave-velocity--(vp)"
    inversion_param%param_ref_name(3)="Swave-velocity--(vs)"
    inversion_param%param_ref_name(4)="epsillon--(ep)"
    inversion_param%param_ref_name(5)="delta--(de)"
    inversion_param%param_ref_name(6)="gamma--(gm)"

    is_selected(:)=.false.
    ipar_inv=0
    inversion_param%NfamilyPar=6

    !! look for wanted parameters
    do ipar=1, inversion_param%NfamilyPar !! loop on all parameters : rho, vp, vs, ep, gm, de

       select case(trim(inversion_param%param_inv_name(ipar)))

       case('rho')
          if (.not. is_selected(1)) then
             ipar_inv=ipar_inv+1
             is_selected(1)=.true.
          endif

       case('vp')
          if (.not. is_selected(2)) then
             ipar_inv=ipar_inv+1
             is_selected(2)=.true.
          endif

       case('vs')
          if (.not. is_selected(3)) then
             ipar_inv=ipar_inv+1
             is_selected(3)=.true.
          endif

       case('ep')
          if (.not. is_selected(4)) then
             ipar_inv=ipar_inv+1
             is_selected(4)=.true.
          endif

       case('gm')
          if (.not. is_selected(5)) then
             ipar_inv=ipar_inv+1
             is_selected(5)=.true.
          endif

       case('de')
          if (.not. is_selected(6)) then
             ipar_inv=ipar_inv+1
             is_selected(6)=.true.
          endif

       end select

    enddo

    !! set wanted parameters in inversion structure
    inversion_param%NinvPar=ipar_inv
    allocate(inversion_param%Index_Invert(inversion_param%NinvPar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 563')
    ipar_inv=0
    do ipar=1, inversion_param%NfamilyPar !! loop on all parameters : rho, vp, vs, ep, gm, de
       if (is_selected(ipar)) then
          ipar_inv=ipar_inv+1
          inversion_param%Index_Invert(ipar_inv) = ipar
          inversion_param%param_inv_name(ipar_inv) = vti_family_name(ipar)
       endif
    enddo



  end subroutine selector_vti_family


!!================================================================================================================================

    subroutine translate_cijkl_gradient_2_vti(inversion_param, ispec, gradient)

       type(inver),                                                  intent(in)      :: inversion_param
       integer,                                                      intent(in)      :: ispec
       real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   allocatable,  intent(inout)   :: gradient

       integer                                                                       :: ipar, index_in_Thomsen
       !! now in only one array
       real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_Thomsen, gradi_Thomsen

       !! first we translate from cijkl -> vti
       call cijkl_2_vti(ispec, model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                        model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))

       !! gradient in Thomsen
       call grad_cijkl_2_vti(ispec, model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                             model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6), &
                             gradi_Thomsen(:,:,:,1), gradi_Thomsen(:,:,:,2), gradi_Thomsen(:,:,:,3), &
                             gradi_Thomsen(:,:,:,4), gradi_Thomsen(:,:,:,5), gradi_Thomsen(:,:,:,6))

       !! We need to get just the inveter parameter and put them in right place
       do ipar = 1, inversion_param%NinvPar
          index_in_Thomsen = inversion_param%Index_Invert(ipar)
          gradient(:,:,:,ipar) = gradi_Thomsen(:,:,:, index_in_Thomsen)
       enddo


    end subroutine translate_cijkl_gradient_2_vti

!!================================================================================================================================

    subroutine translate_from_vti_2_cijkl(inversion_param, ispec, model)

       type(inver),                                                  intent(in)      :: inversion_param
       integer,                                                      intent(in)      :: ispec
       real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   allocatable,  intent(in)      :: model

       integer                                                                       :: ipar, index_in_Thomsen
       !! now in only one array
       real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_Thomsen


       !! first we translate from cijkl -> vti
       call cijkl_2_vti(ispec, model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                        model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))

       !! We need to get just the inveter parameter and put them in right place
       do ipar = 1, inversion_param%NinvPar
          index_in_Thomsen = inversion_param%Index_Invert(ipar)
          model_Thomsen(:,:,:, index_in_Thomsen) = model(:,:,:,ipar)
       enddo

       call  vti_2_cijkl(ispec, model_Thomsen(:,:,:,1) ,model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                         model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))


    end subroutine translate_from_vti_2_cijkl

!!================================================================================================================================

    subroutine translate_from_cijkl_2_vti(inversion_param, ispec, model)

      type(inver),                                                  intent(in)      :: inversion_param
      integer,                                                      intent(in)      :: ispec
      real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   allocatable,  intent(inout)   :: model

      integer                                                                       :: ipar, index_in_Thomsen
      !! now in only one array
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_Thomsen


      !! first we translate from cijkl -> vti
      call cijkl_2_vti(ispec, model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                       model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))

      !! We need to get just the inveter parameter and put them in right place
      do ipar = 1, inversion_param%NinvPar
         index_in_Thomsen = inversion_param%Index_Invert(ipar)
         model(:,:,:, ipar) = model_Thomsen(:,:,:, index_in_Thomsen)
      enddo

    end subroutine translate_from_cijkl_2_vti

!!================================================================================================================================

    subroutine vti_2_cijkl(ispec, rho, vp, vs, ep, gm, de)

      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp, vs, ep, gm, de
      integer,                                                      intent(in)      :: ispec

      !! now put new model directly in specfem database
      rhostore(:,:,:,ispec) =  rho(:,:,:)
      rho_vp(:,:,:,ispec) = rho(:,:,:)* vp(:,:,:)
      rho_vs(:,:,:,ispec) = rho(:,:,:)* vs(:,:,:)

      kappastore(:,:,:,ispec) =   rho(:,:,:)*( vp(:,:,:)**2  - FOUR_THIRDS* vs(:,:,:)**2 )
      mustore(:,:,:,ispec) =      rho(:,:,:) * vs(:,:,:)**2

      c11store(:,:,:,ispec) = (1._CUSTOM_REAL+2._CUSTOM_REAL*ep(:,:,:))* rho(:,:,:) * vp(:,:,:) * vp(:,:,:)

      c12store(:,:,:,ispec) = -2.*(2.*gm(:,:,:)+1.)*vs(:,:,:)*vs(:,:,:)*rho(:,:,:) + &
           (1._CUSTOM_REAL+2._CUSTOM_REAL*ep(:,:,:))* rho(:,:,:) * vp(:,:,:) * vp(:,:,:)


      c13store(:,:,:,ispec) =  rho(:,:,:) * (sqrt( (vp(:,:,:)**2-vs(:,:,:)**2)*&
           ((1.+2.*de(:,:,:))*vp(:,:,:)**2-vs(:,:,:)**2)) - vs(:,:,:)**2)

      c14store(:,:,:,ispec) =  0._CUSTOM_REAL
      c15store(:,:,:,ispec) =  0._CUSTOM_REAL
      c16store(:,:,:,ispec) =  0._CUSTOM_REAL

      c22store(:,:,:,ispec) = (1._CUSTOM_REAL+2._CUSTOM_REAL*ep(:,:,:))* rho(:,:,:) * vp(:,:,:) * vp(:,:,:)

      c23store(:,:,:,ispec) = rho(:,:,:) * ( sqrt((vp(:,:,:)**2-vs(:,:,:)**2)*&
           ((1.+2.*de(:,:,:))*vp(:,:,:)**2-vs(:,:,:)**2)) - vs(:,:,:)**2)

      c24store(:,:,:,ispec) =  0._CUSTOM_REAL
      c25store(:,:,:,ispec) =  0._CUSTOM_REAL
      c26store(:,:,:,ispec) =  0._CUSTOM_REAL

      c33store(:,:,:,ispec) = rho(:,:,:) *vp(:,:,:) *vp(:,:,:)

      c34store(:,:,:,ispec) =  0._CUSTOM_REAL
      c35store(:,:,:,ispec) =  0._CUSTOM_REAL
      c36store(:,:,:,ispec) =  0._CUSTOM_REAL

      c44store(:,:,:,ispec) = rho(:,:,:)*vs(:,:,:)*vs(:,:,:)

      c45store(:,:,:,ispec) =  0._CUSTOM_REAL
      c46store(:,:,:,ispec) =  0._CUSTOM_REAL

      c55store(:,:,:,ispec) = rho(:,:,:)*vs(:,:,:)*vs(:,:,:)

      c56store(:,:,:,ispec) =  0._CUSTOM_REAL

      c66store(:,:,:,ispec) = (2*gm(:,:,:)+1)*vs(:,:,:)*vs(:,:,:)*rho(:,:,:)

    end subroutine vti_2_cijkl

!!================================================================================================================================
    subroutine cijkl_2_vti(ispec, rho, vp, vs, ep, gm, de)

      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: rho, vp, vs, ep, gm, de
      integer,                                                      intent(in)      :: ispec

      rho(:,:,:) = rhostore(:,:,:,ispec)
      vp(:,:,:)=  sqrt(c33store(:,:,:,ispec) / rhostore(:,:,:,ispec))
      vs(:,:,:)=  sqrt( c44store(:,:,:,ispec) / rhostore(:,:,:,ispec))
      ep(:,:,:)= (c11store(:,:,:,ispec) - c33store(:,:,:,ispec)) / (2.*c33store(:,:,:,ispec))
      gm(:,:,:)= (c66store(:,:,:,ispec) - c44store(:,:,:,ispec)) / (2.*c44store(:,:,:,ispec))
      de(:,:,:)= 0.5 * ( (c13store(:,:,:,ispec) + c44store(:,:,:,ispec) )**2 - (c33store(:,:,:,ispec) - c44store(:,:,:,ispec))**2 )&
           / (c33store(:,:,:,ispec)*(c33store(:,:,:,ispec) - c44store(:,:,:,ispec)))

    end subroutine cijkl_2_vti
!!================================================================================================================================


    subroutine grad_cijkl_2_vti(ispec, rho, vp, vs, ep, gm, de, Grho, Gvp, Gvs, Gep, Ggm, Gde)

      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: rho, vp, vs, ep, gm, de
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: Grho, Gvp, Gvs, Gep, Ggm, Gde
      integer,                                                      intent(in)      :: ispec


      !! tmp work arrays
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: w1, w2, w3, w4, w5
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: rh_vp, rh_vs

      !! finalize specfem kernel need to multiply by -1
      rho_kl(:,:,:, ispec) = - rho_kl(:,:,:, ispec)
      cijkl_kl(:,:,:,:,ispec) = - cijkl_kl(:,:,:,:,ispec)


      !! preconpute temp arrays
      rh_vp(:,:,:) = rho(:,:,:)*vp(:,:,:)
      rh_vs(:,:,:) = rho(:,:,:)*vs(:,:,:)

      w1(:,:,:) = sqrt( (vp(:,:,:)**2 - vs(:,:,:)**2) * ( (1.+2.*de(:,:,:))*vp(:,:,:)**2 - vs(:,:,:)**2 ))
      w2(:,:,:) = (1. + 2. *ep(:,:,:))
      w3(:,:,:) = (1. + 2. *gm(:,:,:))
      w4(:,:,:) = (1. + 2.*de(:,:,:)) *vp(:,:,:)**2 - (1.+de(:,:,:))*vs(:,:,:)**2
      w5(:,:,:) = (1. +    de(:,:,:))* vp(:,:,:)**2 - vs(:,:,:)**2

      !! RHO
      Grho(:,:,:) = rho_kl(:,:,:, ispec) + &
      cijkl_kl(1, :,:,:, ispec) * (w2(:,:,:) * vp(:,:,:)**2) + &
      cijkl_kl(2, :,:,:, ispec) * (w2(:,:,:) * vp(:,:,:)**2 - 2.*w3(:,:,:)*vs(:,:,:)**2 ) + &
      cijkl_kl(3, :,:,:, ispec) * (w1(:,:,:) - vs(:,:,:)**2) + &
      cijkl_kl(7, :,:,:, ispec) * (w2(:,:,:) * vp(:,:,:)**2) + &
      cijkl_kl(8, :,:,:, ispec) * (w1(:,:,:) - vs(:,:,:)**2) + &
      cijkl_kl(12, :,:,:, ispec)* vp(:,:,:)**2 + &
      cijkl_kl(16, :,:,:, ispec)* vs(:,:,:)**2 + &
      cijkl_kl(19, :,:,:, ispec)* vs(:,:,:)**2 + &
      cijkl_kl(21, :,:,:, ispec)*(w3(:,:,:)*vs(:,:,:)**2)

      !! VP
      Gvp(:,:,:) =  &
      cijkl_kl(1, :,:,:, ispec) *(2.*w2(:,:,:)*rh_vp(:,:,:)) + &
      cijkl_kl(2, :,:,:, ispec) *(2.*w2(:,:,:)*rh_vp(:,:,:)) + &
      cijkl_kl(3, :,:,:, ispec) *(2. * rh_vp(:,:,:) * w4(:,:,:) / w1(:,:,:) ) + &
      cijkl_kl(7, :,:,:, ispec) *(2.*w2(:,:,:)*rh_vp(:,:,:)) + &
      cijkl_kl(8, :,:,:, ispec) *(2. * rh_vp(:,:,:) * w4(:,:,:) / w1(:,:,:) ) + &
      cijkl_kl(12, :,:,:, ispec)*(2.*rh_vp(:,:,:))

      !! VS
      Gvs(:,:,:) =  &
      cijkl_kl(2, :,:,:, ispec) *(-4.*w3(:,:,:)*rh_vs(:,:,:)) + &
      cijkl_kl(3, :,:,:, ispec) *(-2.*rh_vs(:,:,:) * (1+ w5(:,:,:) / w1(:,:,:))) + &
      cijkl_kl(8, :,:,:, ispec) *(-2.*rh_vs(:,:,:) * (1+ w5(:,:,:) / w1(:,:,:))) + &
      cijkl_kl(16, :,:,:, ispec)*(2.*rh_vs(:,:,:)) + &
      cijkl_kl(19, :,:,:, ispec)*(2.*rh_vs(:,:,:)) + &
      cijkl_kl(21, :,:,:, ispec)*(2.*rh_vs(:,:,:)*w3(:,:,:))

      !! EPSILLON
      Gep(:,:,:) =  &
      cijkl_kl(1, :,:,:, ispec) *(2.*rh_vp(:,:,:)*vp(:,:,:)) + &
      cijkl_kl(2, :,:,:, ispec) *(2.*rh_vp(:,:,:)*vp(:,:,:)) + &
      cijkl_kl(7, :,:,:, ispec) *(2.*rh_vp(:,:,:)*vp(:,:,:))

      !! GAMMA
      Ggm(:,:,:) =  &
      cijkl_kl(2, :,:,:, ispec) *(-4.*rh_vs(:,:,:)*vs(:,:,:)) + &
      cijkl_kl(21, :,:,:, ispec)*( 2.*rh_vs(:,:,:)*vs(:,:,:))

      !! DELTA
      Gde(:,:,:) =  &
      ( cijkl_kl(3, :,:,:, ispec) +  cijkl_kl(8, :,:,:, ispec) ) &
      *( 2.* rho(:,:,:) * (vp(:,:,:)**2) * ( vp(:,:,:)**2 - vs(:,:,:)**2 ) / w1(:,:,:) )

    end subroutine grad_cijkl_2_vti

!!================================================================================================================================
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! ===================================== ACOUSTIC CASE ========================================

    subroutine translate_cijkl_gradient_2_vti_ac(inversion_param, ispec, gradient)

      type(inver),                                                  intent(in)      :: inversion_param
      integer,                                                      intent(in)      :: ispec
      real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   allocatable,  intent(inout)   :: gradient

      integer                                                                       :: ipar, index_in_iso
      !! now in only one array
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_iso, gradi_iso

      !! first we translate from cijkl -> vti
      call cijkl_2_vti_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2))

      !! gradient in Thomsen
      gradi_iso(:,:,:,:) = 0._CUSTOM_REAL
      call grad_cijkl_2_vti_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2), &
                         gradi_iso(:,:,:,1), gradi_iso(:,:,:,2))

      !! We need to get just the inveter parameter and put them in right place
      gradient(:,:,:,:) = 0._CUSTOM_REAL
      do ipar = 1, 2
         index_in_iso = inversion_param%Index_Invert(ipar)
         gradient(:,:,:, ipar) = gradi_iso(:,:,:, index_in_iso)
      enddo

    end subroutine translate_cijkl_gradient_2_vti_ac

!!================================================================================================================================

  subroutine translate_from_vti_2_cijkl_ac(inversion_param, ispec, model)

    type(inver),                                                  intent(in)      :: inversion_param
    integer,                                                      intent(in)      :: ispec
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:),   allocatable,  intent(in)      :: model

    integer                                                                       :: ipar, index_in_iso
    !! full inv parametrization
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_iso

    model_iso(:,:,:,:)=0._CUSTOM_REAL
    !! (modeling -> inv)
    call cijkl_2_vti_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2))

    !! We need to get just the inverter parameter and put them in right place
    do ipar = 1, inversion_param%NinvPar
       index_in_iso = inversion_param%Index_Invert(ipar)
       model_iso(:,:,:, index_in_iso) = model(:,:,:,ipar)
    enddo

    call  vti_2_cijkl_ac(ispec, model_iso(:,:,:,1) ,model_iso(:,:,:,2))


  end subroutine translate_from_vti_2_cijkl_ac

!!================================================================================================================================

  subroutine translate_from_cijkl_2_vti_ac(inversion_param, ispec, model)

      type(inver),                                                  intent(in)      :: inversion_param
      integer,                                                      intent(in)      :: ispec
      real(kind=CUSTOM_REAL),   dimension(:,:,:,:), allocatable,    intent(inout)   :: model

      integer                                                                       :: ipar, index_in_iso
      !! now in only one array
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_iso

      model_iso(:,:,:,:) = 0._CUSTOM_REAL
      !! modeling -> inv
      call cijkl_2_vti_ac(ispec, model_iso(:,:,:,1), model_iso(:,:,:,2))

      !! We need to get just the inverter parameter and put them in right place
      do ipar = 1, inversion_param%NinvPar
         index_in_iso = inversion_param%Index_Invert(ipar)
         model(:,:,:,ipar) = model_iso(:,:,:, index_in_iso)
      enddo

    end subroutine translate_from_cijkl_2_vti_ac
!!================================================================================================================================

    subroutine vti_2_cijkl_ac(ispec, rho, vp)
      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp
      integer,                                                      intent(in)      :: ispec
      rhostore(:,:,:,ispec)  =  rho(:,:,:)
      kappastore(:,:,:,ispec) =   rho(:,:,:)* vp(:,:,:)**2
    end subroutine vti_2_cijkl_ac

!!================================================================================================================================

    subroutine cijkl_2_vti_ac(ispec, rho, vp)
      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp
      integer,                                                      intent(in)      :: ispec
      rho(:,:,:) = rhostore(:,:,:,ispec)
      vp(:,:,:)  = sqrt(kappastore(:,:,:,ispec)  / rhostore(:,:,:,ispec))
    end subroutine cijkl_2_vti_ac

!!================================================================================================================================

    subroutine grad_cijkl_2_vti_ac(ispec, rho, vp, Grho, Gvp)
      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp
      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: Grho, Gvp
      integer,                                                      intent(in)      :: ispec

      !! put specfem kernel in **not** log
      rho_ac_kl(:,:,:,ispec) =  rho_ac_kl(:,:,:,ispec) / rhostore(:,:,:,ispec)
      kappa_ac_kl(:,:,:,ispec) = kappa_ac_kl(:,:,:,ispec) / kappastore(:,:,:,ispec)

      Grho(:,:,:) = rho_ac_kl(:,:,:,ispec) + kappa_ac_kl(:,:,:,ispec)*vp(:,:,:)**2
      Gvp(:,:,:) = 2._CUSTOM_REAL * kappa_ac_kl(:,:,:,ispec)*rho(:,:,:)*vp(:,:,:)

    end subroutine grad_cijkl_2_vti_ac


end module vti_parameters_mod
