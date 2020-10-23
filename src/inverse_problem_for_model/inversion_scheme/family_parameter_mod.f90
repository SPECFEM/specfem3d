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

module family_parameter

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,  ANISOTROPIC_KL, ANISOTROPY, APPROXIMATE_HESS_KL, &
    ELASTIC_SIMULATION,ACOUSTIC_SIMULATION

  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, NSPEC_AB, myrank

  use specfem_par_elastic, only: cijkl_kl, rho_kl, mu_kl, kappa_kl, ispec_is_elastic, hess_rho_kl, &
       hess_mu_kl, hess_kappa_kl

  use specfem_par_acoustic, only: hess_kappa_ac_kl, hess_rho_ac_kl, rho_ac_kl, kappa_ac_kl, ispec_is_acoustic

  !---------------------------------------------------------------------------------------------------------------------------------
  use inverse_problem_par

  !! one module per family parameter for inversion
  use iso_parameter_mod
  use vti_parameters_mod


  implicit none


  private

  !---- locals --
  integer,                private                                     :: ier, ispec
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable  :: wks
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable  :: wks1
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable  :: gradient_wks
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable  :: model_wks
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable  :: model_ref_wks

  public   ::  InvertParam2Specfem,  SpecfemParam2Invert, PrepareArraysfamilyParam, StoreGradientInfamilyParam, &
               mpi_sum_grad_all_to_all_simultaneous_runs, SpecfemPrior2Invert, Modeling2RefInvert

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
!> allocate working array for family paprmeter (TODO: need to be called in initialize_specfem_for_inversion())
!-------------------------------------------------------------------------
  subroutine PrepareArraysfamilyParam(inversion_param)

    type(inver),                                                  intent(inout)    :: inversion_param
    integer                                                                        :: i

    !! temporary array useful for MPI comm
    allocate(wks(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 564')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks in  PrepareArraysfamilyParam subroutine, family_parameter_mod")
    if (ANISOTROPIC_KL) then
       allocate(wks1(21, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 565')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks1 in  PrepareArraysfamilyParam subroutine, family_parameter_mod")
    endif

    !! manage family parameters for inversion
    call choose_inversion_parameters(inversion_param)

    !! temporay arrays used for translation : inversion parmeters <-> modeling parameters
    allocate(gradient_wks(NGLLX, NGLLY, NGLLZ, inversion_param%NinvPar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 566')
    allocate(model_wks(NGLLX, NGLLY, NGLLZ, inversion_param%NinvPar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 567')
    allocate(model_ref_wks(NGLLX, NGLLY, NGLLZ, inversion_param%NinvPar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 568')


    if (myrank == 0) then

       write(INVERSE_LOG_FILE,*) '*************************************************************************'
       write(INVERSE_LOG_FILE,*) ' FAMILY PARAMETER USED ', trim(inversion_param%parameter_family_name)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)

       if (ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)  ' Pure Acoustic : '
       endif
       write(INVERSE_LOG_FILE,*)  ' Nb inverse param :', inversion_param%NinvPar
       write(INVERSE_LOG_FILE,*)  ' Nb param in full family:', inversion_param%NfamilyPar
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '  full family parameters :'
       write(INVERSE_LOG_FILE,*)
       do i=1, inversion_param%NfamilyPar
          write(INVERSE_LOG_FILE,*) ' Param ', i, ' :  ' , trim(inversion_param%param_ref_name(i))
       enddo
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       do i=1, inversion_param%NinvPar
          write(INVERSE_LOG_FILE,*) ' inverse parameter ', i, ' : ',trim(inversion_param%param_inv_name(i))
       enddo

    endif

  end subroutine PrepareArraysfamilyParam

!!=========================================================================================================================

  subroutine choose_inversion_parameters(inversion_param)
    type(inver),               intent(inout)      :: inversion_param

    select case (inversion_param%parameter_family_name)

    case('VTI')
       call selector_vti_family(inversion_param)

    case('ISO')
       call selector_iso_family(inversion_param)

    case('IP')

    case('LAME')

    case default
      write(*,*)  " inversion parameters not known ", inversion_param%parameter_family_name
      stop
   end select

 end subroutine choose_inversion_parameters



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
!conversion from inversion parameters to inversion  specfem parameters
!-------------------------------------------------------------------------

  subroutine InvertParam2Specfem(inversion_param, model, model_ref)
    use input_output, only: create_mass_matrices_Stacey_duplication_routine

    type(inver),                                                  intent(in)      :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(in)      :: model, model_ref


    do ispec = 1, NSPEC_AB  !!


       !! get model in locals arrays memory ---------------------------------------------------------------------------------
       !! and set in physical units
       select case(inversion_param%parameter_metric)
       case(0)  !! directly the parameter P
          model_wks(:,:,:,:) = model(:,:,:,ispec,:)

       case(1)  !! P / Pref
          model_wks(:,:,:,:) = model(:,:,:,ispec,:) * model_ref(:,:,:,ispec,:)

       case(2)  !! log(P)
          model_wks(:,:,:,:) = exp( model(:,:,:,ispec,:))

       case(3) !! log(P/Pref)
          !!??
       end select

       !! todo : call change_metric_2units(ispec, model, model_wks, model_ref_wks)

       ! elastic simulations
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test
                                          !! because we can use both elastic or acoustic elements (also poroelastic)

          if (ANISOTROPIC_KL) then

             select case (trim(inversion_param%parameter_family_name))
             case("VTI")
                call translate_from_vti_2_cijkl(inversion_param, ispec, model_wks)

             case default
                write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                     &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                     &you just need to get parameters from array model.... in &
                     &code  : family_parameter_module.f90 subroutine InvertParam2Specfem'
                stop

             end select

          else

             select case(trim(adjustl(inversion_param%parameter_family_name)))
             case('ISO')
                call translate_from_iso_2_lame(inversion_param, ispec, model_wks)

             case('LAME')

             case('IP')

             case default
                write(*,*) 'Error : unknonwn family for elastic case ', &
                     trim(adjustl(inversion_param%parameter_family_name))
                stop

             end select

          endif

       endif

       if (ispec_is_acoustic(ispec)) then

          select case(trim(adjustl(inversion_param%parameter_family_name)))

          case('VTI')
             call translate_from_vti_2_cijkl_ac(inversion_param, ispec, model_wks)

          case('ISO')
             call translate_from_iso_2_lame_ac(inversion_param, ispec, model_wks)

          case('LAME')

          case ('IP')

          case default
             write(*,*) 'Error : unknown family for acoustic case ', &
                  trim(adjustl(inversion_param%parameter_family_name))
             stop

          end select

       endif

    enddo

    !! re-compute mass matrix  (this routine is in mesh_tool module)
    call create_mass_matrices_Stacey_duplication_routine()

  end subroutine InvertParam2Specfem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
! conversion from specfem parameters to inversion parameters
!-------------------------------------------------------------------------

  subroutine SpecfemParam2Invert(inversion_param, model, model_ref)
    type(inver),                                                  intent(in)      :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(in)      :: model_ref
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: model

    do ispec = 1, NSPEC_AB  !!

       !! set model ref in local array
       model_ref_wks(:,:,:,:) = model_ref(:,:,:,ispec,:)

       ! elastic simulations
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test
                                          !! because we can use both elastic or acoustic elements (also proelastic)

          if (ANISOTROPIC_KL) then

             select case(trim(inversion_param%parameter_family_name))
             case ("VTI")
                call translate_from_cijkl_2_vti(inversion_param, ispec, model_wks)

             case default
                write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                     &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                     &you just need to store it in array model .... in &
                     &code  : family_parameter_module.f90 subroutine SpecfemParam2Invert'
                stop
             end select

          else

             select case (trim(inversion_param%parameter_family_name))
             case('ISO')
                call translate_from_lame_2_iso(inversion_param, ispec, model_wks)

             case('LAME')

             case('IP')

             case default
                write(*,*) 'Error : unknonwn family for elastic case ', &
                     trim(adjustl(inversion_param%parameter_family_name))
                stop
             end select

          endif

       endif

       if (ispec_is_acoustic(ispec)) then

          select case (trim(inversion_param%parameter_family_name))
          case ("VTI")
             call translate_from_cijkl_2_vti_ac(inversion_param, ispec, model_wks)

          case('ISO')
             call translate_from_lame_2_iso_ac(inversion_param, ispec, model_wks)

          case('LAME')

          case('IP')

          case default
             write(*,*) 'Error : unknonwn family for elastic case ', &
                  trim(adjustl(inversion_param%parameter_family_name))
             stop
          end select

       endif


       !! store local array model in global array and set the metric for inversion  ------------------------------------------
       select case(inversion_param%parameter_metric)
       case(0)  !! directly the parameter P
          model(:,:,:,ispec,:) = model_wks(:,:,:,:)
       case(1)  !! P / Pref
         model(:,:,:,ispec,:) = model_wks(:,:,:,:) / model_ref_wks(:,:,:,:)
       case(2)  !! log(P)
          model(:,:,:,ispec,:) = log(model_wks(:,:,:,:))
       case(3)  !! log(P/Pref)
          !! ??
          model(:,:,:,ispec,:) = log(model_wks(:,:,:,:) / model_ref_wks(:,:,:,:))
       end select
       !! todo : call call change_model_units_2metric(ispec, model, model_wks, model_ref_wks)
    enddo

  end subroutine SpecfemParam2Invert

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
! store the current gradient in the inversion family parameter
!-------------------------------------------------------------------------

  subroutine StoreGradientInfamilyParam(inversion_param, gradient, model, model_ref, hess_approxim)

    type(inver),                                                  intent(in)      :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(in)      :: model, model_ref
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: gradient, hess_approxim



    do ispec = 1, NSPEC_AB  !! loop on elements


       !! get model in loacals arrays memory ---------------------------------------------------------------------------------
       !! and set in physical units
       select case(inversion_param%parameter_metric)
       case(0)  !! directly the parameter P
          model_wks(:,:,:,:) = model(:,:,:,ispec,:)
          model_ref_wks(:,:,:,:) = model_ref(:,:,:,ispec,:)
       case(1)  !! P / Pref
          model_wks(:,:,:,:) = model(:,:,:,ispec,:) * model_ref(:,:,:,ispec,:)
          model_ref_wks(:,:,:,:) = model_ref(:,:,:,ispec,:)
       case(2)  !! log(P)
          model_wks(:,:,:,:) = exp( model(:,:,:,ispec,:))
          model_ref_wks(:,:,:,:) = model_ref(:,:,:,ispec,:)
       case(3) !! log(P/Pref)
          !!??
       end select

       ! elastic element -------------------------------------------------------------------------------------------------------
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test

          if (ANISOTROPIC_KL) then  !! Anisotropic

             select case (trim(inversion_param%parameter_family_name))
             case ("VTI")

                call translate_cijkl_gradient_2_vti(inversion_param ,ispec, gradient_wks)

             case default
                write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                     &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                     &you just need to store it in array gradient .... in &
                     &code  : family_parameter_module.f90 subroutine StoreGradientInfamilyParam'
                stop
             end select

          else  !! isotropic

             select case (trim(inversion_param%parameter_family_name))
             case("ISO")
                call translate_lame_gradient_2_iso(inversion_param, ispec, gradient_wks)

             case("LAME")

             case("IP")

             case default
                write(*,*) 'Error : unknonwn family for elastic case ',trim(adjustl(inversion_param%parameter_family_name))
                stop
             end select

             !! TODO store preconditionner kernels  ( Hessian ... may be not necessary )

          endif

       endif

       ! acoustic element ---------------------------------------------------------------------------------------------------
       if (ispec_is_acoustic(ispec)) then

          select case (trim(inversion_param%parameter_family_name))
          case("VTI")
             call translate_cijkl_gradient_2_vti_ac(inversion_param ,ispec, gradient_wks)

          case("ISO")
             call translate_lame_gradient_2_iso_ac(inversion_param, ispec, gradient_wks)

          case("LAME")

          case("IP")

          case default
             write(*,*) hess_approxim(1,1,1,1,1)
             write(*,*) 'Error : unknonwn family for elastic case ',trim(adjustl(inversion_param%parameter_family_name))
             stop
          end select

          !! TODO store preconditionner kernels  ( Hessian ... may be not necessary )

       endif


       !! store local array gradient in global array and set the metric for inversion  ------------------------------------------
       select case(inversion_param%parameter_metric)
       case(0)  !! directly the parameter P
          gradient(:,:,:,ispec,:) = gradient_wks(:,:,:,:)
       case(1)  !! P / Pref
          gradient(:,:,:,ispec,:) = gradient_wks(:,:,:,:) * model_ref_wks(:,:,:,:)
       case(2)  !! log(P)
          gradient(:,:,:,ispec,:) = gradient_wks(:,:,:,:) * model_wks(:,:,:,:)
       case(3)  !! log(P/Pref)
          !! ??
       end select
       ! todo  : call change_gradient_unit_2metric(ispec, gradient, gradient_wks, model_wks,  model_ref_wks)
    enddo   !! loop on elements

  end subroutine StoreGradientInfamilyParam

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
! conversion from specfem parameters to inversion parameters  !! TODO this is currently not working : fix it
!-------------------------------------------------------------------------

  subroutine SpecfemPrior2Invert(inversion_param, model)
    type(inver),                                                  intent(inout)   :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: model
    integer :: ipar
    write(*,*) " ABORT : OPTION READING DIRECLTY PRIOR MODEL NOT TESTED "
    stop
    do ispec = 1, NSPEC_AB  !!
       do ipar =1 ,  inversion_param%NinvPar
          model(:,:,:,ispec, ipar) =  inversion_param%prior_model(:,:,:,ispec, inversion_param%Index_Invert(ipar))
       enddo
    enddo

    !! clear memory since we do not need it any more
    deallocate(inversion_param%prior_model)

  end subroutine SpecfemPrior2Invert

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------
! conversion from specfem parameters to inversion parameters in physical unit metric
!-------------------------------------------------------------------------------------

  subroutine Modeling2RefInvert(inversion_param, model)
    type(inver),                                                  intent(inout)   :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: model
    !integer :: ipar


    do ispec = 1, NSPEC_AB  !!


       ! elastic simulations
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test
                                          !! because we can use both elastic or acoustic elements (also proelastic)
          if (ANISOTROPIC_KL) then

             select case(trim(inversion_param%parameter_family_name))
             case ("VTI")
                call translate_from_cijkl_2_vti(inversion_param, ispec, model_wks)

             case default
                write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                     &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                     &you just need to store it in array model .... in &
                     &code  : family_parameter_module.f90 subroutine SpecfemParam2Invert'
                stop
             end select

          else

             select case (trim(inversion_param%parameter_family_name))
             case('ISO')
                call translate_from_lame_2_iso(inversion_param, ispec, model_wks)

             case('LAME')

             case('IP')

             case default
                write(*,*) 'Error : unknonwn family for elastic case ', &
                     trim(adjustl(inversion_param%parameter_family_name))
                stop
             end select

          endif

       endif

       if (ispec_is_acoustic(ispec)) then

          select case (trim(inversion_param%parameter_family_name))
          case ("VTI")
             call translate_from_cijkl_2_vti_ac(inversion_param, ispec, model_wks)

          case('ISO')
             call translate_from_lame_2_iso_ac(inversion_param, ispec, model_wks)

          case('LAME')

          case('IP')

          case default
             write(*,*) 'Error : unknonwn family for elastic case ', &
                  trim(adjustl(inversion_param%parameter_family_name))
             stop
          end select

       endif

       !! set model in physical unit
       model(:,:,:,ispec,:) = model_wks(:,:,:,:)

    enddo

  end subroutine Modeling2RefInvert


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                  routine puerly for MPI communcation of gradients (simultaneous runs)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine mpi_sum_grad_all_to_all_simultaneous_runs(cost_function)

    real(kind=CUSTOM_REAL), intent(inout) :: cost_function
    real(kind=CUSTOM_REAL)                :: cost_function_tmp

    if ( NUMBER_OF_SIMULTANEOUS_RUNS < 2 ) return !! not need of MPI commutication througth simultaneous sources

    !! communicate kernels from each simulataneaou run
    if (ACOUSTIC_SIMULATION) then

       wks(:,:,:,:)=rho_ac_kl(:,:,:,:)
       call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), rho_ac_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
       wks(:,:,:,:)=kappa_ac_kl(:,:,:,:)
       call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), kappa_ac_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

       if (APPROXIMATE_HESS_KL) then
          wks(:,:,:,:)=hess_rho_ac_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), hess_rho_ac_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          wks(:,:,:,:)=hess_kappa_ac_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), hess_kappa_ac_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
       endif

    endif

    if (ELASTIC_SIMULATION) then

       wks(:,:,:,:)=rho_kl(:,:,:,:)
       call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), rho_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

       if (APPROXIMATE_HESS_KL) then
          wks(:,:,:,:)=hess_rho_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), hess_rho_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
       endif


       if (ANISOTROPIC_KL) then
          wks1(:,:,:,:,:)= cijkl_kl(:,:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks1(1,1,1,1,1), cijkl_kl(1,1,1,1,1), 21*NGLLX*NGLLY*NGLLZ*NSPEC_AB)
       else
          wks(:,:,:,:)=mu_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), mu_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          wks(:,:,:,:)=kappa_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), kappa_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

          if (APPROXIMATE_HESS_KL) then
             wks(:,:,:,:)=hess_mu_kl(:,:,:,:)
             call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), hess_mu_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
             wks(:,:,:,:)=hess_kappa_kl(:,:,:,:)
             call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), hess_kappa_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          endif
       endif

    endif

    !! cost function reduction
    cost_function_tmp=cost_function
    cost_function=0.
    call sum_all_all_cr_for_simulatenous_runs(cost_function_tmp, cost_function, 1)

  end subroutine mpi_sum_grad_all_to_all_simultaneous_runs


!!==========================================================================================================================

end module family_parameter
