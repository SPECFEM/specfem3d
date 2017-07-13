module family_parameter

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,  ANISOTROPIC_KL, ANISOTROPY

  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, NSPEC_AB, myrank, &
                                   rhostore, mustore, kappastore, FOUR_THIRDS

  use specfem_par_elastic, only: cijkl_kl, rho_kl, mu_kl, kappa_kl, ispec_is_elastic, rho_vs, rho_vp, &
                                   c11store,c12store,c13store,c14store,c15store,c16store, &
                                   c22store,c23store,c24store,c25store,c26store,c33store, &
                                   c34store,c35store,c36store,c44store,c45store,c46store, &
                                   c55store,c56store,c66store, &
                                   ispec_is_elastic, ELASTIC_SIMULATION

  use specfem_par_acoustic, only: ispec_is_acoustic, rho_ac_kl, kappa_ac_kl,  ACOUSTIC_SIMULATION

  !---------------------------------------------------------------------------------------------------------------------------------
  use inverse_problem_par

  !---- locals --
  integer,                private                                     :: ierror, ispec
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable  :: wks
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable  :: wks1

  public   ::  InvertParam2Specfem,  SpecfemParam2Invert, PrepareArraysfamilyParam, StoreGradientInfamilyParam, &
               mpi_sum_grad_all_to_all_simultaneous_runs

  !private ::

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
!> allocate working array for family paprmeter (TODO: need to be called in initialize_specfem_for_inversion())
!-------------------------------------------------------------------------
  subroutine PrepareArraysfamilyParam(inversion_param)

    type(inver),                                                  intent(inout)    :: inversion_param

    allocate(wks(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT),stat=ierror)

    if (ierror /= 0) call exit_MPI(myrank,"error allocation wks in  PrepareArraysfamilyParam subroutine, family_parameter_mod")

    if (ANISOTROPIC_KL) then
       allocate(wks1(21, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT),stat=ierror)
       if (ierror /= 0) call exit_MPI(myrank,"error allocation wks1 in  PrepareArraysfamilyParam subroutine, family_parameter_mod")
    endif

    !! manage family parameters for inversion

    !! set the number of parameters in family
    if ( ANISOTROPY) then
       inversion_param%NfamilyPar=21
    else
       inversion_param%NfamilyPar=3
       if (ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION)  inversion_param%NfamilyPar=2
    endif

    !! set the number of parameter to invert in family
    select case (trim(inversion_param%param_family))

    case('vp')
       inversion_param%NinvPar=1

    case('rho_vp', 'vp_vs','rho_kappa')

       inversion_param%NinvPar=2
    case('rho_vp_vs','rho_kappa_mu')
        inversion_param%NinvPar=3
        if (ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION)  inversion_param%NinvPar=2

    case default
       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*) '*****************'
          write(INVERSE_LOG_FILE,*) 'STOP FAMILY ', trim(inversion_param%param_family) , 'NOT KNOWN'
          write(INVERSE_LOG_FILE,*) '*****************'
       endif
       write(*,*) 'STOP FAMILY ', trim(inversion_param%param_family) , 'NOT KNOWN'
       stop

    end select


!!$    if (ANISOTROPIC_KL) then
!!$       !! TODO
!!$       write(*,*) ' ERROR : Families parameters for 21 cijkl coef not yet implemented'
!!$       stop
!!$    else
!!$       allocate(inversion_param%Index_Invert(inversion_param%NinvPar))
!!$
!!$       select case(trim(adjustl(inversion_param%param_family)))
!!$
!!$       case('rho_vp_vs')
!!$
!!$       case('rho_lambda_mu')
!!$
!!$       case('rho_ip_is')
!!$
!!$       case('rho_kappa_mu')
!!$
!!$       case default
!!$          write(*,*) 'Error : unknonwn family ',trim(adjustl(inversion_param%param_family))
!!$          stop
!!$       end select
!!$    endif

    if (myrank == 0) then

       write(INVERSE_LOG_FILE,*) '*************************************************************************'
       write(INVERSE_LOG_FILE,*) ' FAMILY PARAMETER USED ', trim(inversion_param%param_family) 
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) 
       
       if (ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)  ' Pure Acoustic : '
          write(INVERSE_LOG_FILE,*)  ' Nb inverse param :', inversion_param%NinvPar
          write(INVERSE_LOG_FILE,*)  ' Nb param :', inversion_param%NfamilyPar
          write(INVERSE_LOG_FILE,*)
       else
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)  
          write(INVERSE_LOG_FILE,*)  ' Nb inverse param :', inversion_param%NinvPar
          write(INVERSE_LOG_FILE,*)  ' Nb param :', inversion_param%NfamilyPar
       end if
       write(INVERSE_LOG_FILE,*) '*************************************************************************'
       
    endif

  end subroutine PrepareArraysfamilyParam


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------
!conversion from inversion parameters to inversion  specfem parameters
!-------------------------------------------------------------------------
  subroutine InvertParam2Specfem(inversion_param, model)
    use input_output

    type(inver),                                                  intent(in)      :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: model

    do ispec = 1, NSPEC_AB  !!
       ! elastic simulations
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test
                                          !! because we can use both elastic or acoustic elements (also poroelastic)

          if (ANISOTROPIC_KL) then

             write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                  &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                  &you just need to get parameters from array model.... in &
                  &code  : family_parameter_module.f90 subroutine InvertParam2Specfem'
             stop

             !! IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS
             !! FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem
             !! you just need to store it in array gradient ....

          else

             select case(trim(adjustl(inversion_param%param_family)))

             case('vp')
                rho_vp(:,:,:,ispec)  = model(:,:,:,ispec,1) * rhostore(:,:,:,ispec)
                kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * ( model(:,:,:,ispec,1)**2 -  &
                     FOUR_THIRDS*( rho_vs(:,:,:,ispec) / rhostore(:,:,:,ispec) )**2)

             case('rho_vp')
                rhostore(:,:,:,ispec) = model(:,:,:,ispec,1)
                rho_vp(:,:,:,ispec)   = model(:,:,:,ispec,2) * rhostore(:,:,:,ispec)
                kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * ( model(:,:,:,ispec,1)**2 -  &
                     FOUR_THIRDS*( rho_vs(:,:,:,ispec) / rhostore(:,:,:,ispec)**2))

             case('vp_vs')
                rho_vp(:,:,:,ispec)   = model(:,:,:,ispec,1) * rhostore(:,:,:,ispec)
                rho_vs(:,:,:,ispec)   = model(:,:,:,ispec,2) * rhostore(:,:,:,ispec)
                kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * ( model(:,:,:,ispec,1)**2 -  &
                     FOUR_THIRDS* model(:,:,:,ispec,2)**2 )
                mustore(:,:,:,ispec) = rhostore(:,:,:,ispec) * model(:,:,:,ispec,2)**2

             case('rho_vp_vs')
                rhostore(:,:,:,ispec) = model(:,:,:,ispec,1)
                rho_vp(:,:,:,ispec)   = model(:,:,:,ispec,2) * rhostore(:,:,:,ispec)
                rho_vs(:,:,:,ispec)   = model(:,:,:,ispec,3) * rhostore(:,:,:,ispec)
                kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * ( model(:,:,:,ispec,2)**2 -  &
                     FOUR_THIRDS* model(:,:,:,ispec,3)**2 )
                mustore(:,:,:,ispec) = rhostore(:,:,:,ispec) * model(:,:,:,ispec,3)**2

!!$             case('rho_lambda_mu')
!!$
!!$             case ('rho_ip_is')

             case('rho_kappa_mu')
                rhostore(:,:,:,ispec)   = model(:,:,:,ispec,1)
                kappastore(:,:,:,ispec) = model(:,:,:,ispec,2)
                mustore(:,:,:,ispec)    = model(:,:,:,ispec,3)
                rho_vp(:,:,:,ispec)     = sqrt( ( kappastore(:,:,:,ispec) + (4./3.) * mustore(:,:,:,ispec) ) &
                                                / rhostore(:,:,:,ispec) )  * rhostore(:,:,:,ispec)
                rho_vs(:,:,:,ispec)     = sqrt(  mustore(:,:,:,ispec) /rhostore(:,:,:,ispec) ) * rhostore(:,:,:,ispec)

             case default

                write(*,*) 'Error : unknonwn family for elastic case ',trim(adjustl(inversion_param%param_family))
                stop

             end select

          endif

       endif

       if (ispec_is_acoustic(ispec)) then

          select case(trim(adjustl(inversion_param%param_family)))

          case('vp', 'vp_vs')
             kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * model(:,:,:,ispec,1)**2

          case('rho_vp', 'rho_vp_vs')
             rhostore(:,:,:,ispec) = model(:,:,:,ispec,1)
             kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * model(:,:,:,ispec,2)**2

          case('rho_kappa', 'rho_kappa_mu')
             rhostore(:,:,:,ispec) = model(:,:,:,ispec,1)
             kappastore(:,:,:,ispec) = model(:,:,:,ispec,2)

          case default

             write(*,*) 'Error : unknonwn family for acoustic case ',trim(adjustl(inversion_param%param_family))
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
  subroutine SpecfemParam2Invert(inversion_param, model)
    type(inver),                                                  intent(in)      :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: model

    do ispec = 1, NSPEC_AB  !!
       ! elastic simulations
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test
                                          !! because we can use both elastic or acoustic elements (also proelastic)

          if (ANISOTROPIC_KL) then

             write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                  &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                  &you just need to store it in array model .... in &
                  &code  : family_parameter_module.f90 subroutine SpecfemParam2Invert'
             stop

             !! IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS
             !! FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem
             !! you just need to store it in array gradient ....

          else

             select case(trim(adjustl(inversion_param%param_family)))

             case('vp')
                model(:,:,:,ispec,1) = (kappastore(:,:,:,ispec) + (4./3.) * mustore(:,:,:,ispec) ) / rho_vp(:,:,:,ispec)

             case('rho_vp')
                model(:,:,:,ispec,1) = rho_vs(:,:,:,ispec) * rho_vs(:,:,:,ispec) / mustore(:,:,:,ispec)
                model(:,:,:,ispec,2) = (kappastore(:,:,:,ispec) + (4./3.) * mustore(:,:,:,ispec) ) / rho_vp(:,:,:,ispec)

             case('vp_vs')
                model(:,:,:,ispec,1) = rho_vs(:,:,:,ispec) * rho_vs(:,:,:,ispec) / mustore(:,:,:,ispec)
                model(:,:,:,ispec,2) = mustore(:,:,:,ispec) /  rho_vs(:,:,:,ispec)

             case('rho_vp_vs')
                model(:,:,:,ispec,1) = rho_vs(:,:,:,ispec) * rho_vs(:,:,:,ispec) / mustore(:,:,:,ispec)
                model(:,:,:,ispec,2) = (kappastore(:,:,:,ispec) + (4./3.) * mustore(:,:,:,ispec) ) / rho_vp(:,:,:,ispec)
                model(:,:,:,ispec,3) = mustore(:,:,:,ispec) /  rho_vs(:,:,:,ispec)

!!$             case('rho_lambda_mu')
!!$                model(:,:,:,ispec,1) =
!!$                model(:,:,:,ispec,2) =
!!$                model(:,:,:,ispec,3) =

!!$             case ('rho_ip_is')
!!$                model(:,:,:,ispec,1) =
!!$                model(:,:,:,ispec,2) =
!!$                model(:,:,:,ispec,3) =

             case('rho_kappa_mu')
                model(:,:,:,ispec,1) = rhostore(:,:,:,ispec)
                model(:,:,:,ispec,2) = kappastore(:,:,:,ispec)
                model(:,:,:,ispec,3) = mustore(:,:,:,ispec)


             case default

                write(*,*) 'Error : unknonwn family for elastic case ',trim(adjustl(inversion_param%param_family))
                stop

             end select

          endif

       endif

       if (ispec_is_acoustic(ispec)) then

          select case(trim(adjustl(inversion_param%param_family)))

          case('vp', 'vp_vs')
             model(:,:,:,ispec,1) = sqrt(kappastore(:,:,:,ispec) /  rhostore(:,:,:,ispec))

          case('rho_vp', 'rho_vp_vs')
             model(:,:,:,ispec,1) = rhostore(:,:,:,ispec)
             model(:,:,:,ispec,2) = sqrt(kappastore(:,:,:,ispec) /  rhostore(:,:,:,ispec))

          case('rho_kappa', 'rho_kappa_mu')
             model(:,:,:,ispec,1) = rhostore(:,:,:,ispec)
             model(:,:,:,ispec,2) = kappastore(:,:,:,ispec)

          case default

              write(*,*) 'Error : unknonwn family for acoustic case ',trim(adjustl(inversion_param%param_family))
              stop

          end select


       endif

    enddo

  end subroutine SpecfemParam2Invert

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !--------------------------------------------------------------------------
  ! store the current gradient in the inversion family parameter
  !-------------------------------------------------------------------------
  subroutine StoreGradientInfamilyParam(inversion_param, gradient)

    type(inver),                                                  intent(in)      :: inversion_param
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(inout)   :: gradient

    do ispec = 1, NSPEC_AB  !!
       ! elastic simulations
       if (ispec_is_elastic(ispec)) then  !! we need to process element by element because we need to do this test
          !! because we can use both elastic or acoustic elements (also poroelastic)
          !! get kernel in element ispec
          rho_kl(:,:,:,ispec) = - rho_kl(:,:,:,ispec) * &
               rho_vs(:,:,:,ispec) * rho_vs(:,:,:,ispec) /  mustore(:,:,:,ispec)

          if (ANISOTROPIC_KL) then

             write(*,*) ' IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS &
                  &FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem &
                  &you just need to store it in array gradient .... in &
                  &code  : family_parameter_module.f90 subroutine StoreGradientInfamilyParam'
             stop

             !! IF YOU NEED family FROM 21 ELASTIC COEFFICIENTS
             !! FEEL FREE TO DO IT ..... cijkl_kl are already computed by specfem
             !! you just need to store it in array gradient ....

          else
             mu_kl(:,:,:,ispec)  =  - 2._CUSTOM_REAL *  mustore(:,:,:,ispec) * mu_kl(:,:,:,ispec)
             kappa_kl(:,:,:,ispec) = - kappastore(:,:,:,ispec) * kappa_kl(:,:,:,ispec)

             !---------------------put kernel in the choosen family---------------------------------------------------------
             select case(trim(adjustl(inversion_param%param_family)))

             case('vp')

                !! vp
                gradient(:,:,:,ispec,1)= 2._CUSTOM_REAL * (1._CUSTOM_REAL &
                + 4._CUSTOM_REAL * mustore(:,:,:,ispec) / (3._CUSTOM_REAL * kappastore(:,:,:,ispec) ) ) * kappa_kl(:,:,:,ispec)

             case('rho_vp')

                !! rho
                gradient(:,:,:,ispec,1)=rho_kl(:,:,:,ispec) + kappa_kl(:,:,:,ispec) + mu_kl(:,:,:,ispec)

                !! vp
                gradient(:,:,:,ispec,2)=2._CUSTOM_REAL * (1._CUSTOM_REAL &
                + 4._CUSTOM_REAL * mustore(:,:,:,ispec) / (3._CUSTOM_REAL * kappastore(:,:,:,ispec) ) ) * kappa_kl(:,:,:,ispec)

             case('vp_vs')

                !! vp
                gradient(:,:,:,ispec,1)=2._CUSTOM_REAL * (1._CUSTOM_REAL &
                + 4._CUSTOM_REAL * mustore(:,:,:,ispec) / (3._CUSTOM_REAL * kappastore(:,:,:,ispec) ) ) * kappa_kl(:,:,:,ispec)

                !! vs
                gradient(:,:,:,ispec,2)= 2._CUSTOM_REAL * (mu_kl(:,:,:,ispec) &
                - 4._CUSTOM_REAL * mustore(:,:,:,ispec) / (3._CUSTOM_REAL * kappastore(:,:,:,ispec)) * kappa_kl(:,:,:,ispec))

             case('rho_vp_vs')

                !! rho
                gradient(:,:,:,ispec,1)=rho_kl(:,:,:,ispec) + kappa_kl(:,:,:,ispec) + mu_kl(:,:,:,ispec)

                !! vp
                gradient(:,:,:,ispec,2)=2._CUSTOM_REAL * (1._CUSTOM_REAL &
                + 4._CUSTOM_REAL * mustore(:,:,:,ispec) / (3._CUSTOM_REAL * kappastore(:,:,:,ispec) ) ) * kappa_kl(:,:,:,ispec)

                !!vs
                gradient(:,:,:,ispec,3)=2._CUSTOM_REAL * (mu_kl(:,:,:,ispec) &
                - 4._CUSTOM_REAL * mustore(:,:,:,ispec) / (3._CUSTOM_REAL * kappastore(:,:,:,ispec)) * kappa_kl(:,:,:,ispec))


!!$             case('rho_lambda_mu')
!!$                gradient(:,:,:,ispec,1)=
!!$                gradient(:,:,:,ispec,2)=
!!$                gradient(:,:,:,ispec,3)=
!!$
!!$             case ('rho_ip_is')
!!$                gradient(:,:,:,ispec,1)=
!!$                gradient(:,:,:,ispec,2)=
!!$                gradient(:,:,:,ispec,3)=

             case('rho_kappa_mu')
                gradient(:,:,:,ispec,1)=rho_kl(:,:,:,ispec)
                gradient(:,:,:,ispec,2)=kappa_kl(:,:,:,ispec)
                gradient(:,:,:,ispec,3)=mu_kl(:,:,:,ispec)


             case default

                write(*,*) 'Error : unknonwn family for elastic case ',trim(adjustl(inversion_param%param_family))
                stop

             end select

             !----------------------------------------------------------------------------

          endif

       endif

       if (ispec_is_acoustic(ispec)) then

          select case(trim(adjustl(inversion_param%param_family)))

          case('vp', 'vp_vs')

             !! vp
             gradient(:,:,:,ispec,1) = 2._CUSTOM_REAL * kappa_ac_kl(:,:,:,ispec)

          case('rho_vp', 'rho_vp_vs')

             !! rho
             gradient(:,:,:,ispec,1) = rho_ac_kl(:,:,:,ispec) + kappa_ac_kl(:,:,:,ispec)

             !! vp
             gradient(:,:,:,ispec,2) = 2._CUSTOM_REAL * kappa_ac_kl(:,:,:,ispec)

          case('rho_kappa', 'rho_kappa_mu')
             gradient(:,:,:,ispec,1) = rho_ac_kl(:,:,:,ispec)
             gradient(:,:,:,ispec,2) = kappa_ac_kl(:,:,:,ispec)

          case default

             write(*,*) 'Error : unknonwn family for acoustic case ',trim(adjustl(inversion_param%param_family))
             stop

          end select

       endif


    enddo


  end subroutine StoreGradientInfamilyParam



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    endif

    if (ELASTIC_SIMULATION) then

       wks(:,:,:,:)=rho_kl(:,:,:,:)
       call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), rho_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

       if (ANISOTROPIC_KL) then
          wks1(:,:,:,:,:)= cijkl_kl(:,:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks1(1,1,1,1,1), cijkl_kl(1,1,1,1,1), 21*NGLLX*NGLLY*NGLLZ*NSPEC_AB)
       else
          wks(:,:,:,:)=mu_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), mu_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          wks(:,:,:,:)=kappa_kl(:,:,:,:)
          call sum_all_all_cr_for_simulatenous_runs(wks(1,1,1,1), kappa_kl(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
       endif

    endif

    !! cost function reduction
    cost_function_tmp=cost_function
    cost_function=0.
    call sum_all_all_cr_for_simulatenous_runs(cost_function_tmp, cost_function, 1)

  end subroutine mpi_sum_grad_all_to_all_simultaneous_runs


end module family_parameter
