module vti_parameters_mod

  

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
 

  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, NSPEC_AB, myrank, &
                                   rhostore, mustore, kappastore, FOUR_THIRDS

  use specfem_par_elastic, only:   rho_vs, rho_vp, &
                                   c11store,c12store,c13store,c14store,c15store,c16store, &
                                   c22store,c23store,c24store,c25store,c26store,c33store, &
                                   c34store,c35store,c36store,c44store,c45store,c46store, &
                                   c55store,c56store,c66store

  use inverse_problem_par

  implicit none

  
  private
  public : translate_from_vti_2_cijkl, translate_from_cijkl_2_vti, translate_cijkl_gradient_2_vti

  contains 
  
!!========================================================================================================================================================

    subroutine translate_cijkl_gradient_2_vti(inversion_param, ispec, gradient)

       type(inver),                                                  intent(in)      :: inversion_param
       integer,                                                      intent(in)      :: ispec
       real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(in)      :: gradient

       integer                                                                       :: ipar, index_in_Thomsen
       !! now in only one array 
       real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_Thomsen, gradi_Thomsen

       !! first we translate from cijkl -> vti
       call cijkl_2_vti(model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                        model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))
       
       !! gradient in Thomsen
       call grad_cijkl_2_vti(model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                             model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6), &
                             gradi_Thomsen(:,:,:,1), gradi_Thomsen(:,:,:,2), gradi_Thomsen(:,:,:,3), &
                             gradi_Thomsen(:,:,:,4), gradi_Thomsen(:,:,:,5), gradi_Thomsen(:,:,:,6))
       
       !! We need to get just the inveter parameter and put them in right place
       do ipar = 1, inversion_param%NinvPar
          index_in_Thomsen = inversion_param%Index_Invert(ipar)
          gradient(:,:,:,ispec, ipar) = gradi_Thomsen(:,:,:, index_in_Thomsen)
       end do


    end subroutine translate_cijkl_gradient_2_vti

!!========================================================================================================================================================
   
    subroutine translate_from_vti_2_cijkl(inversion_param, ispec, model)
      
       type(inver),                                                  intent(in)      :: inversion_param
       integer,                                                      intent(in)      :: ispec
       real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(in)      :: model

       integer                                                                       :: ipar, index_in_Thomsen
       !! now in only one array 
       real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6)                     :: model_Thomsen
      
      
       !! first we translate from cijkl -> vti
       call cijkl_2_vti(model_Thomsen(:,:,:,1), model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                        model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))
       
       !! We need to get just the inveter parameter and put them in right place
       do ipar = 1, inversion_param%NinvPar

          index_in_Thomsen = inversion_param%Index_Invert(ipar)
          model_Thomsen(:,:,:, index_in_Thomsen) = model(:,:,:,ispec,ipar) 
          
       end do
        
       call  vti_2_cijkl(model_Thomsen(:,:,:,1) ,model_Thomsen(:,:,:,2), model_Thomsen(:,:,:,3), &
                         model_Thomsen(:,:,:,4), model_Thomsen(:,:,:,5), model_Thomsen(:,:,:,6))
     

    end subroutine translate_from_vti_2_cijkl

!!========================================================================================================================================================

    subroutine translate_from_cijkl_2_vti(inversion_param, ispec, model)

      type(inver),                                                  intent(in)      :: inversion_param
      integer,                                                      intent(in)      :: ispec
      real(kind=CUSTOM_REAL),   dimension(:,:,:,:,:), allocatable,  intent(in)      :: model

      integer 
      !! now in only one array 
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, 6), target             :: model_Thomsen
     
 
       !! first we translate from cijkl -> local model 
       call cijkl_2_vti(model_Thomsen)

       !! We need to get just the inveter parameter and put them in right place
       do ipar = 1, inversion_param%NinvPar
          index_in_Thomsen = inversion_param%Index_Invert(ipar)
          model(:,:,:,ispec, ipar) = model_Thomsen(:,:,:, index_in_Thomsen)
       end do

    end subroutine translate_from_cijkl_2_vti

!!========================================================================================================================================================
    
    subroutine vti_2_cijkl(rho, vp, vs, ep, gm, de)

      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ)          :: rho, vp, vs, ep, gm, de 
      
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

!!========================================================================================================================================================
    subroutine cijkl_2_vti(rho, vp, vs, ep, gm, de)

      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: rho, vp, vs, ep, gm, de
      
       
       rho(:,:,:,) = rhostore(:,:,:,ispec)
       vp(:,:,:)=  sqrt(c33store(:,:,:,ispec) / rhostore(:,:,:,ispec))
       vs(:,:,:)=  sqrt( c44store(:,:,:,:) / rhostore(:,:,:,ispec))
       ep(:,:,:)= (c11store(:,:,:,ispec) - c33store(:,:,:,ispec)) / 2*c33store(:,:,:,ispec)
       gm(:,:,:)= (c66store(:,:,:,ispec) - c44store(:,:,:,ispec)) / 2*c66store(:,:,:,ispec)
       de(:,:,:)= 0.5 * ( (c13store(:,:,:,ispec) + c44store(:,:,:,ispec) )**2 - (c33store(:,:,:,ispec) - c44store(:,:,:,ispec))**2 )&
            / (c33store(:,:,:,ispec)*(c33store(:,:,:,ispec) - c44store(:,:,:,ispec))) 

    end subroutine cijkl_2_vti
!!==========================================================================================================================================================


    subroutine grad_cijkl_2vti(rho, vp, vs, ep, gm, de, Grho, Gvp, Gvs, Gep, Ggm, Gde)
      
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: rho, vp, vs, ep, gm, de
      real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)            :: Grho, Gvp, Gvs, Gep, Ggm, Gde
    

      
      Grho(:,:,:) = rho_kl(:,:,:, ispec) + &
      cijkl_kl(1, :,:,:, ispec) *( 1. + 2. * ep(:,:,:)) * vp(:,:,:) **2 + &
      cijkl_kl(2, :,:,:, ispec) *( ( 1. + 2. * ep(:,:,:)) * vp(:,:,:) **2 - 2.*(1.+2.*gm(:,:,:))*vs(:,:,:)**2 ) + &
      cijkl_kl(3, :,:,:, ispec) *( sqrt( (vp(:,:,:)**2 - vs(:,:,:)**2) * ( (1.+2.*de(:,:,:))*vp(:,:,:)**2 - vs(:,:,:)***2 ) ) - vs(:,:,:)**2) + &
      cijkl_kl(7, :,:,:, ispec) *( 1. + 2. * ep(:,:,:)) * vp(:,:,:) **2 + &
      cijkl_kl(8, :,:,:, ispec) *( sqrt( (vp(:,:,:)**2 - vs(:,:,:)**2) * ( (1.+2.*de(:,:,:))*vp(:,:,:)**2 - vs(:,:,:)***2 ) ) - vs(:,:,:)**2) + &
      cijkl_kl(12, :,:,:, ispec)*vp(:,:,:)**2 + &
      cijkl_kl(16, :,:,:, ispec)*vs(:,:,:)**2 + &
      cijkl_kl(19, :,:,:, ispec)*vs(:,:,:)**2 + &
      cijkl_kl(21, :,:,:, ispec)*(1.+2.*de(:,:,:))*vs(:,:,:)**2

      Gvp(:,:,:) =  &
      cijkl_kl(1, :,:,:, ispec) *(2.*(1.+2.*ep(:,:,:)*rho(:,:,:)*vp(:,:,:)) + &
      cijkl_kl(2, :,:,:, ispec) *(2.*(1.+2.*ep(:,:,:)*rho(:,:,:)*vp(:,:,:)) + &
      cijkl_kl(3, :,:,:, ispec) *() + &
      cijkl_kl(7, :,:,:, ispec) *(2.*(1.+2.*ep(:,:,:)*rho(:,:,:)*vp(:,:,:)) + &
      cijkl_kl(8, :,:,:, ispec) *() + &
      cijkl_kl(12, :,:,:, ispec)*(2.*rho(:,:,:)*vp(:,:,:)) + &
     
      Gvs(:,:,:) =  &
      cijkl_kl(2, :,:,:, ispec) *(-4.*(1.+2.*gm(:,:,:))*rho(:,:,:)*vs(:,:,:)) + &
      cijkl_kl(3, :,:,:, ispec) *() + &
      cijkl_kl(8, :,:,:, ispec) *() + &
      cijkl_kl(16, :,:,:, ispec)*(2.*rho(:,:,:)*vs(:,:,:)) + &
      cijkl_kl(19, :,:,:, ispec)*(2.*rho(:,:,:)*vs(:,:,:)) + &
      cijkl_kl(21, :,:,:, ispec)*(2.*rho(:,:,:)*vs(:,:,:)*(1.+2.*gm(:,:,:)))

      Gep(:,:,:) =  &
      cijkl_kl(1, :,:,:, ispec) *(2.*rho(:,:,:)*vp(:,:,:)**2) + &
      cijkl_kl(2, :,:,:, ispec) *(2.*rho(:,:,:)*vp(:,:,:)**2) + &
      cijkl_kl(7, :,:,:, ispec) *(2.*rho(:,:,:)*vp(:,:,:)**2) 

      Ggm(:,:,:) =  &
      cijkl_kl(2, :,:,:, ispec) *(-4.*rho(:,:,:)*vs(:,:,:)**2) + &
      cijkl_kl(21, :,:,:, ispec)*(2.*rho(:,:,:)*vs(:,:,:)**2)


      Gde(:,:,:) =  &
      cijkl_kl(3, :,:,:, ispec) *() + &
      cijkl_kl(8, :,:,:, ispec) *() 

    end subroutine grad_cijkl_2vti

!!==========================================================================================================================================================


end module vti_parameters_mod
