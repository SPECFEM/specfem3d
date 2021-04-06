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

module  inversion_scheme

  !! IMPORT VARIABLES FROM SPECFEM
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, myrank

  use inverse_problem_par

  implicit none

  integer,                private                                        ::  ier
  integer,                private                                        ::  Ninvpar, Mbfgs

  !! size(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar, Mbfgs
  real(kind=CUSTOM_REAL), private,  dimension(:,:,:,:,:,:), allocatable  ::  bfgs_stored_gradient
  real(kind=CUSTOM_REAL), private,  dimension(:,:,:,:,:,:), allocatable  ::  bfgs_stored_model

  !! for l-bfgs working arrays
  real(kind=CUSTOM_REAL), private,  dimension(:,:,:,:,:),   allocatable  ::  wks_1, wks_2
  real(kind=CUSTOM_REAL), private,  dimension(0:1000)                    ::  ak_store, pk_store

  !! for normalization
  real(kind=8), private,  dimension(:,:,:,:,:),   allocatable  ::  wks_1n, wks_2n

  public  :: AllocateArraysForInversion, DeAllocateArraysForInversion, wolfe_rules, ComputeDescentDirection
  private :: L_BFGS_GENERIC


contains


!--------------------------------------------------------------------------
!> allocate arrays for inversion
!-------------------------------------------------------------------------

  subroutine AllocateArraysForInversion(inversion_param)

    type(inver),                           intent(in)      :: inversion_param

    Ninvpar=inversion_param%NinvPar
    Mbfgs=inversion_param%max_history_bfgs

    allocate(bfgs_stored_gradient(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar, 0:Mbfgs),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 557')
    if (ier /= 0) call exit_MPI(myrank,"error allocation bfgs_stored_gradient in AllocateArraysForInversion subroutine")

    allocate(bfgs_stored_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar, 0:Mbfgs),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 558')
    if (ier /= 0) call exit_MPI(myrank,"error allocation bfgs_stored_model in AllocateArraysForInversion subroutine")

    allocate(wks_1(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 559')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_1 in AllocateArraysForInversion subroutine")

    allocate(wks_2(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 560')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_2 in AllocateArraysForInversion subroutine")

    allocate(wks_1n(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 561')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_1n in AllocateArraysForInversion subroutine")

    allocate(wks_2n(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 562')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_2n in AllocateArraysForInversion subroutine")

    bfgs_stored_gradient(:,:,:,:,:,:) = 0._CUSTOM_REAL
    bfgs_stored_model(:,:,:,:,:,:) = 0._CUSTOM_REAL
    wks_1(:,:,:,:,:) = 0._CUSTOM_REAL
    wks_2(:,:,:,:,:) = 0._CUSTOM_REAL
    wks_1n(:,:,:,:,:) = 0._CUSTOM_REAL
    wks_2n(:,:,:,:,:) = 0._CUSTOM_REAL


  end subroutine AllocateArraysForInversion


!----------------------------------------------------------------------------------------------------------------------------------

  subroutine DeAllocateArraysForInversion()

    deallocate(bfgs_stored_gradient)
    deallocate(bfgs_stored_model)
    deallocate(wks_1)
    deallocate(wks_2)

  end subroutine DeAllocateArraysForInversion


!----------------------------------------------------------------------------------------------------------------------------------

  subroutine ComputeDescentDirection(current_iteration, descent_direction, fwi_precond)

    real(kind=CUSTOM_REAL),  dimension(:,:,:,:,:),   allocatable, intent(inout)  :: descent_direction, fwi_precond
    integer,                                                      intent(in)     :: current_iteration
    integer                                                                      :: ilast, ipar

    ilast = mod(current_iteration, Mbfgs)

    if (current_iteration == 0 .or. USE_GRADIENT_OPTIM) then

       do ipar = 1, Ninvpar
          descent_direction(:,:,:,:,ipar) =  -fwi_precond(:,:,:,:,ipar) *  bfgs_stored_gradient(:,:,:,:,ipar,ilast)
       enddo

    else

       call L_BFGS_GENERIC(current_iteration, descent_direction, fwi_precond)

    endif

  end subroutine ComputeDescentDirection


!--------------------------------------------------------------------------
!> l-bfgs generic routine
!-------------------------------------------------------------------------

  subroutine L_BFGS_GENERIC(current_iteration, descent_direction, fwi_precond)

    real(kind=CUSTOM_REAL),  dimension(:,:,:,:,:),   allocatable, intent(inout)  :: descent_direction, fwi_precond
    integer,                                                      intent(in)     :: current_iteration
    integer                                                                      :: k, imin, imax
    real(kind=CUSTOM_REAL)                                                       :: ak, pk, norme_yiter, beta

    !! define boundary loop
    imin = 0
    if (current_iteration >= Mbfgs) then
       imax = Mbfgs-1
    else
       imax = current_iteration - 1
    endif

    if (VERBOSE_MODE .and. myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '  Calling l-bfgs, iter  :', current_iteration
       write(INVERSE_LOG_FILE,*) '  Loops :',imin, imax
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
    endif

    !! initialize L-BFGS
    ak_store(:)=0._CUSTOM_REAL
    pk_store(:)=0._CUSTOM_REAL
    descent_direction(:,:,:,:,:) = bfgs_stored_gradient(:,:,:,:,:,imax+1)

    wks_1(:,:,:,:,:) = bfgs_stored_gradient(:,:,:,:,:,imax+1) -  bfgs_stored_gradient(:,:,:,:,:,imax)
    call Parallel_ComputeL2normSquare(wks_1, Ninvpar, norme_yiter)


    do k = imax, imin, -1

       wks_1(:,:,:,:,:) = bfgs_stored_gradient(:,:,:,:,:,k+1) -  bfgs_stored_gradient(:,:,:,:,:,k)
       wks_2(:,:,:,:,:) = bfgs_stored_model(:,:,:,:,:,k+1)    -  bfgs_stored_model(:,:,:,:,:,k)

       ! Liu and Nocedal 1989 (Mathematical programming)
       ! Diagonal_prec(:,:,:,:,:) =  Diagonal_prec(:,:,:,:,:) +  (wks_1(:,:,:,:,:)*wks_2(:,:,:,:,:)) /(wks_1(:,:,:,:,:)**2)
       !

       call Parallel_ComputeInnerProduct(wks_2, wks_1, Ninvpar, pk)
       pk_store(k) = 1._CUSTOM_REAL / pk

       call Parallel_ComputeInnerProduct(wks_2, descent_direction, Ninvpar, ak)
       ak_store(k) = pk_store(k) * ak

       descent_direction(:,:,:,:,:) = descent_direction(:,:,:,:,:) - ak_store(k) * wks_1(:,:,:,:,:)

    enddo

    !! Nocedal's default preconditionning
    k=imax
    pk = 1._CUSTOM_REAL / (pk_store(k) * norme_yiter)
    descent_direction(:,:,:,:,:) = pk * descent_direction(:,:,:,:,:)

    !! customer preconditionning
    descent_direction(:,:,:,:,:) =  fwi_precond(:,:,:,:,:) * descent_direction(:,:,:,:,:)

    !! diagonal precoditionner
    !!descent_direction(:,:,:,:,:) =  Diagonal_prec(:,:,:,:,:) *  descent_direction(:,:,:,:,:)

    do k = imin, imax

       wks_1(:,:,:,:,:) = bfgs_stored_gradient(:,:,:,:,:,k+1) -  bfgs_stored_gradient(:,:,:,:,:,k)
       wks_2(:,:,:,:,:) = bfgs_stored_model(:,:,:,:,:,k+1)    -  bfgs_stored_model(:,:,:,:,:,k)

       call Parallel_ComputeInnerProduct(wks_1, descent_direction, Ninvpar, beta)
       beta = pk_store(k) * beta

       descent_direction(:,:,:,:,:) = descent_direction(:,:,:,:,:) + &
            (ak_store(k) - beta) *  wks_2(:,:,:,:,:)

    enddo

    descent_direction(:,:,:,:,:) = -1._CUSTOM_REAL * descent_direction(:,:,:,:,:)

  end subroutine L_BFGS_GENERIC


!----------------------------------------------------------------------------------------------------------

  subroutine wolfe_rules(mwl1, mwl2, q0, qt, qp0, vqpt, step_length, td, tg, flag_wolfe)

    implicit none
    real(kind=CUSTOM_REAL), intent(in)    :: mwl1, mwl2, q0, qt, qp0, vqpt
    real(kind=CUSTOM_REAL), intent(inout) :: step_length, td, tg
    logical,                intent(inout) :: flag_wolfe
    real(kind=CUSTOM_REAL)                :: qpt

    if (DEBUG_MODE) write(IIDD,*) ' WOLFE RULES ', mwl1, mwl2, q0, qt, &
         qp0, vqpt, step_length, td, tg, flag_wolfe

    if (myrank == 0)  write(INVERSE_LOG_FILE,*)

    qpt = (qt - q0) / step_length

    if (mwl2*qp0 <= vqpt .and. qpt <= mwl1*qp0 ) then
       flag_wolfe=.true.
       if (myrank == 0) write(INVERSE_LOG_FILE,*) '   --- > Wolfe  rules :  step accepted '
       return
    endif

    if (mwl1*qp0 < qpt) then
       td=step_length
       if (myrank == 0) write(INVERSE_LOG_FILE,*) '   --- > Wolfe rules :  right step update'
    endif

    if (qpt <= mwl1*qp0 .and. vqpt < mwl2 *qp0) then
       tg=step_length
       if (myrank == 0) write(INVERSE_LOG_FILE,*) '   --- > Wolfe rules:  left step update '
    endif

    if (td == 0.) then
       step_length = 2.*step_length
       if (myrank == 0) write(INVERSE_LOG_FILE,*) '   --- > Wolfe rules :  step too small '
    else
       step_length = 0.5*(td+tg)
       if (myrank == 0) write(INVERSE_LOG_FILE,*) '   --- > Wolfe rules :  step too big '
    endif
    if (myrank == 0)  write(INVERSE_LOG_FILE,*)

  end subroutine wolfe_rules


!-------------------------------------------------------------------------------------------------

  subroutine StoreModelAndGradientForLBFGS(models_to_store, gradients_to_store, iteration_to_store)

    integer,                                                   intent(in) :: iteration_to_store
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(in) :: models_to_store, gradients_to_store
    integer                                                               :: k

    if (iteration_to_store <= Mbfgs ) then

       if (VERBOSE_MODE .and. myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '  Storing iterations in l-bfgs memory ', iteration_to_store
          write(INVERSE_LOG_FILE,*) '  Total  iteration to store :', Mbfgs
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
       endif

       bfgs_stored_model(:,:,:,:,:,iteration_to_store) =  models_to_store(:,:,:,:,:)
       bfgs_stored_gradient(:,:,:,:,:,iteration_to_store) = gradients_to_store(:,:,:,:,:)

    else

       if (VERBOSE_MODE .and. myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '  Storing iterations in l-bfgs memory ', iteration_to_store
          write(INVERSE_LOG_FILE,*) '  Shifting previous arrays and store in ', Mbfgs, ' index '
          write(INVERSE_LOG_FILE,*) '  Total  iteration to store :', Mbfgs
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
       endif

       do k = 0, Mbfgs-1
          bfgs_stored_model(:,:,:,:,:,k) = bfgs_stored_model(:,:,:,:,:,k+1)
          bfgs_stored_gradient(:,:,:,:,:,k) = bfgs_stored_gradient(:,:,:,:,:,k+1)
       enddo

       bfgs_stored_model(:,:,:,:,:,Mbfgs) =    models_to_store(:,:,:,:,:)
       bfgs_stored_gradient(:,:,:,:,:,Mbfgs) = gradients_to_store(:,:,:,:,:)

    endif

  end subroutine StoreModelAndGradientForLBFGS


!---------------------------------------------------------------------------------------------

  subroutine Parallel_ComputeInnerProduct(vect1, vect2, Niv, qp)

    use specfem_par, only: NSPEC_AB,  jacobianstore, wxgll, wygll, wzgll

    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(in)    :: vect1, vect2
    real(kind=CUSTOM_REAL),                                    intent(inout) :: qp
    integer,                                                   intent(in)    :: Niv
    real(kind=CUSTOM_REAL)                                                   :: qp_tmp_single
    ! try double precision
    real(kind=8)                                                             :: jacobianl, weight, qp_tmp
    integer                                                                  :: ipar, i, j, k, ispec
    real(kind=CUSTOM_REAL)                                                   :: coeff, coeff_n1, coeff_n2
    real(kind=8)                                                             :: coeff_n1_dp, coeff_n2_dp

    !! try normalization to avoid numerical errors
    !call Parallel_ComputeL2normSquare(vect1 , Niv, coeff_n1)
    !call Parallel_ComputeL2normSquare(vect2 , Niv, coeff_n2)

    coeff=maxval(abs(vect1(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n1)
    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL
    wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1

    coeff=maxval(abs(vect2(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n2)
    if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2=1._CUSTOM_REAL
    wks_2n(:,:,:,:,:) = vect2(:,:,:,:,:) / coeff_n2

    coeff_n1_dp = coeff_n1
    coeff_n2_dp = coeff_n2

    qp_tmp=0._CUSTOM_REAL

    do ipar=1, Niv
       do ispec = 1, NSPEC_AB

          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   weight = wxgll(i)*wygll(j)*wzgll(k)
                   jacobianl = jacobianstore(i,j,k,ispec)
                   qp_tmp = qp_tmp + jacobianl * weight * wks_1n(i,j,k,ispec,ipar) * wks_2n(i,j,k,ispec,ipar)
                   !qp = qp + jacobianl * weight * vect1(i,j,k,ispec,ipar) * vect2(i,j,k,ispec,ipar)
                enddo
             enddo
          enddo
       enddo
    enddo

    qp_tmp_single = qp_tmp * coeff_n1_dp * coeff_n2_dp
    qp=0.
    call sum_all_all_cr(qp_tmp_single, qp)

  end subroutine Parallel_ComputeInnerProduct


!---------------------------------------------------------------------------------------------

   subroutine Parallel_ComputeL2normSquare(vect1 , Niv, qp)

     use specfem_par, only: NSPEC_AB,  jacobianstore, wxgll, wygll, wzgll

     real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(in)    :: vect1
     real(kind=CUSTOM_REAL),                                    intent(inout) :: qp
     integer,                                         intent(in)    :: Niv
     real(kind=CUSTOM_REAL)                                         :: coeff, coeff_n1
     real(kind=CUSTOM_REAL)                                         :: qp_tmp
     real(kind=8)                                                   :: jacobianl, weight, qp_dp, coeff_n1_dp
     integer                                                        :: ipar, i, j, k, ispec

     qp=0.d0
     qp_dp=0.d0

     coeff=maxval(abs(vect1(:,:,:,:,:)))
     call max_all_all_cr(coeff, coeff_n1)

     if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL

     wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1
     coeff_n1_dp=coeff_n1

     do ipar=1,Niv
        do ispec = 1, NSPEC_AB

           do k=1,NGLLZ
              do j=1,NGLLY
                 do i=1,NGLLX
                    weight = wxgll(i)*wygll(j)*wzgll(k)
                    jacobianl = jacobianstore(i,j,k,ispec)
                    qp_dp = qp_dp + jacobianl * weight * wks_1n(i,j,k,ispec,ipar) **2
                 enddo
              enddo
           enddo
        enddo
     enddo

     qp_tmp = qp_dp * coeff_n1_dp * coeff_n1_dp
     qp=0.
     call sum_all_all_cr(qp_tmp, qp)

  end subroutine Parallel_ComputeL2normSquare

end module inversion_scheme
