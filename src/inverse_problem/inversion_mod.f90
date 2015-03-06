module inversion_mod
  use specfem_par
  use specfem_par_elastic
  use project_tomo_grid_mod
  use tv_reg
  implicit none
  integer iter,iter_wolfe
  integer, parameter :: IIOUT_INV=670
  integer, parameter :: maxiter_wolfe=10
  real(KIND=CUSTOM_REAL) , allocatable :: model_inv_0(:,:,:),model_inv_1(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: grad_save1(:,:,:),grad_save2(:,:,:),grad_save3(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: grad_inv(:,:,:),grad_inv_0(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: grad_inv_1(:,:,:),grad_inv_tmp(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: pert_alpha(:,:,:),pert_beta(:,:,:),pert_rho(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: grad_alpha(:,:,:,:),grad_beta(:,:,:,:),grad_rho(:,:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: sum_grad_alpha(:,:,:),sum_grad_beta(:,:,:),sum_grad_rho(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: reg_wks(:,:,:),precond(:,:,:),wk(:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: pert_alpha_in_SEM_mesh(:,:,:,:),pert_beta_in_SEM_mesh(:,:,:,:),pert_rho_in_SEM_mesh(:,:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: alpha_in_SEM_mesh(:,:,:,:),beta_in_SEM_mesh(:,:,:,:), rho_in_SEM_mesh(:,:,:,:)
  real(KIND=CUSTOM_REAL) , allocatable :: alpha_in_SEM_mesh_saved(:,:,:,:),beta_in_SEM_mesh_saved(:,:,:,:), rho_in_SEM_mesh_saved(:,:,:,:)
  real(KIND=CUSTOM_REAL) step_length,q0,qt,qp0,qpt,vqpt,tg,td
  real(KIND=CUSTOM_REAL) fit0,fit1,dfit0,dfit1,deriv_adjust,penalty_value1,penalty_value2,penalty_value3
  logical FLAG_WOLFE,GRADIENT,PERTURB_IN_TOMO_GRID
  real(KIND=CUSTOM_REAL), parameter ::  mwl1=0.01
  real(KIND=CUSTOM_REAL), parameter :: mwl2=0.7
  !integer maxiter_wolfe
  contains

!------------------------------------------------------
  subroutine allocate_inv_wsk(nsrc)
   implicit none
   integer nsrc
   allocate(model_inv_0(nx,ny,nz),model_inv_1(nx,ny,nz))
   allocate(grad_inv_0(nx,ny,nz),grad_inv_1(nx,ny,nz))
   allocate(grad_save1(nx,ny,nz),grad_save2(nx,ny,nz),grad_save3(nx,ny,nz))
   allocate(grad_inv(nx,ny,nz),grad_inv_tmp(nx,ny,nz))
   allocate(pert_alpha(nx,ny,nz),pert_beta(nx,ny,nz),pert_rho(nx,ny,nz))
   allocate(grad_alpha(nx,ny,nz,nsrc),grad_beta(nx,ny,nz,nsrc),grad_rho(nx,ny,nz,nsrc))
   allocate(sum_grad_alpha(nx,ny,nz),sum_grad_beta(nx,ny,nz),sum_grad_rho(nx,ny,nz))
   allocate(reg_wks(nx,ny,nz),precond(nx,ny,nz))
   allocate(wk(nx,ny,nz))


   call allocate_wksp(nx,ny,nz)

   precond(:,:,:)=1._CUSTOM_REAL

   if (myrank==0) open(IIOUT_INV,file='log.adjust')
   !maxiter_wolfe=6
   GRADIENT=.false.
   PERTURB_IN_TOMO_GRID=.false.

   if (.not.PERTURB_IN_TOMO_GRID) then

     allocate(pert_alpha_in_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(pert_beta_in_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(pert_rho_in_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))

     allocate(alpha_in_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(beta_in_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(rho_in_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))

     allocate(alpha_in_SEM_mesh_saved(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(beta_in_SEM_mesh_saved(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(rho_in_SEM_mesh_saved(NGLLX,NGLLY,NGLLZ,NSPEC_AB))


   endif



  end subroutine allocate_inv_wsk
!------------------------------------------------------


  subroutine compute_descent_direction(nsrc)
    implicit none
    real(kind=CUSTOM_REAL) a_k(10000),p_k(10000)
    real(kind=CUSTOM_REAL) beta,pk,norme_yiter
    real(kind=CUSTOM_REAL) maxi_alpha,maxi_beta,maxi_grad
    integer iter_get,k,nsrc,ier
    character(len=5) name_data


    if (myrank==0)  then

      if (use_precond) precond(:,:,:)=sqrt(profondeur_tomo_grid(:,:,:)/1000.)
      write(*,*) 'des dir',iter
      a_k(:)=0._CUSTOM_REAL
      p_k(:)=0._CUSTOM_REAL
      name_data='alpha'
      iter_get=iter

      call get_gradient(iter_get,nsrc,grad_save1,name_data)
      call read_mod(model_inv_0,iter,name_data)
      !call add_penalty(grad_save1,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

      pert_alpha(:,:,:)=grad_save1(:,:,:)
      grad_inv(:,:,:)=grad_save1(:,:,:)

      !write(*,*) 'verif pert : ',pert_alpha(10,10,10), profondeur_tomo_grid(10,10,10)

      !pert_alpha = (profondeur_tomo_grid/1000)* pert_alpha

      if (iter == 0 .or. GRADIENT) then ! 1rst iteration steepest descent
        !! lecture du gradient courrant
        !! lecture du modele courant
        !call read_mod(model_inv_0,iter,name_data)
        !pert_alpha(:,:,:) =-1_CUSTOM_REAL*grad_inv(:,:,:)
        if (use_precond) pert_alpha = precond * pert_alpha

      else ! l-bfgs method

        do k=iter-1,0,-1

           call get_gradient(k,  nsrc,grad_inv_0,name_data)
           call get_gradient(k+1,nsrc,grad_inv_1,name_data)

           call read_mod(model_inv_0,k,name_data)
           call read_mod(model_inv_1,k+1,name_data)

           !call add_penalty(grad_inv_0,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
           !call add_penalty(grad_inv_1,model_inv_1,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)


           if (k == iter - 1) &
           norme_yiter = sum( (grad_inv_1 - grad_inv_0) **2   )

           p_k(k) = sum( (grad_inv_1 - grad_inv_0) * (model_inv_1 - model_inv_0) )
           p_k(k) = 1._CUSTOM_REAL / p_k(k)

           a_k(k) = p_k(k) * sum((model_inv_1 - model_inv_0)*grad_inv)

           grad_inv = grad_inv - a_k(k) * (grad_inv_1 - grad_inv_0)

        enddo

        k=iter-1
        pk = 1. / (p_k(k) * norme_yiter)
        pert_alpha=pk*grad_inv

        !pert_alpha = (profondeur_tomo_grid/1000)* pert_alpha ! preconditionnning
         if (use_precond) pert_alpha = precond * pert_alpha

        do k=0,iter-1

          call get_gradient(k,  nsrc,grad_inv_0,name_data)
          call get_gradient(k+1,nsrc,grad_inv_1,name_data)

          call read_mod(model_inv_0,k,name_data)
          call read_mod(model_inv_1,k+1,name_data)

          !call add_penalty(grad_inv_0,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
          !call add_penalty(grad_inv_1,model_inv_1,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)


          beta = p_k(k) * sum ((grad_inv_1 - grad_inv_0)*pert_alpha)

          pert_alpha = pert_alpha + (a_k(k) - beta ) * (model_inv_1 - model_inv_0)
        enddo


      endif

      call add_tapper(pert_alpha,nx,ny,nz)

      maxi_alpha=maxval(abs(pert_alpha))
      pert_alpha(:,:,:)=-1._CUSTOM_REAL* pert_alpha(:,:,:)
      !deriv_adjust = sum(grad_save(:,:,:)*pert_alpha(:,:,:))

      name_data='beta_'
      call get_gradient(iter_get,nsrc,grad_save2,name_data)
      call read_mod(model_inv_0,iter,name_data)
      !call add_penalty(grad_save2,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
      pert_beta(:,:,:)=grad_save2(:,:,:)
      grad_inv(:,:,:)=grad_save2(:,:,:)
      !if (use_precond) pert_beta = precond * pert_beta

      if (iter == 0 .or. GRADIENT) then ! 1rst iteration steepest descent
        !! lecture du gradient courrant
        !! lecture du modele courant
        !call read_mod(model_inv_0,iter,name_data)
        !pert_beta(:,:,:) =-1_CUSTOM_REAL*grad_inv(:,:,:)
        if (use_precond) pert_beta = precond * pert_beta

      else ! l-bfgs method

        do k=iter-1,0,-1

          call get_gradient(k,  nsrc,grad_inv_0,name_data)
          call get_gradient(k+1,nsrc,grad_inv_1,name_data)

          call read_mod(model_inv_0,k,name_data)
          call read_mod(model_inv_1,k+1,name_data)

          !call add_penalty(grad_inv_0,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
          !call add_penalty(grad_inv_1,model_inv_1,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

          if (k == iter - 1) &
          norme_yiter = sum( (grad_inv_1 - grad_inv_0) **2   )

          p_k(k) = sum( (grad_inv_1 - grad_inv_0) * (model_inv_1 - model_inv_0) )
          p_k(k) = 1._CUSTOM_REAL / p_k(k)

          a_k(k) = p_k(k) * sum((model_inv_1 - model_inv_0)*grad_inv)

          grad_inv = grad_inv - a_k(k) * (grad_inv_1 - grad_inv_0)

        enddo

        k=iter-1
        pk = 1. / (p_k(k) * norme_yiter)
        pert_beta=pk*grad_inv

        !pert_beta = (profondeur_tomo_grid/1000)* pert_beta
        if (use_precond) pert_beta = precond * pert_beta

        do k=0,iter-1

          call get_gradient(k,  nsrc,grad_inv_0,name_data)
          call get_gradient(k+1,nsrc,grad_inv_1,name_data)

          call read_mod(model_inv_0,k,name_data)
          call read_mod(model_inv_1,k+1,name_data)

          !call add_penalty(grad_inv_0,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
          !call add_penalty(grad_inv_1,model_inv_1,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)


          beta = p_k(k) * sum ((grad_inv_1 - grad_inv_0)*pert_beta)

          pert_beta = pert_beta + (a_k(k) - beta ) * (model_inv_1 - model_inv_0)

        enddo


      endif
      !! tapper in pert

      call add_tapper(pert_beta,nx,ny,nz)

      maxi_beta=maxval(abs(pert_beta))
      pert_beta(:,:,:)=-1._CUSTOM_REAL* pert_beta(:,:,:)
      !deriv_adjust = deriv_adjust + sum(grad_save(:,:,:)*pert_beta(:,:,:))
      if (iparam_inv==1) then ! vp,vs
         maxi_grad=max(maxi_alpha,maxi_beta)
      else if (iparam_inv==2) then
         maxi_grad=maxi_alpha
         pert_beta=0.
      endif

      maxi_grad=maxi_grad/0.025 ! 2.5% pert maxi

      pert_alpha=pert_alpha/maxi_grad
      pert_beta=pert_beta/maxi_grad

      deriv_adjust =sum(grad_save1(:,:,:)*pert_alpha(:,:,:)) + sum(grad_save2(:,:,:)*pert_beta(:,:,:))

      write(*,*)
      write(*,*) 'Pert alpha %  :', 100.*minval(pert_alpha),100.*maxval(pert_alpha)
      write(*,*) 'Pert beta  %  :', 100.*minval(pert_beta),100.*maxval(pert_beta)
      write(*,*)

      open(667,file='pert_alpha',access='direct',recl=CUSTOM_REAL*nx*ny*nz)
      write(667,rec=1) pert_alpha
      close(667)

      open(667,file='pert_beta',access='direct',recl=CUSTOM_REAL*nx*ny*nz)
      write(667,rec=1) pert_beta
      close(667)

      open(667,file='precond_bfgs',access='direct',recl=CUSTOM_REAL*nx*ny*nz)
      write(667,rec=1) precond
      close(667)

    endif

    call mpi_bcast(pert_alpha,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(pert_beta,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(deriv_adjust,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine compute_descent_direction


  subroutine compute_descent_direction_1(nsrc)
    implicit none

    real(kind=CUSTOM_REAL) beta,pk,norme_yiter
    real(kind=CUSTOM_REAL) maxi_alpha,maxi_beta,maxi_rho,maxi_grad
    integer k,nsrc,ier
    character(len=5) name_data


    if (myrank==0) then

       if (iter == 0 .or. GRADIENT) then
          name_data='alpha'
          call GRAD(name_data,pert_alpha,grad_save1,nsrc)

          name_data='beta_'
          call GRAD(name_data,pert_beta,grad_save2,nsrc)

          name_data='_rho_'
          call GRAD(name_data,pert_rho,grad_save3,nsrc)

       else
          name_data='alpha'
          call L_BFGS(name_data,pert_alpha,grad_save1,nsrc)

          name_data='beta_'
          call L_BFGS(name_data,pert_beta,grad_save2,nsrc)

          name_data='_rho_'
          call L_BFGS(name_data,pert_rho,grad_save3,nsrc)
       endif

       call add_tapper(pert_alpha,nx,ny,nz)
       call add_tapper(pert_beta,nx,ny,nz)
       call add_tapper(pert_rho,nx,ny,nz)

       maxi_alpha=maxval(abs(pert_alpha))
       maxi_beta=maxval(abs(pert_beta))
       maxi_rho=maxval(abs(pert_rho))

       if (iparam_inv==1) then
          maxi_grad=max(maxi_alpha,maxi_beta)
          !pert_alpha=pert_alpha/maxi_grad
          !pert_beta=pert_beta/maxi_grad
          pert_rho=0._CUSTOM_REAL
       else if (iparam_inv==2) then
          maxi_grad=maxi_alpha
          !pert_alpha=pert_alpha/maxi_grad
          pert_beta=0._CUSTOM_REAL
          pert_rho=0._CUSTOM_REAL

       else if (iparam_inv==3) then
          WRITE(*,*) MAXI_GRAD
          maxi_grad=max(maxi_alpha,maxi_beta)
          WRITE(*,*) MAXI_GRAD
          maxi_grad=max(maxi_grad,maxi_rho)
          WRITE(*,*) MAXI_GRAD

       else
          write(*,*) ' stop: familly parameter not implemented yet '
          stop
       endif




       maxi_grad=maxi_grad/0.025 ! 2.5% pert maxi

       pert_alpha=pert_alpha/maxi_grad
       pert_beta=pert_beta/maxi_grad
       pert_rho=pert_rho/maxi_grad

       deriv_adjust =sum(grad_save1(:,:,:)*pert_alpha(:,:,:)) &
                   + sum(grad_save2(:,:,:)*pert_beta(:,:,:))&
                   + sum(grad_save3(:,:,:)*pert_rho(:,:,:))


       write(*,*)
       write(*,*)  maxi_grad
       write(*,*) 'Pert alpha %  :', 100.*minval(pert_alpha),100.*maxval(pert_alpha)
       write(*,*) 'Pert beta  %  :', 100.*minval(pert_beta),100.*maxval(pert_beta)
       write(*,*) 'Pert rho   %  :', 100.*minval(pert_rho),100.*maxval(pert_rho)
       write(*,*)

       OPEN(666,FILE='PERTa',ACCESS='DIRECT',RECL=4*NX*NY*NZ)
       WRITE(666,REC=1) pert_alpha
       CLOSE(666)


       OPEN(666,FILE='PERTb',ACCESS='DIRECT',RECL=4*NX*NY*NZ)
       WRITE(666,REC=1) pert_beta
       CLOSE(666)


       OPEN(666,FILE='PERTr',ACCESS='DIRECT',RECL=4*NX*NY*NZ)
       WRITE(666,REC=1) pert_rho
       CLOSE(666)

    endif

    call mpi_bcast(pert_alpha,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(pert_beta,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(pert_rho,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(deriv_adjust,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)


  end subroutine compute_descent_direction_1


  subroutine GRAD(name_data,pert_mod,grad_save,nsrc)
    real(kind=CUSTOM_REAL) pert_mod(nx,ny,nz),grad_save(nx,ny,nz)
    integer iter_get,nsrc
    character(len=5) name_data

    iter_get=iter
    if (use_precond) precond(:,:,:)=sqrt(profondeur_tomo_grid(:,:,:)/1000.)
    call get_gradient(iter_get,nsrc,grad_save,name_data)
    pert_mod(:,:,:)=grad_save(:,:,:)
    pert_mod(:,:,:) = -1._CUSTOM_REAL* pert_mod(:,:,:)

  end subroutine GRAD



  subroutine L_BFGS(name_data,pert_mod,grad_save,nsrc)
    real(kind=CUSTOM_REAL) a_k(10000),p_k(10000)
    real(kind=CUSTOM_REAL) pert_mod(nx,ny,nz),grad_save(nx,ny,nz)
    real(kind=CUSTOM_REAL) beta,pk,norme_yiter
    integer iter_get,k,nsrc
    character(len=5) name_data

    iter_get=iter
    a_k(:)=0._CUSTOM_REAL
    p_k(:)=0._CUSTOM_REAL

    if (use_precond) precond(:,:,:)=sqrt(profondeur_tomo_grid(:,:,:)/1000.)
    call get_gradient(iter_get,nsrc,grad_save,name_data)
    call read_mod(model_inv_0,iter,name_data)
    pert_mod(:,:,:)=grad_save(:,:,:)
    grad_inv(:,:,:)=grad_save(:,:,:)


    do k=iter-1,0,-1

       call get_gradient(k,  nsrc,grad_inv_0,name_data)
       call get_gradient(k+1,nsrc,grad_inv_1,name_data)

       call read_mod(model_inv_0,k,name_data)
       call read_mod(model_inv_1,k+1,name_data)

       !call add_penalty(grad_inv_0,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
       !call add_penalty(grad_inv_1,model_inv_1,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)


       if (k == iter - 1) &
            norme_yiter = sum( (grad_inv_1 - grad_inv_0) **2   )

       p_k(k) = sum( (grad_inv_1 - grad_inv_0) * (model_inv_1 - model_inv_0) )
       p_k(k) = 1._CUSTOM_REAL / p_k(k)

       a_k(k) = p_k(k) * sum((model_inv_1 - model_inv_0)*grad_inv)

       grad_inv = grad_inv - a_k(k) * (grad_inv_1 - grad_inv_0)

    enddo

    k=iter-1
    pk = 1. / (p_k(k) * norme_yiter)
    pert_mod=pk*grad_inv

    !pert_alpha = (profondeur_tomo_grid/1000)* pert_alpha ! preconditionnning
    if (use_precond) pert_mod = precond * pert_mod


    do k=0,iter-1

       call get_gradient(k,  nsrc,grad_inv_0,name_data)
       call get_gradient(k+1,nsrc,grad_inv_1,name_data)

       call read_mod(model_inv_0,k,name_data)
       call read_mod(model_inv_1,k+1,name_data)

       !call add_penalty(grad_inv_0,model_inv_0,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
       !call add_penalty(grad_inv_1,model_inv_1,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)


       beta = p_k(k) * sum ((grad_inv_1 - grad_inv_0)*pert_mod)

       pert_mod = pert_mod + (a_k(k) - beta ) * (model_inv_1 - model_inv_0)

    enddo

    pert_mod(:,:,:) = -1._CUSTOM_REAL* pert_mod(:,:,:)

  end subroutine L_BFGS



  subroutine get_gradient(k,nsrc,grad,name_data)
    implicit none
    real(kind=CUSTOM_REAL) grad(nx,ny,nz)
    character(len=5) name_data
    integer isrc,nsrc,k

    grad(:,:,:)=0._CUSTOM_REAL
    !do isrc=1,Nsrc
      !write(*,*) 'read kl',isrc,k
      call read_kl(grad_inv_tmp,isrc,k,name_data)
      grad(:,:,:)=grad(:,:,:)+grad_inv_tmp(:,:,:)
    !enddo

  end subroutine get_gradient

!-------------------------------------------------------

  subroutine compute_new_model()

    implicit none
    character(len=5) name_data

    if (PERTURB_IN_TOMO_GRID) then
      if (myrank==0) then
        name_data='alpha'
        write(*,*)
        write(*,*) 'read mod ',trim(name_data),iter
        write(*,*)
        call read_mod(model_inv_0,iter,name_data)
        alpha_grid(:,:,:)=model_inv_0(:,:,:)*exp(step_length*pert_alpha(:,:,:))

        name_data='beta_'
        call read_mod(model_inv_0,iter,name_data)
        beta_grid(:,:,:)=model_inv_0(:,:,:)*exp(step_length*pert_beta(:,:,:))

        write(*,*)
        write(*,*) 'perturbed alp:',minval(alpha_grid),maxval(alpha_grid)
        write(*,*) 'perturbed bet:',minval(beta_grid),maxval(beta_grid)
        write(*,*)

      endif

      !write(*,*) 'new_model 1',myrank
      call mpi_bcast(alpha_grid,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      !write(*,*) 'new_model 2',myrank

      call mpi_bcast(beta_grid,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      !write(*,*) 'new_model 3',myrank

      !! project new model in SEM grid
      ! projection de rho*vs
      data_tomo(:,:,:) = rho_grid(:,:,:)*beta_grid(:,:,:)
      call project_sem_interp(data_tomo,rho_vs)
      !write(*,*) 'new_model 4',myrank

      ! projection de rho*vp
      data_tomo(:,:,:) = rho_grid(:,:,:)*alpha_grid(:,:,:)
      call project_sem_interp(data_tomo,rho_vp)


      ! projection de mu
      mu_grid(:,:,:) = rho_grid(:,:,:)*beta_grid(:,:,:)*beta_grid(:,:,:)
      call project_sem_interp(mu_grid,mustore)

      ! projection de kappa
      kappa_grid(:,:,:) = rho_grid(:,:,:)*( (alpha_grid(:,:,:)**2) -FOUR_THIRDS*(beta_grid(:,:,:)**2) )
      call project_sem_interp(kappa_grid,kappastore)
   else



      ! in a case when regularization is used, I have to compute this
      if (myrank==0) then

        name_data='alpha'
        call read_mod(model_inv_0,iter,name_data)
        alpha_grid(:,:,:)=model_inv_0(:,:,:)*exp(step_length*pert_alpha(:,:,:))

        name_data='beta_'
        call read_mod(model_inv_0,iter,name_data)
        beta_grid(:,:,:)=model_inv_0(:,:,:)*exp(step_length*pert_beta(:,:,:))

        name_data='_rho_'
        call read_mod(model_inv_0,iter,name_data)
        rho_grid(:,:,:)=model_inv_0(:,:,:)*exp(step_length*pert_rho(:,:,:))

        write(*,*)
        write(*,*) 'perturbed alp:',minval(alpha_grid),maxval(alpha_grid)
        write(*,*) 'perturbed bet:',minval(beta_grid),maxval(beta_grid)
        write(*,*) 'perturbed rho:',minval(rho_grid),maxval(rho_grid)
        write(*,*)

      endif

      !
      call mpi_bcast(alpha_grid,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      call mpi_bcast(beta_grid,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      call mpi_bcast(rho_grid,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

      ! project perturbation in SEM grid
      call mpi_bcast(pert_alpha,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pert_beta,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pert_rho,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)


      data_tomo(:,:,:) = step_length*pert_alpha(:,:,:)
      call project_sem_interp(data_tomo,pert_alpha_in_SEM_mesh)

      data_tomo(:,:,:) = step_length*pert_beta(:,:,:)
      call project_sem_interp(data_tomo,pert_beta_in_SEM_mesh)

      data_tomo(:,:,:) = step_length*pert_rho(:,:,:)
      call project_sem_interp(data_tomo,pert_rho_in_SEM_mesh)

      call new_model_in_SEM_mesh()


   endif

  end subroutine compute_new_model
!------------------------------------------------------
  subroutine save_current_model_SEM_in_memmory()

   implicit none
   integer i,j,k,iglob,ispec
   real(kind=CUSTOM_REAL) :: weight, jacobianl

   alpha_in_SEM_mesh_saved(:,:,:,:)  = ( kappastore(:,:,:,:) + (4./3.) * mustore(:,:,:,:) ) / rho_vp(:,:,:,:)
   beta_in_SEM_mesh_saved(:,:,:,:)   =  mustore(:,:,:,:) / rho_vs(:,:,:,:)
   rho_in_SEM_mesh_saved(:,:,:,:) = rho_vs(:,:,:,:) /  beta_in_SEM_mesh_saved(:,:,:,:)

  end subroutine save_current_model_SEM_in_memmory

  subroutine new_model_in_SEM_mesh()
    implicit none
    integer i,j,k,ispec,iglob
    real(kind=CUSTOM_REAL) weight,jacobianl

   !alpha_in_SEM_mesh(:,:,:,:)  =  alpha_in_SEM_mesh_saved(:,:,:,:)   !( kappastore(:,:,:,:) + (4./3.) * mustore(:,:,:,:) ) / rho_vp(:,:,:,:)
   !beta_in_SEM_mesh(:,:,:,:)   =  beta_in_SEM_mesh_saved(:,:,:,:)    !mustore(:,:,:,:) / rho_vs(:,:,:,:)
   !rho_in_SEM_mesh(:,:,:,:)    =  rho_in_SEM_mesh_saved(:,:,:,:)   !rho_vs(:,:,:,:) /  beta_in_SEM_mesh(:,:,:,:)

   alpha_in_SEM_mesh(:,:,:,:)  =  alpha_in_SEM_mesh_saved(:,:,:,:)  * exp(pert_alpha_in_SEM_mesh(:,:,:,:))
   beta_in_SEM_mesh(:,:,:,:)   =  beta_in_SEM_mesh_saved(:,:,:,:)   * exp(pert_beta_in_SEM_mesh(:,:,:,:))
   rho_in_SEM_mesh(:,:,:,:)    =  rho_in_SEM_mesh_saved(:,:,:,:)    * exp(pert_rho_in_SEM_mesh(:,:,:,:))

   ! new model
   rho_vs(:,:,:,:) = rho_in_SEM_mesh(:,:,:,:) * beta_in_SEM_mesh(:,:,:,:)
   rho_vp(:,:,:,:) = rho_in_SEM_mesh(:,:,:,:) * alpha_in_SEM_mesh(:,:,:,:)
   kappastore(:,:,:,:) = rho_in_SEM_mesh(:,:,:,:) * ( (alpha_in_SEM_mesh(:,:,:,:)**2) -FOUR_THIRDS*(beta_in_SEM_mesh(:,:,:,:)**2))
   mustore(:,:,:,:) = rho_in_SEM_mesh(:,:,:,:) * (beta_in_SEM_mesh(:,:,:,:)**2)

   ! new mass matrix
   rmass(:) = 0._CUSTOM_REAL
   do ispec=1,nspec_ab
      do k=1,NGLLZ
         do j=1,NGLLY
            do i=1,NGLLX
               iglob = ibool(i,j,k,ispec)
               weight = wxgll(i)*wygll(j)*wzgll(k)
               jacobianl = jacobian(i,j,k,ispec)
               rmass(iglob) = rmass(iglob) + &
                    jacobianl * weight * rho_in_SEM_mesh(i,j,k,ispec)
            enddo
         enddo
      enddo
   enddo


  end subroutine

!-------------------------------------------------------
  subroutine compute_deriv()
    implicit none
    deriv_adjust = sum(sum_grad_alpha(:,:,:)*pert_alpha(:,:,:))
    deriv_adjust = deriv_adjust + sum(sum_grad_beta(:,:,:)*pert_beta(:,:,:))
    deriv_adjust = deriv_adjust + sum(sum_grad_rho(:,:,:)*pert_rho(:,:,:))

  end subroutine compute_deriv
!-------------------------------------------------------
  subroutine wolfe_criteria()
    implicit none
    !logical FLAG_WOLFE
    !real(kind=CUSTOM_REAL) qp0,vqpt,qpt

    !write(*,*) 'in wolfe ',myrank,q0,qt,qp0,vqpt
    if (mwl2*qp0<=vqpt .and. qpt<=mwl1*qp0 ) then
      flag_wolfe=.true.
      !write(*,*) 'inwolfe 1',myrank
      return
    endif

    if (mwl1*qp0<qpt) then
      td=step_length
      if (myrank ==0) write(*,*) 'mwl1*qp0<qpt ',myrank, mwl1*qp0,qpt
    endif

    if (qpt<=mwl1*qp0 .and. vqpt<mwl2 *qp0) then
      tg=step_length
      if (myrank==0) write(*,*) 'inwolfe 3',myrank
    endif

    if (td==0.) then
      step_length = 2.*step_length
      if (myrank==0) write(*,*) 'inwolfe 4',myrank
    else
      step_length = 0.5*(td+tg)
      if (myrank==0) write(*,*) 'inwolfe 5',myrank
    endif



  end subroutine wolfe_criteria


  subroutine LA_pena(nrme,m,nx,ny,nz,l,lz)

    implicit none
    integer  nx,ny,nz,ix,iy,iz
    real(kind=CUSTOM_REAL) m(nx,ny,nz),nrme,l,lz
    nrme=0._CUSTOM_REAL
    do iz=2,nz-1
      do iy=2,ny-1
        do ix=2,nx-1
        ! lissage selon x et y
        nrme = nrme+0.5*( -4*l* m(ix,iy,iz) + l*m(ix-1,iy,iz) + l*m(ix+1,iy,iz) + l*m(ix,iy-1,iz) + l*m(ix,iy+1,iz))**2
        ! lissage selon z
        nrme = nrme+0.5*( -2*lz*m(ix,iy,iz) + lz*m(ix,iy,iz-1) + lz*m(ix,iy,iz+1))**2

        enddo
      enddo
    enddo

  end subroutine LA_pena



  subroutine LA_regu(q,m,mw,nx,ny,nz,l,lz)
    implicit none
    integer nx,ny,nz
    integer n,i,ix,iy,iz
    real(kind=CUSTOM_REAL) mw(nx,ny,nz),m(nx,ny,nz),q(nx,ny,nz),l,lz

    mw(:,:,:)=0.
    do iz=2,nz-1
       do iy=2,ny-1
          do ix=2,nx-1
            ! lissage selon x et y
            mw(ix,iy,iz) = -4*l* m(ix,iy,iz) + l*m(ix-1,iy,iz) + l*m(ix+1,iy,iz) + l*m(ix,iy-1,iz) + l*m(ix,iy+1,iz)
            ! lissage selon z
            mw(ix,iy,iz) = mw(ix,iy,iz) -2*lz*m(ix,iy,iz) + lz*m(ix,iy,iz-1) + lz*m(ix,iy,iz+1)

          enddo
       enddo
    enddo

    do iz=2,nz-1
       do iy=2,ny-1
          do ix=2,nx-1
            ! lissage selon x et y
            q(ix,iy,iz) = q(ix,iy,iz) -4*l* mw(ix,iy,iz) + l*mw(ix-1,iy,iz) + l*mw(ix+1,iy,iz) + l*mw(ix,iy-1,iz) + l*mw(ix,iy+1,iz)
            ! lissage selon z
            q(ix,iy,iz) = q(ix,iy,iz) -2*lz*mw(ix,iy,iz) + lz*mw(ix,iy,iz-1) + lz*mw(ix,iy,iz+1)
          enddo
       enddo
    enddo

  end subroutine LA_regu

  subroutine TV_regu(g,w,m,nx,ny,nz,l,wl,hx,hy,hz)
    implicit none
    integer nx,ny,nz
    real(kind=CUSTOM_REAL) m(nx,ny,nz),g(nx,ny,nz),w(nx,ny,nz),l,wl
    real(kind=CUSTOM_REAL) hx,hy,hz
    delta_tv=wl
    call sub_grad_tv(m,w,nx,ny,nz,hx,hy,hz)
    g(:,:,:)=g(:,:,:)-l*w(:,:,:)
  end subroutine TV_regu

  subroutine TV_pena(nrme,m,nx,ny,nz,l,hx,hy,hz)
    implicit none
    integer nx,ny,nz
    real(kind=CUSTOM_REAL) m(nx,ny,nz)
    real(kind=CUSTOM_REAL) hx,hy,hz,l,nrme
    call sub_pena_tv(m,nrme,nx,ny,nz,hx,hy,hz)
    nrme=l*nrme
  end subroutine TV_pena

  subroutine CO_regu(g,w,m,mp,nx,ny,nz,l,wl,hx,hy,hz)
    implicit none
    integer nx,ny,nz
    real(kind=CUSTOM_REAL) m(nx,ny,nz),g(nx,ny,nz)
    real(kind=CUSTOM_REAL) w(nx,ny,nz),mp(nx,ny,nz)
    real(kind=CUSTOM_REAL) l,wl
    real(kind=CUSTOM_REAL) hx,hy,hz

    wk(:,:,:)=m(:,:,:)-mp(:,:,:)
    call inv_exp_l1_cov(wk,w,nx,ny,nz,hx,hy,hz,l,wl)
    g(:,:,:)=g(:,:,:) - scale_pena * w(:,:,:)

  end subroutine CO_regu

  subroutine CO_pena(nrme,m,mp,g,nx,ny,nz)
    implicit none
    integer nx,ny,nz
    real(kind=CUSTOM_REAL) m(nx,ny,nz),mp(nx,ny,nz)
    real(kind=CUSTOM_REAL) g(nx,ny,nz)
    real(kind=CUSTOM_REAL) nrme

    wk(:,:,:)=m(:,:,:)-mp(:,:,:)
    call norm_inv_exp_l1_cov(wk,g,nrme,nx,ny,nz)
    nrme = 0.5 * nrme * scale_pena
  end subroutine CO_pena

  subroutine  add_penalty(gradient,model,model_prior,mw,nx,ny,nz,lambda_penalty,water_level,type_reg,hx,hy,hz)


    implicit none
    integer nx,ny,nz
    real(kind=CUSTOM_REAL) mw(nx,ny,nz),model(nx,ny,nz)
    real(kind=CUSTOM_REAL) gradient(nx,ny,nz),model_prior(nx,ny,nz)
    real(kind=CUSTOM_REAL) lambda_penalty,water_level
    real(kind=CUSTOM_REAL) hx,hy,hz
    character(len=2) type_reg

    select case (type_reg)

    case('LA')
       call LA_regu(gradient,model,mw,nx,ny,nz,lambda_penalty,water_level)
    case('TV')
       call TV_regu(gradient,mw,model,nx,ny,nz,lambda_penalty,water_level,hx,hy,hz)
    case('CO')
       call CO_regu(gradient,mw,model,model_prior,nx,ny,nz,lambda_penalty,water_level,hx,hy,hz)
    case('NO')
      lambda_penalty=0.
    case default
      lambda_penalty=0.
    end select

  end subroutine add_penalty


  subroutine compute_penalty_reg(model,model_prior,grad,nx,ny,nz,norme_penalisation,lambda_penalty,water_level,type_reg,hx,hy,hz)

  implicit none
  integer nx,ny,nz
  real(kind=CUSTOM_REAL) model(nx,ny,nz),model_prior(nx,ny,nz),grad(nx,ny,nz)
  real(kind=CUSTOM_REAL) lambda_penalty,water_level,norme_penalisation
  real(kind=CUSTOM_REAL) hx,hy,hz
  character(len=2) type_reg

  select case (type_reg)

  case('LA')
     call LA_pena(norme_penalisation,model,nx,ny,nz,lambda_penalty,water_level)
  case('TV')
     call TV_pena(norme_penalisation,model,nx,ny,nz,lambda_penalty,hx,hy,hz)
  case ('CO')
     call CO_pena(norme_penalisation,model,model_prior,grad,nx,ny,nz)
  case('NO')
     lambda_penalty=0.
     norme_penalisation=0.
  case default
     lambda_penalty=0.
     norme_penalisation=0.

  end select


 end subroutine compute_penalty_reg


!---

 subroutine addition_penalty(reg_par_1,reg_par_2,reg_type)

   implicit none
   character(len=2) reg_type,reg_type_custom
   real(kind=CUSTOM_REAL) reg_par_1,reg_par_2,hx,hy,hz

   select case (iparam_inv )

      case(1)  ! vp,vs

         call add_penalty(sum_grad_alpha,alpha_grid,alpha_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
         call compute_penalty_reg(alpha_grid,alpha_grid_prior,reg_wks,nx,ny,nz,penalty_value1,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

         call add_penalty(sum_grad_beta,beta_grid,beta_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
         call compute_penalty_reg(beta_grid,beta_grid_prior,reg_wks,nx,ny,nz,penalty_value2,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

         reg_type_custom='NO'
         call add_penalty(sum_grad_rho,rho_grid,rho_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type_custom,hx,hy,hz)
         call compute_penalty_reg(rho_grid,rho_grid_prior,reg_wks,nx,ny,nz,penalty_value3,reg_par_1,reg_par_2,reg_type_custom,hx,hy,hz)

      case(2) ! vp

         call add_penalty(sum_grad_alpha,alpha_grid,alpha_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
         call compute_penalty_reg(alpha_grid,alpha_grid_prior,reg_wks,nx,ny,nz,penalty_value1,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

         reg_type_custom='NO'
         call add_penalty(sum_grad_beta,beta_grid,beta_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type_custom,hx,hy,hz)
         call compute_penalty_reg(beta_grid,beta_grid_prior,reg_wks,nx,ny,nz,penalty_value2,reg_par_1,reg_par_2,reg_type_custom,hx,hy,hz)


         call add_penalty(sum_grad_rho,rho_grid,rho_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type_custom,hx,hy,hz)
         call compute_penalty_reg(rho_grid,rho_grid_prior,reg_wks,nx,ny,nz,penalty_value3,reg_par_1,reg_par_2,reg_type_custom,hx,hy,hz)

      case(3) ! vp, vs, rho

         call add_penalty(sum_grad_alpha,alpha_grid,alpha_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
         call compute_penalty_reg(alpha_grid,alpha_grid_prior,reg_wks,nx,ny,nz,penalty_value1,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

         call add_penalty(sum_grad_beta,beta_grid,beta_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
         call compute_penalty_reg(beta_grid,beta_grid_prior,reg_wks,nx,ny,nz,penalty_value2,reg_par_1,reg_par_2,reg_type,hx,hy,hz)


         call add_penalty(sum_grad_rho,rho_grid,rho_grid_prior,reg_wks,nx,ny,nz,reg_par_1,reg_par_2,reg_type,hx,hy,hz)
         call compute_penalty_reg(rho_grid,rho_grid_prior,reg_wks,nx,ny,nz,penalty_value3,reg_par_1,reg_par_2,reg_type,hx,hy,hz)

      end select


 end subroutine addition_penalty


!---

subroutine add_tapper(gd,nx,ny,nz)
   implicit none
   integer i0,i1,i2,i3
   integer j0,j1,j2,j3
   integer k0,k1,k2,k3
   integer nx,ny,nz
   integer i,j,k
   real(kind=CUSTOM_REAL) tapperx,tappery,tapperz,gd(nx,ny,nz)
   real(kind=CUSTOM_REAL) wx0,px0,wx1,px1
   real(kind=CUSTOM_REAL) wy0,py0,wy1,py1
   real(kind=CUSTOM_REAL) wz0,pz0,wz1,pz1
   integer irec
   !real(kind=CUSTOM_REAL) pi

   ! PI = 3.141592653589793
   i0=2
   i1=4
   i2=nx-i1+1
   i3=nx-i0+1

   wx0=pi/(i0-i1)
   px0=-i1*wx0

   wx1=pi/(i3-i2)
   px1=-i2*wx1


   j0=2
   j1=4
   j2=ny-j1+1
   j3=ny-j0+1

   wy0=pi/(j0-j1)
   py0=-j1*wx0

   wy1=pi/(j3-j2)
   py1=-j2*wx1


   k0=3
   k1=5
   k2=nz-k1+1
   k3=nz-k0+1


   wz0=pi/(k0-k1)
   pz0=-k1*wz0

   wz1=pi/(k3-k2)
   pz1=-k2*wz1
  if (myrank==0) open(667,file='tapper',access='direct',recl=4)
  irec=0
  do k=1,nz

    tapperz = 1._CUSTOM_REAL
    if (k<=k0 .or. k>=k3) tapperz = 0._CUSTOM_REAL
    if (k>k0 .and. k<=k1) tapperz = (0.5+0.5*cos(wz0*k+pz0))
    if (k>k2 .and. k<=k3) tapperz = (0.5+0.5*cos(wz1*k+pz1))
    do j=1,ny
     tappery = 1._CUSTOM_REAL
     if (j<=j0 .or. j>=j3) tappery = 0._CUSTOM_REAL
     if (j>j0 .and. j<=j1) tappery = (0.5+0.5*cos(wy0*j+py0))
     if (j>j2 .and. j<=j3) tappery = (0.5+0.5*cos(wy1*j+py1))
     do i=1,nx
       irec=irec+1
       tapperx =  1._CUSTOM_REAL
       if (i<=i0 .or. i>=i3) tapperx = 0._CUSTOM_REAL
       if (i>i0 .and. i<=i1) tapperx = (0.5+0.5*cos(wx0*i+px0))
       if (i>i2 .and. i<=i3) tapperx = (0.5+0.5*cos(wx1*i+px1))
       gd(i,j,k) = tapperx * tappery * tapperz * gd(i,j,k)
       !write(*,*)  tapperx,tappery,tapperz
       if (myrank==0) write(667,rec=irec) tapperx * tappery * tapperz
     enddo
   enddo
  enddo
  if (myrank==0) close(667)


end subroutine add_tapper






end module inversion_mod



