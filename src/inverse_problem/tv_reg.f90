module tv_reg
  use specfem_par,only : CUSTOM_REAL
  real(KIND=CUSTOM_REAL) delta_tv
  real(KIND=CUSTOM_REAL), allocatable, dimension(:,:,:) :: d1m,d2m,d3m
  real(KIND=CUSTOM_REAL), allocatable, dimension(:,:,:) :: d11m,d22m,d33m
  real(KIND=CUSTOM_REAL), allocatable, dimension(:,:,:) :: d12m,d13m,d23m

contains

  subroutine allocate_wksp(nx,ny,nz)

    implicit none
    integer nx,ny,nz

    allocate(d1m(nx,ny,nz),d2m(nx,ny,nz),d3m(nx,ny,nz))
    allocate(d11m(nx,ny,nz),d22m(nx,ny,nz),d33m(nx,ny,nz))
    allocate(d12m(nx,ny,nz),d13m(nx,ny,nz),d23m(nx,ny,nz))


  end subroutine allocate_wksp




  subroutine sub_grad_tv(m,grad_tv,nx,ny,nz,hx,hy,hz)

    use deriv_DF_mod
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL), intent(in) :: m(nx,ny,nz)
    real(KIND=CUSTOM_REAL) grad_tv(nx,ny,nz)
    real(KIND=CUSTOM_REAL) hx,hy,hz

    call sub_D1(m,d1m,nx,ny,nz,hx)
    call sub_D2(m,d2m,nx,ny,nz,hy)
    call sub_D3(m,d3m,nx,ny,nz,hz)

    d11m(:,:,:) = sqrt(d1m(:,:,:)**2 + d2m(:,:,:)**2 + d3m(:,:,:)**2 ) + delta_tv

!!$    open(666,file='denm.bin',access='direct',recl=4.*nx*ny*nz)
!!$    write(666,rec=1) sqrt(d1m(:,:,:)**2 + d2m(:,:,:)**2 + d3m(:,:,:)**2 )
!!$    close(666)
!!$
!!$    open(666,file='model_used.bin',access='direct',recl=4.*nx*ny*nz)
!!$    write(666,rec=1) m
!!$    close(666)

    grad_tv(:,:,:) =  d1m(:,:,:) /  d11m(:,:,:)
    call sub_D1(grad_tv,d1m,nx,ny,nz,hx)

    grad_tv(:,:,:) =  d2m(:,:,:) /  d11m(:,:,:)
    call sub_D2(grad_tv,d2m,nx,ny,nz,hy)

    grad_tv(:,:,:) =  d3m(:,:,:) /   d11m(:,:,:)
    call sub_D3(grad_tv,d3m,nx,ny,nz,hz)


    grad_tv(:,:,:) = d1m(:,:,:) + d2m(:,:,:) + d3m(:,:,:)
!!$    open(666,file='num.bin',access='direct',recl=4.*nx*ny*nz)
!!$    write(666,rec=1) grad_tv
!!$    close(666)


  end subroutine sub_grad_tv

  subroutine sub_pena_tv(m,pena,nx,ny,nz,hx,hy,hz)

    use deriv_DF_mod
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL) m(nx,ny,nz),pena
    real(KIND=CUSTOM_REAL) hx,hy,hz

    call sub_D1(m,d1m,nx,ny,nz,hx)
    call sub_D2(m,d2m,nx,ny,nz,hy)
    call sub_D3(m,d3m,nx,ny,nz,hz)

    pena=sum(sqrt(d1m(:,:,:)**2 + d2m(:,:,:)**2 + d3m(:,:,:)**2 ))


  end subroutine sub_pena_tv


  subroutine def_grad_mod(gradm,m,nx,ny,nz,hx,hy,hz)
    use deriv_DF_mod
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL) m(nx,ny,nz),gradm(nx,ny,nz)
    real(KIND=CUSTOM_REAL) hx,hy,hz

    call sub_D1(m,d1m,nx,ny,nz,hx)
    call sub_D2(m,d2m,nx,ny,nz,hy)
    call sub_D3(m,d3m,nx,ny,nz,hz)

    gradm(:,:,:) = d1m(:,:,:)**2 + d2m(:,:,:)**2 + d3m(:,:,:)**2


  end subroutine def_grad_mod

  subroutine def_grad_mod_l1(gradm,m,nx,ny,nz,hx,hy,hz)
    use deriv_DF_mod
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL) m(nx,ny,nz),gradm(nx,ny,nz)
    real(KIND=CUSTOM_REAL) hx,hy,hz

    call sub_D1(m,d1m,nx,ny,nz,hx)
    call sub_D2(m,d2m,nx,ny,nz,hy)
    call sub_D3(m,d3m,nx,ny,nz,hz)

    gradm(:,:,:) = sqrt(d1m(:,:,:)**2 + d2m(:,:,:)**2 + d3m(:,:,:)**2)


  end subroutine def_grad_mod_l1

  subroutine def_laplac_mod(lap,m,nx,ny,nz,hx,hy,hz)
    use deriv_DF_mod
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL) m(nx,ny,nz),lap(nx,ny,nz)
    real(KIND=CUSTOM_REAL) hx,hy,hz
    call sub_D11(m,d11m,nx,ny,nz,hx)
    call sub_D22(m,d22m,nx,ny,nz,hy)
    call sub_D33(m,d33m,nx,ny,nz,hz)

    lap(:,:,:)=d11m(:,:,:)+d22m(:,:,:)+d33m(:,:,:)

  end subroutine def_laplac_mod


!! ------------------------------- a priori Tarantola-Valette --------------

  subroutine inv_exp_l1_cov(m,grad,nx,ny,nz,hx,hy,hz,sigma,lambda)
    !
    ! Cm-1 * (m - mprior)
    !
    ! below m means (m - m_prior)
    use deriv_DF_mod
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL), intent(in) :: m(nx,ny,nz),sigma,lambda
    real(KIND=CUSTOM_REAL) grad(nx,ny,nz)
    real(KIND=CUSTOM_REAL) hx,hy,hz
    real(KIND=CUSTOM_REAL) c0,c1,c2,c3

    !
    ! from Tarantola Inverse Problems P313
    !

    c0 = (hx*hy*hz) / (8._CUSTOM_REAL*3.141592653589793*sigma*sigma)
    c1 = 1._CUSTOM_REAL*c0 / (lambda**3)
    c2 = 2._CUSTOM_REAL*c0 / lambda
    c3 = lambda*c0


    call def_laplac_mod(d1m,m,nx,ny,nz,hx,hy,hz)
    call def_laplac_mod(d2m,d1m,nx,ny,nz,hx,hy,hz)

    grad(:,:,:) =  c1*m(:,:,:) - c2*d1m(:,:,:) + c3*d2m(:,:,:)

    open(666,file='grad_pena.bin',access='direct',recl=4*nx*ny*nz)
    write(666,rec=1) grad
    close(666)

  end subroutine inv_exp_l1_cov


  subroutine norm_inv_exp_l1_cov(m,grad,pena,nx,ny,nz)
    !
    ! (m - mprior)t *Cm-1 * (m-mprior)
    !
    implicit none
    integer nx,ny,nz
    real(KIND=CUSTOM_REAL), intent(in) :: m(nx,ny,nz),grad(nx,ny,nz)
    real(KIND=CUSTOM_REAL) :: pena

    pena=sum(m(:,:,:)*grad(:,:,:))


  end subroutine norm_inv_exp_l1_cov


end module tv_reg
