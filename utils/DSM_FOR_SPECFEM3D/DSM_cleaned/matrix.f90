
!! DK DK these subroutines to create the submatrices are called only once and for all
!! DK DK before the main loop on all discrete frequencies, thus it is *not* critical to optimize their performance level

  subroutine comp_submatrix_mod( ngrid_r,grid_r, &
            grid_rho,grid_kappa,grid_mu, &
            grid_Ak,grid_Am,grid_L,grid_N,grid_Fk,grid_Fm, &
            submatrix_I0,submatrix_I1k,submatrix_I3k, &
            submatrix_I3m,submatrix_I4,submatrix_I5k, &
            submatrix_I5m,submatrix_I6,submatrix_I7, &
            submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,xi_Gauss,weight_Gauss,n_structure_zone,boundary_r,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone)

! compute modified submatrices for elastic part of the medium

  implicit none

! variables for input/output
  integer ngrid_r,n_structure_zone
  real(kind=8) grid_r(*)
  real(kind=8) grid_rho(2,*),grid_kappa(2,*),grid_mu(2,*)
  real(kind=8) grid_Ak(2,*),grid_Am(2,*),grid_L(2,*)
  real(kind=8) grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
  real(kind=8) submatrix_I0(4,*),submatrix_I1k(4,*)
  real(kind=8) submatrix_I3k(4,*),submatrix_I3m(4,*),submatrix_I4(4,*)
  real(kind=8) submatrix_I5k(4,*),submatrix_I5m(4,*)
  real(kind=8) submatrix_I6(4,*),submatrix_I7(4,*)
  real(kind=8) submatrix_I3k_mod(6,*),submatrix_I3m_mod(6,*)
  real(kind=8) submatrix_I4_mod(6,*)
  integer boundary_r(n_structure_zone+1),flag,flag1,jjj,izone1
  real(kind=8) rho_structure_zone(4,n_structure_zone)
  real(kind=8) vpv_structure_zone(4,n_structure_zone)
  real(kind=8) vph_structure_zone(4,n_structure_zone)
  real(kind=8) vsv_structure_zone(4,n_structure_zone)
  real(kind=8) vsh_structure_zone(4,n_structure_zone)
  real(kind=8) eta_structure_zone(4,n_structure_zone)

! variables for numerical integrations using Simpson's or Gauss' quadrature rule for mass lumping
  include "integration_points_for_submatrices.h"

  integer ir,iy,istart,iend
  real(kind=8) rho0,Fk0,Fm0,Ak0,Am0,L0,N0,lambda0
  real(kind=8) y,r,dr,rhor(0:2*ns),laminvr(0:2*ns),rhoinv(0:2*ns)
  real(kind=8) Fkr(0:2*ns),Fmr(0:2*ns),Lr(0:2*ns)
  real(kind=8) Ak(0:2*ns),Am(0:2*ns),L(0:2*ns),N(0:2*ns)
  real(kind=8) rhor_ave,laminvr_ave,rhoinv_ave
  real(kind=8) Fkr_ave,Fmr_ave,Lr_ave,Ak_ave,Am_ave,L_ave,N_ave

!! DK DK added this for Gauss integration
  real(kind=8), dimension(ns) :: xi_Gauss,weight_Gauss
  real(kind=8) :: A,coef(4)

  real(kind=8), external :: grid_interpolate

! **********************************************************************
!  compute submatrix elements for each cell
! **********************************************************************

     flag = 1

  do ir=1,ngrid_r-1

  do jjj = 1,n_structure_zone
      if ( ir >= boundary_r(jjj) .and. ir + 1 <= boundary_r(jjj+1) )  then
       izone1 = jjj
       exit
      endif
  enddo

  dr = grid_r(ir+1) - grid_r(ir)

  if ( grid_mu(1,ir)*grid_mu(2,ir) /= 0.d0 ) then

! evaluate the density or elastic constants at each sub-grid
  if(USE_GAUSS_INTEGRATION) then
    istart = 1
    iend = ns
    A = (grid_r(ir+1) - grid_r(ir)) / 2.d0
  else
    istart = 0
    iend = 2*ns
!! DK DK this is unused in the case of Simpson
    A = 0.d0
  endif

  do iy = istart, iend
    if(USE_GAUSS_INTEGRATION) then
!! DK DK map Gauss points, which are in ]-1,+1[, to ]0,1[
      y = (xi_Gauss(iy) + 1.d0) / 2.d0
    else
      y = dble(iy) / dble(2*ns)
    endif
!! DK DK this is the interpolated radius
    r = ( 1.d0 - y ) * grid_r(ir) + y * grid_r(ir+1)
    call grid_interpolate_1(r,rho0,Fk0,Fm0,Ak0,Am0,L0,N0,rho_structure_zone(1:4,izone1),vpv_structure_zone(1:4,izone1),vph_structure_zone(1:4,izone1),vsv_structure_zone(1:4,izone1),vsh_structure_zone(1:4,izone1),eta_structure_zone(1:4,izone1) )
    rhor(iy) = rho0 * r * r
    Fkr(iy) = Fk0 * r
    Fmr(iy) = Fm0 * r
    Lr(iy) = L0 * r
    Ak(iy) = Ak0
    Am(iy) = Am0
    L(iy) = L0
    N(iy) = N0


!    rhor(iy) = grid_interpolate(y,grid_rho(1,ir), grid_rho(2,ir)) * r * r
!    Fkr(iy) = grid_interpolate(y,grid_Fk(1,ir), grid_Fk(2,ir)) * r
!    Fmr(iy) = grid_interpolate(y,grid_Fm(1,ir), grid_Fm(2,ir)) * r
!    Lr(iy) = grid_interpolate(y,grid_L(1,ir), grid_L(2,ir)) * r
!    Ak(iy) = grid_interpolate(y,grid_Ak(1,ir), grid_Ak(2,ir))
!    Am(iy) = grid_interpolate(y,grid_Am(1,ir), grid_Am(2,ir))
!    L(iy) = grid_interpolate(y,grid_L(1,ir), grid_L(2,ir))
!    N(iy) = grid_interpolate(y,grid_N(1,ir), grid_N(2,ir))
  enddo

! compute lumped parameters
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), rhor,rhor_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), Fkr,Fkr_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), Fmr,Fmr_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), Lr,Lr_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), Ak,Ak_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), Am,Am_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), L,L_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), N,N_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )

  rhor_ave = rhor_ave / dr
  Fkr_ave = Fkr_ave / dr
  Fmr_ave = Fmr_ave / dr
  Lr_ave = Lr_ave / dr
  Ak_ave = Ak_ave / dr
  Am_ave = Am_ave / dr
  L_ave = L_ave / dr
  N_ave = N_ave / dr

! -----------------------------------------------------------
! -------- compute the optimally accurate           ---------
! -------- submatrices                              ---------
! -----------------------------------------------------------

!! DK DK the explanation of these 1 / 12. correction coefficients
!! DK DK is in equation (1.7) and in the comments below equation (1.3) of
!! DK DK   Geller, R. J., Takeuchi, N., 1995,
!! DK DK   A new method for computing highly accurate DSM synthetic seismograms,
!! DK DK   Geophys. J. Int., vol. 123, p. 449-470.

! modified submatrix_I0
  submatrix_I0(1,ir) = submatrix_I0(1,ir) + dr / 12.d0 * rhor_ave
  submatrix_I0(2,ir) = submatrix_I0(2,ir) - dr / 12.d0 * rhor_ave
  submatrix_I0(3,ir) = submatrix_I0(3,ir) - dr / 12.d0 * rhor_ave
  submatrix_I0(4,ir) = submatrix_I0(4,ir) + dr / 12.d0 * rhor_ave

! modified submatrix_I3k_mod
  submatrix_I3k_mod(1,ir) = submatrix_I3k(1,ir) - 1.d0 / 12.d0 * Fkr_ave
  submatrix_I3k_mod(2,ir) = submatrix_I3k(2,ir) + 2.d0 / 12.d0 * Fkr_ave
  submatrix_I3k_mod(3,ir) = submatrix_I3k(3,ir) + 1.d0 / 12.d0 * Fkr_ave
  submatrix_I3k_mod(4,ir) = submatrix_I3k(4,ir) - 2.d0 / 12.d0 * Fkr_ave
  submatrix_I3k_mod(5,ir) = - 1.d0 / 12.d0 * Fkr_ave
  submatrix_I3k_mod(6,ir) =   1.d0 / 12.d0 * Fkr_ave

! modified submatrix_I3m_mod
  submatrix_I3m_mod(1,ir) = submatrix_I3m(1,ir) - 1.d0 / 12.d0 * Fmr_ave
  submatrix_I3m_mod(2,ir) = submatrix_I3m(2,ir) + 2.d0 / 12.d0 * Fmr_ave
  submatrix_I3m_mod(3,ir) = submatrix_I3m(3,ir) + 1.d0 / 12.d0 * Fmr_ave
  submatrix_I3m_mod(4,ir) = submatrix_I3m(4,ir) - 2.d0 / 12.d0 * Fmr_ave
  submatrix_I3m_mod(5,ir) = - 1.d0 / 12.d0 * Fmr_ave
  submatrix_I3m_mod(6,ir) =   1.d0 / 12.d0 * Fmr_ave

! modified submatrix_I4_mod
  submatrix_I4_mod(1,ir) = submatrix_I4(1,ir) + 1.d0 / 12.d0 * Lr_ave
  submatrix_I4_mod(2,ir) = submatrix_I4(2,ir) - 2.d0 / 12.d0 * Lr_ave
  submatrix_I4_mod(3,ir) = submatrix_I4(3,ir) - 1.d0 / 12.d0 * Lr_ave
  submatrix_I4_mod(4,ir) = submatrix_I4(4,ir) + 2.d0 / 12.d0 * Lr_ave
  submatrix_I4_mod(5,ir) =   1.d0 / 12.d0 * Lr_ave
  submatrix_I4_mod(6,ir) = - 1.d0 / 12.d0 * Lr_ave

! modified submatrix_I5k
  submatrix_I5k(1,ir) = submatrix_I5k(1,ir) + dr / 12.d0 * Ak_ave
  submatrix_I5k(2,ir) = submatrix_I5k(2,ir) - dr / 12.d0 * Ak_ave
  submatrix_I5k(3,ir) = submatrix_I5k(3,ir) - dr / 12.d0 * Ak_ave
  submatrix_I5k(4,ir) = submatrix_I5k(4,ir) + dr / 12.d0 * Ak_ave

! modified submatrix_I5m
  submatrix_I5m(1,ir) = submatrix_I5m(1,ir) + dr / 12.d0 * Am_ave
  submatrix_I5m(2,ir) = submatrix_I5m(2,ir) - dr / 12.d0 * Am_ave
  submatrix_I5m(3,ir) = submatrix_I5m(3,ir) - dr / 12.d0 * Am_ave
  submatrix_I5m(4,ir) = submatrix_I5m(4,ir) + dr / 12.d0 * Am_ave

! modified submatrix_I6
  submatrix_I6(1,ir) = submatrix_I6(1,ir) + dr / 12.d0 * L_ave
  submatrix_I6(2,ir) = submatrix_I6(2,ir) - dr / 12.d0 * L_ave
  submatrix_I6(3,ir) = submatrix_I6(3,ir) - dr / 12.d0 * L_ave
  submatrix_I6(4,ir) = submatrix_I6(4,ir) + dr / 12.d0 * L_ave

! modified submatrix_I7
  submatrix_I7(1,ir) = submatrix_I7(1,ir) + dr / 12.d0 * N_ave
  submatrix_I7(2,ir) = submatrix_I7(2,ir) - dr / 12.d0 * N_ave
  submatrix_I7(3,ir) = submatrix_I7(3,ir) - dr / 12.d0 * N_ave
  submatrix_I7(4,ir) = submatrix_I7(4,ir) + dr / 12.d0 * N_ave
  else

  if(USE_GAUSS_INTEGRATION) then
    istart = 1
    iend = ns
    A = (grid_r(ir+1) - grid_r(ir)) / 2.d0
  else
    istart = 0
    iend = 2*ns
!! DK DK this is unused in the case of Simpson
    A = 0.d0
  endif

  do iy = istart, iend
    if(USE_GAUSS_INTEGRATION) then
!! DK DK map Gauss points, which are in ]-1,+1[, to ]0,1[
      y = (xi_Gauss(iy) + 1.d0) / 2.d0
    else
      y = dble(iy) / dble(2*ns)
    endif
!! DK DK this is the interpolated radius
    r = ( 1.d0 - y ) * grid_r(ir) + y * grid_r(ir+1)
    call grid_interpolate_3(r,rho0,lambda0,rho_structure_zone(1:4,izone1),vpv_structure_zone(1:4,izone1),vph_structure_zone(1:4,izone1),vsv_structure_zone(1:4,izone1),vsh_structure_zone(1:4,izone1),eta_structure_zone(1:4,izone1) )


    laminvr(iy) = r * r / lambda0
    rhoinv(iy) = 1.d0 / rho0

!    laminvr(iy) = r * r / grid_interpolate(y,grid_kappa(1,ir), grid_kappa(2,ir))
!    rhoinv(iy) = 1.d0 / grid_interpolate(y,grid_rho(1,ir), grid_rho(2,ir))
  enddo

! compute lumped parameters
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), laminvr,laminvr_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), rhoinv,rhoinv_ave,weight_Gauss,A,USE_GAUSS_INTEGRATION )
  laminvr_ave = laminvr_ave / dr
  rhoinv_ave = rhoinv_ave / dr

! -----------------------------------------------------------
! -------- compute the optimally accurate           ---------
! -------- submatrices                              ---------
! -----------------------------------------------------------

! modified submatrix_I0
  submatrix_I0(1,ir) = submatrix_I0(1,ir) + dr / 12.d0 * laminvr_ave
  submatrix_I0(2,ir) = submatrix_I0(2,ir) - dr / 12.d0 * laminvr_ave
  submatrix_I0(3,ir) = submatrix_I0(3,ir) - dr / 12.d0 * laminvr_ave
  submatrix_I0(4,ir) = submatrix_I0(4,ir) + dr / 12.d0 * laminvr_ave

! modified submatrix_I1k
  submatrix_I1k(1,ir) = submatrix_I1k(1,ir) + dr / 12.d0 * rhoinv_ave
  submatrix_I1k(2,ir) = submatrix_I1k(2,ir) - dr / 12.d0 * rhoinv_ave
  submatrix_I1k(3,ir) = submatrix_I1k(3,ir) - dr / 12.d0 * rhoinv_ave
  submatrix_I1k(4,ir) = submatrix_I1k(4,ir) + dr / 12.d0 * rhoinv_ave
  endif

  enddo ! of do ir=1,ngrid_r-1


  end subroutine comp_submatrix_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_submatrix( ngrid_r,grid_r, &
              grid_rho,grid_kappa,grid_mu,grid_Ak,grid_Am, &
              grid_Ck,grid_Cm,grid_L,grid_N,grid_Fk,grid_Fm, &
              submatrix_I0,submatrix_I1k,submatrix_I1m,submatrix_I2, &
              submatrix_I3k,submatrix_I3m,submatrix_I4,submatrix_I5k, &
              submatrix_I5m,submatrix_I6,submatrix_I7,xi_Gauss,weight_Gauss,n_structure_zone,boundary_r,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone )

! compute submatrices for elastic part of the medium

  implicit none

! variables for input/output
  integer ngrid_r,n_structure_zone
  real(kind=8) grid_r(*)
  real(kind=8) grid_rho(2,*),grid_kappa(2,*),grid_mu(2,*)
  real(kind=8) grid_Ak(2,*),grid_Am(2,*),grid_Ck(2,*),grid_Cm(2,*)
  real(kind=8) grid_L(2,*),grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
  real(kind=8) submatrix_I0(4,*),submatrix_I1k(4,*),submatrix_I1m(4,*)
  real(kind=8) submatrix_I2(4,*),submatrix_I3k(4,*),submatrix_I3m(4,*)
  real(kind=8) submatrix_I4(4,*),submatrix_I5k(4,*),submatrix_I5m(4,*)
  real(kind=8) submatrix_I6(4,*),submatrix_I7(4,*)

! variables for numerical integrations using Simpson's or Gauss' quadrature rule for mass lumping
  include "integration_points_for_submatrices.h"

  integer ir,iy,ib,istart,iend
  real(kind=8) y,yw(2),dyw(2),r
  real(kind=8) rho0,lambda0,Ck0,Cm0,L0,Fk0,Fm0,Ak0,Am0,N0
  real(kind=8) val_integrand0(0:2*ns,4)
  real(kind=8) val_integrand1k(0:2*ns,4)
  real(kind=8) val_integrand1m(0:2*ns,4)
  real(kind=8) val_integrand2(0:2*ns,4)
  real(kind=8) val_integrand3k(0:2*ns,4)
  real(kind=8) val_integrand3m(0:2*ns,4)
  real(kind=8) val_integrand4(0:2*ns,4)
  real(kind=8) val_integrand5k(0:2*ns,4)
  real(kind=8) val_integrand5m(0:2*ns,4)
  real(kind=8) val_integrand6(0:2*ns,4)
  real(kind=8) val_integrand7(0:2*ns,4)
  integer boundary_r(n_structure_zone+1),flag,flag1,jjj,izone1
  real(kind=8) rho_structure_zone(4,n_structure_zone)
  real(kind=8) vpv_structure_zone(4,n_structure_zone)
  real(kind=8) vph_structure_zone(4,n_structure_zone)
  real(kind=8) vsv_structure_zone(4,n_structure_zone)
  real(kind=8) vsh_structure_zone(4,n_structure_zone)
  real(kind=8) eta_structure_zone(4,n_structure_zone)

 
!! DK DK added this for Gauss integration
  real(kind=8), dimension(ns) :: xi_Gauss,weight_Gauss
  real(kind=8) :: A

! function for numerical integrations
  real(kind=8), external :: grid_interpolate

! **********************************************************************
!  compute stiffness submatrix elements for each cell
! **********************************************************************


  do ir = 1,ngrid_r-1

  do jjj = 1,n_structure_zone
      if ( ir >= boundary_r(jjj) .and. ir + 1 <= boundary_r(jjj+1) )  then
        izone1 = jjj
        exit
      endif
  enddo

  if ( grid_mu(1,ir)*grid_mu(2,ir) /= 0.d0 ) then

! -----------------------------------------------------------
! --------  evaluate the stiffness and              ---------
! --------  the integrated function for stiffness   ---------
! --------  matrix at the sub-grids                 ---------
! -----------------------------------------------------------

  if(USE_GAUSS_INTEGRATION) then
    istart = 1
    iend = ns
    A = (grid_r(ir+1) - grid_r(ir)) / 2.d0
  else
    istart = 0
    iend = 2*ns
!! DK DK this is unused in the case of Simpson
    A = 0.d0
  endif

  do iy = istart, iend
    if(USE_GAUSS_INTEGRATION) then
!! DK DK map Gauss points, which are in ]-1,+1[, to ]0,1[
      y = (xi_Gauss(iy) + 1.d0) / 2.d0
    else
      y = dble(iy) / dble(2*ns)
    endif
!! DK DK this is the interpolated radius
    r = ( 1.d0 - y ) * grid_r(ir) + y * grid_r(ir+1)
    call grid_interpolate_2(r,rho0,Ck0,Cm0,L0,Fk0,Fm0,Ak0,Am0,N0,rho_structure_zone(1:4,izone1),vpv_structure_zone(1:4,izone1),vph_structure_zone(1:4,izone1),vsv_structure_zone(1,izone1),vsh_structure_zone(1:4,izone1),eta_structure_zone(1:4,izone1) )

!    rho0 = grid_interpolate(y,grid_rho(1,ir),grid_rho(2,ir) )
!    Ck0 = grid_interpolate(y,grid_Ck(1,ir),grid_Ck(2,ir) )
!    Cm0 = grid_interpolate(y,grid_Cm(1,ir),grid_Cm(2,ir) )
!    L0 = grid_interpolate(y,grid_L(1,ir),grid_L(2,ir) )
!    Fk0 = grid_interpolate(y,grid_Fk(1,ir),grid_Fk(2,ir) )
!    Fm0 = grid_interpolate(y,grid_Fm(1,ir),grid_Fm(2,ir) )
!    Ak0 = grid_interpolate(y,grid_Ak(1,ir),grid_Ak(2,ir) )
!    Am0 = grid_interpolate(y,grid_Am(1,ir),grid_Am(2,ir) )
!    N0 = grid_interpolate(y,grid_N(1,ir),grid_N(2,ir) )
    yw(1) = 1.d0 - y
    yw(2) = y
    dyw(1) = - 1.d0 / ( grid_r(ir+1) - grid_r(ir) )
    dyw(2) =   1.d0 / ( grid_r(ir+1) - grid_r(ir) )
    do ib=1,4
      val_integrand0(iy,ib) = sub_integrand0( rho0, yw((ib+1)/2), yw(mod(ib+1,2)+1), r )
      val_integrand1k(iy,ib) = sub_integrand1( Ck0, dyw((ib+1)/2), dyw(mod(ib+1,2)+1), r )
      val_integrand1m(iy,ib) = sub_integrand1( Cm0, dyw((ib+1)/2), dyw(mod(ib+1,2)+1), r )
      val_integrand2(iy,ib) = sub_integrand2( L0, dyw((ib+1)/2), dyw(mod(ib+1,2)+1), r )
      val_integrand3k(iy,ib) = sub_integrand3( Fk0, yw((ib+1)/2), dyw(mod(ib+1,2)+1), r )
      val_integrand3m(iy,ib) = sub_integrand3( Fm0, yw((ib+1)/2), dyw(mod(ib+1,2)+1), r )
      val_integrand4(iy,ib) = sub_integrand4( L0, dyw((ib+1)/2), yw(mod(ib+1,2)+1), r )
      val_integrand5k(iy,ib) = sub_integrand5( Ak0, yw((ib+1)/2), yw(mod(ib+1,2)+1) )
      val_integrand5m(iy,ib) = sub_integrand5( Am0, yw((ib+1)/2), yw(mod(ib+1,2)+1) )
      val_integrand6(iy,ib) = sub_integrand6( L0, yw((ib+1)/2), yw(mod(ib+1,2)+1) )
      val_integrand7(iy,ib) = sub_integrand7( N0, yw((ib+1)/2), yw(mod(ib+1,2)+1) )
    enddo
  enddo

! -----------------------------------------------------------
! -------- integrate function_stiffness to obtain   ---------
! -------- conventional stiffness matrix            ---------
! -----------------------------------------------------------
  do ib = 1,4
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand0(0,ib),submatrix_I0(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand1k(0,ib),submatrix_I1k(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand1m(0,ib),submatrix_I1m(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand2(0,ib),submatrix_I2(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand3k(0,ib),submatrix_I3k(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand3m(0,ib),submatrix_I3m(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand4(0,ib),submatrix_I4(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand5k(0,ib),submatrix_I5k(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand5m(0,ib),submatrix_I5m(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand6(0,ib),submatrix_I6(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand7(0,ib),submatrix_I7(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
  enddo

  else

  if(USE_GAUSS_INTEGRATION) then
    istart = 1
    iend = ns
    A = (grid_r(ir+1) - grid_r(ir)) / 2.d0
  else
    istart = 0
    iend = 2*ns
!! DK DK this is unused in the case of Simpson
    A = 0.d0
  endif

  do iy = istart, iend
    if(USE_GAUSS_INTEGRATION) then
!! DK DK map Gauss points, which are in ]-1,+1[, to ]0,1[
      y = (xi_Gauss(iy) + 1.d0) / 2.d0
    else
      y = dble(iy) / dble(2*ns)
    endif
!! DK DK this is the interpolated radius
    r = ( 1.d0 - y ) * grid_r(ir) + y * grid_r(ir+1)
    call grid_interpolate_3(r,rho0,lambda0,rho_structure_zone(1:4,izone1),vpv_structure_zone(1:4,izone1),vph_structure_zone(1:4,izone1),vsv_structure_zone(1:4,izone1),vsh_structure_zone(1:4,izone1),eta_structure_zone(1:4,izone1) )
!    rho0 = grid_interpolate(y,grid_rho(1,ir), grid_rho(2,ir) )
!    lambda0 = grid_interpolate(y,grid_kappa(1,ir), grid_kappa(2,ir) )
    yw(1) = 1.d0 - y
    yw(2) = y
    dyw(1) = - 1.d0 / ( grid_r(ir+1) - grid_r(ir) )
    dyw(2) =   1.d0 / ( grid_r(ir+1) - grid_r(ir) )
    do ib=1,4
      val_integrand0(iy,ib) = sub_integrandF0( lambda0, yw((ib+1)/2), yw(mod(ib+1,2)+1), r )
      val_integrand1k(iy,ib) = sub_integrandF1( rho0, yw((ib+1)/2), yw(mod(ib+1,2)+1) )
      val_integrand2(iy,ib) = sub_integrandF2( rho0, dyw((ib+1)/2), dyw(mod(ib+1,2)+1), r )
    enddo
  enddo

! -----------------------------------------------------------
! -------- integrate function_stiffness to obtain   ---------
! -------- conventional stiffness matrix            ---------
! -----------------------------------------------------------
  do ib = 1,4
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand0(0,ib),submatrix_I0(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand1k(0,ib),submatrix_I1k(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    submatrix_I1m(ib,ir) = 0.d0
    call simpson_or_gauss( ns,grid_r(ir),grid_r(ir+1), val_integrand2(0,ib),submatrix_I2(ib,ir), &
                           weight_Gauss,A,USE_GAUSS_INTEGRATION )
    submatrix_I3k(ib,ir) = 0.d0
    submatrix_I3m(ib,ir) = 0.d0
    submatrix_I4(ib,ir) = 0.d0
    submatrix_I5k(ib,ir) = 0.d0
    submatrix_I5m(ib,ir) = 0.d0
    submatrix_I6(ib,ir) = 0.d0
    submatrix_I7(ib,ir) = 0.d0
  enddo

  endif

  enddo ! of do ir=1,ngrid_r-1


  contains

! functions for numerical integrations

  real(kind=8) function sub_integrand0(rho0,y1,y2,r)
  implicit none
    real(kind=8) :: rho0,y1,y2,r
    sub_integrand0 = y1 * r * rho0 * y2 * r
  end function sub_integrand0

  real(kind=8) function sub_integrand1(C0,dy1,dy2,r)
  implicit none
    real(kind=8) :: C0,dy1,dy2,r
    sub_integrand1 = dy1 * r * C0 * dy2 * r
  end function sub_integrand1

  real(kind=8) function sub_integrand2(L0,dy1,dy2,r)
  implicit none
    real(kind=8) :: L0,dy1,dy2,r
    sub_integrand2 = dy1 * r * L0 * dy2 * r
  end function sub_integrand2

  real(kind=8) function sub_integrand3(F0,y1,dy2,r)
  implicit none
    real(kind=8) :: F0,y1,dy2,r
    sub_integrand3 = y1 * F0 * dy2 * r
  end function sub_integrand3

  real(kind=8) function sub_integrand4(L0,dy1,y2,r)
  implicit none
    real(kind=8) :: L0,dy1,y2,r
    sub_integrand4 = dy1 * r * L0 * y2
  end function sub_integrand4

  real(kind=8) function sub_integrand5(A0,y1,y2)
  implicit none
    real(kind=8) :: A0,y1,y2
    sub_integrand5 = y1 * A0 * y2
  end function sub_integrand5

  real(kind=8) function sub_integrand6(L0,y1,y2)
  implicit none
    real(kind=8) :: L0,y1,y2
    sub_integrand6 = y1 * L0 * y2
  end function sub_integrand6

  real(kind=8) function sub_integrand7(N0,y1,y2)
  implicit none
    real(kind=8) :: N0,y1,y2
    sub_integrand7 = y1 * N0 * y2
  end function sub_integrand7

  real(kind=8) function sub_integrandF0(lambda0,y1,y2,r)
  implicit none
    real(kind=8) :: lambda0,y1,y2,r
    sub_integrandF0 = y1 * r / lambda0 * y2 * r
  end function sub_integrandF0

  real(kind=8) function sub_integrandF1(rho0,y1,y2)
  implicit none
    real(kind=8) :: rho0,y1,y2
    sub_integrandF1 = y1 / rho0 * y2
  end function sub_integrandF1

  real(kind=8) function sub_integrandF2(rho0,dy1,dy2,r)
  implicit none
    real(kind=8) :: rho0,dy1,dy2,r
    sub_integrandF2 = dy1 * r / rho0 * dy2 * r
  end function sub_integrandF2

  end subroutine comp_submatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function for numerical integrations

  real(kind=8) function grid_interpolate(y,grid_val1,grid_val2)

    implicit none

    real(kind=8) y,grid_val1,grid_val2

    grid_interpolate = grid_val1 + y * ( grid_val2 - grid_val1 )

  end function grid_interpolate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simpson_or_gauss( ns,xs,xe,f,integ,weight_Gauss,A,USE_GAUSS_INTEGRATION )

! integration using Simpson's (also known as Newton-Cotes) numerical quadrature rule based on a set of evenly-spaced points.

!    ns        :  integer   number of grid intervals
!    xs,xe     :  real(kind=8)    start and the end of the integration
!    f(0:2*ns) :  real(kind=8)    function values at the nodes
!    integ     :  real(kind=8)    the integrated value

!                                                June 1992, N. Takeuchi

!! DK DK Dimitri Komatitsch, CNRS, Marseille, France, March 2014: added support for more precise Gauss quadrature

  implicit none

  integer ns
  real(kind=8) xe,xs,f(0:2*ns),integ
  integer i,nn
  real(kind=8) dx,intego,intege

!! DK DK added this for Gauss integration
  logical :: USE_GAUSS_INTEGRATION
  real(kind=8), dimension(ns) :: weight_Gauss
  real(kind=8) :: A

!! DK DK Dimitri Komatitsch, CNRS, Marseille, France, March 2014: added support for Gauss integration instead of Simpson

!! DK DK the tables of Gauss-Legendre points and weights are taken from Pavel Holoborodko
!! DK DK at http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/

! Gauss-Legendre n-point quadrature, exact for all polynomials of degree <= 2n-1
!
! When n is even (which is always the case in this code):
!
! int(f(t), t = a..b) = A * sum(weight(i) * f(A*xi(i) + B), i = 1..ns)
!
!  where A = (b-a) / 2  and  B = (a+b) / 2
!

!! DK DK added this for Gauss integration
 if(USE_GAUSS_INTEGRATION) then

  integ = weight_Gauss(1) * f(1)

  do i = 2,ns
    integ = integ + weight_Gauss(i) * f(i)
  enddo

  integ = integ * A

 else ! Simpson (also known as Newton-Cotes) integration based on a set of evenly-spaced points

  nn=2*ns

  integ = f(0)+f(nn)

  intego = 0.d0
  intege = 0.d0

  do i=1,ns-1
    nn = 2 * i - 1
    intego = intego + f(nn)
    nn = 2 * i
    intege = intege + f(nn)
  enddo

  nn = 2 * ns - 1
  intego = intego + f(nn)
  dx = ( xe - xs ) / ( 2.d0 * ns )
  integ = ( integ + 4.d0 * intego + 2.0 * intege ) * dx / 3.d0

 endif

  end

 subroutine compute_coef_cubic_poly(x,y,coef )
!this subroutine is used for computing the cubic polynomial coefficients for every parameters
 implicit none

 real(kind=8) x(2385),y(2385),coef(4),coef_check(4)
 integer i,j,k
 real(kind=8) rmax,a(0:3),b(0:3,0:3),bb(0:3),denominator,inv(0:3,0:3),rr,param,param1

 rmax = 6371.d0
 a(0:3) = x(1:4) / rmax
 bb(0:3) = y(1:4)
 inv(:,:) = 0.0d0
! write(*,*) 'a and bb for check',a,bb
! inv is the inverse of the 3 order vandermonde's matrix

j = 0
denominator = (a(1)-a(0))*(a(2)-a(0))*(a(3)-a(0)) !j=0 case
b(0,j) = (a(1)*a(2)*a(3)) / denominator
b(1,j) = ((-1)**1) * ( a(1)*a(2) + a(1)*a(3) + a(2)*a(3)  ) / denominator
b(2,j) = ((-1)**2) * ( a(1) + a(2) + a(3)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator

j = 1
denominator = (a(0)-a(1))*(a(2)-a(1))*(a(3)-a(1)) !j=1 case
b(0,j) = (a(0)*a(2)*a(3)) / denominator
b(1,j) = ((-1)**1) * ( a(0)*a(2) + a(0)*a(3) + a(2)*a(3)  ) / denominator
b(2,j) = ((-1)**2) * ( a(0) + a(2) + a(3)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator

j = 2
denominator = (a(0)-a(2))*(a(1)-a(2))*(a(3)-a(2)) !j=2 case
b(0,j) = (a(0)*a(1)*a(3)) / denominator
b(1,j) = ((-1)**1) * ( a(0)*a(1) + a(0)*a(3) + a(1)*a(3)  ) / denominator
b(2,j) = ((-1)**2) * ( a(0) + a(1) + a(3)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator

j = 3
denominator = (a(0)-a(3))*(a(1)-a(3))*(a(2)-a(3)) !j=3 case
b(0,j) = (a(0)*a(1)*a(2)) / denominator
b(1,j) = ((-1)**1) * ( a(0)*a(1) + a(0)*a(2) + a(1)*a(2)  ) / denominator
b(2,j) = ((-1)**2) * ( a(0) + a(1) + a(2)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator


do i=0,3
 do j=0,3
  do k=0,3
  inv(i,j)=inv(i,j) + (a(k)**j) * b(i,k)
  enddo
 enddo
enddo

write(*,*) 'check',inv
write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

coef(:)=0.0d0
do i=0,3
 do k=0,3
  coef(i+1) = coef(i+1) + b(i,k)*bb(k)
 end do
end do

a(0:3) = x(2382:2385) / rmax
bb(0:3) = y(2382:2385)
!write(*,*) 'a and bb for check_2nd',a,bb
! inv is the inverse of the 3 order vandermonde's matrix

j = 0
denominator = (a(1)-a(0))*(a(2)-a(0))*(a(3)-a(0)) !j=0 case
b(0,j) = (a(1)*a(2)*a(3)) / denominator
b(1,j) = ((-1)**1) * ( a(1)*a(2) + a(1)*a(3) + a(2)*a(3)  ) / denominator
b(2,j) = ((-1)**2) * ( a(1) + a(2) + a(3)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator

j = 1
denominator = (a(0)-a(1))*(a(2)-a(1))*(a(3)-a(1)) !j=1 case
b(0,j) = (a(0)*a(2)*a(3)) / denominator
b(1,j) = ((-1)**1) * ( a(0)*a(2) + a(0)*a(3) + a(2)*a(3)  ) / denominator
b(2,j) = ((-1)**2) * ( a(0) + a(2) + a(3)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator

j = 2
denominator = (a(0)-a(2))*(a(1)-a(2))*(a(3)-a(2)) !j=2 case
b(0,j) = (a(0)*a(1)*a(3)) / denominator
b(1,j) = ((-1)**1) * ( a(0)*a(1) + a(0)*a(3) + a(1)*a(3)  ) / denominator
b(2,j) = ((-1)**2) * ( a(0) + a(1) + a(3)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator

j = 3
denominator = (a(0)-a(3))*(a(1)-a(3))*(a(2)-a(3)) !j=3 case
b(0,j) = (a(0)*a(1)*a(2)) / denominator
b(1,j) = ((-1)**1) * ( a(0)*a(1) + a(0)*a(2) + a(1)*a(2)  ) / denominator
b(2,j) = ((-1)**2) * ( a(0) + a(1) + a(2)  ) / denominator
b(3,j) = ((-1)**3) * 1.0d0 / denominator


coef_check(:)=0.0d0
do i=0,3
 do k=0,3
   coef_check(i+1) = coef_check(i+1) + b(i,k)*bb(k)
 end do
end do

write(*,*) 'coef',coef
write(*,*) 'coef_check',coef_check
         
! open(1,file = './data/compare_integration', status = 'unknown', form = 'formatted')
! do k=1,2385
!    rr = x(k) / rmax
!    param = coef(1) + coef(2)*rr + coef(3)*rr**2 + coef(4)*rr**3
!    param1 = coef_check(1) + coef_check(2)*rr + coef_check(3)*rr**2 + coef_check(4)*rr**3
!    write(1,"(f12.6,f14.7,f14.7,f14.7)") x(k),param,param1,y(k)
! ! write(*,*) x(k),param/y(k)-1.0d0,param1/y(k)-1.0d0
! enddo
! close(1)

 end subroutine

subroutine grid_interpolate_1(r,rho0,Fk0,Fm0,Ak0,Am0,L0,N0,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone )
  implicit none
    real(kind=8) rho_structure_zone(4)
    real(kind=8) vpv_structure_zone(4),vph_structure_zone(4)
    real(kind=8) vsv_structure_zone(4),vsh_structure_zone(4)
    real(kind=8) eta_structure_zone(4)
    real(kind=8) r,rho0,Fk0,Fm0,Ak0,Am0,L0,N0
    real(kind=8) rho_r,vpv_r,vph_r,vsv_r,vsh_r,eta_r
    real(kind=8) A_r,C_r,L_r,N_r,F_r
    real(kind=8) kappa_r,mu_r
    real(kind=8), external :: cal_PREM_structure

      rho_r = cal_PREM_structure( r, rho_structure_zone(1) )
      vpv_r = cal_PREM_structure( r, vpv_structure_zone(1) )
      vph_r = cal_PREM_structure( r, vph_structure_zone(1) )
      vsv_r = cal_PREM_structure( r, vsv_structure_zone(1) )
      vsh_r = cal_PREM_structure( r, vsh_structure_zone(1) )
      eta_r = cal_PREM_structure( r, eta_structure_zone(1) )

  A_r = rho_r * vph_r * vph_r
  C_r = rho_r * vpv_r * vpv_r
  L_r = rho_r * vsv_r * vsv_r
  N_r = rho_r * vsh_r * vsh_r
  F_r = eta_r * ( A_r - 2.d0 * L_r )

  kappa_r = ( 4.d0 * A_r + C_r + 4.d0 * F_r - 4.d0 * N_r ) / 9.d0

  mu_r = ( A_r + C_r - 2.d0 * F_r + 5.d0 * N_r + 6.d0 * L_r ) / 15.d0

 rho0 = rho_r

      Ak0 = A_r * kappa_r / ( kappa_r + 4.d0/3.d0 * mu_r )
      Am0 = A_r * 4.d0/3.d0 * mu_r / ( kappa_r + 4.d0/3.d0 * mu_r )
!      Ck0 = C_r * kappa_r / ( kappa_r + 4.d0/3.d0 * mu_r )
!      Cm0 = C_r * 4.d0/3.d0 * mu_r / ( kappa_r + 4.d0/3.d0 * mu_r )
      L0 = L_r
      N0 = N_r
      Fk0 = F_r * kappa_r / ( kappa_r - 2.d0/3.d0 * mu_r )
      Fm0 = F_r * (-2.d0/3.d0) * mu_r / ( kappa_r-2.d0/3.d0 * mu_r )

!      kappa0 = kappa_r
!      mu0 = mu_r
!rho0,Fk0,Fm0,Ak0,Am0,L0,N0

end subroutine

subroutine grid_interpolate_2(r,rho0,Ck0,Cm0,L0,Fk0,Fm0,Ak0,Am0,N0,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone )
  implicit none
    real(kind=8) rho_structure_zone(4)
    real(kind=8) vpv_structure_zone(4),vph_structure_zone(4)
    real(kind=8) vsv_structure_zone(4),vsh_structure_zone(4)
    real(kind=8) eta_structure_zone(4)
    real(kind=8) r,rho0,Ck0,Cm0,L0,Fk0,Fm0,Ak0,Am0,N0
    real(kind=8) rho_r,vpv_r,vph_r,vsv_r,vsh_r,eta_r
    real(kind=8) A_r,C_r,L_r,N_r,F_r
    real(kind=8) kappa_r,mu_r
    real(kind=8), external :: cal_PREM_structure

      rho_r = cal_PREM_structure( r, rho_structure_zone(1) )
      vpv_r = cal_PREM_structure( r, vpv_structure_zone(1) )
      vph_r = cal_PREM_structure( r, vph_structure_zone(1) )
      vsv_r = cal_PREM_structure( r, vsv_structure_zone(1) )
      vsh_r = cal_PREM_structure( r, vsh_structure_zone(1) )
      eta_r = cal_PREM_structure( r, eta_structure_zone(1) )

  A_r = rho_r * vph_r * vph_r
  C_r = rho_r * vpv_r * vpv_r
  L_r = rho_r * vsv_r * vsv_r
  N_r = rho_r * vsh_r * vsh_r
  F_r = eta_r * ( A_r - 2.d0 * L_r )

  kappa_r = ( 4.d0 * A_r + C_r + 4.d0 * F_r - 4.d0 * N_r ) / 9.d0

  mu_r = ( A_r + C_r - 2.d0 * F_r + 5.d0 * N_r + 6.d0 * L_r ) / 15.d0

 rho0 = rho_r

      Ak0 = A_r * kappa_r / ( kappa_r + 4.d0/3.d0 * mu_r )
      Am0 = A_r * 4.d0/3.d0 * mu_r / ( kappa_r + 4.d0/3.d0 * mu_r )
      Ck0 = C_r * kappa_r / ( kappa_r + 4.d0/3.d0 * mu_r )
      Cm0 = C_r * 4.d0/3.d0 * mu_r / ( kappa_r + 4.d0/3.d0 * mu_r )
      L0 = L_r
      N0 = N_r
      Fk0 = F_r * kappa_r / ( kappa_r - 2.d0/3.d0 * mu_r )
      Fm0 = F_r * (-2.d0/3.d0) * mu_r / ( kappa_r-2.d0/3.d0 * mu_r )
!      kappa0 = kappa_r
!      mu0 = mu_r
!rho0,Ck0,Cm0,L0,Fk0,Fm0,Ak0,Am0,N0

end subroutine

subroutine grid_interpolate_3(r,rho0,lambda0,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone )
  implicit none
    real(kind=8) rho_structure_zone(4)
    real(kind=8) vpv_structure_zone(4),vph_structure_zone(4)
    real(kind=8) vsv_structure_zone(4),vsh_structure_zone(4)
    real(kind=8) eta_structure_zone(4)
    real(kind=8) r,rho0,lambda0
    real(kind=8) rho_r,vpv_r,vph_r,vsv_r,vsh_r,eta_r
    real(kind=8) A_r,C_r,L_r,N_r,F_r
    real(kind=8) kappa_r,mu_r
    real(kind=8), external :: cal_PREM_structure

      rho_r = cal_PREM_structure( r, rho_structure_zone(1) )
      vpv_r = cal_PREM_structure( r, vpv_structure_zone(1) )
      vph_r = cal_PREM_structure( r, vph_structure_zone(1) )
      vsv_r = cal_PREM_structure( r, vsv_structure_zone(1) )
      vsh_r = cal_PREM_structure( r, vsh_structure_zone(1) )
      eta_r = cal_PREM_structure( r, eta_structure_zone(1) )

  A_r = rho_r * vph_r * vph_r
  C_r = rho_r * vpv_r * vpv_r
  L_r = rho_r * vsv_r * vsv_r
  N_r = rho_r * vsh_r * vsh_r
  F_r = eta_r * ( A_r - 2.d0 * L_r )

  kappa_r = ( 4.d0 * A_r + C_r + 4.d0 * F_r - 4.d0 * N_r ) / 9.d0

!  mu_r = ( A_r + C_r - 2.d0 * F_r + 5.d0 * N_r + 6.d0 * L_r ) / 15.d0

 rho0 = rho_r

!      Ak0 = A_r * kappa_r / ( kappa_r + 4.d0/3.d0 * mu_r )
!      Am0 = A_r * 4.d0/3.d0 * mu_r / ( kappa_r + 4.d0/3.d0 * mu_r )
!      Ck0 = C_r * kappa_r / ( kappa_r + 4.d0/3.d0 * mu_r )
!      Cm0 = C_r * 4.d0/3.d0 * mu_r / ( kappa_r + 4.d0/3.d0 * mu_r )
!      L0 = L_r
!      N0 = N_r
!      Fk0 = F_r * kappa_r / ( kappa_r - 2.d0/3.d0 * mu_r )
!      Fm0 = F_r * (-2.d0/3.d0) * mu_r / ( kappa_r-2.d0/3.d0 * mu_r )
      lambda0 = kappa_r
!      mu0 = mu_r
!rho0,lambda0

end subroutine


