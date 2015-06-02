cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_submatrix_mod
     &	  ( ngrid_r,grid_r,
     &	    grid_rho,grid_kappa,grid_mu,
     &	    grid_Ak,grid_Am,grid_L,grid_N,grid_Fk,grid_Fm,
     &	    submatrix_I0,submatrix_I1k,submatrix_I3k,
     &	    submatrix_I3m,submatrix_I4,submatrix_I5k,
     &	    submatrix_I5m,submatrix_I6,submatrix_I7,
     &	    submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing modified submatrices for elastic part of the medium.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer ngrid_r
	real*8 grid_r(*)
	real*8 grid_rho(2,*),grid_kappa(2,*),grid_mu(2,*)
	real*8 grid_Ak(2,*),grid_Am(2,*),grid_L(2,*)
	real*8 grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
	real*8 submatrix_I0(4,*),submatrix_I1k(4,*)
	real*8 submatrix_I3k(4,*),submatrix_I3m(4,*),submatrix_I4(4,*)
	real*8 submatrix_I5k(4,*),submatrix_I5m(4,*)
	real*8 submatrix_I6(4,*),submatrix_I7(4,*)
	real*8 submatrix_I3k_mod(6,*),submatrix_I3m_mod(6,*)
	real*8 submatrix_I4_mod(6,*)
c variables for numerical integrations
	integer ns
c	parameter ( ns = 6 )
	parameter ( ns = 128 )
	integer ir,iy
	real*8 y,r,dr,rhor(0:2*ns),laminvr(0:2*ns),rhoinv(0:2*ns)
	real*8 Fkr(0:2*ns),Fmr(0:2*ns),Lr(0:2*ns)
	real*8 Ak(0:2*ns),Am(0:2*ns),L(0:2*ns),N(0:2*ns)
	real*8 rhor_ave,laminvr_ave,rhoinv_ave
	real*8 Fkr_ave,Fmr_ave,Lr_ave,Ak_ave,Am_ave,L_ave,N_ave
c functions for numerical integrations
	real*8 grid_interpolate
	real*8 grid_val1,grid_val2
	grid_interpolate(y,grid_val1,grid_val2)
     &    = grid_val1
     &	    + y * ( grid_val2 - grid_val1 )
c **********************************************************************
c  computing submatrix elements for each cell
c **********************************************************************
	do 990 ir=1,ngrid_r-1
	dr = grid_r(ir+1) - grid_r(ir)
	if ( grid_mu(1,ir)*grid_mu(2,ir).ne.0.d0 ) then
c evaluating the density or elastic constatns at each sub-grid
	  do 100 iy=0,2*ns
	    y = dble(iy) / dble(2*ns)
	    r = ( 1.d0 - y ) * grid_r(ir)
     &	        + y * grid_r(ir+1)
	    rhor(iy)
     &	      = grid_interpolate(y,grid_rho(1,ir),
     &	                           grid_rho(2,ir)  )
     &	        * r * r
	    Fkr(iy)
     &	      = grid_interpolate(y,grid_Fk(1,ir),
     &	                           grid_Fk(2,ir)  )
     &	        * r
	    Fmr(iy)
     &	      = grid_interpolate(y,grid_Fm(1,ir),
     &	                           grid_Fm(2,ir)  )
     &	        * r
	    Lr(iy)
     &	      = grid_interpolate(y,grid_L(1,ir),
     &	                           grid_L(2,ir)  )
     &	        * r
	    Ak(iy)
     &	      = grid_interpolate(y,grid_Ak(1,ir),
     &	                           grid_Ak(2,ir)  )
	    Am(iy)
     &	      = grid_interpolate(y,grid_Am(1,ir),
     &	                           grid_Am(2,ir)  )
	    L(iy)
     &	      = grid_interpolate(y,grid_L(1,ir),
     &	                           grid_L(2,ir)  )
	    N(iy)
     &	      = grid_interpolate(y,grid_N(1,ir),
     &	                           grid_N(2,ir)  )
  100	  continue
c computing lumped parameters
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                rhor,rhor_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                Fkr,Fkr_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                Fmr,Fmr_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                Lr,Lr_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                Ak,Ak_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                Am,Am_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                L,L_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                N,N_ave )
	  rhor_ave = rhor_ave / dr
	  Fkr_ave = Fkr_ave / dr
	  Fmr_ave = Fmr_ave / dr
	  Lr_ave = Lr_ave / dr
	  Ak_ave = Ak_ave / dr
	  Am_ave = Am_ave / dr
	  L_ave = L_ave / dr
	  N_ave = N_ave / dr
c -----------------------------------------------------------
c -------- computing the optimally accurate         ---------
c -------- submatrices                              ---------
c -----------------------------------------------------------
c modified submatrix_I0
	  submatrix_I0(1,ir)
     &	      = submatrix_I0(1,ir) + dr / 12.d0 * rhor_ave
	  submatrix_I0(2,ir)
     &	      = submatrix_I0(2,ir) - dr / 12.d0 * rhor_ave
	  submatrix_I0(3,ir)
     &	      = submatrix_I0(3,ir) - dr / 12.d0 * rhor_ave
	  submatrix_I0(4,ir)
     &	      = submatrix_I0(4,ir) + dr / 12.d0 * rhor_ave
c modified submatrix_I3k_mod
	  submatrix_I3k_mod(1,ir)
     &	      = submatrix_I3k(1,ir) - 1.d0 / 12.d0 * Fkr_ave
	  submatrix_I3k_mod(2,ir)
     &	      = submatrix_I3k(2,ir) + 2.d0 / 12.d0 * Fkr_ave
	  submatrix_I3k_mod(3,ir)
     &	      = submatrix_I3k(3,ir) + 1.d0 / 12.d0 * Fkr_ave
	  submatrix_I3k_mod(4,ir)
     &	      = submatrix_I3k(4,ir) - 2.d0 / 12.d0 * Fkr_ave
	  submatrix_I3k_mod(5,ir) = - 1.d0 / 12.d0 * Fkr_ave
	  submatrix_I3k_mod(6,ir) =   1.d0 / 12.d0 * Fkr_ave
c modified submatrix_I3m_mod
	  submatrix_I3m_mod(1,ir)
     &	      = submatrix_I3m(1,ir) - 1.d0 / 12.d0 * Fmr_ave
	  submatrix_I3m_mod(2,ir)
     &	      = submatrix_I3m(2,ir) + 2.d0 / 12.d0 * Fmr_ave
	  submatrix_I3m_mod(3,ir)
     &	      = submatrix_I3m(3,ir) + 1.d0 / 12.d0 * Fmr_ave
	  submatrix_I3m_mod(4,ir)
     &	      = submatrix_I3m(4,ir) - 2.d0 / 12.d0 * Fmr_ave
	  submatrix_I3m_mod(5,ir) = - 1.d0 / 12.d0 * Fmr_ave
	  submatrix_I3m_mod(6,ir) =   1.d0 / 12.d0 * Fmr_ave
c modified submatrix_I4_mod
	  submatrix_I4_mod(1,ir)
     &	      = submatrix_I4(1,ir) + 1.d0 / 12.d0 * Lr_ave
	  submatrix_I4_mod(2,ir)
     &	      = submatrix_I4(2,ir) - 2.d0 / 12.d0 * Lr_ave
	  submatrix_I4_mod(3,ir)
     &	      = submatrix_I4(3,ir) - 1.d0 / 12.d0 * Lr_ave
	  submatrix_I4_mod(4,ir)
     &	      = submatrix_I4(4,ir) + 2.d0 / 12.d0 * Lr_ave
	  submatrix_I4_mod(5,ir) =   1.d0 / 12.d0 * Lr_ave
	  submatrix_I4_mod(6,ir) = - 1.d0 / 12.d0 * Lr_ave
c modified submatrix_I5k
	  submatrix_I5k(1,ir)
     &	      = submatrix_I5k(1,ir) + dr / 12.d0 * Ak_ave
	  submatrix_I5k(2,ir)
     &	      = submatrix_I5k(2,ir) - dr / 12.d0 * Ak_ave
	  submatrix_I5k(3,ir)
     &	      = submatrix_I5k(3,ir) - dr / 12.d0 * Ak_ave
	  submatrix_I5k(4,ir)
     &	      = submatrix_I5k(4,ir) + dr / 12.d0 * Ak_ave
c modified submatrix_I5m
	  submatrix_I5m(1,ir)
     &	      = submatrix_I5m(1,ir) + dr / 12.d0 * Am_ave
	  submatrix_I5m(2,ir)
     &	      = submatrix_I5m(2,ir) - dr / 12.d0 * Am_ave
	  submatrix_I5m(3,ir)
     &	      = submatrix_I5m(3,ir) - dr / 12.d0 * Am_ave
	  submatrix_I5m(4,ir)
     &	      = submatrix_I5m(4,ir) + dr / 12.d0 * Am_ave
c modified submatrix_I6
	  submatrix_I6(1,ir)
     &	      = submatrix_I6(1,ir) + dr / 12.d0 * L_ave
	  submatrix_I6(2,ir)
     &	      = submatrix_I6(2,ir) - dr / 12.d0 * L_ave
	  submatrix_I6(3,ir)
     &	      = submatrix_I6(3,ir) - dr / 12.d0 * L_ave
	  submatrix_I6(4,ir)
     &	      = submatrix_I6(4,ir) + dr / 12.d0 * L_ave
c modified submatrix_I7
	  submatrix_I7(1,ir)
     &	      = submatrix_I7(1,ir) + dr / 12.d0 * N_ave
	  submatrix_I7(2,ir)
     &	      = submatrix_I7(2,ir) - dr / 12.d0 * N_ave
	  submatrix_I7(3,ir)
     &	      = submatrix_I7(3,ir) - dr / 12.d0 * N_ave
	  submatrix_I7(4,ir)
     &	      = submatrix_I7(4,ir) + dr / 12.d0 * N_ave
	else
	  do 200 iy=0,2*ns
	    y = dble(iy) / dble(2*ns)
	    r = ( 1.d0 - y ) * grid_r(ir)
     &	        + y * grid_r(ir+1)
	    laminvr(iy)
     &	      = r * r
     &	        / grid_interpolate(y,grid_kappa(1,ir),
     &	                             grid_kappa(2,ir)  )
	    rhoinv(iy)
     &	      = 1.d0 / grid_interpolate(y,grid_rho(1,ir),
     &	                                  grid_rho(2,ir)  )
  200	  continue
c computing lumped parameters
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                laminvr,laminvr_ave )
	  call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                rhoinv,rhoinv_ave )
	  laminvr_ave = laminvr_ave / dr
	  rhoinv_ave = rhoinv_ave / dr
c -----------------------------------------------------------
c -------- computing the optimally accurate         ---------
c -------- submatrices                              ---------
c -----------------------------------------------------------
c modified submatrix_I0
	  submatrix_I0(1,ir)
     &	      = submatrix_I0(1,ir) + dr / 12.d0 * laminvr_ave
	  submatrix_I0(2,ir)
     &	      = submatrix_I0(2,ir) - dr / 12.d0 * laminvr_ave
	  submatrix_I0(3,ir)
     &	      = submatrix_I0(3,ir) - dr / 12.d0 * laminvr_ave
	  submatrix_I0(4,ir)
     &	      = submatrix_I0(4,ir) + dr / 12.d0 * laminvr_ave
c modified submatrix_I1k
	  submatrix_I1k(1,ir)
     &	      = submatrix_I1k(1,ir) + dr / 12.d0 * rhoinv_ave
	  submatrix_I1k(2,ir)
     &	      = submatrix_I1k(2,ir) - dr / 12.d0 * rhoinv_ave
	  submatrix_I1k(3,ir)
     &	      = submatrix_I1k(3,ir) - dr / 12.d0 * rhoinv_ave
	  submatrix_I1k(4,ir)
     &	      = submatrix_I1k(4,ir) + dr / 12.d0 * rhoinv_ave
	endif
  990	continue
  	 write(*,*) 'comp_submatrix_mod function is over now'
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine comp_submatrix
     &	    ( ngrid_r,grid_r,
     &	      grid_rho,grid_kappa,grid_mu,grid_Ak,grid_Am,
     &	      grid_Ck,grid_Cm,grid_L,grid_N,grid_Fk,grid_Fm,
     &	      submatrix_I0,submatrix_I1k,submatrix_I1m,submatrix_I2,
     &	      submatrix_I3k,submatrix_I3m,submatrix_I4,submatrix_I5k,
     &	      submatrix_I5m,submatrix_I6,submatrix_I7 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing submatrices for elastic part of the medium.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
	integer ngrid_r
	real*8 grid_r(*)
	real*8 grid_rho(2,*),grid_kappa(2,*),grid_mu(2,*)
	real*8 grid_Ak(2,*),grid_Am(2,*),grid_Ck(2,*),grid_Cm(2,*)
	real*8 grid_L(2,*),grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
	real*8 submatrix_I0(4,*),submatrix_I1k(4,*),submatrix_I1m(4,*)
	real*8 submatrix_I2(4,*),submatrix_I3k(4,*),submatrix_I3m(4,*)
	real*8 submatrix_I4(4,*),submatrix_I5k(4,*),submatrix_I5m(4,*)
	real*8 submatrix_I6(4,*),submatrix_I7(4,*)
c variables for numerical integrations
	integer ns
c	parameter ( ns = 6 )
	parameter ( ns = 128 )
	integer ir,iy,ib
	real*8 y,yw(2),dyw(2),r
	real*8 rho0,lambda0,Ck0,Cm0,L0,Fk0,Fm0,Ak0,Am0,N0
	real*8 val_integrand0(0:2*ns,4)
	real*8 val_integrand1k(0:2*ns,4)
	real*8 val_integrand1m(0:2*ns,4)
	real*8 val_integrand2(0:2*ns,4)
	real*8 val_integrand3k(0:2*ns,4)
	real*8 val_integrand3m(0:2*ns,4)
	real*8 val_integrand4(0:2*ns,4)
	real*8 val_integrand5k(0:2*ns,4)
	real*8 val_integrand5m(0:2*ns,4)
	real*8 val_integrand6(0:2*ns,4)
	real*8 val_integrand7(0:2*ns,4)
c functions for numerical integrations
	real*8 grid_interpolate,C0,F0,A0
	real*8 sub_integrand0,sub_integrand1,sub_integrand2
	real*8 sub_integrand3,sub_integrand4,sub_integrand5
	real*8 sub_integrand6,sub_integrand7
	real*8 sub_integrandF0,sub_integrandF1,sub_integrandF2
	real*8 y1,dy1,y2,dy2
	real*8 grid_val1,grid_val2
	grid_interpolate(y,grid_val1,grid_val2)
     &    = grid_val1
     &	    + y * ( grid_val2 - grid_val1 )
	sub_integrand0(rho0,y1,dy1,y2,dy2,r)
     &	  = y1 * r * rho0 * y2 * r
	sub_integrand1(C0,y1,dy1,y2,dy2,r)
     &	  = dy1 * r * C0 * dy2 * r
	sub_integrand2(L0,y1,dy1,y2,dy2,r)
     &	  = dy1 * r * L0 * dy2 * r
	sub_integrand3(F0,y1,dy1,y2,dy2,r)
     &	  = y1 * F0 * dy2 * r
	sub_integrand4(L0,y1,dy1,y2,dy2,r)
     &	  = dy1 * r * L0 * y2
	sub_integrand5(A0,y1,dy1,y2,dy2,r)
     &	  = y1 * A0 * y2
	sub_integrand6(L0,y1,dy1,y2,dy2,r)
     &	  = y1 * L0 * y2
	sub_integrand7(N0,y1,dy1,y2,dy2,r)
     &	  = y1 * N0 * y2
	sub_integrandF0(lambda0,y1,dy1,y2,dy2,r)
     &	  = y1 * r / lambda0 * y2 * r
	sub_integrandF1(rho0,y1,dy1,y2,dy2,r)
     &	  = y1 / rho0 * y2
	sub_integrandF2(rho0,y1,dy1,y2,dy2,r)
     &	  = dy1 * r / rho0 * dy2 * r
c **********************************************************************
c  computing stiffness submatrix elements for each cell
c **********************************************************************
	do 990 ir=1,ngrid_r-1
	if ( grid_mu(1,ir)*grid_mu(2,ir).ne.0.d0 ) then
c -----------------------------------------------------------
c --------  evaluating the stiffness and            ---------
c --------  the integrated function for stiffness   ---------
c --------  matrix at the sub-grids                 ---------
c -----------------------------------------------------------
	  do 110 iy=0,2*ns
	    y = dble(iy) / dble(2*ns)
	    r = ( 1.d0 - y ) * grid_r(ir)
     &	        + y * grid_r(ir+1)
	    rho0 = grid_interpolate(y,grid_rho(1,ir),grid_rho(2,ir) )
	    Ck0 = grid_interpolate(y,grid_Ck(1,ir),grid_Ck(2,ir) )
	    Cm0 = grid_interpolate(y,grid_Cm(1,ir),grid_Cm(2,ir) )
	    L0 = grid_interpolate(y,grid_L(1,ir),grid_L(2,ir) )
	    Fk0 = grid_interpolate(y,grid_Fk(1,ir),grid_Fk(2,ir) )
	    Fm0 = grid_interpolate(y,grid_Fm(1,ir),grid_Fm(2,ir) )
	    Ak0 = grid_interpolate(y,grid_Ak(1,ir),grid_Ak(2,ir) )
	    Am0 = grid_interpolate(y,grid_Am(1,ir),grid_Am(2,ir) )
	    N0 = grid_interpolate(y,grid_N(1,ir),grid_N(2,ir) )
	    yw(1) = 1.d0 - y
	    yw(2) = y
	    dyw(1) = - 1.d0 / ( grid_r(ir+1)
     &	                        - grid_r(ir) )
	    dyw(2) =   1.d0 / ( grid_r(ir+1)
     &	                        - grid_r(ir) )
	    do 100 ib=1,4
	      val_integrand0(iy,ib)
     &	        = sub_integrand0( rho0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand1k(iy,ib)
     &	        = sub_integrand1( Ck0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand1m(iy,ib)
     &	        = sub_integrand1( Cm0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand2(iy,ib)
     &	        = sub_integrand2( L0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand3k(iy,ib)
     &	        = sub_integrand3( Fk0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand3m(iy,ib)
     &	        = sub_integrand3( Fm0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand4(iy,ib)
     &	        = sub_integrand4( L0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand5k(iy,ib)
     &	        = sub_integrand5( Ak0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand5m(iy,ib)
     &	        = sub_integrand5( Am0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand6(iy,ib)
     &	        = sub_integrand6( L0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
	      val_integrand7(iy,ib)
     &	        = sub_integrand7( N0,
     &	                          yw((ib+1)/2),dyw((ib+1)/2),
     &	                          yw(mod(ib+1,2)+1),
     &	                          dyw(mod(ib+1,2)+1),
     &	                          r )
  100	    continue
  110	  continue
c -----------------------------------------------------------
c -------- integrating function_stiffness to obtain ---------
c -------- conventional stiffness matrix            ---------
c -----------------------------------------------------------
	  do 120 ib=1,4
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand0(0,ib),submatrix_I0(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand1k(0,ib),submatrix_I1k(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand1m(0,ib),submatrix_I1m(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand2(0,ib),submatrix_I2(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand3k(0,ib),submatrix_I3k(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand3m(0,ib),submatrix_I3m(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand4(0,ib),submatrix_I4(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand5k(0,ib),submatrix_I5k(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand5m(0,ib),submatrix_I5m(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand6(0,ib),submatrix_I6(ib,ir) )
c         write(*,*) 'simpson interal over 3',ib
c	    call simpson( ns,grid_r(ir),grid_r(ir+1),
c     &	                  val_integrand7(0,ib),submatrix_I7(ib,ir) )
  120	  continue
         if(mod(ir,100) == 1)   write(*,*) 'index for the grid:',ir
	else
	  do 210 iy=0,2*ns
	    y = dble(iy) / dble(2*ns)
	    r = ( 1.d0 - y ) * grid_r(ir)
     &	        + y * grid_r(ir+1)
	    rho0 = grid_interpolate(y,grid_rho(1,ir),
     &	                              grid_rho(2,ir) )
	    lambda0 = grid_interpolate(y,grid_kappa(1,ir),
     &	                                 grid_kappa(2,ir) )
	    yw(1) = 1.d0 - y
	    yw(2) = y
	    dyw(1) = - 1.d0 / ( grid_r(ir+1) - grid_r(ir) )
	    dyw(2) =   1.d0 / ( grid_r(ir+1) - grid_r(ir) )
	    do 200 ib=1,4
	      val_integrand0(iy,ib)
     &	        = sub_integrandF0( lambda0,
     &	                           yw((ib+1)/2),dyw((ib+1)/2),
     &	                           yw(mod(ib+1,2)+1),
     &	                           dyw(mod(ib+1,2)+1),
     &	                           r )
	      val_integrand1k(iy,ib)
     &	        = sub_integrandF1( rho0,
     &	                           yw((ib+1)/2),dyw((ib+1)/2),
     &	                           yw(mod(ib+1,2)+1),
     &	                           dyw(mod(ib+1,2)+1),
     &	                           r )
	      val_integrand2(iy,ib)
     &	        = sub_integrandF2( rho0,
     &	                           yw((ib+1)/2),dyw((ib+1)/2),
     &	                           yw(mod(ib+1,2)+1),
     &	                           dyw(mod(ib+1,2)+1),
     &	                           r )
  200	    continue
  210	  continue
c -----------------------------------------------------------
c -------- integrating function_stiffness to obtain ---------
c -------- conventional stiffness matrix            ---------
c -----------------------------------------------------------
	  do 220 ib=1,4
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand0(0,ib),submatrix_I0(ib,ir) )
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand1k(0,ib),submatrix_I1k(ib,ir) )
	    submatrix_I1m(ib,ir) = 0.d0
	    call simpson( ns,grid_r(ir),grid_r(ir+1),
     &	                  val_integrand2(0,ib),submatrix_I2(ib,ir) )
	    submatrix_I3k(ib,ir) = 0.d0
	    submatrix_I3m(ib,ir) = 0.d0
	    submatrix_I4(ib,ir) = 0.d0
	    submatrix_I5k(ib,ir) = 0.d0
	    submatrix_I5m(ib,ir) = 0.d0
	    submatrix_I6(ib,ir) = 0.d0
	    submatrix_I7(ib,ir) = 0.d0
  220	  continue
	endif
c
  990	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine simpson( ns,xs,xe,f,integ )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c integralation by using the Simpsons numerical integral.
c   required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine simpson(ns,f,integ)
c    ns        :  integer   number of grid intervals
c    xs,xe     :  real*8    start and the end of the integration
c    f(0:2*ns) :  real*8    function values at the nodes
c    integ     :  real*8    the integrated value
c
c                                                    1992.6  N.Takeuchi
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer ns
	real*8 xe,xs,f(0:2*ns),integ
	integer i,nn
	real*8 dx,sum,sumo,sume
c
	nn=2*ns
	sum = f(0)+f(nn)
	sumo = 0.d0
	sume = 0.d0
	do 100 i=1,ns-1
	  nn = 2 * i - 1
	  sumo = sumo + f(nn)
	  nn = 2 * i
	  sume = sume + f(nn)
  100	continue
	nn = 2 * ns - 1
	sumo = sumo + f(nn)
	dx = ( xe - xs ) / ( 2.d0 * ns )
	sum = ( sum + 4.d0 * sumo + 2.0 * sume )
     &	         * dx / 3.d0
	integ = sum
c
	return
	end
