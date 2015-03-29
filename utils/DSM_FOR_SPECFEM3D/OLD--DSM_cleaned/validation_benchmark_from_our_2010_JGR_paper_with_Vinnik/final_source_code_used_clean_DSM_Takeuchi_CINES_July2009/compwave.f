cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_excitation
     &        ( maxngrid_r,omega,
     &          ngrid_r,grid_r,l,source_r,source_mt,
     &          igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &          grid_qkappas,grid_qmus,
     &          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &          submatrix_I6,submatrix_I7,
     &          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &          idim_rs_sph,idim_rs_tor,
     &          whole_vector_sph,whole_vector_tor )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing a excitation vector for the given frequency.
c    required subroutines: error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,ngrid_r,l,igrid_rs
      integer idim_rs_sph,idim_rs_tor
      real*8 grid_r(*)
      real*8 source_r,source_mt(3,3)
      real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
      real*8 grid_qkappas,grid_qmus
      real*8 submatrix_I0(4,maxngrid_r)
      real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
      real*8 submatrix_I2(4,maxngrid_r)
      real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
      real*8 submatrix_I4(4,maxngrid_r)
      real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
      real*8 submatrix_I6(4,maxngrid_r)
      real*8 submatrix_I7(4,maxngrid_r)
      real*8 submatrix_I3k_mod(6,maxngrid_r)
      real*8 submatrix_I3m_mod(6,maxngrid_r)
      real*8 submatrix_I4_mod(6,maxngrid_r)
      complex*16 omega
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
c other variables
      real*8 b1,b2,lsq2,lsq
      complex*16 D1,D2_p,D2_m,D3_p,D3_m,unelastic_factor
      complex*16 factors_qkappa,factors_qmu
c constant
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      if ( l.le.0 ) call error_handling(51)
c **********************************************************************
c Initialing the whole_vector
c **********************************************************************
      call init_complex_array( 10*maxngrid_r,whole_vector_sph )
      call init_complex_array(  5*maxngrid_r,whole_vector_tor )
c **********************************************************************
c Computing the excitation vector
c **********************************************************************
      b1 = dsqrt( dble(2*l+1) / ( 16.d0 * pi ) )
      b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2) / ( 64.d0 * pi ) )
      factors_qkappa = unelastic_factor( dble(omega), grid_qkappas )
      factors_qmu = unelastic_factor( dble(omega), grid_qmus )
      lsq2 = dble(l) * dble(l+1)
      lsq = dsqrt( lsq2 )
c --- spheroidal excitations due to the traction diccontinuities
      whole_vector_sph(idim_rs_sph,0)
     &        = whole_vector_sph(idim_rs_sph,0)
     &          - b1 * 2.d0
     &            * ( source_mt(2,2) + source_mt(3,3)
     &                - 2.d0 * source_mt(1,1)
     &                  * ( grid_Fks * factors_qkappa
     &                      + grid_Fms * factors_qmu )
     &                  / ( grid_Cks * factors_qkappa
     &                      + grid_Cms * factors_qmu )
     &                ) / source_r
      whole_vector_sph(idim_rs_sph+1,0)
     &        = whole_vector_sph(idim_rs_sph+1,0)
     &          - b1 * lsq
     &            * ( - source_mt(2,2) - source_mt(3,3)
     &                + 2.d0 * source_mt(1,1)
     &                  * ( grid_Fks * factors_qkappa
     &                      + grid_Fms * factors_qmu )
     &                  / ( grid_Cks * factors_qkappa
     &                      + grid_Cms * factors_qmu )
     &                ) / source_r
      whole_vector_sph(idim_rs_sph+1,-2)
     &        = whole_vector_sph(idim_rs_sph+1,-2)
     &          + b2
     &            * dcmplx( - source_mt(2,2) + source_mt(3,3),
     &                      - 2.d0 * source_mt(2,3) )
     &              / source_r
      whole_vector_sph(idim_rs_sph+1,2)
     &        = whole_vector_sph(idim_rs_sph+1,2)
     &          + b2
     &            * dcmplx( - source_mt(2,2) + source_mt(3,3),
     &                        2.d0 * source_mt(2,3) )
     &              / source_r
c --- toroidal excitations due to the traction diccontinuities
      whole_vector_tor(idim_rs_tor,-2)
     &        = whole_vector_tor(idim_rs_tor,-2)
     &          - b2 * ( dcmplx( - 2.d0 * source_mt(2,3),
     &                             source_mt(2,2) - source_mt(3,3) ) )
     &            / source_r
      whole_vector_tor(idim_rs_tor,2)
     &        = whole_vector_tor(idim_rs_tor,2)
     &          - b2 * ( dcmplx( - 2.d0 * source_mt(2,3),
     &                           - source_mt(2,2) + source_mt(3,3) ) )
     &            / source_r
c --- excitations due to the displacement discontinuities
      D1 = b1 * 2.d0 * source_mt(1,1)
     &           / ( source_r * source_r
     &               * ( grid_Cks * factors_qkappa
     &                   + grid_Cms * factors_qmu )
     &              )
      D2_p = b1 * dcmplx( -source_mt(1,2), source_mt(1,3) )
     &           / ( source_r * source_r * grid_Ls * factors_qmu )
      D2_m = b1 * dcmplx( source_mt(1,2), source_mt(1,3) )
     &           / ( source_r * source_r * grid_Ls * factors_qmu )
      D3_p = b1 * dcmplx( source_mt(1,3), source_mt(1,2) )
     &           / ( source_r * source_r * grid_Ls * factors_qmu )
      D3_m = b1 * dcmplx( -source_mt(1,3), source_mt(1,2) )
     &           / ( source_r * source_r * grid_Ls * factors_qmu )
c ---- spheroidal, m=0
      whole_vector_sph(idim_rs_sph,0)
     &        = whole_vector_sph(idim_rs_sph,0)
     &            + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &                + ( submatrix_I1k(1,igrid_rs)
     &                    + 4.d0 * submatrix_I3k(1,igrid_rs)
     &                    + 4.d0 * submatrix_I5k(1,igrid_rs)
     &                  ) * factors_qkappa
     &                + ( submatrix_I1m(1,igrid_rs)
     &                    + 4.d0 * submatrix_I3m(1,igrid_rs)
     &                    + 4.d0 * submatrix_I5m(1,igrid_rs)
     &                    + lsq2 * submatrix_I6(1,igrid_rs)
     &                    - 4.d0 * submatrix_I7(1,igrid_rs)
     &                  ) * factors_qmu
     &              ) * D1
      whole_vector_sph(idim_rs_sph+1,0)
     &        = whole_vector_sph(idim_rs_sph+1,0)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(1,igrid_rs)
     &                      + 2.d0 * submatrix_I5k(1,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(1,igrid_rs)
     &                        - submatrix_I4_mod(1,igrid_rs)
     &                        + 2.d0 * submatrix_I5m(1,igrid_rs)
     &                        + submatrix_I6(1,igrid_rs)
     &                        - 2.d0 * submatrix_I7(1,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D1
      whole_vector_sph(idim_rs_sph+2,0)
     &        = whole_vector_sph(idim_rs_sph+2,0)
     &            + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &                + ( submatrix_I1k(3,igrid_rs)
     &                    + 2.d0 * submatrix_I3k(2,igrid_rs)
     &                    + 2.d0 * submatrix_I3k(3,igrid_rs)
     &                    + 4.d0 * submatrix_I5k(3,igrid_rs)
     &                  ) * factors_qkappa
     &                + ( submatrix_I1m(3,igrid_rs)
     &                    + 2.d0 * submatrix_I3m(2,igrid_rs)
     &                    + 2.d0 * submatrix_I3m(3,igrid_rs)
     &                    + 4.d0 * submatrix_I5m(3,igrid_rs)
     &                    + lsq2 * submatrix_I6(3,igrid_rs)
     &                    - 4.d0 * submatrix_I7(3,igrid_rs)
     &                  ) * factors_qmu
     &              ) * D1
      whole_vector_sph(idim_rs_sph+3,0)
     &        = whole_vector_sph(idim_rs_sph+3,0)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(3,igrid_rs)
     &                      + 2.d0 * submatrix_I5k(3,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(3,igrid_rs)
     &                        - submatrix_I4_mod(3,igrid_rs)
     &                        + 2.d0 * submatrix_I5m(3,igrid_rs)
     &                        + submatrix_I6(3,igrid_rs)
     &                        - 2.d0 * submatrix_I7(3,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D1
c ---- spheroidal, m=-1
      whole_vector_sph(idim_rs_sph,-1)
     &        = whole_vector_sph(idim_rs_sph,-1)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(1,igrid_rs)
     &                      + 2.d0 * submatrix_I5k(1,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(1,igrid_rs)
     &                        - submatrix_I4_mod(1,igrid_rs)
     &                        + 2.d0 * submatrix_I5m(1,igrid_rs)
     &                        + submatrix_I6(1,igrid_rs)
     &                        - 2.d0 * submatrix_I7(1,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D2_m
      whole_vector_sph(idim_rs_sph+1,-1)
     &        = whole_vector_sph(idim_rs_sph+1,-1)
     &            + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &                + ( lsq2 * submatrix_I5k(1,igrid_rs)
     &                  ) * factors_qkappa
     &                  + ( submatrix_I2(1,igrid_rs)
     &                      - 2.d0 * submatrix_I4(1,igrid_rs)
     &                      + lsq2 * submatrix_I5m(1,igrid_rs)
     &                      + submatrix_I6(1,igrid_rs)
     &                      - 2.d0 * submatrix_I7(1,igrid_rs)
     &                  ) * factors_qmu
     &              )  * D2_m
      whole_vector_sph(idim_rs_sph+2,-1)
     &        = whole_vector_sph(idim_rs_sph+2,-1)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(2,igrid_rs)
     &                      + 2.d0 * submatrix_I5k(3,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(2,igrid_rs)
     &                        - submatrix_I4_mod(2,igrid_rs)
     &                        + 2.d0 * submatrix_I5m(3,igrid_rs)
     &                        + submatrix_I6(3,igrid_rs)
     &                        - 2.d0 * submatrix_I7(3,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D2_m
      whole_vector_sph(idim_rs_sph+3,-1)
     &        = whole_vector_sph(idim_rs_sph+3,-1)
     &            + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &                + ( lsq2 * submatrix_I5k(3,igrid_rs)
     &                  ) * factors_qkappa
     &                  + ( submatrix_I2(3,igrid_rs)
     &                      - submatrix_I4(2,igrid_rs)
     &                      - submatrix_I4(3,igrid_rs)
     &                      + lsq2 * submatrix_I5m(3,igrid_rs)
     &                      + submatrix_I6(3,igrid_rs)
     &                      - 2.d0 * submatrix_I7(3,igrid_rs)
     &                  ) * factors_qmu
     &              )  * D2_m
      whole_vector_sph(idim_rs_sph+4,-1)
     &        = whole_vector_sph(idim_rs_sph+4,-1)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(5,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(5,igrid_rs)
     &                        - submatrix_I4_mod(5,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D2_m
c ---- spheroidal, m=+1
      whole_vector_sph(idim_rs_sph,1)
     &        = whole_vector_sph(idim_rs_sph,1)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(1,igrid_rs)
     &                      + 2.d0 * submatrix_I5k(1,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(1,igrid_rs)
     &                        - submatrix_I4_mod(1,igrid_rs)
     &                        + 2.d0 * submatrix_I5m(1,igrid_rs)
     &                        + submatrix_I6(1,igrid_rs)
     &                        - 2.d0 * submatrix_I7(1,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D2_p
      whole_vector_sph(idim_rs_sph+1,1)
     &        = whole_vector_sph(idim_rs_sph+1,1)
     &            + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &                + ( lsq2 * submatrix_I5k(1,igrid_rs)
     &                  ) * factors_qkappa
     &                  + ( submatrix_I2(1,igrid_rs)
     &                      - 2.d0 * submatrix_I4(1,igrid_rs)
     &                      + lsq2 * submatrix_I5m(1,igrid_rs)
     &                      + submatrix_I6(1,igrid_rs)
     &                      - 2.d0 * submatrix_I7(1,igrid_rs)
     &                  ) * factors_qmu
     &              )  * D2_p
      whole_vector_sph(idim_rs_sph+2,1)
     &        = whole_vector_sph(idim_rs_sph+2,1)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(2,igrid_rs)
     &                      + 2.d0 * submatrix_I5k(3,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(2,igrid_rs)
     &                        - submatrix_I4_mod(2,igrid_rs)
     &                        + 2.d0 * submatrix_I5m(3,igrid_rs)
     &                        + submatrix_I6(3,igrid_rs)
     &                        - 2.d0 * submatrix_I7(3,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D2_p
      whole_vector_sph(idim_rs_sph+3,1)
     &        = whole_vector_sph(idim_rs_sph+3,1)
     &            + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &                + ( lsq2 * submatrix_I5k(3,igrid_rs)
     &                  ) * factors_qkappa
     &                  + ( submatrix_I2(3,igrid_rs)
     &                      - submatrix_I4(2,igrid_rs)
     &                      - submatrix_I4(3,igrid_rs)
     &                      + lsq2 * submatrix_I5m(3,igrid_rs)
     &                      + submatrix_I6(3,igrid_rs)
     &                      - 2.d0 * submatrix_I7(3,igrid_rs)
     &                  ) * factors_qmu
     &              )  * D2_p
      whole_vector_sph(idim_rs_sph+4,1)
     &        = whole_vector_sph(idim_rs_sph+4,1)
     &            + ( - lsq * (
     &                    ( submatrix_I3k_mod(5,igrid_rs)
     &                    ) * factors_qkappa
     &                    + ( submatrix_I3m_mod(5,igrid_rs)
     &                        - submatrix_I4_mod(5,igrid_rs)
     &                    ) * factors_qmu
     &              ) ) * D2_p
c ---- toroidal, m=-1
      whole_vector_tor(idim_rs_tor,-1)
     &        = whole_vector_tor(idim_rs_tor,-1)
     &            + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &                + ( submatrix_I2(1,igrid_rs)
     &                    - 2.d0 * submatrix_I4(1,igrid_rs)
     &                    + submatrix_I6(1,igrid_rs)
     &                    - ( lsq2 - 2.d0 ) * submatrix_I7(1,igrid_rs)
     &                   ) * factors_qmu
     &               ) * D3_m
      whole_vector_tor(idim_rs_tor+1,-1)
     &        = whole_vector_tor(idim_rs_tor+1,-1)
     &            + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &                + ( submatrix_I2(3,igrid_rs)
     &                    - submatrix_I4(2,igrid_rs)
     &                    - submatrix_I4(3,igrid_rs)
     &                    + submatrix_I6(3,igrid_rs)
     &                    - ( lsq2 - 2.d0 ) * submatrix_I7(3,igrid_rs)
     &                   ) * factors_qmu
     &               ) * D3_m
c ---- toroidal, m=+1
      whole_vector_tor(idim_rs_tor,1)
     &        = whole_vector_tor(idim_rs_tor,1)
     &            + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &                + ( submatrix_I2(1,igrid_rs)
     &                    - 2.d0 * submatrix_I4(1,igrid_rs)
     &                    + submatrix_I6(1,igrid_rs)
     &                    - ( lsq2 - 2.d0 ) * submatrix_I7(1,igrid_rs)
     &                   ) * factors_qmu
     &               ) * D3_p
      whole_vector_tor(idim_rs_tor+1,1)
     &        = whole_vector_tor(idim_rs_tor+1,1)
     &            + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &                + ( submatrix_I2(3,igrid_rs)
     &                    - submatrix_I4(2,igrid_rs)
     &                    - submatrix_I4(3,igrid_rs)
     &                    + submatrix_I6(3,igrid_rs)
     &                    - ( lsq2 - 2.d0 ) * submatrix_I7(3,igrid_rs)
     &                   ) * factors_qmu
     &               ) * D3_p
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_excitation0
     &        ( maxngrid_r,omega,
     &          ngrid_r,grid_r,l,source_r,source_mt,
     &          igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &          grid_qkappas,grid_qmus,
     &          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &          submatrix_I6,submatrix_I7,
     &          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &          idim_rs_sph,idim_rs_tor,
     &          whole_vector_sph,whole_vector_tor )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing a excitation vector for the given frequency.
c    required subroutines: error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,ngrid_r,l,igrid_rs
      integer idim_rs_sph,idim_rs_tor
      real*8 grid_r(*)
      real*8 source_r,source_mt(3,3)
      real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
      real*8 grid_qkappas,grid_qmus
      real*8 submatrix_I0(4,maxngrid_r)
      real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
      real*8 submatrix_I2(4,maxngrid_r)
      real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
      real*8 submatrix_I4(4,maxngrid_r)
      real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
      real*8 submatrix_I6(4,maxngrid_r)
      real*8 submatrix_I7(4,maxngrid_r)
      real*8 submatrix_I3k_mod(6,maxngrid_r)
      real*8 submatrix_I3m_mod(6,maxngrid_r)
      real*8 submatrix_I4_mod(6,maxngrid_r)
      complex*16 omega
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
c other variables
      real*8 b1,b2,lsq2
      complex*16 D1,D2_p,D2_m,D3_p,D3_m,unelastic_factor
      complex*16 factors_qkappa,factors_qmu
c constant
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      if ( l.ne.0 ) call error_handling(52)
c **********************************************************************
c Initialing the whole_vector
c **********************************************************************
      call init_complex_array( 10*maxngrid_r,whole_vector_sph )
      call init_complex_array(  5*maxngrid_r,whole_vector_tor )
c **********************************************************************
c Computing the excitation vector
c **********************************************************************
      b1 = dsqrt( dble(2*l+1) / ( 16.d0 * pi ) )
c     b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2) / ( 64.d0 * pi ) )
      b2 = 0.d0
c --- excitations due to the traction diccontinuities
      factors_qkappa = unelastic_factor( dble(omega), grid_qkappas )
      factors_qmu = unelastic_factor( dble(omega), grid_qmus )
      lsq2 = dble(l) * dble(l+1)
      whole_vector_sph(idim_rs_sph,0)
     &        = whole_vector_sph(idim_rs_sph,0)
     &          - b1 * 2.d0
     &            * ( source_mt(2,2) + source_mt(3,3)
     &                - 2.d0 * source_mt(1,1)
     &                  * ( grid_Fks * factors_qkappa
     &                      + grid_Fms * factors_qmu
     &                    )
     &                  / ( grid_Cks * factors_qkappa
     &                      + grid_Cms * factors_qmu
     &                    )
     &                ) / source_r
c --- excitations due to the displacement discontinuities
      D1 = b1 * 2.d0 * source_mt(1,1)
     &           / ( source_r * source_r
     &               * ( grid_Cks * factors_qkappa
     &                   + grid_Cms * factors_qmu ) )
      D2_p = dcmplx( 0.d0 )
      D2_m = dcmplx( 0.d0 )
      D3_p = dcmplx( 0.d0 )
      D3_m = dcmplx( 0.d0 )
c ---- spheroidal, m=0
      whole_vector_sph(idim_rs_sph,0)
     &        = whole_vector_sph(idim_rs_sph,0)
     &            + ( - omega * omega * submatrix_I0(1,igrid_rs)
     &                + ( submatrix_I1k(1,igrid_rs)
     &                    + 4.d0 * submatrix_I3k(1,igrid_rs)
     &                    + 4.d0 * submatrix_I5k(1,igrid_rs)
     &                  ) * factors_qkappa
     &                + ( submatrix_I1m(1,igrid_rs)
     &                    + 4.d0 * submatrix_I3m(1,igrid_rs)
     &                    + 4.d0 * submatrix_I5m(1,igrid_rs)
     &                    + lsq2 * submatrix_I6(1,igrid_rs)
     &                    - 4.d0 * submatrix_I7(1,igrid_rs)
     &                  ) * factors_qmu
     &              ) * D1
      whole_vector_sph(idim_rs_sph+1,0)
     &        = whole_vector_sph(idim_rs_sph+1,0)
     &            + ( - omega * omega * submatrix_I0(3,igrid_rs)
     &                + ( submatrix_I1k(3,igrid_rs)
     &                    + 2.d0 * submatrix_I3k(2,igrid_rs)
     &                    + 2.d0 * submatrix_I3k(3,igrid_rs)
     &                    + 4.d0 * submatrix_I5k(3,igrid_rs)
     &                  ) * factors_qkappa
     &                + ( submatrix_I1m(3,igrid_rs)
     &                    + 2.d0 * submatrix_I3m(2,igrid_rs)
     &                    + 2.d0 * submatrix_I3m(3,igrid_rs)
     &                    + 4.d0 * submatrix_I5m(3,igrid_rs)
     &                    + lsq2 * submatrix_I6(3,igrid_rs)
     &                    - 4.d0 * submatrix_I7(3,igrid_rs)
     &                  ) * factors_qmu
     &              ) * D1
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_wavefield
     &        ( maxngrid_r,omega,
     &          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &          submatrix_I6,submatrix_I7,
     &          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &          ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &          idim1_sph0,idim2_sph,idim1_tor0,idim2_tor,
     &          idim0,init_npos_sph,init_npos_tor,
     &          idim_rs_sph,idim_rs_tor,
     &          idim_station_sph,idim_station_tor,
     &          whole_matrix_sph,whole_matrix_tor,
     &          whole_matrix_dr_sph,whole_matrix_dr_tor,
     &          whole_vector_sph,whole_vector_tor,work_vector )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing wavefield for the given frequency.
c    required subroutines: init_complex_array
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,ngrid_r,l
      integer idim1_sph0,idim2_sph,idim1_tor0,idim2_tor
      integer idim0,init_npos_sph,init_npos_tor
      integer idim_rs_sph,idim_rs_tor
      integer idim_station_sph,idim_station_tor
      real*8 submatrix_I0(4,maxngrid_r)
      real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
      real*8 submatrix_I2(4,maxngrid_r)
      real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
      real*8 submatrix_I4(4,maxngrid_r)
      real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
      real*8 submatrix_I6(4,maxngrid_r)
      real*8 submatrix_I7(4,maxngrid_r)
      real*8 submatrix_I3k_mod(6,maxngrid_r)
      real*8 submatrix_I3m_mod(6,maxngrid_r)
      real*8 submatrix_I4_mod(6,maxngrid_r)
      real*8 grid_r(*),grid_mu(2,*)
      real*8 grid_qkappa(*),grid_qmu(*)
      complex*16 omega
      complex*16 whole_matrix_sph(4,*),whole_matrix_tor(2,*)
      complex*16 whole_matrix_dr_sph(*),whole_matrix_dr_tor(*)
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
      complex*16 work_vector(*)
c other variables
      integer idim1_sph,idim1_tor,ir,npos,m,itype_medium
c itype_medium=1: solid, itype_medium=0: liquid
      integer ndim_whole_matrix_sph,ndim_whole_matrix_tor,ier
      integer init_grid,end_grid,ns,nq
      real*8 lsq,lsq2,eps
      complex*16 unelastic_factor,factor_qkappa(2),factor_qmu(2)
      data eps / -1.d0 /
c
      if ( l.le.0 ) call error_handling(53)
c **********************************************************************
c Initialing the whole_matrix
c **********************************************************************
      call init_complex_array( 8*maxngrid_r,whole_matrix_sph )
      call init_complex_array( 2*maxngrid_r,whole_matrix_tor )
      call init_complex_array( 2*maxngrid_r,whole_matrix_dr_sph )
      call init_complex_array(   maxngrid_r,whole_matrix_dr_tor )
      idim1_sph = max0( idim0,idim1_sph0 )
      idim1_tor = max0( idim0,idim1_tor0 )
c **********************************************************************
c **********************************************************************
c Spheroidal Component
c **********************************************************************
c **********************************************************************
c **********************************************************************
c constructing the whole_matrix
c **********************************************************************
      lsq2 = dble(l) * dble(l+1)
      lsq  = dsqrt( lsq2 )
c
      factor_qkappa(2)
     &        = unelastic_factor( dble(omega),grid_qkappa(idim1_sph) )
      factor_qmu(2)
     &        = unelastic_factor( dble(omega),grid_qmu(idim1_sph) )
      if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).eq.0.d0 ) then
        npos = init_npos_sph
        itype_medium = 0
        whole_matrix_sph(4,npos)
     &        = omega * omega / factor_qkappa(2)
     &            * dcmplx( submatrix_I0(1,idim1_sph)
     &                    )
     &            - dcmplx( lsq2 * submatrix_I1k(1,idim1_sph)
     &                        + submatrix_I2(1,idim1_sph)
     &                    )
      else
        npos = init_npos_sph
        itype_medium = 1
        whole_matrix_sph(4,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(1,idim1_sph)
     &                    )
     &            - factor_qkappa(2)
     &              * dcmplx( submatrix_I1k(1,idim1_sph)
     &                        + 4.d0 * submatrix_I3k(1,idim1_sph)
     &                        + 4.d0 * submatrix_I5k(1,idim1_sph)
     &                      )
     &            - factor_qmu(2)
     &              * dcmplx( submatrix_I1m(1,idim1_sph)
     &                        + 4.d0 * submatrix_I3m(1,idim1_sph)
     &                        + 4.d0 * submatrix_I5m(1,idim1_sph)
     &                        + lsq2 * submatrix_I6(1,idim1_sph)
     &                        - 4.d0 * submatrix_I7(1,idim1_sph)
     &                      )
        whole_matrix_sph(3,npos+1)
     &        = dcmplx( lsq )
     &            * (
     &              factor_qkappa(2)
     &              * dcmplx( submatrix_I3k_mod(1,idim1_sph)
     &                        + 2.d0 * submatrix_I5k(1,idim1_sph)
     &                      )
     &            + factor_qmu(2)
     &              * dcmplx( submatrix_I3m_mod(1,idim1_sph)
     &                        - submatrix_I4_mod(1,idim1_sph)
     &                        + 2.d0 * submatrix_I5m(1,idim1_sph)
     &                        + submatrix_I6(1,idim1_sph)
     &                        - 2.d0 * submatrix_I7(1,idim1_sph)
     &                      )
     &             )
        whole_matrix_sph(4,npos+1)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(1,idim1_sph)
     &                    )
     &            - factor_qkappa(2)
     &              * dcmplx( lsq2 * submatrix_I5k(1,idim1_sph)
     &                      )
     &            - factor_qmu(2)
     &              * dcmplx( submatrix_I2(1,idim1_sph)
     &                        - 2.d0 * submatrix_I4(1,idim1_sph)
     &                        + lsq2 * submatrix_I5m(1,idim1_sph)
     &                        + submatrix_I6(1,idim1_sph)
     &                        - 2.d0 * submatrix_I7(1,idim1_sph)
     &                      )
      endif
      do 120 ir=idim1_sph+1,idim2_sph-1
        factor_qkappa(1)
     &          = unelastic_factor( dble(omega),grid_qkappa(ir-1) )
        factor_qmu(1)
     &          = unelastic_factor( dble(omega),grid_qmu(ir-1) )
        factor_qkappa(2)
     &          = unelastic_factor( dble(omega),grid_qkappa(ir) )
        factor_qmu(2)
     &          = unelastic_factor( dble(omega),grid_qmu(ir) )
        if ( grid_mu(1,ir)*grid_mu(2,ir).eq.0.d0 ) then
          if ( itype_medium.eq.1 ) then
            npos = npos + 2
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &                * dcmplx( submatrix_I1k(2,ir-1)
     &                      + 2.d0 * submatrix_I3k(2,ir-1)
     &                      + 2.d0 * submatrix_I3k(3,ir-1)
     &                      + 4.d0 * submatrix_I5k(2,ir-1)
     &                        )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I1m(2,ir-1)
     &                      + 2.d0 * submatrix_I3m(2,ir-1)
     &                      + 2.d0 * submatrix_I3m(3,ir-1)
     &                      + 4.d0 * submatrix_I5m(2,ir-1)
     &                      + lsq2 * submatrix_I6(2,ir-1)
     &                      - 4.d0 * submatrix_I7(2,ir-1)
     &                    )
            whole_matrix_sph(3,npos)
     &            = whole_matrix_sph(3,npos)
     &              + dcmplx( lsq )
     &              * (
     &              factor_qkappa(1)
     &              * dcmplx( submatrix_I3k_mod(2,ir-1)
     &                        + submatrix_I3k_mod(5,ir-1)
     &                        + 2.d0 * submatrix_I5k(2,ir-1)
     &                      )
     &              + factor_qmu(1)
     &                * dcmplx( submatrix_I3m_mod(2,ir-1)
     &                          + submatrix_I3m_mod(5,ir-1)
     &                          - submatrix_I4_mod(2,ir-1)
     &                          - submatrix_I4_mod(5,ir-1)
     &                          + 2.d0 * submatrix_I5m(2,ir-1)
     &                          + submatrix_I6(2,ir-1)
     &                          - 2.d0 * submatrix_I7(2,ir-1)
     &                        )
     &               )
            whole_matrix_sph(4,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &                * dcmplx( submatrix_I1k(4,ir-1)
     &                          + 4.d0 * submatrix_I3k(4,ir-1)
     &                          + 4.d0 * submatrix_I5k(4,ir-1)
     &                        )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I1m(4,ir-1)
     &                          + 4.d0 * submatrix_I3m(4,ir-1)
     &                          + 4.d0 * submatrix_I5m(4,ir-1)
     &                          + lsq2 * submatrix_I6(4,ir-1)
     &                          - 4.d0 * submatrix_I7(4,ir-1)
     &                         )
            whole_matrix_sph(1,npos+1)
     &            = dcmplx( lsq )
     &              * (
     &                factor_qkappa(1)
     &                * dcmplx( submatrix_I3k_mod(3,ir-1)
     &                          + 2.d0 * submatrix_I5k(2,ir-1)
     &                        )
     &                + factor_qmu(1)
     &                  * dcmplx( submatrix_I3m_mod(3,ir-1)
     &                            - submatrix_I4_mod(3,ir-1)
     &                            + 2.d0 * submatrix_I5m(2,ir-1)
     &                            + submatrix_I6(2,ir-1)
     &                            - 2.d0 * submatrix_I7(2,ir-1)
     &                          )
     &                 )
            whole_matrix_sph(2,npos+1)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &                * dcmplx( lsq2 * submatrix_I5k(2,ir-1)
     &                        )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I2(2,ir-1)
     &                          - submatrix_I4(2,ir-1)
     &                          - submatrix_I4(3,ir-1)
     &                          + lsq2 * submatrix_I5m(2,ir-1)
     &                          + submatrix_I6(2,ir-1)
     &                          - 2.d0 * submatrix_I7(2,ir-1)
     &                        )
            whole_matrix_sph(3,npos+1)
     &            = dcmplx( lsq )
     &              * (
     &                factor_qkappa(1)
     &                * dcmplx( submatrix_I3k_mod(4,ir-1)
     &                          + submatrix_I3k_mod(6,ir-1)
     &                          + 2.d0 * submatrix_I5k(4,ir-1)
     &                        )
     &                + factor_qmu(1)
     &                  * dcmplx( submatrix_I3m_mod(4,ir-1)
     &                            + submatrix_I3m_mod(6,ir-1)
     &                            - submatrix_I4_mod(4,ir-1)
     &                            - submatrix_I4_mod(6,ir-1)
     &                            + 2.d0 * submatrix_I5m(4,ir-1)
     &                            + submatrix_I6(4,ir-1)
     &                            - 2.d0 * submatrix_I7(4,ir-1)
     &                          )
     &                )
            whole_matrix_sph(4,npos+1)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                      )
     &                - factor_qkappa(1)
     &                  * dcmplx( lsq2 * submatrix_I5k(4,ir-1)
     &                          )
     &                - factor_qmu(1)
     &                  * dcmplx( submatrix_I2(4,ir-1)
     &                            - 2.d0 * submatrix_I4(4,ir-1)
     &                            + lsq2 * submatrix_I5m(4,ir-1)
     &                            + submatrix_I6(4,ir-1)
     &                            - 2.d0 * submatrix_I7(4,ir-1)
     &                          )
            npos = npos + 2
            itype_medium = 0
            whole_matrix_sph(2,npos)
     &            = omega * grid_r(ir) * grid_r(ir)
            whole_matrix_sph(4,npos)
     &              = omega * omega / factor_qkappa(2)
     &              * dcmplx( submatrix_I0(1,ir)
     &                      )
     &              - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &                        + submatrix_I2(1,ir)
     &                      )
          else
            npos = npos + 1
            itype_medium = 0
            whole_matrix_sph(3,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1) ) / factor_qkappa(1)
     &              - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &                        + submatrix_I2(2,ir-1)
     &                      )
            whole_matrix_sph(4,npos)
     &             = omega * omega
     &             * ( dcmplx( submatrix_I0(4,ir-1) ) / factor_qkappa(1)
     &                 + dcmplx( submatrix_I0(1,ir) ) / factor_qkappa(2)
     &                )
     &              - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &                        + submatrix_I2(4,ir-1)
     &                      )
     &              - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &                        + submatrix_I2(1,ir)
     &                      )
          endif
        else
          if ( itype_medium.eq.0 ) then
            npos = npos + 1
            whole_matrix_sph(3,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(2,ir-1) )
     &              - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &                        + submatrix_I2(2,ir-1)
     &                      )
            whole_matrix_sph(4,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                      )
     &              - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &                        + submatrix_I2(4,ir-1)
     &                      )
            npos = npos + 1
            itype_medium = 1
            whole_matrix_sph(3,npos)
     &            = - omega * grid_r(ir) * grid_r(ir)
            whole_matrix_sph(4,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(1,ir)
     &                      )
     &              - factor_qkappa(2)
     &              * dcmplx( submatrix_I1k(1,ir)
     &                        + 4.d0 * submatrix_I3k(1,ir)
     &                        + 4.d0 * submatrix_I5k(1,ir)
     &                      )
     &              - factor_qmu(2)
     &              * dcmplx( submatrix_I1m(1,ir)
     &                        + 4.d0 * submatrix_I3m(1,ir)
     &                        + 4.d0 * submatrix_I5m(1,ir)
     &                        + lsq2 * submatrix_I6(1,ir)
     &                        - 4.d0 * submatrix_I7(1,ir)
     &                      )
            whole_matrix_sph(3,npos+1)
     &            = dcmplx( lsq )
     &            * (
     &                factor_qkappa(2)
     &              * dcmplx( submatrix_I3k_mod(1,ir)
     &                        + 2.d0 * submatrix_I5k(1,ir)
     &                      )
     &              + factor_qmu(2)
     &              * dcmplx( submatrix_I3m_mod(1,ir)
     &                        - submatrix_I4_mod(1,ir)
     &                        + 2.d0 * submatrix_I5m(1,ir)
     &                        + submatrix_I6(1,ir)
     &                        - 2.d0 * submatrix_I7(1,ir)
     &                      )
     &             )
            whole_matrix_sph(4,npos+1)
     &            = omega * omega
     &            * dcmplx( submatrix_I0(1,ir)
     &                    )
     &            - factor_qkappa(2)
     &              * dcmplx( lsq2 * submatrix_I5k(1,ir)
     &                      )
     &            - factor_qmu(2)
     &              * dcmplx( submatrix_I2(1,ir)
     &                        - 2.d0 * submatrix_I4(1,ir)
     &                        + lsq2 * submatrix_I5m(1,ir)
     &                        + submatrix_I6(1,ir)
     &                        - 2.d0 * submatrix_I7(1,ir)
     &                      )
          else
            npos = npos + 2
            itype_medium = 1
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &              * dcmplx( submatrix_I1k(2,ir-1)
     &                        + 2.d0 * submatrix_I3k(2,ir-1)
     &                        + 2.d0 * submatrix_I3k(3,ir-1)
     &                        + 4.d0 * submatrix_I5k(2,ir-1)
     &                      )
     &              - factor_qmu(1)
     &              * dcmplx( submatrix_I1m(2,ir-1)
     &                        + 2.d0 * submatrix_I3m(2,ir-1)
     &                        + 2.d0 * submatrix_I3m(3,ir-1)
     &                        + 4.d0 * submatrix_I5m(2,ir-1)
     &                        + lsq2 * submatrix_I6(2,ir-1)
     &                        - 4.d0 * submatrix_I7(2,ir-1)
     &                      )
            whole_matrix_sph(3,npos)
     &            = whole_matrix_sph(3,npos)
     &            + dcmplx( lsq )
     &            * (
     &              factor_qkappa(1)
     &              * dcmplx( submatrix_I3k_mod(2,ir-1)
     &                        + 2.d0 * submatrix_I5k(2,ir-1)
     &                      )
     &              + factor_qmu(1)
     &              * dcmplx( submatrix_I3m_mod(2,ir-1)
     &                        - submatrix_I4_mod(2,ir-1)
     &                        + 2.d0 * submatrix_I5m(2,ir-1)
     &                        + submatrix_I6(2,ir-1)
     &                        - 2.d0 * submatrix_I7(2,ir-1)
     &                      )
     &             )
            whole_matrix_sph(4,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                        + submatrix_I0(1,ir)
     &                      )
     &              - factor_qkappa(1)
     &              * dcmplx( submatrix_I1k(4,ir-1)
     &                        + 4.d0 * submatrix_I3k(4,ir-1)
     &                        + 4.d0 * submatrix_I5k(4,ir-1)
     &                      )
     &              - factor_qmu(1)
     &              * dcmplx( submatrix_I1m(4,ir-1)
     &                        + 4.d0 * submatrix_I3m(4,ir-1)
     &                        + 4.d0 * submatrix_I5m(4,ir-1)
     &                        + lsq2 * submatrix_I6(4,ir-1)
     &                        - 4.d0 * submatrix_I7(4,ir-1)
     &                      )
     &              - factor_qkappa(2)
     &              * dcmplx( submatrix_I1k(1,ir)
     &                        + 4.d0 * submatrix_I3k(1,ir)
     &                        + 4.d0 * submatrix_I5k(1,ir)
     &                      )
     &              - factor_qmu(2)
     &              * dcmplx( submatrix_I1m(1,ir)
     &                        + 4.d0 * submatrix_I3m(1,ir)
     &                        + 4.d0 * submatrix_I5m(1,ir)
     &                        + lsq2 * submatrix_I6(1,ir)
     &                        - 4.d0 * submatrix_I7(1,ir)
     &                      )
            whole_matrix_sph(1,npos+1)
     &            = dcmplx( lsq )
     &            * (
     &              factor_qkappa(1)
     &              * dcmplx( submatrix_I3k_mod(3,ir-1)
     &                        + 2.d0 * submatrix_I5k(2,ir-1)
     &                      )
     &              + factor_qmu(1)
     &              * dcmplx( submatrix_I3m_mod(3,ir-1)
     &                        - submatrix_I4_mod(3,ir-1)
     &                        + 2.d0 * submatrix_I5m(2,ir-1)
     &                        + submatrix_I6(2,ir-1)
     &                        - 2.d0 * submatrix_I7(2,ir-1)
     &                      )
     &               )
            whole_matrix_sph(2,npos+1)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                    )
     &              - factor_qkappa(1)
     &              * dcmplx( lsq2 * submatrix_I5k(2,ir-1)
     &                      )
     &              - factor_qmu(1)
     &              * dcmplx( submatrix_I2(2,ir-1)
     &                        - submatrix_I4(2,ir-1)
     &                        - submatrix_I4(3,ir-1)
     &                        + lsq2 * submatrix_I5m(2,ir-1)
     &                        + submatrix_I6(2,ir-1)
     &                        - 2.d0 * submatrix_I7(2,ir-1)
     &                      )
            whole_matrix_sph(3,npos+1)
     &            = dcmplx( lsq )
     &              * (
     &              factor_qkappa(1)
     &              * dcmplx( submatrix_I3k_mod(4,ir-1)
     &                        + 2.d0 * submatrix_I5k(4,ir-1)
     &                      )
     &              + factor_qmu(1)
     &              * dcmplx( submatrix_I3m_mod(4,ir-1)
     &                        - submatrix_I4_mod(4,ir-1)
     &                        + 2.d0 * submatrix_I5m(4,ir-1)
     &                        + submatrix_I6(4,ir-1)
     &                        - 2.d0 * submatrix_I7(4,ir-1)
     &                      )
     &              + factor_qkappa(2)
     &              * dcmplx( submatrix_I3k_mod(1,ir)
     &                        + 2.d0 * submatrix_I5k(1,ir)
     &                      )
     &              + factor_qmu(2)
     &              * dcmplx( submatrix_I3m_mod(1,ir)
     &                        - submatrix_I4_mod(1,ir)
     &                        + 2.d0 * submatrix_I5m(1,ir)
     &                        + submatrix_I6(1,ir)
     &                        - 2.d0 * submatrix_I7(1,ir)
     &                      )
     &               )
            whole_matrix_sph(4,npos+1)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                      + submatrix_I0(1,ir)
     &                    )
     &              - factor_qkappa(1)
     &              * dcmplx( lsq2 * submatrix_I5k(4,ir-1)
     &                      )
     &              - factor_qmu(1)
     &              * dcmplx( submatrix_I2(4,ir-1)
     &                        - 2.d0 * submatrix_I4(4,ir-1)
     &                        + lsq2 * submatrix_I5m(4,ir-1)
     &                        + submatrix_I6(4,ir-1)
     &                        - 2.d0 * submatrix_I7(4,ir-1)
     &                      )
     &              - factor_qkappa(2)
     &              * dcmplx( lsq2 * submatrix_I5k(1,ir)
     &                      )
     &              - factor_qmu(2)
     &              * dcmplx( submatrix_I2(1,ir)
     &                        - 2.d0 * submatrix_I4(1,ir)
     &                        + lsq2 * submatrix_I5m(1,ir)
     &                        + submatrix_I6(1,ir)
     &                        - 2.d0 * submatrix_I7(1,ir)
     &                      )
            whole_matrix_sph(1,npos+2)
     &            = dcmplx( lsq )
     &            * (
     &              factor_qkappa(1)
     &              * dcmplx( submatrix_I3k_mod(5,ir-1)
     &                      )
     &              + factor_qmu(1)
     &              * dcmplx( submatrix_I3m_mod(5,ir-1)
     &                        - submatrix_I4_mod(5,ir-1)
     &                      )
     &             )
            whole_matrix_sph(3,npos+2)
     &            = dcmplx( lsq )
     &            * (
     &              factor_qkappa(1)
     &              * dcmplx( submatrix_I3k_mod(6,ir-1)
     &                      )
     &              + factor_qmu(1)
     &              * dcmplx( submatrix_I3m_mod(6,ir-1)
     &                        - submatrix_I4_mod(6,ir-1)
     &                      )
     &             )
          endif
        endif
  120 continue
      factor_qkappa(1)
     &        = unelastic_factor( dble(omega),grid_qkappa(idim2_sph-1) )
      factor_qmu(1)
     &        = unelastic_factor( dble(omega),grid_qmu(idim2_sph-1) )
      if ( itype_medium.eq.0 ) then
        npos = npos + 1
        whole_matrix_sph(3,npos)
     &          = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(2,idim2_sph-1)
     &                      )
     &            - dcmplx( lsq2 * submatrix_I1k(2,idim2_sph-1)
     &                        + submatrix_I2(2,idim2_sph-1)
     &                    )
          whole_matrix_sph(4,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(4,idim2_sph-1)
     &                      )
     &            - dcmplx( lsq2 * submatrix_I1k(4,idim2_sph-1)
     &                        + submatrix_I2(4,idim2_sph-1)
     &                      )
        ndim_whole_matrix_sph = npos
      else
        npos = npos + 2
        whole_matrix_sph(2,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(2,idim2_sph-1)
     &                    )
     &          - factor_qkappa(1)
     &            * dcmplx( submatrix_I1k(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3k(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3k(3,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5k(2,idim2_sph-1)
     &                    )
     &          - factor_qmu(1)
     &            * dcmplx( submatrix_I1m(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3m(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3m(3,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5m(2,idim2_sph-1)
     &                      + lsq2 * submatrix_I6(2,idim2_sph-1)
     &                      - 4.d0 * submatrix_I7(2,idim2_sph-1)
     &                    )
        whole_matrix_sph(3,npos)
     &        = whole_matrix_sph(3,npos)
     &          + dcmplx( lsq )
     &          * (
     &            factor_qkappa(1)
     &            * dcmplx( submatrix_I3k_mod(2,idim2_sph-1)
     &                      + submatrix_I3k_mod(5,idim2_sph-1)
     &                      + 2.d0 * submatrix_I5k(2,idim2_sph-1)
     &                    )
     &          + factor_qmu(1)
     &            * dcmplx( submatrix_I3m_mod(2,idim2_sph-1)
     &                      + submatrix_I3m_mod(5,idim2_sph-1)
     &                      - submatrix_I4_mod(2,idim2_sph-1)
     &                      - submatrix_I4_mod(5,idim2_sph-1)
     &                      + 2.d0 * submatrix_I5m(2,idim2_sph-1)
     &                      + submatrix_I6(2,idim2_sph-1)
     &                      - 2.d0 * submatrix_I7(2,idim2_sph-1)
     &                    )
     &           )
        whole_matrix_sph(4,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(4,idim2_sph-1)
     &                    )
     &          - factor_qkappa(1)
     &            * dcmplx( submatrix_I1k(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I3k(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5k(4,idim2_sph-1)
     &                    )
     &          - factor_qmu(1)
     &            * dcmplx( submatrix_I1m(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I3m(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5m(4,idim2_sph-1)
     &                      + lsq2 * submatrix_I6(4,idim2_sph-1)
     &                      - 4.d0 * submatrix_I7(4,idim2_sph-1)
     &                    )
        whole_matrix_sph(1,npos+1)
     &        = dcmplx( lsq )
     &          * (
     &            factor_qkappa(1)
     &            * dcmplx( submatrix_I3k_mod(3,idim2_sph-1)
     &                      + 2.d0 * submatrix_I5k(2,idim2_sph-1)
     &                    )
     &          + factor_qmu(1)
     &            * dcmplx( submatrix_I3m_mod(3,idim2_sph-1)
     &                      - submatrix_I4_mod(3,idim2_sph-1)
     &                      + 2.d0 * submatrix_I5m(2,idim2_sph-1)
     &                      + submatrix_I6(2,idim2_sph-1)
     &                      - 2.d0 * submatrix_I7(2,idim2_sph-1)
     &                    )
     &           )
        whole_matrix_sph(2,npos+1)
     &        = omega * omega
     &          * dcmplx( submatrix_I0(2,idim2_sph-1)
     &                  )
     &          - factor_qkappa(1)
     &            * dcmplx( lsq2 * submatrix_I5k(2,idim2_sph-1)
     &                    )
     &          - factor_qmu(1)
     &            * dcmplx( submatrix_I2(2,idim2_sph-1)
     &                      - submatrix_I4(2,idim2_sph-1)
     &                      - submatrix_I4(3,idim2_sph-1)
     &                      + lsq2 * submatrix_I5m(2,idim2_sph-1)
     &                      + submatrix_I6(2,idim2_sph-1)
     &                      - 2.d0 * submatrix_I7(2,idim2_sph-1)
     &                    )
        whole_matrix_sph(3,npos+1)
     &        = dcmplx( lsq )
     &          * (
     &            factor_qkappa(1)
     &            * dcmplx( submatrix_I3k_mod(4,idim2_sph-1)
     &                      + submatrix_I3k_mod(6,idim2_sph-1)
     &                      + 2.d0 * submatrix_I5k(4,idim2_sph-1)
     &                    )
     &          + factor_qmu(1)
     &            * dcmplx( submatrix_I3m_mod(4,idim2_sph-1)
     &                      + submatrix_I3m_mod(6,idim2_sph-1)
     &                      - submatrix_I4_mod(4,idim2_sph-1)
     &                      - submatrix_I4_mod(6,idim2_sph-1)
     &                      + 2.d0 * submatrix_I5m(4,idim2_sph-1)
     &                      + submatrix_I6(4,idim2_sph-1)
     &                      - 2.d0 * submatrix_I7(4,idim2_sph-1)
     &                    )
     &           )
        whole_matrix_sph(4,npos+1)
     &        = omega * omega
     &          * dcmplx( submatrix_I0(4,idim2_sph-1)
     &                  )
     &          - factor_qkappa(1)
     &            * dcmplx( lsq2 * submatrix_I5k(4,idim2_sph-1)
     &                    )
     &          - factor_qmu(1)
     &            * dcmplx( submatrix_I2(4,idim2_sph-1)
     &                      - 2.d0 * submatrix_I4(4,idim2_sph-1)
     &                      + lsq2 * submatrix_I5m(4,idim2_sph-1)
     &                      + submatrix_I6(4,idim2_sph-1)
     &                      - 2.d0 * submatrix_I7(4,idim2_sph-1)
     &                    )
        ndim_whole_matrix_sph = npos+1
      endif
c **********************************************************************
c computing the wavefield
c **********************************************************************
c imposing fixed boundary conditions at r=0 and free surface boundary
c conditions at the surface
      if ( ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).ne.0.d0 )
     &           .and.( grid_r(idim1_sph).eq.0.d0 ) ) then
        init_grid = max0(init_npos_sph,3)
      else
        init_grid = init_npos_sph
      endif
      if ( grid_mu(1,idim2_sph-1)*grid_mu(2,idim2_sph-1).eq.0.d0 )
     &        then
        end_grid = ndim_whole_matrix_sph - 1
      else
        end_grid = ndim_whole_matrix_sph
      endif
      m = max0(-l,-2)
      ns = idim_rs_sph - init_grid + 1
      if ( mod(l,100).eq.0 ) then
        nq = end_grid - init_grid + 1
      else
        nq = min0( end_grid - idim_station_sph + 1,
     &                   end_grid - init_grid + 1 )
      endif
      call dclisb(whole_matrix_sph(1,init_grid),
     &                  end_grid - init_grid + 1,
     &                  3,4,ns,nq,whole_vector_sph(init_grid,m),
     &                  eps,whole_matrix_dr_sph(init_grid),
     &                  work_vector(init_grid),ier)
      do 200 m=max0(-l,-2)+1,min0(l,2)
        call dcsbsub(whole_matrix_sph(1,init_grid),
     &                     end_grid - init_grid + 1,
     &                     3,4,ns,nq,whole_vector_sph(init_grid,m),
     &                     eps,whole_matrix_dr_sph(init_grid),
     &                     work_vector(init_grid),ier)
  200 continue
c
c **********************************************************************
c **********************************************************************
c Toroidal Component
c **********************************************************************
c **********************************************************************
c **********************************************************************
c constructing the whole_matrix
c **********************************************************************
      factor_qmu(2)
     &        = unelastic_factor( dble(omega),grid_qmu(idim1_tor) )
      npos = init_npos_tor
      whole_matrix_tor(2,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(1,idim1_tor)
     &                    )
     &            - factor_qmu(2)
     &              * dcmplx( submatrix_I2(1,idim1_tor)
     &              - 2.d0 * submatrix_I4(1,idim1_tor)
     &                     + submatrix_I6(1,idim1_tor)
     &                     + ( lsq2 - 2.d0 ) * submatrix_I7(1,idim1_tor)
     &                      )
      do 320 ir=idim1_tor+1,idim2_tor-1
        factor_qmu(1)
     &          = unelastic_factor( dble(omega),grid_qmu(ir-1) )
        factor_qmu(2)
     &          = unelastic_factor( dble(omega),grid_qmu(ir) )
        npos = npos + 1
        whole_matrix_tor(1,npos)
     &          = omega * omega
     &              * dcmplx(  submatrix_I0(2,ir-1)
     &                      )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I2(2,ir-1)
     &                            - submatrix_I4(2,ir-1)
     &                            - submatrix_I4(3,ir-1)
     &                            + submatrix_I6(2,ir-1)
     &                            + ( lsq2 - 2.d0 )
     &                              * submatrix_I7(2,ir-1)
     &                        )
        whole_matrix_tor(2,npos)
     &          = omega * omega
     &              * dcmplx(   submatrix_I0(4,ir-1)
     &                        + submatrix_I0(1,ir)
     &                      )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I2(4,ir-1)
     &                            - 2.d0 * submatrix_I4(4,ir-1)
     &                            + submatrix_I6(4,ir-1)
     &                            + ( lsq2 - 2.d0 )
     &                              * submatrix_I7(4,ir-1)
     &                        )
     &              - factor_qmu(2)
     &                * dcmplx( submatrix_I2(1,ir)
     &                            - 2.d0 * submatrix_I4(1,ir)
     &                            + submatrix_I6(1,ir)
     &                            + ( lsq2 - 2.d0 )
     &                              * submatrix_I7(1,ir)
     &                        )
  320 continue
      factor_qmu(1)
     &        = unelastic_factor( dble(omega),
     &                            grid_qmu(idim2_tor-1) )
      npos = npos + 1
      whole_matrix_tor(1,npos)
     &        = omega * omega
     &            * dcmplx(   submatrix_I0(2,idim2_tor-1)
     &                    )
     &            - factor_qmu(1)
     &              * dcmplx( submatrix_I2(2,idim2_tor-1)
     &                        - submatrix_I4(2,idim2_tor-1)
     &                        - submatrix_I4(3,idim2_tor-1)
     &                        + submatrix_I6(2,idim2_tor-1)
     &                        + ( lsq2 - 2.d0 )
     &                          * submatrix_I7(2,idim2_tor-1) )
      whole_matrix_tor(2,npos)
     &        = omega * omega
     &            * dcmplx(  submatrix_I0(4,idim2_tor-1)
     &                    )
     &            - factor_qmu(1)
     &              * dcmplx( submatrix_I2(4,idim2_tor-1)
     &                        - 2.d0 * submatrix_I4(4,idim2_tor-1)
     &                        + submatrix_I6(4,idim2_tor-1)
     &                        + ( lsq2 - 2.d0 )
     &                          * submatrix_I7(4,idim2_tor-1) )
      ndim_whole_matrix_tor = npos
c **********************************************************************
c computing the wavefield
c **********************************************************************
      m = max0(-l,-2)
      init_grid = init_npos_tor
      end_grid = ndim_whole_matrix_tor
      ns = idim_rs_tor - init_grid + 1
      if ( mod(l,100).eq.0 ) then
        nq = end_grid - init_grid + 1
      else
        nq = min0( end_grid - idim_station_tor + 1,
     &                   end_grid - init_grid + 1 )
      endif
      call dclisb(whole_matrix_tor(1,init_grid),
     &                  end_grid - init_grid + 1,
     &                  1,2,ns,nq,whole_vector_tor(init_grid,m),
     &                  eps,whole_matrix_dr_tor(init_grid),
     &                  work_vector(init_grid),ier)
      do 400 m=max0(-l,-2)+1,min0(l,2)
        call dcsbsub(whole_matrix_tor(1,init_grid),
     &                     end_grid - init_grid + 1,
     &                     1,2,ns,nq,whole_vector_tor(init_grid,m),
     &                     eps,whole_matrix_dr_tor(init_grid),
     &                     work_vector(init_grid),ier)
  400 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_wavefield0
     &        ( maxngrid_r,omega,
     &          submatrix_I0,submatrix_I1k,submatrix_I1m,
     &          submatrix_I2,submatrix_I3k,submatrix_I3m,
     &          submatrix_I4,submatrix_I5k,submatrix_I5m,
     &          submatrix_I6,submatrix_I7,
     &          submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &          ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &          idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &          idim0,init_npos_sph,init_npos_tor,
     &          idim_rs_sph,idim_rs_tor,
     &          idim_station_sph,idim_station_tor,
     &          whole_matrix_sph,whole_matrix_tor,
     &          whole_matrix_dr_sph,whole_matrix_dr_tor,
     &          whole_vector_sph,whole_vector_tor,work_vector )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing wavefield for the given frequency.
c    required subroutines: init_complex_array,error_handling
c    required functions: unelastic_factor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,ngrid_r,l
      integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
      integer idim0,init_npos_sph,init_npos_tor
      integer idim_rs_sph,idim_rs_tor
      integer idim_station_sph,idim_station_tor
      real*8 submatrix_I0(4,maxngrid_r)
      real*8 submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
      real*8 submatrix_I2(4,maxngrid_r)
      real*8 submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
      real*8 submatrix_I4(4,maxngrid_r)
      real*8 submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
      real*8 submatrix_I6(4,maxngrid_r)
      real*8 submatrix_I7(4,maxngrid_r)
      real*8 submatrix_I3k_mod(6,maxngrid_r)
      real*8 submatrix_I3m_mod(6,maxngrid_r)
      real*8 submatrix_I4_mod(6,maxngrid_r)
      real*8 grid_r(*),grid_mu(2,*)
      real*8 grid_qkappa(*),grid_qmu(*)
      complex*16 omega
      complex*16 whole_matrix_sph(4,*),whole_matrix_tor(2,*)
      complex*16 whole_matrix_dr_sph(*),whole_matrix_dr_tor(*)
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
      complex*16 work_vector(*)
c other variables
      integer ir,npos,m,itype_medium
c itype_medium=1: solid, itype_medium=0: liquid
      integer ndim_whole_matrix_sph,ndim_whole_matrix_tor,ier
      integer init_grid,end_grid,ns,nq
      real*8 lsq,lsq2,eps
      complex*16 unelastic_factor,factor_qkappa(2),factor_qmu(2)
      data eps / -1.d0 /
c
      if ( l.ne.0 ) call error_handling(54)
c **********************************************************************
c Initialing the whole_matrix
c **********************************************************************
      call init_complex_array( 8*maxngrid_r,whole_matrix_sph )
      call init_complex_array( 2*maxngrid_r,whole_matrix_dr_sph )
      idim0 = 1
      init_npos_sph = 1
      init_npos_tor = 1
c **********************************************************************
c **********************************************************************
c Spheroidal Component
c **********************************************************************
c **********************************************************************
c **********************************************************************
c constructing the whole_matrix
c **********************************************************************
      lsq2 = dble(l) * dble(l+1)
      lsq  = dsqrt( lsq2 )
c
      factor_qkappa(2)
     &        = unelastic_factor( dble(omega),grid_qkappa(idim1_sph) )
      factor_qmu(2)
     &        = unelastic_factor( dble(omega),grid_qmu(idim1_sph) )
      if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).eq.0.d0 ) then
        npos = 1
        itype_medium = 0
        whole_matrix_sph(2,npos)
     &        = omega * omega / factor_qkappa(2)
     &            * dcmplx( submatrix_I0(1,idim1_sph)
     &                    )
     &          - dcmplx( lsq2 * submatrix_I1k(1,idim1_sph)
     &                        + submatrix_I2(1,idim1_sph)
     &                  )
      else
        npos = 1
        itype_medium = 1
        whole_matrix_sph(2,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(1,idim1_sph)
     &                    )
     &            - factor_qkappa(2)
     &              * dcmplx( submatrix_I1k(1,idim1_sph)
     &                        + 4.d0 * submatrix_I3k(1,idim1_sph)
     &                        + 4.d0 * submatrix_I5k(1,idim1_sph)
     &                      )
     &            - factor_qmu(2)
     &              * dcmplx( submatrix_I1m(1,idim1_sph)
     &                        + 4.d0 * submatrix_I3m(1,idim1_sph)
     &                        + 4.d0 * submatrix_I5m(1,idim1_sph)
     &                        + lsq2 * submatrix_I6(1,idim1_sph)
     &                        - 4.d0 * submatrix_I7(1,idim1_sph)
     &                      )
      endif
      do 120 ir=idim1_sph+1,idim2_sph-1
        factor_qkappa(1)
     &            = unelastic_factor( dble(omega),grid_qkappa(ir-1) )
        factor_qmu(1)
     &            = unelastic_factor( dble(omega),grid_qmu(ir-1) )
        factor_qkappa(2)
     &          = unelastic_factor( dble(omega),grid_qkappa(ir) )
        factor_qmu(2)
     &          = unelastic_factor( dble(omega),grid_qmu(ir) )
        if ( grid_mu(1,ir)*grid_mu(2,ir).eq.0.d0 ) then
          if ( itype_medium.eq.1 ) then
            npos = npos + 1
            whole_matrix_sph(1,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &                * dcmplx( submatrix_I1k(2,ir-1)
     &                      + 2.d0 * submatrix_I3k(2,ir-1)
     &                      + 2.d0 * submatrix_I3k(3,ir-1)
     &                      + 4.d0 * submatrix_I5k(2,ir-1)
     &                        )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I1m(2,ir-1)
     &                      + 2.d0 * submatrix_I3m(2,ir-1)
     &                      + 2.d0 * submatrix_I3m(3,ir-1)
     &                      + 4.d0 * submatrix_I5m(2,ir-1)
     &                      + lsq2 * submatrix_I6(2,ir-1)
     &                      - 4.d0 * submatrix_I7(2,ir-1)
     &                    )
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &                * dcmplx( submatrix_I1k(4,ir-1)
     &                          + 4.d0 * submatrix_I3k(4,ir-1)
     &                          + 4.d0 * submatrix_I5k(4,ir-1)
     &                        )
     &              - factor_qmu(1)
     &                * dcmplx( submatrix_I1m(4,ir-1)
     &                          + 4.d0 * submatrix_I3m(4,ir-1)
     &                          + 4.d0 * submatrix_I5m(4,ir-1)
     &                          + lsq2 * submatrix_I6(4,ir-1)
     &                          - 4.d0 * submatrix_I7(4,ir-1)
     &                         )
            npos = npos + 1
            itype_medium = 0
            whole_matrix_sph(1,npos)
     &            = omega * grid_r(ir) * grid_r(ir)
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * (
     &                dcmplx( submatrix_I0(1,ir)   ) / factor_qkappa(2)
     &                )
     &            - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &                        + submatrix_I2(1,ir)
     &                    )
          else
            npos = npos + 1
            whole_matrix_sph(1,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(2,ir-1) )
     &              - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &                        + submatrix_I2(2,ir-1)
     &                      )
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * (
     &                dcmplx( submatrix_I0(4,ir-1) ) / factor_qkappa(1)
     &              + dcmplx( submatrix_I0(1,ir)   ) / factor_qkappa(2)
     &                )
     &            - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &                        + submatrix_I2(4,ir-1)
     &                    )
     &            - dcmplx( lsq2 * submatrix_I1k(1,ir)
     &                        + submatrix_I2(1,ir)
     &                    )
          endif
        else
          if ( itype_medium.eq.0 ) then
            npos = npos + 1
            whole_matrix_sph(1,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                      )
     &              - dcmplx( lsq2 * submatrix_I1k(2,ir-1)
     &                        + submatrix_I2(2,ir-1)
     &                      )
            whole_matrix_sph(2,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                      )
     &            - dcmplx( lsq2 * submatrix_I1k(4,ir-1)
     &                        + submatrix_I2(4,ir-1)
     &                    )
            npos = npos + 1
            itype_medium = 1
            whole_matrix_sph(1,npos)
     &            = - omega * grid_r(ir) * grid_r(ir)
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(1,ir)
     &                      )
     &              - factor_qkappa(2)
     &              * dcmplx( submatrix_I1k(1,ir)
     &                        + 4.d0 * submatrix_I3k(1,ir)
     &                        + 4.d0 * submatrix_I5k(1,ir)
     &                      )
     &              - factor_qmu(2)
     &              * dcmplx( submatrix_I1m(1,ir)
     &                        + 4.d0 * submatrix_I3m(1,ir)
     &                        + 4.d0 * submatrix_I5m(1,ir)
     &                        + lsq2 * submatrix_I6(1,ir)
     &                        - 4.d0 * submatrix_I7(1,ir)
     &                      )
          else
            npos = npos + 1
            itype_medium = 1
            whole_matrix_sph(1,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(2,ir-1)
     &                      )
     &              - factor_qkappa(1)
     &              * dcmplx( submatrix_I1k(2,ir-1)
     &                        + 2.d0 * submatrix_I3k(2,ir-1)
     &                        + 2.d0 * submatrix_I3k(3,ir-1)
     &                        + 4.d0 * submatrix_I5k(2,ir-1)
     &                      )
     &              - factor_qmu(1)
     &              * dcmplx( submatrix_I1m(2,ir-1)
     &                        + 2.d0 * submatrix_I3m(2,ir-1)
     &                        + 2.d0 * submatrix_I3m(3,ir-1)
     &                        + 4.d0 * submatrix_I5m(2,ir-1)
     &                        + lsq2 * submatrix_I6(2,ir-1)
     &                        - 4.d0 * submatrix_I7(2,ir-1)
     &                      )
            whole_matrix_sph(2,npos)
     &            = omega * omega
     &              * dcmplx( submatrix_I0(4,ir-1)
     &                        + submatrix_I0(1,ir)
     &                      )
     &              - factor_qkappa(1)
     &              * dcmplx( submatrix_I1k(4,ir-1)
     &                        + 4.d0 * submatrix_I3k(4,ir-1)
     &                        + 4.d0 * submatrix_I5k(4,ir-1)
     &                      )
     &              - factor_qmu(1)
     &              * dcmplx( submatrix_I1m(4,ir-1)
     &                        + 4.d0 * submatrix_I3m(4,ir-1)
     &                        + 4.d0 * submatrix_I5m(4,ir-1)
     &                        + lsq2 * submatrix_I6(4,ir-1)
     &                        - 4.d0 * submatrix_I7(4,ir-1)
     &                      )
     &              - factor_qkappa(2)
     &              * dcmplx( submatrix_I1k(1,ir)
     &                        + 4.d0 * submatrix_I3k(1,ir)
     &                        + 4.d0 * submatrix_I5k(1,ir)
     &                      )
     &              - factor_qmu(2)
     &              * dcmplx( submatrix_I1m(1,ir)
     &                        + 4.d0 * submatrix_I3m(1,ir)
     &                        + 4.d0 * submatrix_I5m(1,ir)
     &                        + lsq2 * submatrix_I6(1,ir)
     &                        - 4.d0 * submatrix_I7(1,ir)
     &                      )
          endif
        endif
  120 continue
      factor_qkappa(1)
     &        = unelastic_factor( dble(omega),grid_qkappa(idim2_sph-1) )
      factor_qmu(1)
     &        = unelastic_factor( dble(omega),grid_qmu(idim2_sph-1) )
      if ( itype_medium.eq.0 ) then
        npos = npos + 1
        whole_matrix_sph(1,npos)
     &          = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(2,idim2_sph-1)
     &                      )
     &            - dcmplx( lsq2 * submatrix_I1k(2,idim2_sph-1)
     &                        + submatrix_I2(2,idim2_sph-1)
     &                    )
        whole_matrix_sph(2,npos)
     &            = omega * omega / factor_qkappa(1)
     &              * dcmplx( submatrix_I0(4,idim2_sph-1)
     &                      )
     &            - dcmplx( lsq2 * submatrix_I1k(4,idim2_sph-1)
     &                      + submatrix_I2(4,idim2_sph-1)
     &                    )
        ndim_whole_matrix_sph = npos
      else
        npos = npos + 1
        whole_matrix_sph(1,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(2,idim2_sph-1)
     &                    )
     &          - factor_qkappa(1)
     &            * dcmplx( submatrix_I1k(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3k(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3k(3,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5k(2,idim2_sph-1)
     &                    )
     &          - factor_qmu(1)
     &            * dcmplx( submatrix_I1m(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3m(2,idim2_sph-1)
     &                      + 2.d0 * submatrix_I3m(3,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5m(2,idim2_sph-1)
     &                      + lsq2 * submatrix_I6(2,idim2_sph-1)
     &                      - 4.d0 * submatrix_I7(2,idim2_sph-1)
     &                    )
        whole_matrix_sph(2,npos)
     &        = omega * omega
     &            * dcmplx( submatrix_I0(4,idim2_sph-1)
     &                    )
     &          - factor_qkappa(1)
     &            * dcmplx( submatrix_I1k(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I3k(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5k(4,idim2_sph-1)
     &                    )
     &          - factor_qmu(1)
     &            * dcmplx( submatrix_I1m(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I3m(4,idim2_sph-1)
     &                      + 4.d0 * submatrix_I5m(4,idim2_sph-1)
     &                      + lsq2 * submatrix_I6(4,idim2_sph-1)
     &                      - 4.d0 * submatrix_I7(4,idim2_sph-1)
     &                    )
        ndim_whole_matrix_sph = npos
      endif
c **********************************************************************
c computing the wavefield
c **********************************************************************
c imposing fixed boundary conditions at r=0 and free surface boundary
c conditions at the surface
      if ( ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph).ne.0.d0 )
     &           .and.( grid_r(idim1_sph).eq.0.d0 ) ) then
        init_grid = 2
      else
        init_grid = 1
      endif
      if ( grid_mu(1,idim2_sph-1)*grid_mu(2,idim2_sph-1).eq.0.d0 )
     &        then
        end_grid = ndim_whole_matrix_sph - 1
      else
        end_grid = ndim_whole_matrix_sph
      endif
      m = 0
      ns = idim_rs_sph - init_grid + 1
      if ( mod(l,100).eq.0 ) then
        nq = end_grid - init_grid + 1
      else
        nq = end_grid - idim_station_sph + 1
      endif
      call dclisb(whole_matrix_sph(1,init_grid),
     &                  end_grid - init_grid + 1,
     &                  1,4,ns,nq,whole_vector_sph(init_grid,m),
     &                  eps,whole_matrix_dr_sph(init_grid),
     &                  work_vector(init_grid),ier)
c
      ndim_whole_matrix_tor = 0
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function unelastic_factor( omega,qmu )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing the unelastic factor (ratio between complex mu and
c real mu) for the given quality factor, qmu.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c input/output variables
      real*8 omega,qmu
c other variables
      real*8 vr,vi
c constants
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
c **********************************************************************
c computing the complex velocity for the given qmu
c **********************************************************************
      if ( (omega.eq.0.d0).or.(qmu.lt.0.d0) ) then
        vr = 1.d0
        vi = 0.d0
      else
        vr = 1.d0
     &             + dlog( omega / ( 2.d0 * pi ) ) / ( pi * qmu )
        vi = 1.d0 / ( 2.d0 * qmu )
      endif
c
c **********************************************************************
c computing the unelastic factor
c **********************************************************************
      unelastic_factor = dcmplx( vr, vi ) * dcmplx( vr, vi )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dclisb(a, n, nud, n1, np, nq, b, eps, dr, z, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
*  simultaneous linear equations with real symmetric positive definite *
*      band matrix by cholesky method.                                 *
*  parameters                                                          *
*    (1) a : 2-dim. array containing the matrix.                       *
*    (2) n : order of the matrix.                                      *
*    (3) nud : size of band's half width.                              *
*    (4) n1 : row size of the array a in the 'dimension' statement.    *
*    (5) b : 1-dim. array containing the right hand side vector.       *
*    (6) eps : parameter to check singurarity off the matrix           *
*              standard value = 1.0d-14                                *
*    (7) dr : 1-dim. working array.                                    *
*    (8) z : 1-dim. working array.                                     *
*    (9) ier : error code.                                             *
*  copy right   t. oguni   july 30 1989   version 1.0                  *
************************************************************************
      complex*16 a(n1,n), b(n), dr(n), z(n)
      real*8 eps
      integer n, nud, n1, np, nq, ier
      complex*16 xx, s, sum, au, t
      real*8 eps1
      integer i ,m, j, k1, mj, i1, k, j1
c  check the input data
      ier = 0
      eps1 = 1.0d-14
      m = nud + 1
      if ((n .le. 0) .or. (nud .le. 0 ) .or. (n1 .lt. m)) then
       ier = 2
       write(*,*) '(subr. lisb) invalid argument. ', n, nud, n1
       return
      endif
      if (eps .le. 0.0) eps = eps1
c  modified cholesky decomposition
      j = 1
      if (cdabs(a(m,1)) .le. eps) then
       ier = 1
       write(*,*) '(subr. lisb) singular at step # ', j
       return
      endif
      dr(1) = dcmplx(1.0d0) / a(m,1)
      xx = a(m-1,2)
      a(m-1,2) = a(m-1,2) * dr(1)
      s = a(m,2) - xx * a(m-1,2)
      j = 2
      if (cdabs(s) .le. eps) then
       ier = 1
       write(*,*) '(subr. lisb) singular at step # ', j
       return
      endif
      dr(2) = dcmplx(1.0d0) / s
      if (m .lt. 3) then
       do 5 j=3,n
        xx = a(1,j)
        a(1,j) = xx * dr(j-1)
        s = a(2,j) - xx * a(1,j)
        if (cdabs(s) .le. eps) then
         ier = 1
         write(*,*) ' (subr. lisb) singular at step # ', j
         return
        endif
        dr(j) = dcmplx(1.0d0) / s
    5  continue
      else
       do 30 j=3,n
        k1 = 1
        if (j .ge. m) k1 = j - m + 1
        mj = m - j
        do 20 i=k1+1,j-1
         sum = dcmplx(0.0d0)
         do 10 k=k1,i-1
   10     sum = sum + a(m-i+k,i) * a(mj+k,j)
         a(mj+i,j) = a(mj+i,j) - sum
   20   continue
        sum = dcmplx(0.0d0)
        do 25 i=k1,j-1
         xx = a(mj+i,j)
         au = xx * dr(i)
         sum = sum + xx *au
         a(mj+i,j) = au
   25   continue
        t = a(m,j) - sum
        if (cdabs(t) .le. eps) then
         ier = 1
         write(*,*) ' (subr. lisb) singular at step # ', j
         return
        endif
        dr(j) = dcmplx(1.0d0) / t
   30  continue
      endif
c subtitution
      entry dcsbsub(a, n, nud, n1, np, nq, b, eps, dr, z, ier)
c  forward substitution
      m = nud + 1
      if (m .lt. 3) then
        z(np) = b(np)
        do 40 j=np+1,n
   40   z(j) = b(j) - a(1,j) * z(j-1)
        do 45 j=1,np-1
   45   z(j) = dcmplx(0.d0)
        do 50 j=np,n
   50   z(j) = z(j) * dr(j)
        b(n) = z(n)
c       do 60 j=1,n-1
        do 60 j=1,nq-1
   60   b(n-j) = z(n-j) - a(1,n-j+1) * b(n-j+1)
      else
        z(np) = b(np)
        z(np+1) = b(np+1) - a(m-1,np+1) * z(np)
        do 80 j=np+2,n
          if (j .gt. np-1+m) then
            i1 = 1
          else
            i1 = np-1+m - j + 1
          endif
          sum = dcmplx(0.0d0)
          do 70 k=i1,m-1
   70       sum = sum + a(k,j) * z(j-m+k)
   80   z(j) = b(j) - sum
        do 85 j=1,np-1
   85   z(j) = dcmplx(0.d0)
        do 90 j=np,n
   90   z(j) = z(j) * dr(j)
c
        b(n) = z(n)
        b(n-1) = z(n-1) - a(m-1,n) * z(n)
        do 110 j=3,nq
          j1 = n - j + 1
          i1 = 1
          if (j .lt. m) i1 = m - j + 1
          sum = dcmplx(0.0d0)
          do 100 k=i1,m-1
  100     sum = sum + a(k,m-k+j1) * b(m-k+j1)
          b(j1) = z(j1) - sum
  110   continue
      endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_complex_array
     &                   ( n_station,station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initializing the acculuated displacement at the station.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer n_station
      complex*16 station_displacement(*)
c other variables
      integer i_station
c
      do 100 i_station=1,n_station
        station_displacement(i_station) = dcmplx(0.d0)
  100 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_displacement_station
     &          ( maxngrid_r,whole_vector_tor,whole_vector_sph,
     &            l,n_station,station_theta,station_phi,
     &            idim_station_sph,idim_station_tor,
     &            station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c accumulating the displacement at the station.
c    required subroutines: error_handling,comp_vecsph_tor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,l,n_station
      integer idim_station_sph,idim_station_tor
      real*8 station_theta(*),station_phi(*)
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
      complex*16 station_displacement(3,*)
c other variables
      integer m,i_station,icomp
      real*8 lsq
      complex*16 vecsph_sph1(3),vecsph_sph2(3),vecsph_tor(3)
c
      if ( l.le.0 ) call error_handling(55)
      lsq = dsqrt( dble(l) * dble(l+1) )
      do 200 i_station=1,n_station
c ---- computing the value of the trial functions at the station
        do 120 m=max0(-l,-2),min0(l,2)
c -------- horizontal dependent part
          call comp_vecsph(l,m,
     &                           station_theta(i_station),
     &                           station_phi(i_station),
     &                           vecsph_sph1,vecsph_sph2,vecsph_tor)
c ---- computing the displacement at the station
          do 110 icomp=1,3
            station_displacement(icomp,i_station)
     &            = station_displacement(icomp,i_station)
     &              + whole_vector_sph(idim_station_sph,m)
     &                * vecsph_sph1(icomp)
     &              + whole_vector_sph(idim_station_sph+1,m)
     &                * vecsph_sph2(icomp) / dcmplx( lsq )
     &              + whole_vector_tor(idim_station_tor,m)
     &                * vecsph_tor(icomp) / dcmplx( lsq )
  110     continue
  120   continue
  200 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_displacement_station0
     &          ( maxngrid_r,whole_vector_tor,whole_vector_sph,
     &            l,n_station,station_theta,station_phi,
     &            idim_station_sph,idim_station_tor,
     &            station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c accumulating the displacement at the station.
c    required subroutines: error_handling,comp_vecsph_tor
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,l,n_station
      integer idim_station_sph,idim_station_tor
      real*8 station_theta(*),station_phi(*)
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
      complex*16 station_displacement(3,*)
c other variables
      integer m,i_station,icomp
      real*8 lsq
      complex*16 vecsph_sph1(3),vecsph_sph2(3),vecsph_tor(3)
c
      if ( l.ne.0 ) call error_handling(56)
      lsq = dsqrt( dble(l) * dble(l+1) )
      do 200 i_station=1,n_station
c ---- computing the value of the trial functions at the station
        do 120 m=max0(-l,-2),min0(l,2)
c -------- horizontal dependent part
          call comp_vecsph(l,m,
     &                           station_theta(i_station),
     &                           station_phi(i_station),
     &                           vecsph_sph1,vecsph_sph2,vecsph_tor)
c ---- computing the displacement at the station
          do 110 icomp=1,3
            station_displacement(icomp,i_station)
     &            = station_displacement(icomp,i_station)
     &              + whole_vector_sph(idim_station_sph,m)
     &                * vecsph_sph1(icomp)
  110     continue
  120   continue
  200 continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_vecsph(l,m,theta,phi,
     &                             vecsph_sph1,vecsph_sph2,vecsph_tor)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing the vector spherical harmonics.
c   required subroutines: error_handling
c   required functions: plgndr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer l,m
      real*8 theta,phi
      complex*16 vecsph_sph1(3),vecsph_sph2(3),vecsph_tor(3)
c other variables
      integer i,m0
      real*8 factor,x,plgndr
      complex*16 expimp
c constant
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
c **********************************************************************
c checking the argumants
c **********************************************************************
      m0 = iabs(m)
      if ( (l.lt.0).or.(m0.gt.l) ) call error_handling(41)
c **********************************************************************
c computing the normalization factor (including the sign)
c **********************************************************************
      factor = 1.d0
      do 100 i=l-m0+1,l+m0
        factor = factor * dble(i)
  100 continue
      factor = dsqrt( dble(2*l+1)/(4.d0*pi) / factor )
      if ( ( m0.ne.m ).and.( mod(m0,2).eq.1 ) ) factor = - factor
c **********************************************************************
c computing each component of the vector spherical harmonics
c **********************************************************************
      x = dcos(theta)
      expimp = cdexp( dcmplx( 0.d0, dble(m)*phi ) )
c
      vecsph_sph1(1) = factor
     &                       * plgndr(l,m0,dcos(theta))
     &                       * expimp
      vecsph_sph1(2) = 0.d0
      vecsph_sph1(3) = 0.d0
c
      vecsph_sph2(1) = 0.d0
      if ( l.ge.m0+1 ) then
        vecsph_sph2(2) = factor
     &                        * (   dble(m0) * x / dsin(theta)
     &                              * plgndr(l,m0,x)
     &                            + plgndr(l,m0+1,x) )
     &                        * expimp
      else
        vecsph_sph2(2) = factor
     &                        * (   dble(m0) * x / dsin(theta)
     &                              * plgndr(l,m0,x) )
     &                        * expimp
      endif
      vecsph_sph2(3) = factor
     &                      * dcmplx( 0.d0, dble(m) ) / dsin(theta)
     &                      * plgndr(l,m0,dcos(theta))
     &                      * expimp
c
      vecsph_tor(1) = dcmplx(0.d0)
      vecsph_tor(2) = vecsph_sph2(3)
      vecsph_tor(3) = - vecsph_sph2(2)
c     vecsph_tor(2) = factor
c     &                     * dcmplx( 0.d0, dble(m0) ) / dsin(theta)
c     &                     * plgndr(l,m0,dcos(theta))
c     &                     * expimp
c     if ( l.ge.m0+1 ) then
c       vecsph_tor(3) = - factor
c     &                       * (   dble(m0) * x / dsin(theta)
c     &                             * plgndr(l,m0,x)
c     &                           + plgndr(l,m0+1,x) )
c     &                       * expimp
c     else
c       vecsph_tor(3) = - factor
c     &                       * (   dble(m) * x / dsin(theta)
c     &                             * plgndr(l,m0,x) )
c     &                       * expimp
c     endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function plgndr(l,m,x)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing the associated Legendre polynominal
c   (from Numerical Recipies).
c   required subroutines: error_handling
c   required functions: plgndr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer l,m
      real*8 x
      integer i,ll
      real*8 fact,pll,pmm,pmmp1,somx2
c
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) call error_handling(42)
      pmm=1.d0
      if(m.gt.0) then
        somx2=dsqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*dble(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
c
      return
      end
