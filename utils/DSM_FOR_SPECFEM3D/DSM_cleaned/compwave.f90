
  subroutine comp_excitation( maxngrid_r,omega, &
            l,source_r,source_mt, &
            igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms, &
            grid_qkappas,grid_qmus, &
            submatrix_I0,submatrix_I1k,submatrix_I1m, &
            submatrix_I2,submatrix_I3k,submatrix_I3m, &
            submatrix_I4,submatrix_I5k,submatrix_I5m, &
            submatrix_I6,submatrix_I7, &
            submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod, &
            idim_rs_sph,idim_rs_tor, &
            whole_vector_sph,whole_vector_tor )

! compute a excitation vector for the given frequency

  implicit none

! variables for input/output
  integer maxngrid_r,l,igrid_rs
  integer idim_rs_sph,idim_rs_tor
  real(kind=8) source_r,source_mt(3,3)
  real(kind=8) grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
  real(kind=8) grid_qkappas,grid_qmus
  real(kind=8) submatrix_I0(4,maxngrid_r)
  real(kind=8) submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
  real(kind=8) submatrix_I2(4,maxngrid_r)
  real(kind=8) submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
  real(kind=8) submatrix_I4(4,maxngrid_r)
  real(kind=8) submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
  real(kind=8) submatrix_I6(4,maxngrid_r)
  real(kind=8) submatrix_I7(4,maxngrid_r)
  real(kind=8) submatrix_I3k_mod(6,maxngrid_r)
  real(kind=8) submatrix_I3m_mod(6,maxngrid_r)
  real(kind=8) submatrix_I4_mod(6,maxngrid_r)
  complex(kind=8) omega
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) whole_vector_tor(maxngrid_r,-2:2)

! other variables
  real(kind=8) b1,b2,lsq2,lsq
  complex(kind=8) D1,D2_p,D2_m,D3_p,D3_m,anelastic_factor
  complex(kind=8) factors_qkappa,factors_qmu

! constant
  real(kind=8), parameter :: pi=3.1415926535897932d0

  if ( l <= 0 ) call error_handling(51)

! **********************************************************************
! initialize the whole_vector
! **********************************************************************
  call init_complex_array( 10*maxngrid_r,whole_vector_sph )
  call init_complex_array(  5*maxngrid_r,whole_vector_tor )

! **********************************************************************
! compute the excitation vector
! **********************************************************************
  b1 = dsqrt( dble(2*l+1) / ( 16.d0 * pi ) )
  b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2) / ( 64.d0 * pi ) )
  factors_qkappa = anelastic_factor( dble(omega), grid_qkappas )
  factors_qmu = anelastic_factor( dble(omega), grid_qmus )
  lsq2 = dble(l) * dble(l+1)
  lsq = dsqrt( lsq2 )

! --- spheroidal excitations due to the traction discontinuities
  whole_vector_sph(idim_rs_sph,0) = whole_vector_sph(idim_rs_sph,0) &
            - b1 * 2.d0 * ( source_mt(2,2) + source_mt(3,3) &
                  - 2.d0 * source_mt(1,1) &
                    * ( grid_Fks * factors_qkappa &
                        + grid_Fms * factors_qmu ) &
                    / ( grid_Cks * factors_qkappa &
                        + grid_Cms * factors_qmu ) ) / source_r

  whole_vector_sph(idim_rs_sph+1,0) = whole_vector_sph(idim_rs_sph+1,0) &
            - b1 * lsq * ( - source_mt(2,2) - source_mt(3,3) &
                  + 2.d0 * source_mt(1,1) &
                    * ( grid_Fks * factors_qkappa &
                        + grid_Fms * factors_qmu ) &
                    / ( grid_Cks * factors_qkappa &
                        + grid_Cms * factors_qmu ) &
                  ) / source_r

  whole_vector_sph(idim_rs_sph+1,-2) = whole_vector_sph(idim_rs_sph+1,-2) &
            + b2 * cmplx( - source_mt(2,2) + source_mt(3,3), &
                        - 2.d0 * source_mt(2,3) ) / source_r

  whole_vector_sph(idim_rs_sph+1,2) = whole_vector_sph(idim_rs_sph+1,2) &
            + b2 * cmplx( - source_mt(2,2) + source_mt(3,3), &
                          2.d0 * source_mt(2,3) ) / source_r

! --- toroidal excitations due to the traction discontinuities

  whole_vector_tor(idim_rs_tor,-2) = whole_vector_tor(idim_rs_tor,-2) &
            - b2 * ( cmplx( - 2.d0 * source_mt(2,3), &
                               source_mt(2,2) - source_mt(3,3) ) ) / source_r

  whole_vector_tor(idim_rs_tor,2) = whole_vector_tor(idim_rs_tor,2) &
            - b2 * ( cmplx( - 2.d0 * source_mt(2,3), &
                             - source_mt(2,2) + source_mt(3,3) ) ) / source_r

! --- excitations due to the displacement discontinuities

  D1 = b1 * 2.d0 * source_mt(1,1) / ( source_r * source_r * ( grid_Cks * factors_qkappa + grid_Cms * factors_qmu ) )
  D2_p = b1 * cmplx( -source_mt(1,2), source_mt(1,3) ) &
             / ( source_r * source_r * grid_Ls * factors_qmu )
  D2_m = b1 * cmplx( source_mt(1,2), source_mt(1,3) ) &
             / ( source_r * source_r * grid_Ls * factors_qmu )
  D3_p = b1 * cmplx( source_mt(1,3), source_mt(1,2) ) &
             / ( source_r * source_r * grid_Ls * factors_qmu )
  D3_m = b1 * cmplx( -source_mt(1,3), source_mt(1,2) ) &
             / ( source_r * source_r * grid_Ls * factors_qmu )

! ---- spheroidal, m=0

  whole_vector_sph(idim_rs_sph,0) = whole_vector_sph(idim_rs_sph,0) &
              + ( - omega * omega * submatrix_I0(1,igrid_rs) &
                  + ( submatrix_I1k(1,igrid_rs) &
                      + 4.d0 * submatrix_I3k(1,igrid_rs) &
                      + 4.d0 * submatrix_I5k(1,igrid_rs) &
                    ) * factors_qkappa &
                  + ( submatrix_I1m(1,igrid_rs) &
                      + 4.d0 * submatrix_I3m(1,igrid_rs) &
                      + 4.d0 * submatrix_I5m(1,igrid_rs) &
                      + lsq2 * submatrix_I6(1,igrid_rs) &
                      - 4.d0 * submatrix_I7(1,igrid_rs) &
                    ) * factors_qmu ) * D1

  whole_vector_sph(idim_rs_sph+1,0) = whole_vector_sph(idim_rs_sph+1,0) &
              + ( - lsq * ( ( submatrix_I3k_mod(1,igrid_rs) &
                        + 2.d0 * submatrix_I5k(1,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(1,igrid_rs) &
                          - submatrix_I4_mod(1,igrid_rs) &
                          + 2.d0 * submatrix_I5m(1,igrid_rs) &
                          + submatrix_I6(1,igrid_rs) &
                          - 2.d0 * submatrix_I7(1,igrid_rs) &
                      ) * factors_qmu ) ) * D1

  whole_vector_sph(idim_rs_sph+2,0) = whole_vector_sph(idim_rs_sph+2,0) &
              + ( - omega * omega * submatrix_I0(3,igrid_rs) &
                  + ( submatrix_I1k(3,igrid_rs) &
                      + 2.d0 * submatrix_I3k(2,igrid_rs) &
                      + 2.d0 * submatrix_I3k(3,igrid_rs) &
                      + 4.d0 * submatrix_I5k(3,igrid_rs) &
                    ) * factors_qkappa &
                  + ( submatrix_I1m(3,igrid_rs) &
                      + 2.d0 * submatrix_I3m(2,igrid_rs) &
                      + 2.d0 * submatrix_I3m(3,igrid_rs) &
                      + 4.d0 * submatrix_I5m(3,igrid_rs) &
                      + lsq2 * submatrix_I6(3,igrid_rs) &
                      - 4.d0 * submatrix_I7(3,igrid_rs) &
                    ) * factors_qmu ) * D1

  whole_vector_sph(idim_rs_sph+3,0) = whole_vector_sph(idim_rs_sph+3,0) &
              + ( - lsq * ( ( submatrix_I3k_mod(3,igrid_rs) &
                        + 2.d0 * submatrix_I5k(3,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(3,igrid_rs) &
                          - submatrix_I4_mod(3,igrid_rs) &
                          + 2.d0 * submatrix_I5m(3,igrid_rs) &
                          + submatrix_I6(3,igrid_rs) &
                          - 2.d0 * submatrix_I7(3,igrid_rs) &
                      ) * factors_qmu ) ) * D1

! ---- spheroidal, m=-1

  whole_vector_sph(idim_rs_sph,-1) = whole_vector_sph(idim_rs_sph,-1) &
              + ( - lsq * ( ( submatrix_I3k_mod(1,igrid_rs) &
                        + 2.d0 * submatrix_I5k(1,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(1,igrid_rs) &
                          - submatrix_I4_mod(1,igrid_rs) &
                          + 2.d0 * submatrix_I5m(1,igrid_rs) &
                          + submatrix_I6(1,igrid_rs) &
                          - 2.d0 * submatrix_I7(1,igrid_rs) &
                      ) * factors_qmu ) ) * D2_m

  whole_vector_sph(idim_rs_sph+1,-1) = whole_vector_sph(idim_rs_sph+1,-1) &
              + ( - omega * omega * submatrix_I0(1,igrid_rs) &
                  + ( lsq2 * submatrix_I5k(1,igrid_rs) &
                    ) * factors_qkappa &
                    + ( submatrix_I2(1,igrid_rs) &
                        - 2.d0 * submatrix_I4(1,igrid_rs) &
                        + lsq2 * submatrix_I5m(1,igrid_rs) &
                        + submatrix_I6(1,igrid_rs) &
                        - 2.d0 * submatrix_I7(1,igrid_rs) &
                    ) * factors_qmu )  * D2_m

  whole_vector_sph(idim_rs_sph+2,-1) = whole_vector_sph(idim_rs_sph+2,-1) &
              + ( - lsq * ( ( submatrix_I3k_mod(2,igrid_rs) &
                        + 2.d0 * submatrix_I5k(3,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(2,igrid_rs) &
                          - submatrix_I4_mod(2,igrid_rs) &
                          + 2.d0 * submatrix_I5m(3,igrid_rs) &
                          + submatrix_I6(3,igrid_rs) &
                          - 2.d0 * submatrix_I7(3,igrid_rs) &
                      ) * factors_qmu ) ) * D2_m

  whole_vector_sph(idim_rs_sph+3,-1) = whole_vector_sph(idim_rs_sph+3,-1) &
              + ( - omega * omega * submatrix_I0(3,igrid_rs) &
                  + ( lsq2 * submatrix_I5k(3,igrid_rs) &
                    ) * factors_qkappa &
                    + ( submatrix_I2(3,igrid_rs) &
                        - submatrix_I4(2,igrid_rs) &
                        - submatrix_I4(3,igrid_rs) &
                        + lsq2 * submatrix_I5m(3,igrid_rs) &
                        + submatrix_I6(3,igrid_rs) &
                        - 2.d0 * submatrix_I7(3,igrid_rs) &
                    ) * factors_qmu )  * D2_m

  whole_vector_sph(idim_rs_sph+4,-1) = whole_vector_sph(idim_rs_sph+4,-1) &
              + ( - lsq * ( ( submatrix_I3k_mod(5,igrid_rs) ) * factors_qkappa &
                      + ( submatrix_I3m_mod(5,igrid_rs) &
                          - submatrix_I4_mod(5,igrid_rs) ) * factors_qmu ) ) * D2_m

! ---- spheroidal, m=+1

  whole_vector_sph(idim_rs_sph,1) = whole_vector_sph(idim_rs_sph,1) &
              + ( - lsq * ( ( submatrix_I3k_mod(1,igrid_rs) &
                        + 2.d0 * submatrix_I5k(1,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(1,igrid_rs) &
                          - submatrix_I4_mod(1,igrid_rs) &
                          + 2.d0 * submatrix_I5m(1,igrid_rs) &
                          + submatrix_I6(1,igrid_rs) &
                          - 2.d0 * submatrix_I7(1,igrid_rs) &
                      ) * factors_qmu ) ) * D2_p

  whole_vector_sph(idim_rs_sph+1,1) = whole_vector_sph(idim_rs_sph+1,1) &
              + ( - omega * omega * submatrix_I0(1,igrid_rs) &
                  + ( lsq2 * submatrix_I5k(1,igrid_rs) &
                    ) * factors_qkappa &
                    + ( submatrix_I2(1,igrid_rs) &
                        - 2.d0 * submatrix_I4(1,igrid_rs) &
                        + lsq2 * submatrix_I5m(1,igrid_rs) &
                        + submatrix_I6(1,igrid_rs) &
                        - 2.d0 * submatrix_I7(1,igrid_rs) &
                    ) * factors_qmu )  * D2_p

  whole_vector_sph(idim_rs_sph+2,1) = whole_vector_sph(idim_rs_sph+2,1) &
              + ( - lsq * ( ( submatrix_I3k_mod(2,igrid_rs) &
                        + 2.d0 * submatrix_I5k(3,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(2,igrid_rs) &
                          - submatrix_I4_mod(2,igrid_rs) &
                          + 2.d0 * submatrix_I5m(3,igrid_rs) &
                          + submatrix_I6(3,igrid_rs) &
                          - 2.d0 * submatrix_I7(3,igrid_rs) &
                      ) * factors_qmu ) ) * D2_p

  whole_vector_sph(idim_rs_sph+3,1) = whole_vector_sph(idim_rs_sph+3,1) &
              + ( - omega * omega * submatrix_I0(3,igrid_rs) &
                  + ( lsq2 * submatrix_I5k(3,igrid_rs) &
                    ) * factors_qkappa &
                    + ( submatrix_I2(3,igrid_rs) &
                        - submatrix_I4(2,igrid_rs) &
                        - submatrix_I4(3,igrid_rs) &
                        + lsq2 * submatrix_I5m(3,igrid_rs) &
                        + submatrix_I6(3,igrid_rs) &
                        - 2.d0 * submatrix_I7(3,igrid_rs) &
                    ) * factors_qmu )  * D2_p

  whole_vector_sph(idim_rs_sph+4,1) = whole_vector_sph(idim_rs_sph+4,1) &
              + ( - lsq * ( ( submatrix_I3k_mod(5,igrid_rs) &
                      ) * factors_qkappa &
                      + ( submatrix_I3m_mod(5,igrid_rs) &
                          - submatrix_I4_mod(5,igrid_rs) &
                      ) * factors_qmu ) ) * D2_p

! ---- toroidal, m=-1

  whole_vector_tor(idim_rs_tor,-1) = whole_vector_tor(idim_rs_tor,-1) &
              + ( - omega * omega * submatrix_I0(1,igrid_rs) &
                  + ( submatrix_I2(1,igrid_rs) &
                      - 2.d0 * submatrix_I4(1,igrid_rs) &
                      + submatrix_I6(1,igrid_rs) &
                      - ( lsq2 - 2.d0 ) * submatrix_I7(1,igrid_rs) &
                     ) * factors_qmu ) * D3_m

  whole_vector_tor(idim_rs_tor+1,-1) = whole_vector_tor(idim_rs_tor+1,-1) &
              + ( - omega * omega * submatrix_I0(3,igrid_rs) &
                  + ( submatrix_I2(3,igrid_rs) &
                      - submatrix_I4(2,igrid_rs) &
                      - submatrix_I4(3,igrid_rs) &
                      + submatrix_I6(3,igrid_rs) &
                      - ( lsq2 - 2.d0 ) * submatrix_I7(3,igrid_rs) &
                     ) * factors_qmu ) * D3_m

! ---- toroidal, m=+1

  whole_vector_tor(idim_rs_tor,1) = whole_vector_tor(idim_rs_tor,1) &
              + ( - omega * omega * submatrix_I0(1,igrid_rs) &
                  + ( submatrix_I2(1,igrid_rs) &
                      - 2.d0 * submatrix_I4(1,igrid_rs) &
                      + submatrix_I6(1,igrid_rs) &
                      - ( lsq2 - 2.d0 ) * submatrix_I7(1,igrid_rs) &
                     ) * factors_qmu ) * D3_p

  whole_vector_tor(idim_rs_tor+1,1) = whole_vector_tor(idim_rs_tor+1,1) &
              + ( - omega * omega * submatrix_I0(3,igrid_rs) &
                  + ( submatrix_I2(3,igrid_rs) &
                      - submatrix_I4(2,igrid_rs) &
                      - submatrix_I4(3,igrid_rs) &
                      + submatrix_I6(3,igrid_rs) &
                      - ( lsq2 - 2.d0 ) * submatrix_I7(3,igrid_rs) &
                     ) * factors_qmu ) * D3_p

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_excitation0( maxngrid_r,omega, &
            l,source_r,source_mt, &
            igrid_rs,grid_Cks,grid_Cms,grid_Fks,grid_Fms, &
            grid_qkappas,grid_qmus, &
            submatrix_I0,submatrix_I1k,submatrix_I1m, &
            submatrix_I3k,submatrix_I3m, &
            submatrix_I5k,submatrix_I5m, &
            submatrix_I6,submatrix_I7, &
            idim_rs_sph, &
            whole_vector_sph,whole_vector_tor )

! compute a excitation vector for the given frequency

  implicit none

! variables for input/output
  integer maxngrid_r,l,igrid_rs
  integer idim_rs_sph
  real(kind=8) source_r,source_mt(3,3)
  real(kind=8) grid_Cks,grid_Cms,grid_Fks,grid_Fms
  real(kind=8) grid_qkappas,grid_qmus
  real(kind=8) submatrix_I0(4,maxngrid_r)
  real(kind=8) submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
  real(kind=8) submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
  real(kind=8) submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
  real(kind=8) submatrix_I6(4,maxngrid_r)
  real(kind=8) submatrix_I7(4,maxngrid_r)
  complex(kind=8) omega
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) whole_vector_tor(maxngrid_r,-2:2)

! other variables
  real(kind=8) b1,b2,lsq2
  complex(kind=8) D1,D2_p,D2_m,D3_p,D3_m,anelastic_factor
  complex(kind=8) factors_qkappa,factors_qmu

! constant
  real(kind=8), parameter :: pi=3.1415926535897932d0

  if ( l/=0 ) call error_handling(52)

! **********************************************************************
! initialize the whole_vector
! **********************************************************************
  call init_complex_array( 10*maxngrid_r,whole_vector_sph )
  call init_complex_array(  5*maxngrid_r,whole_vector_tor )

! **********************************************************************
! compute the excitation vector
! **********************************************************************
  b1 = dsqrt( dble(2*l+1) / ( 16.d0 * pi ) )
! b2 = dsqrt( dble(2*l+1)*dble(l-1)*dble(l+2) / ( 64.d0 * pi ) )
  b2 = 0.d0

! --- excitations due to the traction discontinuities
  factors_qkappa = anelastic_factor( dble(omega), grid_qkappas )
  factors_qmu = anelastic_factor( dble(omega), grid_qmus )
  lsq2 = dble(l) * dble(l+1)
  whole_vector_sph(idim_rs_sph,0) = whole_vector_sph(idim_rs_sph,0) &
            - b1 * 2.d0 * ( source_mt(2,2) + source_mt(3,3) &
                  - 2.d0 * source_mt(1,1) * ( grid_Fks * factors_qkappa &
                        + grid_Fms * factors_qmu ) / ( grid_Cks * factors_qkappa &
                        + grid_Cms * factors_qmu ) ) / source_r

! --- excitations due to the displacement discontinuities
  D1 = b1 * 2.d0 * source_mt(1,1) / ( source_r * source_r * ( grid_Cks * factors_qkappa + grid_Cms * factors_qmu ) )
  D2_p = cmplx( 0.d0 )
  D2_m = cmplx( 0.d0 )
  D3_p = cmplx( 0.d0 )
  D3_m = cmplx( 0.d0 )

! ---- spheroidal, m=0

  whole_vector_sph(idim_rs_sph,0) = whole_vector_sph(idim_rs_sph,0) &
              + ( - omega * omega * submatrix_I0(1,igrid_rs) &
                  + ( submatrix_I1k(1,igrid_rs) &
                      + 4.d0 * submatrix_I3k(1,igrid_rs) &
                      + 4.d0 * submatrix_I5k(1,igrid_rs) &
                    ) * factors_qkappa &
                  + ( submatrix_I1m(1,igrid_rs) &
                      + 4.d0 * submatrix_I3m(1,igrid_rs) &
                      + 4.d0 * submatrix_I5m(1,igrid_rs) &
                      + lsq2 * submatrix_I6(1,igrid_rs) &
                      - 4.d0 * submatrix_I7(1,igrid_rs) &
                    ) * factors_qmu ) * D1

  whole_vector_sph(idim_rs_sph+1,0) = whole_vector_sph(idim_rs_sph+1,0) &
              + ( - omega * omega * submatrix_I0(3,igrid_rs) &
                  + ( submatrix_I1k(3,igrid_rs) &
                      + 2.d0 * submatrix_I3k(2,igrid_rs) &
                      + 2.d0 * submatrix_I3k(3,igrid_rs) &
                      + 4.d0 * submatrix_I5k(3,igrid_rs) &
                    ) * factors_qkappa &
                  + ( submatrix_I1m(3,igrid_rs) &
                      + 2.d0 * submatrix_I3m(2,igrid_rs) &
                      + 2.d0 * submatrix_I3m(3,igrid_rs) &
                      + 4.d0 * submatrix_I5m(3,igrid_rs) &
                      + lsq2 * submatrix_I6(3,igrid_rs) &
                      - 4.d0 * submatrix_I7(3,igrid_rs) &
                    ) * factors_qmu ) * D1

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_wavefield( maxngrid_r,omega, &
            submatrix_I0,submatrix_I1k,submatrix_I1m, &
            submatrix_I2,submatrix_I3k,submatrix_I3m, &
            submatrix_I4,submatrix_I5k,submatrix_I5m, &
            submatrix_I6,submatrix_I7, &
            submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod, &
            grid_r,grid_mu,grid_qkappa,grid_qmu,l, &
            idim1_sph0,idim2_sph,idim1_tor0,idim2_tor, &
            idim0,init_npos_sph,init_npos_tor, &
            idim_rs_sph,idim_rs_tor, &
            idim_station_sph,idim_station_tor, &
            whole_matrix_sph,whole_matrix_tor, &
            whole_matrix_dr_sph,whole_matrix_dr_tor, &
            whole_vector_sph,whole_vector_tor,work_vector )

! compute wavefield for the given frequency

  implicit none

! variables for input/output
  integer maxngrid_r,l
  integer idim1_sph0,idim2_sph,idim1_tor0,idim2_tor
  integer idim0,init_npos_sph,init_npos_tor
  integer idim_rs_sph,idim_rs_tor
  integer idim_station_sph,idim_station_tor
  real(kind=8) submatrix_I0(4,maxngrid_r)
  real(kind=8) submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
  real(kind=8) submatrix_I2(4,maxngrid_r)
  real(kind=8) submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
  real(kind=8) submatrix_I4(4,maxngrid_r)
  real(kind=8) submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
  real(kind=8) submatrix_I6(4,maxngrid_r)
  real(kind=8) submatrix_I7(4,maxngrid_r)
  real(kind=8) submatrix_I3k_mod(6,maxngrid_r)
  real(kind=8) submatrix_I3m_mod(6,maxngrid_r)
  real(kind=8) submatrix_I4_mod(6,maxngrid_r)
  real(kind=8) grid_r(*),grid_mu(2,*)
  real(kind=8) grid_qkappa(*),grid_qmu(*)
  complex(kind=8) omega
  complex(kind=8) whole_matrix_sph(4,*),whole_matrix_tor(2,*)
  complex(kind=8) whole_matrix_dr_sph(*),whole_matrix_dr_tor(*)
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) whole_vector_tor(maxngrid_r,-2:2)
  complex(kind=8) work_vector(*)

! other variables
  integer idim1_sph,idim1_tor,ir,npos,m,itype_medium

!==============================================
! itype_medium=1: solid, itype_medium=0: liquid
!==============================================

  integer ndim_whole_matrix_sph,ndim_whole_matrix_tor,ier
  integer init_grid,end_grid,ns,nq
  real(kind=8) lsq,lsq2,eps
  complex(kind=8) anelastic_factor,factor_qkappa(2),factor_qmu(2)

  eps = -1.d0

  if ( l<=0 ) call error_handling(53)

! **********************************************************************
! initialize the whole_matrix
! **********************************************************************
  call init_complex_array( 8*maxngrid_r,whole_matrix_sph )
  call init_complex_array( 2*maxngrid_r,whole_matrix_tor )
  call init_complex_array( 2*maxngrid_r,whole_matrix_dr_sph )
  call init_complex_array(   maxngrid_r,whole_matrix_dr_tor )
  idim1_sph = max0( idim0,idim1_sph0 )
  idim1_tor = max0( idim0,idim1_tor0 )

! **********************************************************************
! **********************************************************************
! Spheroidal Component
! **********************************************************************
! **********************************************************************

! **********************************************************************
! constructing the whole_matrix
! **********************************************************************
  lsq2 = dble(l) * dble(l+1)
  lsq  = dsqrt( lsq2 )

  factor_qkappa(2) = anelastic_factor( dble(omega),grid_qkappa(idim1_sph) )
  factor_qmu(2) = anelastic_factor( dble(omega),grid_qmu(idim1_sph) )

!! DK DK it seems that this test is their way of seeing if that layer is fluid or solid
  if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph) == 0.d0 ) then
!! DK DK fluid case?
    npos = init_npos_sph
    itype_medium = 0
    whole_matrix_sph(4,npos) = omega * omega / factor_qkappa(2) &
              * cmplx( submatrix_I0(1,idim1_sph) ) &
              - cmplx( lsq2 * submatrix_I1k(1,idim1_sph) + submatrix_I2(1,idim1_sph) )
  else
!! DK DK solid case?
    npos = init_npos_sph
    itype_medium = 1
    whole_matrix_sph(4,npos) = omega * omega * cmplx( submatrix_I0(1,idim1_sph) ) &
              - factor_qkappa(2) * cmplx( submatrix_I1k(1,idim1_sph) &
                          + 4.d0 * submatrix_I3k(1,idim1_sph) &
                          + 4.d0 * submatrix_I5k(1,idim1_sph) ) &
              - factor_qmu(2) * cmplx( submatrix_I1m(1,idim1_sph) &
                          + 4.d0 * submatrix_I3m(1,idim1_sph) &
                          + 4.d0 * submatrix_I5m(1,idim1_sph) &
                          + lsq2 * submatrix_I6(1,idim1_sph) &
                          - 4.d0 * submatrix_I7(1,idim1_sph) )

    whole_matrix_sph(3,npos+1) = cmplx( lsq ) * ( factor_qkappa(2) &
                * cmplx( submatrix_I3k_mod(1,idim1_sph) + 2.d0 * submatrix_I5k(1,idim1_sph) ) &
              + factor_qmu(2) * cmplx( submatrix_I3m_mod(1,idim1_sph) &
                          - submatrix_I4_mod(1,idim1_sph) &
                          + 2.d0 * submatrix_I5m(1,idim1_sph) &
                          + submatrix_I6(1,idim1_sph) &
                          - 2.d0 * submatrix_I7(1,idim1_sph) ) )

    whole_matrix_sph(4,npos+1) = omega * omega * cmplx( submatrix_I0(1,idim1_sph) ) &
              - factor_qkappa(2) * cmplx( lsq2 * submatrix_I5k(1,idim1_sph) ) &
              - factor_qmu(2) * cmplx( submatrix_I2(1,idim1_sph) &
                          - 2.d0 * submatrix_I4(1,idim1_sph) &
                          + lsq2 * submatrix_I5m(1,idim1_sph) &
                          + submatrix_I6(1,idim1_sph) &
                          - 2.d0 * submatrix_I7(1,idim1_sph) )
  endif

  do ir=idim1_sph+1,idim2_sph-1

  factor_qkappa(1) = anelastic_factor( dble(omega),grid_qkappa(ir-1) )
  factor_qmu(1) = anelastic_factor( dble(omega),grid_qmu(ir-1) )
  factor_qkappa(2) = anelastic_factor( dble(omega),grid_qkappa(ir) )
  factor_qmu(2) = anelastic_factor( dble(omega),grid_qmu(ir) )

!! DK DK it seems that this test is their way of seeing if that layer is fluid or solid
  if ( grid_mu(1,ir)*grid_mu(2,ir)==0.d0 ) then

!! DK DK solid case?
    if ( itype_medium==1 ) then
      npos = npos + 2
      whole_matrix_sph(2,npos) = omega * omega * cmplx( submatrix_I0(2,ir-1) ) &
                - factor_qkappa(1) * cmplx( submatrix_I1k(2,ir-1) &
                        + 2.d0 * submatrix_I3k(2,ir-1) &
                        + 2.d0 * submatrix_I3k(3,ir-1) &
                        + 4.d0 * submatrix_I5k(2,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I1m(2,ir-1) &
                        + 2.d0 * submatrix_I3m(2,ir-1) &
                        + 2.d0 * submatrix_I3m(3,ir-1) &
                        + 4.d0 * submatrix_I5m(2,ir-1) &
                        + lsq2 * submatrix_I6(2,ir-1) &
                        - 4.d0 * submatrix_I7(2,ir-1) )

      whole_matrix_sph(3,npos) = whole_matrix_sph(3,npos) &
                + cmplx( lsq ) * ( factor_qkappa(1) &
                * cmplx( submatrix_I3k_mod(2,ir-1) &
                          + submatrix_I3k_mod(5,ir-1) &
                          + 2.d0 * submatrix_I5k(2,ir-1) ) &
                + factor_qmu(1) * cmplx( submatrix_I3m_mod(2,ir-1) &
                            + submatrix_I3m_mod(5,ir-1) &
                            - submatrix_I4_mod(2,ir-1) &
                            - submatrix_I4_mod(5,ir-1) &
                            + 2.d0 * submatrix_I5m(2,ir-1) &
                            + submatrix_I6(2,ir-1) &
                            - 2.d0 * submatrix_I7(2,ir-1) ) )

      whole_matrix_sph(4,npos) = omega * omega * cmplx( submatrix_I0(4,ir-1) ) &
                - factor_qkappa(1) * cmplx( submatrix_I1k(4,ir-1) &
                            + 4.d0 * submatrix_I3k(4,ir-1) &
                            + 4.d0 * submatrix_I5k(4,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I1m(4,ir-1) &
                            + 4.d0 * submatrix_I3m(4,ir-1) &
                            + 4.d0 * submatrix_I5m(4,ir-1) &
                            + lsq2 * submatrix_I6(4,ir-1) &
                            - 4.d0 * submatrix_I7(4,ir-1) )

      whole_matrix_sph(1,npos+1) = cmplx( lsq ) * ( factor_qkappa(1) &
                  * cmplx( submatrix_I3k_mod(3,ir-1) &
                            + 2.d0 * submatrix_I5k(2,ir-1) ) &
                  + factor_qmu(1) * cmplx( submatrix_I3m_mod(3,ir-1) &
                              - submatrix_I4_mod(3,ir-1) &
                              + 2.d0 * submatrix_I5m(2,ir-1) &
                              + submatrix_I6(2,ir-1) &
                              - 2.d0 * submatrix_I7(2,ir-1) ) )

      whole_matrix_sph(2,npos+1) = omega * omega * cmplx( submatrix_I0(2,ir-1) ) &
                - factor_qkappa(1) * cmplx( lsq2 * submatrix_I5k(2,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I2(2,ir-1) &
                            - submatrix_I4(2,ir-1) &
                            - submatrix_I4(3,ir-1) &
                            + lsq2 * submatrix_I5m(2,ir-1) &
                            + submatrix_I6(2,ir-1) &
                            - 2.d0 * submatrix_I7(2,ir-1) )

      whole_matrix_sph(3,npos+1) = cmplx( lsq ) * ( factor_qkappa(1) &
                  * cmplx( submatrix_I3k_mod(4,ir-1) &
                            + submatrix_I3k_mod(6,ir-1) &
                            + 2.d0 * submatrix_I5k(4,ir-1) ) &
                  + factor_qmu(1) * cmplx( submatrix_I3m_mod(4,ir-1) &
                              + submatrix_I3m_mod(6,ir-1) &
                              - submatrix_I4_mod(4,ir-1) &
                              - submatrix_I4_mod(6,ir-1) &
                              + 2.d0 * submatrix_I5m(4,ir-1) &
                              + submatrix_I6(4,ir-1) &
                              - 2.d0 * submatrix_I7(4,ir-1) ) )

      whole_matrix_sph(4,npos+1) = omega * omega * cmplx( submatrix_I0(4,ir-1) ) &
                  - factor_qkappa(1) * cmplx( lsq2 * submatrix_I5k(4,ir-1) ) &
                  - factor_qmu(1) * cmplx( submatrix_I2(4,ir-1) &
                              - 2.d0 * submatrix_I4(4,ir-1) &
                              + lsq2 * submatrix_I5m(4,ir-1) &
                              + submatrix_I6(4,ir-1) &
                              - 2.d0 * submatrix_I7(4,ir-1) )

      npos = npos + 2
      itype_medium = 0

      whole_matrix_sph(2,npos) = omega * grid_r(ir) * grid_r(ir)

      whole_matrix_sph(4,npos) = omega * omega / factor_qkappa(2) &
                * cmplx( submatrix_I0(1,ir) ) &
                - cmplx( lsq2 * submatrix_I1k(1,ir) &
                          + submatrix_I2(1,ir) )

!! DK DK fluid case?
    else

      npos = npos + 1
      itype_medium = 0

      whole_matrix_sph(3,npos) = omega * omega &
                * cmplx( submatrix_I0(2,ir-1) ) / factor_qkappa(1) &
                - cmplx( lsq2 * submatrix_I1k(2,ir-1) + submatrix_I2(2,ir-1) )

      whole_matrix_sph(4,npos) = omega * omega &
                * ( cmplx( submatrix_I0(4,ir-1) ) / factor_qkappa(1) &
                    + cmplx( submatrix_I0(1,ir) ) / factor_qkappa(2) ) &
                - cmplx( lsq2 * submatrix_I1k(4,ir-1) + submatrix_I2(4,ir-1) ) &
                - cmplx( lsq2 * submatrix_I1k(1,ir) + submatrix_I2(1,ir) )
    endif

  else

    if ( itype_medium==0 ) then
      npos = npos + 1

      whole_matrix_sph(3,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(2,ir-1) ) &
                - cmplx( lsq2 * submatrix_I1k(2,ir-1) + submatrix_I2(2,ir-1) )

      whole_matrix_sph(4,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(4,ir-1) ) &
                - cmplx( lsq2 * submatrix_I1k(4,ir-1) + submatrix_I2(4,ir-1) )

      npos = npos + 1
      itype_medium = 1

      whole_matrix_sph(3,npos) = - omega * grid_r(ir) * grid_r(ir)

      whole_matrix_sph(4,npos) = omega * omega &
                * cmplx( submatrix_I0(1,ir) ) &
                - factor_qkappa(2) * cmplx( submatrix_I1k(1,ir) + 4.d0 * submatrix_I3k(1,ir) + 4.d0 * submatrix_I5k(1,ir) ) &
                - factor_qmu(2) * cmplx( submatrix_I1m(1,ir) &
                          + 4.d0 * submatrix_I3m(1,ir) &
                          + 4.d0 * submatrix_I5m(1,ir) &
                          + lsq2 * submatrix_I6(1,ir) &
                          - 4.d0 * submatrix_I7(1,ir) )

      whole_matrix_sph(3,npos+1) = cmplx( lsq ) * ( factor_qkappa(2) &
                * cmplx( submatrix_I3k_mod(1,ir) + 2.d0 * submatrix_I5k(1,ir) ) &
                + factor_qmu(2) * cmplx( submatrix_I3m_mod(1,ir) &
                          - submatrix_I4_mod(1,ir) &
                          + 2.d0 * submatrix_I5m(1,ir) &
                          + submatrix_I6(1,ir) &
                          - 2.d0 * submatrix_I7(1,ir) ) )

      whole_matrix_sph(4,npos+1) = omega * omega &
              * cmplx( submatrix_I0(1,ir) ) &
              - factor_qkappa(2) * cmplx( lsq2 * submatrix_I5k(1,ir) ) &
              - factor_qmu(2) * cmplx( submatrix_I2(1,ir) &
                          - 2.d0 * submatrix_I4(1,ir) &
                          + lsq2 * submatrix_I5m(1,ir) &
                          + submatrix_I6(1,ir) &
                          - 2.d0 * submatrix_I7(1,ir) )

    else

      npos = npos + 2
      itype_medium = 1

      whole_matrix_sph(2,npos) = omega * omega &
                * cmplx( submatrix_I0(2,ir-1) ) &
                - factor_qkappa(1) * cmplx( submatrix_I1k(2,ir-1) &
                          + 2.d0 * submatrix_I3k(2,ir-1) &
                          + 2.d0 * submatrix_I3k(3,ir-1) &
                          + 4.d0 * submatrix_I5k(2,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I1m(2,ir-1) &
                          + 2.d0 * submatrix_I3m(2,ir-1) &
                          + 2.d0 * submatrix_I3m(3,ir-1) &
                          + 4.d0 * submatrix_I5m(2,ir-1) &
                          + lsq2 * submatrix_I6(2,ir-1) &
                          - 4.d0 * submatrix_I7(2,ir-1) )

      whole_matrix_sph(3,npos) = whole_matrix_sph(3,npos) &
              + cmplx( lsq ) * ( factor_qkappa(1) * cmplx( submatrix_I3k_mod(2,ir-1) &
                          + 2.d0 * submatrix_I5k(2,ir-1) ) &
                + factor_qmu(1) * cmplx( submatrix_I3m_mod(2,ir-1) &
                          - submatrix_I4_mod(2,ir-1) &
                          + 2.d0 * submatrix_I5m(2,ir-1) &
                          + submatrix_I6(2,ir-1) &
                          - 2.d0 * submatrix_I7(2,ir-1) ) )

      whole_matrix_sph(4,npos) = omega * omega &
                * cmplx( submatrix_I0(4,ir-1) + submatrix_I0(1,ir) ) &
                - factor_qkappa(1) * cmplx( submatrix_I1k(4,ir-1) &
                          + 4.d0 * submatrix_I3k(4,ir-1) &
                          + 4.d0 * submatrix_I5k(4,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I1m(4,ir-1) &
                          + 4.d0 * submatrix_I3m(4,ir-1) &
                          + 4.d0 * submatrix_I5m(4,ir-1) &
                          + lsq2 * submatrix_I6(4,ir-1) &
                          - 4.d0 * submatrix_I7(4,ir-1) ) &
                - factor_qkappa(2) * cmplx( submatrix_I1k(1,ir) &
                          + 4.d0 * submatrix_I3k(1,ir) &
                          + 4.d0 * submatrix_I5k(1,ir) ) &
                - factor_qmu(2) * cmplx( submatrix_I1m(1,ir) &
                          + 4.d0 * submatrix_I3m(1,ir) &
                          + 4.d0 * submatrix_I5m(1,ir) &
                          + lsq2 * submatrix_I6(1,ir) &
                          - 4.d0 * submatrix_I7(1,ir) )

      whole_matrix_sph(1,npos+1) = cmplx( lsq ) * ( &
                factor_qkappa(1) * cmplx( submatrix_I3k_mod(3,ir-1) &
                          + 2.d0 * submatrix_I5k(2,ir-1) ) &
                + factor_qmu(1) * cmplx( submatrix_I3m_mod(3,ir-1) &
                          - submatrix_I4_mod(3,ir-1) &
                          + 2.d0 * submatrix_I5m(2,ir-1) &
                          + submatrix_I6(2,ir-1) &
                          - 2.d0 * submatrix_I7(2,ir-1) ) )

      whole_matrix_sph(2,npos+1) = omega * omega &
                * cmplx( submatrix_I0(2,ir-1) ) &
                - factor_qkappa(1) * cmplx( lsq2 * submatrix_I5k(2,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I2(2,ir-1) &
                          - submatrix_I4(2,ir-1) &
                          - submatrix_I4(3,ir-1) &
                          + lsq2 * submatrix_I5m(2,ir-1) &
                          + submatrix_I6(2,ir-1) &
                          - 2.d0 * submatrix_I7(2,ir-1) )

      whole_matrix_sph(3,npos+1) = cmplx( lsq ) * ( factor_qkappa(1) &
                * cmplx( submatrix_I3k_mod(4,ir-1) &
                          + 2.d0 * submatrix_I5k(4,ir-1) ) &
                + factor_qmu(1) * cmplx( submatrix_I3m_mod(4,ir-1) &
                          - submatrix_I4_mod(4,ir-1) &
                          + 2.d0 * submatrix_I5m(4,ir-1) &
                          + submatrix_I6(4,ir-1) &
                          - 2.d0 * submatrix_I7(4,ir-1) ) &
                + factor_qkappa(2) * cmplx( submatrix_I3k_mod(1,ir) &
                          + 2.d0 * submatrix_I5k(1,ir) ) &
                + factor_qmu(2) * cmplx( submatrix_I3m_mod(1,ir) &
                          - submatrix_I4_mod(1,ir) &
                          + 2.d0 * submatrix_I5m(1,ir) &
                          + submatrix_I6(1,ir) &
                          - 2.d0 * submatrix_I7(1,ir) ) )

      whole_matrix_sph(4,npos+1) = omega * omega &
                * cmplx( submatrix_I0(4,ir-1) &
                        + submatrix_I0(1,ir) ) &
                - factor_qkappa(1) &
                * cmplx( lsq2 * submatrix_I5k(4,ir-1) ) &
                - factor_qmu(1) &
                * cmplx( submatrix_I2(4,ir-1) &
                          - 2.d0 * submatrix_I4(4,ir-1) &
                          + lsq2 * submatrix_I5m(4,ir-1) &
                          + submatrix_I6(4,ir-1) &
                          - 2.d0 * submatrix_I7(4,ir-1) ) &
                - factor_qkappa(2) &
                * cmplx( lsq2 * submatrix_I5k(1,ir) ) &
                - factor_qmu(2) &
                * cmplx( submatrix_I2(1,ir) &
                          - 2.d0 * submatrix_I4(1,ir) &
                          + lsq2 * submatrix_I5m(1,ir) &
                          + submatrix_I6(1,ir) &
                          - 2.d0 * submatrix_I7(1,ir) )

      whole_matrix_sph(1,npos+2) = cmplx( lsq ) * ( factor_qkappa(1) &
                * cmplx( submatrix_I3k_mod(5,ir-1) ) &
                + factor_qmu(1) * cmplx( submatrix_I3m_mod(5,ir-1) - submatrix_I4_mod(5,ir-1) ) )

      whole_matrix_sph(3,npos+2) = cmplx( lsq ) * ( factor_qkappa(1) * cmplx( submatrix_I3k_mod(6,ir-1) ) &
                + factor_qmu(1) * cmplx( submatrix_I3m_mod(6,ir-1) - submatrix_I4_mod(6,ir-1) ) )

    endif
  endif
  enddo ! of do ir=idim1_sph+1,idim2_sph-1

  factor_qkappa(1) = anelastic_factor( dble(omega),grid_qkappa(idim2_sph-1) )
  factor_qmu(1) = anelastic_factor( dble(omega),grid_qmu(idim2_sph-1) )

  if ( itype_medium==0 ) then

    npos = npos + 1

    whole_matrix_sph(3,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(2,idim2_sph-1) ) &
              - cmplx( lsq2 * submatrix_I1k(2,idim2_sph-1) &
                          + submatrix_I2(2,idim2_sph-1) )

    whole_matrix_sph(4,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(4,idim2_sph-1) ) &
              - cmplx( lsq2 * submatrix_I1k(4,idim2_sph-1) &
                          + submatrix_I2(4,idim2_sph-1) )

    ndim_whole_matrix_sph = npos

  else

    npos = npos + 2

    whole_matrix_sph(2,npos) = omega * omega &
              * cmplx( submatrix_I0(2,idim2_sph-1) ) &
            - factor_qkappa(1) &
              * cmplx( submatrix_I1k(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3k(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3k(3,idim2_sph-1) &
                        + 4.d0 * submatrix_I5k(2,idim2_sph-1) ) &
            - factor_qmu(1) &
              * cmplx( submatrix_I1m(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3m(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3m(3,idim2_sph-1) &
                        + 4.d0 * submatrix_I5m(2,idim2_sph-1) &
                        + lsq2 * submatrix_I6(2,idim2_sph-1) &
                        - 4.d0 * submatrix_I7(2,idim2_sph-1) )

    whole_matrix_sph(3,npos) = whole_matrix_sph(3,npos) &
            + cmplx( lsq ) * ( factor_qkappa(1) &
              * cmplx( submatrix_I3k_mod(2,idim2_sph-1) &
                        + submatrix_I3k_mod(5,idim2_sph-1) &
                        + 2.d0 * submatrix_I5k(2,idim2_sph-1) ) &
            + factor_qmu(1) &
              * cmplx( submatrix_I3m_mod(2,idim2_sph-1) &
                        + submatrix_I3m_mod(5,idim2_sph-1) &
                        - submatrix_I4_mod(2,idim2_sph-1) &
                        - submatrix_I4_mod(5,idim2_sph-1) &
                        + 2.d0 * submatrix_I5m(2,idim2_sph-1) &
                        + submatrix_I6(2,idim2_sph-1) &
                        - 2.d0 * submatrix_I7(2,idim2_sph-1) ) )

    whole_matrix_sph(4,npos) = omega * omega &
              * cmplx( submatrix_I0(4,idim2_sph-1) ) &
            - factor_qkappa(1) &
              * cmplx( submatrix_I1k(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I3k(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I5k(4,idim2_sph-1) ) &
            - factor_qmu(1) &
              * cmplx( submatrix_I1m(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I3m(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I5m(4,idim2_sph-1) &
                        + lsq2 * submatrix_I6(4,idim2_sph-1) &
                        - 4.d0 * submatrix_I7(4,idim2_sph-1) )

    whole_matrix_sph(1,npos+1) = cmplx( lsq ) * ( &
              factor_qkappa(1) &
              * cmplx( submatrix_I3k_mod(3,idim2_sph-1) &
                        + 2.d0 * submatrix_I5k(2,idim2_sph-1) &
                      ) &
            + factor_qmu(1) &
              * cmplx( submatrix_I3m_mod(3,idim2_sph-1) &
                        - submatrix_I4_mod(3,idim2_sph-1) &
                        + 2.d0 * submatrix_I5m(2,idim2_sph-1) &
                        + submatrix_I6(2,idim2_sph-1) &
                        - 2.d0 * submatrix_I7(2,idim2_sph-1) ) )

    whole_matrix_sph(2,npos+1) = omega * omega &
            * cmplx( submatrix_I0(2,idim2_sph-1) ) &
            - factor_qkappa(1) &
              * cmplx( lsq2 * submatrix_I5k(2,idim2_sph-1) ) &
            - factor_qmu(1) &
              * cmplx( submatrix_I2(2,idim2_sph-1) &
                        - submatrix_I4(2,idim2_sph-1) &
                        - submatrix_I4(3,idim2_sph-1) &
                        + lsq2 * submatrix_I5m(2,idim2_sph-1) &
                        + submatrix_I6(2,idim2_sph-1) &
                        - 2.d0 * submatrix_I7(2,idim2_sph-1) )

    whole_matrix_sph(3,npos+1) = cmplx( lsq ) * ( &
              factor_qkappa(1) &
              * cmplx( submatrix_I3k_mod(4,idim2_sph-1) &
                        + submatrix_I3k_mod(6,idim2_sph-1) &
                        + 2.d0 * submatrix_I5k(4,idim2_sph-1) ) &
            + factor_qmu(1) &
              * cmplx( submatrix_I3m_mod(4,idim2_sph-1) &
                        + submatrix_I3m_mod(6,idim2_sph-1) &
                        - submatrix_I4_mod(4,idim2_sph-1) &
                        - submatrix_I4_mod(6,idim2_sph-1) &
                        + 2.d0 * submatrix_I5m(4,idim2_sph-1) &
                        + submatrix_I6(4,idim2_sph-1) &
                        - 2.d0 * submatrix_I7(4,idim2_sph-1) ) )

    whole_matrix_sph(4,npos+1) = omega * omega &
            * cmplx( submatrix_I0(4,idim2_sph-1) ) &
            - factor_qkappa(1) &
              * cmplx( lsq2 * submatrix_I5k(4,idim2_sph-1) ) &
            - factor_qmu(1) &
              * cmplx( submatrix_I2(4,idim2_sph-1) &
                        - 2.d0 * submatrix_I4(4,idim2_sph-1) &
                        + lsq2 * submatrix_I5m(4,idim2_sph-1) &
                        + submatrix_I6(4,idim2_sph-1) &
                        - 2.d0 * submatrix_I7(4,idim2_sph-1) )

    ndim_whole_matrix_sph = npos+1

  endif

! **********************************************************************
! compute the wavefield
! **********************************************************************
! imposing fixed boundary conditions at r=0 and free surface boundary
! conditions at the surface

  if ( ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph)/=0.d0 ) .and.( grid_r(idim1_sph)==0.d0 ) ) then
    init_grid = max0(init_npos_sph,3)
  else
    init_grid = init_npos_sph
  endif

  if ( grid_mu(1,idim2_sph-1)*grid_mu(2,idim2_sph-1)==0.d0 ) then
    end_grid = ndim_whole_matrix_sph - 1
  else
    end_grid = ndim_whole_matrix_sph
  endif

  m = max0(-l,-2)
  ns = idim_rs_sph - init_grid + 1

  if ( mod(l,100)==0 ) then
    nq = end_grid - init_grid + 1
  else
    nq = min0( end_grid - idim_station_sph + 1, end_grid - init_grid + 1 )
  endif

  call dclisb(whole_matrix_sph(1,init_grid), &
                    end_grid - init_grid + 1, &
                    3,4,ns,nq,whole_vector_sph(init_grid,m), &
                    eps,whole_matrix_dr_sph(init_grid), &
                    work_vector(init_grid),ier)

  do m=max0(-l,-2)+1,min0(l,2)
    call dcsbsub(whole_matrix_sph(1,init_grid), &
                       end_grid - init_grid + 1, &
                       3,4,ns,nq,whole_vector_sph(init_grid,m), &
                       whole_matrix_dr_sph(init_grid), &
                       work_vector(init_grid))
  enddo

! **********************************************************************
! **********************************************************************
! Toroidal Component
! **********************************************************************
! **********************************************************************

! **********************************************************************
! constructing the whole_matrix
! **********************************************************************

  factor_qmu(2) = anelastic_factor( dble(omega),grid_qmu(idim1_tor) )
  npos = init_npos_tor
  whole_matrix_tor(2,npos) = omega * omega &
              * cmplx( submatrix_I0(1,idim1_tor) ) &
              - factor_qmu(2) &
                * cmplx( submatrix_I2(1,idim1_tor) &
                          - 2.d0 * submatrix_I4(1,idim1_tor) &
                          + submatrix_I6(1,idim1_tor) &
                          + ( lsq2 - 2.d0 ) * submatrix_I7(1,idim1_tor) )

  do ir=idim1_tor+1,idim2_tor-1

  factor_qmu(1) = anelastic_factor( dble(omega),grid_qmu(ir-1) )
  factor_qmu(2) = anelastic_factor( dble(omega),grid_qmu(ir) )

  npos = npos + 1

  whole_matrix_tor(1,npos) = omega * omega &
                * cmplx(  submatrix_I0(2,ir-1) ) &
                - factor_qmu(1) &
                  * cmplx( submatrix_I2(2,ir-1) &
                              - submatrix_I4(2,ir-1) &
                              - submatrix_I4(3,ir-1) &
                              + submatrix_I6(2,ir-1) &
                              + ( lsq2 - 2.d0 ) * submatrix_I7(2,ir-1) )

  whole_matrix_tor(2,npos) = omega * omega &
                * cmplx(   submatrix_I0(4,ir-1) &
                          + submatrix_I0(1,ir) ) &
                - factor_qmu(1) &
                  * cmplx( submatrix_I2(4,ir-1) &
                              - 2.d0 * submatrix_I4(4,ir-1) &
                              + submatrix_I6(4,ir-1) &
                              + ( lsq2 - 2.d0 ) * submatrix_I7(4,ir-1) ) &
                - factor_qmu(2) &
                  * cmplx( submatrix_I2(1,ir) &
                              - 2.d0 * submatrix_I4(1,ir) &
                              + submatrix_I6(1,ir) &
                              + ( lsq2 - 2.d0 ) * submatrix_I7(1,ir) )

  enddo

  factor_qmu(1) = anelastic_factor( dble(omega), grid_qmu(idim2_tor-1) )

  npos = npos + 1

  whole_matrix_tor(1,npos) = omega * omega &
              * cmplx(   submatrix_I0(2,idim2_tor-1) ) &
              - factor_qmu(1) &
                * cmplx( submatrix_I2(2,idim2_tor-1) &
                          - submatrix_I4(2,idim2_tor-1) &
                          - submatrix_I4(3,idim2_tor-1) &
                          + submatrix_I6(2,idim2_tor-1) &
                          + ( lsq2 - 2.d0 ) * submatrix_I7(2,idim2_tor-1) )

  whole_matrix_tor(2,npos) = omega * omega &
              * cmplx(  submatrix_I0(4,idim2_tor-1) ) &
              - factor_qmu(1) &
                * cmplx( submatrix_I2(4,idim2_tor-1) &
                          - 2.d0 * submatrix_I4(4,idim2_tor-1) &
                          + submatrix_I6(4,idim2_tor-1) &
                          + ( lsq2 - 2.d0 ) &
                            * submatrix_I7(4,idim2_tor-1) )

  ndim_whole_matrix_tor = npos

! **********************************************************************
! compute the wavefield
! **********************************************************************

  m = max0(-l,-2)
  init_grid = init_npos_tor
  end_grid = ndim_whole_matrix_tor
  ns = idim_rs_tor - init_grid + 1

  if ( mod(l,100)==0 ) then
    nq = end_grid - init_grid + 1
  else
    nq = min0( end_grid - idim_station_tor + 1, end_grid - init_grid + 1 )
  endif

  call dclisb(whole_matrix_tor(1,init_grid), end_grid - init_grid + 1, &
                    1,2,ns,nq,whole_vector_tor(init_grid,m), &
                    eps,whole_matrix_dr_tor(init_grid), &
                    work_vector(init_grid),ier)

  do m=max0(-l,-2)+1,min0(l,2)
    call dcsbsub(whole_matrix_tor(1,init_grid), &
                       end_grid - init_grid + 1, &
                       1,2,ns,nq,whole_vector_tor(init_grid,m), &
                       whole_matrix_dr_tor(init_grid), &
                       work_vector(init_grid))
  enddo

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_wavefield0( maxngrid_r,omega, &
            submatrix_I0,submatrix_I1k,submatrix_I1m, &
            submatrix_I2,submatrix_I3k,submatrix_I3m, &
            submatrix_I5k,submatrix_I5m, &
            submatrix_I6,submatrix_I7, &
            grid_r,grid_mu,grid_qkappa,grid_qmu,l, &
            idim1_sph,idim2_sph, &
            idim0,init_npos_sph,init_npos_tor, &
            idim_rs_sph, &
            idim_station_sph, &
            whole_matrix_sph, &
            whole_matrix_dr_sph, &
            whole_vector_sph,work_vector )

! compute wavefield for the given frequency

  implicit none

! variables for input/output
  integer maxngrid_r,l
  integer idim1_sph,idim2_sph
  integer idim0,init_npos_sph,init_npos_tor
  integer idim_rs_sph
  integer idim_station_sph
  real(kind=8) submatrix_I0(4,maxngrid_r)
  real(kind=8) submatrix_I1k(4,maxngrid_r),submatrix_I1m(4,maxngrid_r)
  real(kind=8) submatrix_I2(4,maxngrid_r)
  real(kind=8) submatrix_I3k(4,maxngrid_r),submatrix_I3m(4,maxngrid_r)
  real(kind=8) submatrix_I5k(4,maxngrid_r),submatrix_I5m(4,maxngrid_r)
  real(kind=8) submatrix_I6(4,maxngrid_r)
  real(kind=8) submatrix_I7(4,maxngrid_r)
  real(kind=8) grid_r(*),grid_mu(2,*)
  real(kind=8) grid_qkappa(*),grid_qmu(*)
  complex(kind=8) omega
  complex(kind=8) whole_matrix_sph(4,*)
  complex(kind=8) whole_matrix_dr_sph(*)
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) work_vector(*)

! other variables
  integer ir,npos,m,itype_medium

!==============================================
! itype_medium=1: solid, itype_medium=0: liquid
!==============================================

  integer ndim_whole_matrix_sph,ndim_whole_matrix_tor,ier
  integer init_grid,end_grid,ns,nq
  real(kind=8) lsq,lsq2,eps
  complex(kind=8) anelastic_factor,factor_qkappa(2),factor_qmu(2)

  eps = -1.d0

  if ( l /= 0 ) call error_handling(54)

! **********************************************************************
! initialize the whole_matrix
! **********************************************************************
  call init_complex_array( 8*maxngrid_r,whole_matrix_sph )
  call init_complex_array( 2*maxngrid_r,whole_matrix_dr_sph )

  idim0 = 1
  init_npos_sph = 1
  init_npos_tor = 1

! **********************************************************************
! **********************************************************************
! Spheroidal Component
! **********************************************************************
! **********************************************************************

! **********************************************************************
! constructing the whole_matrix
! **********************************************************************
  lsq2 = dble(l) * dble(l+1)
  lsq  = dsqrt( lsq2 )

  factor_qkappa(2) = anelastic_factor( dble(omega),grid_qkappa(idim1_sph) )
  factor_qmu(2) = anelastic_factor( dble(omega),grid_qmu(idim1_sph) )

  if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph)==0.d0 ) then
    npos = 1
    itype_medium = 0

    whole_matrix_sph(2,npos) = omega * omega / factor_qkappa(2) &
              * cmplx( submatrix_I0(1,idim1_sph) ) &
            - cmplx( lsq2 * submatrix_I1k(1,idim1_sph) &
                          + submatrix_I2(1,idim1_sph) )
  else
    npos = 1
    itype_medium = 1

    whole_matrix_sph(2,npos) = omega * omega &
              * cmplx( submatrix_I0(1,idim1_sph) ) &
              - factor_qkappa(2) &
                * cmplx( submatrix_I1k(1,idim1_sph) &
                          + 4.d0 * submatrix_I3k(1,idim1_sph) &
                          + 4.d0 * submatrix_I5k(1,idim1_sph) ) &
              - factor_qmu(2) &
                * cmplx( submatrix_I1m(1,idim1_sph) &
                          + 4.d0 * submatrix_I3m(1,idim1_sph) &
                          + 4.d0 * submatrix_I5m(1,idim1_sph) &
                          + lsq2 * submatrix_I6(1,idim1_sph) &
                          - 4.d0 * submatrix_I7(1,idim1_sph) )
  endif

  do ir=idim1_sph+1,idim2_sph-1

  factor_qkappa(1) = anelastic_factor( dble(omega),grid_qkappa(ir-1) )
  factor_qmu(1) = anelastic_factor( dble(omega),grid_qmu(ir-1) )
  factor_qkappa(2) = anelastic_factor( dble(omega),grid_qkappa(ir) )
  factor_qmu(2) = anelastic_factor( dble(omega),grid_qmu(ir) )

  if ( grid_mu(1,ir)*grid_mu(2,ir) == 0.d0 ) then

    if ( itype_medium == 1 ) then

      npos = npos + 1

      whole_matrix_sph(1,npos) = omega * omega &
                * cmplx( submatrix_I0(2,ir-1) ) &
                - factor_qkappa(1) &
                  * cmplx( submatrix_I1k(2,ir-1) &
                        + 2.d0 * submatrix_I3k(2,ir-1) &
                        + 2.d0 * submatrix_I3k(3,ir-1) &
                        + 4.d0 * submatrix_I5k(2,ir-1) ) &
                - factor_qmu(1) &
                  * cmplx( submatrix_I1m(2,ir-1) &
                        + 2.d0 * submatrix_I3m(2,ir-1) &
                        + 2.d0 * submatrix_I3m(3,ir-1) &
                        + 4.d0 * submatrix_I5m(2,ir-1) &
                        + lsq2 * submatrix_I6(2,ir-1) &
                        - 4.d0 * submatrix_I7(2,ir-1) )

      whole_matrix_sph(2,npos) = omega * omega &
                * cmplx( submatrix_I0(4,ir-1) ) &
                - factor_qkappa(1) &
                  * cmplx( submatrix_I1k(4,ir-1) &
                            + 4.d0 * submatrix_I3k(4,ir-1) &
                            + 4.d0 * submatrix_I5k(4,ir-1) ) &
                - factor_qmu(1) &
                  * cmplx( submatrix_I1m(4,ir-1) &
                            + 4.d0 * submatrix_I3m(4,ir-1) &
                            + 4.d0 * submatrix_I5m(4,ir-1) &
                            + lsq2 * submatrix_I6(4,ir-1) &
                            - 4.d0 * submatrix_I7(4,ir-1) )

      npos = npos + 1
      itype_medium = 0

      whole_matrix_sph(1,npos) = omega * grid_r(ir) * grid_r(ir)

      whole_matrix_sph(2,npos) = omega * omega * ( cmplx( submatrix_I0(1,ir)   ) / factor_qkappa(2) ) &
              - cmplx( lsq2 * submatrix_I1k(1,ir) + submatrix_I2(1,ir) )
    else

      npos = npos + 1

      whole_matrix_sph(1,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(2,ir-1) ) - cmplx( lsq2 * submatrix_I1k(2,ir-1) + submatrix_I2(2,ir-1) )

      whole_matrix_sph(2,npos) = omega * omega * ( cmplx( submatrix_I0(4,ir-1) ) / factor_qkappa(1) &
                  + cmplx( submatrix_I0(1,ir)   ) / factor_qkappa(2) ) &
              - cmplx( lsq2 * submatrix_I1k(4,ir-1) + submatrix_I2(4,ir-1) ) &
              - cmplx( lsq2 * submatrix_I1k(1,ir) + submatrix_I2(1,ir) )
    endif

  else

    if ( itype_medium==0 ) then

      npos = npos + 1

      whole_matrix_sph(1,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(2,ir-1) ) - cmplx( lsq2 * submatrix_I1k(2,ir-1) + submatrix_I2(2,ir-1) )

      whole_matrix_sph(2,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(4,ir-1) ) - cmplx( lsq2 * submatrix_I1k(4,ir-1) + submatrix_I2(4,ir-1) )

      npos = npos + 1
      itype_medium = 1

      whole_matrix_sph(1,npos) = - omega * grid_r(ir) * grid_r(ir)

      whole_matrix_sph(2,npos) = omega * omega * cmplx( submatrix_I0(1,ir) ) &
                - factor_qkappa(2) * cmplx( submatrix_I1k(1,ir) &
                          + 4.d0 * submatrix_I3k(1,ir) &
                          + 4.d0 * submatrix_I5k(1,ir) ) &
                - factor_qmu(2) * cmplx( submatrix_I1m(1,ir) &
                          + 4.d0 * submatrix_I3m(1,ir) &
                          + 4.d0 * submatrix_I5m(1,ir) &
                          + lsq2 * submatrix_I6(1,ir) &
                          - 4.d0 * submatrix_I7(1,ir) )

    else

      npos = npos + 1
      itype_medium = 1

      whole_matrix_sph(1,npos) = omega * omega * cmplx( submatrix_I0(2,ir-1) ) &
                - factor_qkappa(1) * cmplx( submatrix_I1k(2,ir-1) &
                          + 2.d0 * submatrix_I3k(2,ir-1) &
                          + 2.d0 * submatrix_I3k(3,ir-1) &
                          + 4.d0 * submatrix_I5k(2,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I1m(2,ir-1) &
                          + 2.d0 * submatrix_I3m(2,ir-1) &
                          + 2.d0 * submatrix_I3m(3,ir-1) &
                          + 4.d0 * submatrix_I5m(2,ir-1) &
                          + lsq2 * submatrix_I6(2,ir-1) &
                          - 4.d0 * submatrix_I7(2,ir-1) )

      whole_matrix_sph(2,npos) = omega * omega * cmplx( submatrix_I0(4,ir-1) &
                          + submatrix_I0(1,ir) ) &
                - factor_qkappa(1) * cmplx( submatrix_I1k(4,ir-1) &
                          + 4.d0 * submatrix_I3k(4,ir-1) &
                          + 4.d0 * submatrix_I5k(4,ir-1) ) &
                - factor_qmu(1) * cmplx( submatrix_I1m(4,ir-1) &
                          + 4.d0 * submatrix_I3m(4,ir-1) &
                          + 4.d0 * submatrix_I5m(4,ir-1) &
                          + lsq2 * submatrix_I6(4,ir-1) &
                          - 4.d0 * submatrix_I7(4,ir-1) ) &
                - factor_qkappa(2) * cmplx( submatrix_I1k(1,ir) &
                          + 4.d0 * submatrix_I3k(1,ir) &
                          + 4.d0 * submatrix_I5k(1,ir) ) &
                - factor_qmu(2) * cmplx( submatrix_I1m(1,ir) &
                          + 4.d0 * submatrix_I3m(1,ir) &
                          + 4.d0 * submatrix_I5m(1,ir) &
                          + lsq2 * submatrix_I6(1,ir) &
                          - 4.d0 * submatrix_I7(1,ir) )

    endif
  endif
  enddo ! of do ir=idim1_sph+1,idim2_sph-1

  factor_qkappa(1) = anelastic_factor( dble(omega),grid_qkappa(idim2_sph-1) )
  factor_qmu(1) = anelastic_factor( dble(omega),grid_qmu(idim2_sph-1) )

  if ( itype_medium==0 ) then

    npos = npos + 1

    whole_matrix_sph(1,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(2,idim2_sph-1) ) &
              - cmplx( lsq2 * submatrix_I1k(2,idim2_sph-1) + submatrix_I2(2,idim2_sph-1) )

    whole_matrix_sph(2,npos) = omega * omega / factor_qkappa(1) &
                * cmplx( submatrix_I0(4,idim2_sph-1) ) &
              - cmplx( lsq2 * submatrix_I1k(4,idim2_sph-1) + submatrix_I2(4,idim2_sph-1) )

    ndim_whole_matrix_sph = npos

  else

    npos = npos + 1

    whole_matrix_sph(1,npos) = omega * omega * cmplx( submatrix_I0(2,idim2_sph-1) ) &
            - factor_qkappa(1) * cmplx( submatrix_I1k(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3k(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3k(3,idim2_sph-1) &
                        + 4.d0 * submatrix_I5k(2,idim2_sph-1) ) &
            - factor_qmu(1) * cmplx( submatrix_I1m(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3m(2,idim2_sph-1) &
                        + 2.d0 * submatrix_I3m(3,idim2_sph-1) &
                        + 4.d0 * submatrix_I5m(2,idim2_sph-1) &
                        + lsq2 * submatrix_I6(2,idim2_sph-1) &
                        - 4.d0 * submatrix_I7(2,idim2_sph-1) )

    whole_matrix_sph(2,npos) = omega * omega * cmplx( submatrix_I0(4,idim2_sph-1) ) &
            - factor_qkappa(1) * cmplx( submatrix_I1k(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I3k(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I5k(4,idim2_sph-1) ) &
            - factor_qmu(1) * cmplx( submatrix_I1m(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I3m(4,idim2_sph-1) &
                        + 4.d0 * submatrix_I5m(4,idim2_sph-1) &
                        + lsq2 * submatrix_I6(4,idim2_sph-1) &
                        - 4.d0 * submatrix_I7(4,idim2_sph-1) )

    ndim_whole_matrix_sph = npos

  endif

! **********************************************************************
! compute the wavefield
! **********************************************************************

! imposing fixed boundary conditions at r=0 and free surface boundary conditions at the surface
  if ( grid_mu(1,idim1_sph)*grid_mu(2,idim1_sph) /= 0.d0 .and. grid_r(idim1_sph) == 0.d0 ) then
    init_grid = 2
  else
    init_grid = 1
  endif

  if ( grid_mu(1,idim2_sph-1)*grid_mu(2,idim2_sph-1) == 0.d0 ) then
    end_grid = ndim_whole_matrix_sph - 1
  else
    end_grid = ndim_whole_matrix_sph
  endif

  m = 0
  ns = idim_rs_sph - init_grid + 1

  if ( mod(l,100) == 0 ) then
    nq = end_grid - init_grid + 1
  else
    nq = end_grid - idim_station_sph + 1
  endif

  call dclisb(whole_matrix_sph(1,init_grid), end_grid - init_grid + 1, &
                    1,4,ns,nq,whole_vector_sph(init_grid,m), &
                    eps,whole_matrix_dr_sph(init_grid), &
                    work_vector(init_grid),ier)

  ndim_whole_matrix_tor = 0

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=8) function anelastic_factor( omega,qmu )

! compute the anelastic factor (ratio between complex mu and
! real mu) for the given quality factor, qmu.

  implicit none

! input/output variables
  real(kind=8) omega,qmu

! other variables
  real(kind=8) vr,vi

! constants
  real(kind=8), parameter :: pi=3.1415926535897932d0

! **********************************************************************
! compute the complex velocity for the given qmu
! **********************************************************************
  if ( omega == 0.d0 .or. qmu < 0.d0 ) then
!! DK DK no attenuation
    vr = 1.d0
    vi = 0.d0
  else
    vr = 1.d0 + dlog( omega / ( 2.d0 * pi ) ) / ( pi * qmu )
    vi = 1.d0 / ( 2.d0 * qmu )
  endif

! **********************************************************************
! compute the anelastic factor
! **********************************************************************
  anelastic_factor = cmplx( vr, vi ) * cmplx( vr, vi )

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dclisb(a, n, nud, n1, np, nq, b, eps, dr, z, ier)

!***********************************************************************
!  simultaneous linear equations with real symmetric positive definite *
!      band matrix by the Cholesky method.                             *
!  parameters                                                          *
!    (1) a : 2-dim. array containing the matrix.                       *
!    (2) n : order of the matrix.                                      *
!    (3) nud : size of band's half width.                              *
!    (4) n1 : row size of the array a in the 'dimension' statement.    *
!    (5) b : 1-dim. array containing the right hand side vector.       *
!    (6) eps : parameter to check singurarity off the matrix           *
!              standard value = 1.0d-14                                *
!    (7) dr : 1-dim. working array.                                    *
!    (8) z : 1-dim. working array.                                     *
!    (9) ier : error code.                                             *
!  copyright   T. Oguni   July 30, 1989   version 1.0                  *
!***********************************************************************

!! DK DK this solver apparently works for any "nud" value of the half bandwidth of the matrix,
!! DK DK thus currently they use it to solve a tridiagonal system (i.e. nud = 1 I suppose) but
!! DK DK we could easily use it for pentadiagonal (nud = 2?), heptadiagonal (nud = 3?) or anything
!! DK DK of higher order.

  implicit none

  integer n, nud, n1, np, nq, ier
  complex(kind=8) a(n1,n), b(n), dr(n), z(n)
  real(kind=8) eps
  complex(kind=8) xx, s, sumval, au, t
  real(kind=8) eps1
  integer i ,m, j, k1, mj, k

! check the input data
  ier = 0
  eps1 = 1.0d-14
  m = nud + 1

  if (n <= 0 .or. nud <= 0 .or. n1 < m) then
    ier = 2
    write(*,*) '(subr. lisb) invalid argument. ', n, nud, n1
    return
  endif

  if (eps <= 0.0) eps = eps1

! modified Cholesky decomposition
  j = 1

  if (abs(a(m,1)) <= eps) then
    ier = 1
    write(*,*) '(subr. lisb) singular at step # ', j
    return
  endif

  dr(1) = cmplx(1.0d0) / a(m,1)
  xx = a(m-1,2)
  a(m-1,2) = a(m-1,2) * dr(1)
  s = a(m,2) - xx * a(m-1,2)
  j = 2

  if (abs(s) <= eps) then
    ier = 1
    write(*,*) '(subr. lisb) singular at step # ', j
    return
  endif

  dr(2) = cmplx(1.0d0) / s

  if (m < 3) then

  do j=3,n

    xx = a(1,j)
    a(1,j) = xx * dr(j-1)
    s = a(2,j) - xx * a(1,j)

    if (abs(s) <= eps) then
      ier = 1
      write(*,*) ' (subr. lisb) singular at step # ', j
      return
    endif

    dr(j) = cmplx(1.0d0) / s
  enddo

  else

 do j=3,n

  k1 = 1
  if (j >= m) k1 = j - m + 1
  mj = m - j

  do i=k1+1,j-1
    sumval = cmplx(0.0d0)
    do k=k1,i-1
      sumval = sumval + a(m-i+k,i) * a(mj+k,j)
    enddo
    a(mj+i,j) = a(mj+i,j) - sumval
  enddo

  sumval = cmplx(0.0d0)

  do i=k1,j-1
    xx = a(mj+i,j)
    au = xx * dr(i)
    sumval = sumval + xx *au
    a(mj+i,j) = au
  enddo

  t = a(m,j) - sumval

  if (abs(t) <= eps) then
    ier = 1
    write(*,*) ' (subr. lisb) singular at step # ', j
    return
  endif

  dr(j) = cmplx(1.0d0) / t

   enddo
  endif

! perform the substitution for Cholesky
  call dcsbsub(a, n, nud, n1, np, nq, b, dr, z)

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! perform the substitution for Cholesky
  subroutine dcsbsub(a, n, nud, n1, np, nq, b, dr, z)

  implicit none

  integer n, nud, n1, np, nq
  complex(kind=8) a(n1,n), b(n), dr(n), z(n)
  complex(kind=8) sumval
  integer m, j, i1, k, j1

! forward substitution
  m = nud + 1

  if (m < 3) then

  z(np) = b(np)

  do j=np+1,n
    z(j) = b(j) - a(1,j) * z(j-1)
  enddo

  do j=1,np-1
    z(j) = cmplx(0.d0)
  enddo

  do j=np,n
    z(j) = z(j) * dr(j)
  enddo

  b(n) = z(n)

  do j=1,nq-1
    b(n-j) = z(n-j) - a(1,n-j+1) * b(n-j+1)
  enddo

  else

  z(np) = b(np)
  z(np+1) = b(np+1) - a(m-1,np+1) * z(np)

  do j=np+2,n
    if (j > np-1+m) then
      i1 = 1
    else
      i1 = np-1+m - j + 1
    endif
    sumval = cmplx(0.0d0)
    do k=i1,m-1
      sumval = sumval + a(k,j) * z(j-m+k)
    enddo
    z(j) = b(j) - sumval
  enddo

  do j=1,np-1
    z(j) = cmplx(0.d0)
  enddo

  do j=np,n
    z(j) = z(j) * dr(j)
  enddo

  b(n) = z(n)
  b(n-1) = z(n-1) - a(m-1,n) * z(n)

  do j=3,nq
    j1 = n - j + 1
    i1 = 1
    if (j < m) i1 = m - j + 1
    sumval = cmplx(0.0d0)
    do k=i1,m-1
      sumval = sumval + a(k,m-k+j1) * b(m-k+j1)
    enddo
    b(j1) = z(j1) - sumval
  enddo

  endif

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_complex_array( n,A )

! initializing the accumulated displacement at the station

  implicit none

! variables for input/output
  integer n
  complex(kind=8) A(n)

  A(:) = cmplx(0.d0)

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_displacement_station( maxngrid_r,maxlmax,whole_vector_tor,whole_vector_sph, &
              l,n_station,idim_station_sph,idim_station_tor, &
              vecsph_sph1,vecsph_sph2,vecsph_tor, &
              station_displacement )

! accumulating the displacement at the station

  implicit none

! variables for input/output
  integer maxngrid_r,maxlmax,l,n_station
  integer idim_station_sph,idim_station_tor
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) whole_vector_tor(maxngrid_r,-2:2)
  complex(kind=8) vecsph_sph1(3,0:maxlmax,-2:2,*)
  complex(kind=8) vecsph_sph2(3,0:maxlmax,-2:2,*)
  complex(kind=8) vecsph_tor(3,0:maxlmax,-2:2,*)
  complex(kind=8) station_displacement(3,*)

! other variables
  integer m,i_station,icomp
  real(kind=8) lsq

  if ( l <= 0 ) call error_handling(55)

  lsq = dsqrt( dble(l) * dble(l+1) )

  do i_station = 1,n_station

! ---- compute the value of the trial functions at the station
  do m = max0(-l,-2),min0(l,2)

! -------- horizontal dependent part
!         call comp_vecsph(l,m,
!     &                          station_theta(i_station),
!     &                          station_phi(i_station),
!     &                          vecsph_sph1,vecsph_sph2,vecsph_tor)

! ---- compute the displacement at the station
    do icomp = 1,3
      station_displacement(icomp,i_station) &
              = station_displacement(icomp,i_station) &
                + whole_vector_sph(idim_station_sph,m) &
                  * vecsph_sph1(icomp,l,m,i_station) &
                + whole_vector_sph(idim_station_sph+1,m) &
                  * vecsph_sph2(icomp,l,m,i_station) / cmplx( lsq ) &
                + whole_vector_tor(idim_station_tor,m) &
                  * vecsph_tor(icomp,l,m,i_station) / cmplx( lsq )
      enddo
    enddo
  enddo

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_displacement_station0( maxngrid_r,maxlmax,whole_vector_sph, &
              l,n_station,idim_station_sph,vecsph_sph1,station_displacement )

! accumulating the displacement at the station

  implicit none

! variables for input/output
  integer maxngrid_r,maxlmax,l,n_station
  integer idim_station_sph
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) vecsph_sph1(3,0:maxlmax,-2:2,*)
  complex(kind=8) station_displacement(3,*)

! other variables
  integer m,i_station,icomp
  real(kind=8) lsq

  if ( l /= 0 ) call error_handling(56)

  lsq = dsqrt( dble(l) * dble(l+1) )

  do i_station = 1,n_station

! ---- compute the value of the trial functions at the station
  do m = max0(-l,-2),min0(l,2)

! -------- horizontal dependent part
!         call comp_vecsph(l,m,
!     &                          station_theta(i_station),
!     &                          station_phi(i_station),
!     &                          vecsph_sph1,vecsph_sph2,vecsph_tor)

! ---- compute the displacement at the station
    do icomp = 1,3
      station_displacement(icomp,i_station) = station_displacement(icomp,i_station) &
                + whole_vector_sph(idim_station_sph,m) * vecsph_sph1(icomp,l,m,i_station)
      enddo
    enddo
  enddo

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_vecsph_station( maxlmax,lmax,n_station,station_theta,station_phi, &
            vecsph_sph1,vecsph_sph2,vecsph_tor )

  implicit none

! variables for input/output
  integer maxlmax,lmax,n_station
  real(kind=8) station_theta(*),station_phi(*)
  complex(kind=8) vecsph_sph1(3,0:maxlmax,-2:2,*)
  complex(kind=8) vecsph_sph2(3,0:maxlmax,-2:2,*)
  complex(kind=8) vecsph_tor(3,0:maxlmax,-2:2,*)

! other variables
  integer i_station

  if ( lmax>maxlmax ) call error_handling(40)

  do i_station=1,n_station
    call  comp_vecsph_all( maxlmax,lmax, station_theta(i_station),station_phi(i_station), &
            vecsph_sph1(1,0,-2,i_station), vecsph_sph2(1,0,-2,i_station), vecsph_tor(1,0,-2,i_station) )
  enddo

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_vecsph_all(maxlmax,lmax,theta,phi,vecsph_sph1,vecsph_sph2,vecsph_tor)

! compute the vector spherical harmonics

  implicit none

! variables for input/output
  integer maxlmax,lmax
  real(kind=8) theta,phi
  complex(kind=8) vecsph_sph1(3,0:maxlmax,-2:2)
  complex(kind=8) vecsph_sph2(3,0:maxlmax,-2:2)
  complex(kind=8) vecsph_tor(3,0:maxlmax,-2:2)

! other variables
  integer i,l,m,m0
  real(kind=8) factor,x,compute_Legendre_polynomial
  complex(kind=8) expimp

! constant
  real(kind=8), parameter :: pi=3.1415926535897932d0

  do l = 0,lmax
    do m = max0(-l,-2),min0(l,2)

! **********************************************************************
! checking the arguments
! **********************************************************************

  m0 = iabs(m)
  if ( l < 0 .or. m0 > l ) call error_handling(41)

! **********************************************************************
! compute the normalization factor (including the sign)
! **********************************************************************

  factor = 1.d0

  do i = l-m0+1,l+m0
    factor = factor * dble(i)
  enddo

  factor = dsqrt( dble(2*l+1)/(4.d0*pi) / factor )
  if ( m0 /= m .and. mod(m0,2) == 1 ) factor = - factor

! **********************************************************************
! compute each component of the vector spherical harmonics
! **********************************************************************

  x = dcos(theta)
  expimp = exp( cmplx( 0.d0, dble(m)*phi ) )

  vecsph_sph1(1,l,m) = factor * compute_Legendre_polynomial(l,m0,dcos(theta)) * expimp
  vecsph_sph1(2,l,m) = 0.d0
  vecsph_sph1(3,l,m) = 0.d0

  vecsph_sph2(1,l,m) = 0.d0

  if ( l >= m0+1 ) then
    vecsph_sph2(2,l,m) = factor * (   dble(m0) * x / dsin(theta) * compute_Legendre_polynomial(l,m0,x) &
                                  + compute_Legendre_polynomial(l,m0+1,x) ) * expimp
  else
    vecsph_sph2(2,l,m) = factor * (   dble(m0) * x / dsin(theta) &
                                     * compute_Legendre_polynomial(l,m0,x) ) * expimp
  endif

  vecsph_sph2(3,l,m) = factor * cmplx( 0.d0, dble(m) ) / dsin(theta) &
                             * compute_Legendre_polynomial(l,m0,dcos(theta)) * expimp

  vecsph_tor(1,l,m) = cmplx(0.d0)
  vecsph_tor(2,l,m) = vecsph_sph2(3,l,m)
  vecsph_tor(3,l,m) = - vecsph_sph2(2,l,m)

    enddo ! of do  m=max0(-l,-2),min0(l,2)
  enddo ! of do l=0,lmax

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function compute_Legendre_polynomial(vall,valm,valx)

! compute Legendre polynominal P(vall,valm) at point valx

  implicit none

! input parameters
  integer, intent(in) :: vall,valm
  real(kind=8), intent(in) :: valx

! local variables
  integer indexval,lother
  real(kind=8) local_factor,partial_result,partial_m,partial_m_one,square_root_of_product

! error if a parameter is not correct
  if(valm < 0 .or. valm > vall .or. abs(valx) > 1.d0) call error_handling(42)

  partial_m = 1.d0

  if(valm > 0) then
    square_root_of_product = dsqrt((1.d0 - valx) * (1.d0 + valx))
    local_factor = 1.d0

    do indexval = 1, valm
      partial_m = - partial_m * local_factor * square_root_of_product
      local_factor = local_factor + 2.d0
    enddo
  endif

  if(vall == valm) then

    compute_Legendre_polynomial = partial_m

  else

    partial_m_one = valx * dble(2*valm + 1) * partial_m
    if(vall == valm + 1) then
      compute_Legendre_polynomial = partial_m_one
    else
      do lother = valm + 2, vall
        partial_result = (valx * dble(2*lother - 1) * partial_m_one - dble(lother + valm - 1) * partial_m) / dble(lother - valm)
        partial_m = partial_m_one
        partial_m_one = partial_result
      enddo
      compute_Legendre_polynomial = partial_result
    endif

  endif

  end

