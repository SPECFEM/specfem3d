      program dsmti
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computation of synthetic seismograms for spherical symmetric
c TI media by using the Direct Solution Method.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c constant
      integer maxngrid_r
      integer max_nstation,maxn_structure_zone,maxnfreq
      parameter ( maxngrid_r = 100000 )
ccccccccc DK DK      parameter ( max_nstation = 100 )
      parameter ( max_nstation = 20 )
ccccccccc DK DK      parameter ( maxn_structure_zone = 15 )
      parameter ( maxn_structure_zone = 7 )
      parameter ( maxnfreq = 8192 )
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c variables for time series
      integer n_frequency,i_frequency
      real*8 time_series_length,omega_imag
c variables for numerical grids
      integer ngrid_r0,lmin0,lmax0
      integer ngrid_r_estimated,ngrid_r,lmin,lmax,l
      integer n_frequency_band,i_frequency_band
      integer i_frequency_min,i_frequency_max
      real*8 grid_r(maxngrid_r)
      real*8 grid_rho(2,maxngrid_r)
      real*8 grid_Ak(2,maxngrid_r),grid_Am(2,maxngrid_r)
      real*8 grid_Ck(2,maxngrid_r),grid_Cm(2,maxngrid_r)
      real*8 grid_L(2,maxngrid_r),grid_N(2,maxngrid_r)
      real*8 grid_Fk(2,maxngrid_r),grid_Fm(2,maxngrid_r)
      real*8 grid_kappa(2,maxngrid_r)
      real*8 grid_mu(2,maxngrid_r)
      real*8 grid_qkappa(maxngrid_r)
      real*8 grid_qmu(maxngrid_r)
c variables for the structure
      integer n_structure_zone
      integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
      real*8 rmin_structure_zone(maxn_structure_zone)
      real*8 rmax_structure_zone(maxn_structure_zone)
      real*8 rho_structure_zone(4,maxn_structure_zone)
      real*8 vpv_structure_zone(4,maxn_structure_zone)
      real*8 vph_structure_zone(4,maxn_structure_zone)
      real*8 vsv_structure_zone(4,maxn_structure_zone)
      real*8 vsh_structure_zone(4,maxn_structure_zone)
      real*8 eta_structure_zone(4,maxn_structure_zone)
      real*8 qkappa_structure_zone(maxn_structure_zone)
      real*8 qmu_structure_zone(maxn_structure_zone)
c variables for a source
      integer igrid_rs
      integer idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0
      real*8 source_r,source_mt(3,3)
      real*8 source_depth,source_lat,source_lon
      real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
      real*8 grid_qkappas,grid_qmus
c variables for stations
      integer n_station
      integer idim_station_sph,idim_station_tor,
     &              idim_station_sph0,idim_station_tor0
      real*8 station_lat(max_nstation),station_lon(max_nstation)
      real*8 station_theta(max_nstation),station_phi(max_nstation)
      complex*16 station_displacement(3,max_nstation,maxnfreq)
      character*80 sac_file(max_nstation)
c variables for matrix elements
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
c variables for wavefield
      integer i_significance,idim0,init_npos_sph,init_npos_tor
      real*8 amp_max
      complex*16 omega
      complex*16 whole_matrix_sph(4,2*maxngrid_r)
      complex*16 whole_matrix_tor(2,maxngrid_r)
      complex*16 whole_matrix_dr_sph(2*maxngrid_r)
      complex*16 whole_matrix_dr_tor(maxngrid_r)
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
      complex*16 work_vector(2*maxngrid_r)
      real work_time(32*maxnfreq)
      complex*16 work_spc(16*maxnfreq)
      equivalence ( work_vector,work_spc )
c MPI variables
      include 'mpicom.h'
      integer info,irank,nsnd
c functions
      omega(time_series_length,i_frequency,omega_imag)
     &     = dcmplx( 2.d0*pi*dble(i_frequency)/dble(time_series_length),
     &      -omega_imag )
      ngrid_r_estimated
     &        (i_frequency_band,n_frequency_band,ngrid_r0)
     &        = int(
     &            dnint (
     &              ngrid_r0
     &                * dble(i_frequency_band) / dble(n_frequency_band)
     &            )
     &          )
      lmin
     &        (i_frequency_band,n_frequency_band,lmin0 )
     &    = lmin0
      lmax
     &        (i_frequency_band,n_frequency_band,lmax0 )
     &        = int(
     &            dnint (
     &              lmax0
     &                * dble(i_frequency_band) / dble(n_frequency_band)
     &            )
     &          )
      i_frequency_min
     &        (i_frequency_band,n_frequency_band,n_frequency)
     &        = int(
     &            dnint (
     &              n_frequency
     &            * dble(i_frequency_band-1) / dble(n_frequency_band)
     &            )
     &          ) + 1 - int(1.d0/dble(i_frequency_band))
      i_frequency_max
     &        (i_frequency_band,n_frequency_band,n_frequency)
     &        = int(
     &            dnint (
     &              n_frequency
     &                * dble(i_frequency_band) / dble(n_frequency_band)
     &            )
     &          )
c
      data n_frequency_band / 4 /
c **********************************************************************
c inputing parameters
c **********************************************************************
c initialize MPI
      call minit
c
      if ( rank.eq.0 ) then
        call param_input
     &       ( maxngrid_r,max_nstation,maxn_structure_zone,
     &         time_series_length,n_frequency,omega_imag,
     &         ngrid_r0,lmin0,lmax0,
     &         n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &         rho_structure_zone,vpv_structure_zone,vph_structure_zone,
     &         vsv_structure_zone,vsh_structure_zone,eta_structure_zone,
     &         qkappa_structure_zone,qmu_structure_zone,
     &         source_r,source_mt,source_depth,source_lat,source_lon,
     &         n_station,station_lat,station_lon,
     &         station_theta,station_phi,sac_file )
        call msndi
        call mputi( n_frequency,1,info )
        call mputi( ngrid_r0,1,info )
        call mputi( lmin0,1,info )
        call mputi( lmax0,1,info )
        call mputi( n_structure_zone,1,info )
        call mputi( n_station,1,info )
        call mputd( time_series_length,1,info )
        call mputd( omega_imag,1,info )
        call mputd( rmin_structure_zone,maxn_structure_zone,info )
        call mputd( rmax_structure_zone,maxn_structure_zone,info )
        call mputd( rho_structure_zone,4*maxn_structure_zone,info )
        call mputd( vpv_structure_zone,4*maxn_structure_zone,info )
        call mputd( vph_structure_zone,4*maxn_structure_zone,info )
        call mputd( vsv_structure_zone,4*maxn_structure_zone,info )
        call mputd( vsh_structure_zone,4*maxn_structure_zone,info )
        call mputd( eta_structure_zone,4*maxn_structure_zone,info )
        call mputd( qkappa_structure_zone,maxn_structure_zone,info )
        call mputd( qmu_structure_zone,maxn_structure_zone,info )
        call mputd( source_r,1,info )
        call mputd( source_mt,9,info )
        call mputd( station_theta,max_nstation,info )
        call mputd( station_phi,max_nstation,info )
        call msend( -1,10,info )
      else
        call mrecv( 10,0,info )
        call mgeti( n_frequency,1,info )
        call mgeti( ngrid_r0,1,info )
        call mgeti( lmin0,1,info )
        call mgeti( lmax0,1,info )
        call mgeti( n_structure_zone,1,info )
        call mgeti( n_station,1,info )
        call mgetd( time_series_length,1,info )
        call mgetd( omega_imag,1,info )
        call mgetd( rmin_structure_zone,maxn_structure_zone,info )
        call mgetd( rmax_structure_zone,maxn_structure_zone,info )
        call mgetd( rho_structure_zone,4*maxn_structure_zone,info )
        call mgetd( vpv_structure_zone,4*maxn_structure_zone,info )
        call mgetd( vph_structure_zone,4*maxn_structure_zone,info )
        call mgetd( vsv_structure_zone,4*maxn_structure_zone,info )
        call mgetd( vsh_structure_zone,4*maxn_structure_zone,info )
        call mgetd( eta_structure_zone,4*maxn_structure_zone,info )
        call mgetd( qkappa_structure_zone,maxn_structure_zone,info )
        call mgetd( qmu_structure_zone,maxn_structure_zone,info )
        call mgetd( source_r,1,info )
        call mgetd( source_mt,9,info )
        call mgetd( station_theta,max_nstation,info )
        call mgetd( station_phi,max_nstation,info )
      endif
      write(*,*) rank,'mpi interface setting is OK'
      do 1000 i_frequency_band=1,n_frequency_band
c **********************************************************************
c generating numerical grids
c **********************************************************************
        call grid_generation
     &        ( maxngrid_r,ngrid_r_estimated
     &        (i_frequency_band,n_frequency_band,ngrid_r0),
     &         n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &         vsv_structure_zone,vsh_structure_zone,
     &         grid_r,source_r,ngrid_r )
        call assign_structure
     &       ( ngrid_r,grid_r,
     &         n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &         rho_structure_zone,
     &         vpv_structure_zone,vph_structure_zone,
     &         vsv_structure_zone,vsh_structure_zone,
     &         eta_structure_zone,
     &         qkappa_structure_zone,qmu_structure_zone,
     &         grid_rho,grid_Ak,grid_Am,grid_Ck,grid_Cm,
     &         grid_L,grid_N,grid_Fk,grid_Fm,
     &         grid_kappa,grid_mu,grid_qkappa,grid_qmu,
     &         idim1_sph,idim2_sph,idim1_tor,idim2_tor )
        call assign_source
     &          ( ngrid_r,grid_r,source_r,
     &            grid_L,grid_Ck,grid_Cm,grid_Fk,grid_Fm,
     &            grid_qkappa,grid_qmu,grid_mu,
     &            igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &            grid_qkappas,grid_qmus,
     &            idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0 )
        call assign_station
     &          ( ngrid_r,grid_r,grid_mu,
     &            idim_station_sph,idim_station_tor,
     &            idim_station_sph0,idim_station_tor0 )
       write(*,*) rank,'assign parameter running is over'
       write(*,*) rank,'number of grid:',ngrid_r
c **********************************************************************
c computing submatrices for the elastic part of the medium
c **********************************************************************
        call comp_submatrix
     &          ( ngrid_r,grid_r,
     &            grid_rho,grid_kappa,grid_mu,grid_Ak,grid_Am,
     &            grid_Ck,grid_Cm,grid_L,grid_N,grid_Fk,grid_Fm,
     &            submatrix_I0,submatrix_I1k,submatrix_I1m,submatrix_I2,
     &            submatrix_I3k,submatrix_I3m,submatrix_I4,submatrix_I5k,
     &            submatrix_I5m,submatrix_I6,submatrix_I7 )
       write(*,*) rank,'comp_submatrix function is calling over'
        call comp_submatrix_mod
     &          ( ngrid_r,grid_r,
     &            grid_rho,grid_kappa,grid_mu,
     &            grid_Ak,grid_Am,grid_L,grid_N,grid_Fk,grid_Fm,
     &            submatrix_I0,submatrix_I1k,submatrix_I3k,
     &            submatrix_I3m,submatrix_I4,submatrix_I5k,
     &            submatrix_I5m,submatrix_I6,submatrix_I7,
     &            submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod )
c **********************************************************************
c **********************************************************************
c **********************************************************************
c computing wavefield for each frequency (for complex omega)
c **********************************************************************
c **********************************************************************
c **********************************************************************
        do 200 i_frequency=i_frequency_min(i_frequency_band,
     &                                           n_frequency_band,
     &                                           n_frequency)
     &                           +rank,
     &                           i_frequency_max(i_frequency_band,
     &                                           n_frequency_band,
     &                                           n_frequency),
     &                           size
c         write(6,*) i_frequency,i_frequency_band,
c     &                    omega(time_series_length,i_frequency,omega_imag)
cc DK DK
       print *,'ifreq, imin, imax = ',i_frequency,
     &                           i_frequency_min(i_frequency_band,
     &                                           n_frequency_band,
     &                                           n_frequency),
     &                           i_frequency_max(i_frequency_band,
     &                                           n_frequency_band,
     &                                           n_frequency)
cc DK DK
         if ( i_frequency.ne.0 ) then
          if ( i_frequency.gt.maxnfreq ) stop 'Too large i_freqnecy.'
          call init_complex_array
     &         ( 3*max_nstation,station_displacement(1,1,i_frequency) )
c
          do 110 l=lmin(i_frequency_band,n_frequency_band,lmin0 ),
     &                   lmax(i_frequency_band,n_frequency_band,lmax0 )
c **********************************************************************
c computing the displacement and the traction at the boundary
c of the heterogeneous region.
c **********************************************************************
            if ( l.eq.0 ) then
              call comp_excitation0
     &           ( maxngrid_r,
     &             omega(time_series_length,i_frequency,omega_imag),
     &             ngrid_r,grid_r,l,source_r,source_mt,
     &             igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &             grid_qkappas,grid_qmus,
     &             submatrix_I0,submatrix_I1k,submatrix_I1m,
     &             submatrix_I2,submatrix_I3k,submatrix_I3m,
     &             submatrix_I4,submatrix_I5k,submatrix_I5m,
     &             submatrix_I6,submatrix_I7,
     &             submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &             idim_rs_sph0,idim_rs_tor0,
     &             whole_vector_sph,whole_vector_tor )
              call comp_wavefield0
     &           ( maxngrid_r,
     &             omega(time_series_length,i_frequency,omega_imag),
     &             submatrix_I0,submatrix_I1k,submatrix_I1m,
     &             submatrix_I2,submatrix_I3k,submatrix_I3m,
     &             submatrix_I4,submatrix_I5k,submatrix_I5m,
     &             submatrix_I6,submatrix_I7,
     &             submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &             ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &             idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &             idim0,init_npos_sph,init_npos_tor,
     &             idim_rs_sph0,idim_rs_tor0,
     &             idim_station_sph0,idim_station_tor0,
     &             whole_matrix_sph,whole_matrix_tor,
     &             whole_matrix_dr_sph,whole_matrix_dr_tor,
     &             whole_vector_sph,whole_vector_tor,work_vector )
              call check_amp_significance
     &              ( maxngrid_r,l,amp_max,i_significance,idim0,
     &                init_npos_sph,init_npos_tor,
     &                idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &                grid_mu,idim_station_sph0,idim_station_tor0,
     &                whole_vector_sph,whole_vector_tor,work_vector )
              call comp_displacement_station0
     &              ( maxngrid_r,whole_vector_tor,whole_vector_sph,
     &                l,n_station,station_theta,station_phi,
     &                idim_station_sph0,idim_station_tor0,
     &                station_displacement(1,1,i_frequency) )
            else
              if ( i_significance.eq.1 ) then
              call comp_excitation
     &           ( maxngrid_r,
     &             omega(time_series_length,i_frequency,omega_imag),
     &             ngrid_r,grid_r,l,source_r,source_mt,
     &             igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &             grid_qkappas,grid_qmus,
     &             submatrix_I0,submatrix_I1k,submatrix_I1m,
     &             submatrix_I2,submatrix_I3k,submatrix_I3m,
     &             submatrix_I4,submatrix_I5k,submatrix_I5m,
     &             submatrix_I6,submatrix_I7,
     &             submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &             idim_rs_sph,idim_rs_tor,
     &             whole_vector_sph,whole_vector_tor )
              call comp_wavefield
     &           ( maxngrid_r,
     &             omega(time_series_length,i_frequency,omega_imag),
     &             submatrix_I0,submatrix_I1k,submatrix_I1m,
     &             submatrix_I2,submatrix_I3k,submatrix_I3m,
     &             submatrix_I4,submatrix_I5k,submatrix_I5m,
     &             submatrix_I6,submatrix_I7,
     &             submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,
     &             ngrid_r,grid_r,grid_mu,grid_qkappa,grid_qmu,l,
     &             idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &             idim0,init_npos_sph,init_npos_tor,
     &             idim_rs_sph,idim_rs_tor,
     &             idim_station_sph,idim_station_tor,
     &             whole_matrix_sph,whole_matrix_tor,
     &             whole_matrix_dr_sph,whole_matrix_dr_tor,
     &             whole_vector_sph,whole_vector_tor,work_vector )
              call check_amp_significance
     &              ( maxngrid_r,l,amp_max,i_significance,idim0,
     &                init_npos_sph,init_npos_tor,
     &                idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &                grid_mu,idim_station_sph,idim_station_tor,
     &                whole_vector_sph,whole_vector_tor,work_vector )
              call comp_displacement_station
     &              ( maxngrid_r,whole_vector_tor,whole_vector_sph,
     &                l,n_station,station_theta,station_phi,
     &                idim_station_sph,idim_station_tor,
     &                station_displacement(1,1,i_frequency) )
              endif
            endif
  110     continue
         endif
c **********************************************************************
c writing the result to the spectrum files
c **********************************************************************
c        call write_sac_file
c     &           ( i_frequency,n_station,sac_file,
c     &             station_displacement )
  200   continue
c
 1000 continue
c
        nsnd = 2 * 3 * max_nstation
        if ( rank.ne.0 ) then
          call msndi
          do 1030 i_frequency_band=1,n_frequency_band
          do 1020 i_frequency=i_frequency_min(i_frequency_band,
     &                                             n_frequency_band,
     &                                             n_frequency)
     &                             +rank,
     &                             i_frequency_max(i_frequency_band,
     &                                             n_frequency_band,
     &                                             n_frequency),
     &                             size
            call mputd( station_displacement(1,1,i_frequency),
     &                        nsnd,info )
 1020     continue
 1030     continue
          call msend( 0,20,info )
        else
          do 1060 irank=1,size-1
            call mrecv( 20,irank,info )
            do 1050 i_frequency_band=1,n_frequency_band
            do 1040 i_frequency=i_frequency_min(i_frequency_band,
     &                                             n_frequency_band,
     &                                             n_frequency)
     &                             +irank,
     &                             i_frequency_max(i_frequency_band,
     &                                             n_frequency_band,
     &                                             n_frequency),
     &                             size
              call mgetd( station_displacement(1,1,i_frequency),
     &                          nsnd,info )
 1040       continue
 1050       continue
 1060     continue
        endif
 1100 continue
      if ( rank.eq.0 ) then
        call convspc
     &     ( max_nstation,maxnfreq,n_station,
     &       time_series_length,n_frequency,omega_imag,
     &       station_displacement,
     &       source_depth,source_lat,source_lon,station_lat,station_lon,
     &       work_spc,work_time,sac_file )
      endif
c
      call mpi_finalize( info )
c
      end
