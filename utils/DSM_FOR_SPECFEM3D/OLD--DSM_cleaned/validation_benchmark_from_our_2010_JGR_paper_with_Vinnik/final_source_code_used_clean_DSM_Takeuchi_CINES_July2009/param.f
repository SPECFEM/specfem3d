cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine param_input
     &        (maxngrid_r,max_nstation,maxn_structure_zone,
     &         time_series_length,n_frequency,omega_imag,
     &         ngrid_r,lmin,lmax,
     &         n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &         rho_structure_zone,vpv_structure_zone,vph_structure_zone,
     &         vsv_structure_zone,vsh_structure_zone,eta_structure_zone,
     &         qkappa_structure_zone,qmu_structure_zone,
     &         source_r,source_mt,source_depth,source_lat,source_lon,
     &         n_station,station_lat,station_lon,
     &         station_theta,station_phi,sac_file )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c inputting parameters
c    required subroutines: error_handling,distaz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,max_nstation,maxn_structure_zone
      integer n_frequency,ngrid_r,lmin,lmax
      integer n_structure_zone,n_station
      real*8 time_series_length,omega_imag
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
      real*8 source_r,source_mt(3,3)
      real*8 source_depth,source_lat,source_lon
      real*8 station_lat(max_nstation),station_lon(max_nstation)
      real*8 station_theta(max_nstation),station_phi(max_nstation)
      character*80 sac_file(max_nstation)
c other variables
      integer i,nexp,nerr
      real dist,az,baz,xdeg
      character*80 dummy,tmpfile
c constants
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      data tmpfile / 'work' /
c
c **********************************************************************
c reading input file from the standard input and writing out to the
c temporary file
c **********************************************************************
c opening a temporary file
      open( unit=11, file=tmpfile, status='unknown' )
c writing to the temporary file
  100 continue
        read(5,110) dummy
  110   format(a80)
        if ( dummy(1:1).eq.'c' ) goto 100
        if ( dummy(1:3).eq.'end' ) goto 120
        write(11,110) dummy
        goto 100
  120 continue
c closing the temporary file
      close(11)
c
c **********************************************************************
c reading the parameter from the temporary file
c **********************************************************************
c opening the temporary file
      open( unit=11, file=tmpfile, status='unknown' )
c reading the parameters
c ---- parameters for time series ---
      read(11,*) time_series_length,n_frequency
      read(11,*) omega_imag
c ---- parameters for numerical grids ----
      read(11,*) ngrid_r,lmin,lmax
      if ( ngrid_r.gt.maxngrid_r )
     &        call error_handling(1)
c ---- parameters for structure ---
      read(11,*) n_structure_zone
      if ( n_structure_zone.gt.maxn_structure_zone )
     &        call error_handling(2)
      do 130 i=1,n_structure_zone
        read(11,*) rmin_structure_zone(i),rmax_structure_zone(i),
     &                  rho_structure_zone(1,i),rho_structure_zone(2,i),
     &                   rho_structure_zone(3,i),rho_structure_zone(4,i)
        read(11,*) vpv_structure_zone(1,i),vpv_structure_zone(2,i),
     &                   vpv_structure_zone(3,i),vpv_structure_zone(4,i)
        read(11,*) vph_structure_zone(1,i),vph_structure_zone(2,i),
     &                   vph_structure_zone(3,i),vph_structure_zone(4,i)
        read(11,*) vsv_structure_zone(1,i),vsv_structure_zone(2,i),
     &                   vsv_structure_zone(3,i),vsv_structure_zone(4,i)
        read(11,*) vsh_structure_zone(1,i),vsh_structure_zone(2,i),
     &                   vsh_structure_zone(3,i),vsh_structure_zone(4,i)
        read(11,*) eta_structure_zone(1,i),eta_structure_zone(2,i),
     &                  eta_structure_zone(3,i),eta_structure_zone(4,i),
     &                   qmu_structure_zone(i),qkappa_structure_zone(i)
  130 continue
c ---- parameters for a source ---
      read(11,*) source_depth,source_lat,source_lon
      read(11,*) nexp,
     &                 source_mt(1,1),source_mt(2,2),source_mt(3,3),
     &                 source_mt(1,2),source_mt(1,3),source_mt(2,3)
      source_r = rmax_structure_zone(n_structure_zone)
     &                  - source_depth
      source_mt(1,1) = source_mt(1,1) * ( 10.d0**(nexp-25) )
      source_mt(2,2) = source_mt(2,2) * ( 10.d0**(nexp-25) )
      source_mt(3,3) = source_mt(3,3) * ( 10.d0**(nexp-25) )
      source_mt(1,2) = source_mt(1,2) * ( 10.d0**(nexp-25) )
      source_mt(1,3) = source_mt(1,3) * ( 10.d0**(nexp-25) )
      source_mt(2,3) = source_mt(2,3) * ( 10.d0**(nexp-25) )
      source_mt(2,1) = source_mt(1,2)
      source_mt(3,1) = source_mt(1,3)
      source_mt(3,2) = source_mt(2,3)
c --- parameters for stations ---
      read(11,*) n_station
      if ( n_station.gt.max_nstation ) call error_handling(5)
      do 140 i=1,n_station
        read(11,*) station_lat(i),station_lon(i)
        call distaz(real(source_lat),real(source_lon),
     &                    real(station_lat(i)),real(station_lon(i)),
     &                    1,dist,az,baz,xdeg,nerr)
        station_theta(i) = dble(xdeg) * pi / 180.d0
        station_phi(i) = ( 180.d0 - dble(az) ) * pi / 180.d0
  140 continue
      do 150 i=1,n_station
        read(11,'(a80)') sac_file(i)
  150 continue
c closing the temporary file
      close(11)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_generation
     &        ( maxngrid_r,ngrid_r0,
     &         n_structure_zone,rmin_structure_zone,rmax_structure_zone,
     &          vsv_structure_zone,vsh_structure_zone,
     &          grid_r,source_r,ngrid_r )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the number and the location of grid points.
c    required subroutines: error_handling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,ngrid_r0,n_structure_zone,ngrid_r
      real*8 rmin_structure_zone(*),rmax_structure_zone(*)
      real*8 vsv_structure_zone(4,*)
      real*8 vsh_structure_zone(4,*)
      real*8 grid_r(*),source_r
c other variables
      integer ngrid,m_structure_zone
      integer i_zone,itmp,i
      real*8 rmin,rmax,rh
c
c computing the layer number of the top solid layer
      m_structure_zone = 0
      do 100 i_zone=1,n_structure_zone
        if ( ( ( vsv_structure_zone(1,i_zone).eq.0.d0 ).and.
     &         ( vsv_structure_zone(2,i_zone).eq.0.d0 ).and.
     &         ( vsv_structure_zone(3,i_zone).eq.0.d0 ).and.
     &         ( vsv_structure_zone(4,i_zone).eq.0.d0 )      ).or.
     &       ( ( vsh_structure_zone(1,i_zone).eq.0.d0 ).and.
     &         ( vsh_structure_zone(2,i_zone).eq.0.d0 ).and.
     &         ( vsh_structure_zone(3,i_zone).eq.0.d0 ).and.
     &         ( vsh_structure_zone(4,i_zone).eq.0.d0 )      ) ) then
          continue
        else
          m_structure_zone = i_zone
        endif
  100 continue
c
c computing the number and the location of the grid points
      rmin = rmin_structure_zone(1)
      rmax = rmax_structure_zone(n_structure_zone)
      grid_r(1) = rmin
      itmp = 1
      do 210 i_zone=1,n_structure_zone
        rh = rmax_structure_zone(i_zone) - rmin_structure_zone(i_zone)
        ngrid = int( dble(ngrid_r0) * rh / ( rmax - rmin ) ) + 1
        do 200 i=1,ngrid
          itmp = itmp + 1
          if ( itmp.gt.maxngrid_r ) call error_handling(11)
          grid_r(itmp) = rmin_structure_zone(i_zone)
     &                                  + rh * dble(i) / dble( ngrid )
c ------ putting an additional node on the source location
          if ( itmp.gt.2 ) then
            if ( ( grid_r(itmp-1).lt.source_r ).and.
     &                 ( source_r.lt.grid_r(itmp) ) ) then
              if ( itmp+1.gt.maxngrid_r ) call error_handling(11)
              grid_r(itmp+1) = grid_r(itmp)
              grid_r(itmp) = source_r
              itmp = itmp + 1
              if ( ( i_zone.eq.m_structure_zone ).and.
     &                   ( i.eq.ngrid )                      ) then
                if ( itmp+1.gt.maxngrid_r ) call error_handling(11)
                grid_r(itmp+1) = grid_r(itmp)
                grid_r(itmp) = ( grid_r(itmp) + source_r ) / 2.d0
                itmp = itmp + 1
              endif
            else
            if ( ( i_zone.eq.m_structure_zone ).and.
     &                   ( i.eq.ngrid ).and.
     &                   ( source_r.eq.grid_r(itmp) )        ) then
                if ( itmp+2.gt.maxngrid_r ) call error_handling(11)
                source_r = grid_r(itmp-1) * 0.001d0
     &                             + grid_r(itmp) * 0.999d0
                grid_r(itmp+2) = grid_r(itmp)
                grid_r(itmp+1) = ( grid_r(itmp) + source_r ) / 2.d0
                grid_r(itmp) = source_r
                itmp = itmp + 2
            endif
            endif
          endif
c ------
  200   continue
  210 continue
c recouting the total number of grid points
      ngrid_r = itmp
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine assign_structure
     &        ( ngrid_r,grid_r,
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c assigning the input model to the numerical grids
c    required subroutines: error_handling
c    required function: cal_PREM_structure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc DK DK added these flags
      include "flags.h"

c variables for input/output
      integer ngrid_r,n_structure_zone
      integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
      real*8 grid_r(*)
      real*8 rmin_structure_zone(*),rmax_structure_zone(*)
      real*8 rho_structure_zone(4,*)
      real*8 vpv_structure_zone(4,*),vph_structure_zone(4,*)
      real*8 vsv_structure_zone(4,*),vsh_structure_zone(4,*)
      real*8 eta_structure_zone(4,*)
      real*8 qkappa_structure_zone(*),qmu_structure_zone(*)
      real*8 grid_rho(2,*)
      real*8 grid_Ak(2,*),grid_Am(2,*),grid_Ck(2,*),grid_Cm(2,*)
      real*8 grid_L(2,*),grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
      real*8 grid_kappa(2,*),grid_mu(2,*)
      real*8 grid_qkappa(*),grid_qmu(*)
c other variables
      integer ir,izone,izone1,izone2
      real*8 rho_r(2),vpv_r(2),vph_r(2),vsv_r(2),vsh_r(2),eta_r(2)
      real*8 A_r(2),C_r(2),L_r(2),N_r(2),F_r(2)
      real*8 kappa_r(2),mu_r(2)
      real*8 cal_PREM_structure
c
      do 500 ir=1,ngrid_r-1
        izone1 = 0
        izone2 = 0
        do 150 izone=1,n_structure_zone
          if ( rmin_structure_zone(izone).le.grid_r(ir) )
     &            izone1 = izone
  150   continue
        do 160 izone=n_structure_zone,1,-1
          if ( rmax_structure_zone(izone).ge.grid_r(ir+1) )
     &            izone2 = izone
  160   continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        rho_r(1)
     &          = cal_PREM_structure( grid_r(ir),
     &                     rho_structure_zone(1,izone1),FLAG_RHO,izone1)
        rho_r(2)
     &          = cal_PREM_structure( grid_r(ir+1),
     &                     rho_structure_zone(1,izone2),FLAG_RHO,izone2)

cc DK DK added this to output and check the smoothed model
c     print *,sngl(grid_r(ir)),sngl(rho_r(1))
cc DK DK added this to output and check the smoothed model

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        vpv_r(1)
     &          = cal_PREM_structure( grid_r(ir),
     &                      vpv_structure_zone(1,izone1),FLAG_VP,izone1)
        vpv_r(2)
     &          = cal_PREM_structure( grid_r(ir+1),
     &                      vpv_structure_zone(1,izone2),FLAG_VP,izone2)
        vph_r(1)
     &          = cal_PREM_structure( grid_r(ir),
     &                      vph_structure_zone(1,izone1),FLAG_VP,izone1)
        vph_r(2)
     &          = cal_PREM_structure( grid_r(ir+1),
     &                      vph_structure_zone(1,izone2),FLAG_VP,izone2)

cc DK DK added this to output and check the smoothed model
c     print *,sngl(grid_r(ir)),sngl(vpv_r(1))
cc DK DK added this to output and check the smoothed model

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        vsv_r(1)
     &          = cal_PREM_structure( grid_r(ir),
     &                      vsv_structure_zone(1,izone1),FLAG_VS,izone1)
        vsv_r(2)
     &          = cal_PREM_structure( grid_r(ir+1),
     &                      vsv_structure_zone(1,izone2),FLAG_VS,izone2)
        vsh_r(1)
     &          = cal_PREM_structure( grid_r(ir),
     &                      vsh_structure_zone(1,izone1),FLAG_VS,izone1)
        vsh_r(2)
     &          = cal_PREM_structure( grid_r(ir+1),
     &                      vsh_structure_zone(1,izone2),FLAG_VS,izone2)

cc DK DK added this to output and check the smoothed model
c     print *,sngl(grid_r(ir)),sngl(vsv_r(1))
cc DK DK added this to output and check the smoothed model

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        eta_r(1)
     &          = cal_PREM_structure( grid_r(ir),
     &                     eta_structure_zone(1,izone1),FLAG_ETA,izone1)
        eta_r(2)
     &          = cal_PREM_structure( grid_r(ir+1),
     &                     eta_structure_zone(1,izone2),FLAG_ETA,izone2)

cc DK DK added this to output and check the smoothed model
c     print *,sngl(grid_r(ir)),sngl(eta_r(1))
cc DK DK added this to output and check the smoothed model

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        A_r(1) = rho_r(1) * vph_r(1) * vph_r(1)
        A_r(2) = rho_r(2) * vph_r(2) * vph_r(2)
        C_r(1) = rho_r(1) * vpv_r(1) * vpv_r(1)
        C_r(2) = rho_r(2) * vpv_r(2) * vpv_r(2)
        L_r(1) = rho_r(1) * vsv_r(1) * vsv_r(1)
        L_r(2) = rho_r(2) * vsv_r(2) * vsv_r(2)
        N_r(1) = rho_r(1) * vsh_r(1) * vsh_r(1)
        N_r(2) = rho_r(2) * vsh_r(2) * vsh_r(2)
        F_r(1) = eta_r(1) * ( A_r(1) - 2.d0 * L_r(1) )
        F_r(2) = eta_r(2) * ( A_r(2) - 2.d0 * L_r(2) )

        kappa_r(1) = ( 4.d0 * A_r(1) + C_r(1)
     &                         + 4.d0 * F_r(1) - 4.d0 * N_r(1) ) / 9.d0
        kappa_r(2) = ( 4.d0 * A_r(2) + C_r(2)
     &                         + 4.d0 * F_r(2) - 4.d0 * N_r(2) ) / 9.d0
        mu_r(1) = ( A_r(1) + C_r(1) - 2.d0 * F_r(1)
     &                      + 5.d0 * N_r(1) + 6.d0 * L_r(1) ) / 15.d0
        mu_r(2) = ( A_r(2) + C_r(2) - 2.d0 * F_r(2)
     &                      + 5.d0 * N_r(2) + 6.d0 * L_r(2) ) / 15.d0
        grid_rho(1,ir) = rho_r(1)
        grid_rho(2,ir) = rho_r(2)
        grid_Ak(1,ir) = A_r(1)
     &                          * kappa_r(1)
     &                            / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
        grid_Ak(2,ir) = A_r(2)
     &                          * kappa_r(2)
     &                            / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
        grid_Am(1,ir) = A_r(1)
     &                          * 4.d0/3.d0*mu_r(1)
     &                            / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
        grid_Am(2,ir) = A_r(2)
     &                          * 4.d0/3.d0*mu_r(2)
     &                            / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
        grid_Ck(1,ir) = C_r(1)
     &                          * kappa_r(1)
     &                            / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
        grid_Ck(2,ir) = C_r(2)
     &                          * kappa_r(2)
     &                            / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
        grid_Cm(1,ir) = C_r(1)
     &                          * 4.d0/3.d0*mu_r(1)
     &                            / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
        grid_Cm(2,ir) = C_r(2)
     &                          * 4.d0/3.d0*mu_r(2)
     &                            / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
        grid_L(1,ir) = L_r(1)
        grid_L(2,ir) = L_r(2)
        grid_N(1,ir) = N_r(1)
        grid_N(2,ir) = N_r(2)
        grid_Fk(1,ir) = F_r(1)
     &                          * kappa_r(1)
     &                            / ( kappa_r(1)-2.d0/3.d0*mu_r(1) )
        grid_Fk(2,ir) = F_r(2)
     &                          * kappa_r(2)
     &                            / ( kappa_r(2)-2.d0/3.d0*mu_r(2) )
        grid_Fm(1,ir) = F_r(1)
     &                          * (-2.d0/3.d0)*mu_r(1)
     &                            / ( kappa_r(1)-2.d0/3.d0*mu_r(1) )
        grid_Fm(2,ir) = F_r(2)
     &                          * (-2.d0/3.d0)*mu_r(2)
     &                            / ( kappa_r(2)-2.d0/3.d0*mu_r(2) )
        grid_kappa(1,ir) = kappa_r(1)
        grid_kappa(2,ir) = kappa_r(2)
        grid_mu(1,ir) = mu_r(1)
        grid_mu(2,ir) = mu_r(2)
        grid_qkappa(ir)
     &          = 2.d0/(1.d0/qkappa_structure_zone(izone1)
     &                  +1.d0/qkappa_structure_zone(izone2))
        grid_qmu(ir)
     &          = 2.d0/(1.d0/qmu_structure_zone(izone1)
     &                  +1.d0/qmu_structure_zone(izone2))
  500 continue
c --- computing positions of in solution vectors
      idim1_sph = 1
      idim2_sph = ngrid_r
      idim1_tor = 0
      idim2_tor = 0
      do 600 ir=1,ngrid_r-1
        if ( grid_mu(1,ir)*grid_mu(2,ir).ne.0.d0 ) idim2_tor = ir+1
  600 continue
      if ( idim2_tor.eq.0 ) call error_handling(16)
      do 700 ir=1,idim2_tor-1
        if ( grid_mu(1,ir)*grid_mu(2,ir).eq.0.d0 ) then
          idim1_tor = ir
        endif
  700 continue
      idim1_tor = idim1_tor + 1
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine assign_source
     &          ( ngrid_r,grid_r,source_r,
     &            grid_L,grid_Ck,grid_Cm,grid_Fk,grid_Fm,
     &            grid_qkappa,grid_qmu,grid_mu,
     &            igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms,
     &            grid_qkappas,grid_qmus,
     &            idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing intermediate parameters used in assing source parameters
c to the numerical grids
c    required subroutines: error_handling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer ngrid_r,igrid_rs
      integer idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0
      real*8 grid_r(*),source_r
      real*8 grid_L(2,*),grid_Ck(2,*),grid_Cm(2,*)
      real*8 grid_Fk(2,*),grid_Fm(2,*)
      real*8 grid_qkappa(*),grid_qmu(*),grid_mu(2,*)
      real*8 grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
      real*8 grid_qkappas,grid_qmus
c other variables
      integer i,ipos1,ipos2,ipos3
c
c ---- serching the grid on the source
      igrid_rs = 0
      do 100 i=1,ngrid_r
        if ( grid_r(i).eq.source_r ) igrid_rs = i
  100 continue
      if ( igrid_rs.eq.0 ) call error_handling(17)
c ---- computing the elastic constants at the source
      if ( igrid_rs.eq.ngrid_r ) then
        grid_Ls = grid_L(2,igrid_rs-1)
        grid_Cks = grid_Ck(2,igrid_rs-1)
        grid_Cms = grid_Cm(2,igrid_rs-1)
        grid_Fks = grid_Fk(2,igrid_rs-1)
        grid_Fms = grid_Fm(2,igrid_rs-1)
        grid_qkappas = grid_qkappa(igrid_rs-1)
        grid_qmus = grid_qmu(igrid_rs-1)
      else
        grid_Ls = grid_L(1,igrid_rs)
        grid_Cks = grid_Ck(1,igrid_rs)
        grid_Cms = grid_Cm(1,igrid_rs)
        grid_Fks = grid_Fk(1,igrid_rs)
        grid_Fms = grid_Fm(1,igrid_rs)
        grid_qkappas = grid_qkappa(igrid_rs)
        grid_qmus = grid_qmu(igrid_rs)
      endif
c --- computing positions of non-zero elements in excitation vectors
      ipos1 = 0
      ipos3 = 0
      ipos2 = 0
      if ( igrid_rs.gt.1 ) then
        do 200 i=1,igrid_rs-1
          if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
            ipos1 = ipos1 + 1
            ipos3 = ipos3 + 1
            ipos2 = 0
          else
            ipos1 = ipos1 + 2
            ipos3 = ipos3 + 1
            ipos2 = ipos2 + 1
          endif
          if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &                 ( grid_mu(1,i+1).ne.0.d0 )    ) then
            ipos1 = ipos1 + 1
            ipos3 = ipos3 + 1
            ipos2 = 0
          endif
          if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &                 ( grid_mu(1,i+1).eq.0.d0 )    ) then
            ipos1 = ipos1 + 2
            ipos3 = ipos3 + 1
            ipos2 = 0
          endif
  200   continue
      endif
      idim_rs_sph = ipos1 + 1
      idim_rs_tor = ipos2 + 1
      idim_rs_sph0 = ipos3 + 1
      idim_rs_tor0 = 0
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine assign_station
     &          ( ngrid_r,grid_r,grid_mu,
     &            idim_station_sph,idim_station_tor,
     &            idim_station_sph0,idim_station_tor0 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c computing intermediate parameters used in assing station parameters
c    required subroutines: none
c    required function: error_handling
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer ngrid_r
      integer idim_station_sph,idim_station_tor,
     &              idim_station_sph0,idim_station_tor0
      real*8 grid_r(*),grid_mu(2,*)
c other variables
      integer ir,idim,i,ipos1,ipos2,ipos3
c
      idim = 0
      do 600 ir=1,ngrid_r-1
        if ( grid_mu(1,ir)*grid_mu(2,ir).ne.0.d0 ) idim = ir
  600 continue
      if ( idim.eq.0 ) call error_handling(18)
c
c --- computing positions of non-zero elements in excitation vectors
      ipos1 = 0
      ipos3 = 0
      ipos2 = 0
      if ( idim.gt.1 ) then
        do 200 i=1,idim-1
          if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
            ipos1 = ipos1 + 1
            ipos3 = ipos3 + 1
            ipos2 = 0
          else
            ipos1 = ipos1 + 2
            ipos3 = ipos3 + 1
            ipos2 = ipos2 + 1
          endif
          if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &                 ( grid_mu(1,i+1).ne.0.d0 )    ) then
            ipos1 = ipos1 + 1
            ipos3 = ipos3 + 1
            ipos2 = 0
          endif
          if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &                 ( grid_mu(1,i+1).eq.0.d0 )    ) then
            ipos1 = ipos1 + 2
            ipos3 = ipos3 + 1
            ipos2 = 0
          endif
  200   continue
      endif
      idim_station_sph = ipos1 + 1 + 2
      idim_station_tor = ipos2 + 1 + 1
      idim_station_sph0 = ipos3 + 1 + 1
      idim_station_tor0 = 0
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function cal_PREM_structure( r,param,iflag ,izone)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 r,param(4)
      real*8 rmax,a

cc DK DK added these flags
      include "flags.h"
      integer NPOIN,i,iflag,izone
      parameter (NPOIN = 1193)
      double precision x1
      double precision rho(NPOIN)
      double precision vp(NPOIN)
      double precision vs(NPOIN)
      double precision deltax
      parameter (deltax = 1.d0 / 10000.d0)

      rmax = 6371.d0
      a = r / rmax
      cal_PREM_structure
     &        = param(1)
     &          + param(2) * a
     &          + param(3) * a * a
     &          + param(4) * a * a * a

cccccccccc DK DK better/safer test      if(a .ge. x1) then
cccccccccc DK DK replace only in the last layer (layer number 5)
      if(izone .eq. 5) then

cc DK DK added this to smooth the upper mantle
      include "values_smoothed_model.dat"

        i = nint((a - x1) / deltax) + 1

cc DK DK avoid edge effects
        if(i .lt. 1) i = 1
        if(i .gt. NPOIN) i = NPOIN

        if(iflag .eq. FLAG_RHO) then
          cal_PREM_structure = rho(i)

        else if(iflag .eq. FLAG_VP) then
          cal_PREM_structure = vp(i)

        else if(iflag .eq. FLAG_VS) then
          cal_PREM_structure = vs(i)

        else if(iflag .eq. FLAG_ETA) then
          cal_PREM_structure = 1.d0

        else
          stop 'wrong flag in smoothed model by DK'
        endif

      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_spectrum_file
     &          ( i_frequency,n_station,spectrum_file,
     &            station_displacement )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c writing the displacement at each station to the spectrum files.
c   required subroutine: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer i_frequency,n_station
      character*80 spectrum_file(*)
      complex*16 station_displacement(3,*)
c other variables
      integer i_station
c
      if ( i_frequency.eq.0 ) then
        do 100 i_station=1,n_station
          open( unit=11,file=spectrum_file(i_station),
     &                status='unknown' )
          write(11,*) i_frequency,
     &                      station_displacement(1,i_station),
     &                      station_displacement(2,i_station),
     &                      station_displacement(3,i_station)
          close(11)
  100   continue
      else
        do 110 i_station=1,n_station
          open( unit=11,file=spectrum_file(i_station),
     &                access='append',status='old' )
          write(11,*) i_frequency,
     &                      station_displacement(1,i_station),
     &                      station_displacement(2,i_station),
     &                      station_displacement(3,i_station)
          close(11)
  110   continue
      endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine distaz(the,phe,ths,phs,ns,dist,az,baz,xdeg,nerr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*=====================================================================
* PURPOSE:  To compute the distance and azimuth between locations.
*=====================================================================
* INPUT ARGUMENTS:
*    THE:     Event latitude in decimal degrees, North positive. [r]
*    PHE:     Event longitude, East positive. [r]
*    THS:     Array of station latitudes. [r]
*    PHS:     Array of station longitudes. [r]
*    NS:      Length of THS and PHS. [i]
*=====================================================================
* OUTPUT ARGUMENTS:
*    DIST:    Array of epicentral distances in km. [r]
*    AZ:      Array of azimuths in degrees. [r]
*    BAZ:     Array of back azimuths. [r]
*    XDEG:    Array of great circle arc lengths. [r]
*    NERR:    Error flag:
*             =    0   No error.
*             = 0904   Calculation failed internal consistency checks.
*=====================================================================
* MODULE/LEVEL:  DFM/4
*=====================================================================
* GLOBAL INPUT:
*    MACH:
*=====================================================================
* SUBROUTINES CALLED:
*    SACLIB:  SETMSG, APCMSG
*=====================================================================
* LOCAL VARIABLES:
*=====================================================================
* KNOWN ERRORS:
* - Problem with equation for distance. See discussion below.
*=====================================================================

      real pi,todeg,torad
      parameter (PI=3.141592654)
      parameter (TODEG=57.29577950)
      parameter (TORAD=1./TODEG)
c
      real ths(ns), phs(ns)
      real dist(ns), az(ns), baz(ns), xdeg(ns)
      integer ns,nerr,i
      logical laz,lbaz,ldist,lxdeg
      real the,phe,ec2,onemec2,eps,temp,therad,pherad,thg
      real d,e,f,c,a,b,g,h,thsrad,phsrad
      real d1,e1,f1,c1,a1,b1,g1,h1,sc,sd,ss,t1,p1,t2,p2,el
      real costhi,costhk,sinthi,sinthk,tanthi,tanthk,al,dl
      real a12top,a12bot,a12,cosa12,sina12,e2,e3,c0,c2,c4
      real v1,v2,z1,z2,x2,y2
      real e1p1,sqrte1p1,u1bot,u1,u2top,u2bot,u2,b0,du,pdist
      real rad,fl,twopideg,degtokm
      real c00,c01,c02,c03,c21,c22,c23,c42,c43
* PROCEDURE:

* - Calculations are based upon the reference spheroid of 1968 and
*   are defined by the major radius (RAD) and the flattening (FL).

      data rad/6378.160/,fl/0.00335293/
      data twopideg/360./
      data c00,c01,c02,c03/1.,0.25,-4.6875e-02,1.953125e-02/
      data c21,c22,c23/-0.125,3.125e-02,-1.46484375e-02/
      data c42,c43/-3.90625e-03,2.9296875e-03/
      data degtokm/111.3199/

* - Initialize.

      nerr=0
      ec2=2.*fl-fl*fl
      onemec2=1.-ec2
      eps=1.+ec2/onemec2

* - Check which output items are required.

      laz=.true.
c      if(az(1).lt.0.)laz=.false.
      lbaz=.true.
c      if(baz(1).lt.0.)lbaz=.false.
      ldist=.true.
c      if(dist(1).lt.0.)ldist=.false.
      lxdeg=.true.
c      if(xdeg(1).lt.0.)lxdeg=.false.

* - Convert event location to radians.
*   (Equations are unstable for latidudes of exactly 0 degrees.)

      temp=the
      if(temp.eq.0.)temp=1.0e-08
      therad=torad*temp
      pherad=torad*phe

* - Must convert from geographic to geocentric coordinates in order
*   to use the spherical trig equations.  This requires a latitude
*   correction given by: 1-EC2=1-2*FL+FL*FL

      thg=atan(onemec2*tan(therad))
      d=sin(pherad)
      e=-cos(pherad)
      f=-cos(thg)
      c=sin(thg)
      a= f*e
      b=-f*d
      g=-c*e
      h=c*d

* - Loop on stations:

      do 5000 i=1,ns

* -- Convert to radians.
        temp=ths(i)
        if(temp.eq.0.)temp=1.0e-08
        thsrad=torad*temp
        phsrad=torad*phs(i)

* -- Calculate some trig constants.
        thg=atan(onemec2*tan(thsrad))
        d1=sin(phsrad)
        e1=-cos(phsrad)
        f1=-cos(thg)
        c1=sin(thg)
        a1=f1*e1
        b1=-f1*d1
        g1=-c1*e1
        h1=c1*d1
        sc=a*a1+b*b1+c*c1

* - Spherical trig relationships used to compute angles.

        if(lxdeg)then
          sd=0.5*sqrt(((a-a1)**2+(b-b1)**2+(c-c1)**2)*((a+a1)**2
     #       +(b+b1)**2+(c+c1)**2))
          xdeg(i)=atan2(sd,sc)*todeg
          if(xdeg(i).lt.0.)xdeg(i)=xdeg(i)+twopideg
        endif
        if(laz)then
          ss = ((a1-d)**2+(b1-e)**2+c1**2-2.)
          sc = ((a1-g)**2+(b1-h)**2+(c1-f)**2-2.)
          az(i)=atan2(ss,sc)*todeg
          if(az(i).lt.0.)az(i)=az(i)+twopideg
        endif
        if(lbaz)then
          ss=((a-d1)**2+(b-e1)**2+c**2-2.)
          sc=((a-g1)**2+(b-h1)**2+(c-f1)**2-2.)
          baz(i)=atan2(ss,sc)*todeg
          if(baz(i).lt.0.)baz(i)=baz(i)+twopideg
        endif

* - Now compute the distance between the two points using Rudoe's
*   formula given in GEODESY, section 2.15(b).
*   (There is some numerical problem with the following formulae.
*   If the station is in the southern hemisphere and the event in
*   in the northern, these equations give the longer, not the
*   shorter distance between the two locations.  Since the equations
*   are fairly messy, the simplist solution is to reverse the
*   meanings of the two locations for this case.)
        if(ldist)then
          if(thsrad.gt.0.)then
            t1=thsrad
            p1=phsrad
            t2=therad
            p2=pherad
          else
            t1=therad
            p1=pherad
            t2=thsrad
            p2=phsrad
          endif
          el=ec2/onemec2
          e1=1.+el
          costhi=cos(t1)
          costhk=cos(t2)
          sinthi=sin(t1)
          sinthk=sin(t2)
          tanthi=sinthi/costhi
          tanthk=sinthk/costhk
          al=tanthi/(e1*tanthk)+
     #       ec2*sqrt((e1+tanthi**2)/(e1+tanthk**2))
          dl=p1-p2
          a12top=sin(dl)
          a12bot=(al-cos(dl))*sinthk
          a12=atan2(a12top,a12bot)
          cosa12=cos(a12)
          sina12=sin(a12)
          e1=el*((costhk*cosa12)**2+sinthk**2)
          e2=e1*e1
          e3=e1*e2
          c0=c00+c01*e1+c02*e2+c03*e3
          c2=c21*e1+c22*e2+c23*e3
          c4=c42*e2+c43*e3
          v1=rad/sqrt(1.-ec2*sinthk**2)
          v2=rad/sqrt(1.-ec2*sinthi**2)
          z1=v1*(1.-ec2)*sinthk
          z2=v2*(1.-ec2)*sinthi
          x2=v2*costhi*cos(dl)
          y2=v2*costhi*sin(dl)
          e1p1=e1+1.
          sqrte1p1=sqrt(e1p1)
          u1bot=sqrte1p1*cosa12
          u1=atan2(tanthk,u1bot)
          u2top=v1*sinthk+e1p1*(z2-z1)
          u2bot=sqrte1p1*(x2*cosa12-y2*sinthk*sina12)
          u2=atan2(u2top,u2bot)
          b0=v1*sqrt(1.+el*(costhk*cosa12)**2)/e1p1
          du=u2 -u1
          pdist=b0*(c2*(sin(2.*u2)-sin(2.*u1))+
     #       c4*(sin(4.*u2)-sin(4.*u1)))
          dist(i)=abs(b0*c0*du+pdist)
          if(lxdeg .and. (abs(dist(i)-degtokm*xdeg(i))).gt.100.)then
            nerr=0904
c            call setmsg('ERROR',nerr)
c            call apimsg(i)
          endif
        endif
 5000   continue

 8888 return

*=====================================================================
* MODIFICATION HISTORY:
*    830603:  Fixed bug with negative station latiudes.
*    810000:  Original version.
*=====================================================================

      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine error_handling(id)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c error handling.
c    required subroutines: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer id
      integer lnblnk
      character*72 message(99)
c
      message(1) = 'ngrid_r is too large (parameter_input).'
      message(2) = 'n_structure_zone is too large (parameter_input).'
      message(5) = 'n_station is too large (parameter_input).'
      message(11) = 'ngrid_r is too large
     & (grid_generation).'
      message(16) = 'Something is wrong (assign_structure).'
      message(17) = 'Something is wrong (assign_source).'
      message(18) = 'Something is wrong (assign_station).'
      message(41) = 'Bad arguments (comp_vecsph_tor)'
      message(42) = 'Bad arguments (plgndr)'
      message(51) = 'Bad arguments (comp_excitation)'
      message(51) = 'Bad arguments (comp_excitation0)'
      message(53) = 'Bad arguments (comp_wavefield)'
      message(54) = 'Bad arguments (comp_wavefield0)'
      message(55) = 'Bad arguments (comp_displacement_station)'
      message(56) = 'Bad arguments (comp_displacement_station0)'
c
      write(6,*) '******************** ERROR ********************'
      write(6,*) message(id)(1:lnblnk(message(id)))
      write(6,*) '***********************************************'
      stop
c
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_amp_significance
     &        ( maxngrid_r,l,amp_max,i_significance,idim0,
     &          init_npos_sph,init_npos_tor,
     &          idim1_sph,idim2_sph,idim1_tor,idim2_tor,
     &          grid_mu,idim_station_sph,idim_station_tor,
     &          whole_vector_sph,whole_vector_tor,work_vector )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c variables for input/output
      integer maxngrid_r,l,i_significance,idim0
      integer init_npos_sph,init_npos_tor
      integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
      integer idim_station_sph,idim_station_tor
      real*8 amp_max,grid_mu(2,*)
      complex*16 whole_vector_sph(2*maxngrid_r,-2:2)
      complex*16 whole_vector_tor(maxngrid_r,-2:2)
      real*8 work_vector(maxngrid_r,2)
c other variables
      integer i,ipos1,ipos2,itmp,idim0_sph,idim0_tor
      real*8 amp,eps,amaxr_sph,amaxr_tor,fac
c
      data eps / 1.d-12 /
c
      if ( mod(l,100).eq.0 ) then
        if ( l.eq.0 ) then
          amp_max = cdabs( whole_vector_sph(idim_station_sph,0) )**2
          i_significance = 1
c
          ipos1 = 0
          amaxr_sph = 0.d0
          amaxr_tor = 0.d0
          call init_complex_array( maxngrid_r,work_vector )
          if ( idim2_sph.gt.idim1_sph ) then
c introducing 'fac' to avoid underflow exceptions
            fac = 1.d100
            do 100 i=idim1_sph,idim2_sph-1
              if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
                ipos1 = ipos1 + 1
              else
                ipos1 = ipos1 + 1
                work_vector(i,1)
     &                  = cdabs( whole_vector_sph(ipos1,0)*fac )**2
                if ( work_vector(i,1).gt.amaxr_sph )
     &                  amaxr_sph = work_vector(i,1)
              endif
              if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &                     ( grid_mu(1,i+1).ne.0.d0 )    ) then
                ipos1 = ipos1 + 1
              endif
              if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &                     ( grid_mu(1,i+1).eq.0.d0 )    ) then
                ipos1 = ipos1 + 1
              endif
  100       continue
            itmp = idim2_sph
            do 110 i=idim2_sph-1,idim1_sph,-1
              if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
                if ( itmp.eq.i+1 ) itmp = i
              else
c               if ( work_vector(i,1).gt.eps*amaxr_sph ) itmp=i
                if ( work_vector(i,1).ge.eps*amaxr_sph ) itmp=i
              endif
  110       continue
            if ( itmp.ne.idim2_sph ) idim0_sph = itmp
            if ( idim0_sph.gt.2 ) idim0 = idim0_sph
          endif
c
        else
          amp =
     &              cdabs( whole_vector_sph(idim_station_sph,  -2) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph+1,-2) )**2
     &            + cdabs( whole_vector_tor(idim_station_tor,  -2) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph,  -1) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph+1,-1) )**2
     &            + cdabs( whole_vector_tor(idim_station_tor,  -1) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph,   0) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph+1, 0) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph,   1) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph+1, 1) )**2
     &            + cdabs( whole_vector_tor(idim_station_tor,   1) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph,   2) )**2
     &            + cdabs( whole_vector_sph(idim_station_sph+1, 2) )**2
     &            + cdabs( whole_vector_tor(idim_station_tor,   2) )**2
          if ( amp.gt.amp_max ) amp_max = amp
          if ( amp.lt.eps*amp_max ) i_significance = 0
c
          ipos1 = 0
          ipos2 = 0
          amaxr_sph = 0.d0
          amaxr_tor = 0.d0
          call init_complex_array( maxngrid_r,work_vector )
          ipos1 = 0
          ipos2 = 0
          if ( idim2_sph.gt.idim1_sph ) then
c introducing 'fac' to avoid underflow exceptions
            fac = 1.d100
            do 200 i=idim1_sph,idim2_sph-1
              if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
                ipos1 = ipos1 + 1
              else
                ipos1 = ipos1 + 2
                work_vector(i,1)
     &                  =   cdabs( whole_vector_sph(ipos1,  -2)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1+1,-2)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1,  -1)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1+1,-1)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1,   0)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1+1, 0)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1,   1)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1+1, 1)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1,   2)*fac )**2
     &                    + cdabs( whole_vector_sph(ipos1+1, 2)*fac )**2
                if ( work_vector(i,1).gt.amaxr_sph )
     &                  amaxr_sph = work_vector(i,1)
              endif
              if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &                     ( grid_mu(1,i+1).ne.0.d0 )    ) then
                ipos1 = ipos1 + 1
              endif
              if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &                     ( grid_mu(1,i+1).eq.0.d0 )    ) then
                ipos1 = ipos1 + 2
              endif
  200       continue
c
            do 210 i=idim1_tor,idim2_tor-1
              ipos2 = ipos2 + 1
              work_vector(i,2)
     &                =   cdabs( whole_vector_tor(ipos2,-2)*fac )**2
     &                  + cdabs( whole_vector_tor(ipos2,-1)*fac )**2
     &                  + cdabs( whole_vector_tor(ipos2, 1)*fac )**2
     &                  + cdabs( whole_vector_tor(ipos2, 2)*fac )**2
              if ( work_vector(i,2).gt.amaxr_tor )
     &                amaxr_tor = work_vector(i,2)
  210       continue
c
            itmp = idim2_sph
            do 220 i=idim2_sph-1,idim1_sph,-1
              if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
                if ( itmp.eq.i+1 ) itmp = i
              else
c               if ( work_vector(i,1).gt.eps*amaxr_sph ) itmp=i
                if ( work_vector(i,1).ge.eps*amaxr_sph ) itmp=i
              endif
  220       continue
            if ( itmp.ne.idim2_sph ) then
              idim0_sph = itmp
            else
              idim0_sph = idim1_sph
            endif
            itmp = idim2_tor
            do 230 i=idim2_tor-1,idim1_tor,-1
c             if ( work_vector(i,2).gt.eps*amaxr_tor ) itmp=i
              if ( work_vector(i,2).ge.eps*amaxr_tor ) itmp=i
  230       continue
            if ( itmp.ne.idim2_tor ) then
              idim0_tor = itmp
            else
              idim0_tor = idim1_tor
            endif
            if ( ( idim0.lt.min0(idim0_sph,idim0_tor) ).and.
     &                 ( min0(idim0_sph,idim0_tor).gt.2 ) )
     &              idim0 = min0(idim0_sph,idim0_tor)
          endif
        endif
c
        ipos1 = 0
        ipos2 = 0
        if ( idim0.gt.1 ) then
          do 300 i=1,idim0-1
            if ( grid_mu(1,i)*grid_mu(2,i).eq.0.d0 ) then
              ipos1 = ipos1 + 1
              ipos2 = 0
            else
              ipos1 = ipos1 + 2
              ipos2 = ipos2 + 1
            endif
            if ( ( grid_mu(2,i).eq.0.d0 ).and.
     &                   ( grid_mu(1,i+1).ne.0.d0 )    ) then
              ipos1 = ipos1 + 1
              ipos2 = 0
            endif
            if ( ( grid_mu(2,i).ne.0.d0 ).and.
     &                   ( grid_mu(1,i+1).eq.0.d0 )    ) then
              ipos1 = ipos1 + 2
              ipos2 = 0
            endif
  300     continue
        endif
        init_npos_sph = ipos1 + 1
        if ( idim0.ge.idim1_tor ) then
          init_npos_tor = ipos2 + 1
        else
          init_npos_tor = 1
        endif
c
      endif
c
      end
