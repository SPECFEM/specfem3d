
  subroutine param_input ( maxngrid_r,maxlmax,max_nstation,maxn_structure_zone, &
            time_series_length,n_frequency,omega_imag, &
            ngrid_r,lmin,lmax, &
            n_structure_zone,rmin_structure_zone,rmax_structure_zone, &
            rho_structure_zone,vpv_structure_zone,vph_structure_zone, &
            vsv_structure_zone,vsh_structure_zone,eta_structure_zone, &
            qkappa_structure_zone,qmu_structure_zone, &
            source_r,source_mt,source_depth,source_lat,source_lon, &
            n_station,station_lat,station_lon, &
            station_theta,station_phi,sac_file )

  implicit none

! input parameters

! variables for input/output
  integer maxngrid_r,maxlmax,max_nstation,maxn_structure_zone
  integer n_frequency,ngrid_r,lmin,lmax
  integer n_structure_zone,n_station
  real(kind=8) time_series_length,omega_imag
  real(kind=8) rmin_structure_zone(maxn_structure_zone)
  real(kind=8) rmax_structure_zone(maxn_structure_zone)
  real(kind=8) rho_structure_zone(4,maxn_structure_zone)
  real(kind=8) vpv_structure_zone(4,maxn_structure_zone)
  real(kind=8) vph_structure_zone(4,maxn_structure_zone)
  real(kind=8) vsv_structure_zone(4,maxn_structure_zone)
  real(kind=8) vsh_structure_zone(4,maxn_structure_zone)
  real(kind=8) eta_structure_zone(4,maxn_structure_zone)
  real(kind=8) qkappa_structure_zone(maxn_structure_zone)
  real(kind=8) qmu_structure_zone(maxn_structure_zone)
  real(kind=8) source_r,source_mt(3,3)
  real(kind=8) source_depth,source_lat,source_lon
  real(kind=8) station_lat(max_nstation),station_lon(max_nstation)
  real(kind=8) station_theta(max_nstation),station_phi(max_nstation)
  character(len=80) sac_file(max_nstation)

! other variables
  integer i,nexp,nerr
  real dist,az,baz,xdeg
  character(len=80) dummy,tmpfile

! constants
  real(kind=8), parameter :: pi=3.1415926535897932d0

  tmpfile = 'work'

! **********************************************************************
! reading input file from the standard input and writing out to the
! temporary file
! **********************************************************************

! opening a temporary file
  open( unit=11, file=tmpfile, status='unknown' )
! writing to the temporary file
  100 continue
  read(5,110) dummy
  110   format(a80)
  if ( dummy(1:1) == 'c' ) goto 100
  if ( dummy(1:3) == 'end' ) goto 120
  write(11,110) dummy
  goto 100
  120 continue
! closing the temporary file
  close(11)

! **********************************************************************
! reading the parameter from the temporary file
! **********************************************************************

! opening the temporary file
  open( unit=11, file=tmpfile, status='unknown' )

! reading the parameters

! ---- parameters for time series ---
  read(11,*) time_series_length,n_frequency
  read(11,*) omega_imag

! ---- parameters for numerical grids ----
  read(11,*) ngrid_r,lmin,lmax
  if ( ngrid_r > maxngrid_r ) call error_handling(1)
  if ( lmax > maxlmax ) call error_handling(3)

! ---- parameters for structure ---
  read(11,*) n_structure_zone
  if ( n_structure_zone > maxn_structure_zone ) call error_handling(2)

  do i=1,n_structure_zone
    read(11,*) rmin_structure_zone(i),rmax_structure_zone(i), &
                     rho_structure_zone(1,i),rho_structure_zone(2,i), &
                     rho_structure_zone(3,i),rho_structure_zone(4,i)
    read(11,*) vpv_structure_zone(1,i),vpv_structure_zone(2,i), &
                     vpv_structure_zone(3,i),vpv_structure_zone(4,i)
    read(11,*) vph_structure_zone(1,i),vph_structure_zone(2,i), &
                     vph_structure_zone(3,i),vph_structure_zone(4,i)
    read(11,*) vsv_structure_zone(1,i),vsv_structure_zone(2,i), &
                     vsv_structure_zone(3,i),vsv_structure_zone(4,i)
    read(11,*) vsh_structure_zone(1,i),vsh_structure_zone(2,i), &
                     vsh_structure_zone(3,i),vsh_structure_zone(4,i)
    read(11,*) eta_structure_zone(1,i),eta_structure_zone(2,i), &
                     eta_structure_zone(3,i),eta_structure_zone(4,i), &
                     qmu_structure_zone(i),qkappa_structure_zone(i)
  enddo

! ---- parameters for a source ---
  read(11,*) source_depth,source_lat,source_lon
  read(11,*) nexp, source_mt(1,1),source_mt(2,2),source_mt(3,3), &
                   source_mt(1,2),source_mt(1,3),source_mt(2,3)
  source_r = rmax_structure_zone(n_structure_zone) - source_depth
  source_mt(1,1) = source_mt(1,1) * ( 10.d0**(nexp-25) )
  source_mt(2,2) = source_mt(2,2) * ( 10.d0**(nexp-25) )
  source_mt(3,3) = source_mt(3,3) * ( 10.d0**(nexp-25) )
  source_mt(1,2) = source_mt(1,2) * ( 10.d0**(nexp-25) )
  source_mt(1,3) = source_mt(1,3) * ( 10.d0**(nexp-25) )
  source_mt(2,3) = source_mt(2,3) * ( 10.d0**(nexp-25) )
  source_mt(2,1) = source_mt(1,2)
  source_mt(3,1) = source_mt(1,3)
  source_mt(3,2) = source_mt(2,3)

! --- parameters for stations ---
  read(11,*) n_station
  if ( n_station > max_nstation ) call error_handling(5)

  do i=1,n_station
    read(11,*) station_lat(i),station_lon(i)
    call distaz(real(source_lat),real(source_lon),real(station_lat(i)),real(station_lon(i)),dist,az,baz,xdeg,nerr)
    station_theta(i) = dble(xdeg) * pi / 180.d0
    station_phi(i) = ( 180.d0 - dble(az) ) * pi / 180.d0
  enddo

  do i=1,n_station
    read(11,'(a80)') sac_file(i)
  enddo

! close the temporary file
  close(11)

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine grid_generation ( maxngrid_r,ngrid_r0, &
            n_structure_zone,rmin_structure_zone,rmax_structure_zone, &
            vsv_structure_zone,vsh_structure_zone, &
            grid_r,source_r,ngrid_r )

! compute the number and the location of grid points

  implicit none

! variables for input/output
  integer maxngrid_r,ngrid_r0,n_structure_zone,ngrid_r
  real(kind=8) rmin_structure_zone(*),rmax_structure_zone(*)
  real(kind=8) vsv_structure_zone(4,*)
  real(kind=8) vsh_structure_zone(4,*)
  real(kind=8) grid_r(*),source_r

! other variables
  integer ngrid,m_structure_zone
  integer i_zone,itmp,i
  real(kind=8) rmin,rmax,rh

! compute the layer number of the top solid layer
  m_structure_zone = 0

  do i_zone=1,n_structure_zone
    if ( ( ( vsv_structure_zone(1,i_zone)==0.d0 ).and. &
                 ( vsv_structure_zone(2,i_zone)==0.d0 ).and. &
                 ( vsv_structure_zone(3,i_zone)==0.d0 ).and. &
                 ( vsv_structure_zone(4,i_zone)==0.d0 )      ).or. &
               ( ( vsh_structure_zone(1,i_zone)==0.d0 ).and. &
                 ( vsh_structure_zone(2,i_zone)==0.d0 ).and. &
                 ( vsh_structure_zone(3,i_zone)==0.d0 ).and. &
                 ( vsh_structure_zone(4,i_zone)==0.d0 )      ) ) then
      continue !! DK DK i.e. do nothing
    else
      m_structure_zone = i_zone
    endif
  enddo

! compute the number and the location of the grid points
  rmin = rmin_structure_zone(1)
  rmax = rmax_structure_zone(n_structure_zone)
  grid_r(1) = rmin
  itmp = 1

  do i_zone=1,n_structure_zone

  rh = rmax_structure_zone(i_zone) - rmin_structure_zone(i_zone)
  ngrid = int( dble(ngrid_r0) * rh / ( rmax - rmin ) ) + 1

  do i=1,ngrid

    itmp = itmp + 1
    if ( itmp>maxngrid_r ) call error_handling(11)
    grid_r(itmp) = rmin_structure_zone(i_zone) + rh * dble(i) / dble( ngrid )

!! DK DK it could be that the source discontinuity is already handled as suggested by Vadim
!! DK DK (following eq 17 of Friederich_Dalkolmo_GEMINI_GJI_1995.pdf)
!! DK DK because the source code of Takeuchi contains this,
!! DK DK but I am not sure if this just adds a point at the source location (which is *not* sufficient)
!! DK DK or if there are more code lines in the rest of the code that handle the explicit treatment
!! DK DK of the source discontinuity. It will be important to check that.
!
! ------ putting an additional node on the source location
    if ( itmp>2 ) then
      if ( ( grid_r(itmp-1)<source_r ).and. ( source_r<grid_r(itmp) ) ) then
        if ( itmp+1>maxngrid_r ) call error_handling(11)
        grid_r(itmp+1) = grid_r(itmp)
        grid_r(itmp) = source_r
        itmp = itmp + 1
        if ( ( i_zone==m_structure_zone ).and. ( i==ngrid ) ) then
          if ( itmp+1>maxngrid_r ) call error_handling(11)
          grid_r(itmp+1) = grid_r(itmp)
          grid_r(itmp) = ( grid_r(itmp) + source_r ) / 2.d0
          itmp = itmp + 1
        endif
      else
      if ( ( i_zone==m_structure_zone ).and. ( i==ngrid ).and. ( source_r==grid_r(itmp) ) ) then
          if ( itmp+2>maxngrid_r ) call error_handling(11)
          source_r = grid_r(itmp-1) * 0.001d0 + grid_r(itmp) * 0.999d0
          grid_r(itmp+2) = grid_r(itmp)
          grid_r(itmp+1) = ( grid_r(itmp) + source_r ) / 2.d0
          grid_r(itmp) = source_r
          itmp = itmp + 2
      endif
      endif
    endif

    enddo ! of do i=1,ngrid
  enddo ! of do i_zone=1,n_structure_zone

! recouting the total number of grid points
  ngrid_r = itmp

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assign_structure ( ngrid_r,grid_r, &
            n_structure_zone,rmin_structure_zone,rmax_structure_zone, &
            rho_structure_zone, &
            vpv_structure_zone,vph_structure_zone, &
            vsv_structure_zone,vsh_structure_zone, &
            eta_structure_zone, &
            qkappa_structure_zone,qmu_structure_zone, &
            grid_rho,grid_Ak,grid_Am,grid_Ck,grid_Cm, &
            grid_L,grid_N,grid_Fk,grid_Fm, &
            grid_kappa,grid_mu,grid_qkappa,grid_qmu, &
            idim1_sph,idim2_sph,idim1_tor,idim2_tor )

! assigning the input model to the numerical grids

  implicit none

! variables for input/output
  integer ngrid_r,n_structure_zone
  integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
  real(kind=8) grid_r(*)
  real(kind=8) rmin_structure_zone(*),rmax_structure_zone(*)
  real(kind=8) rho_structure_zone(4,*)
  real(kind=8) vpv_structure_zone(4,*),vph_structure_zone(4,*)
  real(kind=8) vsv_structure_zone(4,*),vsh_structure_zone(4,*)
  real(kind=8) eta_structure_zone(4,*)
  real(kind=8) qkappa_structure_zone(*),qmu_structure_zone(*)
  real(kind=8) grid_rho(2,*)
  real(kind=8) grid_Ak(2,*),grid_Am(2,*),grid_Ck(2,*),grid_Cm(2,*)
  real(kind=8) grid_L(2,*),grid_N(2,*),grid_Fk(2,*),grid_Fm(2,*)
  real(kind=8) grid_kappa(2,*),grid_mu(2,*)
  real(kind=8) grid_qkappa(*),grid_qmu(*)

! other variables
  integer ir,izone,izone1,izone2
  real(kind=8) rho_r(2),vpv_r(2),vph_r(2),vsv_r(2),vsh_r(2),eta_r(2)
  real(kind=8) A_r(2),C_r(2),L_r(2),N_r(2),F_r(2)
  real(kind=8) kappa_r(2),mu_r(2)
  real(kind=8), external :: cal_PREM_structure

  do ir=1,ngrid_r-1

  izone1 = 0
  izone2 = 0

  do izone=1,n_structure_zone
    if ( rmin_structure_zone(izone)<=grid_r(ir) ) izone1 = izone
  enddo

  do izone=n_structure_zone,1,-1
    if ( rmax_structure_zone(izone)>=grid_r(ir+1) ) izone2 = izone
  enddo

  rho_r(1) = cal_PREM_structure( grid_r(ir), rho_structure_zone(1,izone1) )
  rho_r(2) = cal_PREM_structure( grid_r(ir+1), rho_structure_zone(1,izone2) )
  vpv_r(1) = cal_PREM_structure( grid_r(ir), vpv_structure_zone(1,izone1) )
  vpv_r(2) = cal_PREM_structure( grid_r(ir+1), vpv_structure_zone(1,izone2) )
  vph_r(1) = cal_PREM_structure( grid_r(ir), vph_structure_zone(1,izone1) )
  vph_r(2) = cal_PREM_structure( grid_r(ir+1), vph_structure_zone(1,izone2) )
  vsv_r(1) = cal_PREM_structure( grid_r(ir), vsv_structure_zone(1,izone1) )
  vsv_r(2) = cal_PREM_structure( grid_r(ir+1), vsv_structure_zone(1,izone2) )
  vsh_r(1) = cal_PREM_structure( grid_r(ir), vsh_structure_zone(1,izone1) )
  vsh_r(2) = cal_PREM_structure( grid_r(ir+1), vsh_structure_zone(1,izone2) )
  eta_r(1) = cal_PREM_structure( grid_r(ir), eta_structure_zone(1,izone1) )
  eta_r(2) = cal_PREM_structure( grid_r(ir+1), eta_structure_zone(1,izone2) )

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

  kappa_r(1) = ( 4.d0 * A_r(1) + C_r(1) + 4.d0 * F_r(1) - 4.d0 * N_r(1) ) / 9.d0
  kappa_r(2) = ( 4.d0 * A_r(2) + C_r(2) + 4.d0 * F_r(2) - 4.d0 * N_r(2) ) / 9.d0

  mu_r(1) = ( A_r(1) + C_r(1) - 2.d0 * F_r(1) + 5.d0 * N_r(1) + 6.d0 * L_r(1) ) / 15.d0
  mu_r(2) = ( A_r(2) + C_r(2) - 2.d0 * F_r(2) + 5.d0 * N_r(2) + 6.d0 * L_r(2) ) / 15.d0

  grid_rho(1,ir) = rho_r(1)
  grid_rho(2,ir) = rho_r(2)

  grid_Ak(1,ir) = A_r(1) * kappa_r(1) / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
  grid_Ak(2,ir) = A_r(2) * kappa_r(2) / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
  grid_Am(1,ir) = A_r(1) * 4.d0/3.d0*mu_r(1) / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
  grid_Am(2,ir) = A_r(2) * 4.d0/3.d0*mu_r(2) / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
  grid_Ck(1,ir) = C_r(1) * kappa_r(1) / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
  grid_Ck(2,ir) = C_r(2) * kappa_r(2) / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
  grid_Cm(1,ir) = C_r(1) * 4.d0/3.d0*mu_r(1) / ( kappa_r(1)+4.d0/3.d0*mu_r(1) )
  grid_Cm(2,ir) = C_r(2) * 4.d0/3.d0*mu_r(2) / ( kappa_r(2)+4.d0/3.d0*mu_r(2) )
  grid_L(1,ir) = L_r(1)
  grid_L(2,ir) = L_r(2)
  grid_N(1,ir) = N_r(1)
  grid_N(2,ir) = N_r(2)
  grid_Fk(1,ir) = F_r(1) * kappa_r(1) / ( kappa_r(1)-2.d0/3.d0*mu_r(1) )
  grid_Fk(2,ir) = F_r(2) * kappa_r(2) / ( kappa_r(2)-2.d0/3.d0*mu_r(2) )
  grid_Fm(1,ir) = F_r(1) * (-2.d0/3.d0)*mu_r(1) / ( kappa_r(1)-2.d0/3.d0*mu_r(1) )
  grid_Fm(2,ir) = F_r(2) * (-2.d0/3.d0)*mu_r(2) / ( kappa_r(2)-2.d0/3.d0*mu_r(2) )
  grid_kappa(1,ir) = kappa_r(1)
  grid_kappa(2,ir) = kappa_r(2)
  grid_mu(1,ir) = mu_r(1)
  grid_mu(2,ir) = mu_r(2)
  grid_qkappa(ir) = 2.d0/(1.d0/qkappa_structure_zone(izone1) +1.d0/qkappa_structure_zone(izone2))
  grid_qmu(ir) = 2.d0/(1.d0/qmu_structure_zone(izone1) +1.d0/qmu_structure_zone(izone2))

  enddo ! of do ir=1,ngrid_r-1

! --- compute positions of in solution vectors
  idim1_sph = 1
  idim2_sph = ngrid_r
  idim1_tor = 0
  idim2_tor = 0

  do ir=1,ngrid_r-1
    if ( grid_mu(1,ir)*grid_mu(2,ir) /= 0.d0 ) idim2_tor = ir+1
  enddo

  if ( idim2_tor == 0 ) call error_handling(16)

  do ir=1,idim2_tor-1
    if ( grid_mu(1,ir)*grid_mu(2,ir) == 0.d0 ) idim1_tor = ir
  enddo

  idim1_tor = idim1_tor + 1

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assign_source ( ngrid_r,grid_r,source_r, &
              grid_L,grid_Ck,grid_Cm,grid_Fk,grid_Fm, &
              grid_qkappa,grid_qmu,grid_mu, &
              igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms, &
              grid_qkappas,grid_qmus, &
              idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0 )

! compute intermediate parameters used in assing source parameters
! to the numerical grids

  implicit none

! variables for input/output
  integer ngrid_r,igrid_rs
  integer idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0
  real(kind=8) grid_r(*),source_r
  real(kind=8) grid_L(2,*),grid_Ck(2,*),grid_Cm(2,*)
  real(kind=8) grid_Fk(2,*),grid_Fm(2,*)
  real(kind=8) grid_qkappa(*),grid_qmu(*),grid_mu(2,*)
  real(kind=8) grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
  real(kind=8) grid_qkappas,grid_qmus

! other variables
  integer i,ipos1,ipos2,ipos3

! ---- serching the grid on the source
  igrid_rs = 0

  do i=1,ngrid_r
    if ( grid_r(i) == source_r ) igrid_rs = i
  enddo

  if ( igrid_rs == 0 ) call error_handling(17)

! ---- compute the elastic constants at the source
  if ( igrid_rs == ngrid_r ) then
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

! --- compute positions of non-zero elements in excitation vectors
  ipos1 = 0
  ipos3 = 0
  ipos2 = 0
  if ( igrid_rs > 1 ) then

  do i=1,igrid_rs-1
    if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
      ipos1 = ipos1 + 1
      ipos3 = ipos3 + 1
      ipos2 = 0
    else
      ipos1 = ipos1 + 2
      ipos3 = ipos3 + 1
      ipos2 = ipos2 + 1
    endif
    if ( ( grid_mu(2,i)==0.d0 ).and. ( grid_mu(1,i+1)/=0.d0 ) ) then
      ipos1 = ipos1 + 1
      ipos3 = ipos3 + 1
      ipos2 = 0
    endif
    if ( ( grid_mu(2,i)/=0.d0 ).and. ( grid_mu(1,i+1)==0.d0 ) ) then
      ipos1 = ipos1 + 2
      ipos3 = ipos3 + 1
      ipos2 = 0
    endif
  enddo

  endif

  idim_rs_sph = ipos1 + 1
  idim_rs_tor = ipos2 + 1
  idim_rs_sph0 = ipos3 + 1
  idim_rs_tor0 = 0

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assign_station ( ngrid_r,grid_mu, &
              idim_station_sph,idim_station_tor, &
              idim_station_sph0,idim_station_tor0 )

! compute intermediate parameters used in assing station parameters

  implicit none

! variables for input/output
  integer ngrid_r
  integer idim_station_sph,idim_station_tor, &
                idim_station_sph0,idim_station_tor0
  real(kind=8) grid_mu(2,*)

! other variables
  integer ir,idim,i,ipos1,ipos2,ipos3

  idim = 0

  do ir=1,ngrid_r-1
    if ( grid_mu(1,ir)*grid_mu(2,ir)/=0.d0 ) idim = ir
  enddo

  if ( idim==0 ) call error_handling(18)

! --- compute positions of non-zero elements in excitation vectors
  ipos1 = 0
  ipos3 = 0
  ipos2 = 0

  if ( idim>1 ) then

  do i=1,idim-1
    if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
      ipos1 = ipos1 + 1
      ipos3 = ipos3 + 1
      ipos2 = 0
    else
      ipos1 = ipos1 + 2
      ipos3 = ipos3 + 1
      ipos2 = ipos2 + 1
    endif
    if ( ( grid_mu(2,i)==0.d0 ).and. ( grid_mu(1,i+1)/=0.d0 ) ) then
      ipos1 = ipos1 + 1
      ipos3 = ipos3 + 1
      ipos2 = 0
    endif
    if ( ( grid_mu(2,i)/=0.d0 ).and. ( grid_mu(1,i+1)==0.d0 ) ) then
      ipos1 = ipos1 + 2
      ipos3 = ipos3 + 1
      ipos2 = 0
    endif
  enddo

  endif

  idim_station_sph = ipos1 + 1 + 2
  idim_station_tor = ipos2 + 1 + 1
  idim_station_sph0 = ipos3 + 1 + 1
  idim_station_tor0 = 0

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=8) function cal_PREM_structure( r,param )

  real(kind=8) r,param(4)
  real(kind=8) rmax,a

  rmax = 6371.d0
  a = r / rmax
  cal_PREM_structure = param(1) + param(2) * a + param(3) * a * a + param(4) * a * a * a

  end function cal_PREM_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_spectrum_file ( i_frequency,n_station,spectrum_file, station_displacement )

! writing the displacement at each station to the spectrum files

  implicit none

! variables for input/output
  integer i_frequency,n_station
  character(len=80) spectrum_file(*)
  complex(kind=8) station_displacement(3,*)

! other variables
  integer i_station

  if ( i_frequency==0 ) then

  do i_station=1,n_station
    open( unit=11,file=spectrum_file(i_station), status='unknown' )
    write(11,*) i_frequency, station_displacement(1,i_station), &
                             station_displacement(2,i_station), &
                             station_displacement(3,i_station)
    close(11)
  enddo

  else
  do i_station=1,n_station
    open( unit=11,file=spectrum_file(i_station), position='append',status='old' )
    write(11,*) i_frequency, station_displacement(1,i_station), &
                             station_displacement(2,i_station), &
                             station_displacement(3,i_station)
    close(11)
  enddo
  endif

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine distaz(the,phe,ths,phs,dist,az,baz,xdeg,nerr)

!=====================================================================
! PURPOSE:  compute the distance and azimuth between locations.
!=====================================================================
! INPUT ARGUMENTS:
!    THE:     Event latitude in decimal degrees, North positive. [r]
!    PHE:     Event longitude, East positive. [r]
!    THS:     Array of station latitudes. [r]
!    PHS:     Array of station longitudes. [r]
!    NS:      Length of THS and PHS. [i]
!=====================================================================
! OUTPUT ARGUMENTS:
!    DIST:    Array of epicentral distances in km. [r]
!    AZ:      Array of azimuths in degrees. [r]
!    BAZ:     Array of back azimuths. [r]
!    XDEG:    Array of great circle arc lengths. [r]
!    NERR:    Error flag:
!             =    0   No error.
!             = 0904   Calculation failed internal consistency checks.
!=====================================================================
! MODULE/LEVEL:  DFM/4
!=====================================================================
! GLOBAL INPUT:
!    MACH:
!=====================================================================
! SUBROUTINES CALLED:
!    SACLIB:  SETMSG, APCMSG
!=====================================================================
! LOCAL VARIABLES:
!=====================================================================
! KNOWN ERRORS:
! - Problem with equation for distance. See discussion below.
!=====================================================================

!=====================================================================
! MODIFICATION HISTORY:
!   1983/06/03:  Fixed bug with negative station latiudes.
!   1981/00/00:  Original version.
!=====================================================================

  implicit none

  real, parameter :: PI=3.141592654
  real, parameter :: TODEG=57.29577950
  real, parameter :: TORAD=1./TODEG

  integer nerr

  real :: ths, phs, dist, az, baz, xdeg

  logical laz,lbaz,ldist,lxdeg
  real the,phe,ec2,onemec2,eps,temp,therad,pherad,thg
  real d,e,f,c,a,b,g,h,thsrad,phsrad
  real d1,e1,f1,c1,a1,b1,g1,h1,sc,sd,ss,t1,p1,t2,p2,el
  real costhi,costhk,sinthi,sinthk,tanthi,tanthk,al,dl
  real a12top,a12bot,a12,cosa12,sina12,e2,e3,c0,c2,c4,v1,v2,z1,z2,x2,y2
  real e1p1,sqrte1p1,u1bot,u1,u2top,u2bot,u2,b0,du,pdist
  real rad,fl,twopideg,degtokm
  real c00,c01,c02,c03,c21,c22,c23,c42,c43

! - Calculations are based upon the reference spheroid of 1968 and
!   are defined by the major radius (RAD) and the flattening (FL).

  rad = 6378.160
  fl = 0.00335293
  twopideg = 360.

  c00 = 1.
  c01 = 0.25
  c02 = -4.6875e-02
  c03 = 1.953125e-02

  c21 = -0.125
  c22 = 3.125e-02
  c23 = -1.46484375e-02

  c42 = -3.90625e-03
  c43 = 2.9296875e-03

  degtokm = 111.3199

  nerr=0
  ec2=2.*fl-fl*fl
  onemec2=1.-ec2
  eps=1.+ec2/onemec2

! - Check which output items are required.

  laz=.true.
  lbaz=.true.
  ldist=.true.
  lxdeg=.true.

! - Convert event location to radians.
!   (Equations are unstable for latidudes of exactly 0 degrees.)

  temp=the
  if(temp==0.)temp=1.0e-08
  therad=torad*temp
  pherad=torad*phe

! - Must convert from geographic to geocentric coordinates in order
!   to use the spherical trig equations.  This requires a latitude
!   correction given by: 1-EC2=1-2*FL+FL*FL

  thg=atan(onemec2*tan(therad))
  d=sin(pherad)
  e=-cos(pherad)
  f=-cos(thg)
  c=sin(thg)
  a= f*e
  b=-f*d
  g=-c*e
  h=c*d

! - Loop on stations:

! -- Convert to radians.
  temp=ths
  if(temp==0.)temp=1.0e-08
  thsrad=torad*temp
  phsrad=torad*phs

! -- Calculate some trig constants.
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

! - Spherical trig relationships used to compute angles.

  if(lxdeg)then
    sd=0.5*sqrt(((a-a1)**2+(b-b1)**2+(c-c1)**2)*((a+a1)**2 +(b+b1)**2+(c+c1)**2))
    xdeg=atan2(sd,sc)*todeg
    if(xdeg<0.)xdeg=xdeg+twopideg
  endif
  if(laz)then
    ss = ((a1-d)**2+(b1-e)**2+c1**2-2.)
    sc = ((a1-g)**2+(b1-h)**2+(c1-f)**2-2.)
    az=atan2(ss,sc)*todeg
    if(az<0.)az=az+twopideg
  endif
  if(lbaz)then
    ss=((a-d1)**2+(b-e1)**2+c**2-2.)
    sc=((a-g1)**2+(b-h1)**2+(c-f1)**2-2.)
    baz=atan2(ss,sc)*todeg
    if(baz<0.)baz=baz+twopideg
  endif

! - Now compute the distance between the two points using Rudoe's
!   formula given in GEODESY, section 2.15(b).
!   (There is some numerical problem with the following formulae.
!   If the station is in the southern hemisphere and the event in
!   in the northern, these equations give the longer, not the
!   shorter distance between the two locations.  Since the equations
!   are fairly messy, the simplist solution is to reverse the
!   meanings of the two locations for this case.)
  if(ldist)then
    if(thsrad>0.)then
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
    al=tanthi/(e1*tanthk)+ ec2*sqrt((e1+tanthi**2)/(e1+tanthk**2))
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
    pdist=b0*(c2*(sin(2.*u2)-sin(2.*u1))+ c4*(sin(4.*u2)-sin(4.*u1)))
    dist=abs(b0*c0*du+pdist)
    if(lxdeg .and. (abs(dist-degtokm*xdeg))>100.)then
      nerr=0904
      stop 'error in distaz'
!     call setmsg('ERROR',nerr)
!     call apimsg
    endif
  endif

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine error_handling(id)

! error handling

  implicit none

  integer id
  character(len=72) message(99)

  message(1) = 'ngrid_r is too large (parameter_input).'
  message(2) = 'n_structure_zone is too large (parameter_input).'
  message(3) = 'lmax is too large (parameter_input).'
  message(5) = 'n_station is too large (parameter_input).'
  message(11) = 'ngrid_r is too large (grid_generation).'
  message(16) = 'Something is wrong (assign_structure).'
  message(17) = 'Something is wrong (assign_source).'
  message(18) = 'Something is wrong (assign_station).'
  message(40) = 'Bad arguments (comp_vecsph_station)'
  message(41) = 'Bad arguments (comp_vecsph_all)'
  message(42) = 'Bad arguments (compute_Legendre_polynomial)'
  message(51) = 'Bad arguments (comp_excitation)'
  message(51) = 'Bad arguments (comp_excitation0)'
  message(53) = 'Bad arguments (comp_wavefield)'
  message(54) = 'Bad arguments (comp_wavefield0)'
  message(55) = 'Bad arguments (comp_displacement_station)'
  message(56) = 'Bad arguments (comp_displacement_station0)'

  write(6,*) '******************** ERROR ********************'
  write(6,*) message(id)(1:len_trim(message(id)))
  write(6,*) '***********************************************'
  stop 'error detected'

  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_amp_significance ( maxngrid_r,l,amp_max,i_significance,idim0, &
            init_npos_sph,init_npos_tor, &
            idim1_sph,idim2_sph,idim1_tor,idim2_tor, &
            grid_mu,idim_station_sph,idim_station_tor, &
            whole_vector_sph,whole_vector_tor,work_vector )

  implicit none

! variables for input/output
  integer maxngrid_r,l,i_significance,idim0
  integer init_npos_sph,init_npos_tor
  integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
  integer idim_station_sph,idim_station_tor
  real(kind=8) amp_max,grid_mu(2,*)
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) whole_vector_tor(maxngrid_r,-2:2)
  real(kind=8) work_vector(maxngrid_r,2)

! other variables
  integer i,ipos1,ipos2,itmp,idim0_sph,idim0_tor
  real(kind=8) amp,eps,amaxr_sph,amaxr_tor,fac

  eps = 1.d-12

  if ( mod(l,100)==0 ) then
  if ( l==0 ) then
    amp_max = abs( whole_vector_sph(idim_station_sph,0) )**2
    i_significance = 1

    ipos1 = 0
    amaxr_sph = 0.d0
    amaxr_tor = 0.d0

!! DK DK to avoid a warning by the compiler because here the array type is real, not complex
!! DK DK    call init_complex_array( maxngrid_r,work_vector )
    work_vector(:,:) = 0.d0

    if ( idim2_sph>idim1_sph ) then
! introducing 'fac' to avoid underflow exceptions
      fac = 1.d100

      do i=idim1_sph,idim2_sph-1
        if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
          ipos1 = ipos1 + 1
        else
          ipos1 = ipos1 + 1
          work_vector(i,1) = abs( whole_vector_sph(ipos1,0)*fac )**2
          if ( work_vector(i,1)>amaxr_sph ) amaxr_sph = work_vector(i,1)
        endif
        if ( ( grid_mu(2,i)==0.d0 ).and. ( grid_mu(1,i+1)/=0.d0 ) ) then
          ipos1 = ipos1 + 1
        endif
        if ( ( grid_mu(2,i)/=0.d0 ).and. ( grid_mu(1,i+1)==0.d0 ) ) then
          ipos1 = ipos1 + 1
        endif
      enddo

      itmp = idim2_sph

      do i=idim2_sph-1,idim1_sph,-1
        if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
          if ( itmp==i+1 ) itmp = i
        else
!         if ( work_vector(i,1)>eps*amaxr_sph ) itmp=i
          if ( work_vector(i,1)>=eps*amaxr_sph ) itmp=i
        endif
      enddo

      if ( itmp/=idim2_sph ) idim0_sph = itmp
      if ( idim0_sph>2 ) idim0 = idim0_sph
    endif

  else
    amp =       abs( whole_vector_sph(idim_station_sph,  -2) )**2 &
              + abs( whole_vector_sph(idim_station_sph+1,-2) )**2 &
              + abs( whole_vector_tor(idim_station_tor,  -2) )**2 &
              + abs( whole_vector_sph(idim_station_sph,  -1) )**2 &
              + abs( whole_vector_sph(idim_station_sph+1,-1) )**2 &
              + abs( whole_vector_tor(idim_station_tor,  -1) )**2 &
              + abs( whole_vector_sph(idim_station_sph,   0) )**2 &
              + abs( whole_vector_sph(idim_station_sph+1, 0) )**2 &
              + abs( whole_vector_sph(idim_station_sph,   1) )**2 &
              + abs( whole_vector_sph(idim_station_sph+1, 1) )**2 &
              + abs( whole_vector_tor(idim_station_tor,   1) )**2 &
              + abs( whole_vector_sph(idim_station_sph,   2) )**2 &
              + abs( whole_vector_sph(idim_station_sph+1, 2) )**2 &
              + abs( whole_vector_tor(idim_station_tor,   2) )**2
    if ( amp>amp_max ) amp_max = amp
    if ( amp<eps*amp_max ) i_significance = 0

    ipos1 = 0
    ipos2 = 0
    amaxr_sph = 0.d0
    amaxr_tor = 0.d0

!! DK DK to avoid a warning by the compiler because here the array type is real, not complex
!! DK DK    call init_complex_array( maxngrid_r,work_vector )
    work_vector(:,:) = 0.d0

    ipos1 = 0
    ipos2 = 0
    if ( idim2_sph>idim1_sph ) then
! introducing 'fac' to avoid underflow exceptions
      fac = 1.d100

      do i=idim1_sph,idim2_sph-1
        if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
          ipos1 = ipos1 + 1
        else
          ipos1 = ipos1 + 2
          work_vector(i,1) =   abs( whole_vector_sph(ipos1,  -2)*fac )**2 &
                             + abs( whole_vector_sph(ipos1+1,-2)*fac )**2 &
                             + abs( whole_vector_sph(ipos1,  -1)*fac )**2 &
                             + abs( whole_vector_sph(ipos1+1,-1)*fac )**2 &
                             + abs( whole_vector_sph(ipos1,   0)*fac )**2 &
                             + abs( whole_vector_sph(ipos1+1, 0)*fac )**2 &
                             + abs( whole_vector_sph(ipos1,   1)*fac )**2 &
                             + abs( whole_vector_sph(ipos1+1, 1)*fac )**2 &
                             + abs( whole_vector_sph(ipos1,   2)*fac )**2 &
                             + abs( whole_vector_sph(ipos1+1, 2)*fac )**2
          if ( work_vector(i,1)>amaxr_sph ) amaxr_sph = work_vector(i,1)
        endif
        if ( ( grid_mu(2,i)==0.d0 ).and. ( grid_mu(1,i+1)/=0.d0 ) ) then
          ipos1 = ipos1 + 1
        endif
        if ( ( grid_mu(2,i)/=0.d0 ).and. ( grid_mu(1,i+1)==0.d0 ) ) then
          ipos1 = ipos1 + 2
        endif
      enddo

      do i=idim1_tor,idim2_tor-1
        ipos2 = ipos2 + 1
        work_vector(i,2) =   abs( whole_vector_tor(ipos2,-2)*fac )**2 &
                           + abs( whole_vector_tor(ipos2,-1)*fac )**2 &
                           + abs( whole_vector_tor(ipos2, 1)*fac )**2 &
                           + abs( whole_vector_tor(ipos2, 2)*fac )**2
        if ( work_vector(i,2)>amaxr_tor ) amaxr_tor = work_vector(i,2)
      enddo

      itmp = idim2_sph

      do i=idim2_sph-1,idim1_sph,-1
        if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
          if ( itmp==i+1 ) itmp = i
        else
!         if ( work_vector(i,1)>eps*amaxr_sph ) itmp=i
          if ( work_vector(i,1)>=eps*amaxr_sph ) itmp=i
        endif
      enddo

      if ( itmp/=idim2_sph ) then
        idim0_sph = itmp
      else
        idim0_sph = idim1_sph
      endif
      itmp = idim2_tor

      do i=idim2_tor-1,idim1_tor,-1
!       if ( work_vector(i,2)>eps*amaxr_tor ) itmp=i
        if ( work_vector(i,2)>=eps*amaxr_tor ) itmp=i
      enddo

      if ( itmp/=idim2_tor ) then
        idim0_tor = itmp
      else
        idim0_tor = idim1_tor
      endif
      if ( amaxr_tor/=0.d0 ) then
        if ( ( idim0<min0(idim0_sph,idim0_tor) ).and. ( min0(idim0_sph,idim0_tor)>2 ) ) idim0 = min0(idim0_sph,idim0_tor)
      else
        if ( idim0_sph>2 ) idim0 = idim0_sph
      endif
    endif
  endif

  ipos1 = 0
  ipos2 = 0

  if ( idim0>1 ) then

    do i=1,idim0-1
      if ( grid_mu(1,i)*grid_mu(2,i)==0.d0 ) then
        ipos1 = ipos1 + 1
        ipos2 = 0
      else
        ipos1 = ipos1 + 2
        ipos2 = ipos2 + 1
      endif
      if ( ( grid_mu(2,i)==0.d0 ).and. ( grid_mu(1,i+1)/=0.d0 ) ) then
        ipos1 = ipos1 + 1
        ipos2 = 0
      endif
      if ( ( grid_mu(2,i)/=0.d0 ).and. ( grid_mu(1,i+1)==0.d0 ) ) then
        ipos1 = ipos1 + 2
        ipos2 = 0
      endif
    enddo

  endif

  init_npos_sph = ipos1 + 1

  if ( idim0>=idim1_tor ) then
    init_npos_tor = ipos2 + 1
  else
    init_npos_tor = 1
  endif

  endif

  end

