
! DK DK determine if we compile in serial or in MPI mode
#ifndef USE_SERIAL
#define USE_MPI
#endif

  program dsmti

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computation of synthetic seismograms for spherically-symmetric
! transversely isotropic (TI) media by using the Direct Solution Method.

! Converted to Fortran90 and cleaned by Dimitri Komatitsch, CNRS, Marseille, France, March 2014

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

!
! Important remark from Nozomu Takeuchi:
!
! The code computes synthetics in the frequency domain and incorporates
! a small (constant) imaginary part in the angular frequency to reduce wraprounds.
! The obtained seismograms are therefore less accurate for lower frequencies.
! You might have to apply high-pass filters to remove low-frequency noise.
!

  use mpi

  implicit none

!
! Another important remark from Nozomu Takeuchi and Dimitri Komatitsch:
!
! The code computes velocity Green functions for Heaviside function source time history
! (or displacement Green functions for delta function source time history).
! To obtain displacement seismograms, we should integrate the obtained seismograms.
! The value below is used to select that:
!
  integer, parameter :: ITYPE_SEISMOGRAMS = 2   !   1 = displacement,  2 = velocity

! constants
  integer, parameter :: maxngrid_r = 100000
  integer, parameter :: maxlmax = 10000
  integer, parameter :: max_nstation = 70
  integer, parameter :: maxn_structure_zone = 15
  integer, parameter :: maxnfreq = 8192

  real(kind=8), parameter :: pi=3.1415926535897932d0

! variables for time series
  integer n_frequency,i_frequency
  real(kind=8) time_series_length,omega_imag

! variables for numerical grids
  integer ngrid_r0,lmin0,lmax0
  integer ngrid_r,l
  integer n_frequency_band,i_frequency_band
  integer boundary_r(maxn_structure_zone)
  real(kind=8) grid_r(maxngrid_r)
  real(kind=8) grid_rho(2,maxngrid_r)
  real(kind=8) grid_Ak(2,maxngrid_r),grid_Am(2,maxngrid_r)
  real(kind=8) grid_Ck(2,maxngrid_r),grid_Cm(2,maxngrid_r)
  real(kind=8) grid_L(2,maxngrid_r),grid_N(2,maxngrid_r)
  real(kind=8) grid_Fk(2,maxngrid_r),grid_Fm(2,maxngrid_r)
  real(kind=8) grid_kappa(2,maxngrid_r)
  real(kind=8) grid_mu(2,maxngrid_r)
  real(kind=8) grid_qkappa(maxngrid_r)
  real(kind=8) grid_qmu(maxngrid_r)

! variables for the structure
  integer n_structure_zone
  integer idim1_sph,idim2_sph,idim1_tor,idim2_tor
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

! variables for a source
  integer igrid_rs
  integer idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0
  real(kind=8) source_r,source_mt(3,3)
  real(kind=8) source_depth,source_lat,source_lon
  real(kind=8) grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms
  real(kind=8) grid_qkappas,grid_qmus

! variables for stations
  integer n_station
  integer idim_station_sph,idim_station_tor, &
                idim_station_sph0,idim_station_tor0
  real(kind=8) station_lat(max_nstation),station_lon(max_nstation)
  real(kind=8) station_theta(max_nstation),station_phi(max_nstation)
  complex(kind=8) vecsph_sph1(3,0:maxlmax,-2:2,max_nstation)
  complex(kind=8) vecsph_sph2(3,0:maxlmax,-2:2,max_nstation)
  complex(kind=8) vecsph_tor(3,0:maxlmax,-2:2,max_nstation)
  complex(kind=8) station_displacement(3,max_nstation,maxnfreq)
  character(len=80) sac_file(max_nstation),coutfile

! variables for matrix elements
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

! variables for wavefield
  integer i_significance,idim0,init_npos_sph,init_npos_tor
  real(kind=8) amp_max
  complex(kind=8) whole_matrix_sph(4,2*maxngrid_r)
  complex(kind=8) whole_matrix_tor(2,maxngrid_r)
  complex(kind=8) whole_matrix_dr_sph(2*maxngrid_r)
  complex(kind=8) whole_matrix_dr_tor(maxngrid_r)
  complex(kind=8) whole_vector_sph(2*maxngrid_r,-2:2)
  complex(kind=8) whole_vector_tor(maxngrid_r,-2:2)
  complex(kind=8) work_vector(2*maxngrid_r)
  complex(kind=8) work_spc(16*maxnfreq)
!! DK DK added this to avoid a incompatible type warning when calling routine four1() later
  real(kind=8) work_spc_real(2*16*maxnfreq)
  equivalence ( work_vector,work_spc,work_spc_real )
  real work_time(32*maxnfreq)

!! DK DK added this to define the Gauss points and weights
! variables for numerical integrations using Simpson's or Gauss' quadrature rule for mass lumping
  include "integration_points_for_submatrices.h"
  real(kind=8), dimension(:), allocatable :: xi_Gauss,weight_Gauss

  integer :: rank,size

!! DK DK for MPI
!! DK DK we could/should consider making these two constants a variable instead
!! DK DK and then use dynamic memory allocation instead of static allocation
!! DK DK for the size of the MPI buffers below
#ifdef USE_MPI
  integer, parameter :: mbfi=20000,mbfd=16000000
  integer, dimension(4) :: param
  integer, dimension(mbfi) :: ibuff
  real(kind=8), dimension(mbfd) :: dbuff
  integer :: ipos,dpos,irank,nsnd,ierr,info,jjj
  complex(kind(0e0)) :: displacementsngl(3,max_nstation)
  integer i,j
#endif

! **********************************************************************
! input parameters
! **********************************************************************

! initialize MPI
#ifdef USE_MPI
  call mpi_init( ierr )
  call mpi_comm_rank( mpi_comm_world,rank,ierr )
  call mpi_comm_size( mpi_comm_world,size,ierr )
#else
  rank = 0
  size = 1
#endif

  n_frequency_band = 4

#ifdef USE_MPI
  ipos = 0
  dpos = 0

  if ( rank==0 ) then
#endif

  call param_input( maxngrid_r,maxlmax,max_nstation,maxn_structure_zone, &
            time_series_length,n_frequency,omega_imag, &
            ngrid_r0,lmin0,lmax0, &
            n_structure_zone,rmin_structure_zone,rmax_structure_zone, &
            rho_structure_zone,vpv_structure_zone,vph_structure_zone, &
            vsv_structure_zone,vsh_structure_zone,eta_structure_zone, &
            qkappa_structure_zone,qmu_structure_zone, &
            source_r,source_mt,source_depth,source_lat,source_lon, &
            n_station,station_lat,station_lon, &
            station_theta,station_phi,sac_file )

#ifdef USE_MPI
  call mputi( n_frequency,1,info ,mbfi,ipos,ibuff)
  call mputi( ngrid_r0,1,info ,mbfi,ipos,ibuff)
  call mputi( lmin0,1,info ,mbfi,ipos,ibuff)
  call mputi( lmax0,1,info ,mbfi,ipos,ibuff)
  call mputi( n_structure_zone,1,info ,mbfi,ipos,ibuff)
  call mputi( n_station,1,info ,mbfi,ipos,ibuff)

  call mputd( time_series_length,1,info ,mbfd,dpos,dbuff)
  call mputd( omega_imag,1,info ,mbfd,dpos,dbuff)
  call mputd( rmin_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( rmax_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( rho_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( vpv_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( vph_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( vsv_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( vsh_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( eta_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( qkappa_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( qmu_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mputd( source_r,1,info ,mbfd,dpos,dbuff)
  call mputd( source_mt,9,info ,mbfd,dpos,dbuff)
  call mputd( station_theta,max_nstation,info ,mbfd,dpos,dbuff)
  call mputd( station_phi,max_nstation,info ,mbfd,dpos,dbuff)

  call msend( -1,10,info ,mbfi,mbfd,size,ipos,dpos,param,ibuff,dbuff)

  else

  call mrecv( 10,0,info ,mbfi,mbfd,ipos,dpos,param,ibuff,dbuff)

  call mgeti( n_frequency,1,info ,mbfi,ipos,ibuff)
  call mgeti( ngrid_r0,1,info ,mbfi,ipos,ibuff)
  call mgeti( lmin0,1,info ,mbfi,ipos,ibuff)
  call mgeti( lmax0,1,info ,mbfi,ipos,ibuff)
  call mgeti( n_structure_zone,1,info ,mbfi,ipos,ibuff)
  call mgeti( n_station,1,info ,mbfi,ipos,ibuff)

  call mgetd( time_series_length,1,info ,mbfd,dpos,dbuff)
  call mgetd( omega_imag,1,info ,mbfd,dpos,dbuff)
  call mgetd( rmin_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( rmax_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( rho_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( vpv_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( vph_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( vsv_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( vsh_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( eta_structure_zone,4*maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( qkappa_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( qmu_structure_zone,maxn_structure_zone,info ,mbfd,dpos,dbuff)
  call mgetd( source_r,1,info ,mbfd,dpos,dbuff)
  call mgetd( source_mt,9,info ,mbfd,dpos,dbuff)
  call mgetd( station_theta,max_nstation,info ,mbfd,dpos,dbuff)
  call mgetd( station_phi,max_nstation,info ,mbfd,dpos,dbuff)

  endif
#endif

  call comp_vecsph_station( maxlmax,lmax0,n_station,station_theta,station_phi, &
            vecsph_sph1,vecsph_sph2,vecsph_tor )

  if(USE_GAUSS_INTEGRATION) then
    allocate(xi_Gauss(ns))
    allocate(weight_Gauss(ns))
    call define_Gauss_points_and_weights(xi_Gauss,weight_Gauss,ns)
  else
!! DK DK allocate dummy values just to be able to use the arrays in subroutine calls
    allocate(xi_Gauss(1))
    allocate(weight_Gauss(1))
  endif

! **********************************************************************
! loop on frequency bands
! **********************************************************************
  do i_frequency_band=1,n_frequency_band

! **********************************************************************
! generating numerical grids
! **********************************************************************

  call grid_generation( maxngrid_r,ngrid_r_estimated &
!                (i_frequency_band,n_frequency_band,ngrid_r0), &
                (n_frequency_band,n_frequency_band,ngrid_r0), &
              n_structure_zone,rmin_structure_zone,rmax_structure_zone, &
              vsv_structure_zone,vsh_structure_zone, &
              grid_r,source_r,ngrid_r,boundary_r )
!write(*,*) 'rank,i_frequency_band,n_frequency_band,ngrid_r,ngrid_r0',rank,i_frequency_band,n_frequency_band,ngrid_r,ngrid_r0
!write(*,*) rank,'boundary of the every layer:',boundary_r
!if (rank == 0) then
!  do jjj =1,n_structure_zone + 1
!   write(*,*) 'boundary of every layer:',jjj,boundary_r(jjj),grid_r(boundary_r(jjj))
!  enddo
!endif

  call assign_structure( ngrid_r,grid_r, &
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
!write(*,*) rank,'grid_r(1),grid_r(ngrid_r:ngrid_r+1)',grid_r(1),grid_r(ngrid_r:ngrid_r+1),'grid_rho',grid_rho(1,1),grid_rho(1,ngrid_r-1:ngrid_r)

  call assign_source( ngrid_r,grid_r,source_r, &
              grid_L,grid_Ck,grid_Cm,grid_Fk,grid_Fm, &
              grid_qkappa,grid_qmu,grid_mu, &
              igrid_rs,grid_Ls,grid_Cks,grid_Cms,grid_Fks,grid_Fms, &
              grid_qkappas,grid_qmus, &
              idim_rs_sph,idim_rs_tor,idim_rs_sph0,idim_rs_tor0 )

  call assign_station( ngrid_r,grid_mu, &
              idim_station_sph,idim_station_tor, &
              idim_station_sph0,idim_station_tor0 )

! **********************************************************************
! compute submatrices for the elastic part of the medium
! **********************************************************************
  call comp_submatrix( ngrid_r,grid_r, &
              grid_rho,grid_kappa,grid_mu,grid_Ak,grid_Am, &
              grid_Ck,grid_Cm,grid_L,grid_N,grid_Fk,grid_Fm, &
              submatrix_I0,submatrix_I1k,submatrix_I1m,submatrix_I2, &
              submatrix_I3k,submatrix_I3m,submatrix_I4,submatrix_I5k, &
              submatrix_I5m,submatrix_I6,submatrix_I7,xi_Gauss,weight_Gauss,n_structure_zone,boundary_r,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone )

  call comp_submatrix_mod( ngrid_r,grid_r, &
              grid_rho,grid_kappa,grid_mu, &
              grid_Ak,grid_Am,grid_L,grid_N,grid_Fk,grid_Fm, &
              submatrix_I0,submatrix_I1k,submatrix_I3k, &
              submatrix_I3m,submatrix_I4,submatrix_I5k, &
              submatrix_I5m,submatrix_I6,submatrix_I7, &
              submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod,xi_Gauss,weight_Gauss,n_structure_zone,boundary_r,rho_structure_zone,vpv_structure_zone,vph_structure_zone,vsv_structure_zone,vsh_structure_zone,eta_structure_zone )

! **********************************************************************
! compute wavefield for each frequency (for complex omega)
! **********************************************************************

! **********************************************************************
! loop on the frequencies
! **********************************************************************

  do i_frequency=i_frequency_min(i_frequency_band, n_frequency_band, n_frequency) +rank, &
                             i_frequency_max(i_frequency_band, n_frequency_band, n_frequency), size

  write(*,*) i_frequency,i_frequency_band,omega(time_series_length,i_frequency,omega_imag)

   if ( i_frequency/=0 ) then

    if ( i_frequency>maxnfreq ) stop 'error: i_frequency is too large.'

    call init_complex_array( 3*max_nstation,station_displacement(1,1,i_frequency) )

    do l=lmin0, lmax(i_frequency_band,n_frequency_band,lmax0 )

! **********************************************************************
! compute the displacement and the traction at the boundary
! of the heterogeneous region.
! **********************************************************************
      if ( l == 0 ) then

        call comp_excitation0( maxngrid_r, &
                  omega(time_series_length,i_frequency,omega_imag), &
                  l,source_r,source_mt, &
                  igrid_rs,grid_Cks,grid_Cms,grid_Fks,grid_Fms, &
                  grid_qkappas,grid_qmus, &
                  submatrix_I0,submatrix_I1k,submatrix_I1m, &
                  submatrix_I3k,submatrix_I3m, &
                  submatrix_I5k,submatrix_I5m, &
                  submatrix_I6,submatrix_I7, &
                  idim_rs_sph0, &
                  whole_vector_sph,whole_vector_tor )

        call comp_wavefield0( maxngrid_r, &
                  omega(time_series_length,i_frequency,omega_imag), &
                  submatrix_I0,submatrix_I1k,submatrix_I1m, &
                  submatrix_I2,submatrix_I3k,submatrix_I3m, &
                  submatrix_I5k,submatrix_I5m, &
                  submatrix_I6,submatrix_I7, &
                  grid_r,grid_mu,grid_qkappa,grid_qmu,l, &
                  idim1_sph,idim2_sph, &
                  idim0,init_npos_sph,init_npos_tor, &
                  idim_rs_sph0, &
                  idim_station_sph0, &
                  whole_matrix_sph, &
                  whole_matrix_dr_sph, &
                  whole_vector_sph,work_vector )

        call check_amp_significance( maxngrid_r,l,amp_max,i_significance,idim0, &
                  init_npos_sph,init_npos_tor, &
                  idim1_sph,idim2_sph,idim1_tor,idim2_tor, &
                  grid_mu,idim_station_sph0,idim_station_tor0, &
                  whole_vector_sph,whole_vector_tor,work_vector )

        call comp_displacement_station0( maxngrid_r,maxlmax, &
                  whole_vector_sph, &
                  l,n_station, &
                  idim_station_sph0, &
                  vecsph_sph1, &
                  station_displacement(1,1,i_frequency) )
      else
        if ( i_significance==1 ) then
        call comp_excitation( maxngrid_r, &
                  omega(time_series_length,i_frequency,omega_imag), &
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

        call comp_wavefield( maxngrid_r, &
                  omega(time_series_length,i_frequency,omega_imag), &
                  submatrix_I0,submatrix_I1k,submatrix_I1m, &
                  submatrix_I2,submatrix_I3k,submatrix_I3m, &
                  submatrix_I4,submatrix_I5k,submatrix_I5m, &
                  submatrix_I6,submatrix_I7, &
                  submatrix_I3k_mod,submatrix_I3m_mod,submatrix_I4_mod, &
                  grid_r,grid_mu,grid_qkappa,grid_qmu,l, &
                  idim1_sph,idim2_sph,idim1_tor,idim2_tor, &
                  idim0,init_npos_sph,init_npos_tor, &
                  idim_rs_sph,idim_rs_tor, &
                  idim_station_sph,idim_station_tor, &
                  whole_matrix_sph,whole_matrix_tor, &
                  whole_matrix_dr_sph,whole_matrix_dr_tor, &
                  whole_vector_sph,whole_vector_tor,work_vector )

        call check_amp_significance( maxngrid_r,l,amp_max,i_significance,idim0, &
                  init_npos_sph,init_npos_tor, &
                  idim1_sph,idim2_sph,idim1_tor,idim2_tor, &
                  grid_mu,idim_station_sph,idim_station_tor, &
                  whole_vector_sph,whole_vector_tor,work_vector )

        call comp_displacement_station( maxngrid_r,maxlmax, &
                  whole_vector_tor,whole_vector_sph, &
                  l,n_station, &
                  idim_station_sph,idim_station_tor, &
                  vecsph_sph1,vecsph_sph2,vecsph_tor, &
                  station_displacement(1,1,i_frequency) )
        endif
      endif
    enddo ! of do l=lmin0, lmax(i_frequency_band,n_frequency_band,lmax0 )
   endif

! **********************************************************************
! write the result to the spectrum files
! **********************************************************************
!   call write_sac_file( i_frequency,n_station,sac_file,station_displacement )

  enddo ! of do i_frequency=i_frequency_min(i_frequency_band, n_frequency_band, n_frequency) +rank

 enddo ! of do i_frequency_band=1,n_frequency_band

#ifdef USE_MPI
  nsnd = 2 * 3 * max_nstation

  if ( rank /= 0 ) then

    ipos = 0
    dpos = 0

    do i_frequency_band=1,n_frequency_band
      do i_frequency=i_frequency_min(i_frequency_band, n_frequency_band, n_frequency) +rank, &
                               i_frequency_max(i_frequency_band, n_frequency_band, n_frequency), size
      call mputd( station_displacement(1,1,i_frequency), nsnd,info ,mbfd,dpos,dbuff)
      enddo
    enddo

    call msend( 0,20,info ,mbfi,mbfd,size,ipos,dpos,param,ibuff,dbuff)
  else
    do irank=1,size-1
      call mrecv( 20,irank,info ,mbfi,mbfd,ipos,dpos,param,ibuff,dbuff)
      do i_frequency_band=1,n_frequency_band
        do i_frequency=i_frequency_min(i_frequency_band, n_frequency_band, n_frequency) +irank, &
                               i_frequency_max(i_frequency_band, n_frequency_band, n_frequency), size
        call mgetd( station_displacement(1,1,i_frequency), nsnd,info ,mbfd,dpos,dbuff)

        enddo
      enddo
    enddo
  endif
#endif

  if ( rank == 0 ) then
    call convspc( max_nstation,maxnfreq,n_station, time_series_length,n_frequency,omega_imag, &
           station_displacement, source_depth,source_lat,source_lon,station_lat,station_lon, &
           work_spc,work_spc_real,work_time,sac_file,ITYPE_SEISMOGRAMS )
    
!    coutfile = "./data/test_velocity"

!    do i = 0,n_frequency
!       displacementsngl(:,:) = cmplx(0.e0)
!       if (i /= 0) then
!         do j = 1,n_station
!          displacementsngl(1:3,j) =  station_displacement(1:3,j,i)
!         end do
!       endif

!      open(1,file=coutfile,status='unknown',form='unformatted',access = 'direct', recl=2*3*kind(0e0)*n_station)
!      write(1,rec=i+1)displacementsngl(1:3,1:n_station)
!      close(1)
!     enddo


  endif

#ifdef USE_MPI
  call mpi_finalize( info )
#endif

! functions
  contains

  complex(kind=8) function omega(time_series_length,i_frequency,omega_imag)
    implicit none
    integer :: i_frequency
    real(kind=8) :: time_series_length,omega_imag
    omega = cmplx( 2.d0*pi*dble(i_frequency)/dble(time_series_length), -omega_imag )
  end function omega

  integer function ngrid_r_estimated(i_frequency_band,n_frequency_band,ngrid_r0)
    implicit none
    integer :: i_frequency_band,n_frequency_band,ngrid_r0
    ngrid_r_estimated = int( dnint ( ngrid_r0 * dble(i_frequency_band) / dble(n_frequency_band) ) )
  end function ngrid_r_estimated

  integer function lmax(i_frequency_band,n_frequency_band,lmax0)
    implicit none
    integer :: i_frequency_band,n_frequency_band,lmax0
    lmax = int( dnint ( lmax0 * dble(i_frequency_band) / dble(n_frequency_band) ) )
  end function lmax

  integer function i_frequency_min(i_frequency_band,n_frequency_band,n_frequency)
    implicit none
    integer :: i_frequency_band,n_frequency_band,n_frequency
    i_frequency_min = int( dnint ( n_frequency &
                  * dble(i_frequency_band-1) / dble(n_frequency_band) ) ) + 1 - int(1.d0/dble(i_frequency_band))
  end function i_frequency_min

  integer function i_frequency_max(i_frequency_band,n_frequency_band,n_frequency)
    implicit none
    integer :: i_frequency_band,n_frequency_band,n_frequency
    i_frequency_max = int( dnint ( n_frequency * dble(i_frequency_band) / dble(n_frequency_band) ) )
  end function i_frequency_max

  end program dsmti

