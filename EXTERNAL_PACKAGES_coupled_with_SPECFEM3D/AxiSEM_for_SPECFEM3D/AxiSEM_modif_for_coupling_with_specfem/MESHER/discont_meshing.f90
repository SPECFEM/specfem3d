!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
module discont_meshing

  use global_parameters
  use data_bkgrdmodel
  use data_spec
  use model_discontinuities, only: define_discont
  use background_models
  use data_diag

  implicit none

  public :: create_subregions
  private

  contains

!-----------------------------------------------------------------------------------------
subroutine create_subregions

  real(kind=dp), allocatable                 :: rdisc_top(:), rdisc_bot(:)
  real(kind=dp), dimension(:), allocatable   :: ds_glob, radius_arr
  real(kind=dp), dimension(:), allocatable   :: vp_arr, vs_arr
  real(kind=dp)     :: ds, minh, maxh, aveh, dz, current_radius
  integer           :: idom, ic, icount_glob, iz_glob
  integer           :: ns_ref,  ns_ref_icb1, ns_ref_icb, ns_ref_surf
  logical           :: memorydz, current
  integer, dimension(1)             :: iloc1, iloc2
  real(kind=realkind), dimension(1) :: rad1, rad2
  character(len=32) :: fmtstring

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! constraining factors/relations for mesh architecture:
  ! 1) number of coarsening levels:
  !    ncoars = log( vs(cmb)/vs(surf)*r(surf)/r(cmb) ) / log(2) = 3.2 (for PREM)
  ! 2) period/resolution: nlambda = period * min(vs/h)
  ! 3) time step:  dt < courant * min(h/vp)
  ! 3),4) == => keep
  ! 4) resulting number of lateral elements ns(r) = nlambda * r / ( vp * period)
  ! = => first "derive" global ns(r), ds(r) and from that define nz(r),dz(r)
  !     for each subdomain such that dz(r) varies less than ds(r)
  ! 5) additionally: number of linear central elements
  ! 6) number of lateral processor domains (cakepieces)
  ! 7) coarsening levels need to be between discontinuities
  ! = => total number of subdomains: ndisc+ncoars
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Define background model discontinuities (e.g. PREM)
  ! & model specific boolean "solid" (fluid layers?) for each subdomain

   have_fluid = .false.
   have_solid = .true.
   nfluidregions = 0

   call define_discont

   allocate(rdisc_top(ndisc), rdisc_bot(ndisc), solid_domain(ndisc))
   allocate(idom_fluid(ndisc))
   idom_fluid=0

   solid_domain(1:ndisc) = .true.
   do idom = 1, ndisc
      if (idom < ndisc) then
         rdisc_top(idom) = discont(idom)
         rdisc_bot(idom) = discont(idom+1)
      else
         rdisc_top(idom) = discont(idom)
         rdisc_bot(idom) = zero
      endif

      if (vs(idom,1) == 0.0d0 .and. vs(idom,2) == 0.0d0) then
         nfluidregions = nfluidregions + 1
         have_fluid = .true.
         solid_domain(idom) = .false.
         idom_fluid(nfluidregions) = idom
      else if ( vs(idom,1) == 0.0d0 .and. vs(idom,2) /= 0.0d0 .or. &
              vs(idom,1) /= 0.0d0 .and. vs(idom,2) == 0.0d0) then
         write(*,*) 'ERROR in background model region:'
         write(*,*) 'Cannot have one region with fluid and solid parts...'
         write(*,*) 'upper radius/vs:', rdisc_top(idom), vs(idom,1)
         write(*,*) 'lower radius/vs:', rdisc_bot(idom), vs(idom,2)
      endif

      write(*,*) '#######################################################################'
      fmtstring = '("  ", A, I12, F12.2)'
      write(*,fmtstring)'discontinuities:    ', idom,real(discont(idom))
      fmtstring = '("  ", A, L12, I12)'
      write(*,fmtstring)'solid/fluid domain: ', solid_domain(idom),idom_fluid(idom)
      fmtstring = '("  ", A, F12.2, F12.2)'
      write(*,fmtstring)'upper/lower radius: ', real(rdisc_top(idom)),real(rdisc_bot(idom))
      write(*,fmtstring)'vs jump:            ', real(vs(idom,1)),real(vs(idom,2))
      write(*,*) '#######################################################################'
      write(*,*)

   enddo

   if (nfluidregions == ndisc) then
      write(*,*) 'COMPLETELY acoustic domain!'
      have_solid = .false.
   endif

   write(*,"(10x,'Number of discontinuities/regions:     ',i3)") ndisc
   write(*,"(10x,'Number of fluid regions:               ',i3)") nfluidregions

   print *
   write(*,*)'Constructing the mesh....'
   print *

  ! Loop over discontinuities/subregions
  icount_glob = 0   ! element depth levels
  ic = 0            ! coarsening levels

  ! nc expected for PREM :
  ! nc = int( log( (r_surface/min_velocity_surface (S))/ &
  !                 (r_icb/min_velocity_icb (P)) ) )

  ! take surface/crustal values as constraint on ns resolution
  if (solid_domain(1)) then
     ns_ref_surf = estimate_ns(pts_wavelngth, rdisc_top(1), vs(1,1), period)
  else ! top layer is a fluid
     ns_ref_surf = estimate_ns(pts_wavelngth, rdisc_top(1), vp(1,1), period)
  endif

  if (dump_mesh_info_screen) &
        write(*,*) 'ns_ref initial estimate from crust   :', ns_ref_surf

  ! take ICB value as constraint on ns resolution
  if (solid_domain(ndisc)) then
     ns_ref_icb1 = estimate_ns(pts_wavelngth, rdisc_top(ndisc), vs(ndisc,1), period) &
                        * 2 ** nc_init
  else ! bottom layer is a fluid
     ns_ref_icb1 = estimate_ns(pts_wavelngth, rdisc_top(ndisc), vp(ndisc,1), period) &
                        * 2 ** nc_init
  endif

  if (dump_mesh_info_screen) &
        write(*,*) 'ns_ref initial estimate from icb     :', ns_ref_icb1

  ! need to resolve the crust and icb
  ns_ref = max(ns_ref_icb1, ns_ref_surf)

  ! need to make sure that ns_ref is defined such that ns is even below
  ! last coarsening level (fix to make sure nel and lnodes counts are correct)
  if ( mod(ns_ref, 2 ** nc_init * nthetaslices*2) /= 0 ) &
        ns_ref = (2 ** nc_init * nthetaslices * 2) &
                    * (ns_ref / (2 ** nc_init * nthetaslices * 2) + 1)

  if (dump_mesh_info_screen) &
        write(*,*) 'ns_ref fixed with procs & coarsenings:', ns_ref

  ns_glob = ns_ref

  ! DETERMINE rmin such that elements in central region are not too large.
  ! Assumptions taken here:
  !  - velocities in lowermost layer very low such that they constrain the
  !         amount of lateral elements needed at the surface
  !  - almost constant, non-decreasing velocities in the inner core

  if (solid_domain(ndisc)) then
     maxh_icb = period * vs(ndisc,1) / (pts_wavelngth * max_spacing(npol))
  else
     maxh_icb = period * vp(ndisc,1) / (pts_wavelngth * max_spacing(npol))
  endif

  !rmin = maxh_icb * (dble(ns_ref / dble(2. * dble(2**nc_init))) + 1.d0)
  rmin = maxh_icb * (ns_ref / (2.d0 * 2**nc_init) + 1)

  if (dump_mesh_info_screen) &
        write(*,*) 'actual ds at innermost discontinuity [km] :', &
                    0.5 * pi * rdisc_top(ndisc) / real(ns_ref) * real(2**nc_init) / 1000.
  if (dump_mesh_info_screen) &
        write(*,*) 'maximal dz at innermost discontinuity [km]:', maxh_icb / 1000.

  if (rmin > rdisc_top(ndisc) - maxh_icb * .9) then
     ! at least the ICB....
     rmin = rdisc_top(ndisc) - maxh_icb * .9
  endif

  if (dump_mesh_info_screen) then
     write(*,*) 'CALCULATED RMIN=', rmin
     write(*,*) 'MAXH_ICB=', maxh_icb
     write(*,*) '# central region elements (incl. buffer, e.g. along axis):', &
            int(ns_ref / (2.* 2**nc_init) + 1.)
  endif

  rdisc_bot(ndisc) = rmin

  ! trial loop to calculate global amount of radial layers icount_glob and
  ! coarsening levels ic
  do idom =1, ndisc
     current_radius = rdisc_top(idom)
     memorydz = .false.
     do while (current_radius > rdisc_bot(idom) )
        call compute_dz_nz(idom, rdisc_bot, current_radius, dz, ds, current, memorydz, &
                           icount_glob, ic, ns_ref)
     enddo
  enddo

  write(*,*)
  write(*,"(10x,'Spherical shell part of the mesh:')")
  write(*,"(10x,'Total number of depth levels:        ',i6)") icount_glob
  write(*,"(10x,'Actual number of coarsening levels:     ',i3)") ic
  write(*,"(10x,'Anticipated number of coarsening levels:',i3)") nc_init

  !TODO: MvD: can't this be automatized?
  if (nc_init /= ic) then
     write(*,*) ' a bit of rethinking is needed for nc_init!', &
                  'Check your calculus'
     stop
  endif

  ! assign global values after counting 'em
  nz_glob = icount_glob
  nc_glob = ic
  allocate(dz_glob(nz_glob))
  allocate(iclev_glob(1:nc_glob))
  allocate(ds_glob(nz_glob), radius_arr(nz_glob))
  allocate(vp_arr(nz_glob), vs_arr(nz_glob))

  ! Final output needed to feed the skeleton (data_bkgrdmodel.f90):
  ! ns_glob, nz_glob counting all depth levels
  ! dz_glob(1:nz_glob)
  ! nc_glob number of coarsening levels
  ! iclev_glob(1:nc_glob)

  ns_ref = max(ns_ref_icb1, ns_ref_surf)

  ! need to make sure that ns_ref is defined such that ns is even below
  ! last coarsening level (fix to make sure nel and lnodes counts are correct)
  if (mod(ns_ref, 2**nc_init*nthetaslices*2) /= 0 ) &
        ns_ref = (2**nc_init*nthetaslices*2)* ( ns_ref/(2**nc_init*nthetaslices*2) + 1 )

  if (dump_mesh_info_screen) &
        write(*,*) 'ns_ref fixed with # processors/coarsenings:' , ns_ref

  ns_glob = ns_ref

  ndisc = ndisc
  icount_glob = 0
  ic = 0

  do idom=1, ndisc
     current_radius = rdisc_top(idom)
     memorydz = .false.
     do while (current_radius > rdisc_bot(idom))
        call compute_dz_nz(idom, rdisc_bot, current_radius, dz, ds, current, memorydz, &
                           icount_glob, ic, ns_ref)
        ! Storing radial info into global arrays
        if (current) iclev_glob(ic) = nz_glob - icount_glob + 1
        dz_glob(icount_glob) = dz
        ds_glob(icount_glob) = ds
        radius_arr(icount_glob) = current_radius
        vp_arr(icount_glob) = velocity(current_radius, 'v_p', idom, bkgrdmodel, &
                                       lfbkgrdmodel)
        vs_arr(icount_glob) = velocity(current_radius, 'v_s', idom, bkgrdmodel, &
                                       lfbkgrdmodel)
        if (vs_arr(icount_glob) < 0.1d0 * vs_arr(1)) &
                vs_arr(icount_glob) = vp_arr(icount_glob)
     enddo
  enddo

  ! scale for smallest grid spacing in GLL clustering
  call gll_spacing(npol,minh,maxh,aveh)

  dt = courant * min(minval(ds_glob/vp_arr), minval(dz_glob/vp_arr)) &
            / dble(npol) * minh/ aveh

  if (dump_mesh_info_screen) write(*,*) 'TIME STEP:', dt

  if (dump_mesh_info_files) then
     open(unit=667,file=diagpath(1:lfdiag)//'/period_courant_pts_wavelength.txt')
     write(667,11)period,courant,pts_wavelngth,dt
     write(667,11)minh,maxh,aveh,real(npol)
     open(unit=668,file=diagpath(1:lfdiag)//'/ds_dz.txt')
     iloc1=minloc(ds_glob/vp_arr) ; rad1=radius_arr(iloc1)
     iloc2=minloc(dz_glob/vp_arr) ; rad2=radius_arr(iloc2)
     write(668,11)minval(ds_glob/vp_arr),rad1, minval(dz_glob/vp_arr),rad2

    ! %%%%%%%%%%% for MATLAB %%%%%%%%%%%%%
     open(unit=666,file=diagpath(1:lfdiag)//'/ds_dz_matlab.txt')
     do iz_glob = 1, nz_glob
        write(666,11) radius_arr(iz_glob),dz_glob(iz_glob),ds_glob(iz_glob), &
             vs_arr(iz_glob)*period,vp_arr(iz_glob)/min(ds_glob(iz_glob), &
             dz_glob(iz_glob))*dt*real(npol)/minh*aveh
     enddo
11   format(40(1pe12.5,2x))
    ! %%%%%%%%%%% end MATLAB %%%%%%%%%%%%%

     close(666)
     close(667)
     close(668)
  endif

  if (dump_mesh_info_files) then
     open(unit=30,file=diagpath(1:lfdiag)//'/coarsening_radii.dat')
     do ic=1,nc_glob
        write(30,*)iclev_glob(ic),radius_arr(iclev_glob(ic))
     enddo
     close(30)
  endif

  ! INNER CORE

  ns_ref_icb = ns_ref
  minh_ic = vp(ndisc,1) * min(minval(ds_glob / vp_arr), minval(dz_glob / vp_arr)) &
              * minh / aveh

  if (solid_domain(ndisc)) then
     ns_ref = estimate_ns(pts_wavelngth, discont(ndisc), vs(ndisc,1), period)
     maxh_ic = vs(ndisc,1) * max(maxval(ds_glob / vs_arr), maxval(dz_glob / vs_arr)) &
                    * maxh / aveh
  else
     ns_ref = estimate_ns(pts_wavelngth, discont(ndisc), vp(ndisc,1), period)
     maxh_ic = vp(ndisc,1) * max(maxval(ds_glob / vs_arr), maxval(dz_glob / vs_arr)) &
                    * maxh / aveh
  endif

  if (dump_mesh_info_screen) then
     write(*,*) ' NS_REF at ICB + from meshing ', ns_ref_icb
     write(*,*) ' WHAT WE WANT ', ns_ref
     write(*,*)
     write(*,*) 'MESH EFFICIENCY: smallest/largest spacing etc. '
     write(*,*) 'min (h_el/vp):   ', &
          min(minval(ds_glob/vp_arr),minval(dz_glob/vp_arr))*minh/aveh
     write(*,*) 'dt/courant*npol: ',dt/courant*real(npol)
     write(*,*)
     write(*,*) 'max (h_el/vs):         ', &
          max(maxval(ds_glob/vs_arr),maxval(dz_glob/vs_arr))*maxh/aveh
     write(*,*) 'period/pts_Wavelength: ',period/pts_wavelngth
     write(*,*)
  endif

  if (dump_mesh_info_files) then
     open(30,file=diagpath(1:lfdiag)//'/test_innercore_hminmax.dat')
     write(30,*) minh_ic, maxh_ic
     close(30)
  endif

  if (dump_mesh_info_screen) then
     write(*,*) 'Inner Core element sizes:'
     write(*,*) 'r_min=', rmin
     write(*,*) 'max h r/ns(icb):', pi * 0.5 * discont(ndisc) / real(ns_ref_icb)
     write(*,*) 'precalculated max h:', maxh_icb
     write(*,*) 'max h:', maxh_ic
     write(*,*) 'min h:', minh_ic
  endif

  ! for mesh_params.h
  minhvp = min(minval(ds_glob / vp_arr), minval(dz_glob / vp_arr)) * minh / aveh
  maxhvs = max(maxval(ds_glob / vs_arr), maxval(dz_glob / vs_arr)) * maxh / aveh
  maxhnsicb = pi * 0.5d0 * discont(ndisc) / dble(ns_ref_icb)

end subroutine create_subregions
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_dz_nz(idom, rdisc_bot, current_radius, dz, ds, current, memorydz, &
                         icount_glob, ic, ns_ref)

  use data_grid, only: fluidfac

  integer, intent(in)           :: idom
  real(kind=dp), intent(in)     :: rdisc_bot(ndisc)
  real(kind=dp), intent(inout)  :: current_radius, dz, ds
  logical, intent(inout)        :: current, memorydz
  integer, intent(inout)        :: icount_glob, ic
  integer, intent(inout)        :: ns_ref

  real(kind=dp)                 :: dz_trial
  real(kind=dp)                 :: velo
  integer                       :: nz_trial,ns_trial

  current = .false.

  if (solid_domain(idom)) then
     velo = velocity(current_radius,'v_s',idom,bkgrdmodel,lfbkgrdmodel)
  else
     ! add a prefactor here to make smaller elements in outer core and
     ! reduce dispersion error in outer core!
     ! Outer Core P-Waves see a lot more dispersion error than in the
     ! mantle, due to beeing on the resolution edge. But as the fluid is
     ! really cheap, we could just make the elements smaller...
     velo = fluidfac * velocity(current_radius,'v_p',idom,bkgrdmodel,lfbkgrdmodel)
  endif
  ns_trial = estimate_ns(pts_wavelngth,current_radius,velo,period)
  icount_glob = icount_glob+1

  if (ns_trial < ns_ref/2 .and. (ic < nc_init) ) then
     ! Is coarsening possible within this subregion
     ! ( <- are there at least two elemental
     ! layers between the actual layer and the bottom of the subregion?)

     dz_trial = .5d0* pi * current_radius / dble(ns_ref)
     nz_trial = max(ceiling((current_radius-rdisc_bot(idom))/dz_trial),1)
     if (nz_trial >= 3) then
        ns_trial = ns_ref
        ns_ref = ns_ref / 2
        ic = ic + 1
        memorydz = .true.
        current = .true.
     endif
  endif
  dz_trial = .5d0* pi * current_radius / dble(ns_trial)
  if (memorydz .and. .not. current) then
     dz_trial = .5d0* pi * current_radius / dble(2*ns_trial)
     memorydz = .false.
  endif
  nz_trial = max(ceiling((current_radius-rdisc_bot(idom))/dz_trial),1)
  dz = (current_radius-rdisc_bot(idom))/dble(nz_trial)
  ds = .5d0* pi * current_radius / dble(ns_ref)
  current_radius = current_radius -dz

end subroutine compute_dz_nz
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function estimate_ns(el_per_lambda,r,v,period)

  real(kind=dp) ,intent(in) :: el_per_lambda,r
  real(kind=dp)             :: v, period
  real(kind=dp)             :: minh, maxh, aveh

  ! scaling for irregular GLL spacing
  call gll_spacing(npol, minh, maxh, aveh)
  estimate_ns=ceiling(el_per_lambda * .5d0 * pi * r / (v * period) * maxh / aveh)

end function estimate_ns
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine spacing_info(npol)

  use splib

  integer, intent(in) :: npol
  real(kind=dp)       :: minh, maxh, aveh
  real(kind=dp),allocatable,dimension(:) :: eta, xi_k, dxi
  real(kind=dp),allocatable,dimension(:) :: spacing_eta, spacing_xi

  integer :: i

  allocate(eta(0:npol), xi_k(0:npol), dxi(0:npol))
  allocate(spacing_eta(1:npol), spacing_xi(1:npol))

  call ZELEGL(npol, eta, dxi)
  call zemngl2(npol, xi_k)

  ! spacing within [0,1]
  do i=0, npol-1
     spacing_eta(i+1) = dabs(eta(i) - eta(i+1)) / 2.d0
     spacing_xi(i+1) = dabs(xi_k(i) - xi_k(i+1)) / 2.d0
  enddo

  minh = min(minval(spacing_eta), minval(spacing_xi))
  maxh = max(maxval(spacing_eta), maxval(spacing_xi))
  aveh = 1.d0 / dble(npol)

  write(*,*) ' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '
  write(*,*) 'GLL SPACING for polynomial order', npol
  write(*,*) '  AVERAGE SPACING:       ', aveh
  write(*,*) '  MINIMAL SPACING:       ', minh
  write(*,*) '  MAXIMAL SPACING:       ', maxh
  write(*,*) '  MIN/AVE, MAX/AVE:      ', minh / aveh, maxh / aveh
  write(*,*) '  MIN GLL,GLJ SPACING:   ', minval(spacing_eta), minval(spacing_xi)
  write(*,*) '  MAX GLL,GLJ SPACING:   ', maxval(spacing_eta), maxval(spacing_xi)
  write(*,*) '  MINGLL/AVE,MINGLJ/AVE: ', minval(spacing_eta) / aveh, &
                                          minval(spacing_xi) / aveh
  write(*,*) '  MAXGLL/AVE,MAXGLJ/AVE: ', maxval(spacing_eta) / aveh, &
                                          maxval(spacing_xi) / aveh
  write(*,*) '  VALUES GLL, GLJ in [0,1]:'

  do i=0, npol
     write(*,*) i, eta(i), xi_k(i)
  enddo
  write(*,*)' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '

end subroutine spacing_info
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gll_spacing(npol, minh, maxh, aveh)

  use splib

  integer, intent(in)           :: npol
  real(kind=dp), intent(out)    :: minh, maxh, aveh
  real(kind=dp) :: spacing_eta(npol), spacing_xi(npol)
  real(kind=dp) :: eta(0:npol), xi_k(0:npol), dxi(0:npol)

  integer :: i

  call zelegl(npol, eta, dxi)
  call zemngl2(npol, xi_k)

  ! spacing within [0,1]
  do i=0, npol-1
     spacing_eta(i+1) = dabs(eta(i) - eta(i+1)) / 2.d0
     spacing_xi(i+1) = dabs(xi_k(i) - xi_k(i+1)) / 2.d0
  enddo

  minh = min(minval(spacing_eta), minval(spacing_xi))
  maxh = max(maxval(spacing_eta), maxval(spacing_xi))
  aveh = 1.d0 / dble(npol)

end subroutine gll_spacing
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function max_spacing(npol)

  use splib

  integer, intent(in) :: npol
  real(kind=dp) :: minh,maxh,aveh
  real(kind=dp) :: spacing_eta(npol), spacing_xi(npol)
  real(kind=dp) :: eta(0:npol), xi_k(0:npol), dxi(0:npol)

  integer :: i

  call ZELEGL(npol,eta,dxi)
  call zemngl2(npol,xi_k)

  ! spacing within [0,1]
  do i=0, npol-1
     spacing_eta(i+1) = dabs(eta(i) - eta(i+1)) / 2.d0
     spacing_xi(i+1) = dabs(xi_k(i) - xi_k(i+1)) / 2.d0
  enddo

  minh = min(minval(spacing_eta), minval(spacing_xi))
  maxh = max(maxval(spacing_eta), maxval(spacing_xi))
  aveh = 1.d0 / dble(npol)

  max_spacing = maxh / aveh
end function max_spacing
!-----------------------------------------------------------------------------------------

end module discont_meshing
!=========================================================================================
