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
module data_bkgrdmodel

  use global_parameters, only: dp, sp
  implicit none

  integer                    :: ndisc, nfluidregions
  integer, allocatable       :: idom_fluid(:)
  real(kind=dp), allocatable :: discont(:)
  real(kind=dp), allocatable :: vp(:,:), vs(:,:), rho(:,:)
  logical, allocatable       :: solid_domain(:)
  integer                    :: lfbkgrdmodel
  character(len=100)         :: bkgrdmodel
  logical                    :: have_fluid, have_solid
  real(kind=dp)              :: pts_wavelngth
  real(kind=dp)              :: period, courant
  real(kind=dp)              :: dt
  integer                    :: nc_init, nthetaslices, nradialslices

  ! the sole quantities to be created in create_subregions
  ! that are needed by the rest of the mesher
  integer                    :: nz_glob, ns_glob, nc_glob
  integer, allocatable       :: iclev_glob(:)
  real(kind=dp), allocatable :: dz_glob(:)
  real(kind=dp)              :: rmin, minh_ic, maxh_ic, maxh_icb
  real(kind=dp)              :: minhvp, maxhvs, maxhnsicb

  ! The following variables are only needed by external models
  character(len=100)         :: fnam_ext_model

end module data_bkgrdmodel
!=========================================================================================
