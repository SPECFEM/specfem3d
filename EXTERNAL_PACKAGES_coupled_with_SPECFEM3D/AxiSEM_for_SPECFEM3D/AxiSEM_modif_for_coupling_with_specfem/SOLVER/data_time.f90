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
!> Various variables around timing
module data_time

  use global_parameters
  implicit none
  public

  real(kind=dp)       :: enforced_dt        ! < Enforced time step in inparam
  real(kind=dp)       :: enforced_period    ! < Enforced source period in inparam
  character(len=8)    :: time_scheme        ! < time extrapolation scheme, better be:
                                            !! symplec4, sympqua4, symp_5_4, newmark2
  real(kind=dp)       :: period             ! < Dominant source period,  given by mesher
  real(kind=dp)       :: courant            ! < Courant number,  given by mesher
  real(kind=dp)       :: t                  ! < current time
  real(kind=dp)       :: deltat             ! < Time step size
  real(kind=dp)       :: half_dt            ! < Half of time step
  real(kind=dp)       :: half_dt_sq         ! < Half of squared time step
  real(kind=dp)       :: deltat_strain      ! < Time step for strain/kernel dumps
  real(kind=dp)       :: deltat_coarse      ! < Time step of slowest dump
  real(kind=dp)       :: seislength_t       ! < seismogram length in seconds
  integer             :: niter              ! < Number of iterations
  integer             :: iclockold, idold   ! < tick labels for timer
  integer             :: iclockcomm, idcomm ! < tick labels for comm timer
  integer             :: iclockmpi, idmpi   ! < tick labels for MPI timer
  integer             :: iclockmpiws, idmpiws   ! < tick labels for solid MPI wait timer
  integer             :: iclockmpiwf, idmpiwf   ! < tick labels for fluid MPI wait timer
  integer             :: iclockanelst, idanelst ! < tick labels for anelastic stiffness timer
  integer             :: iclockanelts, idanelts ! < tick labels for anelastic time step timer
  integer             :: iclockstiff, idstiff   ! < tick labels for stiffness timer
  integer             :: iclockdump, iddump ! < tick labels for dump timer
  integer             :: iclocknbio, idnbio ! < tick labels for non blocking IO timer
  real(kind=dp)       :: seis_dt            ! < seismogram sampling rate in seconds
  integer             :: seis_it            ! < seismogram sampling rate in time steps
  real(kind=dp)       :: snap_dt            ! < time interval between snaps in seconds
  integer             :: snap_it            ! < equivalent as snap_dt in time steps
  integer             :: strain_it          ! < strain/kernel dump interval in time steps
  integer             :: check_it           ! < Checkpointing of seismograms in NetCDF output
  integer             :: nstages            ! < number of substages in symplectic schemes
  real(kind=realkind) :: decay, shift_fact

end module data_time
!=========================================================================================
