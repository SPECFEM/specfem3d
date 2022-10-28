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
!> Miscellaneous variables relevant to any read/write process such as
!! paths, logicals describing what to save, sampling rate of dumps
module data_io

  use global_parameters
  implicit none
  public

  character(len=200) :: datapath,infopath
  integer           :: lfdata,lfinfo
  logical           :: dump_energy
  logical           :: dump_snaps_glob
  logical           :: dump_vtk
  logical           :: dump_xdmf
  logical           :: dump_wavefields
  logical           :: checkpointing
  logical           :: diagfiles ! < Write diagnostic files (seismograms at antipodes,
                                 !! list of surface elements, blabla), default: false

  logical           :: need_fluid_displ
  real(kind=dp)     :: strain_samp
  integer           :: nseismo  ! < Number of seismogram samples
  integer           :: nstrain  ! < Number of wavefield dumps for kernels
  integer           :: nsnap    ! < Number of wavefield snapshots for movies
  character(len=12) :: dump_type
  character(len=8)  :: rec_file_type
  logical           :: sum_seis, sum_fields
  logical           :: add_hetero, file_exists, use_netcdf
  character(len=6)  :: output_format  ! < netcdf or binary
  integer           :: ncid_out
  logical           :: do_anel
  integer           :: verbose

  integer           :: deflate_level  ! < Level of deflate compression in NetCDF. Only used
                                      !! for the XDMF visualization so far.

  ! indices to limit dumping to select contiguous range of GLL points:
  ! 0 <= ibeg <= iend <= npol
  ! For the time being: dump the same in xeta and eta directions
  integer           :: ibeg, iend
  integer           :: jbeg, jend
  ! ndumppts_el = (iend-ibeg+1) * (jend-jbeg+1)
  integer           :: ndumppts_el

  ! for xdmf dumps
  integer               :: i_n_xdmf, j_n_xdmf
  integer, allocatable  :: i_arr_xdmf(:), j_arr_xdmf(:)
  real(kind=dp)         :: xdmf_rmin, xdmf_rmax, xdmf_thetamin, xdmf_thetamax
  real(kind=dp)         :: kwf_rmin, kwf_rmax, kwf_thetamin, kwf_thetamax

  ! rotations
  real(kind=dp)    :: rot_mat(3,3),trans_rot_mat(3,3)
  real(kind=dp), allocatable, dimension(:,:) :: recfac

  character(len=80), dimension(:), allocatable :: fname_rec_seis
  character(len=80), dimension(:), allocatable :: fname_rec_velo

  ! Flag for coupling (SB) Don't know where to put
  logical :: coupling

contains

!-----------------------------------------------------------------------------------------
subroutine define_io_appendix(app, iproc)
  ! Defines the 4 digit character string appended to any
  ! data or io file related to process myid.

  implicit none
  integer, intent(in)           :: iproc
  character(len=4), intent(out) :: app

  write(app,"(I4.4)") iproc

end subroutine define_io_appendix
!-----------------------------------------------------------------------------------------

end module data_io
!=========================================================================================
