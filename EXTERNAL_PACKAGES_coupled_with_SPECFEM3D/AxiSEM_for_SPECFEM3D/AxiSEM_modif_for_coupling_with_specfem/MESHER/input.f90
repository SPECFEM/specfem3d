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
module input

  use data_grid
  use data_diag
  use data_coarse
  use data_bkgrdmodel

  implicit none

  public :: read_params
  private

  contains

!-----------------------------------------------------------------------------------------
subroutine read_params

  use global_parameters
  use data_mesh
  use data_spec
  use background_models, only: override_ext_q

  character(len=100)    :: keyword, keyvalue, line
  integer               :: iinparam_mesh = 500, ioerr

  keyword = ' '
  keyvalue = ' '
  line = ' '

  ! Default values
  bkgrdmodel      = 'UNDEFINED'
  override_ext_q  = 'none'
  nthetaslices    = -1
  nradialslices   = 1
  period          = -1
  dump_mesh_vtk   = .true.
  dump_1dmodel    = .true.
  nc_init         = 3
  npol            = 4
  pts_wavelngth   = 1.5
  courant         = 0.6
  router          = 6.371e6
  axisfac         = 0.7
  fluidfac        = 0.9
  dump_mesh_info_files = .false.
  dump_mesh_info_screen = .false.
  only_suggest_ntheta = .false.
  diagpath        = 'Diags'
  lfdiag          = index(diagpath,' ') - 1


  write(*, '(A)', advance='no') 'Reading inparam_mesh...'
  open(unit=iinparam_mesh, file='./inparam_mesh', status='old', action='read', &
       iostat=ioerr)
  if (ioerr /= 0) stop 'Check input file ''inparam_mesh''! Is it still there?'

  do
      read(iinparam_mesh,fmt='(a100)',iostat=ioerr) line
      if (ioerr < 0) exit
      if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

      read(line,*) keyword, keyvalue
      parameter_to_read : select case(trim(keyword))

      case('BACKGROUND_MODEL')
          bkgrdmodel = keyvalue
          lfbkgrdmodel = index(bkgrdmodel,' ') - 1

      case('EXT_MODEL')
          fnam_ext_model = keyvalue

      case('OVERRIDE_EXT_Q')
          override_ext_q = keyvalue

      case('DOMINANT_PERIOD')
          read(keyvalue, *) period

      case('NTHETA_SLICES')
          read(keyvalue, *) nthetaslices

      case('NRADIAL_SLICES')
          read(keyvalue, *) nradialslices

      case('ONLY_SUGGEST_NTHETA')
          read(keyvalue, *) only_suggest_ntheta

      case('WRITE_VTK')
          read(keyvalue, *) dump_mesh_vtk

      case('WRITE_1DMODEL')
          read(keyvalue, *) dump_1dmodel

      case('COARSENING_LAYERS')
          read(keyvalue, *) nc_init

      case('NPOL')
          read(keyvalue, *) npol

      case('EL_PER_LAMBDA')
          read(keyvalue, *) pts_wavelngth

      case('COURANT_NR')
          read(keyvalue, *) courant

      case('RADIUS')
          read(keyvalue, *) router

      case('AXIS_SHRINKING_FACTOR')
          read(keyvalue, *) axisfac

      case('FLUID_SHRINKING_FACTOR')
          read(keyvalue, *) fluidfac

      case('SAVE_DEBUG_FILES')
          read(keyvalue, *) dump_mesh_info_files

      case('VERBOSE')
          read(keyvalue, *) dump_mesh_info_screen

      end select parameter_to_read
  enddo

  if (trim(bkgrdmodel) == 'UNDEFINED') then
      write(*,20) 'BACKGROUND_MODEL'
      stop
  endif

  if (nthetaslices == -1) then
      write(*,20) 'NTHETA_SLICES'
      stop
  endif

  if (only_suggest_ntheta) nthetaslices = 4

  if (period == -1) then
      write(*,20) 'DOMINANT_PERIOD'
      stop
  endif
20 format('ERROR: Parameter ', A, ' not set in inparam_mesh')
  write(*,*) 'done'


  print *
  write(*,*) 'PREDEFINED MODEL/SIMULATION PARAMETERS'
  write(*,*) 'Background model                 : ',bkgrdmodel(1:lfbkgrdmodel)
  write(*,*) 'Dominant period [s]              : ',period
  write(*,*) 'Elements per dominant wavelength : ',pts_wavelngth
  write(*,*) 'Courant number                   : ',courant
  write(*,*) 'coarsening levels                : ',nc_init
  write(*,*) 'processors used in solver        : ',nthetaslices
  write(*,*) 'outer radius [m]                 : ',router
  write(*,*) 'save mesh info files?            : ',dump_mesh_info_files
  write(*,*) 'print mesh info to screen?       : ',dump_mesh_info_screen
  write(*,*) 'path to dump output files        : ',trim(diagpath)
  write(*,*)
  call flush(6)

end subroutine read_params
!-----------------------------------------------------------------------------------------

end module input
!=========================================================================================
