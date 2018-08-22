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
!> This module contains routines that compute and dump the respective meshes
!! underlying the actual wavefields to be dumped in the time loop
!! which is done in wavefields_io. This module needs pre-loop variables such
!! as mesh coordinates and is therefore cut off from the dumping module.
module meshes_io

  use global_parameters
  use data_mesh
  use data_proc
  use data_io

  use utlity, only: scoord, zcoord, rcoord, thetacoord

  implicit none

  private
  public :: finish_xdmf_xml
  public :: dump_wavefields_mesh_1d
  public :: dump_glob_grid_midpoint
  public :: dump_xdmf_grid
  public :: prepare_mesh_memoryvar_vtk
  public :: build_kwf_grid
  public :: dump_kwf_midpoint_xdmf
  public :: dump_kwf_fem_xdmf
  public :: dump_kwf_sem_xdmf
  public :: dump_kwf_gll_xdmf

contains

!-----------------------------------------------------------------------------------------
!> Dumps the mesh (s,z) [m] in ASCII format as needed to visualize global
!! snapshots, and additionally the constant factors preceding the displacement
!! in the fluid, namely rho^{-1} and (rho s)^{-1}.
!! When reading the fluid wavefield, one therefore needs to multiply all
!! components with inv_rho_fluid and the phi component with one/scoord!
!! Convention for order in the file: First the fluid, then the solid domain.
subroutine dump_glob_grid_midpoint(ibeg,iend,jbeg,jend)

  use data_pointwise, only: inv_rho_fluid
  use data_mesh, only: npol, nel_fluid, nel_solid

  integer, intent(in) :: ibeg,iend,jbeg,jend
  integer             :: iel,ipol,jpol

  open(unit=25000+mynum,file=datapath(1:lfdata)//'/glob_grid_'&
                             //appmynum//'.dat')

  if (have_fluid) then
     open(unit=26000+mynum,file=datapath(1:lfdata)// &
                                '/inv_rho_s_fluid_globsnaps_' &
                                //appmynum//'.dat')
     do iel=1,nel_fluid
        do jpol=0,npol,npol/2
           do ipol=0,npol,npol/2

           if ( axis_fluid(iel) .and. ipol == 0 ) then
               ! Axis s=0! write 1 instead of 1/s and then multiply
               ! with the correct factor dsdchi, obtained by L'Hospital's rule
               ! (see routine glob_snapshot in wavefields_io).
               write(26000+mynum,*)inv_rho_fluid(ipol,jpol,iel),one
           else
               write(26000+mynum,*)inv_rho_fluid(ipol,jpol,iel), &
                                   one/scoord(ipol,jpol,ielfluid(iel))
           endif

           write(25000+mynum,*)scoord(ipol,jpol,ielfluid(iel)), &
                               zcoord(ipol,jpol,ielfluid(iel))
           enddo
        enddo
     enddo
     close(26000+mynum)
  endif ! have_fluid

  do iel=1,nel_solid
     do jpol=0,npol,npol/2
        do ipol=0,npol,npol/2
           write(25000+mynum,*)scoord(ipol,jpol,ielsolid(iel)), &
                               zcoord(ipol,jpol,ielsolid(iel))
           enddo
        enddo
  enddo
  close(25000+mynum)

end subroutine dump_glob_grid_midpoint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_xdmf_grid()

  use nc_snapshots, only: nc_dump_snap_points, nc_dump_snap_grid, nc_make_snapfile
  use data_mesh

  integer               :: iel, ipol, jpol, ipol1, jpol1, i, j, ct, ipt, idest
  real(sp), allocatable :: points(:,:)
  integer, allocatable  :: grid(:,:), mapping(:)
  logical, allocatable  :: check(:), mask_tp_elem(:)
  character(len=120)    :: fname

  allocate(mask_tp_elem(nelem))
  mask_tp_elem = .false.

  ct = 0

  do iel=1, nel_fluid
      if (min(min(rcoord(0,0,ielfluid(iel)), rcoord(0,npol,ielfluid(iel))), &
              min(rcoord(npol,0,ielfluid(iel)), rcoord(npol,npol,ielfluid(iel)))) < xdmf_rmax &
          .and. &
          max(max(rcoord(0,0,ielfluid(iel)), rcoord(0,npol,ielfluid(iel))), &
              max(rcoord(npol,0,ielfluid(iel)), rcoord(npol,npol,ielfluid(iel)))) > xdmf_rmin &
          .and. &
          min(min(thetacoord(0,0,ielfluid(iel)), thetacoord(0,npol,ielfluid(iel))), &
              min(thetacoord(npol,0,ielfluid(iel)), thetacoord(npol,npol,ielfluid(iel)))) < xdmf_thetamax &
          .and. &
          max(max(thetacoord(0,0,ielfluid(iel)), thetacoord(0,npol,ielfluid(iel))), &
              max(thetacoord(npol,0,ielfluid(iel)), thetacoord(npol,npol,ielfluid(iel)))) > xdmf_thetamin) &
          then
          ct = ct + 1
          mask_tp_elem(iel) = .true.
      endif
  enddo

  do iel=1, nel_solid
      if (min(min(rcoord(0,0,ielsolid(iel)), rcoord(0,npol,ielsolid(iel))), &
              min(rcoord(npol,0,ielsolid(iel)), rcoord(npol,npol,ielsolid(iel)))) < xdmf_rmax &
          .and. &
          max(max(rcoord(0,0,ielsolid(iel)), rcoord(0,npol,ielsolid(iel))), &
              max(rcoord(npol,0,ielsolid(iel)), rcoord(npol,npol,ielsolid(iel)))) > xdmf_rmin &
          .and. &
          min(min(thetacoord(0,0,ielsolid(iel)), thetacoord(0,npol,ielsolid(iel))), &
              min(thetacoord(npol,0,ielsolid(iel)), thetacoord(npol,npol,ielsolid(iel)))) < xdmf_thetamax &
          .and. &
          max(max(thetacoord(0,0,ielsolid(iel)), thetacoord(0,npol,ielsolid(iel))), &
              max(thetacoord(npol,0,ielsolid(iel)), thetacoord(npol,npol,ielsolid(iel)))) > xdmf_thetamin) &
          then
          ct = ct + 1
          mask_tp_elem(iel + nel_fluid) = .true.
      endif
  enddo

  nelem_plot = ct * (i_n_xdmf - 1) * (j_n_xdmf - 1)

  allocate(check(nglob_fluid + nglob_solid))
  allocate(mapping(nglob_fluid + nglob_solid))
  allocate(mapping_ijel_iplot(i_n_xdmf, j_n_xdmf, nelem))
  allocate(plotting_mask(i_n_xdmf, j_n_xdmf, nelem))

  check = .false.
  plotting_mask = .false.

  ct = 0

  if (lpr) write(*,*) '   construction of mapping for xdmf plotting...'
  if (lpr) write(*,*) '   ...fluid part...'

  do iel=1, nel_fluid
      if (.not. mask_tp_elem(iel)) cycle
      do i=1, i_n_xdmf
          ipol = i_arr_xdmf(i)
          do j=1, j_n_xdmf
              jpol = j_arr_xdmf(j)

              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_fluid(ipt)

              if (.not. check(idest)) then
                  ct = ct + 1
                  check(idest) = .true.
                  mapping(idest) = ct
                  plotting_mask(i,j,iel) = .true.
              endif
              mapping_ijel_iplot(i,j,iel) = mapping(idest)
          enddo
      enddo
  enddo

  if (lpr) write(*,*) '   ...solid part...'

  do iel=1, nel_solid
      if (.not. mask_tp_elem(iel + nel_fluid)) cycle
      do i=1, i_n_xdmf
          ipol = i_arr_xdmf(i)
          do j=1, j_n_xdmf
              jpol = j_arr_xdmf(j)

              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_solid(ipt) + nglob_fluid

              if (.not. check(idest)) then
                  ct = ct + 1
                  check(idest) = .true.
                  mapping(idest) = ct
                  plotting_mask(i,j,iel + nel_fluid) = .true.
              endif
              mapping_ijel_iplot(i,j,iel + nel_fluid) = mapping(idest)
          enddo
      enddo
  enddo

  deallocate(check, mapping)
  npoint_plot = ct

  allocate(points(1:2,1:npoint_plot))

  if (lpr) write(*,*) '   ...collecting coordinates...'

  points = 0.

  do iel=1, nel_fluid

      do i=1, i_n_xdmf - 1
          ipol = i_arr_xdmf(i)
          ipol1 = i_arr_xdmf(i+1)

          do j=1, j_n_xdmf - 1
              jpol = j_arr_xdmf(j)
              jpol1 = j_arr_xdmf(j+1)

              if (plotting_mask(i,j,iel)) then
                  ct = mapping_ijel_iplot(i,j,iel)
                  points(1,ct) = scoord(ipol,jpol,ielfluid(iel))
                  points(2,ct) = zcoord(ipol,jpol,ielfluid(iel))
              endif

              if (plotting_mask(i+1,j,iel)) then
                  ct = mapping_ijel_iplot(i+1,j,iel)
                  points(1,ct) = scoord(ipol1,jpol,ielfluid(iel))
                  points(2,ct) = zcoord(ipol1,jpol,ielfluid(iel))
              endif

              if (plotting_mask(i+1,j+1,iel)) then
                  ct = mapping_ijel_iplot(i+1,j+1,iel)
                  points(1,ct) = scoord(ipol1,jpol1,ielfluid(iel))
                  points(2,ct) = zcoord(ipol1,jpol1,ielfluid(iel))
              endif

              if (plotting_mask(i,j+1,iel)) then
                  ct = mapping_ijel_iplot(i,j+1,iel)
                  points(1,ct) = scoord(ipol,jpol1,ielfluid(iel))
                  points(2,ct) = zcoord(ipol,jpol1,ielfluid(iel))
              endif
          enddo
      enddo
  enddo

  do iel=1, nel_solid

      do i=1, i_n_xdmf - 1
          ipol = i_arr_xdmf(i)
          ipol1 = i_arr_xdmf(i+1)

          do j=1, j_n_xdmf - 1
              jpol = j_arr_xdmf(j)
              jpol1 = j_arr_xdmf(j+1)

              if (plotting_mask(i,j,iel + nel_fluid)) then
                  ct = mapping_ijel_iplot(i,j,iel + nel_fluid)
                  points(1,ct) = scoord(ipol,jpol,ielsolid(iel))
                  points(2,ct) = zcoord(ipol,jpol,ielsolid(iel))
              endif

              if (plotting_mask(i+1,j,iel + nel_fluid)) then
                  ct = mapping_ijel_iplot(i+1,j,iel + nel_fluid)
                  points(1,ct) = scoord(ipol1,jpol,ielsolid(iel))
                  points(2,ct) = zcoord(ipol1,jpol,ielsolid(iel))
              endif

              if (plotting_mask(i+1,j+1,iel + nel_fluid)) then
                  ct = mapping_ijel_iplot(i+1,j+1,iel + nel_fluid)
                  points(1,ct) = scoord(ipol1,jpol1,ielsolid(iel))
                  points(2,ct) = zcoord(ipol1,jpol1,ielsolid(iel))
              endif

              if (plotting_mask(i,j+1,iel + nel_fluid)) then
                  ct = mapping_ijel_iplot(i,j+1,iel + nel_fluid)
                  points(1,ct) = scoord(ipol,jpol1,ielsolid(iel))
                  points(2,ct) = zcoord(ipol,jpol1,ielsolid(iel))
              endif
          enddo
      enddo
  enddo

  if (lpr) write(*,*) '   .... finished construction of mapping for xdmf plotting'

  if (use_netcdf) then
      call nc_make_snapfile
      call nc_dump_snap_points(points)
  else
      fname = datapath(1:lfdata) // '/xdmf_points_' // appmynum // '.dat'

#if defined(_CRAYFTN)
      open(110, file=trim(fname), access='stream', status='replace', &
          form='unformatted')
#else
      open(110, file=trim(fname), access='stream', status='replace', &
          form='unformatted', convert='big_endian')
#endif

      write(110) points
      close(110)
  endif

  deallocate(points)

  allocate(grid(1:4, 1:nelem_plot))

  if (lpr) write(*,*) '   .... constructing grid for xdmf plotting'

  ct = 1

  do iel=1, nel_fluid
      if (.not. mask_tp_elem(iel)) cycle
      do i=1, i_n_xdmf - 1
          do j=1, j_n_xdmf - 1
              grid(1,ct) = mapping_ijel_iplot(i,j,iel) - 1
              grid(2,ct) = mapping_ijel_iplot(i+1,j,iel) - 1
              grid(3,ct) = mapping_ijel_iplot(i+1,j+1,iel) - 1
              grid(4,ct) = mapping_ijel_iplot(i,j+1,iel) - 1
              ct = ct + 1
          enddo
      enddo
  enddo

  do iel=1, nel_solid
      if (.not. mask_tp_elem(iel + nel_fluid)) cycle
      do i=1, i_n_xdmf - 1
          do j=1, j_n_xdmf - 1
              grid(1,ct) = mapping_ijel_iplot(i,j,iel + nel_fluid) - 1
              grid(2,ct) = mapping_ijel_iplot(i+1,j,iel + nel_fluid) - 1
              grid(3,ct) = mapping_ijel_iplot(i+1,j+1,iel + nel_fluid) - 1
              grid(4,ct) = mapping_ijel_iplot(i,j+1,iel + nel_fluid) - 1
              ct = ct + 1
          enddo
      enddo
  enddo

  if (lpr) write(*,*) '   .... writing grid + header of xdmf to file'

  if (use_netcdf) then
      call nc_dump_snap_grid(grid)
  else
      fname = datapath(1:lfdata) // '/xdmf_grid_' // appmynum // '.dat'

#if defined(_CRAYFTN)
      open(110, file=trim(fname), access='stream', status='replace', &
          form='unformatted')
#else
      open(110, file=trim(fname), access='stream', status='replace', &
          form='unformatted', convert='big_endian')
#endif

      write(110) grid
      close(110)
  endif

  fname = datapath(1:lfdata) // '/xdmf_meshonly_' // appmynum // '.xdmf'
  open(110, file=trim(fname))
  if (use_netcdf) then
      write(110, 732) nelem_plot, nelem_plot, 'hdf', 'netcdf_snap_'//appmynum//'.nc:/grid', &
                      npoint_plot, 'hdf', 'netcdf_snap_'//appmynum//'.nc:/points'
  else
      write(110, 732) nelem_plot, nelem_plot, 'binary', 'xdmf_grid_' // appmynum // '.dat', &
                      npoint_plot, 'binary', 'xdmf_points_' // appmynum // '.dat'
  endif
  close(110)

732 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/&
    '<Grid GridType="Uniform">',/&
    '<Time Value="0.000" />',/&
    '<Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '<DataItem Dimensions="',i10,' 4" NumberType="Int" Format="',A,'" Endian="Big">',/&
    '        ',A,/&
    '</DataItem>',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Dimensions="',i10,' 2" NumberType="Float" Format="',A,'" Endian="Big">',/&
    '        ',A/&
    '</DataItem>',/&
    '</Geometry>',/&
    '</Grid>',/&
    '</Grid>',/&
    '</Domain>',/&
    '</Xdmf>')


  fname = datapath(1:lfdata) // '/xdmf_xml_' // appmynum // '.xdmf'
  open(110, file=trim(fname))
  if (use_netcdf) then
     write(110, 733) nelem_plot, 'hdf', 'netcdf_snap_'//appmynum//'.nc:/grid', &
                     npoint_plot, 'hdf', 'netcdf_snap_'//appmynum//'.nc:/points'
  else
     write(110, 733) nelem_plot, 'binary', 'xdmf_grid_' // appmynum // '.dat', &
                     npoint_plot, 'binary', 'xdmf_points_' // appmynum // '.dat'
  endif
  close(110)

733 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="grid" Dimensions="',i10,' 4" NumberType="Int" Format="',A,'" Endian="Big">',/&
    '  ', A,/&
    '</DataItem>',/&
    '<DataItem Name="points" Dimensions="',i10,' 2" NumberType="Float" Format="',A,'" Endian="Big">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/)


end subroutine dump_xdmf_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine finish_xdmf_xml()

  use data_source, only: src_type

  character(len=120) :: fname

  fname = datapath(1:lfdata) // '/xdmf_xml_' // appmynum // '.xdmf'
  !open(110, file=trim(fname), access='append')
  open(110, file=trim(fname), position='append')
  write(110, 736)

736 format(&
    '</Grid>',/&
    '</Domain>',/&
    '</Xdmf>')

  close(110)
  close(13100)
  if (.not. src_type(1) == 'monopole') close(13101)
  close(13102)
  close(13103)
  close(13104)

end subroutine finish_xdmf_xml
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_mesh_memoryvar_vtk()

  use data_matr, only: points_solid

  integer :: iel, ipol,jpol

  allocate(points_solid(0:npol,0:npol,nel_solid,2))

  do iel=1, nel_solid
     do jpol=0, npol
        do ipol=0, npol
           points_solid(ipol,jpol,iel,1) = scoord(ipol,jpol,ielsolid(iel))
           points_solid(ipol,jpol,iel,2) = zcoord(ipol,jpol,ielsolid(iel))
        enddo
     enddo
  enddo

end subroutine prepare_mesh_memoryvar_vtk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine build_kwf_grid()

  use nc_routines, only: set_npoints
  use data_mesh
  use data_io, only: ibeg, iend, jbeg, jend

  integer               :: iel, ipol, jpol, ct, ipt, idest
  integer, allocatable  :: mapping(:)
  logical, allocatable  :: check(:), mask_tp_elem(:)

  allocate(mask_tp_elem(nelem))
  mask_tp_elem = .false.

  ct = 0

  do iel=1, nel_solid
      if (min(min(rcoord(0,0,ielsolid(iel)), rcoord(0,npol,ielsolid(iel))), &
              min(rcoord(npol,0,ielsolid(iel)), rcoord(npol,npol,ielsolid(iel)))) < kwf_rmax &
          .and. &
          max(max(rcoord(0,0,ielsolid(iel)), rcoord(0,npol,ielsolid(iel))), &
              max(rcoord(npol,0,ielsolid(iel)), rcoord(npol,npol,ielsolid(iel)))) > kwf_rmin &
          .and. &
          min(min(thetacoord(0,0,ielsolid(iel)), thetacoord(0,npol,ielsolid(iel))), &
              min(thetacoord(npol,0,ielsolid(iel)), thetacoord(npol,npol,ielsolid(iel)))) < kwf_thetamax &
          .and. &
          max(max(thetacoord(0,0,ielsolid(iel)), thetacoord(0,npol,ielsolid(iel))), &
              max(thetacoord(npol,0,ielsolid(iel)), thetacoord(npol,npol,ielsolid(iel)))) > kwf_thetamin) &
          then
          ct = ct + 1
          mask_tp_elem(iel) = .true.
      endif
  enddo

  do iel=1, nel_fluid
      if (min(min(rcoord(0,0,ielfluid(iel)), rcoord(0,npol,ielfluid(iel))), &
              min(rcoord(npol,0,ielfluid(iel)), rcoord(npol,npol,ielfluid(iel)))) < kwf_rmax &
          .and. &
          max(max(rcoord(0,0,ielfluid(iel)), rcoord(0,npol,ielfluid(iel))), &
              max(rcoord(npol,0,ielfluid(iel)), rcoord(npol,npol,ielfluid(iel)))) > kwf_rmin &
          .and. &
          min(min(thetacoord(0,0,ielfluid(iel)), thetacoord(0,npol,ielfluid(iel))), &
              min(thetacoord(npol,0,ielfluid(iel)), thetacoord(npol,npol,ielfluid(iel)))) < kwf_thetamax &
          .and. &
          max(max(thetacoord(0,0,ielfluid(iel)), thetacoord(0,npol,ielfluid(iel))), &
              max(thetacoord(npol,0,ielfluid(iel)), thetacoord(npol,npol,ielfluid(iel)))) > kwf_thetamin) &
          then
          ct = ct + 1
          mask_tp_elem(iel + nel_solid) = .true.
      endif
  enddo

  nelem_kwf = ct

  allocate(check(nglob_fluid + nglob_solid))
  allocate(mapping(nglob_fluid + nglob_solid))
  allocate(mapping_ijel_ikwf(0:npol, 0:npol, nelem))
  allocate(kwf_mask(0:npol, 0:npol, nelem))

  check = .false.
  kwf_mask = .false.

  ct = 0

  if (lpr) write(*,*) '   construction of mapping for kwf output...'
  if (lpr) write(*,*) '   ...solid part...'

  do iel=1, nel_solid
      if (.not. mask_tp_elem(iel)) cycle
      do jpol=jbeg, jend
          do ipol=ibeg, iend

              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_solid(ipt) + nglob_fluid

              if (.not. check(idest)) then
                  ct = ct + 1
                  check(idest) = .true.
                  mapping(idest) = ct
                  kwf_mask(ipol,jpol,iel) = .true.
              endif
              mapping_ijel_ikwf(ipol,jpol,iel) = mapping(idest)
          enddo

          ! add GLL points on the axis
          if (axis_solid(iel)) then
              ipol = 0

              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_solid(ipt) + nglob_fluid

              if (.not. check(idest)) then
                  ct = ct + 1
                  check(idest) = .true.
                  mapping(idest) = ct
                  kwf_mask(ipol,jpol,iel) = .true.
              endif
              mapping_ijel_ikwf(ipol,jpol,iel) = mapping(idest)
          endif
      enddo
  enddo

  npoint_solid_kwf = ct

  if (lpr) write(*,*) '   ...fluid part...'

  do iel=1, nel_fluid
      if (.not. mask_tp_elem(iel + nel_solid)) cycle
      do jpol=jbeg, jend
          do ipol=ibeg, iend

              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_fluid(ipt)

              if (.not. check(idest)) then
                  ct = ct + 1
                  check(idest) = .true.
                  mapping(idest) = ct
                  kwf_mask(ipol,jpol,iel + nel_solid) = .true.
              endif
              mapping_ijel_ikwf(ipol,jpol,iel + nel_solid) = mapping(idest)
          enddo

          ! add GLL points on the axis
          if (axis_fluid(iel)) then
              ipol = 0

              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              idest = igloc_fluid(ipt)

              if (.not. check(idest)) then
                  ct = ct + 1
                  check(idest) = .true.
                  mapping(idest) = ct
                  kwf_mask(ipol,jpol,iel + nel_solid) = .true.
              endif
              mapping_ijel_ikwf(ipol,jpol,iel + nel_solid) = mapping(idest)
          endif
      enddo
  enddo

  deallocate(check, mapping)
  npoint_kwf = ct
  npoint_fluid_kwf = ct - npoint_solid_kwf

  call set_npoints(npoint_kwf)

  if (lpr) then
     write(*,*) 'local point number:        ', nelem_kwf * (iend - ibeg + 1) * (jend - jbeg + 1)
     write(*,*) 'after removing duplicates: ', npoint_kwf
     write(*,*) 'compression:               ', &
                 real(npoint_kwf) / real(nelem_kwf * (npol + 1)**2)
  endif

  if (trim(dump_type) == 'displ_only') then
     allocate(midpoint_mesh_kwf(1:nelem_kwf))

     if (lpr) write(*,*) '   .... constructing midpoint grid for kwf output'

     ct = 1

     do iel=1, nel_solid
         if (.not. mask_tp_elem(iel)) cycle
         midpoint_mesh_kwf(ct) = mapping_ijel_ikwf(npol/2,npol/2,iel) - 1
         ct = ct + 1
     enddo

     do iel=1, nel_fluid
         if (.not. mask_tp_elem(iel + nel_solid)) cycle
         midpoint_mesh_kwf(ct) = mapping_ijel_ikwf(npol/2,npol/2,iel + nel_solid) - 1
         ct = ct + 1
     enddo

     allocate(eltype_kwf(1:nelem_kwf))

     eltype_kwf(:) = -2

     ct = 1

     do iel=1, nel_solid
         if (.not. mask_tp_elem(iel)) cycle
         if (eltype(ielsolid(iel)) == 'curved') then
            eltype_kwf(ct) = 0
         else if (eltype(ielsolid(iel)) == 'linear') then
            eltype_kwf(ct) = 1
         else if (eltype(ielsolid(iel)) == 'semino') then
            eltype_kwf(ct) = 2
         else if (eltype(ielsolid(iel)) == 'semiso') then
            eltype_kwf(ct) = 3
         else
            eltype_kwf(ct) = -1
         endif
         ct = ct + 1
     enddo

     do iel=1, nel_fluid
         if (.not. mask_tp_elem(iel + nel_solid)) cycle
         if (eltype(ielfluid(iel)) == 'curved') then
            eltype_kwf(ct) = 0
         else if (eltype(ielfluid(iel)) == 'linear') then
            eltype_kwf(ct) = 1
         else if (eltype(ielfluid(iel)) == 'semino') then
            eltype_kwf(ct) = 2
         else if (eltype(ielfluid(iel)) == 'semiso') then
            eltype_kwf(ct) = 3
         else
            eltype_kwf(ct) = -1
         endif
         ct = ct + 1
     enddo

     allocate(axis_kwf(1:nelem_kwf))

     axis_kwf(:) = -1

     ct = 1

     do iel=1, nel_solid
         if (.not. mask_tp_elem(iel)) cycle
         if (axis_solid(iel)) then
            axis_kwf(ct) = 1
         else
            axis_kwf(ct) = 0
         endif
         ct = ct + 1
     enddo

     do iel=1, nel_fluid
         if (.not. mask_tp_elem(iel + nel_solid)) cycle
         if (axis_fluid(iel)) then
            axis_kwf(ct) = 1
         else
            axis_kwf(ct) = 0
         endif
         ct = ct + 1
     enddo

     allocate(fem_mesh_kwf(1:4, 1:nelem_kwf))

     if (lpr) write(*,*) '   .... constructing finite element grid for kwf output'

     ct = 1

     do iel=1, nel_solid
         if (.not. mask_tp_elem(iel)) cycle
         fem_mesh_kwf(1,ct) = mapping_ijel_ikwf(   0,   0,iel) - 1
         fem_mesh_kwf(2,ct) = mapping_ijel_ikwf(npol,   0,iel) - 1
         fem_mesh_kwf(3,ct) = mapping_ijel_ikwf(npol,npol,iel) - 1
         fem_mesh_kwf(4,ct) = mapping_ijel_ikwf(   0,npol,iel) - 1
         ct = ct + 1
     enddo

     do iel=1, nel_fluid
         if (.not. mask_tp_elem(iel + nel_solid)) cycle
         fem_mesh_kwf(1,ct) = mapping_ijel_ikwf(   0,   0,iel + nel_solid) - 1
         fem_mesh_kwf(2,ct) = mapping_ijel_ikwf(npol,   0,iel + nel_solid) - 1
         fem_mesh_kwf(3,ct) = mapping_ijel_ikwf(npol,npol,iel + nel_solid) - 1
         fem_mesh_kwf(4,ct) = mapping_ijel_ikwf(   0,npol,iel + nel_solid) - 1
         ct = ct + 1
     enddo
  endif

  allocate(sem_mesh_kwf(ibeg:iend, jbeg:jend, 1:nelem_kwf))

  if (lpr) write(*,*) '   .... constructing spectral element grid for kwf output'

  ct = 1

  do iel=1, nel_solid
      if (.not. mask_tp_elem(iel)) cycle
      do ipol=ibeg, iend
         do jpol=jbeg, jend
            sem_mesh_kwf(ipol,jpol,ct) = mapping_ijel_ikwf(ipol,jpol,iel) - 1
         enddo
      enddo
      ct = ct + 1
  enddo

  do iel=1, nel_fluid
      if (.not. mask_tp_elem(iel + nel_solid)) cycle
      do ipol=ibeg, iend
         do jpol=jbeg, jend
            sem_mesh_kwf(ipol,jpol,ct) = mapping_ijel_ikwf(ipol,jpol,iel + nel_solid) - 1
         enddo
      enddo
      ct = ct + 1
  enddo

  if (lpr) write(*,*) '   .... finished construction of mapping for kwf output'

end subroutine build_kwf_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_kwf_grid()

  use nc_routines, only: nc_dump_mesh_kwf, nc_dump_mesh_mp_kwf
  use data_mesh

  integer               :: iel, ipol, jpol, ct
  real(sp), allocatable :: points(:,:)
  real(sp), allocatable :: points_mp(:,:)

  if (lpr) write(*,*) '   ...collecting coordinates...'

  allocate(points(1:npoint_kwf, 2))

  points = 0.

  do iel=1, nel_solid
      do ipol=0, npol
          do jpol=0, npol
              if (kwf_mask(ipol,jpol,iel)) then
                  ct = mapping_ijel_ikwf(ipol,jpol,iel)
                  points(ct,1) = scoord(ipol,jpol,ielsolid(iel))
                  points(ct,2) = zcoord(ipol,jpol,ielsolid(iel))
              endif
          enddo
      enddo
  enddo

  do iel=1, nel_fluid
      do ipol=0, npol
          do jpol=0, npol
              if (kwf_mask(ipol,jpol,iel + nel_solid)) then
                  ct = mapping_ijel_ikwf(ipol,jpol,iel + nel_solid)
                  points(ct,1) = scoord(ipol,jpol,ielfluid(iel))
                  points(ct,2) = zcoord(ipol,jpol,ielfluid(iel))
              endif
          enddo
      enddo
  enddo

  allocate(points_mp(1:nelem_kwf, 2))

  points_mp = 0.

  ct = 1
  do iel=1, nel_solid
     if (kwf_mask(npol/2,jpol/2,iel)) then
        points_mp(ct,1) = scoord(npol/2,npol/2,ielsolid(iel))
        points_mp(ct,2) = zcoord(npol/2,npol/2,ielsolid(iel))
        ct = ct + 1
     endif
  enddo

  do iel=1, nel_fluid
     if (kwf_mask(npol/2,jpol/2,iel + nel_solid)) then
        points_mp(ct,1) = scoord(npol/2,npol/2,ielfluid(iel))
        points_mp(ct,2) = zcoord(npol/2,npol/2,ielfluid(iel))
        ct = ct + 1
     endif
  enddo

  if (use_netcdf) then
      call nc_dump_mesh_kwf(points, npoint_solid_kwf, npoint_fluid_kwf)
      call nc_dump_mesh_mp_kwf(points_mp, nelem_kwf)
  else
     write(*,*) 'ERROR: binary output for non-duplicate mesh not implemented'
     call abort()
  endif

end subroutine dump_kwf_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_kwf_midpoint_xdmf(filename, npoints, nelem)
  character(len=*), intent(in)      :: filename
  integer, intent(in)               :: npoints, nelem

  integer                           :: iinput_xdmf
  character(len=512)                :: filename_np


  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'_mp.xdmf')
  write(iinput_xdmf, 733) npoints, npoints, trim(filename_np), npoints, trim(filename_np)

  write(iinput_xdmf, 734) nelem, nelem, trim(filename_np), "'", "'"

  close(iinput_xdmf)

733 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="points" ItemType="function" function="join($0, $1)" Dimensions="', i10, ' 2">',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_S',/&
    '</DataItem>',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_Z',/&
    '</DataItem>',/&
    '</DataItem>',/,/)

734 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Topology TopologyType="Polyvertex" NumberOfElements="',i10,'">',/&
    '<DataItem ItemType="Uniform" Name="points" DataType="Int" Dimensions="', i10, '" Format="HDF">',/&
    '        ', A, ':/Mesh/midpoint_mesh',/&
    '</DataItem>',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')


  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'_mp_vect.xdmf')
  write(iinput_xdmf, 743) nelem, nelem, trim(filename_np), nelem, trim(filename_np)

  write(iinput_xdmf, 744) nelem, "'", "'"

  close(iinput_xdmf)

743 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="points" ItemType="function" function="join($0, $1)" Dimensions="', i10, ' 2">',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mp_mesh_S',/&
    '</DataItem>',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mp_mesh_Z',/&
    '</DataItem>',/&
    '</DataItem>',/,/)

744 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Topology TopologyType="Polyvertex" NumberOfElements="',i10,'">',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_kwf_fem_xdmf(filename, npoints, nelem)
  character(len=*), intent(in)      :: filename
  integer, intent(in)               :: npoints, nelem

  integer                           :: iinput_xdmf
  character(len=512)                :: filename_np


  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'_fem.xdmf')
  write(iinput_xdmf, 733) npoints, npoints, trim(filename_np), npoints, trim(filename_np)

  write(iinput_xdmf, 734) nelem, nelem, trim(filename_np), "'", "'", &
                          nelem, trim(filename_np), nelem, trim(filename_np)


  close(iinput_xdmf)

733 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="points" ItemType="function" function="join($0, $1)" Dimensions="', i10, ' 2">',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_S',/&
    '</DataItem>',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_Z',/&
    '</DataItem>',/&
    '</DataItem>',/,/)

734 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Topology TopologyType="Quadrilateral" NumberOfElements="', i10, '">',/&
    '<DataItem ItemType="Uniform" Name="points" DataType="Int" Dimensions="', i10, ' 4" Format="HDF">',/&
    '        ', A, ':/Mesh/fem_mesh',/&
    '</DataItem>',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '<Attribute Name="eltype" AttributeType="Scalar" Center="Cell">',/&
    '<DataItem ItemType="Uniform" Name="points" DataType="Int" Dimensions="', i10, '" Format="HDF">',/&
    '        ', A, ':/Mesh/eltype',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="axis" AttributeType="Scalar" Center="Cell">',/&
    '<DataItem ItemType="Uniform" Name="points" DataType="Int" Dimensions="', i10, '" Format="HDF">',/&
    '        ', A, ':/Mesh/axis',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_kwf_sem_xdmf(filename, npoints, nelem)
  character(len=*), intent(in)      :: filename
  integer, intent(in)               :: npoints, nelem

  integer                           :: iinput_xdmf
  character(len=512)                :: filename_np


  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'_sem.xdmf')
  write(iinput_xdmf, 733) npoints, npoints, trim(filename_np), npoints, trim(filename_np)

  write(iinput_xdmf, 734) nelem, (npol+1)**2, nelem, npol+1, npol+1, &
                          trim(filename_np), "'", "'"

  close(iinput_xdmf)

733 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="points" ItemType="function" function="join($0, $1)" Dimensions="', i10, ' 2">',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_S',/&
    '</DataItem>',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_Z',/&
    '</DataItem>',/&
    '</DataItem>',/,/)

734 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Topology TopologyType="Polygon" NumberOfElements="',i10,'" NodesPerElement="', i5, '">',/&
    '<DataItem ItemType="Uniform" Name="points" DataType="Int" Dimensions="', i10, i3, i3, '" Format="HDF">',/&
    '        ', A, ':/Mesh/sem_mesh',/&
    '</DataItem>',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_kwf_gll_xdmf(filename, npoints)
  character(len=*), intent(in)      :: filename
  integer, intent(in)               :: npoints

  integer                           :: iinput_xdmf
  character(len=512)                :: filename_np


  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'_gll.xdmf')
  write(iinput_xdmf, 733) npoints, npoints, trim(filename_np), npoints, trim(filename_np)

  write(iinput_xdmf, 734) npoints, "'", "'"

  close(iinput_xdmf)

733 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="points" ItemType="function" function="join($0, $1)" Dimensions="', i10, ' 2">',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_S',/&
    '</DataItem>',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_Z',/&
    '</DataItem>',/&
    '</DataItem>',/,/)

734 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Topology TopologyType="Polyvertex" NumberOfElements="',i10'">',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Dumps the mesh (s,z) [m] and related constant fields in binary format as
!! needed to compute waveform kernels from the strain and velocity fields.
!! The distinction between different dumping methods is honored here,
!! and influences the amount of additional dumpsters (prefactors, derivatives).
!! In a nutshell, the end-member dumping methods constitute
!! 1) computing strain and velocity on-the-fly, i.e. only dumping the mesh here;
!! 2) only dumping the already-known fields (displacement, potential) on-the-fly
!!    and dump a number of constant fields here.
!! The latter choice is more memory- and CPU-efficient, but requires
!! significant post-processing AND dumping the entire SEM mesh.
!! See compute_strain in time_evol_wave.f90 for more info.
subroutine dump_wavefields_mesh_1d

  use data_mesh
  use data_io, only: ibeg, iend, jbeg, jend, ndumppts_el
  use data_spec, only: G1T, G2T, G2
  use data_pointwise
  use nc_routines, only: nc_dump_mesh_sol, nc_dump_mesh_flu

  real(kind=dp), dimension(:,:,:), allocatable :: ssol, zsol
  real(kind=dp), dimension(:,:,:), allocatable :: sflu, zflu

  integer :: iel, ipol, jpol

  ! Dump subset of GLL point grid for 'fullfields'
  if (dump_type == 'fullfields') then

     if (lpr) then
        write(*,*)'  set strain dumping GLL boundaries to:'
        write(*,*)'    ipol=', ibeg, iend
        write(*,*)'    jpol=', jbeg, jend
     endif

     allocate(ssol(ibeg:iend,jbeg:jend,nel_solid))
     allocate(zsol(ibeg:iend,jbeg:jend,nel_solid))
     allocate(sflu(ibeg:iend,jbeg:jend,nel_fluid))
     allocate(zflu(ibeg:iend,jbeg:jend,nel_fluid))

  endif

  !!! SB
  if (dump_type == 'coupling' .or. dump_type == 'coupling_box') then

     ibeg = 0
     iend = 4
     jbeg = 0
     jend = 4
     if (lpr) then
        write(*,*)'  Coupling : FORCE GLL boundaries to:'
        write(*,*)'    ipol=', ibeg, iend
        write(*,*)'    jpol=', jbeg, jend
     endif

     allocate(ssol(ibeg:iend,jbeg:jend,nel_solid))
     allocate(zsol(ibeg:iend,jbeg:jend,nel_solid))
     allocate(sflu(ibeg:iend,jbeg:jend,nel_fluid))
     allocate(zflu(ibeg:iend,jbeg:jend,nel_fluid))

!     have_fluid = .false. !!!! SB for debug...

  endif
  !!! END SB

  if (dump_type == 'displ_only' .or. dump_type == 'strain_only') then
     call dump_kwf_grid()
  else
     ! compute solid grid
     do iel=1,nel_solid
         do jpol=jbeg,jend
             do ipol=ibeg,iend
                 ssol(ipol,jpol,iel) = scoord(ipol,jpol,ielsolid(iel))
                 zsol(ipol,jpol,iel) = zcoord(ipol,jpol,ielsolid(iel))
             enddo
        enddo
     enddo

     if (lpr) write(*,*)'  dumping solid submesh for kernel wavefields...'
     if (use_netcdf) then
         call nc_dump_mesh_sol(real(ssol(ibeg:iend,jbeg:jend,:)), &
                               real(zsol(ibeg:iend,jbeg:jend,:)))

     else
         open(unit=2500+mynum,file=datapath(1:lfdata)//'/strain_mesh_sol_'&
                                   //appmynum//'.dat', &
                                   FORM="UNFORMATTED",STATUS="REPLACE")

         write(2500+mynum)ssol(ibeg:iend,jbeg:jend,:),zsol(ibeg:iend,jbeg:jend,:)
         close(2500+mynum)
     endif
     deallocate(ssol,zsol)

     ! compute fluid grid
     if (have_fluid) then
         do iel=1,nel_fluid
             do jpol=jbeg,jend
                 do ipol=ibeg,iend
                     sflu(ipol,jpol,iel) = scoord(ipol,jpol,ielfluid(iel))
                     zflu(ipol,jpol,iel) = zcoord(ipol,jpol,ielfluid(iel))
                 enddo
              enddo
         enddo
         if (lpr) write(*,*)'  dumping fluid submesh for kernel wavefields...'
         if (use_netcdf) then
             call nc_dump_mesh_flu(real(sflu(ibeg:iend,jbeg:jend,:)), &
                                   real(zflu(ibeg:iend,jbeg:jend,:)))
         else
             open(unit=2600+mynum,file=datapath(1:lfdata)//'/strain_mesh_flu_'&
                  //appmynum//'.dat', &
                  FORM="UNFORMATTED",STATUS="REPLACE")
             write(2600+mynum)sflu(ibeg:iend,jbeg:jend,:),zflu(ibeg:iend,jbeg:jend,:)
             close(2600+mynum)
         endif
         deallocate(sflu,zflu)
     endif ! have_fluid
  endif


  ! In the following: Only dumping additional arrays if displacements only
  ! are dumped as wavefields to reconstruct the strains.

  select case (dump_type)
  case ('displ_only')
     if (lpr) then
        write(*,*)'  strain dump: only elementwise displacement'
     endif
  case ('strain_only')
     if (lpr) then
        write(*,*)'  strain dump: only pointwise strain '
     endif
  case ('displ_velo')
     if (lpr) then
        write(*,*)'  strain dump: only displacement/velocity, potentials'
        write(*,*)'  ...now dumping global pointwise deriv. terms, etc....'
     endif

     ! Dump pointwise derivative matrices in solid
     open(unit=2600+mynum,file=datapath(1:lfdata)//'/pointwise_deriv_sol_'&
                               //appmynum//'.dat', &
                               FORM="UNFORMATTED",STATUS="REPLACE")
     write(2600+mynum)DzDeta_over_J_sol,DzDxi_over_J_sol, &
                      DsDeta_over_J_sol,DsDxi_over_J_sol
     close(2600+mynum)

     if (have_fluid) then
     ! Dump pointwise derivative matrices in fluid
        open(unit=2600+mynum,file=datapath(1:lfdata)//'/pointwise_deriv_flu_'&
             //appmynum//'.dat', &
             FORM="UNFORMATTED",STATUS="REPLACE")
        write(2600+mynum)DzDeta_over_J_flu,DzDxi_over_J_flu, &
             DsDeta_over_J_flu,DsDxi_over_J_flu
        close(2600+mynum)

        ! Dump inverse density inside fluid
        open(unit=2600+mynum,file=datapath(1:lfdata)//'/inv_rho_fluid_'&
             //appmynum//'.dat', &
             FORM="UNFORMATTED",STATUS="REPLACE")
        write(2600+mynum) inv_rho_fluid
        close(2600+mynum)
     endif

     ! Dump Lagrange interpolant derivatives
     open(unit=2600+mynum,file=datapath(1:lfdata)//'/lagrange_derivs_'&
                               //appmynum//'.dat', &
                               FORM="UNFORMATTED",STATUS="REPLACE")
     write(2600+mynum) G1T, G2T, G2
     close(2600+mynum)

     write(*,*)'  ...dumped it all.'

  case ('fullfields') !Hardcoded choice in parameters.f90:110
     if (lpr) then
        write(*,*)'  strain dump: Global strain tensor and velocity fields'
        write(*,*)'  ....no need to dump anything else.'
     endif

  !!! SB
  case ('coupling')
     if (lpr) then
        write(*,*)'  strain dump: Global strain tensor and velocity fields'
        write(*,*)'  ....no need to dump anything else.'
     endif

  case ('coupling_box')
     if (lpr) then
        write(*,*)'  strain dump: Global strain tensor and velocity fields'
        write(*,*)'  ....no need to dump anything else.'
     endif
  !!! END SB
  case default
     if (lpr) then
        write(*,*)'  wavefield dumping type',dump_type,' unknown!'
        write(*,*)'  select from 1) displ_only, 2) displ_velo, 2) fullfields'
     endif
     stop
  end select

end subroutine dump_wavefields_mesh_1d
!-----------------------------------------------------------------------------------------

end module meshes_io
!=========================================================================================
