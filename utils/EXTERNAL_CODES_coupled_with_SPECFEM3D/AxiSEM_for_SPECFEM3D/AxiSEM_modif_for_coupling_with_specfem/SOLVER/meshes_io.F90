!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
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
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!> This module contains routines that compute and dump the respective meshes
!! underlying the actual wavefields to be dumped in the time loop
!! which is done in wavefields_io. This module needs pre-loop variables such
!! as mesh coordinates and is therefore cut off from the dumping module.
module meshes_io

  use global_parameters
  use data_mesh
  use data_proc
  use data_io
  
  use utlity, ONLY : scoord, zcoord, rcoord, thetacoord
  
  implicit none
  
  private
  public :: fldout_cyl2
  public :: finish_xdmf_xml
  public :: dump_wavefields_mesh_1d
  public :: dump_glob_grid_midpoint
  public :: dump_xdmf_grid
  public :: dump_solid_grid
  public :: dump_fluid_grid
  public :: prepare_mesh_memoryvar_vtk

contains

!-----------------------------------------------------------------------------
!> Dumps the mesh (s,z) [m] in ASCII format as needed to visualize global 
!! snapshots, and additionally the constant factors preceding the displacement
!! in the fluid, namely rho^{-1} and (rho s)^{-1}.
!! When reading the fluid wavefield, one therefore needs to multiply all 
!! components with inv_rho_fluid and the phi component with one/scoord!
!! Convention for order in the file: First the fluid, then the solid domain.
subroutine dump_glob_grid_midpoint(ibeg,iend,jbeg,jend)
  
  use data_pointwise, ONLY : inv_rho_fluid
  use data_mesh,      only : npol, nel_fluid, nel_solid
  
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

           if ( axis_fluid(iel) .and. ipol==0 ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
subroutine dump_xdmf_grid()

  use nc_routines,      only: nc_dump_snap_points, nc_dump_snap_grid, nc_make_snapfile
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

  if (lpr) write(6,*) '   construction of mapping for xdmf plotting...'
  if (lpr) write(6,*) '   ...fluid part...'

  do iel=1, nel_fluid
      if (.not.  mask_tp_elem(iel)) cycle
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
  
  if (lpr) write(6,*) '   ...solid part...'
  
  do iel=1, nel_solid
      if (.not.  mask_tp_elem(iel + nel_fluid)) cycle
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
  
  if (lpr) write(6,*) '   ...collecting coordinates...'

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
  
  if (lpr) write(6,*) '   .... finished construction of mapping for xdmf plotting'

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
  end if

  deallocate(points)

  allocate(grid(1:4, 1:nelem_plot))
  
  if (lpr) write(6,*) '   .... constructing grid for xdmf plotting'
  
  ct = 1
  
  do iel=1, nel_fluid
      if (.not.  mask_tp_elem(iel)) cycle
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
      if (.not.  mask_tp_elem(iel + nel_fluid)) cycle
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
  
  if (lpr) write(6,*) '   .... writing grid + header of xdmf to file'
  
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
  end if
  
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
    '  <Grid GridType="Uniform">',/&
    '    <Time Value="0.000" />',/&
    '    <Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '      <DataItem Dimensions="',i10,' 4" NumberType="Int" Format="',A,'" Endian="Big">',/&
    '        ',A,/&
    '      </DataItem>',/&
    '    </Topology>',/&
    '    <Geometry GeometryType="XY">',/&
    '      <DataItem Dimensions="',i10,' 2" NumberType="Float" Format="',A,'" Endian="Big">',/&
    '        ',A/&
    '      </DataItem>',/&
    '    </Geometry>',/&
    '  </Grid>',/&
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
!=============================================================================

!-----------------------------------------------------------------------------
subroutine finish_xdmf_xml()

  use data_source, ONLY : src_type

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
  if (.not. src_type(1)=='monopole') close(13101)
  close(13102)
  close(13103)
  close(13104)

end subroutine finish_xdmf_xml
!=============================================================================

!-----------------------------------------------------------------------------
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
!=============================================================================

!-----------------------------------------------------------------------------
!> Dumps the mesh (s,z) [m] in ASCII format as needed to visualize snapshots 
!! in the solid region only.
!! Convention for order in the file: First the fluid, then the solid domain.
subroutine dump_solid_grid(ibeg,iend,jbeg,jend)

  
  integer, intent(in) :: ibeg,iend,jbeg,jend 
  integer             :: iel, ipol,jpol

  open(unit=2500+mynum,file=datapath(1:lfdata)//'/solid_grid_'&
                            //appmynum//'.dat')
  do iel=1,nel_solid
     do jpol=jbeg,jend
        do ipol=ibeg,iend
           write(2500+mynum,*)scoord(ipol,jpol,ielsolid(iel)), &
                              zcoord(ipol,jpol,ielsolid(iel))
        enddo
     enddo
  enddo
  close(2500+mynum)

end subroutine dump_solid_grid
!=============================================================================

!-----------------------------------------------------------------------------
!> Dumps the mesh (s,z) [m] in ASCII format as needed to visualize snapshots 
!! in the fluid region only, and additionally the constant factors preceding 
!! the displacement in the fluid, namely rho^{-1} and (rho s)^{-1}.
!! When reading the fluid wavefield, one therefore needs to multiply all 
!! components with inv_rho_fluid and the phi component with one/scoord!
!! Convention for order in the file: First the fluid, then the solid domain.
subroutine dump_fluid_grid(ibeg,iend,jbeg,jend)

  use data_pointwise, ONLY : inv_rho_fluid
  
  
  integer, intent(in) :: ibeg,iend,jbeg,jend
  integer             :: iel, ipol,jpol
  
  ! When reading the fluid wavefield, one needs to multiply all components 
  ! with inv_rho_fluid and the phi component with one/scoord!!

  open(unit=2500+mynum,file=datapath(1:lfdata)//&
                            '/fluid_grid_'//appmynum//'.dat')
  open(unit=2600+mynum,file=datapath(1:lfdata)//&
                            '/inv_rho_scoord_fluid_flusnaps_'&
                            //appmynum//'.dat', STATUS="REPLACE")
  do iel=1,nel_fluid
     do jpol=jbeg,jend
        do ipol=ibeg,iend
           write(2500+mynum,*)scoord(ipol,jpol,ielfluid(iel)), &
                              zcoord(ipol,jpol,ielfluid(iel))
           if ( axis_fluid(iel) .and. ipol==0 ) then
              ! Axis s=0! write 1 instead of 1/s and then multiply 
              ! with the correct factor dsdchi, obtained by L'Hospital's rule 
              ! (see routine fluid_snapshot below).
              write(2600+mynum,*)inv_rho_fluid(ipol,jpol,iel),one
           else  
              write(2600+mynum,*)inv_rho_fluid(ipol,jpol,iel), &
                                 one/scoord(ipol,jpol,ielfluid(iel))
           endif
        enddo
     enddo
  enddo
  close(2500+mynum)
  close(2600+mynum)

end subroutine dump_fluid_grid
!=============================================================================

!-----------------------------------------------------------------------------
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
!! 
!! CURRENTLY HARDCODED TO dump_type=='fullfields'
subroutine dump_wavefields_mesh_1d

  use data_mesh
  use data_io, ONLY : ibeg,iend,ndumppts_el
  use data_spec, ONLY : G1T,G2T,G2
  use data_pointwise
  use nc_routines, ONLY: nc_dump_mesh_sol, nc_dump_mesh_flu
  
  real(kind=dp)   , dimension(:,:,:), allocatable :: ssol, zsol
  real(kind=dp)   , dimension(:,:,:), allocatable :: sflu, zflu
  
  integer :: iel, ipol, jpol
  
  ! Dump entire (including duplicate) GLL point grid for displ_only
  if (dump_type=='displ_only') then
     
  else ! Free choice for other dumping method

     if (lpr) then
        write(6,*)'  set strain dumping GLL boundaries to:'
        write(6,*)'    ipol=',ibeg,iend
     endif

     allocate(ssol(ibeg:iend,ibeg:iend,nel_solid))
     allocate(zsol(ibeg:iend,ibeg:iend,nel_solid))
     allocate(sflu(ibeg:iend,ibeg:iend,nel_fluid))
     allocate(zflu(ibeg:iend,ibeg:iend,nel_fluid))

  endif


  ! compute solid grid
  do iel=1,nel_solid
      do jpol=ibeg,iend
          do ipol=ibeg,iend
              ssol(ipol,jpol,iel) = scoord(ipol,jpol,ielsolid(iel))
              zsol(ipol,jpol,iel) = zcoord(ipol,jpol,ielsolid(iel))
          enddo
     enddo
  enddo


  if (lpr) write(6,*)'  dumping solid submesh for kernel wavefields...'
  if (use_netcdf) then
      call nc_dump_mesh_sol(real(ssol(ibeg:iend,ibeg:iend,:)), &
                            real(zsol(ibeg:iend,ibeg:iend,:)))

  else
      open(unit=2500+mynum,file=datapath(1:lfdata)//'/strain_mesh_sol_'&
                                //appmynum//'.dat', &
                                FORM="UNFORMATTED",STATUS="REPLACE")

      write(2500+mynum)ssol(ibeg:iend,ibeg:iend,:),zsol(ibeg:iend,ibeg:iend,:)
      close(2500+mynum)
  end if
  deallocate(ssol,zsol)

  ! compute fluid grid
  if (have_fluid) then
      do iel=1,nel_fluid
          do jpol=ibeg,iend
              do ipol=ibeg,iend
                  sflu(ipol,jpol,iel) = scoord(ipol,jpol,ielfluid(iel))
                  zflu(ipol,jpol,iel) = zcoord(ipol,jpol,ielfluid(iel))
              enddo
           enddo
      enddo
      if (lpr) write(6,*)'  dumping fluid submesh for kernel wavefields...'
      if (use_netcdf) then
          call nc_dump_mesh_flu(real(sflu(ibeg:iend,ibeg:iend,:)),&
                                real(zflu(ibeg:iend,ibeg:iend,:)))
      else
          open(unit=2600+mynum,file=datapath(1:lfdata)//'/strain_mesh_flu_'&
               //appmynum//'.dat', &
               FORM="UNFORMATTED",STATUS="REPLACE")
          write(2600+mynum)sflu(ibeg:iend,ibeg:iend,:),zflu(ibeg:iend,ibeg:iend,:)
          close(2600+mynum)
      end if
      deallocate(sflu,zflu)
  endif ! have_fluid


  ! In the following: Only dumping additional arrays if displacements only 
  ! are dumped as wavefields to reconstruct the strains.

  select case (dump_type)
  case ('displ_only')
     if (lpr) then     
        write(6,*)'  strain dump: only displacement/velocity, potentials'
        write(6,*)'  ...now dumping global pointwise deriv. terms, etc....'
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

     write(6,*)'  ...dumped it all.'

 case ('fullfields') !Hardcoded choice in parameters.f90:110
     if (lpr) then
        write(6,*)'  strain dump: Global strain tensor and velocity fields'
        write(6,*)'  ....no need to dump anything else.'
     endif
  case default
     if (lpr) then 
        write(6,*)'  wavefield dumping type',dump_type,' unknown!'
        write(6,*)'  select from 1) displ_only, 2) fullfields'
     endif
     stop
  end select

end subroutine dump_wavefields_mesh_1d
!=============================================================================

!-----------------------------------------------------------------------------
!> Dumps the mesh (s,z) [m] along with corresponding field f in ASCII format.
!! At this point only used in computing the valence. 
!! Note that for higher-frequency meshes these files can be very large.
!! flag_norm is to be set to one if the output is to be normalized.
subroutine fldout_cyl2(fname,nel,f,ibeg,iend,jbeg,jend,flag_norm,domain)

  use utlity, ONLY : compute_coordinates
  
  
  character(len=80), intent(in)   :: fname
  character(len=5), intent(in)    :: domain
  integer, intent(in)             :: flag_norm,ibeg,iend,jbeg,jend,nel
  real(kind=realkind), intent(in) :: f(ibeg:,jbeg:,:)
  integer                         :: lf,ielem, ipol,jpol,iel
  real(kind=dp)                   :: r, theta, s, z
  real(kind=realkind)             :: fnr,afnr
  
  lf=index(fname,' ')-1
  open(unit=10000+mynum,file=infopath(1:lfinfo)//'/'//fname(1:lf)//'_'&
                             //appmynum//'.dat')

  fnr = 1.
  if (flag_norm == 1) then
     fnr = zero
     do ielem = 1, nel
        do jpol = jbeg, jend
           do ipol = ibeg, iend
              fnr = max(fnr,abs(f(ielem,ipol,jpol)))
           end do
        end do
     end do   
  end if
  afnr = one/fnr
  do ielem = 1, nel
     if (domain=='total') iel=ielem
     if (domain=='solid') iel=ielsolid(ielem)
     if (domain=='fluid') iel=ielfluid(ielem)
     do jpol = jbeg, jend !,npol/2
        do ipol = ibeg, iend !,npol/2
           call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
           write(10000+mynum,*) s/router,z/router,afnr*f(ipol,jpol,ielem)
        end do
     end do
  end do
  close(10000+mynum)

end subroutine fldout_cyl2
!=============================================================================

end module meshes_io
