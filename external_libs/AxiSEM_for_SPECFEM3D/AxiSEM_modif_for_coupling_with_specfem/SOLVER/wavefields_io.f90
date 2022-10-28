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
!> Contains all routines that dump entire wavefields during the time loop.
!! Optimization of I/O therefore happens here and nowhere else.
!! The corresponding meshes are dumped in meshes_io.
module wavefields_io

  use global_parameters
  use data_mesh
  use data_proc
  use data_io
  use nc_routines
  use nc_snapshots

  implicit none

  private

  public :: dump_field_1d
  public :: glob_snapshot_xdmf
  public :: glob_snapshot_midpoint
  public :: dump_velo_global
  public :: dump_disp_global
  public :: dump_disp
  public :: dump_velo_dchi

contains

!-----------------------------------------------------------------------------------------
!> Dumps the global displacement snapshots [m] in binary  format
!! When reading the fluid wavefield, one needs to multiply all
!! components with inv_rho_fluid and the phi component with one/scoord
!! as dumped by the corresponding routine dump_glob_grid!
!! Convention for order in the file: First the fluid, then the solid domain.
!! MvD: loop increment npol/2 -> what if npol odd?
subroutine glob_snapshot_midpoint(f_sol, chi, ibeg, iend, jbeg, jend, isnap)

   use data_source, only: src_type
   use pointwise_derivatives, only: axisym_gradient_fluid, dsdf_fluid_axis
   use data_mesh, only: npol, nel_solid, nel_fluid

   integer, intent(in)             :: ibeg, iend, jbeg, jend, isnap
   real(kind=realkind), intent(in) :: f_sol(0:,0:,:,:)
   real(kind=realkind), intent(in) :: chi(0:,0:,:)

   real(kind=realkind)             :: usz_fluid(0:npol,0:npol,1:nel_fluid,2)
   character(len=4)                :: appisnap
   integer                         :: iel, ipol, jpol, iidim
   real(kind=realkind)             :: dsdchi, prefac

   ! When reading the fluid wavefield, one needs to multiply all components
   ! with inv_rho_fluid and the phi component with one/scoord!!

   if (src_type(1) == 'monopole') prefac = zero
   if (src_type(1) == 'dipole')   prefac = one
   if (src_type(1) == 'quadpole') prefac = two

   call define_io_appendix(appisnap, isnap)

   open(unit=2500+mynum, &
        file=datapath(1:lfdata)//'/snap_'//appmynum//'_'//appisnap//'.dat', &
        FORM="UNFORMATTED", STATUS="REPLACE")

   if (have_fluid) then
      call axisym_gradient_fluid(chi, usz_fluid)
      do iel=1, nel_fluid
         do jpol=0, npol, npol/2
            do ipol=0, npol, npol/2

               if (axis_fluid(iel)) then
                  call dsdf_fluid_axis(chi(:,:,iel), iel, jpol, dsdchi)
                  write(2500+mynum) usz_fluid(ipol,jpol,iel,1), &
                                    0, &
                                    usz_fluid(ipol,jpol,iel,2)
               else
                  write(2500+mynum) usz_fluid(ipol,jpol,iel,1), &
                                    0, &
                                    usz_fluid(ipol,jpol,iel,2)
               endif
            enddo
         enddo
      enddo
   endif ! have_fluid

   do iel=1, nel_solid
      do jpol=0, npol, npol/2
         do ipol=0, npol, npol/2
            write(2500+mynum) (f_sol(ipol,jpol,iel,iidim), iidim=1,3)
         enddo
      enddo
   enddo
   close(2500+mynum)

end subroutine glob_snapshot_midpoint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Dumps the global displacement snapshots in binary plus XDMF descriptor
subroutine glob_snapshot_xdmf(f_sol, chi, t, isnap)

   use data_source, only: src_type
   use data_pointwise, only: inv_rho_fluid
   use pointwise_derivatives, only: axisym_gradient_fluid, dsdf_fluid_axis
   use nc_snapshots, only: nc_dump_snapshot
   use data_mesh, only: npol, nel_solid, nel_fluid

   real(kind=realkind), intent(in) :: f_sol(0:,0:,:,:)
   real(kind=realkind), intent(in) :: chi(0:,0:,:)
   real(dp), intent(in)            :: t
   integer, intent(in)             :: isnap

   character(len=4)                :: appisnap
   integer                         :: n_xdmf_fl, n_xdmf_sol
   real(sp), allocatable           :: u(:,:), usz_fl(:,:,:,:), u_fl(:,:,:,:)
   real(sp), allocatable           :: straintrace(:,:,:,:), straintrace_mask(:,:)
   real(sp), allocatable           :: curlinplane(:,:,:,:), curlinplane_mask(:,:)
   real(kind=realkind)             :: f_sol_spz(0:npol,0:npol,1:nel_solid,3)
   character(len=120)              :: fname


   allocate(u(1:3,1:npoint_plot))
   allocate(straintrace_mask(1, npoint_plot))
   allocate(straintrace(0:npol, 0:npol, nel_fluid + nel_solid,1))
   allocate(curlinplane_mask(1, npoint_plot))
   allocate(curlinplane(0:npol, 0:npol, nel_fluid + nel_solid,1))

   ! convert +- to sp in case of dipole
   if (src_type(1) == 'dipole') then
      f_sol_spz(:,:,:,1) = f_sol(:,:,:,1) + f_sol(:,:,:,2)
      f_sol_spz(:,:,:,2) = f_sol(:,:,:,1) - f_sol(:,:,:,2)
      f_sol_spz(:,:,:,3) = f_sol(:,:,:,3)
   else
      f_sol_spz = f_sol
   endif

   call define_io_appendix(appisnap, isnap)


   ! Collect displacements in variable U
   if (have_fluid) then
      allocate(usz_fl(0:npol,0:npol,1:nel_fluid,2))
      allocate(u_fl(0:npol,0:npol,1:nel_fluid,3))
      call axisym_gradient_fluid(chi, usz_fl)

      n_xdmf_fl = count(plotting_mask(:,:,1:nel_fluid))

      u_fl(:,:,:,1) = usz_fl(:,:,:,1) * inv_rho_fluid
      u_fl(:,:,:,2) = 0
      u_fl(:,:,:,3) = usz_fl(:,:,:,2) * inv_rho_fluid

      call xdmf_mapping(u_fl, mapping_ijel_iplot(:,:,1:nel_fluid), plotting_mask(:,:,1:nel_fluid), &
                        i_arr_xdmf, j_arr_xdmf, u(:,1:n_xdmf_fl))

      deallocate(usz_fl, u_fl)
   endif

   n_xdmf_sol = count(plotting_mask(:,:,nel_fluid+1:))
   call xdmf_mapping(f_sol_spz, mapping_ijel_iplot(:,:,nel_fluid+1:), &
                     plotting_mask(:,:,nel_fluid+1:), &
                     i_arr_xdmf, j_arr_xdmf, u(:,:))

   call calc_straintrace(f_sol, chi, straintrace)
   call xdmf_mapping(straintrace, mapping_ijel_iplot(:,:,:), &
                     plotting_mask(:,:,:), &
                     i_arr_xdmf, j_arr_xdmf, straintrace_mask)

   call calc_curlinplane(f_sol, curlinplane)
   call xdmf_mapping(curlinplane, mapping_ijel_iplot(:,:,:), &
                     plotting_mask(:,:,:), &
                     i_arr_xdmf, j_arr_xdmf, curlinplane_mask)
   ! Write variable u to respective file (binary or netcdf)
   if (use_netcdf) then
       call nc_dump_snapshot(u, straintrace_mask, curlinplane_mask, isnap)
   else
       write(13100) u(1,:)
       if (src_type(1) /= 'monopole') write(13101) u(2,:)
       write(13102) u(3,:)
       write(13103) straintrace_mask
       write(13104) curlinplane_mask
   endif

   deallocate(u)


   ! Write header into XDMF (text) file
   fname = datapath(1:lfdata) // '/xdmf_xml_' // appmynum // '.xdmf'
   !open(110, file=trim(fname), access='append')
   open(110, file=trim(fname), position='append')

   if (use_netcdf) then
       if (src_type(1) == 'monopole') then
           write(110, 736) appisnap, t, nelem_plot, "'", "'", "'", "'", &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, appisnap, appisnap, &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc'
       else
           write(110, 737) appisnap, t, nelem_plot, "'", "'", "'", "'", &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, appisnap, appisnap, appisnap, &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc', &
                       npoint_plot, isnap-1, npoint_plot, &
                       nsnap, npoint_plot, &
                       'netcdf_snap_'//appmynum//'.nc'
       endif !monopole

   else
       if (src_type(1) == 'monopole') then
           write(110, 734) appisnap, t, nelem_plot, "'", "'", "'", "'", &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_s_'//appmynum//'.dat', &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_z_'//appmynum//'.dat', &
                       npoint_plot, appisnap, appisnap, &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_trace_'//appmynum//'.dat', &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_curlip_'//appmynum//'.dat'
       else
           write(110, 735) appisnap, t, nelem_plot, "'", "'", "'", "'", &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_s_'//appmynum//'.dat', &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_p_'//appmynum//'.dat', &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_z_'//appmynum//'.dat', &
                       npoint_plot, appisnap, appisnap, appisnap, &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_trace_'//appmynum//'.dat', &
                       npoint_plot, isnap-1, npoint_plot, nsnap, npoint_plot, &
                       'xdmf_snap_curlip_'//appmynum//'.dat'
       endif !monopole
   endif !use_netcdf

734 format(&
    '<Grid Name="', A,'" GridType="Uniform">',/&
    '<Time Value="',F8.2,'" />',/&
    '<Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '<Attribute Name="u_s" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="u_z" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="abs" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="function" function="sqrt($0 * $0 + $1 * $1)" Dimensions="', I10,'">',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_s"]/DataItem[1]',/&
    '</DataItem>',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_z"]/DataItem[1]',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="straintrace" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="curlinplane" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '</Grid>',/)

735 format(&
    '<Grid Name="', A,'" GridType="Uniform">',/&
    '<Time Value="',F8.2,'" />',/&
    '<Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '<Attribute Name="u_s" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="u_p" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="u_z" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="abs" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="function" function="sqrt($0 * $0 + $1 * $1 + $2 * $2)" Dimensions="', I10,'">',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_s"]/DataItem[1]',/&
    '</DataItem>',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_p"]/DataItem[1]',/&
    '</DataItem>',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_z"]/DataItem[1]',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="straintrace" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="curlinplane" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary" Endian="Big">',/&
    '                   ', A,/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '</Grid>',/)

736 format(&
    '<Grid Name="', A,'" GridType="Uniform">',/&
    '<Time Value="',F8.2,'" />',/&
    '<Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '<Attribute Name="u_s" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 3" Format="XML">',/&
    '                      ', i10,  '         0          0',/&
    '                               1         1          1',/&
    '                               1',   i10,'          1',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, ' 2" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/displacement',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="u_z" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 3" Format="XML">',/&
    '                      ', i10,  '         0          1',/&
    '                               1         1          1',/&
    '                               1',   i10,'          1',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, ' 2" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/displacement',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="abs" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="function" function="sqrt($0 * $0 + $1 * $1)" Dimensions="', I10,'">',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_s"]/DataItem[1]',/&
    '</DataItem>',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_z"]/DataItem[1]',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="divergence(u)" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                      ', i10,  '         0  ',/&
    '                               1         1  ',/&
    '                               1',   i10,'  ',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/straintrace',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="curl(u)" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                      ', i10,  '         0  ',/&
    '                               1         1  ',/&
    '                               1',   i10,'  ',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/curlinplane',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '</Grid>',/)

737 format(&
    '<Grid Name="', A,'" GridType="Uniform">',/&
    '<Time Value="',F8.2,'" />',/&
    '<Topology TopologyType="Quadrilateral" NumberOfElements="',i10,'">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>',/&
    '<Attribute Name="u_s" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 3" Format="XML">',/&
    '                      ', i10,'         0          0',/&
    '                               1         1          1',/&
    '                               1',  i10,'          1',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, ' 3" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/displacement',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="u_p" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 3" Format="XML">',/&
    '                      ', i10,'         0          1',/&
    '                               1         1          1',/&
    '                               1',  i10,'          1',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, ' 3" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/displacement',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="u_z" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 3" Format="XML">',/&
    '                      ', i10,'         0          2',/&
    '                               1         1          1',/&
    '                               1',   i10,'          1',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, ' 3" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/displacement',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="abs" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="function" function="sqrt($0 * $0 + $1 * $1 + $2 * $2)" Dimensions="', I10,'">',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_s"]/DataItem[1]',/&
    '</DataItem>',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_p"]/DataItem[1]',/&
    '</DataItem>',/&
    '<DataItem Reference="XML">',/&
    '                    /Xdmf/Domain/Grid[@Name="CellsTime"]/Grid[@Name="', A,'"]/Attribute[@Name="u_z"]/DataItem[1]',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="divergence(u)" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                      ', i10,  '         0  ',/&
    '                               1         1  ',/&
    '                               1',   i10,'  ',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/straintrace',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '<Attribute Name="curl(u)" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                      ', i10,  '         0  ',/&
    '                               1         1  ',/&
    '                               1',   i10,'  ',/&
    '</DataItem>',/&
    '<DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/curlinplane',/&
    '</DataItem>',/&
    '</DataItem>',/&
    '</Attribute>',/&
    '</Grid>',/)

    close(110)

end subroutine glob_snapshot_xdmf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine calc_curlinplane(f_sol, curlinplane)

   use data_pointwise, only: inv_rho_fluid
   use pointwise_derivatives, only: axisym_gradient_solid, axisym_gradient_solid_add
   use pointwise_derivatives, only: axisym_gradient_fluid, axisym_gradient_fluid_add
   use pointwise_derivatives, only: f_over_s_solid, f_over_s_fluid
   use data_source, only: src_type

   real(kind=realkind), intent(in)  :: f_sol(0:,0:,:,:)
   real(kind=realkind), intent(out) :: curlinplane(0:,0:,:,:)

   real(kind=realkind)              :: grad_sol_s(0:npol,0:npol,nel_solid,2)
   real(kind=realkind)              :: grad_sol_z(0:npol,0:npol,nel_solid,2)

   curlinplane = 0
   if (src_type(1) == 'dipole') then
      call axisym_gradient_solid(f_sol(:,:,:,1) + f_sol(:,:,:,2), grad_sol_s)
   else
      call axisym_gradient_solid(f_sol(:,:,:,1), grad_sol_s) ! 1: dsus, 2: dzus
   endif

   call axisym_gradient_solid(f_sol(:,:,:,3), grad_sol_z) ! 1:dsuz 2:dzuz

   curlinplane(:,:,nel_fluid+1:nel_fluid+nel_solid,1) &
              = grad_sol_s(:,:,:,2) - grad_sol_z(:,:,:,1)


end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine calc_straintrace(f_sol,chi,straintrace)
! < Calculate strain trace (for P-wave visualisation)

   use data_pointwise, only: inv_rho_fluid
   use pointwise_derivatives, only: axisym_gradient_solid, axisym_gradient_solid_add
   use pointwise_derivatives, only: axisym_gradient_fluid, axisym_gradient_fluid_add
   use pointwise_derivatives, only: f_over_s_solid, f_over_s_fluid
   use data_source, only: src_type

   real(kind=realkind), intent(in)  :: f_sol(0:,0:,:,:), chi(0:,0:,:)
   real(kind=realkind), intent(out) :: straintrace(0:,0:,:,:)

   real(kind=realkind)              :: grad_sol(0:npol,0:npol,nel_solid,2)
   real(kind=realkind)              :: buff_solid(0:npol,0:npol,nel_solid)
   real(kind=realkind)              :: usz_fluid(0:npol,0:npol,nel_fluid,2)
   real(kind=realkind)              :: grad_flu(0:npol,0:npol,nel_fluid,2)
   real(kind=realkind)              :: two_rk = 2


   if (src_type(1) == 'dipole') then
      call axisym_gradient_solid(f_sol(:,:,:,1) + f_sol(:,:,:,2), grad_sol)
   else
      call axisym_gradient_solid(f_sol(:,:,:,1), grad_sol) ! 1: dsus, 2: dzus
   endif

   call axisym_gradient_solid_add(f_sol(:,:,:,3), grad_sol) ! 1:dsuz+dzus,2:dzuz+dsus

   if (src_type(1) == 'monopole') then
      buff_solid = f_over_s_solid(f_sol(:,:,:,1))
   else if (src_type(1) == 'dipole') then
      buff_solid = two_rk * f_over_s_solid(f_sol(:,:,:,2))
   else if (src_type(1) == 'quadpole') then
      buff_solid = f_over_s_solid(f_sol(:,:,:,1) - two_rk * f_sol(:,:,:,2))
   endif
   straintrace(:,:,nel_fluid+1:nel_fluid+nel_solid,1) = buff_solid + grad_sol(:,:,:,2)

   if (have_fluid) then
      ! construct displacements in the fluid
      call axisym_gradient_fluid(chi, usz_fluid)
      usz_fluid(:,:,:,1) = usz_fluid(:,:,:,1) * inv_rho_fluid
      usz_fluid(:,:,:,2) = usz_fluid(:,:,:,2) * inv_rho_fluid

      ! gradient of s component
      call axisym_gradient_fluid(usz_fluid(:,:,:,1), grad_flu)   ! 1:dsus, 2:dzus

      ! gradient of z component added to s-comp gradient for strain trace and E13
      call axisym_gradient_fluid_add(usz_fluid(:,:,:,2), grad_flu)   !1:dsuz+dzus
                                                                     !2:dzuz+dsus

      straintrace(:,:,1:nel_fluid,1) = f_over_s_fluid(usz_fluid(:,:,:,1)) &
                                        + grad_flu(:,:,:,2)
   endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine xdmf_mapping(u_in, mapping_ijel_iplot, plotting_mask, i_arr_xdmf, j_arr_xdmf, &
                        u_out)

   real(kind=realkind), intent(in)    :: u_in(0:,0:,:,:)
   integer,             intent(in)    :: mapping_ijel_iplot(:,:,:)
   logical,             intent(in)    :: plotting_mask(:,:,:)
   integer,             intent(in)    :: i_arr_xdmf(:), j_arr_xdmf(:)
   real(kind=realkind), intent(inout) :: u_out(:,:)

   integer                            :: i_n_xdmf, j_n_xdmf, i, j
   integer                            :: nelem
   integer                            :: iel, ipol, ipol1, jpol, jpol1, ct

   nelem = size(u_in,3)

   i_n_xdmf = size(mapping_ijel_iplot,1)
   j_n_xdmf = size(mapping_ijel_iplot,2)
   do iel=1, nelem
      do i=1, i_n_xdmf - 1
         ipol = i_arr_xdmf(i)
         ipol1 = i_arr_xdmf(i+1)

         do j=1, j_n_xdmf - 1
            jpol = j_arr_xdmf(j)
            jpol1 = j_arr_xdmf(j+1)

            if (plotting_mask(i,j,iel)) then
               ct = mapping_ijel_iplot(i,j,iel)
               u_out(:, ct) = u_in(ipol,jpol,iel,:)
            endif

            if (plotting_mask(i+1,j,iel)) then
               ct = mapping_ijel_iplot(i+1,j,iel)
               u_out(:, ct) = u_in(ipol1,jpol,iel,:)
            endif

            if (plotting_mask(i+1,j+1,iel)) then
               ct = mapping_ijel_iplot(i+1,j+1,iel)
               u_out(:, ct) = u_in(ipol1,jpol1,iel,:)
            endif

            if (plotting_mask(i,j+1,iel)) then
               ct = mapping_ijel_iplot(i,j+1,iel)
               u_out(:, ct) = u_in(ipol,jpol1,iel,:)
            endif
         enddo
      enddo
   enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function kwf_mapping_sol(u_in)

   real(kind=realkind), intent(in)    :: u_in(0:npol,0:npol,nel_solid)
   real(kind=realkind)                :: kwf_mapping_sol(npoint_solid_kwf)

   integer                            :: iel, ipol, jpol, ct

   do iel=1, nel_solid
      do ipol=0, npol
         do jpol=0, npol
            if (kwf_mask(ipol,jpol,iel)) then
               ct = mapping_ijel_ikwf(ipol,jpol,iel)
               kwf_mapping_sol(ct) = u_in(ipol,jpol,iel)
            endif
         enddo
      enddo
   enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function kwf_mapping_flu(u_in)

   real(kind=realkind), intent(in)    :: u_in(0:npol,0:npol,nel_fluid)
   real(kind=realkind)                :: kwf_mapping_flu(npoint_fluid_kwf)

   integer                            :: iel, ipol, jpol, ct

   do iel=1, nel_fluid
      do ipol=0, npol
         do jpol=0, npol
            if (kwf_mask(ipol,jpol,iel + nel_solid)) then
               ct = mapping_ijel_ikwf(ipol,jpol,iel + nel_solid) - npoint_solid_kwf
               kwf_mapping_flu(ct) = u_in(ipol,jpol,iel)
            endif
         enddo
      enddo
   enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_field_1d(f, filename, appisnap, n)

   use data_source, only: have_src, src_dump_type
   use data_mesh, only: nel_solid, nel_fluid

   integer, intent(in)                 :: n
   real(kind=realkind),intent(in)      :: f(0:,0:,:)
   character(len=16), intent(in)       :: filename
   character(len=4), intent(in)        :: appisnap

   real(kind=realkind)                 :: floc(0:size(f,1)-1, 0:size(f,2)-1, 1:size(f,3))

   floc = f

   if (src_dump_type == 'mask' .and. n == nel_solid) &
        call eradicate_src_elem_values(floc)

   if (use_netcdf) then

       if (n == nel_solid) then
          if (dump_type == 'strain_only') then
             if (npoint_solid_kwf > 0) &
                call nc_dump_field_solid(kwf_mapping_sol(floc), filename(2:))
          else
             call nc_dump_field_solid(pack(floc(ibeg:iend,jbeg:jend,:), .true.), &
                                    filename(2:))
          endif

       else if (n == nel_fluid) then
          if (dump_type == 'strain_only') then
             if (npoint_fluid_kwf > 0) &
                call nc_dump_field_fluid(kwf_mapping_flu(floc), filename(2:))
          else
             call nc_dump_field_fluid(pack(floc(ibeg:iend,jbeg:jend,:), .true.), &
                                    filename(2:))
          endif

       else
           write(*,*) 'Neither solid nor fluid. What''s wrong here?'
           stop 2
       endif
   else
      open(unit=25000+mynum, file=datapath(1:lfdata)//filename//'_' &
                                  //appmynum//'_'//appisnap//'.bindat', &
           FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
      write(25000+mynum) pack(floc(ibeg:iend,jbeg:jend,:), .true.)
      close(25000+mynum)
   endif

end subroutine dump_field_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_disp(u, chi, istrain)

   use data_source, only: src_type,src_dump_type

   real(kind=realkind), intent(in) :: u(0:,0:,:,:)
   real(kind=realkind), intent(in) :: chi(0:,0:,:)
   integer,             intent(in) :: istrain

   integer                         :: i
   character(len=4)                :: appisnap
   real(kind=realkind)             :: f(0:size(u,1)-1, 0:size(u,1)-1,size(u,3),3)

   call define_io_appendix(appisnap,istrain)

   f = u

   if (src_dump_type == 'mask') then
      call eradicate_src_elem_vec_values(f)
   endif

   ! Dump solid displacement
   open(unit=75000+mynum,file=datapath(1:lfdata)//'/disp_sol_'&
                             //appmynum//'_'//appisnap//'.bindat', &
                             FORM="UNFORMATTED",STATUS="REPLACE")

   if (src_type(1) /= 'monopole') then
      write(75000+mynum) (f(ibeg:iend,jbeg:jend,:,i),i=1,3)
   else
      write(75000+mynum) f(ibeg:iend,jbeg:jend,:,1:3:2)
   endif
   close(75000+mynum)

   ! Dump fluid potential
   if (have_fluid) then
      open(unit=76000+mynum,file=datapath(1:lfdata)//'/chi_flu_'&
                                //appmynum//'_'//appisnap//'.bindat', &
                                FORM="UNFORMATTED",STATUS="REPLACE")

      write(76000+mynum) chi
      close(76000+mynum)
   endif

end subroutine dump_disp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_velo_dchi(v, dchi, istrain)

   use data_source, only: src_type, src_dump_type

   real(kind=realkind),intent(in) :: v(0:,0:,:,:)
   real(kind=realkind),intent(in) :: dchi(0:,0:,:)
   integer,            intent(in) :: istrain

   integer                        :: i
   character(len=4)               :: appisnap
   real(kind=realkind)            :: f(0:size(v,1)-1, 0:size(v,1)-1,size(v,3),3)

   call define_io_appendix(appisnap,istrain)

   f = v

   if (src_dump_type == 'mask') then
      call eradicate_src_elem_vec_values(f)
   endif

   ! Dump solid velocity vector
   open(unit=85000+mynum,file=datapath(1:lfdata)//'/velo_sol_'&
                              //appmynum//'_'//appisnap//'.bindat', &
                              FORM="UNFORMATTED",STATUS="REPLACE")

   if (src_type(1) /= 'monopole') then
      write(85000+mynum) (f(ibeg:iend,jbeg:jend,:,i), i=1,3)
   else
      write(85000+mynum) f(ibeg:iend,jbeg:jend,:,1), f(ibeg:iend,jbeg:jend,:,3)
   endif
   close(85000+mynum)

   ! Dump fluid potential 1st derivative
   if (have_fluid) then
      open(unit=86000+mynum,file=datapath(1:lfdata)//'/dchi_flu_'&
                                 //appmynum//'_'//appisnap//'.bindat', &
                                 FORM="UNFORMATTED",STATUS="REPLACE")

      write(86000+mynum) dchi
      close(86000+mynum)
   endif

end subroutine dump_velo_dchi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_velo_global(v, dchi, istrain)

   use data_pointwise, only: inv_rho_fluid
   use data_source, only: src_type, src_dump_type
   use pointwise_derivatives, only: axisym_gradient_fluid, dsdf_fluid_allaxis
   use data_mesh, only: npol, nel_solid, nel_fluid

   real(kind=realkind), intent(in) :: v(:,:,:,:)
   real(kind=realkind), intent(in) :: dchi(:,:,:)
   integer,             intent(in) :: istrain

   integer                         :: i
   character(len=4)                :: appisnap
   real(kind=realkind)             :: f(0:npol,0:npol,1:nel_solid,3)
   real(kind=realkind)             :: fflu(0:npol,0:npol,1:nel_fluid,3)
   real(kind=realkind)             :: usz_fluid(0:npol,0:npol,1:nel_fluid,2)

   call define_io_appendix(appisnap,istrain)

   ! sssssssssssss dump velocity vector inside solid ssssssssssssssssssssssssssss

   if (src_type(1) == 'dipole') then
      ! tranform from +/- coordinates to s/phi
      f(:,:,:,1) = v(:,:,:,1) + v(:,:,:,2)
      f(:,:,:,2) = v(:,:,:,1) - v(:,:,:,2)
      f(:,:,:,3) = v(:,:,:,3)
   else
      f = v
   endif

   if (src_dump_type == 'mask') then
      call eradicate_src_elem_vec_values(f)
   endif

   if (use_netcdf) then
      if (src_type(1) /= 'monopole') then
         call nc_dump_field_solid(pack(f(ibeg:iend,jbeg:jend,:,2),.true.), 'velo_sol_p')
      endif
      call nc_dump_field_solid(pack(f(ibeg:iend,jbeg:jend,:,1),.true.), 'velo_sol_s')
      call nc_dump_field_solid(pack(f(ibeg:iend,jbeg:jend,:,3),.true.), 'velo_sol_z')
   else
      open(unit=95000+mynum,file=datapath(1:lfdata)//'/velo_sol_'&
                                 //appmynum//'_'//appisnap//'.bindat', &
                                 FORM="UNFORMATTED",STATUS="REPLACE")
      if (src_type(1) /= 'monopole') then
         write(95000+mynum) (f(ibeg:iend,jbeg:jend,:,i), i=1,3)
      else
         write(95000+mynum) f(ibeg:iend,jbeg:jend,:,1), &
                            f(ibeg:iend,jbeg:jend,:,3)
      endif
      close(95000+mynum)
   endif

   ! ffffffff fluid region ffffffffffffffffffffffffffffffffffffffffffffffffffffff

   if (have_fluid) then
      ! compute velocity vector inside fluid
     call axisym_gradient_fluid(dchi, usz_fluid)

     call define_io_appendix(appisnap,istrain)
     fflu(ibeg:iend,jbeg:jend,:,1) = inv_rho_fluid(ibeg:iend,jbeg:jend,:) * &
                                     usz_fluid(jbeg:jend,jbeg:jend,:,1)
     fflu(ibeg:iend,jbeg:jend,:,2) = 0
     fflu(ibeg:iend,jbeg:jend,:,3) = inv_rho_fluid(ibeg:iend,jbeg:jend,:) * &
                                     usz_fluid(ibeg:iend,jbeg:jend,:,2)

     ! dump velocity vector inside fluid
     if (use_netcdf) then
        call nc_dump_field_fluid(pack(fflu(ibeg:iend,jbeg:jend,:,1), .true.), 'velo_flu_s')
        call nc_dump_field_fluid(pack(fflu(ibeg:iend,jbeg:jend,:,3), .true.), 'velo_flu_z')
        if (src_type(1) /= 'monopole') then
          call nc_dump_field_fluid(pack(fflu(ibeg:iend,jbeg:jend,:,2), .true.), 'velo_flu_p')
        endif
     else
        open(unit=960000+mynum,file=datapath(1:lfdata)//'/velo_flu_'&
                                  //appmynum//'_'//appisnap//'.bindat', &
                                   FORM="UNFORMATTED",STATUS="REPLACE")

        write(960000+mynum) (fflu(ibeg:iend,jbeg:jend,:,i), i=1,3)
        close(960000+mynum)
     endif ! netcdf
   endif ! have_fluid

end subroutine dump_velo_global
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_disp_global(u, chi, istrain)

   use data_pointwise, only: inv_rho_fluid
   use data_source, only: src_type, src_dump_type
   use pointwise_derivatives, only: axisym_gradient_fluid, dsdf_fluid_allaxis
   use data_mesh, only: npol, nel_solid, nel_fluid

   real(kind=realkind), intent(in) :: u(:,:,:,:)
   real(kind=realkind), intent(in) :: chi(:,:,:)
   integer,             intent(in) :: istrain

   integer                         :: i
   character(len=4)                :: appisnap
   real(kind=realkind)             :: f(0:npol,0:npol,1:nel_solid,3)
   real(kind=realkind)             :: fflu(0:npol,0:npol,1:nel_fluid,3)
   real(kind=realkind)             :: usz_fluid(0:npol,0:npol,1:nel_fluid,2)

   call define_io_appendix(appisnap, istrain)

   ! sssssssssssss dump disp vector inside solid ssssssssssssssssssssssssssss

   if (src_type(1) == 'dipole') then
      ! tranform from +/- coordinates to s/phi
      f(:,:,:,1) = u(:,:,:,1) + u(:,:,:,2)
      f(:,:,:,2) = u(:,:,:,1) - u(:,:,:,2)
      f(:,:,:,3) = u(:,:,:,3)
   else
      f = u
   endif

   if (src_dump_type == 'mask') then
      call eradicate_src_elem_vec_values(f)
   endif

   if (use_netcdf .and. dump_type == 'displ_only') then
      if (src_type(1) /= 'monopole') then
         call nc_dump_field_solid(kwf_mapping_sol(f(:,:,:,2)), 'disp_sol_p')
      endif
      call nc_dump_field_solid(kwf_mapping_sol(f(:,:,:,1)), 'disp_sol_s')
      call nc_dump_field_solid(kwf_mapping_sol(f(:,:,:,3)), 'disp_sol_z')

   else if (use_netcdf .and. dump_type /= 'displ_only') then
      if (src_type(1) /= 'monopole') then
         call nc_dump_field_solid(pack(f(ibeg:iend,jbeg:jend,:,2),.true.), 'disp_sol_p')
      endif
      call nc_dump_field_solid(pack(f(ibeg:iend,jbeg:jend,:,1),.true.), 'disp_sol_s')
      call nc_dump_field_solid(pack(f(ibeg:iend,jbeg:jend,:,3),.true.), 'disp_sol_z')

   else
      open(unit=95000+mynum,file=datapath(1:lfdata)//'/disp_sol_'&
                                 //appmynum//'_'//appisnap//'.bindat', &
                                 FORM="UNFORMATTED",STATUS="REPLACE")
      if (src_type(1) /= 'monopole') then
         write(95000+mynum) (f(ibeg:iend,jbeg:jend,:,i), i=1,3)
      else
         write(95000+mynum) f(ibeg:iend,jbeg:jend,:,1), &
                            f(ibeg:iend,jbeg:jend,:,3)
      endif
      close(95000+mynum)
   endif

   ! ffffffff fluid region ffffffffffffffffffffffffffffffffffffffffffffffffffffff

   if (have_fluid) then
      ! compute velocity vector inside fluid
     call axisym_gradient_fluid(chi, usz_fluid)

     call define_io_appendix(appisnap,istrain)
     fflu(:,:,:,1) = inv_rho_fluid(:,:,:) * usz_fluid(:,:,:,1)
     fflu(:,:,:,2) = 0
     fflu(:,:,:,3) = inv_rho_fluid(:,:,:) * usz_fluid(:,:,:,2)

     ! dump displacement vector inside fluid
     if (use_netcdf .and. dump_type == 'displ_only') then
        if (src_type(1) /= 'monopole') then
           call nc_dump_field_fluid(kwf_mapping_flu(fflu(:,:,:,2)), 'disp_flu_p')
        endif
        call nc_dump_field_fluid(kwf_mapping_flu(fflu(:,:,:,1)), 'disp_flu_s')
        call nc_dump_field_fluid(kwf_mapping_flu(fflu(:,:,:,3)), 'disp_flu_z')

     else if (use_netcdf .and. dump_type /= 'displ_only') then
        call nc_dump_field_fluid(pack(fflu(ibeg:iend,jbeg:jend,:,1), .true.), 'disp_flu_s')
        call nc_dump_field_fluid(pack(fflu(ibeg:iend,jbeg:jend,:,3), .true.), 'disp_flu_z')
        if (src_type(1) /= 'monopole') then
          call nc_dump_field_fluid(pack(fflu(ibeg:iend,jbeg:jend,:,2), .true.), 'disp_flu_p')
        endif
     else
        open(unit=960000+mynum,file=datapath(1:lfdata)//'/disp_flu_'&
                                  //appmynum//'_'//appisnap//'.bindat', &
                                   FORM="UNFORMATTED",STATUS="REPLACE")

        write(960000+mynum) (fflu(ibeg:iend,jbeg:jend,:,i), i=1,3)
        close(960000+mynum)
     endif ! netcdf
   endif ! have_fluid

end subroutine dump_disp_global
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Deletes all entries to vector field u on ALL GLL points inside
!! elements that have a non-zero source term (i.e. including all
!! assembled neighboring elements)
!! This is a preliminary test for the wavefield dumps.
pure subroutine eradicate_src_elem_vec_values(u)

   use data_source, only: nelsrc, ielsrc

   real(kind=realkind), intent(inout) :: u(0:,0:,:,:)
   integer                            :: iel

   do iel = 1, nelsrc
      u(:,:,ielsrc(iel),:) = 0.0
   enddo

end subroutine eradicate_src_elem_vec_values
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Deletes all entries to scalar field u on ALL GLL points inside
!! elements that have a non-zero source term (i.e. including all
!! assembled neighboring elements)
!! This is a preliminary test for the wavefield dumps.
pure subroutine eradicate_src_elem_values(u)

   use data_source, only: nelsrc, ielsrc

   real(kind=realkind), intent(inout) :: u(0:,0:,:)
   integer                            :: iel

   do iel = 1, nelsrc
      u(:,:,ielsrc(iel)) = 0.0
   enddo

end subroutine eradicate_src_elem_values
!-----------------------------------------------------------------------------------------

end module wavefields_io
!=========================================================================================
