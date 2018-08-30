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
!> Contains some routines to ease interaction with netcdf files
module nc_snapshots

#ifdef enable_netcdf
    use netcdf
#endif
    use nc_helpers
    use global_parameters
    use data_io, only: verbose, deflate_level
    use data_proc, only: mynum, lpr

    implicit none
    private
    save

    !! Variables for dumping of wavefields for plotting purposes
    integer             :: nc_snap_disp_varid, nc_coord_dimid
    integer             :: nc_snap_point_varid, nc_snap_grid_varid
    integer             :: nc_snap_pwave_varid, nc_snap_swave_varid
    integer             :: ncid_out_snap
    integer             :: ndim_disp ! < 2 for monopole, 3 for rest


    public :: nc_dump_snap_points
    public :: nc_dump_snap_grid
    public :: nc_make_snapfile
    public :: nc_close_snapfile
    public :: nc_dump_snapshot

contains

!-----------------------------------------------------------------------------------------
subroutine nc_make_snapfile

    use data_mesh, only: npoint_plot, nelem_plot
    use data_proc, only: appmynum
    use data_io, only: datapath, lfdata, nsnap
    use data_source, only: src_type

#ifdef enable_netcdf
    integer              :: nmode, nc_snappoint_dimid, nc_snapelem_dimid, nc_snapdim_dimid
    integer              :: nc_snaptime_dimid, nc_snapconnec_dimid
    character(len=120)   :: fname

    if (lpr .and. verbose > 1) write(*,*) '   .... preparing xdmf nc file'

    if (src_type(1) == 'monopole') then
        ndim_disp = 2
    else
        ndim_disp = 3
    endif

    fname = datapath(1:lfdata) // '/netcdf_snap_' // appmynum // '.nc'
    nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
    call check(nf90_create(path=fname, cmode=nmode, ncid=ncid_out_snap) )

    call check(nf90_def_dim(ncid_out_snap, 'points', npoint_plot, nc_snappoint_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'elements', nelem_plot, nc_snapelem_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'dimensions', ndim_disp , nc_snapdim_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 's-z-coordinate', 2 , nc_coord_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'connections', 4 , nc_snapconnec_dimid) )
    call check(nf90_def_dim(ncid_out_snap, 'timesteps', nsnap , nc_snaptime_dimid) )

    call flush(6)
    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'displacement', &
                            xtype  = NF90_FLOAT, &
                            dimids = [nc_snapdim_dimid, nc_snappoint_dimid, &
                                      nc_snaptime_dimid], &
                            chunksizes = [ndim_disp, npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_disp_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'straintrace', &
                            xtype  = NF90_FLOAT, &
                            dimids = [nc_snappoint_dimid, nc_snaptime_dimid], &
                            chunksizes = [npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_pwave_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'curlinplane', &
                            xtype  = NF90_FLOAT, &
                            dimids = [nc_snappoint_dimid, nc_snaptime_dimid], &
                            chunksizes = [npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_swave_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'points', &
                            xtype  = NF90_FLOAT, &
                            dimids = [nc_coord_dimid, nc_snappoint_dimid], &
                            varid  = nc_snap_point_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'grid', &
                            xtype  = NF90_INT, &
                            dimids = [nc_snapconnec_dimid, nc_snapelem_dimid], &
                            varid  = nc_snap_grid_varid) )

    call check(nf90_enddef(ncid    = ncid_out_snap))

    if (lpr .and. verbose > 1) write(*,*) '   .... DONE'

#endif

end subroutine nc_make_snapfile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snap_points(points)

    use data_mesh, only: npoint_plot
    real(sp), dimension(2,npoint_plot), intent(in)       :: points

#ifdef enable_netcdf
    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_point_varid, &
                            count  = [2, npoint_plot], &
                            values = points) )
#endif

end subroutine nc_dump_snap_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snap_grid(grid)

    use data_mesh, only: nelem_plot
    integer, dimension(4, nelem_plot), intent(in)       :: grid

#ifdef enable_netcdf
    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_grid_varid, &
                            count  = [4, nelem_plot], &
                            values = grid) )
#endif
end subroutine nc_dump_snap_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snapshot(u, straintrace, curlinplane, isnap)

    use data_mesh, only: npoint_plot
    use data_io, only: nsnap
    use data_source, only: src_type

    real(kind=realkind), dimension(3,npoint_plot), intent(in)  :: u
    real(kind=realkind), dimension(1,npoint_plot), intent(in)  :: straintrace
    real(kind=realkind), dimension(1,npoint_plot), intent(in)  :: curlinplane
    integer,                                       intent(in)  :: isnap

#ifdef enable_netcdf
    if (src_type(1) == 'monopole') then
       call putvar_real3d(ncid   = ncid_out_snap, &
                          varid  = nc_snap_disp_varid, &
                          start  = [1, 1, isnap], &
                          count  = [1, npoint_plot, 1], &
                          values = reshape(u(1,:), [1, npoint_plot,1]) )

       call putvar_real3d(ncid   = ncid_out_snap, &
                          varid  = nc_snap_disp_varid, &
                          start  = [2, 1, isnap], &
                          count  = [1, npoint_plot, 1], &
                          values = reshape(u(3,:), [1, npoint_plot,1]) )
    else
       call putvar_real3d(ncid   = ncid_out_snap, &
                          varid  = nc_snap_disp_varid, &
                          start  = [1, 1, isnap], &
                          count  = [1, npoint_plot, 1], &
                          values = reshape(u(1,:), [1, npoint_plot,1]) )

       call putvar_real3d(ncid   = ncid_out_snap, &
                          varid  = nc_snap_disp_varid, &
                          start  = [2, 1, isnap], &
                          count  = [1, npoint_plot, 1], &
                          values = reshape(u(2,:), [1, npoint_plot,1]) )

       call putvar_real3d(ncid   = ncid_out_snap, &
                          varid  = nc_snap_disp_varid, &
                          start  = [3, 1, isnap], &
                          count  = [1, npoint_plot, 1], &
                          values = reshape(u(3,:), [1, npoint_plot,1]) )
    endif

    call check(nf90_put_var(ncid   = ncid_out_snap, &
                       varid  = nc_snap_pwave_varid, &
                       start  = [1, isnap], &
                       count  = [npoint_plot, 1], &
                       values = straintrace(1,:)) )

    call check(nf90_put_var(ncid   = ncid_out_snap, &
                       varid  = nc_snap_swave_varid, &
                       start  = [1, isnap], &
                       count  = [npoint_plot, 1], &
                       values = curlinplane(1,:)) )
#endif

end subroutine nc_dump_snapshot
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_close_snapfile

#ifdef enable_netcdf
        call check(nf90_close(ncid_out_snap))
#endif

end subroutine
!-----------------------------------------------------------------------------------------

end module nc_snapshots
!=========================================================================================
