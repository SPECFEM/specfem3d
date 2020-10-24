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
program field_transformation

#ifdef unc
    use netcdf
#endif
    implicit none

#ifdef unc
    include 'netcdf.inc'
#endif

#ifdef unc
    integer                         :: nvar, ivar
    integer                         :: nsnap, ngll, ngllread, npol, nelem
    integer                         :: nmode

    integer                         :: ncin_id, ncin_snap_grpid, ncin_mesh_grpid
    integer                         :: ncout_id, ncout_fields_grpid, ncout_gll_dimid
    integer                         :: ncout_snap_dimid, ncout_mesh_grpid
    integer                         :: ncout_mesh_varids(255), ncin_mesh_varids(255)
    integer                         :: ncout_mesh_mps_varid, ncin_mesh_mps_varid
    integer                         :: ncout_mesh_mpz_varid, ncin_mesh_mpz_varid
    integer                         :: ncout_surf_grpid, ncout_surfstrain_dimid, ncout_surf_dimid
    integer                         :: ncout_comp_dimid
    integer                         :: ncout_surf_varids(7), ncin_surf_grpid, ncin_surf_varids(7)

    integer                         :: nc_mesh_npol_dimid, nc_mesh_cntrlpts_dimid, &
                                       nc_mesh_elem_dimid
    integer                         :: ncout_mesh_sem_varid, ncout_mesh_fem_varid, &
                                       ncout_mesh_midpoint_varid
    integer                         :: ncin_mesh_sem_varid, ncin_mesh_fem_varid, &
                                       ncin_mesh_midpoint_varid

    integer                         :: ncin_mesh_G0_varid, ncout_mesh_G0_varid
    integer                         :: ncin_mesh_G1_varid, ncout_mesh_G1_varid
    integer                         :: ncin_mesh_G2_varid, ncout_mesh_G2_varid
    integer                         :: ncin_mesh_gll_varid, ncout_mesh_gll_varid
    integer                         :: ncin_mesh_glj_varid, ncout_mesh_glj_varid

    integer                         :: ncout_mesh_eltype_varid, ncout_mesh_axis_varid
    integer                         :: ncin_mesh_eltype_varid, ncin_mesh_axis_varid

    integer                         :: nsurfelem, ncomp, nstraincomp
    character(len=8)                :: sourcetype
    character(len=12)               :: dump_type
    integer                         :: dimids(2)
    character(len=16), allocatable  :: varnamelist(:), varname_surf(:)
    integer, dimension(9)           :: ncin_field_varid
    integer, dimension(9)           :: ncout_field_varid

    integer                         :: isfinalized, status, percent
    integer                         :: nstep, nvars_mesh
    integer                         :: attnum, nf_att_stat
    character(len=80)               :: attname, varname

    real, allocatable               :: data_mesh(:), data_surf_1d(:), data_surf_3d(:,:,:)
    integer, allocatable            :: int_data_1d(:), int_data_2d(:,:), int_data_3d(:,:,:)
    double precision, allocatable   :: dp_data_1d(:), dp_data_2d(:,:)

    real(kind=8), dimension(:,:), allocatable       :: datat, datat_t

    double precision                :: time_fft, time_i, time_o, tick, tack
    double precision                :: space_i, space_o

    logical                         :: verbose = .true.
    integer                         :: npointsperstep, cache_size, narg
    integer                         :: chunk_gll

    character(len=32)               :: cache_size_char
                                    ! < Contains chunk size in GLL points. Should be system
                                    !! int(disk_block_size / nsnap)
    integer, parameter              :: disk_block_size = 8192

    ! lossy compression: setting the fields to numerical zero before p-wave
    ! arrival helps bzip a lot for the compression. The fields are truncated to
    ! the significant digits below the maximum of each time trace.
    ! the idea goes back to point 8) in
    ! http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    ! lossy deflation does not make sense for ffted fields (no blocks of small
    ! numbers expected)
    logical, parameter              :: deflate = .true.
    integer, parameter              :: deflate_level = 2
    logical, parameter              :: deflate_lossy = .false.
    integer, parameter              :: sigdigits =  5       ! significant digits
                                                            ! below max of time trace



    narg = command_argument_count()
    if (narg < 1) then
        print *, 'Warning: Argument "cache size" is missing, default: 1024 (MB)'
        cache_size = 1024
    else
        call get_command_argument(1, cache_size_char)
        read(cache_size_char, *) cache_size
        print '(A,I6,A)', 'Using ', cache_size, ' MB of memory for field_transformation'
    endif

    ! initialize timer
    time_fft = 0
    time_i = 0
    time_o = 0

    space_i = 0
    space_o = 0

    ! open input netcdf file
    call check( nf90_open(path="./Data/axisem_output.nc4", &
                          mode=NF90_NOWRITE, ncid=ncin_id) )

    status = nf90_get_att(ncin_id, NF90_GLOBAL, 'finalized', isfinalized)

    do while (isfinalized /= 1 .or. status /= NF90_NOERR)
      percent = 0
      status = nf90_get_att(ncin_id, NF90_GLOBAL, 'percent completed', percent)
      print "('Solver run not yet finished (at ', I3, '%). Waiting for 10s')", percent
      call check( nf90_close(ncin_id))
      call sleep(10)

      call check( nf90_open(path="./Data/axisem_output.nc4", &
                            mode=NF90_NOWRITE, ncid=ncin_id) )
      call check( nf90_get_att(ncin_id, NF90_GLOBAL, 'finalized', isfinalized) )
    enddo


    ! get Snapshots group id
    call check( nf90_inq_grp_ncid(ncin_id, "Snapshots", ncin_snap_grpid) )

    ! get excitation type (monopole or multipole?)
    call check( nf90_get_att(ncin_id, NF90_GLOBAL, "excitation type", sourcetype))

    call check( nf90_get_att(ncin_id, NF90_GLOBAL, "npol", npol))

    if (verbose) &
        print *, 'source type  ', sourcetype

    ! get dump type
    call check( nf90_get_att(ncin_id, NF90_GLOBAL, "dump type (displ_only, displ_velo, fullfields)", dump_type))

    if (verbose) &
        print *, 'dump type    ', dump_type

    if (trim(dump_type) == 'displ_only') then
       call check( nf90_get_att(ncin_id, NF90_GLOBAL, "nelem_kwf_global", nelem))
       if (sourcetype == 'monopole') then
           nvar = 2
           allocate(varnamelist(nvar))
           varnamelist = ['disp_s     ', 'disp_z     ']
       else
           nvar = 3
           allocate(varnamelist(nvar))
           varnamelist = ['disp_s     ', 'disp_p     ', 'disp_z     ']
       endif

    else if (trim(dump_type) == 'strain_only') then
       if (sourcetype == 'monopole') then
           nvar = 4
           allocate(varnamelist(nvar))
           varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                           'straintrace'/)
       else
           nvar = 6
           allocate(varnamelist(nvar))
           varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                           'strain_dsup', 'strain_dzup', 'straintrace'/)
       endif

    else if (trim(dump_type) == 'fullfields') then
       if (sourcetype == 'monopole') then
           nvar = 6
           allocate(varnamelist(nvar))
           varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                           'straintrace', 'velo_s     ', 'velo_z     '/)
       else
           nvar = 9
           allocate(varnamelist(nvar))
           varnamelist = (/'strain_dsus', 'strain_dsuz', 'strain_dpup', &
                           'strain_dsup', 'strain_dzup', 'straintrace', &
                           'velo_s     ', 'velo_p     ', 'velo_z     '/)
       endif
    else
       stop
    endif

    ! get variable ids of the fields
    do ivar=1, nvar
        call check( nf90_inq_varid(ncin_snap_grpid, varnamelist(ivar), &
                                   varid=ncin_field_varid(ivar)) )
    enddo

    ! get dimension ids (same for all fields)
    ivar = 1
    call check( nf90_inquire_variable(ncin_snap_grpid, ncin_field_varid(ivar), &
                                      dimids=dimids(:)) )

    ! get number of snapshots (same for all fields)
    call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(2), len=nsnap) )
    if (verbose) &
        print *, 'nsnap = ', nsnap

    ! get number of GLL points (same for all fields)
    call check( nf90_inquire_dimension(ncin_snap_grpid, dimids(1), len=ngll) )
    if (verbose) &
        print *, 'ngll  = ', ngll


    !! Create output file
    if (verbose) print *, 'Creating output file'
    nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
    call check( nf90_create(path="./Data/ordered_output.nc4", cmode=nmode, ncid=ncout_id))

    ! create group for timedomain fields
    call check( nf90_def_grp(ncout_id, "Snapshots", ncout_fields_grpid) )
    print *, 'Defined group "Snapshots"'

    ! copy attributes
    call check( nf90_copy_att( ncin_snap_grpid, NF90_GLOBAL, 'nstrain', &
                               ncout_fields_grpid, NF90_GLOBAL) )
    print *, 'Copyied Attributes'

    ! create dimensions
    call check( nf90_def_dim(ncid=ncout_id, name="gllpoints_all", &
                             len=ngll, dimid=ncout_gll_dimid) )

    call check( nf90_def_dim(ncid=ncout_id, name="snapshots", len=nsnap, &
                             dimid=ncout_snap_dimid) )

    print *, 'Defined snapshots dimension'

    chunk_gll = max(disk_block_size / nsnap, 1)
    print *, 'Chunksize: [', chunk_gll, ',', nsnap, ']'

    ! create variables
    do ivar=1, nvar
        call check( nf90_def_var(ncid       = ncout_fields_grpid, &
                                 name       = trim(varnamelist(ivar)), &
                                 xtype      = NF90_FLOAT, &
                                 dimids     = [ncout_gll_dimid, ncout_snap_dimid], &
                                 varid      = ncout_field_varid(ivar), &
                                 chunksizes = [chunk_gll, nsnap]) )

        call check( nf90_def_var_fill(ncid=ncout_fields_grpid, &
                                      varid=ncout_field_varid(ivar), &
                                      no_fill=1, fill=0) )

        call check( nf90_def_var_fletcher32(ncid=ncout_fields_grpid, &
                                            varid=ncout_field_varid(ivar), &
                                            fletcher32=1) )

        if (deflate) then
            call check( nf90_def_var_deflate(ncid=ncout_fields_grpid, &
                                             varid=ncout_field_varid(ivar), &
                                             shuffle=1, deflate=1, &
                                             deflate_level=deflate_level) )
        endif
    enddo

    ! Create mesh variables
    print *, 'Creating mesh variables'
    call check( nf90_inq_grp_ncid(ncin_id, "Mesh", ncin_mesh_grpid) )
    call check( nf90_inq_varids(ncid   = ncin_mesh_grpid, &
                                nvars  = nvars_mesh, &
                                varids = ncin_mesh_varids) )

    call check( nf90_def_grp(ncout_id, "Mesh", ncout_mesh_grpid) )

    if (trim(dump_type) == 'displ_only') then

       call check( nf90_def_dim( ncid   = ncout_mesh_grpid, &
                                 name   = 'elements', &
                                 len    = nelem, &
                                 dimid  = nc_mesh_elem_dimid) )

       call check( nf90_def_dim( ncid   = ncout_mesh_grpid, &
                                 name   = 'control_points', &
                                 len    = 4, &
                                 dimid  = nc_mesh_cntrlpts_dimid) )

       call check( nf90_def_dim( ncid   = ncout_mesh_grpid, &
                                 name   = 'npol', &
                                 len    = npol+1, &
                                 dimid  = nc_mesh_npol_dimid) )

       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'midpoint_mesh', &
                                 xtype  = NF90_INT, &
                                 dimids = nc_mesh_elem_dimid, &
                                 varid  = ncout_mesh_midpoint_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_midpoint_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'eltype', &
                                 xtype  = NF90_INT, &
                                 dimids = nc_mesh_elem_dimid, &
                                 varid  = ncout_mesh_eltype_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_eltype_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'axis', &
                                 xtype  = NF90_INT, &
                                 dimids = nc_mesh_elem_dimid, &
                                 varid  = ncout_mesh_axis_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_axis_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'fem_mesh', &
                                 xtype  = NF90_INT, &
                                 dimids = [nc_mesh_cntrlpts_dimid, &
                                           nc_mesh_elem_dimid], &
                                 varid  = ncout_mesh_fem_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_fem_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'sem_mesh', &
                                 xtype  = NF90_INT, &
                                 dimids = [nc_mesh_npol_dimid, &
                                           nc_mesh_npol_dimid, &
                                           nc_mesh_elem_dimid], &
                                 varid  = ncout_mesh_sem_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_sem_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid       = ncout_mesh_grpid, &
                                 name       = 'mp_mesh_S', &
                                 xtype      = NF90_FLOAT, &
                                 dimids     = [nc_mesh_elem_dimid], &
                                 varid      = ncout_mesh_mps_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_mps_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid       = ncout_mesh_grpid, &
                                 name       = 'mp_mesh_Z', &
                                 xtype      = NF90_FLOAT, &
                                 dimids     = [nc_mesh_elem_dimid], &
                                 varid      = ncout_mesh_mpz_varid) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                           varid = ncout_mesh_mpz_varid, &
                                           fletcher32=1) )

       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'G0', &
                                 xtype  = NF90_DOUBLE, &
                                 dimids = nc_mesh_npol_dimid, &
                                 varid  = ncout_mesh_G0_varid) )
       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'G1', &
                                 xtype  = NF90_DOUBLE, &
                                 dimids = [nc_mesh_npol_dimid, &
                                           nc_mesh_npol_dimid], &
                                 varid  = ncout_mesh_G1_varid) )
       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'G2', &
                                 xtype  = NF90_DOUBLE, &
                                 dimids = [nc_mesh_npol_dimid, &
                                           nc_mesh_npol_dimid], &
                                 varid  = ncout_mesh_G2_varid) )
       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'gll', &
                                 xtype  = NF90_DOUBLE, &
                                 dimids = nc_mesh_npol_dimid, &
                                 varid  = ncout_mesh_gll_varid) )
       call check( nf90_def_var( ncid   = ncout_mesh_grpid, &
                                 name   = 'glj', &
                                 xtype  = NF90_DOUBLE, &
                                 dimids = nc_mesh_npol_dimid, &
                                 varid  = ncout_mesh_glj_varid) )
    endif


    do ivar = 1, nvars_mesh
        call check( nf90_inquire_variable(ncid  = ncin_mesh_grpid, &
                                          varid = ncin_mesh_varids(ivar), &
                                          name  = varname ))

        !print *, 'Found mesh variable: ', trim(varname)

        if (varname(1:5) /= 'mesh_') cycle

        call check( nf90_def_var( ncid       = ncout_mesh_grpid, &
                                  name       = varname, &
                                  xtype      = NF90_FLOAT, &
                                  dimids     = [ncout_gll_dimid], &
                                  !chunksizes = [ngll], &
                                  varid      = ncout_mesh_varids(ivar)) )

        !call check( nf90_def_var_deflate( ncid    = ncout_mesh_grpid, &
        !                                  varid   = ncout_mesh_varids(ivar), &
        !                                  shuffle = 1, deflate = 1, &
        !                                  deflate_level = deflate_level) )

        call check( nf90_def_var_fletcher32(ncid  = ncout_mesh_grpid, &
                                            varid = ncout_mesh_varids(ivar), &
                                            fletcher32=1) )
    enddo

    ! Create Surface variables
    print *, 'Creating surface variables'
    nstraincomp = 6
    ncomp       = 3
    call check( nf90_inq_grp_ncid(ncin_id, "Surface", ncin_surf_grpid) )
    call check( nf90_def_grp(ncout_id, "Surface", ncout_surf_grpid) )
    call check( nf90_def_dim( ncid  = ncout_surf_grpid, &
                              name  = "straincomponents", &
                              len   = nstraincomp, &
                              dimid = ncout_surfstrain_dimid) )
    call check( nf90_def_dim( ncid  = ncout_surf_grpid, &
                              name  = "components", &
                              len   = ncomp, &
                              dimid = ncout_comp_dimid) )
    call check( nf90_get_att( ncid   = ncin_surf_grpid, &
                              name   = 'nsurfelem', &
                              varid  = NF90_GLOBAL, &
                              values = nsurfelem) )
    call check( nf90_def_dim( ncid   = ncout_surf_grpid, &
                              name   = "surf_elems", &
                              len    = nsurfelem, &
                              dimid  = ncout_surf_dimid) )

    allocate(varname_surf(7))
    varname_surf = ['elem_theta      ', 'displacement    ', 'velocity        ', &
                    'disp_src        ', 'strain          ', 'stf_dump        ', &
                    'stf_d_dump      ']
    do ivar = 1, size(varname_surf)
        call check( nf90_inq_varid(ncid  = ncin_surf_grpid, &
                                   varid = ncin_surf_varids(ivar), &
                                   name  = varname_surf(ivar) ))
        !print *, 'Found surface variable: ', trim(varname_surf(ivar))
    enddo

    call check( nf90_def_var(      ncid       = ncout_surf_grpid, &
                                   name       = varname_surf(1), &
                                   xtype      = NF90_FLOAT, &
                                   dimids     = [ncout_surf_dimid], &
                                   chunksizes = [nsurfelem], &
                                   varid      = ncout_surf_varids(1)) )

    call check( nf90_def_var_deflate( ncid    = ncout_surf_grpid, &
                                      varid   = ncout_surf_varids(1), &
                                      shuffle = 1, deflate = 1, &
                                      deflate_level = deflate_level) )

    call check( nf90_def_var_fletcher32(ncid  = ncout_surf_grpid, &
                                        varid = ncout_surf_varids(1), &
                                        fletcher32=1) )

    do ivar = 2, 4
        call check( nf90_def_var(ncid       = ncout_surf_grpid, &
                                 name       = varname_surf(ivar), &
                                 xtype      = NF90_FLOAT, &
                                 dimids     = [ncout_snap_dimid, ncout_comp_dimid, ncout_surf_dimid], &
                                 chunksizes = [nsnap, ncomp, 1], &
                                 varid      = ncout_surf_varids(ivar)) )

        call check( nf90_def_var_deflate( ncid    = ncout_surf_grpid, &
                                          varid   = ncout_surf_varids(ivar), &
                                          shuffle = 1, deflate = 1, &
                                          deflate_level = deflate_level) )

        call check( nf90_def_var_fletcher32(ncid  = ncout_surf_grpid, &
                                            varid = ncout_surf_varids(ivar), &
                                            fletcher32=1) )
    enddo

    call check( nf90_def_var(ncid       = ncout_surf_grpid, &
                             name       = varname_surf(5), &
                             xtype      = NF90_FLOAT, &
                             dimids     = [ncout_snap_dimid, ncout_surfstrain_dimid, ncout_surf_dimid], &
                             chunksizes = [nsnap, nstraincomp, 1], &
                             varid      = ncout_surf_varids(5)) )

    call check( nf90_def_var_deflate( ncid    = ncout_surf_grpid, &
                                      varid   = ncout_surf_varids(5), &
                                      shuffle = 1, deflate = 1, &
                                      deflate_level = deflate_level) )

    call check( nf90_def_var_fletcher32(ncid  = ncout_surf_grpid, &
                                        varid = ncout_surf_varids(5), &
                                        fletcher32=1) )

    do ivar = 6, 7
       call check( nf90_def_var(ncid       = ncout_surf_grpid, &
                                name       = varname_surf(ivar), &
                                xtype      = NF90_FLOAT, &
                                dimids     = [ncout_snap_dimid], &
                                chunksizes = [nsnap], &
                                varid      = ncout_surf_varids(ivar)) )

       call check( nf90_def_var_deflate( ncid    = ncout_surf_grpid, &
                                         varid   = ncout_surf_varids(ivar), &
                                         shuffle = 1, deflate = 1, &
                                         deflate_level = deflate_level) )

       call check( nf90_def_var_fletcher32(ncid  = ncout_surf_grpid, &
                                           varid = ncout_surf_varids(ivar), &
                                           fletcher32=1) )
    enddo

    print *, 'Surface variables defined'

    ! Copy all attributes
    nf_att_stat = NF90_NOERR
    attnum = 0
    do
        attnum = attnum + 1
        nf_att_stat = nf90_inq_attname(ncin_id, NF90_GLOBAL, attnum, attname)

        if (nf_att_stat /= NF90_NOERR) exit

        !print *, 'Copying attribute: ', attname
        call check(nf90_copy_att( ncid_in   = ncin_id, &
                                  varid_in  = NF90_GLOBAL, &
                                  name      = attname, &
                                  ncid_out  = ncout_id, &
                                  varid_out = NF90_GLOBAL))
    enddo

    call check( nf90_enddef(ncout_id))


    ! Copy mesh variables
    print *, 'Copying mesh variables'
    allocate(data_mesh(ngll))
    do ivar = 1, nvars_mesh
        call check( nf90_inquire_variable(ncid  = ncin_mesh_grpid, &
                                          varid = ncin_mesh_varids(ivar), &
                                          name  = varname ))
        if (varname(1:5) /= 'mesh_') cycle
        call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                  varid  = ncin_mesh_varids(ivar), &
                                  start  = [1], &
                                  count  = [ngll], &
                                  values = data_mesh) )
        call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                  varid  = ncout_mesh_varids(ivar), &
                                  start  = [1], &
                                  count  = [ngll], &
                                  values = data_mesh))
    enddo
    deallocate(data_mesh)

    if (trim(dump_type) == 'displ_only') then

       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_midpoint_varid, &
                                   name  = 'midpoint_mesh' ))

       allocate(int_data_1d(nelem))
       call check(nf90_get_var ( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_midpoint_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = int_data_1d))
       call check(nf90_put_var ( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_midpoint_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = int_data_1d))


       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_eltype_varid, &
                                   name  = 'eltype' ))

       call check(nf90_get_var ( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_eltype_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = int_data_1d))
       call check(nf90_put_var ( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_eltype_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = int_data_1d))


       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_axis_varid, &
                                   name  = 'axis' ))

       call check(nf90_get_var ( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_axis_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = int_data_1d))
       call check(nf90_put_var ( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_axis_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = int_data_1d))
       deallocate(int_data_1d)


       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_fem_varid, &
                                   name  = 'fem_mesh' ))

       allocate(int_data_2d(4, nelem))
       call check(nf90_get_var ( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_fem_varid, &
                                 start  = [1,1], &
                                 count  = [4,nelem], &
                                 values = int_data_2d))
       call check(nf90_put_var ( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_fem_varid, &
                                 start  = [1,1], &
                                 count  = [4,nelem], &
                                 values = int_data_2d))
       deallocate(int_data_2d)


       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_sem_varid, &
                                   name  = 'sem_mesh' ))

       allocate(int_data_3d(npol+1, npol+1, nelem))
       call check(nf90_get_var ( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_sem_varid, &
                                 start  = [1,1,1], &
                                 count  = [npol+1,npol+1,nelem], &
                                 values = int_data_3d))
       call check(nf90_put_var ( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_sem_varid, &
                                 start  = [1,1,1], &
                                 count  = [npol+1,npol+1,nelem], &
                                 values = int_data_3d))
       deallocate(int_data_3d)


       allocate(data_mesh(nelem))
       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_mps_varid, &
                                   name  = 'mp_mesh_S' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_mps_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = data_mesh) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_mps_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = data_mesh))

       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_mpz_varid, &
                                   name  = 'mp_mesh_Z' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_mpz_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = data_mesh) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_mpz_varid, &
                                 start  = [1], &
                                 count  = [nelem], &
                                 values = data_mesh))
       deallocate(data_mesh)

       allocate(dp_data_1d(0:npol))

       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_gll_varid, &
                                   name  = 'gll' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_gll_varid, &
                                 start  = [1], &
                                 count  = [npol+1], &
                                 values = dp_data_1d) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_gll_varid, &
                                 start  = [1], &
                                 count  = [npol+1], &
                                 values = dp_data_1d))

       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_glj_varid, &
                                   name  = 'glj' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_glj_varid, &
                                 start  = [1], &
                                 count  = [npol+1], &
                                 values = dp_data_1d) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_glj_varid, &
                                 start  = [1], &
                                 count  = [npol+1], &
                                 values = dp_data_1d))

       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_G0_varid, &
                                   name  = 'G0' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_G0_varid, &
                                 start  = [1], &
                                 count  = [npol+1], &
                                 values = dp_data_1d) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_G0_varid, &
                                 start  = [1], &
                                 count  = [npol+1], &
                                 values = dp_data_1d))

       deallocate(dp_data_1d)

       allocate(dp_data_2d(0:npol,0:npol))
       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_G1_varid, &
                                   name  = 'G1' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_G1_varid, &
                                 start  = [1,1], &
                                 count  = [npol+1, npol+1], &
                                 values = dp_data_2d) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_G1_varid, &
                                 start  = [1,1], &
                                 count  = [npol+1, npol+1], &
                                 values = dp_data_2d))

       call check( nf90_inq_varid( ncid  = ncin_mesh_grpid, &
                                   varid = ncin_mesh_G2_varid, &
                                   name  = 'G2' ))

       call check( nf90_get_var( ncid   = ncin_mesh_grpid, &
                                 varid  = ncin_mesh_G2_varid, &
                                 start  = [1,1], &
                                 count  = [npol+1, npol+1], &
                                 values = dp_data_2d) )

       call check( nf90_put_var( ncid   = ncout_mesh_grpid, &
                                 varid  = ncout_mesh_G2_varid, &
                                 start  = [1,1], &
                                 count  = [npol+1, npol+1], &
                                 values = dp_data_2d))
       deallocate(dp_data_2d)
    endif

    ! Done with the mesh

    ! Copy surface variables
    print *, 'Copying surface variables'
    allocate(data_surf_1d(nsurfelem))
    call check( nf90_get_var( ncid   = ncin_surf_grpid, &
                              varid  = ncin_surf_varids(1), &
                              start  = [1], &
                              count  = [nsurfelem], &
                              values = data_surf_1d) )
    call check( nf90_put_var( ncid   = ncout_surf_grpid, &
                              varid  = ncout_surf_varids(1), &
                              start  = [1], &
                              count  = [nsurfelem], &
                              values = data_surf_1d))
    deallocate(data_surf_1d)

    allocate(data_surf_3d(nsnap, ncomp, nsurfelem))
    do ivar = 2, 4
        call check( nf90_get_var( ncid   = ncin_surf_grpid, &
                                  varid  = ncin_surf_varids(ivar), &
                                  start  = [1,1,1], &
                                  count  = [nsnap, ncomp, nsurfelem], &
                                  values = data_surf_3d) )
        call check( nf90_put_var( ncid   = ncout_surf_grpid, &
                                  varid  = ncout_surf_varids(ivar), &
                                  start  = [1,1,1], &
                                  count  = [nsnap, ncomp, nsurfelem], &
                                  values = data_surf_3d))
    enddo
    deallocate(data_surf_3d)

    allocate(data_surf_3d(nsnap, nstraincomp, nsurfelem))
    call check( nf90_get_var( ncid   = ncin_surf_grpid, &
                              varid  = ncin_surf_varids(5), &
                              start  = [1,1,1], &
                              count  = [nsnap, nstraincomp, nsurfelem], &
                              values = data_surf_3d) )
    call check( nf90_put_var( ncid   = ncout_surf_grpid, &
                              varid  = ncout_surf_varids(5), &
                              start  = [1,1,1], &
                              count  = [nsnap, nstraincomp, nsurfelem], &
                              values = data_surf_3d))
    deallocate(data_surf_3d)

    allocate(data_surf_1d(nsnap))
    do ivar = 6, 7
       call check( nf90_get_var( ncid   = ncin_surf_grpid, &
                                 varid  = ncin_surf_varids(ivar), &
                                 start  = [1], &
                                 count  = [nsnap], &
                                 values = data_surf_1d) )
       call check( nf90_put_var( ncid   = ncout_surf_grpid, &
                                 varid  = ncout_surf_varids(ivar), &
                                 start  = [1], &
                                 count  = [nsnap], &
                                 values = data_surf_1d))
    enddo
    deallocate(data_surf_1d)



    ! loop over fields
    do ivar=1, nvar
        if (verbose) &
            print *, varnamelist(ivar)
        ! loop over subsets of the GLL points
        nstep = 0
        do while (nstep + 1 < ngll)

            npointsperstep = cache_size * 1048576 / 4 / nsnap

            ngllread = min(npointsperstep, ngll - nstep)

            allocate(datat_t(1:ngllread, 1:nsnap))

            !print *, 'npointsperstep: ', npointsperstep, ', ngllread: ', ngllread
            ! read a chunk of data
            call cpu_time(tick)
            call check( nf90_get_var(ncin_snap_grpid, ncin_field_varid(ivar), &
                                     values=datat_t(:, :), &
                                     start=(/nstep+1, 1/), &
                                     count=(/ngllread, nsnap/)) )

            !datat(1:nsnap,:) = transpose(datat_t(:,1:nsnap))

            call cpu_time(tack)
            time_i = time_i + tack - tick
            space_i = space_i + ngllread * nsnap * 4 / 1048576.
            if (verbose) &
                print "('read  ', F12.2, ' MB in ', F7.2, ' s => ', F7.2, 'MB/s' )", &
                    real(ngllread) * nsnap * 4 / 1048576., tack-tick, &
                    real(ngllread) * nsnap * 4 / 1048576. / (tack-tick)

            ! trunkate for better compression
            !if (deflate .and. deflate_lossy) then
            !    call cpu_time(tick)
            !    call truncate(datat_t, sigdigits)
            !    call cpu_time(tack)
            !    time_fft = time_fft + tack - tick
            !endif

            ! write transposed data to output file
            call cpu_time(tick)
            call check( nf90_put_var(ncout_fields_grpid, ncout_field_varid(ivar), &
                                     values = datat_t(:, :), &
                                     start  = [nstep+1, 1], &
                                     count  = [ngllread, nsnap]) )
            call cpu_time(tack)
            time_o = time_o + tack - tick
            space_o = space_o + ngllread * nsnap * 4 / 1048576.
            if (verbose) &
                print "('wrote ', F12.2, ' MB in ', F7.2, ' s => ', F7.2, 'MB/s' )", &
                    real(ngllread) * nsnap * 4 / 1048576., tack-tick, &
                    real(ngllread) * nsnap * 4 / 1048576. / (tack-tick)

            deallocate(datat_t)

            nstep = nstep + ngllread
        enddo
    enddo

    call check( nf90_close(ncin_id))
    call check( nf90_close(ncout_id))

    call dump_mesh_data_xdmf(filename = 'Data/ordered_output.nc4', varname='Snapshots/straintrace', &
                             npoints=ngll, nsnap=nsnap)

    print '(A, F8.2, A)', 'Time spent for compression/fft: ', time_fft, ' s'
    print '(A, F8.2, A)', 'Time spent for I:               ', time_i,   ' s'
    print '(A, F8.2, A)', 'Time spent for O:               ', time_o,   ' s'
    print *
    print '(3(A, F8.2))', 'MB I: ', space_i, ' MB, av. speed: ', space_i / time_i, ' MB/s'
    print '(3(A, F8.2))', 'MB O: ', space_o, ' MB, av. speed: ', space_o / time_o, ' MB/s'

#else

    stop 'This program can only be run with NetCDF enabled'
#endif

contains

!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    implicit none
    integer, intent ( in) :: status ! < Error code
#ifdef unc
    if (status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        call abort()
    endif
#endif
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine truncate(dataIO, sigdigits)
    implicit none
    real(kind=8), intent(inout)     :: dataIO(:,:)
    integer, intent(in)             :: sigdigits
    integer                         :: bits
    double precision                :: scaleit, maxt
    integer                         :: n

    ! find next power of 2 corresponding to sigdigits
    bits = 1
    do while (log10(2.) * bits < real(sigdigits))
        bits = bits + 1
        scaleit = real(2**bits)
    enddo

    ! trunkate the data time series wise (each GLL point separately)
    do n = lbound(dataIO,2), ubound(dataIO,2)
        maxt = maxval(dataIO(:,n))
        dataIO(:,n) = real(nint(scaleit * dataIO(:,n) / maxt) / scaleit) * maxt
    enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_data_xdmf(filename, varname, npoints, nsnap)
  character(len=*), intent(in)      :: filename, varname
  integer, intent(in)               :: npoints, nsnap

  integer                           :: iinput_xdmf
  integer                           :: i
  character(len=512)                :: filename_np


  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 733) npoints, npoints, trim(filename_np), npoints, trim(filename_np)

  do i=1, nsnap
     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) dble(i), npoints, "'", "'"

     ! write attribute
     write(iinput_xdmf, 7342) varname, npoints, i-1, npoints, nsnap, npoints, &
                              trim(filename_np), trim(varname)

     write(iinput_xdmf, 7343)
  enddo

  ! finish xdmf file
  write(iinput_xdmf, 736)
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
    '</DataItem>',/,/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/)

7341 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Time Value="',F8.2,'" />',/&
    '<Topology TopologyType="Polyvertex" NumberOfElements="',i10,'">',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>')

7342 format(&
    '<Attribute Name="', A,'" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem DataType="Float" Precision="8" Dimensions="', i10, i10, '" Format="HDF">',/&
    '                    ', A, ':/', A, /&
    '</DataItem>',/,/&
    '</DataItem>',/&
    '</Attribute>')

7343 format(&
    '</Grid>',/)

736 format(&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

end program
!=========================================================================================
