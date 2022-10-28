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
!> Contains all the routines for NetCDF handling.
module nc_routines

#ifdef enable_netcdf
    use netcdf
#endif
    use nc_helpers
    use data_io, only: verbose, deflate_level, ncid_out
    use data_proc, only: mynum, nproc, lpr
    use global_parameters
    use commun, only: barrier, comm_elem_number

    implicit none
    save
    private

    !> Buffer variable for recorder
    real(sp), allocatable   :: recdumpvar(:,:,:)
    !> Buffer variable for displacement at surface
    real(sp), allocatable   :: surfdumpvar_disp(:,:,:)
    !> Buffer variable for velocity at surface
    real(sp), allocatable   :: surfdumpvar_velo(:,:,:)
    !> Buffer variable for strain at surface
    real(sp), allocatable   :: surfdumpvar_strain(:,:,:)
    !> Buffer variable for source displacement at surface
    real(sp), allocatable   :: surfdumpvar_srcdisp(:,:,:)

    !> Buffer variable for everything dumped in nc_dump_field_1d
    real(sp), allocatable   :: oneddumpvar(:,:,:)
    real(sp), allocatable   :: scoord1d(:), zcoord1d(:)
    real(sp), allocatable   :: scoord1d_mp(:), zcoord1d_mp(:)
    real(sp), allocatable   :: rho1d(:), mu1d(:), lambda1d(:)
    real(sp), allocatable   :: vp1d(:), vs1d(:)
    real(sp), allocatable   :: xi1d(:), phi1d(:), eta1d(:)
    real(sp), allocatable   :: Q_mu1d(:), Q_kappa1d(:)

    !> Number of steps before kernel specific stuff is dumped
    integer             :: dumpstepsnap
    !> Number of GLL points per element
    integer             :: gllperelem
    !> When is this processor supposed to dump.
    integer             :: outputplan
    !> How many steps since last dump?
    integer             :: stepstodump
    !> Global variables, so that we do not have to pass data to the C subroutine
    integer             :: isnap_global
    !> dito
    integer             :: ndumps
    !> Will any processor dump at this value of isnap?
    logical,allocatable :: dumpposition(:)
    !> Number of GLL points to plot for this processor
    integer             :: npoints
    !> Number of GLL points to plot for all processors
    integer             :: npoints_global
    !> Mapping of this processors GLL points to the global mesh
    integer             :: npoints_myfirst, npoints_mylast
    integer             :: nelem_myfirst, nelem_mylast
    !> Number of GLL points to plot in solid/fluid domain
    integer             :: npts_sol, npts_flu
    !> Number of GLL points to plot in solid domain for all processors
    integer             :: npts_sol_global
    !> Number of GLL points to plot in fluid domain for all processors
    integer             :: npts_flu_global
    !> Mapping of local solid points to global
    integer             :: npts_sol_myfirst, npts_sol_mylast
    !> Mapping of local fluid points to global
    integer             :: npts_flu_myfirst, npts_flu_mylast

    integer            :: ncid_recout, ncid_snapout, ncid_surfout, ncid_meshout
    integer            :: nc_snap_dimid, nc_proc_dimid, nc_rec_dimid, nc_recproc_dimid
    integer            :: nc_times_dimid, nc_comp_dimid, nc_disp_varid, nc_stf_seis_varid, nc_stf_d_seis_varid
    integer            :: nc_time_varid, nc_iter_dimid, nc_stf_iter_varid, nc_stf_d_iter_varid

    integer            :: nc_strcomp_dimid
    integer            :: nc_surfelem_disp_varid, nc_surfelem_velo_varid
    integer            :: nc_surfelem_strain_varid, nc_surfelem_disp_src_varid
    integer            :: nc_mesh_sol_varid, nc_mesh_flu_varid, nc_stf_dump_varid, nc_stf_d_dump_varid
    integer            :: nc_point_dimid, nc_pt_sol_dimid, nc_pt_flu_dimid
    integer            :: nc_szcoord_dimid
    integer            :: nc_snaptime_varid, nc_elem_dom_varid, nc_surfelem_theta_varid
    integer,allocatable :: nc_field_varid(:)
    character(len=16), allocatable  :: varnamelist(:)
    character(len=12), allocatable  :: nc_varnamelist(:)
    integer             :: nvar = -1

    !! Buffer variables to hand over to the dumping thread
    real(kind=sp), allocatable, dimension(:,:,:)  :: copy_oneddumpvar
    real(kind=sp), allocatable, dimension(:,:,:)  :: copy_surfdumpvar_disp
    real(kind=sp), allocatable, dimension(:,:,:)  :: copy_surfdumpvar_strain
    real(kind=sp), allocatable, dimension(:,:,:)  :: copy_surfdumpvar_velo
    real(kind=sp), allocatable, dimension(:,:,:)  :: copy_surfdumpvar_srcdisp


    !! Buffer variables for the STF.
    real(kind=sp), allocatable :: stf_dump_dumpvar(:)
    real(kind=sp), allocatable :: stf_seis_dumpvar(:)
    real(kind=sp), allocatable :: stf_dumpvar(:)

    !! Buffer variables for time derivative of the STF.
    real(kind=sp), allocatable :: stf_d_dump_dumpvar(:)
    real(kind=sp), allocatable :: stf_d_seis_dumpvar(:)
    real(kind=sp), allocatable :: stf_d_dumpvar(:)

    !> How many snaps should be buffered in RAM?
    integer             :: nc_dumpbuffersize

    !> chunking
    logical             :: nc_chunk_time_traces
    integer, parameter  :: disk_block_size = 8192

    public              :: nc_dump_strain, nc_dump_rec, nc_dump_surface
    public              :: nc_dump_field_solid, nc_dump_field_fluid
    public              :: nc_define_outputfile, nc_finish_prepare, nc_end_output
    public              :: nc_dump_strain_to_disk, nc_dump_mesh_sol, nc_dump_mesh_flu
    public              :: nc_dump_mesh_kwf
    public              :: nc_dump_mesh_mp_kwf
    public              :: nc_dump_elastic_parameters
    public              :: nc_dump_stf, nc_rec_checkpoint

    public              :: nc_dumpbuffersize, nc_chunk_time_traces
    public              :: set_npoints
contains

!-----------------------------------------------------------------------------------------
subroutine set_npoints(n)
  integer, intent(in) :: n
  npoints = n
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_data_xdmf(nc_filename_in, xdmf_filename_in, varname, npoints, nsnap)
  character(len=*), intent(in)      :: nc_filename_in, xdmf_filename_in, varname
  integer, intent(in)               :: npoints, nsnap

  integer                           :: iinput_xdmf
  integer                           :: i
  character(len=128)                :: xdmf_filename, nc_filename


  xdmf_filename = trim(nc_filename_in(:index(nc_filename_in, '/', back=.true.))) &
                    // xdmf_filename_in
  nc_filename = trim(nc_filename_in(index(nc_filename_in, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(xdmf_filename))
  write(iinput_xdmf, 733) npoints, npoints, trim(nc_filename), npoints, trim(nc_filename)

  do i=1, nsnap
     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) dble(i), npoints, "'", "'"

     ! write attribute
     write(iinput_xdmf, 7342) varname, npoints, i-1, npoints, nsnap, npoints, &
                              trim(nc_filename), trim(varname)

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

!-----------------------------------------------------------------------------------------
!> Routine to dump the wavefield variables for the Kerner. Collects input in
!! oneddumpvar_sol and oneddumpvar_flu until dumping condition is fulfilled.
subroutine nc_dump_field_solid(f, varname)

    real(kind=realkind), intent(in)   :: f(:)     ! < Data to dump, size should be npts_sol
    character(len=*), intent(in)      :: varname  ! < Internal name of data to dump.
    !! Is used to identify the NetCDF Variable in question
#ifdef enable_netcdf
    integer                           :: ivar

    do ivar=1, nvar
        ! < Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar)) == trim(varname)) exit
    enddo

    if (ivar > nvar/2) then
        write(*,*) 'nc_dump_field_solid: Trying to access variable: ', trim(varname), &
            ' which is a fluid variable. Contact a developer and shout at him!'
        stop 1
    endif

    oneddumpvar(1:npts_sol,stepstodump+1,ivar) = f  !processor specific dump variable
#endif
end subroutine nc_dump_field_solid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Routine to dump the wavefield variables for the Kerner. Collects input in
!! oneddumpvar_sol and oneddumpvar_flu until dumping condition is fulfilled.
subroutine nc_dump_field_fluid(f, varname)

    real(kind=realkind), intent(in)   :: f(:)  ! < Data to dump, size should be npts_flu
    character(len=*), intent(in)      :: varname  ! < Internal name of data to dump.
    !! Is used to identify the NetCDF Variable in question
#ifdef enable_netcdf
    integer                           :: ivar

    do ivar=1, nvar
        ! < Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar)) == trim(varname)) exit
    enddo

    if (ivar <= nvar/2) then
        write(*,*) 'nc_dump_field_fluid: Trying to access variable: ', trim(varname), &
            ' which is a solid variable. Contact a developer and shout at him!'
        stop 1
    endif

    oneddumpvar(npts_sol+1:npoints,stepstodump+1,ivar-nvar/2) = f
#endif
end subroutine nc_dump_field_fluid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain(isnap_loc)

    use data_io, only: nstrain
    use clocks_mod, only: tick
    use data_time, only: iclocknbio, idnbio
    use data_mesh, only: maxind

    ! explicit interfaces to the c functions to avoide the underscore issues
    ! (Fortran 2003 standard)
    interface
        subroutine c_spawn_dumpthread(stepstodump) bind(c, name='c_spawn_dumpthread')
            use, intrinsic      :: iso_c_binding, only: c_int
            integer (c_int)     :: stepstodump
        end subroutine
    end interface

    interface
        subroutine c_wait_for_io() bind(c, name='c_wait_for_io')
        end subroutine
    end interface

    integer, intent(in) :: isnap_loc

#ifdef enable_netcdf
    integer             :: iproc
    real                :: tickl, tackl

    stepstodump = stepstodump + 1
    if (isnap_loc == 0) return

    if (dumpposition(mod(isnap_loc, dumpstepsnap))) then

#ifndef enable_parallel_netcdf

        ! wait for other processes to finish writing, measure waiting time and
        ! issue warning in case waiting longer then .5 sec
        ! MvD: I am not sure if cpu_time is a correct measure here, as we idle
        !      until IO is finished.
        !      Therefor testing system_clock, from the clocks module.
        if (mod(isnap_loc, dumpstepsnap) == outputplan .and. (npoints > 0 .or. maxind > 0)) then
            if (verbose > 1) then
                write(*,"('   Proc ', I4, ': Would like to dump data and waits for his turn')") mynum
                call flush(6)
            endif
            call cpu_time(tickl)
        endif

        iclocknbio = tick()
        do iproc=0, nproc-1
            if (iproc == mynum .and. (npoints > 0 .or. maxind > 0)) then
                call c_wait_for_io()
                call flush(6)
            endif
            call barrier
        enddo
        iclocknbio = tick(id=idnbio, since=iclocknbio)

        ! non blocking write
        if (mod(isnap_loc, dumpstepsnap) == outputplan .and. (npoints > 0 .or. maxind > 0)) then
            call cpu_time(tackl)
            if ((tackl-tickl) > 0.5 .and. verbose > 0) then
                write(*,"('WARNING: Computation was halted for ', F7.2, ' s to wait for ', &
                         & 'dumping processor. Consider adapting netCDF output variables', &
                         & '(disable compression, increase dumpstepsnap)')") tackl-tickl
            endif
#endif

            isnap_global = isnap_loc
            ndumps = stepstodump

            allocate(copy_oneddumpvar(1:npoints,1:ndumps,1:nvar/2))
            allocate(copy_surfdumpvar_disp(1:ndumps, 1:3, 1:maxind))
            allocate(copy_surfdumpvar_strain(1:ndumps, 1:6, 1:maxind))
            allocate(copy_surfdumpvar_velo(1:ndumps, 1:3, 1:maxind))
            allocate(copy_surfdumpvar_srcdisp(1:ndumps, 1:3, 1:maxind))

            copy_oneddumpvar          = oneddumpvar(1:npoints,1:ndumps,1:nvar/2)
            copy_surfdumpvar_disp     = surfdumpvar_disp(1:ndumps, 1:3, 1:maxind)
            copy_surfdumpvar_strain   = surfdumpvar_strain(1:ndumps, 1:6, 1:maxind)
            copy_surfdumpvar_velo     = surfdumpvar_velo(1:ndumps, 1:3, 1:maxind)
            copy_surfdumpvar_srcdisp  = surfdumpvar_srcdisp(1:ndumps, 1:3, 1:maxind)

#ifndef enable_parallel_netcdf

            call c_spawn_dumpthread(stepstodump)
            stepstodump = 0
        endif
#else
            call nc_dump_strain_to_disk
            stepstodump = 0
#endif
    endif

    ! Final and last dump of all remaining
    ! not done in unblocking fashion using a thread, as there is nothing to
    ! compute anymore
    ! Make sure, nobody is accessing the output file anymore
    if (isnap_loc == nstrain) then

#ifndef enable_parallel_netcdf
        do iproc=0, nproc-1
            if (iproc == mynum) then
                call c_wait_for_io()
                call barrier
            endif
        enddo
        do iproc=0,nproc-1
            if (iproc == mynum .and. (npoints > 0 .or. maxind > 0)) then
#endif
                isnap_global = nstrain
                ndumps = stepstodump

                allocate(copy_oneddumpvar(1:npoints,1:ndumps,1:nvar/2))
                allocate(copy_surfdumpvar_disp(1:ndumps, 1:3, 1:maxind))
                allocate(copy_surfdumpvar_strain(1:ndumps, 1:6, 1:maxind))
                allocate(copy_surfdumpvar_velo(1:ndumps, 1:3, 1:maxind))
                allocate(copy_surfdumpvar_srcdisp(1:ndumps, 1:3, 1:maxind))

                copy_oneddumpvar          = oneddumpvar(1:npoints,1:ndumps,1:nvar/2)
                copy_surfdumpvar_disp     = surfdumpvar_disp(1:ndumps, 1:3, 1:maxind)
                copy_surfdumpvar_strain   = surfdumpvar_strain(1:ndumps, 1:6, 1:maxind)
                copy_surfdumpvar_velo     = surfdumpvar_velo(1:ndumps, 1:3, 1:maxind)
                copy_surfdumpvar_srcdisp  = surfdumpvar_srcdisp(1:ndumps, 1:3, 1:maxind)

                call nc_dump_strain_to_disk()
                if (verbose > 1) write(*,*) mynum, 'finished dumping strain'
#ifndef enable_parallel_netcdf
            endif
            call barrier
        enddo
#endif
    endif

#endif
end subroutine nc_dump_strain
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain_to_disk() bind(c, name="nc_dump_strain_to_disk")
#ifdef enable_netcdf

    use data_io
    use data_source, only: src_type
    use global_parameters, only: realkind
    use data_mesh, only: loc2globrec, maxind, ind_first

    integer                    :: ivar, isnap_loc
    real                       :: tick, tack, dumpsize_MB
    integer                    :: dumpsize

    dumpsize = 0
    call cpu_time(tick)

#ifndef enable_parallel_netcdf
    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                          mode=NF90_WRITE, ncid=ncid_out) )
#else
    call barrier
#endif

    isnap_loc = isnap_global
    if (verbose > 1) then
        if (ndumps == 0) then
            write(*,"('   Proc ', I4, ': in dump routine, isnap =', I5, &
                    & ', nothing to dump, returning...')") mynum, isnap_loc
            return
        else
            write(*,"('   Proc ', I4, ': in dump routine, isnap =', I5, &
                    & ', stepstodump = ', I4)") mynum, isnap_loc, ndumps
        endif
    endif

    call getgrpid(ncid_out, "Snapshots", ncid_snapout)
    do ivar=1, nvar/2
        call putvar_real2d(ncid   = ncid_snapout, &
                           varid  = nc_field_varid(ivar), &
                           start  = [npoints_myfirst, isnap_loc-ndumps+1], &
                           count  = [npoints, ndumps], &
                           values = copy_oneddumpvar(1:npoints,1:ndumps,ivar))
        dumpsize = dumpsize + npoints * ndumps
    enddo

    !> Surface dumps
    call putvar_real3d(ncid_surfout, nc_surfelem_disp_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 3, maxind], &
                values = copy_surfdumpvar_disp(1:ndumps, 1:3, 1:maxind))
    dumpsize = dumpsize + 3 * maxind * ndumps

    call putvar_real3d(ncid_surfout, nc_surfelem_velo_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 3, maxind], &
                values = copy_surfdumpvar_velo(1:ndumps, 1:3, 1:maxind))
    dumpsize = dumpsize + 3 * maxind * ndumps

    call putvar_real3d(ncid_surfout, nc_surfelem_strain_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 6, maxind], &
                values = copy_surfdumpvar_strain(1:ndumps, 1:6, 1:maxind))
    dumpsize = dumpsize + 6 * maxind * ndumps

    call putvar_real3d(ncid_surfout, nc_surfelem_disp_src_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 3, maxind], &
                values = copy_surfdumpvar_srcdisp(1:ndumps, 1:3, 1:maxind))
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_put_att(ncid_out, NF90_GLOBAL, 'percent completed', &
                             isnap_loc*100/nstrain) )

#ifndef enable_parallel_netcdf
    call check( nf90_close(ncid_out) )
#endif

    call cpu_time(tack)

    deallocate(copy_oneddumpvar)
    deallocate(copy_surfdumpvar_disp)
    deallocate(copy_surfdumpvar_strain)
    deallocate(copy_surfdumpvar_velo)
    deallocate(copy_surfdumpvar_srcdisp)

    if (verbose > 1) then
        dumpsize_MB = real(dumpsize) * 4 / 1048576
        write(*,"('   Proc', I5,': Wrote ', F8.2, ' MB in ', F6.2, 's (', F8.2, ' MB/s)')") &
            mynum, dumpsize_MB, tack-tick, dumpsize_MB / (tack - tick)
        call flush(6)
    endif

#endif
end subroutine nc_dump_strain_to_disk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_stf(stf)
    use data_io, only: nseismo, nstrain, dump_wavefields
    use data_time, only: seis_it, strain_it, niter, deltat
    real(kind=sp), intent(in), dimension(:) :: stf

#ifdef enable_netcdf
    integer                                 :: it_s, it_d, i

    allocate(stf_dumpvar(niter))
    allocate(stf_seis_dumpvar(nseismo))
    allocate(stf_dump_dumpvar(nstrain))

    allocate(stf_d_dumpvar(niter))
    allocate(stf_d_seis_dumpvar(nseismo))
    allocate(stf_d_dump_dumpvar(nstrain))

    stf_d_dumpvar = 0

    stf_seis_dumpvar = 0
    stf_dump_dumpvar = 0

    stf_d_seis_dumpvar = 0
    stf_d_dump_dumpvar = 0

    ! compute derivative using central differences
    do i = 2, niter-1
       stf_d_dumpvar(i) = (stf(i+1) - stf(i-1)) / (2 * deltat)
    enddo
    ! compute derivative using fwd/bwd difference at the first and last sample
    stf_d_dumpvar(1) = (stf(2) - stf(1)) / deltat
    stf_d_dumpvar(niter) = (stf(niter) - stf(niter-1)) / deltat

    it_s = 1
    it_d = 1

    do i = 1, niter
        ! Dumping the STF in the fine time stepping of the seismogram output
        if ( mod(i,seis_it) == 0) then
           it_s = it_s + 1
           stf_seis_dumpvar(it_s) = stf(i)
           stf_d_seis_dumpvar(it_s) = stf_d_dumpvar(i)
        endif

        if (dump_wavefields) then
            ! Dumping the STF in the coarse time stepping of the strain (KERNER) output
            if ( mod(i,strain_it) == 0) then
               it_d = it_d + 1
               stf_dump_dumpvar(it_d) = stf(i)
               stf_d_dump_dumpvar(it_d) = stf_d_dumpvar(i)
            endif
        endif
    enddo
    stf_dumpvar = stf

#endif
end subroutine nc_dump_stf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Dump receiver specific stuff, especially displacement and velocity
!! N.B.: Works with global indices.
subroutine nc_dump_rec(recfield, iseismo)
    use data_mesh, only: num_rec
    real(sp), intent(in), dimension(3,num_rec) :: recfield
    integer, intent(in)                        :: iseismo

#ifdef enable_netcdf
    recdumpvar(iseismo,:,:) = recfield(:,:)
#endif
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_rec_to_disk
#ifdef enable_netcdf
    use data_mesh, only: loc2globrec, num_rec
    use data_io, only: datapath, lfdata, nseismo

    real              :: tick, tack
    integer           :: irec, dumpsize, icomp

    call cpu_time(tick)

#ifndef enable_parallel_netcdf
    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                          mode=NF90_WRITE, ncid=ncid_out) )
#endif

    call getgrpid( ncid_out, "Seismograms", ncid_recout)
    call getvarid( ncid_recout, "displacement", nc_disp_varid )

    do irec = 1, num_rec
        do icomp = 1, 3
            call putvar_real3d(ncid   = ncid_recout, &
                               varid  = nc_disp_varid, &
                               start  = [1, icomp, loc2globrec(irec)], &
                               count  = [nseismo, 1, 1], &
                               values = reshape(recdumpvar(:,icomp,irec), [nseismo, 1, 1]) )
        enddo
    enddo

    dumpsize = nseismo * 3 * num_rec

#ifndef enable_parallel_netcdf
    call check( nf90_close(ncid=ncid_out))
#endif

    call cpu_time(tack)

    if (verbose > 2) then
        write(*,"(I3,': Receiver data, Wrote ', F8.3, ' MB in ', F6.2, 's')") &
            mynum, real(dumpsize) * 4. / 1048576., tack-tick
    endif
#endif
end subroutine nc_dump_rec_to_disk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_rec_checkpoint
    use data_mesh, only: loc2globrec, num_rec

#ifdef enable_netcdf
#ifndef enable_parallel_netcdf

    interface
        subroutine c_wait_for_io() bind(c, name='c_wait_for_io')
        end subroutine
    end interface

    integer         :: iproc

    call c_wait_for_io()
    do iproc=0, nproc-1
        call barrier
        if (iproc == mynum) then
#endif
            if (num_rec > 0) then
                if (verbose > 2) &
                    write(*,"('   Proc ', I3, ' will dump receiver seismograms')") mynum
                call nc_dump_rec_to_disk()
                call flush(6)
            else
                if (verbose > 2) &
                    write(*,"('   Proc ', I3, ' has no receivers and just waits for the others')") mynum
            endif
#ifndef enable_parallel_netcdf
        endif

    enddo
#endif
    call barrier
#endif
end subroutine nc_rec_checkpoint
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
!> Dump stuff along surface
subroutine nc_dump_surface(surffield, disporvelo)

    real(kind=realkind), intent(in), dimension(:,:) :: surffield
    character(len=4), intent(in)                    :: disporvelo
#ifdef enable_netcdf

    select case(disporvelo)
        case('disp')
            surfdumpvar_disp(stepstodump+1,:,:) = transpose(surffield(:,:))
        case('velo')
            surfdumpvar_velo(stepstodump+1,:,:) = transpose(surffield(:,:))
        case('stra')
            surfdumpvar_strain(stepstodump+1,:,:) = transpose(surffield(:,:))
        case('srcd')
            surfdumpvar_srcdisp(stepstodump+1,:,:) = transpose(surffield(:,:))
    end select

#endif
end subroutine nc_dump_surface
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_sol(scoord_sol, zcoord_sol)

    use data_io, only: ndumppts_el
    use data_mesh, only: nelem
    real(sp), intent(in) :: scoord_sol(:,:,:)
    real(sp), intent(in) :: zcoord_sol(size(scoord_sol,1), size(scoord_sol,2), &
                                       size(scoord_sol,3))
#ifdef enable_netcdf

    npts_sol = size(scoord_sol)

    npoints = ndumppts_el * nelem
    allocate(scoord1d(npoints))
    allocate(zcoord1d(npoints))

    zcoord1d(1:npts_sol) = pack(zcoord_sol,.true.)
    scoord1d(1:npts_sol) = pack(scoord_sol,.true.)


#endif
end subroutine nc_dump_mesh_sol
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_flu(scoord_flu, zcoord_flu)
! can only be called after calling nc_dump_mesh_sol

    real(sp), intent(in) :: scoord_flu(:,:,:)
    real(sp), intent(in) :: zcoord_flu(size(scoord_flu,1), size(scoord_flu,2), &
                                       size(scoord_flu,3))
#ifdef enable_netcdf

    npts_flu = size(scoord_flu)

    zcoord1d(npts_sol+1:) = pack(zcoord_flu, .true.)
    scoord1d(npts_sol+1:) = pack(scoord_flu, .true.)

#endif
end subroutine nc_dump_mesh_flu
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_kwf(coords, nsol, nflu)

    real(sp), intent(in) :: coords(:,:)
    integer, intent(in)  :: nsol, nflu
#ifdef enable_netcdf

    npts_sol = nsol
    npts_flu = nflu

    npoints = nsol + nflu

    if (size(coords, 1) /= npoints) then
       write(*,*) 'ERROR: inconsistent point numbers'
       call abort()
    endif

    allocate(scoord1d(npoints))
    allocate(zcoord1d(npoints))

    scoord1d(:) = coords(:,1)
    zcoord1d(:) = coords(:,2)

#endif
end subroutine nc_dump_mesh_kwf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_mp_kwf(coords, nel)

    real(sp), intent(in) :: coords(:,:)
    integer, intent(in)  :: nel
#ifdef enable_netcdf

    if (size(coords, 1) /= nel) then
       write(*,*) 'ERROR: inconsistent element numbers'
       call abort()
    endif

    allocate(scoord1d_mp(nel))
    allocate(zcoord1d_mp(nel))

    scoord1d_mp(:) = coords(:,1)
    zcoord1d_mp(:) = coords(:,2)

#endif
end subroutine nc_dump_mesh_mp_kwf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_elastic_parameters(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                      fa_ani_theta, fa_ani_phi, Q_mu, Q_kappa)

    use data_io, only: ibeg, iend, jbeg, jend, dump_type
    use data_mesh, only: mapping_ijel_ikwf, ielsolid, ielfluid, nel_solid, nel_fluid, &
                            npol, kwf_mask

    real(kind=dp), dimension(0:,0:,:), intent(in)       :: rho, lambda, mu, xi_ani
    real(kind=dp), dimension(0:,0:,:), intent(in)       :: phi_ani, eta_ani
    real(kind=dp), dimension(0:,0:,:), intent(in)       :: fa_ani_theta, fa_ani_phi
    real(kind=dp), dimension(0:,0:,:), intent(in), optional   :: Q_mu, Q_kappa

    integer :: size1d
    integer :: iel, ipol, jpol, ct

    if (dump_type == 'displ_only' .or. dump_type == 'strain_only') then
       allocate(rho1d(npoints))
       allocate(lambda1d(npoints))
       allocate(mu1d(npoints))
       allocate(vp1d(npoints))
       allocate(vs1d(npoints))
       allocate(xi1d(npoints))
       allocate(phi1d(npoints))
       allocate(eta1d(npoints))
       if (present(Q_mu) .and. present(Q_kappa)) then
           allocate(Q_mu1d(npoints))
           allocate(Q_kappa1d(npoints))
       endif

       do iel=1, nel_solid
           do ipol=0, npol
               do jpol=0, npol
                   if (kwf_mask(ipol,jpol,iel)) then
                       ct = mapping_ijel_ikwf(ipol,jpol,iel)

                       rho1d(ct) = rho(ipol,jpol,ielsolid(iel))
                       lambda1d(ct) = lambda(ipol,jpol,ielsolid(iel))
                       mu1d(ct) = mu(ipol,jpol,ielsolid(iel))
                       xi1d(ct) = xi_ani(ipol,jpol,ielsolid(iel))
                       phi1d(ct) = phi_ani(ipol,jpol,ielsolid(iel))
                       eta1d(ct) = eta_ani(ipol,jpol,ielsolid(iel))

                       if (present(Q_mu) .and. present(Q_kappa)) then
                           Q_mu1d(ct) = Q_mu(ipol,jpol,ielsolid(iel))
                           Q_kappa1d(ct) = Q_kappa(ipol,jpol,ielsolid(iel))
                       endif
                   endif
               enddo
           enddo
       enddo

       do iel=1, nel_fluid
           do ipol=0, npol
               do jpol=0, npol
                   if (kwf_mask(ipol,jpol,iel + nel_solid)) then
                       ct = mapping_ijel_ikwf(ipol,jpol,iel + nel_solid)

                       rho1d(ct) = rho(ipol,jpol,ielfluid(iel))
                       lambda1d(ct) = lambda(ipol,jpol,ielfluid(iel))
                       mu1d(ct) = mu(ipol,jpol,ielfluid(iel))
                       xi1d(ct) = xi_ani(ipol,jpol,ielfluid(iel))
                       phi1d(ct) = phi_ani(ipol,jpol,ielfluid(iel))
                       eta1d(ct) = eta_ani(ipol,jpol,ielfluid(iel))

                       if (present(Q_mu) .and. present(Q_kappa)) then
                           Q_mu1d(ct) = Q_mu(ipol,jpol,ielfluid(iel))
                           Q_kappa1d(ct) = Q_kappa(ipol,jpol,ielfluid(iel))
                       endif
                   endif
               enddo
           enddo
       enddo

       vp1d = sqrt((lambda1d + 2 * mu1d ) / rho1d)
       vs1d = sqrt(mu1d  / rho1d)

    else
       size1d = size(rho(ibeg:iend, jbeg:jend, :))
       print *, ' NetCDF: Mesh elastic parameter variables have size:', size1d
       allocate(rho1d(size1d))
       allocate(lambda1d(size1d))
       allocate(mu1d(size1d))
       allocate(vp1d(size1d))
       allocate(vs1d(size1d))
       allocate(xi1d(size1d))
       allocate(phi1d(size1d))
       allocate(eta1d(size1d))

       rho1d = real(pack(rho(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
       lambda1d = real(pack(lambda(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
       mu1d = real(pack(mu(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
       vp1d = sqrt((lambda1d + 2 * mu1d ) / rho1d)
       vs1d = sqrt(mu1d / rho1d)

       xi1d = real(pack(xi_ani(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
       phi1d = real(pack(phi_ani(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
       eta1d = real(pack(eta_ani(ibeg:iend, jbeg:jend, :), .true.), kind=sp)

       if (present(Q_mu) .and. present(Q_kappa)) then
           allocate(Q_mu1d(size1d))
           allocate(Q_kappa1d(size1d))
           Q_mu1d = real(pack(Q_mu(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
           Q_kappa1d = real(pack(Q_kappa(ibeg:iend, jbeg:jend, :), .true.), kind=sp)
       endif
    endif

end subroutine nc_dump_elastic_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Define the output file variables and dimensions
!! and allocate buffer variables.
subroutine nc_define_outputfile(nrec, rec_names, rec_th, rec_th_req, rec_ph, rec_proc)

    use data_io, only: nseismo, nstrain, nseismo, ibeg, iend, jbeg, jend, &
                           dump_wavefields, dump_type
    use data_io, only: datapath, lfdata, strain_samp
    use data_mesh, only: maxind, num_rec, discont, nelem, nel_solid, nel_fluid, &
                           ndisc, maxind_glob, nelem_kwf_global, npoint_kwf, npoint_solid_kwf, &
                           npoint_fluid_kwf, npol, nelem_kwf, npoint_kwf_global, anel_true
    use data_source, only: src_type, t_0
    use data_time, only: deltat, niter

! using MPI here is fine, because enable_parallel_netcdf can only be set if
! compiled with mpi
#ifdef enable_parallel_netcdf
    use mpi
#endif

    integer, intent(in)                  :: nrec              ! < Number of receivers
    character(len=40),intent(in)         :: rec_names(nrec)   ! < Receiver names
    real(dp), dimension(nrec),intent(in) :: rec_th            ! < Receiver theta
    real(dp), dimension(nrec),intent(in) :: rec_th_req        ! < Requested receiver theta
    real(dp), dimension(nrec),intent(in) :: rec_ph            ! < Receiver phi
    integer, dimension(nrec),intent(in)  :: rec_proc          ! < Receiver processor
#ifdef enable_netcdf
    real(dp), dimension(nstrain)         :: time_strain
    real(dp), dimension(nseismo)         :: time_seis
    character(len=16), allocatable       :: varname(:)
    character(len=256)                   :: nc_fnam
    integer                              :: ivar, i
    integer                              :: irec, iproc, nmode
    integer                              :: nc_ph_varid
    integer                              :: nc_thr_varid, nc_th_varid
    integer                              :: nc_proc_varid, nc_recnam_dimid
    integer                              :: nc_recnam_varid, nc_surf_dimid
    integer                              :: nc_pt_dimid
    integer                              :: nc_mesh_s_varid, nc_mesh_z_varid
    integer                              :: nc_mesh_s_mp_varid, nc_mesh_z_mp_varid
    integer                              :: nc_mesh_vs_varid, nc_mesh_vp_varid
    integer                              :: nc_mesh_mu_varid, nc_mesh_rho_varid
    integer                              :: nc_mesh_lambda_varid
    integer                              :: nc_mesh_xi_varid
    integer                              :: nc_mesh_phi_varid
    integer                              :: nc_mesh_eta_varid
    integer                              :: nc_mesh_qmu_varid
    integer                              :: nc_mesh_qka_varid
    integer                              :: nc_mesh_midpoint_varid
    integer                              :: nc_mesh_fem_varid
    integer                              :: nc_mesh_eltype_varid
    integer                              :: nc_mesh_axis_varid
    integer                              :: nc_mesh_sem_varid
    integer                              :: nc_mesh_G0_varid
    integer                              :: nc_mesh_G1_varid
    integer                              :: nc_mesh_G2_varid
    integer                              :: nc_mesh_gll_varid
    integer                              :: nc_mesh_glj_varid
    integer                              :: nc_mesh_elem_dimid, nc_mesh_npol_dimid
    integer                              :: nc_mesh_cntrlpts_dimid

    integer                              :: chunk_pt, chunk_time
                                    ! < Contains chunk size in GLL points. Should be system
                                    !! int(disk_block_size / nsnap)

    if ((mynum == 0) .and. (verbose > 1)) then
        write(*,*)
        write(*,*) '************************************************************************'
        write(*,*) '**** Producing netcdf output file with its variables and dimensions ****'
        write(*,*) '************************************************************************'
        write(*,*)
    endif

    nc_fnam = datapath(1:lfdata)//"/axisem_output.nc4"

    call barrier

#ifndef enable_parallel_netcdf
    if (nc_dumpbuffersize < nproc) then
        write(*,*) 'ERROR: NETCDF_DUMP_BUFFER needs to be larger then the number of processors'
        stop
    endif

    dumpstepsnap = int(nc_dumpbuffersize / nproc) * nproc
        ! Will later be reduced to nstrain, if this is smaller than value given here

    if (lpr .and. verbose > 1) write(*,*) '  Dumping NetCDF file to disk every', dumpstepsnap, ' snaps'

    call barrier ! for nicer output only

    outputplan = mynum * (dumpstepsnap / nproc)

    allocate(dumpposition(0:dumpstepsnap-1))
    dumpposition = .false.
    do iproc=0, nproc-1
        dumpposition(iproc*(dumpstepsnap/nproc)) = .true.
        if ((iproc == mynum) .and. (verbose > 1)) then
            write(*,"(' Proc ', I4, ' will dump at position ', I4)") mynum, outputplan
            call flush(6)
        endif
        call barrier ! for nicer output only
    enddo
    if (nstrain <= dumpstepsnap) dumpstepsnap = nstrain
#else
    ! in parallel IO, always everybody dumps collectively
    dumpstepsnap = min(nstrain, nc_dumpbuffersize)
    allocate(dumpposition(0:dumpstepsnap-1))
    dumpposition(:) = .false.
    dumpposition(0) = .true.
    outputplan = 0
    if ((mynum == 0) .and. (verbose > 1)) &
            write(*,"(' all Procs will dump every ', I4, ' steps')") dumpstepsnap
#endif


    if (dump_wavefields) then
        select case (trim(dump_type))
           case ('displ_only')
              if (src_type(1) == 'monopole') then
                  nvar = 4
              else
                  nvar = 6
              endif
              allocate(varname(nvar))
              allocate(varnamelist(nvar))
              allocate(nc_varnamelist(nvar/2))
              allocate(nc_field_varid(nvar/2))

              if (src_type(1) == 'monopole') then
                  varnamelist =    ['disp_sol_s     ', 'disp_sol_z     ', &
                                    'disp_flu_s     ', 'disp_flu_z     ']

                  nc_varnamelist = ['disp_s     ', 'disp_z     ']
              else
                  varnamelist =    ['disp_sol_s     ', 'disp_sol_p     ', 'disp_sol_z     ', &
                                    'disp_flu_s     ', 'disp_flu_p     ', 'disp_flu_z     ']

                  nc_varnamelist = ['disp_s     ', 'disp_p     ', 'disp_z     ']
              endif

              npoints = npoint_kwf

              npts_sol = npoint_solid_kwf
              npts_flu = npoint_fluid_kwf

              call comm_elem_number(nelem_kwf, nelem_kwf_global, nelem_myfirst, nelem_mylast)

           case ('strain_only')
              if (src_type(1) == 'monopole') then
                  nvar = 8
              else
                  nvar = 12
              endif
              allocate(varname(nvar))
              allocate(varnamelist(nvar))
              allocate(nc_varnamelist(nvar/2))
              allocate(nc_field_varid(nvar/2))

              if (src_type(1) == 'monopole') then
                  varnamelist =    ['strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                                    'straintrace_sol', &
                                    'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                                    'straintrace_flu']

                  nc_varnamelist = ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
                                    'straintrace']
              else
                  varnamelist =    ['strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                                    'strain_dsup_sol', 'strain_dzup_sol', 'straintrace_sol', &
                                    'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                                    'strain_dsup_flu', 'strain_dzup_flu', 'straintrace_flu']

                  nc_varnamelist = ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
                                    'strain_dsup', 'strain_dzup', 'straintrace']
              endif

              npoints = npoint_kwf


              npts_sol = npoint_solid_kwf
              npts_flu = npoint_fluid_kwf


           case ('displ_velo')
              write(*,*) 'ERROR: not yet implemented with netcdf'
              stop 2

           case ('fullfields')
              if (src_type(1) == 'monopole') then
                  nvar = 12
              else
                  nvar = 18
              endif
              allocate(varname(nvar))
              allocate(varnamelist(nvar))
              allocate(nc_varnamelist(nvar/2))
              allocate(nc_field_varid(nvar/2))

              if (src_type(1) == 'monopole') then
                  varnamelist =    ['strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                                    'straintrace_sol', 'velo_sol_s     ', 'velo_sol_z     ', &
                                    'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                                    'straintrace_flu', 'velo_flu_s     ', 'velo_flu_z     ']

                  nc_varnamelist = ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
                                    'straintrace', 'velo_s     ', 'velo_z     ']
              else
                  varnamelist =    ['strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                                    'strain_dsup_sol', 'strain_dzup_sol', 'straintrace_sol', &
                                    'velo_sol_s     ', 'velo_sol_p     ', 'velo_sol_z     ', &
                                    'strain_dsus_flu', 'strain_dsuz_flu', 'strain_dpup_flu', &
                                    'strain_dsup_flu', 'strain_dzup_flu', 'straintrace_flu', &
                                    'velo_flu_s     ', 'velo_flu_p     ', 'velo_flu_z     ']

                  nc_varnamelist = ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
                                    'strain_dsup', 'strain_dzup', 'straintrace', &
                                    'velo_s     ', 'velo_p     ', 'velo_z     ']
              endif

              gllperelem = (iend - ibeg + 1) * (jend - jbeg + 1)
              npoints = nelem * gllperelem

              npts_sol = nel_solid * gllperelem
              npts_flu = nel_fluid * gllperelem

        end select

        call comm_elem_number(npoints, npoints_global, npoints_myfirst, npoints_mylast)
        npoint_kwf_global = npoints_global
        call comm_elem_number(npts_sol, npts_sol_global, npts_sol_myfirst, npts_sol_mylast)
        call comm_elem_number(npts_flu, npts_flu_global, npts_flu_myfirst, npts_flu_mylast)

        if (lpr) then
          do ivar=1, nvar/2 ! The big snapshot variables for the kerner.
              call dump_mesh_data_xdmf(trim(nc_fnam), trim(nc_varnamelist(ivar))//'.xdmf', &
                                       'Snapshots/'//trim(nc_varnamelist(ivar)), &
                                       npts_sol_global + npts_flu_global,  nstrain)
          enddo
        endif

    endif ! dump_wavefields

    if (mynum == 0) then
        if (verbose > 1) write(*,*) ' Preparing netcdf file for ', nproc, ' processors'
        nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
        call check( nf90_create(path=nc_fnam, cmode=nmode, ncid=ncid_out) )
        if (verbose > 1) write(*,*) ' Netcdf file with ID ', ncid_out, ' produced.'


        if (verbose > 1) write(*,*) '  Producing groups for Seismograms and Snapshots'

        call check( nf90_def_grp(ncid_out, "Seismograms", ncid_recout) )
        call check( nf90_def_grp(ncid_out, "Snapshots", ncid_snapout) )
        call check( nf90_def_grp(ncid_out, "Surface", ncid_surfout) )
        call check( nf90_def_grp(ncid_out, "Mesh", ncid_meshout) )

        if (verbose > 1) write(*,*) 'Define dimensions in ''Seismograms'' group of NetCDF output file'

        call check( nf90_def_dim(ncid_out, "seis_timesteps", nseismo, nc_times_dimid) )
        call check( nf90_def_dim(ncid_out, "sim_timesteps", niter, nc_iter_dimid) )
        call check( nf90_def_dim(ncid_recout, "receivers", nrec, nc_rec_dimid) )
        call check( nf90_def_dim(ncid_out, "components", 3, nc_comp_dimid) )
        call check( nf90_def_dim(ncid_recout, "recnamlength", 40, nc_recnam_dimid) )

        if (verbose > 1) write(*,*) 'Define variables in ''Seismograms'' group of NetCDF output file'
        call flush(6)

        call check( nf90_def_var(ncid=ncid_recout, name="displacement", xtype=NF90_FLOAT, &
                                 dimids=[nc_times_dimid, nc_comp_dimid, nc_rec_dimid], &
                                 contiguous = .false., chunksizes=[nseismo,3,1], &
                                 varid=nc_disp_varid) )

        call check( nf90_put_att(ncid_recout, nc_disp_varid, 'units', 'meters') )
        call check( nf90_put_att(ncid_recout, nc_disp_varid, '_FillValue', 0.0) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_seis", xtype=NF90_FLOAT, &
                                 dimids=[nc_times_dimid], varid=nc_stf_seis_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_d_seis", xtype=NF90_FLOAT, &
                                 dimids=[nc_times_dimid], varid=nc_stf_d_seis_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_iter", xtype=NF90_FLOAT, &
                                 dimids=[nc_iter_dimid], varid=nc_stf_iter_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_d_iter", xtype=NF90_FLOAT, &
                                 dimids=[nc_iter_dimid], varid=nc_stf_d_iter_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="time", xtype=NF90_DOUBLE, &
                                 dimids=[nc_times_dimid], varid=nc_time_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="phi", xtype=NF90_FLOAT, &
                                 dimids=[nc_rec_dimid], varid=nc_ph_varid) )
        call check( nf90_def_var(ncid_recout, "theta_requested", NF90_FLOAT, &
                                 [nc_rec_dimid], nc_thr_varid) )
        call check( nf90_def_var(ncid_recout, "theta", NF90_FLOAT, &
                                 [nc_rec_dimid], nc_th_varid) )
        call check( nf90_def_var(ncid_recout, "processor_of_receiver", NF90_INT, &
                                 [nc_rec_dimid], nc_proc_varid) )
        call check( nf90_def_var(ncid_recout, "receiver_name", NF90_CHAR, &
                                 [nc_rec_dimid, nc_recnam_dimid], nc_recnam_varid) )

        if (dump_wavefields) then
            ! Wavefields group of output file N.B: Snapshots for kernel calculation
            if (verbose > 1) write(*,*) 'Define variables in ''Snapshots'' group of NetCDF output file', &
                                        '  awaiting', nstrain, ' snapshots'

            call check( nf90_def_dim( ncid   = ncid_out, &
                                      name   = 'snapshots', &
                                      len    = nstrain, &
                                      dimid  = nc_snap_dimid) )
            call check( nf90_put_att( ncid   = ncid_snapout, &
                                      varid  = NF90_GLOBAL, &
                                      name   = 'nstrain', &
                                      values = nstrain) )

            call check( nf90_def_dim( ncid   = ncid_out, &
                                      name   = 'gllpoints_all', &
                                      len    = npoints_global, &
                                      dimid  = nc_pt_dimid) )
            call check( nf90_put_att( ncid   = ncid_out, &
                                      varid  = NF90_GLOBAL, &
                                      name   = 'npoints', &
                                      values = npoints_global) )

            call check( nf90_def_var( ncid   = ncid_out, &
                                      name   = 'snapshot_times', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_snap_dimid, &
                                      varid  = nc_snaptime_varid) )

            if (trim(dump_type) == 'displ_only') then
               call check( nf90_def_dim( ncid   = ncid_meshout, &
                                         name   = 'elements', &
                                         len    = nelem_kwf_global, &
                                         dimid  = nc_mesh_elem_dimid) )
               call check( nf90_put_att( ncid   = ncid_out, &
                                         varid  = NF90_GLOBAL, &
                                         name   = 'nelem_kwf_global', &
                                         values = nelem_kwf_global) )

               call check( nf90_def_dim( ncid   = ncid_meshout, &
                                         name   = 'control_points', &
                                         len    = 4, &
                                         dimid  = nc_mesh_cntrlpts_dimid) )

               call check( nf90_def_dim( ncid   = ncid_meshout, &
                                         name   = 'npol', &
                                         len    = npol+1, &
                                         dimid  = nc_mesh_npol_dimid) )
            endif

            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_S', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_s_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_Z', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_z_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_vp', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_vp_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_vs', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_vs_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   ='mesh_rho', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_rho_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_lambda', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_lambda_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_mu', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_mu_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_xi', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_xi_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_phi', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_phi_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_eta', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid, &
                                      varid  = nc_mesh_eta_varid) )

            if (anel_true) then

                call check( nf90_def_var( ncid   = ncid_meshout, &
                                          name   = 'mesh_Qmu', &
                                          xtype  = NF90_FLOAT, &
                                          dimids = nc_pt_dimid, &
                                          varid  = nc_mesh_Qmu_varid) )
                call check( nf90_def_var( ncid   = ncid_meshout, &
                                          name   = 'mesh_Qka', &
                                          xtype  = NF90_FLOAT, &
                                          dimids = nc_pt_dimid, &
                                          varid  = nc_mesh_Qka_varid) )
            endif

            if (trim(dump_type) == 'displ_only') then
               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'midpoint_mesh', &
                                         xtype  = NF90_INT, &
                                         dimids = nc_mesh_elem_dimid, &
                                         varid  = nc_mesh_midpoint_varid) )

               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'eltype', &
                                         xtype  = NF90_INT, &
                                         dimids = nc_mesh_elem_dimid, &
                                         varid  = nc_mesh_eltype_varid) )

               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'axis', &
                                         xtype  = NF90_INT, &
                                         dimids = nc_mesh_elem_dimid, &
                                         varid  = nc_mesh_axis_varid) )

               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'fem_mesh', &
                                         xtype  = NF90_INT, &
                                         dimids = [nc_mesh_cntrlpts_dimid, &
                                                   nc_mesh_elem_dimid], &
                                         varid  = nc_mesh_fem_varid) )

               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'sem_mesh', &
                                         xtype  = NF90_INT, &
                                         dimids = [nc_mesh_npol_dimid, &
                                                   nc_mesh_npol_dimid, &
                                                   nc_mesh_elem_dimid], &
                                         varid  = nc_mesh_sem_varid) )

               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'mp_mesh_S', &
                                         xtype  = NF90_FLOAT, &
                                         dimids = nc_mesh_elem_dimid, &
                                         varid  = nc_mesh_s_mp_varid) )
               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'mp_mesh_Z', &
                                         xtype  = NF90_FLOAT, &
                                         dimids = nc_mesh_elem_dimid, &
                                         varid  = nc_mesh_z_mp_varid) )

               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'G0', &
                                         xtype  = NF90_DOUBLE, &
                                         dimids = nc_mesh_npol_dimid, &
                                         varid  = nc_mesh_G0_varid) )
               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'G1', &
                                         xtype  = NF90_DOUBLE, &
                                         dimids = [nc_mesh_npol_dimid, &
                                                   nc_mesh_npol_dimid], &
                                         varid  = nc_mesh_G1_varid) )
               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'G2', &
                                         xtype  = NF90_DOUBLE, &
                                         dimids = [nc_mesh_npol_dimid, &
                                                   nc_mesh_npol_dimid], &
                                         varid  = nc_mesh_G2_varid) )
               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'gll', &
                                         xtype  = NF90_DOUBLE, &
                                         dimids = nc_mesh_npol_dimid, &
                                         varid  = nc_mesh_gll_varid) )
               call check( nf90_def_var( ncid   = ncid_meshout, &
                                         name   = 'glj', &
                                         xtype  = NF90_DOUBLE, &
                                         dimids = nc_mesh_npol_dimid, &
                                         varid  = nc_mesh_glj_varid) )
            endif

            if (nc_chunk_time_traces) then
                chunk_pt = max(1, disk_block_size / nstrain)
                chunk_time = nstrain
            else
                chunk_pt = npoints_global
                chunk_time = 1
            endif

            do ivar=1, nvar/2 ! The big snapshot variables for the kerner.

                call check( nf90_def_var(ncid       = ncid_snapout, &
                                         name       = trim(nc_varnamelist(ivar)), &
                                         xtype      = NF90_FLOAT, &
                                         dimids     = [nc_pt_dimid, nc_snap_dimid], &
                                         varid      = nc_field_varid(ivar), &
                                         chunksizes = [chunk_pt, chunk_time] ))

                call check( nf90_def_var_fill(ncid    = ncid_snapout, &
                                              varid   = nc_field_varid(ivar), &
                                              no_fill = 1, &
                                              fill    = 0) )
            enddo


            ! Surface group in output file
            if (verbose > 2) write(*,*) 'Define variables in ''Surface'' group of NetCDF output file'
            call check( nf90_put_att( ncid   = ncid_surfout, &
                                      name   = 'nstrain', &
                                      varid  = NF90_GLOBAL, &
                                      values = nstrain) )
            call check( nf90_def_dim( ncid_surfout, "straincomponents", len=6, &
                                      dimid=nc_strcomp_dimid) )

            call check( nf90_def_dim( ncid_surfout, "surf_elems", maxind_glob, nc_surf_dimid) )
            call check( nf90_put_att( ncid   = ncid_surfout, &
                                      name   = 'nsurfelem', &
                                      varid  = NF90_GLOBAL, &
                                      values = maxind_glob) )

            call check( nf90_def_var( ncid_surfout, "elem_theta", NF90_FLOAT, &
                                      [nc_surf_dimid ], nc_surfelem_theta_varid) )
            call check( nf90_put_att(ncid_surfout, nc_surfelem_theta_varid, 'units', 'degrees'))

            call check( nf90_def_var( ncid_surfout, "displacement", NF90_FLOAT, &
                                      [nc_snap_dimid, nc_comp_dimid, nc_surf_dimid ], &
                                      nc_surfelem_disp_varid) )
            call check( nf90_put_att(ncid_surfout, nc_surfelem_disp_varid, 'units', 'meters'))

            call check( nf90_def_var( ncid_surfout, "velocity", NF90_FLOAT, &
                                      [nc_snap_dimid, nc_comp_dimid, nc_surf_dimid ], &
                                      nc_surfelem_velo_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_velo_varid, 'units', &
                                      'meters per second') )

            call check( nf90_def_var( ncid_surfout, "disp_src", NF90_FLOAT, &
                                      [nc_snap_dimid, nc_comp_dimid, nc_surf_dimid ], &
                                      nc_surfelem_disp_src_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_disp_src_varid, 'units', &
                                      'meters') )

            call check( nf90_def_var( ncid_surfout, "strain", NF90_FLOAT, &
                                      [nc_snap_dimid, nc_strcomp_dimid, nc_surf_dimid ], &
                                      nc_surfelem_strain_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_strain_varid, 'units', &
                                      ' ') )

            call check( nf90_def_var( ncid_surfout, "stf_dump", NF90_FLOAT, &
                                      [nc_snap_dimid], nc_stf_dump_varid) )

            call check( nf90_def_var( ncid_surfout, "stf_d_dump", NF90_FLOAT, &
                                      [nc_snap_dimid], nc_stf_d_dump_varid) )
        endif

        if (verbose > 1) write(*,'(I6,a/)') mynum, 'NetCDF variables defined'


        call check( nf90_enddef(ncid_out))


! in case of parallel IO, only the first rank writes
!#ifdef enable_parallel_netcdf
!    if (mynum == 0) then
!#endif
        if (verbose > 1) then
            write(*,*) 'Writing station info into NetCDF file...'
            call flush(6)
        endif
        call check( nf90_put_var( ncid_recout, nc_th_varid,   values = rec_th) )
        call check( nf90_put_var( ncid_recout, nc_ph_varid,   values = rec_ph) )
        call check( nf90_put_var( ncid_recout, nc_thr_varid,  values = rec_th_req) )
        call check( nf90_put_var( ncid_recout, nc_proc_varid, values = rec_proc) )

        do irec=1,nrec
            call check( nf90_put_var( ncid_recout, nc_recnam_varid, start = [irec, 1], &
                                      count = [1, 40], values = (rec_names(irec))) )
        enddo

        ! Write out seismogram dump times
        time_seis = dble([ (i, i = 0, nseismo-1) ]) * deltat
        call check( nf90_put_var( ncid_recout, nc_time_varid, values = time_seis ) )
        if (verbose > 1) then
            write(*,*) '...done'
            call flush(6)
        endif

        ! Write out STFs
        if (verbose > 1) then
            write(*,*) 'Writing stf into NetCDF file...'
            call flush(6)
        endif

        call check( nf90_put_var(ncid   = ncid_recout, &
                                 varid  = nc_stf_iter_varid, &
                                 values = stf_dumpvar) )
        call check( nf90_put_var(ncid   = ncid_recout, &
                                 varid  = nc_stf_d_iter_varid, &
                                 values = stf_d_dumpvar) )
        call check( nf90_put_var(ncid   = ncid_recout, &
                                 varid  = nc_stf_seis_varid, &
                                 values = stf_seis_dumpvar) )
        call check( nf90_put_var(ncid   = ncid_recout, &
                                 varid  = nc_stf_d_seis_varid, &
                                 values = stf_d_seis_dumpvar) )

        if (dump_wavefields) then
            ! Write out strain dump times
            if (verbose > 1) then
                write(*,*) 'Writing strain dump times into NetCDF file...'
                call flush(6)
            endif
            time_strain = dble([ (i, i = 1, nstrain) ])
            time_strain = time_strain * t_0 / strain_samp
            call check( nf90_put_var(ncid   = ncid_out, &
                                     varid  = nc_snaptime_varid, &
                                     values = time_strain ) )

            ! Write out STF values at kernel dump points
            if (verbose > 1) then
                write(*,*) 'Writing STF in strain dumps'
                call flush(6)
            endif
            call check( nf90_put_var(ncid   = ncid_surfout, &
                                     varid  = nc_stf_dump_varid, &
                                     values = stf_dump_dumpvar) )
            call check( nf90_put_var(ncid   = ncid_surfout, &
                                     varid  = nc_stf_d_dump_varid, &
                                     values = stf_d_dump_dumpvar) )

            if (verbose > 1) then
                write(*,*) '...done'
                call flush(6)
            endif
        endif

    endif ! (mynum == 0)


! Allocation of Dump buffer variables. Done on all procs
90  format(' Allocated ', A20, ', uses ',F10.3,'MB')
    allocate(recdumpvar(nseismo,3,num_rec))
    recdumpvar = 0.0

    if (dump_wavefields) then
        stepstodump = 0
        allocate(surfdumpvar_disp(dumpstepsnap,3,maxind))
        allocate(surfdumpvar_velo(dumpstepsnap,3,maxind))
        allocate(surfdumpvar_strain(dumpstepsnap,6,maxind))
        allocate(surfdumpvar_srcdisp(dumpstepsnap,3,maxind))

        if (src_type(1) == 'monopole') then
            allocate(oneddumpvar(npoints, dumpstepsnap, 6))
        else
            allocate(oneddumpvar(npoints, dumpstepsnap, 9))
        endif

        surfdumpvar_disp = 0.0
        surfdumpvar_velo = 0.0
        surfdumpvar_strain = 0.0
        surfdumpvar_srcdisp = 0.0

        if (mynum == 0 .and. verbose > 2) then
            write(*,*)  'Allocating NetCDF buffer variables'
            write(*,90) 'recdumpvar', real(size(recdumpvar))/262144.
            write(*,90) 'surfdumpvar_disp', real(size(surfdumpvar_disp))/262144.
            write(*,90) 'surfdumpvar_velo', real(size(surfdumpvar_velo))/262144.
            write(*,90) 'surfdumpvar_strain', real(size(surfdumpvar_strain))/262144.
            write(*,90) 'surfdumpvar_srcdisp', real(size(surfdumpvar_srcdisp))/262144.
            write(*,90) 'oneddumpvar', real(size(oneddumpvar))/262144.
            write(*,*)
            write(*,*)  '*********************************************************************'
            write(*,*)  '**** NetCDF output file produced, buffer variables all allocated ****'
            write(*,*)  '*********************************************************************'
            write(*,*)
        endif
    endif

#ifdef enable_parallel_netcdf
    if (mynum == 0) call check( nf90_close(ncid = ncid_out))
    call barrier()

    if (verbose > 0 .and. mynum == 0) then
        write(*,*) ' Opening netcdf file for parallel IO'
        call flush(6)
    endif
    nmode = ior(NF90_WRITE, NF90_NETCDF4)
    nmode = ior(nmode, NF90_MPIIO)
    call check( nf90_open(path=nc_fnam, mode=nmode, ncid=ncid_out, &
                          comm=MPI_COMM_WORLD, info=MPI_INFO_NULL) )
    call barrier()
    if (verbose > 0 .and. mynum == 0) then
        write(*,*) ' Netcdf file with ID ', ncid_out, ' opened on all procs.'
        call flush(6)
    endif
#endif
    call barrier()

#endif
end subroutine nc_define_outputfile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Open the NetCDF output file, check for variable IDs and dump meshes.
subroutine nc_finish_prepare
#ifdef enable_netcdf
    use data_io, only: datapath, lfdata, dump_wavefields, dump_type
    use data_mesh, only: maxind, surfcoord, ind_first, ind_last, &
                           midpoint_mesh_kwf, sem_mesh_kwf, fem_mesh_kwf, nelem_kwf, &
                           nelem_kwf_global, npol, eltype_kwf, axis_kwf, num_rec
    use data_spec, only: G0, G1, G2, xi_k, eta

    integer             :: ivar, nmode, iproc
    integer             :: nc_mesh_s_varid, nc_mesh_z_varid
    integer             :: nc_mesh_s_mp_varid, nc_mesh_z_mp_varid
    integer             :: nc_mesh_vs_varid, nc_mesh_vp_varid
    integer             :: nc_mesh_mu_varid, nc_mesh_rho_varid
    integer             :: nc_mesh_lambda_varid
    integer             :: nc_mesh_xi_varid
    integer             :: nc_mesh_phi_varid
    integer             :: nc_mesh_eta_varid
    integer             :: nc_mesh_qmu_varid
    integer             :: nc_mesh_qka_varid

    integer             :: nc_mesh_midpoint_varid
    integer             :: nc_mesh_eltype_varid
    integer             :: nc_mesh_axis_varid
    integer             :: nc_mesh_fem_varid
    integer             :: nc_mesh_sem_varid
    integer             :: nc_mesh_G0_varid
    integer             :: nc_mesh_G1_varid
    integer             :: nc_mesh_G2_varid
    integer             :: nc_mesh_gll_varid
    integer             :: nc_mesh_glj_varid

! in case of parallel IO, we keep it open on all procs
#ifndef enable_parallel_netcdf
    if (mynum == 0) then
        call check(nf90_close(ncid_out))
        if (verbose > 1) then
           write(*,*) '  Root process closed netCDF file, waiting for all procs to'
           write(*,*) '  arrive here and then open it to retrieve IDs'
           if (dump_wavefields) write(*,*) '  and dump mesh coordinates.'
        endif
    endif
    call barrier

    do iproc = 0, nproc
        call barrier
        if (iproc == mynum .and. (npoints > 0 .or. maxind > 0 .or. num_rec > 0 .or. mynum == 0)) then
            if (verbose > 1) then
               write(*,*) '  Processor ', iproc, ' opened the output file and will dump '
               write(*,*) '  his part of the mesh.'
            endif
            nmode = ior(NF90_WRITE, NF90_NETCDF4)
            call check( nf90_open( path = datapath(1:lfdata)//"/axisem_output.nc4", &
                                   mode = nmode, &
                                   ncid = ncid_out) )
            print '(A,I5,A)', '   ', iproc, ': opened file'
#endif


            call getgrpid(ncid_out, "Seismograms", ncid_recout)
            call getgrpid(ncid_out, "Surface", ncid_surfout)
            call getgrpid(ncid_out, "Mesh", ncid_meshout)
            !print '(A,I5,A)', '   ', iproc, ': inquired dimension IDs'
            call getvarid( ncid_recout, "displacement", nc_disp_varid )

            if (dump_wavefields) then

                ! first get all IDs
                call getgrpid(ncid_out, "Snapshots", ncid_snapout)
                do ivar=1, nvar/2
                   call getvarid(ncid_snapout, nc_varnamelist(ivar), nc_field_varid(ivar))
                enddo

                call getvarid(ncid_surfout, "elem_theta",   nc_surfelem_theta_varid)
                call getvarid(ncid_surfout, "displacement", nc_surfelem_disp_varid)
                call getvarid(ncid_surfout, "velocity",     nc_surfelem_velo_varid)
                call getvarid(ncid_surfout, "strain",       nc_surfelem_strain_varid)
                call getvarid(ncid_surfout, "disp_src",     nc_surfelem_disp_src_varid)

                call getvarid(ncid_meshout, "mesh_S",       nc_mesh_s_varid)
                call getvarid(ncid_meshout, "mesh_Z",       nc_mesh_z_varid)
                call getvarid(ncid_meshout, "mesh_vp",      nc_mesh_vp_varid)
                call getvarid(ncid_meshout, "mesh_vs",      nc_mesh_vs_varid)
                call getvarid(ncid_meshout, "mesh_rho",     nc_mesh_rho_varid)
                call getvarid(ncid_meshout, "mesh_lambda",  nc_mesh_lambda_varid)
                call getvarid(ncid_meshout, "mesh_mu",      nc_mesh_mu_varid)
                call getvarid(ncid_meshout, "mesh_phi",     nc_mesh_phi_varid)
                call getvarid(ncid_meshout, "mesh_xi",      nc_mesh_xi_varid)
                call getvarid(ncid_meshout, "mesh_eta",     nc_mesh_eta_varid)

                if (allocated(Q_mu1d) .and. allocated(Q_kappa1d)) then
                    call getvarid(ncid_meshout, "mesh_Qmu", nc_mesh_Qmu_varid)
                    call getvarid(ncid_meshout, "mesh_Qka", nc_mesh_Qka_varid)
                endif

                if (trim(dump_type) == 'displ_only') then
                   call getvarid(ncid_meshout, "midpoint_mesh", nc_mesh_midpoint_varid)
                   call getvarid(ncid_meshout, "eltype",        nc_mesh_eltype_varid)
                   call getvarid(ncid_meshout, "axis",          nc_mesh_axis_varid)
                   call getvarid(ncid_meshout, "fem_mesh",      nc_mesh_fem_varid)
                   call getvarid(ncid_meshout, "sem_mesh",      nc_mesh_sem_varid)
                   call getvarid(ncid_meshout, "mp_mesh_S",     nc_mesh_s_mp_varid)
                   call getvarid(ncid_meshout, "mp_mesh_Z",     nc_mesh_z_mp_varid)
                   call getvarid(ncid_meshout, "gll",           nc_mesh_gll_varid)
                   call getvarid(ncid_meshout, "glj",           nc_mesh_glj_varid)
                   call getvarid(ncid_meshout, "G0",            nc_mesh_G0_varid)
                   call getvarid(ncid_meshout, "G1",            nc_mesh_G1_varid)
                   call getvarid(ncid_meshout, "G2",            nc_mesh_G2_varid)
                endif

                !print '(A,I5,A)', '   ', iproc, ': inquired variable IDs'

#ifdef enable_parallel_netcdf
                ! enable collective IO in case of parallel IO
                do ivar=1, nvar/2
                   call check(nf90_var_par_access(ncid_snapout, nc_field_varid(ivar), &
                                                  NF90_COLLECTIVE))
                enddo
                call check(nf90_var_par_access(ncid_meshout, nc_surfelem_theta_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_surfelem_disp_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_surfelem_velo_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_surfelem_strain_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_surfelem_disp_src_varid, &
                                               NF90_COLLECTIVE))

                call check(nf90_var_par_access(ncid_meshout, nc_mesh_s_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_z_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_vp_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_vs_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_rho_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_lambda_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_mu_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_phi_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_xi_varid, &
                                               NF90_COLLECTIVE))
                call check(nf90_var_par_access(ncid_meshout, nc_mesh_eta_varid, &
                                               NF90_COLLECTIVE))

                if (allocated(Q_mu1d) .and. allocated(Q_kappa1d)) then
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_Qmu_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_Qka_varid, &
                                                  NF90_COLLECTIVE))
                endif

                if (trim(dump_type) == 'displ_only') then
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_midpoint_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_eltype_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_axis_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_fem_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_sem_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_s_mp_varid, &
                                                  NF90_COLLECTIVE))
                   call check(nf90_var_par_access(ncid_meshout, nc_mesh_z_mp_varid, &
                                                  NF90_COLLECTIVE))
                endif
#endif

                ! start writing to file

                call putvar_real1d(ncid   = ncid_surfout, &
                                   varid  = nc_surfelem_theta_varid, &
                                   values = surfcoord, &
                                   start  = ind_first, &
                                   count  = maxind  ) ! count = 0 should be valid?

                ! S-Coordinate
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_s_varid, &
                                    values = scoord1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Z-Coordinate
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_z_varid, &
                                    values = zcoord1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Vp
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_vp_varid, &
                                    values = vp1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Vs
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_vs_varid, &
                                    values = vs1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Rho
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_rho_varid, &
                                    values = rho1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Lambda
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_lambda_varid, &
                                    values = lambda1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Mu
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_mu_varid, &
                                    values = mu1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Anisotropic parameters
                ! Phi
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_phi_varid, &
                                    values = phi1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Xi
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_xi_varid, &
                                    values = xi1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Eta
                call putvar_real1d( ncid   = ncid_meshout, &
                                    varid  = nc_mesh_eta_varid, &
                                    values = eta1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Anelastic parameters
                if (allocated(Q_mu1d) .and. allocated(Q_kappa1d)) then

                    ! Q_mu
                    call putvar_real1d( ncid   = ncid_meshout, &
                                        varid  = nc_mesh_Qmu_varid, &
                                        values = Q_mu1d, &
                                        start  = npoints_myfirst, &
                                        count  = npoints )
                    ! Q_kappa
                    call putvar_real1d( ncid   = ncid_meshout, &
                                        varid  = nc_mesh_Qka_varid, &
                                        values = Q_kappa1d, &
                                        start  = npoints_myfirst, &
                                        count  = npoints )
                endif

                if (trim(dump_type) == 'displ_only') then
                   call check(nf90_put_var ( ncid   = ncid_meshout, &
                                             varid  = nc_mesh_midpoint_varid, &
                                             start  = [nelem_myfirst], &
                                             count  = [nelem_kwf], &
                                             values = midpoint_mesh_kwf + npoints_myfirst - 1))

                   call check(nf90_put_var ( ncid   = ncid_meshout, &
                                             varid  = nc_mesh_eltype_varid, &
                                             start  = [nelem_myfirst], &
                                             count  = [nelem_kwf], &
                                             values = eltype_kwf))

                   call check(nf90_put_var ( ncid   = ncid_meshout, &
                                             varid  = nc_mesh_axis_varid, &
                                             start  = [nelem_myfirst], &
                                             count  = [nelem_kwf], &
                                             values = axis_kwf))

                   call check(nf90_put_var ( ncid   = ncid_meshout, &
                                             varid  = nc_mesh_fem_varid, &
                                             start  = [1, nelem_myfirst], &
                                             count  = [4, nelem_kwf], &
                                             values = fem_mesh_kwf + npoints_myfirst - 1))

                   call check(nf90_put_var ( ncid   = ncid_meshout, &
                                             varid  = nc_mesh_sem_varid, &
                                             start  = [1, 1, nelem_myfirst], &
                                             count  = [npol+1, npol+1, nelem_kwf], &
                                             values = sem_mesh_kwf + npoints_myfirst - 1))

                   ! S-Coordinate
                   call putvar_real1d( ncid   = ncid_meshout, &
                                       varid  = nc_mesh_s_mp_varid, &
                                       values = scoord1d_mp, &
                                       start  = nelem_myfirst, &
                                       count  = nelem_kwf )

                   ! Z-Coordinate
                   call putvar_real1d( ncid   = ncid_meshout, &
                                       varid  = nc_mesh_z_mp_varid, &
                                       values = zcoord1d_mp, &
                                       start  = nelem_myfirst, &
                                       count  = nelem_kwf )

                   ! SEM stuff
                   if (mynum == 0) then ! these are not processorwise, so only rank zero
                                        ! writes
                      call check(nf90_put_var ( ncid   = ncid_meshout, &
                                                varid  = nc_mesh_gll_varid, &
                                                start  = [1], &
                                                count  = [npol+1], &
                                                values = eta))

                      call check(nf90_put_var ( ncid   = ncid_meshout, &
                                                varid  = nc_mesh_glj_varid, &
                                                start  = [1], &
                                                count  = [npol+1], &
                                                values = xi_k))

                      call check(nf90_put_var ( ncid   = ncid_meshout, &
                                                varid  = nc_mesh_G0_varid, &
                                                start  = [1], &
                                                count  = [npol+1], &
                                                values = G0))

                      call check(nf90_put_var ( ncid   = ncid_meshout, &
                                                varid  = nc_mesh_G1_varid, &
                                                start  = [1, 1], &
                                                count  = [npol+1, npol+1], &
                                                values = G1))

                      call check(nf90_put_var ( ncid   = ncid_meshout, &
                                                varid  = nc_mesh_G2_varid, &
                                                start  = [1, 1], &
                                                count  = [npol+1, npol+1], &
                                                values = G2))
                   endif
                endif

                !print '(A,I5,A)', '   ', iproc, ': dumped mesh'

            endif !dump_wavefields
            if (verbose > 1) &
                write(*,"('  Proc ', I3, ' dumped its mesh and is ready to rupture')") &
                    mynum

! in case of parallel IO, we keep it open on all procs
#ifndef enable_parallel_netcdf
            call check( nf90_close( ncid_out))
        endif !mynum == iproc
    enddo
#endif

#endif
end subroutine nc_finish_prepare
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Final dumps to netCDF file. In the moment contains only dump of
!! receiver seismograms.
subroutine nc_end_output
#ifdef enable_netcdf
    use data_mesh, only: num_rec
    use data_io, only: dump_xdmf
    integer           :: iproc

    call flush(6)
#ifndef enable_parallel_netcdf
    call barrier
    do iproc=0, nproc-1
        if (iproc == mynum) then
            if (num_rec > 0) then
                if (verbose > 1) write(*,"('   Proc ', I3, ' will dump receiver seismograms')") mynum
                call nc_dump_rec_to_disk()
            else
                if (verbose > 1) write(*,"('   Proc ', I3, ' has no receivers and just waits for the others')") mynum
            endif
        endif
        call barrier
    enddo
#else
    if (num_rec > 0) then
        if (verbose > 2) write(*,"('   Proc ', I3, ' will dump receiver seismograms')") mynum
        call nc_dump_rec_to_disk()
    endif
#endif

    !Set the finalized flag to true in the output file
#ifndef enable_parallel_netcdf
    if (mynum == 0) &
#endif
        call nc_finalize()

#endif
end subroutine nc_end_output
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_finalize
! < Set the finalized flag to true in the output file
#ifdef enable_netcdf
    use data_io, only: datapath, lfdata

#ifndef enable_parallel_netcdf
    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                          mode=NF90_WRITE, ncid=ncid_out) )
#endif
    call check( nf90_redef(ncid_out))

    call check( nf90_put_att(ncid_out, NF90_GLOBAL, &
                             'finalized', 1) )

    call check( nf90_close(ncid = ncid_out))
#endif
end subroutine nc_finalize
!-----------------------------------------------------------------------------------------

end module nc_routines
!=========================================================================================
