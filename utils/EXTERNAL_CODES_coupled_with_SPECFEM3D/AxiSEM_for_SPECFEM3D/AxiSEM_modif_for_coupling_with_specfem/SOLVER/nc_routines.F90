!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
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

!> Contains all the routines for NetCDF handling.
module nc_routines

#ifdef unc
    use netcdf
#endif
    use data_io,    only : verbose, deflate_level, nseismo
    use data_proc,  only : mynum, nproc, lpr
    use global_parameters
    use commun,     only : barrier, comm_elem_number

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
!    real, allocatable   :: oneddumpvar_flu(:,:,:)    !< Buffer variable for everything fluid dumped in nc_dump_field_1d
!    real, allocatable   :: oneddumpvar_sol(:,:,:)    !< Buffer variable for everything solid dumped in nc_dump_field_1d

    !> Buffer variable for everything dumped in nc_dump_field_1d
    real(sp), allocatable   :: oneddumpvar(:,:,:)
    real(sp), allocatable   :: scoord1d(:), zcoord1d(:)
    real(sp), allocatable   :: rho1d(:), mu1d(:), lambda1d(:)
    real(sp), allocatable   :: vp1d(:), vs1d(:)

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

    integer             :: ncid_out, ncid_recout, ncid_snapout, ncid_surfout, ncid_meshout
    integer             :: nc_snap_dimid, nc_proc_dimid, nc_rec_dimid, nc_recproc_dimid
    integer             :: nc_times_dimid, nc_comp_dimid, nc_disp_varid, nc_stf_seis_varid
    integer             :: nc_time_varid, nc_iter_dimid, nc_stf_iter_varid

    integer             :: nc_strcomp_dimid
    integer             :: nc_surfelem_disp_varid, nc_surfelem_velo_varid
    integer             :: nc_surfelem_strain_varid, nc_surfelem_disp_src_varid
    integer             :: nc_mesh_sol_varid, nc_mesh_flu_varid, nc_stf_dump_varid
    integer             :: nc_point_dimid, nc_pt_sol_dimid, nc_pt_flu_dimid
    integer             :: nc_szcoord_dimid
    integer             :: nc_snaptime_varid, nc_elem_dom_varid, nc_surfelem_theta_varid
    integer,allocatable :: nc_field_varid(:)
    character(len=16), allocatable  :: varnamelist(:)
    character(len=12), allocatable  :: nc_varnamelist(:)
    integer             :: nvar = -1

    !! Variables for dumping of wavefields for plotting purposes
    integer             :: nc_snap_disp_varid, nc_coord_dimid
    integer             :: nc_snap_point_varid, nc_snap_grid_varid
    integer             :: nc_snap_pwave_varid, nc_snap_swave_varid
    integer             :: ncid_out_snap
    integer             :: ndim_disp !< 2 for monopole, 3 for rest


    !! Buffer variables for the STF.
    real(kind=sp), allocatable :: stf_dump_dumpvar(:)
    real(kind=sp), allocatable :: stf_seis_dumpvar(:)
    real(kind=sp), allocatable :: stf_dumpvar(:)

    !! @todo These parameters should move to a input file soon
    !> How many snaps should be buffered in RAM?
    integer             :: dumpbuffersize = 256

    public              :: nc_dump_strain, nc_dump_rec, nc_dump_surface
    public              :: nc_dump_field_solid, nc_dump_field_fluid
    public              :: nc_write_att_char, nc_write_att_real, nc_write_att_int
    public              :: nc_write_att_dble
    public              :: nc_define_outputfile, nc_finish_prepare, nc_end_output
    public              :: nc_dump_strain_to_disk, nc_dump_mesh_sol, nc_dump_mesh_flu
    public              :: nc_write_el_domains, nc_dump_elastic_parameters
    public              :: nc_dump_snapshot, nc_dump_snap_points, nc_dump_snap_grid
    public              :: nc_make_snapfile, nc_dump_stf, nc_rec_checkpoint
contains

!-----------------------------------------------------------------------------------------
!> Write out the model domains for each element
subroutine nc_write_el_domains(idom)
    integer, dimension(:), intent(in) :: idom

#ifdef unc
    call check( nf90_put_var(ncid   = ncid_snapout, &
                             varid  = nc_elem_dom_varid, &
                             values = idom) )
#endif
end subroutine

!-----------------------------------------------------------------------------------------
!> Routine to dump the wavefield variables for the Kerner. Collects input in
!! oneddumpvar_sol and oneddumpvar_flu until dumping condition is fulfilled.
subroutine nc_dump_field_solid(f, varname)

    real(kind=realkind), intent(in)   :: f(:)     !< Data to dump, size should be npts_sol
    character(len=*), intent(in)      :: varname  !< Internal name of data to dump.
    !! Is used to identify the NetCDF Variable in question
#ifdef unc
    integer                           :: ivar

    do ivar=1, nvar
        !< Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar)) == trim(varname)) exit
    enddo

    if (ivar > nvar/2) then
        write(6,*) 'nc_dump_field_solid: Trying to access variable: ', trim(varname), &
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

    real(kind=realkind), intent(in)   :: f(:)  !< Data to dump, size should be npts_flu
    character(len=*), intent(in)      :: varname  !< Internal name of data to dump.
    !! Is used to identify the NetCDF Variable in question
#ifdef unc
    integer                           :: ivar


    do ivar=1, nvar
        !< Check whether this Variable actually exists in file ncid_out
        if (trim(varnamelist(ivar)) == trim(varname)) exit
    enddo

    if (ivar <= nvar/2) then
        write(6,*) 'nc_dump_field_fluid: Trying to access variable: ', trim(varname), &
            ' which is a solid variable. Contact a developer and shout at him!'
        stop 1
    endif

    oneddumpvar(npts_sol+1:npoints,stepstodump+1,ivar-nvar/2) = f
#endif
end subroutine nc_dump_field_fluid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain(isnap_loc)

    use data_io, ONLY    : nstrain
    use clocks_mod, only : tick
    use data_time, only  : iclocknbio, idnbio

    ! explicit interfaces to the c functions to avoide the underscore issues
    ! (fortran 2003 standard)
    interface
        subroutine c_spawn_dumpthread(stepstodump) bind(c, name='c_spawn_dumpthread')
            use, intrinsic      :: iso_c_binding, only : c_int
            integer (c_int)     :: stepstodump
        end subroutine
    end interface

    interface
        subroutine c_wait_for_io() bind(c, name='c_wait_for_io')
        end subroutine
    end interface

    integer, intent(in) :: isnap_loc
#ifdef unc
    integer             :: iproc
    real                :: tickl, tackl

    stepstodump = stepstodump + 1
    if (isnap_loc == 0) return

    if (dumpposition(mod(isnap_loc, dumpstepsnap))) then

        ! wait for other processes to finish writing, measure waiting time and
        ! issue warning in case waiting longer then .5 sec
        ! MvD: I am not sure if cpu_time is a correct measure here, as we idle
        !      until IO is finished.
        !      Therefor testing system_clock, from the clocks module.

        if (mynum == 0) call cpu_time(tickl)

        iclocknbio = tick()
        do iproc=0, nproc-1
            if (iproc == mynum) then
                call c_wait_for_io()
                call flush(6)
            endif
            call barrier
        enddo
        iclocknbio = tick(id=idnbio, since=iclocknbio)

        if (mynum == 0) then
            call cpu_time(tackl)
            if ((tackl-tickl) > 0.5 .and. verbose > 0) then
                write(6,"('WARNING: Computation was halted for ', F7.2, ' s to wait for ',&
                         & 'dumping processor. Consider adapting netCDF output variables',&
                         & '(disable compression, increase dumpstepsnap)')") tackl-tickl
            endif
        endif

        ! non blocking write
        if (mod(isnap_loc, dumpstepsnap) == outputplan) then
            call c_wait_for_io()
            isnap_global = isnap_loc
            ndumps = stepstodump
            call c_spawn_dumpthread(stepstodump)
            stepstodump = 0
        endif
    endif

    ! Final and last dump of all remaining
    ! not done in unblocking fashion using a thread, as there is nothing to
    ! compute anymore
    ! Make sure, nobody is accessing the output file anymore
    if (isnap_loc == nstrain) then
        do iproc=0, nproc-1
            if (iproc == mynum) then
                call c_wait_for_io()
                call barrier
            endif
        enddo
        do iproc=0,nproc-1
            if (iproc == mynum) then
                isnap_global = nstrain
                ndumps = stepstodump
                call nc_dump_strain_to_disk()
                if (verbose > 1) write(6,*) mynum, 'finished dumping strain'
            endif
            call barrier
        enddo
    endif

#endif
end subroutine nc_dump_strain
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_strain_to_disk() bind(c, name="nc_dump_strain_to_disk")
#ifdef unc

    use data_io
    use data_source,       ONLY: src_type
    use global_parameters, ONLY: realkind
    use data_mesh,         ONLY: loc2globrec, maxind, ind_first

    integer                           :: ivar, flen, isnap_loc
    real                              :: tick, tack
    integer                           :: dumpsize

    dumpsize = 0
    call cpu_time(tick)

    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                          mode=NF90_WRITE, ncid=ncid_out) )

    isnap_loc = isnap_global
    if (verbose > 1) then
        if (ndumps == 0) then
            write(6,"('  Proc ', I4, ' in dump routine, isnap =', I5, &
                    & ', nothing to dump, returning...')") mynum, isnap_loc
            return
        else
            write(6,"('  Proc ', I4, ' in dump routine, isnap =', I5, &
                    & ', stepstodump = ', I4)") mynum, isnap_loc, ndumps
        endif
    endif


    do ivar=1, nvar/2
        call check( nf90_put_var(ncid=ncid_snapout, varid=nc_field_varid(ivar), &
                                 start=[npoints_myfirst, isnap_loc-ndumps+1], &
                                 count=[npoints, ndumps], &
                                 values=oneddumpvar(1:npoints,1:ndumps,ivar)) )
        dumpsize = dumpsize + npoints * ndumps
    enddo

    dumpsize = dumpsize + flen * ndumps

    !> Surface dumps
    call check( nf90_put_var(ncid_surfout, nc_surfelem_disp_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 3, maxind], &
                values = surfdumpvar_disp(:,:,:)) )
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_put_var(ncid_surfout, nc_surfelem_velo_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 3, maxind], &
                values = surfdumpvar_velo(:,:,:)) )
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_put_var(ncid_surfout, nc_surfelem_strain_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 6, maxind], &
                values = surfdumpvar_strain(:,:,:)) )
    dumpsize = dumpsize + 6 * maxind * ndumps

    call check( nf90_put_var(ncid_surfout, nc_surfelem_disp_src_varid, &
                start = [isnap_loc-ndumps+1, 1, ind_first], &
                count = [ndumps, 3, maxind], &
                values = surfdumpvar_srcdisp(:,:,:)) )
    dumpsize = dumpsize + 3 * maxind * ndumps

    call check( nf90_close(ncid_out) )
    call cpu_time(tack)

    if (verbose > 1) then
        write(6,"('  Proc', I5,': Wrote ', F8.2, ' MB in ', F6.2, 's')") &
            mynum, real(dumpsize) * realkind / 1048576., tack-tick
        call flush(6)
    endif

#endif
end subroutine nc_dump_strain_to_disk
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine nc_dump_stf(stf)
    use data_io,  only                       : nseismo, nstrain
    use data_time, only                      : seis_it, strain_it, niter
    real(kind=sp), intent(in), dimension(:) :: stf
#ifdef unc
    integer                                 :: it_s, it_d, i

    allocate(stf_dumpvar(niter))
    allocate(stf_seis_dumpvar(nseismo))
    allocate(stf_dump_dumpvar(nstrain))
    stf_seis_dumpvar = 0.0
    it_s = 1
    it_d = 1
    do i = 1, niter
        if ( mod(i,seis_it) == 0) then
           stf_seis_dumpvar(it_s) = stf(i)
           it_s = it_s + 1
        endif
        if ( mod(i,strain_it) == 0) then
           stf_dump_dumpvar(it_d) = stf(i)
           it_d = it_d + 1
        endif
    enddo
    stf_dumpvar = stf

#endif
end subroutine nc_dump_stf
!----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Dump receiver specific stuff, especially displacement and velocity
!! N.B.: Works with global indices.
subroutine nc_dump_rec(recfield)
    use data_mesh, ONLY: num_rec
    use data_io,   ONLY: iseismo
    real(sp), intent(in), dimension(3,num_rec) :: recfield
#ifdef unc

    recdumpvar(iseismo,:,:) = 0.0
    where(abs(recfield)>epsi) recdumpvar(iseismo,:,:) = recfield(:,:)

#endif
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_rec_to_disk
#ifdef unc
    use data_mesh, ONLY: loc2globrec, num_rec
    use data_io,   ONLY: datapath, lfdata, nseismo

    real                              :: tick, tack
    integer                           :: irec, dumpsize, icomp

    call cpu_time(tick)

    call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                          mode=NF90_WRITE, ncid=ncid_out) )
    call getvarid( ncid_recout, "displacement", nc_disp_varid )

    dumpsize = 0
    do irec = 1, num_rec
        do icomp = 1, 3
            call check( nf90_put_var(ncid   = ncid_recout, varid=nc_disp_varid, &
                                     start  = [1, icomp, loc2globrec(irec)], &
                                     count  = [nseismo, 1, 1], &
                                     values = recdumpvar(:,icomp,irec)) )
        enddo
    enddo

    dumpsize = nseismo * 3 * num_rec
    call check( nf90_close(ncid=ncid_out))
    call cpu_time(tack)

    if (verbose > 1) then
        write(6,"(I3,': Receiver data, Wrote ', F8.3, ' MB in ', F6.2, 's')") &
            mynum, real(dumpsize) * realkind / 1048576., tack-tick
    endif
#endif
end subroutine nc_dump_rec_to_disk
!----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_rec_checkpoint
    use data_mesh, ONLY: loc2globrec, num_rec
#ifdef unc
    integer         :: iproc
    do iproc=0, nproc-1
        call barrier
        if (iproc == mynum) then
            if (num_rec>0) then
                if (verbose > 1) write(6,"('Proc ', I3, ' will dump receiver seismograms')") mynum
                call nc_dump_rec_to_disk()
            else
                if (verbose > 1) write(6,"('Proc ', I3, ' has no receivers and just waits for the others')") mynum
            endif
        endif

    enddo
    call barrier
#endif
end subroutine

!----------------------------------------------------------------------------------------
!> Dump stuff along surface
subroutine nc_dump_surface(surffield, disporvelo)!, nrec, dim2)
    use data_mesh, ONLY: maxind

    !integer, intent(in)                          :: nrec, dim2
    real(kind=realkind), intent(in), dimension(:,:) :: surffield
    character(len=4), intent(in)                 :: disporvelo
#ifdef unc

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
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_sol(scoord_sol, zcoord_sol)

    use data_io,   ONLY : ndumppts_el
    use data_mesh, ONLY : nelem
    real(sp), intent(in), dimension(:,:,:) :: scoord_sol, zcoord_sol
#ifdef unc

    npts_sol = size(scoord_sol)
    if (npts_sol/=size(scoord_sol)) then
        write(6,*) 'Inconsistency in mesh size in nc_dump_mesh_solid'
        write(6,*) 'npts_sol=', npts_sol, ', size(scoord_sol)=', size(scoord_sol)
        stop 2
    endif

    npoints = ndumppts_el * nelem
    allocate(scoord1d(npoints))
    allocate(zcoord1d(npoints))

    zcoord1d(1:npts_sol) = pack(zcoord_sol,.true.)
    scoord1d(1:npts_sol) = pack(scoord_sol,.true.)


#endif
end subroutine nc_dump_mesh_sol
!-----------------------------------------------------------------------------------------
subroutine nc_dump_elastic_parameters(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                      fa_ani_theta, fa_ani_phi, Q_mu, Q_kappa)

    use data_io,                                   only: ibeg, iend
    real(kind=dp), dimension(:,:,:), intent(in)       :: rho, lambda, mu, xi_ani
    real(kind=dp), dimension(:,:,:), intent(in)       :: phi_ani, eta_ani
    real(kind=dp), dimension(:,:,:), intent(in)       :: fa_ani_theta, fa_ani_phi
    real(kind=sp), dimension(:), intent(in), optional :: Q_mu, Q_kappa
    integer                                           :: size1d

    !print *, 'Processor', mynum,' has been here'
    !allocate(rho1d(npoints))
    !allocate(lambda1d(npoints))
    !allocate(mu1d(npoints))
    !allocate(vp1d(npoints))
    !allocate(vs1d(npoints))
    size1d = size(rho(ibeg:iend, ibeg:iend, :))
    print *, ' NetCDF: Mesh elastic parameter variables have size:', size1d
    allocate(rho1d(size1d))
    allocate(lambda1d(size1d))
    allocate(mu1d(size1d))
    allocate(vp1d(size1d))
    allocate(vs1d(size1d))

    rho1d     = real(pack(rho(ibeg:iend, ibeg:iend, :)    ,.true.), kind=sp)
    lambda1d  = real(pack(lambda(ibeg:iend, ibeg:iend, :) ,.true.), kind=sp)
    mu1d      = real(pack(mu(ibeg:iend, ibeg:iend, :)     ,.true.), kind=sp)
    vp1d      = sqrt( (lambda1d + 2.*mu1d ) / rho1d  )
    vs1d      = sqrt( mu1d  / rho1d )

    !if (present(Q_mu) .and. present(Q_kappa)) then

end subroutine nc_dump_elastic_parameters
!-----------------------------------------------------------------------------------------
subroutine nc_dump_mesh_flu(scoord_flu, zcoord_flu)

    real(sp), intent(in), dimension(:,:,:) :: scoord_flu, zcoord_flu
#ifdef unc

    npts_flu = size(scoord_flu)
    if (npts_flu/=size(scoord_flu))  then
        write(6,*) 'Inconsistency in mesh size in nc_dump_mesh_flu'
        stop 2
    endif

    zcoord1d(npts_sol+1:) = pack(zcoord_flu, .true.)
    scoord1d(npts_sol+1:) = pack(scoord_flu, .true.)


#endif
end subroutine nc_dump_mesh_flu
!-----------------------------------------------------------------------------------------
!subroutine nc_dump_mesh_to_disk()
!    use data_io, only : datapath, lfdata
!    integer iproc, nc_mesh_s_varid, nc_mesh_z_varid
!
!    do iproc=0, nproc-1
!        if (mynum == iproc) then
!            call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
!                                  mode=NF90_WRITE, ncid=ncid_out) )
!
!            call check( nf90_inq_varid( ncid_snapout, "mesh S", nc_mesh_s_varid ) )
!            call check( nf90_put_var(ncid_snapout, varid = nc_mesh_s_varid, &
!                                     values = scoord1d, &
!                                     start  = (/1,mynum*npoints+1/), &
!                                     count  = (/1,npoints/) ))
!            call check( nf90_inq_varid( ncid_snapout, "mesh Z", nc_mesh_z_varid ) )
!            call check( nf90_put_var(ncid_snapout, varid = nc_mesh_z_varid, &
!                                     values = zcoord1d, &
!                                     start  = (/1,mynum*npoints+1/), &
!                                     count  = (/1,npoints/) ))
!
!            call check( nf90_close(ncid_out))
!        endif
!        call barrier
!    enddo
!
!
!end subroutine nc_dump_mesh_to_disk

!-----------------------------------------------------------------------------------------
!> Define the output file variables and dimensions
!! and allocate buffer variables.
subroutine nc_define_outputfile(nrec, rec_names, rec_th, rec_th_req, rec_ph, rec_proc)

    use data_io,     ONLY: nseismo, nstrain, nseismo, ibeg, iend, dump_wavefields
    use data_io,     ONLY: datapath, lfdata, strain_samp
    use data_mesh,   ONLY: maxind, num_rec, discont, nelem, nel_solid, nel_fluid, &
                           ndisc, maxind_glob
    use data_source, ONLY: src_type, t_0
    use data_time,   ONLY: deltat, niter


    integer, intent(in)                  :: nrec              !< Number of receivers
    character(len=40),intent(in)         :: rec_names(nrec)   !< Receiver names
    real(dp), dimension(nrec),intent(in) :: rec_th            !< Receiver theta
    real(dp), dimension(nrec),intent(in) :: rec_th_req        !< Requested receiver theta
    real(dp), dimension(nrec),intent(in) :: rec_ph            !< Receiver phi
    integer, dimension(nrec),intent(in)  :: rec_proc          !< Receiver processor
#ifdef unc
    real(dp), dimension(nstrain)         :: time_strain
    real(dp), dimension(nseismo)         :: time_seis
    character(len=16), allocatable       :: varname(:)
    integer                              :: ivar, i
    integer                              :: irec, iproc, nmode
    integer                              :: nc_latr_varid, nc_lon_varid
    integer                              :: nc_lat_varid, nc_ph_varid
    integer                              :: nc_thr_varid, nc_th_varid
    integer                              :: nc_proc_varid, nc_recnam_dimid
    integer                              :: nc_recnam_varid, nc_surf_dimid
    integer                              :: nc_pt_dimid
    integer                              :: nc_mesh_s_varid, nc_mesh_z_varid
    integer                              :: nc_mesh_vs_varid, nc_mesh_vp_varid
    integer                              :: nc_mesh_mu_varid, nc_mesh_rho_varid
    integer                              :: nc_mesh_lambda_varid
    integer                              :: nc_disc_dimid, nc_disc_varid

    if ((mynum == 0) .and. (verbose > 1)) then
        write(6,*)
        write(6,*) '************************************************************************'
        write(6,*) '**** Producing netcdf output file with its variables and dimensions ****'
        write(6,*) '************************************************************************'
        write(6,*)
    endif

    call barrier
    dumpstepsnap = int(dumpbuffersize / nproc) * nproc ! Will later be reduced to nstrain, if this is smaller
                                                       ! than value given here

    if (lpr .and. verbose > 1) write(6,*) '  Dumping NetCDF file to disk every', dumpstepsnap, ' snaps'

    call barrier ! for nicer output only

    outputplan = mynum * (dumpstepsnap / nproc)

    allocate(dumpposition(0:dumpstepsnap-1))
    dumpposition = .false.
    do iproc=0, nproc-1
        dumpposition(iproc*(dumpstepsnap/nproc)) = .true.
        if ((iproc == mynum) .and. (verbose > 1)) then
            write(6,"('Proc ', I4, ' will dump at position ', I4)") mynum, outputplan
            call flush(6)
        endif
        call barrier ! for nicer output only
    enddo


    if (mynum == 0) then
        if (verbose > 1) write (6,*) 'Preparing netcdf file for ', nproc, ' processors'
        nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
        call check( nf90_create(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                                cmode=nmode, ncid=ncid_out) )
        if (verbose > 1) write(6,*) ' Netcdf file with ID ', ncid_out, ' produced.'
    endif


    if (dump_wavefields) then
        if (src_type(1) == 'monopole') then
            nvar = 12
        else
            nvar = 18
        endif
        allocate(varname(nvar))
        allocate(varnamelist(nvar))
        allocate(nc_varnamelist(nvar/2))
        allocate(nc_field_varid(nvar/2))

        if (src_type(1)  ==  'monopole') then
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

        gllperelem = (iend - ibeg + 1)**2
        npoints = nelem * gllperelem

        call comm_elem_number(npoints, npoints_global, npoints_myfirst, npoints_mylast)

        npts_sol = nel_solid * gllperelem
        npts_flu = nel_fluid * gllperelem

        call comm_elem_number(npts_sol, npts_sol_global, npts_sol_myfirst, npts_sol_mylast)
        call comm_elem_number(npts_flu, npts_flu_global, npts_flu_myfirst, npts_flu_mylast)

        if (nstrain <= dumpstepsnap) dumpstepsnap = nstrain
    endif ! dump_wavefields

    if (mynum == 0) then
        if (verbose > 1) write(6,*) '  Producing groups for Seismograms and Snapshots'

        call check( nf90_def_grp(ncid_out, "Seismograms", ncid_recout) )
        call check( nf90_def_grp(ncid_out, "Snapshots", ncid_snapout) )
        call check( nf90_def_grp(ncid_out, "Surface", ncid_surfout) )
        call check( nf90_def_grp(ncid_out, "Mesh", ncid_meshout) )
        if (verbose > 1) write(6,*) '  Seismograms group has ID', ncid_recout
        if (verbose > 1) write(6,*) '  Snapshots group has ID', ncid_snapout
        if (verbose > 1) write(6,*) '  Surface group has ID', ncid_surfout
        if (verbose > 1) write(6,*) '  Mesh group has ID', ncid_meshout

        if (verbose > 1) write(6,*) 'Define dimensions in ''Seismograms'' group of NetCDF output file'

        if (verbose > 1) write(6,*) '  ''Seismograms'' group has ID ', ncid_recout
110     format(' Dimension ', A20, ' with length ', I8, ' and ID', I6)
        call check( nf90_def_dim(ncid_out, "seis_timesteps", nseismo, nc_times_dimid) )
        if (verbose > 1) write(6,110) "seis_timesteps", nseismo, nc_times_dimid

        call check( nf90_def_dim(ncid_out, "sim_timesteps", niter, nc_iter_dimid) )
        if (verbose > 1) write(6,110) "sim_timesteps", niter, nc_iter_dimid

        call check( nf90_def_dim(ncid_recout, "receivers", nrec, nc_rec_dimid) )
        if (verbose > 1) write(6,110) "receivers", nrec, nc_rec_dimid

        call check( nf90_def_dim(ncid_out, "components", 3, nc_comp_dimid) )
        if (verbose > 1) write(6,110) "components", 3, nc_comp_dimid

        call check( nf90_def_dim(ncid_recout, "recnamlength", 40, nc_recnam_dimid) )
        if (verbose > 1) write(6,110) "recnamlength", 40, nc_recnam_dimid

        if (verbose > 1) write(6,*) 'NetCDF dimensions defined'

        if (verbose > 1) write(6,*) 'Define variables in ''Seismograms'' group of NetCDF output file'
        call flush(6)

        call check( nf90_def_var(ncid=ncid_recout, name="displacement", xtype=NF90_FLOAT,&
                                 dimids=(/nc_times_dimid, nc_comp_dimid, nc_rec_dimid/), &
                                 !storage = NF90_CHUNKED, chunksizes=(/nseismo,3,1/), &
                                 contiguous = .false., chunksizes=(/nseismo,3,1/), &
                                 !deflate_level = deflate_level, &
                                 varid=nc_disp_varid) )

        call check( nf90_put_att(ncid_recout, nc_disp_varid, 'units', 'meters') )
        call check( nf90_put_att(ncid_recout, nc_disp_varid, '_FillValue', 0.0) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_seis", xtype=NF90_FLOAT,&
                                 dimids=(/nc_times_dimid/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_stf_seis_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="stf_iter", xtype=NF90_FLOAT,&
                                 dimids=(/nc_iter_dimid/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_stf_iter_varid) )

        call check( nf90_def_var(ncid=ncid_recout, name="time", xtype=NF90_DOUBLE,&
                                 dimids=(/nc_times_dimid/), &
                                 deflate_level = deflate_level, &
                                 varid=nc_time_varid) )

        !call check( nf90_def_var(ncid_recout, "Lat_req", NF90_FLOAT, (/nc_rec_dimid/), &
        !                         nc_latr_varid) )
        !call check( nf90_def_var(ncid_recout, "Lon", NF90_FLOAT, (/nc_rec_dimid/), &
        !                         nc_lon_varid) )
        !call check( nf90_def_var(ncid_recout, "Lat", NF90_FLOAT, (/nc_rec_dimid/), &
        !                         nc_lat_varid) )
        call check( nf90_def_var(ncid_recout, "phi", NF90_FLOAT, &
                                 [nc_rec_dimid], nc_ph_varid) )
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
            if (verbose > 1) write(6,*) 'Define variables in ''Snapshots'' group of NetCDF output file', &
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
                                      dimids = nc_snap_dimid,&
                                      varid  = nc_snaptime_varid) )

            call check( nf90_def_dim( ncid   = ncid_meshout, &
                                      name   = 'discontinuities', &
                                      len    = ndisc, &
                                      dimid  = nc_disc_dimid) )
            call check( nf90_put_att( ncid   = ncid_meshout, &
                                      varid  = NF90_GLOBAL, &
                                      name   = 'ndisc', &
                                      values = ndisc) )

            call check( nf90_def_var( ncid   = ncid_meshout,  &
                                      name   = 'mesh_S', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_s_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_Z', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_z_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout,  &
                                      name   = 'mesh_vp', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_vp_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_vs', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_vs_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout,  &
                                      name   ='mesh_rho', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_rho_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_lambda', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_lambda_varid) )
            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'mesh_mu', &
                                      xtype  = NF90_FLOAT, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_mesh_mu_varid) )

            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'model_domain', &
                                      xtype  = NF90_BYTE, &
                                      dimids = nc_pt_dimid,&
                                      varid  = nc_elem_dom_varid) )

            call check( nf90_def_var( ncid   = ncid_meshout, &
                                      name   = 'disc_depths', &
                                      xtype  = NF90_DOUBLE, &
                                      dimids = nc_disc_dimid,&
                                      varid  = nc_disc_varid) )

            do ivar=1, nvar/2 ! The big snapshot variables for the kerner.

                call check( nf90_def_var(ncid=ncid_snapout, name=trim(nc_varnamelist(ivar)), &
                                         xtype = NF90_FLOAT, &
                                         dimids = (/nc_pt_dimid, nc_snap_dimid/),&
                                         varid = nc_field_varid(ivar), &
                                         chunksizes = (/npoints_global, 1/) ))
                call check( nf90_def_var_fill(ncid=ncid_snapout, varid=nc_field_varid(ivar), &
                                              no_fill=1, fill=0) )
                if (verbose > 1) write(6,"(' Netcdf variable ', A16,' with ID ', I3, ' and length', &
                        & I8, ' produced.')") &
                          trim(nc_varnamelist(ivar)), nc_field_varid(ivar), npoints_global
            enddo

            ! Surface group in output file
            if (verbose > 1) write(6,*) 'Define variables in ''Surface'' group of NetCDF output file'
            !call check( nf90_def_dim( ncid_surfout, name="snapshots", len=nstrain, &
            !                          dimid=nc_snap_dimid) )
            call check( nf90_put_att( ncid   = ncid_surfout, &
                                      name   = 'nstrain', &
                                      varid  = NF90_GLOBAL, &
                                      values = nstrain) )
            call check( nf90_def_dim( ncid_surfout, "straincomponents", len=6, &
                                      dimid=nc_strcomp_dimid) )
            if (verbose > 1) write(6,110) "straincomponents", 6, nc_strcomp_dimid

            call check( nf90_def_dim( ncid_surfout, "surf_elems", maxind_glob, nc_surf_dimid) )
            call check( nf90_put_att( ncid   = ncid_surfout, &
                                      name   = 'nsurfelem', &
                                      varid  = NF90_GLOBAL, &
                                      values = maxind_glob) )
            if (verbose > 1) write(6,110) "surf_elems", maxind_glob, nc_surf_dimid

            call check( nf90_def_var( ncid_surfout, "elem_theta", NF90_FLOAT, &
                                      (/nc_surf_dimid /), &
                                      nc_surfelem_theta_varid) )
            call check( nf90_put_att(ncid_surfout, nc_surfelem_theta_varid, 'units', 'degrees'))

            call check( nf90_def_var( ncid_surfout, "displacement", NF90_FLOAT, &
                                      (/nc_snap_dimid, nc_comp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_disp_varid) )
            call check( nf90_put_att(ncid_surfout, nc_surfelem_disp_varid, 'units', 'meters'))

            call check( nf90_def_var( ncid_surfout, "velocity", NF90_FLOAT, &
                                      (/nc_snap_dimid, nc_comp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_velo_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_velo_varid, 'units', &
                                      'meters per second') )

            call check( nf90_def_var( ncid_surfout, "strain", NF90_FLOAT, &
                                      (/nc_snap_dimid, nc_strcomp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_strain_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_strain_varid, 'units', &
                                      ' ') )

            call check( nf90_def_var( ncid_surfout, "disp_src", NF90_FLOAT, &
                                      (/nc_snap_dimid, nc_comp_dimid, nc_surf_dimid /), &
                                      nc_surfelem_disp_src_varid) )
            call check( nf90_put_att( ncid_surfout, nc_surfelem_disp_src_varid, 'units', &
                                      'meters') )

            call check( nf90_def_var( ncid_surfout, "stf_dump", NF90_FLOAT, &
                                      [nc_snap_dimid], &
                                      nc_stf_dump_varid) )
        endif


        if (verbose > 1) write(6,'(a/)') 'NetCDF variables defined'
        ! Leave definition mode
        call check( nf90_enddef(ncid_out))

        if (verbose > 1) write(6,*) 'Writing station info into NetCDF file...'
        call check( nf90_put_var( ncid_recout, nc_th_varid,   values = rec_th) )
        call check( nf90_put_var( ncid_recout, nc_ph_varid,   values = rec_ph) )
        call check( nf90_put_var( ncid_recout, nc_thr_varid,  values = rec_th_req) )
        call check( nf90_put_var( ncid_recout, nc_proc_varid, values = rec_proc) )

        do irec=1,nrec
            call check( nf90_put_var( ncid_recout, nc_recnam_varid, start = (/irec, 1/), &
                                      count = (/1, 40/), values = (rec_names(irec))) )
        enddo

        ! Write out seismogram dump times
        time_seis = dble([ (i, i = 1, nseismo) ]) * deltat
        call check( nf90_put_var( ncid_recout, nc_time_varid, values = time_seis ) )
        if (verbose > 1) write(6,*) '...done'

        ! Write out STFs
        if (verbose > 1) write(6,*) 'Writing stf into NetCDF file...'
        call check( nf90_put_var(ncid   = ncid_recout, &
                                 varid  = nc_stf_iter_varid, &
                                 values = stf_dumpvar) )
        call check( nf90_put_var(ncid   = ncid_recout, &
                                 varid  = nc_stf_seis_varid, &
                                 values = stf_seis_dumpvar) )

        if (dump_wavefields) then
            ! Write out strain dump times
            if (verbose > 1) write(6,*) 'Writing strain dump times into NetCDF file...'
            time_strain = dble((/ (i, i = 1, nstrain) /))
            time_strain = time_strain * t_0 / strain_samp
            call check( nf90_put_var(ncid   = ncid_out, &
                                     varid  = nc_snaptime_varid, &
                                     values = time_strain ) )
            ! Write out discontinuity depths
            if (verbose > 1) write(6,*) 'Writing discontinuity depths into NetCDF file...'
            call check( nf90_put_var(ncid   = ncid_meshout, &
                                     varid  = nc_disc_varid, &
                                     values = discont) )

            ! Write out STF values at kernel dump points
            if (verbose > 1) write(6,*) 'Writing STF in strain dumps'
            call check( nf90_put_var(ncid   = ncid_surfout, &
                                     varid  = nc_stf_dump_varid, &
                                     values = stf_dump_dumpvar) )

            if (verbose > 1) write(6,*) '...done'
        endif

    endif ! (mynum == 0)


! Allocation of Dump buffer variables. Done on all procs
90  format(' Allocated ', A20, ', uses ',F10.3,'MB')
    allocate(recdumpvar(nseismo,3,num_rec))
    recdumpvar = 0.0

    if (dump_wavefields) then
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

        if (mynum == 0 .and. verbose > 1) then
            write(6,*)  'Allocating NetCDF buffer variables'
            write(6,90) 'recdumpvar', real(size(recdumpvar))/262144.
            write(6,90) 'surfdumpvar_disp', real(size(surfdumpvar_disp))/262144.
            write(6,90) 'surfdumpvar_velo', real(size(surfdumpvar_velo))/262144.
            write(6,90) 'surfdumpvar_strain', real(size(surfdumpvar_strain))/262144.
            write(6,90) 'surfdumpvar_srcdisp', real(size(surfdumpvar_srcdisp))/262144.
            write(6,90) 'oneddumpvar', real(size(oneddumpvar))/262144.
            write(6,*)
            write(6,*)  '*********************************************************************'
            write(6,*)  '**** NetCDF output file produced, buffer variables all allocated ****'
            write(6,*)  '*********************************************************************'
            write(6,*)
        endif
    endif

#endif
end subroutine nc_define_outputfile
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Character
subroutine nc_write_att_char(attribute_value, attribute_name)
    character(len=*), intent(in)    :: attribute_name, attribute_value

#ifdef unc
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_char
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Real
subroutine nc_write_att_real(attribute_value, attribute_name)
  character(len=*),  intent(in)      :: attribute_name
  real(sp), intent(in)               :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_real
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Double
subroutine nc_write_att_dble(attribute_value, attribute_name)
  character(len=*),  intent(in)      :: attribute_name
  real(dp), intent(in)               :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_dble
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Integer
subroutine nc_write_att_int(attribute_value, attribute_name)
  character(len=*),  intent(in)     :: attribute_name
  integer, intent(in)               :: attribute_value

#ifdef unc
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Open the NetCDF output file, check for variable IDs and dump meshes.
subroutine nc_finish_prepare
#ifdef unc
    use data_io,   ONLY  : datapath, lfdata, dump_wavefields
    use data_mesh, ONLY  : maxind, surfcoord, ind_first, ind_last
    integer             :: status, ivar, nmode, iproc
    integer             :: nc_mesh_s_varid, nc_mesh_z_varid
    integer             :: nc_mesh_vs_varid, nc_mesh_vp_varid
    integer             :: nc_mesh_mu_varid, nc_mesh_rho_varid
    integer             :: nc_mesh_lambda_varid

    if (mynum == 0) then
        call check(nf90_close(ncid_out))
        if (verbose > 1) then
           write(6,*) '  Root process closed netCDF file, waiting for all procs to'
           write(6,*) '  arrive here and then open it to retrieve IDs'
           if (dump_wavefields) write(6,*) '  and dump mesh coordinates.'
        endif
    endif
    call barrier

    do iproc = 0, nproc
        call barrier
        if (iproc==mynum) then
            if (verbose>1) then
               write(6,*) '  Processor ', iproc, ' opened the output file and will dump '
               write(6,*) '  his part of the mesh.'
            endif
            nmode = ior(NF90_WRITE, NF90_NETCDF4)
            call check( nf90_open(path=datapath(1:lfdata)//"/axisem_output.nc4", &
                                  mode=nmode, ncid=ncid_out) )

            print '(A,I5,A)', '   ', iproc, ': opened file'
            call check( nf90_inq_grp_ncid(ncid_out, "Seismograms", ncid_recout) )
            call check( nf90_inq_grp_ncid(ncid_out, "Surface", ncid_surfout) )
            call check( nf90_inq_grp_ncid(ncid_out, "Mesh", ncid_meshout) )
            print '(A,I5,A)', '   ', iproc, ': inquired dimension IDs'
            call getvarid( ncid_recout, "displacement", nc_disp_varid )

            if (dump_wavefields) then
                call check( nf90_inq_grp_ncid(ncid_out, "Snapshots", ncid_snapout) )
                do ivar=1, nvar/2
                    call getvarid( ncid_snapout, nc_varnamelist(ivar), &
                                nc_field_varid(ivar))
                enddo

                call getvarid( ncid_surfout, "elem_theta", &
                               nc_surfelem_theta_varid)
                call putvar_real1d(ncid   = ncid_surfout, &
                                   varid  = nc_surfelem_theta_varid, &
                                   values = surfcoord, &
                                   start  = ind_first, &
                                   count  = maxind  )


                call getvarid( ncid_surfout, "displacement", &
                               nc_surfelem_disp_varid)

                call getvarid( ncid_surfout, "velocity", &
                               nc_surfelem_velo_varid )

                call getvarid( ncid_surfout, "strain", &
                               nc_surfelem_strain_varid )

                call getvarid( ncid_surfout, "disp_src", &
                               nc_surfelem_disp_src_varid )
                print '(A,I5,A)', '   ', iproc, ': inquired variable IDs'

                ! S-Coordinate
                call getvarid( ncid_meshout, "mesh_S", nc_mesh_s_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_s_varid, &
                                    values = scoord1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Z-Coordinate
                call getvarid( ncid_meshout, "mesh_Z", nc_mesh_z_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_z_varid, &
                                    values = zcoord1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Vp
                call getvarid( ncid_meshout, "mesh_vp", nc_mesh_vp_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_vp_varid, &
                                    values = vp1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Vs
                call getvarid( ncid_meshout, "mesh_vs", nc_mesh_vs_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_vs_varid, &
                                    values = vs1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Rho
                call getvarid( ncid_meshout, "mesh_rho", nc_mesh_rho_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_rho_varid, &
                                    values = rho1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Lambda
                call getvarid( ncid_meshout, "mesh_lambda", &
                                           nc_mesh_lambda_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_lambda_varid, &
                                    values = lambda1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                ! Mu
                call getvarid( ncid_meshout, "mesh_mu", nc_mesh_mu_varid )
                call putvar_real1d( ncid_meshout, varid = nc_mesh_mu_varid, &
                                    values = mu1d, &
                                    start  = npoints_myfirst, &
                                    count  = npoints )

                print '(A,I5,A)', '   ', iproc, ': dumped mesh'

            endif !dump_wavefields
            call check( nf90_close( ncid_out))
            if (verbose > 1) write(6,"('Proc ', I3, ' dumped its mesh and is ready to rupture')") mynum
        endif !mynum==iproc
    enddo
#endif
end subroutine nc_finish_prepare
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Final dumps to netCDF file. In the moment contains only dump of
!! receiver seismograms.
subroutine nc_end_output
#ifdef unc
    use data_mesh, only: num_rec
    use data_io,   only: dump_xdmf
    integer           :: iproc

    call barrier
    do iproc=0, nproc-1
        if (iproc == mynum) then
            if (num_rec>0) then
                if (verbose > 1) write(6,"('Proc ', I3, ' will dump receiver seismograms')") mynum
                call nc_dump_rec_to_disk()
            else
                if (verbose > 1) write(6,"('Proc ', I3, ' has no receivers and just waits for the others')") mynum
            endif
            if (dump_xdmf) then
                call check(nf90_close(ncid_out_snap))
            endif
        endif

        call barrier
    enddo
    !call check( nf90_close(ncid_out) )

#endif
end subroutine nc_end_output
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_make_snapfile

    use data_mesh,    only: npoint_plot, nelem_plot
    use data_proc,    only: appmynum
    use data_io,      only: datapath, lfdata, nsnap
    use data_source,  only: src_type

#ifdef unc
    integer              :: nmode, nc_snappoint_dimid, nc_snapelem_dimid, nc_snapdim_dimid
    integer              :: nc_snaptime_dimid, nc_snapconnec_dimid
    character(len=120)   :: fname

    if (lpr .and. verbose > 1) write(6,*) '   .... preparing xdmf nc file'

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
                            name   = 'displacement',  &
                            xtype  = NF90_FLOAT,     &
                            dimids = [nc_snapdim_dimid, nc_snappoint_dimid, &
                                      nc_snaptime_dimid], &
                            chunksizes = [ndim_disp, npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_disp_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'straintrace',  &
                            xtype  = NF90_FLOAT,     &
                            dimids = [nc_snappoint_dimid, nc_snaptime_dimid], &
                            chunksizes = [npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_pwave_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'curlinplane',  &
                            xtype  = NF90_FLOAT,     &
                            dimids = [nc_snappoint_dimid, nc_snaptime_dimid], &
                            chunksizes = [npoint_plot, 1], &
                            deflate_level = deflate_level, &
                            varid  = nc_snap_swave_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'points',  &
                            xtype  = NF90_FLOAT,     &
                            dimids = [nc_coord_dimid, nc_snappoint_dimid], &
                            varid  = nc_snap_point_varid) )

    call check(nf90_def_var(ncid   = ncid_out_snap, &
                            name   = 'grid',  &
                            xtype  = NF90_INT,     &
                            dimids = [nc_snapconnec_dimid, nc_snapelem_dimid], &
                            varid  = nc_snap_grid_varid) )

    call check(nf90_enddef(ncid    = ncid_out_snap))

    if (lpr .and. verbose > 1) write(6,*) '   .... DONE'

#endif

end subroutine nc_make_snapfile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snapshot(u, straintrace, curlinplane)

    use data_mesh,      only: npoint_plot
    use data_io,        only: isnap, nsnap
    use data_source,    only: src_type

    real(kind=realkind), dimension(3,npoint_plot), intent(in)  :: u
    real(kind=realkind), dimension(1,npoint_plot), intent(in)  :: straintrace
    real(kind=realkind), dimension(1,npoint_plot), intent(in)  :: curlinplane

#ifdef unc
    if (src_type(1) == 'monopole') then
       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [1, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(1,:)) )

       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [2, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(3,:)) )
    else
       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [1, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(1,:)) )

       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [2, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(2,:)) )

       call check(nf90_put_var(ncid   = ncid_out_snap, &
                               varid  = nc_snap_disp_varid, &
                               start  = [3, 1, isnap], &
                               count  = [1, npoint_plot, 1], &
                               values = u(3,:)) )
    endif

    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_pwave_varid, &
                            start  = [1, isnap], &
                            count  = [npoint_plot, 1], &
                            values = straintrace) )

    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_swave_varid, &
                            start  = [1, isnap], &
                            count  = [npoint_plot, 1], &
                            values = curlinplane) )
#endif

end subroutine nc_dump_snapshot
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_dump_snap_points(points)

    use data_mesh, only: npoint_plot
    real(sp), dimension(2,npoint_plot), intent(in)       :: points

#ifdef unc
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

#ifdef unc
    call check(nf90_put_var(ncid   = ncid_out_snap, &
                            varid  = nc_snap_grid_varid, &
                            count  = [4, nelem_plot], &
                            values = grid) )
#endif
end subroutine nc_dump_snap_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvarid(ncid, name, varid)
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out)         :: varid
#ifdef unc
    integer                      :: status

    status = nf90_inq_varid( ncid  = ncid, &
                             name  = name, &
                             varid = varid )
    if (status/=NF90_NOERR) then
        write(6,100) mynum, trim(name), ncid
        stop
    else if (verbose>1) then
        write(6,101) trim(name), ncid, varid
        call flush(6)
    endif
100 format('ERROR: CPU ', I4, ' could not find variable: ''', A, ''' in NCID', I7)
101 format('Variable ''', A, ''' found in NCID', I7, ', has ID:', I7)
#else
    varid = 0
#endif
end subroutine getvarid
!-----------------------------------------------------------------------------------------

subroutine putvar_real1d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   real, intent(in)             :: values(:)

#ifdef unc
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   if (size(values)/=count) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname )
       write(*,100) mynum, trim(varname), varid, ncid, size(values), count
       stop
   endif

   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = [start],        &
                         count  = [count] )


   if (status/=NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)
       if (ndims/=1) then
           write(*,101) mynum, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       endif
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid,     &
                                       dimid = dimid(1), &
                                       name  = dimname,  &
                                       len   = dimsize )
       if (start + count > dimsize) then
           write(*,102) mynum, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))
           stop
       endif

       write(*,103) mynum, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status))
       stop

   else if (verbose>1) then
        write(6,200) mynum, ncid, varid
        call flush(6)
   endif

100 format('ERROR: CPU ', I4, ' could not write variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not write variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not write variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('CPU ', I5, ' wrote variable found in NCID', I7, ', with ID:', I7)
#endif
end subroutine putvar_real1d

!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    integer, intent ( in) :: status !< Error code
#ifdef unc
    if (status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        call abort()
    endif
#endif
end subroutine check
!-----------------------------------------------------------------------------------------


end module nc_routines
