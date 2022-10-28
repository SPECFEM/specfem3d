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
!> Read parameters for the general solver (i.e. NOT mesh, sources, receivers);
!! compute other parameters for the simulation;
!! write out summaries of all relevant simulation settings.
module parameters

    use global_parameters
    use data_proc
    use data_time
    use data_source
    use data_io
    use utlity
    use commun
    use coupling_mod

    implicit none

    character(len=100)  :: hostname, username, git_hash
    character(len=100)  :: compiler, compilerversion
    character(len=255)  :: fflags, cflags, ldflags
    character(len=3)    :: openmp
    character(len=16)   :: simtype

    public :: open_local_output_file, readin_parameters, read_inparam_basic_verbosity
    public :: compute_numerical_parameters, write_parameters
    private

contains

!-----------------------------------------------------------------------------------------
!> Open processor-specific output files
subroutine open_local_output_file

    if (verbose > 1) then
       open(unit=69,file='output_proc'//appmynum//'.dat')
       write(69,*)
       write(69,*)'********** This is the OUTPUT for ',procstrg, &
                  '***********************'
       write(69,*)
       write(69,*)' CONTAINS all information on the mesh, background model,'
       write(69,*)' source, receivers, precomputed matrices, and various tests '
       write(69,*)' such as surface areas, volume, valence, discontinuities,'
       write(69,*)' resolution test, axial masking, solid-fluid boundary copying,'
       write(69,*)' Lagrange interpolants, integration weights, message-passing.'
       write(69,*)
       write(69,*)'****************************************************************'
       write(69,*)
    endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Routine that reads in simulation parameters that are relevant at the
!! stage of the solver, i.e. number of time steps, integration scheme,
!! data paths, specification of wavefield dumping etc.
subroutine readin_parameters

  use data_mesh, only: do_mesh_tests

  call read_inparam_basic

  call read_inparam_advanced

  call get_runinfo

  ! now pre-set. Most of these are to be considered in the post processing stage now.
  sum_seis = .false.
  sum_fields = .false.

  !!! SB MUST BE READ SOMEWHERE
  if (dump_wavefields) then
     if (dump_type == 'coupling' .or. dump_type == 'coupling_box') coupling = .true.
  endif
  !dump_snaps_solflu = .false. !! old axisem
  !dump_type = 'fullfields' !! VM VM
!!!!!!!!!!! SB !!  dump_type = 'coupling' !! VM VM  hardcoded this for test of coupling method
!!!!!!!!!!! SB !!  coupling = .true. !! VM VM  hardcoded too this for test of coupling method
  !!! num_simul = 1 old axisem
  !! SB coupling

  ! netcdf format
  output_format = 'binary'
  if (use_netcdf) output_format='netcdf'

  call barrier
  if (lpr) then
     write(*,20)
     write(*,21) datapath, infopath, seislength_t, enforced_dt, &
                 enforced_period, trim(simtype), rec_file_type, &
                 sum_seis, sum_fields, time_scheme, seis_dt, &
                 dump_energy, dump_vtk, dump_wavefields, &
                 dump_type, ibeg, iend, jbeg, jend, strain_samp, src_dump_type, &
                 add_hetero, do_mesh_tests, output_format

20 format(/&
   08x,'///////////////////////////////////////////////////////////////',/  &
   08x,'//                                                           //',/  &
   08x,'//                  A   x   i   S   E   M                    //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//         Parallel spectral-element solution to             //',/  &
   08x,'//                                                           //',/  &
   08x,'//           3-D seismic wave propagation for                //',/  &
   08x,'//                                                           //',/  &
   08x,'//          axially symmetric background models              //',/  &
   08x,'//                                                           //',/  &
   08x,'//               in a spherical 2-D domain                   //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//  Authors : Tarje Nissen-Meyer (Oxford University)         //',/  &
   08x,'//              Martin van Driel (ETH Zurich)                //',/  &
   08x,'//                 Simon Stahler (LMU Munich)                //',/  &
   08x,'//                Kasra Hosseini (LMU Munich)                //',/  &
   08x,'//               Stefanie Hempel (University of Muenster)    //',/  &
   08x,'//            Alexandre Fournier (IPG Paris)                 //',/  &
   08x,'//                   Tony Dahlen (Princeton University)      //',/  &
   08x,'//                                                           //',/  &
   08x,'//   Contact:     info@axisem.info                           //',/  &
   08x,'//   Information: www.axisem.info                            //',/  &
   08x,'//                                                           //')

21 format(/&
   08x,'//                                                           //',/  &
   08x,'//     If you are publishing results obtained with this      //',/  &
   08x,'//          code, please cite this paper:                    //',/  &
   08x,'//                                                           //',/  &
   08x,'// (1) T. Nissen-Meyer, M. van Driel, S. C. Staehler,        //',/  &
   08x,'//     K. Hosseini, S. Hempel, L. Auer, A. Colombi           //',/  &
   08x,'//     and A. Fournier:                                      //',/  &
   08x,'//     "AxiSEM: broadband 3-D seismic wavefields in          //',/  &
   08x,'//              axisymmetric media"                          //',/  &
   08x,'//     Solid Earth, 5, 425-445, 2014                         //',/  &
   08x,'//     doi:10.5194/se-5-425-2014                             //',/  &
   08x,'//                                                           //',/  &
   08x,'//       Comprehensive description of the underlying         //',/  &
   08x,'//           numerical analysis can be found in:             //',/  &
   08x,'//                                                           //',/  &
   08x,'// (2) Tarje Nissen-Meyer, F. A. Dahlen, A Fournier (2007)   //',/  &
   08x,'//     "Spherical-earth Frechet sensitivity kernels"         //',/  &
   08x,'//     Geophysical Journal International 168(3),1051-1066.   //',/  &
   08x,'//     doi:10.1111/j.1365-246X.2006.03123.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (3) Tarje Nissen-Meyer, A Fournier, F. A. Dahlen (2007)   //',/  &
   08x,'//     "A two-dimensional spectral-element method for        //',/  &
   08x,'//        computing spherical-earth seismograms -            //',/  &
   08x,'//        I. Moment-tensor source"                           //',/  &
   08x,'//     Geophysical Journal International 168(3), 1067-1092.  //',/  &
   08x,'//     doi:10.1111/j.1365-246X.2006.03121.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (4) Tarje Nissen-Meyer, A Fournier, F. A. Dahlen (2007)   //',/  &
   08x,'//     "A two-dimensional spectral-element method for        //',/  &
   08x,'//        computing spherical-earth seismograms -            //',/  &
   08x,'//        II.  Waves in solid-fluid media"                   //',/  &
   08x,'//     Geophysical Journal International 174(3), 873-888.    //',/  &
   08x,'//     doi:10.1111/j.1365-246X.2008.03813.x                  //',/  &
   08x,'//                                                           //',/  &
   08x,'// (5) Martin van Driel and Tarje Nissen-Meyer (2014)        //',/  &
   08x,'//     "Seismic wave propagation in fully anisotropic        //',/  &
   08x,'//        axisymmetric media"                                //',/  &
   08x,'//      Geophysical Journal International 199 (2): 880-893.  //',/  &
   08x,'//      doi: 10.1093/gji/ggu269                              //',/  &
   08x,'//                                                           //',/  &
   08x,'// (6) Martin van Driel and Tarje Nissen-Meyer (2014)        //',/  &
   08x,'//     "Optimized visco-elastic wave propagation for         //',/  &
   08x,'//        weakly dissipative media"                          //',/  &
   08x,'//      Geophysical Journal International 199 (2): 1078-1093.//',/  &
   08x,'//      doi: 10.1093/gji/ggu314                              //',/  &
   08x,'//                                                           //',/  &
   08x,'//                                                           //',/  &
   08x,'//  January 2015: version 1.2                                //',/  &
   08x,'//                                                           //',/  &
   08x,'///////////////////////////////////////////////////////////////',// &
   08x,' ==  ==  ==  ==  ==  == =  I N P U T    P A R A M E T E R S ==  ==  ==  ==  ==  ==  == =',/  &
   12x,'Data I/O path:                      ',a20,/                         &
   12x,'Info I/O path:                      ',a20,/                         &
   12x,'Simulation length [s]:              ',f9.3,/                        &
   12x,'Enforced time step [s]:             ',f7.3,/                        &
   12x,'Enforced source period [s]:         ',f7.3,/                        &
   12x,'Simulation type:                    ',a,/                          &
   12x,'Receiver file type:                 ',a8,/                          &
   12x,'Sum seismograms?                    ',l2,/                          &
   12x,'Sum wavefields?                     ',l2,/                          &
   12x,'Time extrapolation scheme:          ',a8,/                          &
   12x,'Seismogram sampling rate [s]:       ',f7.3,/                        &
   12x,'Dump kin./pot. energy?              ',l2,/                          &
   12x,'Dump global snaps?                  ',l2,/                          &
   12x,'Dump strain?                        ',l2,/                          &
   12x,'Wavefield dumping type:             ',a12,/                         &
   12x,'First GLL to save in strains:       ',i2,/                          &
   12x,'Last GLL to save in strains:        ',i2,/                          &
   12x,'First GLL to save in strains:       ',i2,/                          &
   12x,'Last GLL to save in strains:        ',i2,/                          &
   12x,'Samples per period for strains:     ',f7.3,/                        &
   12x,'Source dumping type:                ',a4,/                          &
   12x,'Add heterogeneous region?           ',l2,/                          &
   12x,'Perform extensive mesh tests?       ',l2,/                          &
   12x,'Output format (seism., wavefields): ',a6,/                          &
   08x,' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == =')

     if (verbose > 1) write(*,'(a/a)') &
           'Processor-specific output is written to: output_proc < PROC ID > .dat', &
           'All potential error messages will appear here...'
  endif !lpr

  call check_basic_parameters

  ! Need to decide here since this boolean is needed in def_precomp_terms
  need_fluid_displ = .false.
  if (dump_vtk .or. dump_xdmf .or. dump_energy .or. dump_wavefields .and. &
        (dump_type == 'fullfields' .or. dump_type == 'displ_only' &
        .or. dump_type == 'strain_only' .or. dump_type == 'coupling' .or. dump_type == 'coupling_box'  )) then !!SB
     ! Need to add this for each new type of wavefield dumping method that
     ! requires the fluid displacement/velocities
     need_fluid_displ = .true.
  endif


  if (lpr .and. verbose > 1) write(*,*)'     small value is:',smallval

end subroutine readin_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read file inparam_basic
subroutine read_inparam_basic
  use data_mesh, only: meshname
  use commun, only: broadcast_int, broadcast_log, broadcast_char, broadcast_dble

  integer             :: iinparam_basic=500, ioerr
  character(len=256)  :: line
  character(len=256)  :: keyword, keyvalue

  ! Default values
  seislength_t = 1800.0
  add_hetero = .false.
  rec_file_type = 'stations'
  do_anel = .false.
  dump_snaps_glob = .false.

  ! These values have to be set
  simtype = 'undefined'
  meshname = 'undefined'

  ! only rank 0 reads the file
  if (mynum == 0) then
     keyword = ' '
     keyvalue = ' '

     if (verbose > 1) write(*,'(A)', advance='no') '    Reading inparam_basic...'
     open(unit=iinparam_basic, file='inparam_basic', status='old', action='read',  iostat=ioerr)
     if (ioerr /= 0) stop 'Check input file ''inparam_basic''! Is it still there?'

     do
       read(iinparam_basic, fmt='(a256)', iostat=ioerr) line
       if (ioerr < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

       read(line,*) keyword, keyvalue

       parameter_to_read : select case(trim(keyword))

       case('SIMULATION_TYPE')
         read(keyvalue, *) simtype

         if (to_lower(trim(simtype)) /= 'single' .and. to_lower(trim(simtype)) &
 /= 'force' .and. to_lower(trim(simtype)) /= 'moment') then
            stop 'SIMULATION_TYPE in inparam_basic has invalid value'
         endif

       case('RECFILE_TYPE')
          rec_file_type = to_lower(keyvalue)

       case('SEISMOGRAM_LENGTH')
          read(keyvalue, *) seislength_t

       case('MESHNAME')
          meshname = keyvalue

       case('LAT_HETEROGENEITY')
          read(keyvalue, *) add_hetero

       case('ATTENUATION')
          read(keyvalue,*) do_anel

       case('SAVE_SNAPSHOTS')
          read(keyvalue, *) dump_snaps_glob

       end select parameter_to_read

     enddo

     close(iinparam_basic)
     if (verbose > 1) print *, 'done'
  endif

  ! broadcast values to other processors
  call broadcast_char(simtype, 0)
  call broadcast_char(rec_file_type, 0)
  call broadcast_dble(seislength_t, 0)
  call broadcast_char(meshname, 0)
  call broadcast_log(add_hetero, 0)
  call broadcast_log(do_anel, 0)
  call broadcast_log(dump_snaps_glob, 0)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> get verbosity in the very beginning
subroutine read_inparam_basic_verbosity

  use commun, only: broadcast_int
  integer             :: iinparam_advanced=500, ioerr
  character(len=256)  :: line
  character(len=256)  :: keyword, keyvalue

  ! Default value
  verbose = 1

  ! only rank 0 reads the file
  if (mynum == 0) then

     keyword = ' '
     keyvalue = ' '

     open(unit=iinparam_advanced, file='inparam_advanced', status='old', action='read',  iostat=ioerr)
     if (ioerr /= 0) stop 'Check input file ''inparam_advanced''! Is it still there?'
     do
         read(iinparam_advanced, fmt='(a256)', iostat=ioerr) line
         if (ioerr < 0) exit
         if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

         read(line,*) keyword, keyvalue

         if (trim(keyword) == 'VERBOSITY') read(keyvalue, *) verbose
     enddo
     close(iinparam_advanced)
  endif

  ! broadcast verbosity to other processors
  call broadcast_int(verbose, 0)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read file inparam_advanced
subroutine read_inparam_advanced

  use nc_routines, only: nc_dumpbuffersize, nc_chunk_time_traces
  use data_mesh, only: naxel, meshname, do_mesh_tests
  use commun, only: broadcast_int, broadcast_log, broadcast_char, broadcast_dble

  integer              :: iinparam_advanced=500, ioerr
  integer              :: npol_max = 12
  character(len=256)   :: line
  character(len=256)   :: keyword, keyvalue

  ! Default values
  seis_dt = 0.0
  enforced_period = 0.0
  enforced_dt = 0.0
  stf_type = 'gauss_0'
  time_scheme = 'newmark2'

  datapath = './Data'
  infopath = './Info'

  diagfiles = .false.
  do_mesh_tests = .false.

  dump_wavefields = .false.
  dump_type = 'fullfields'
  strain_samp = 8
  src_dump_type = 'mask'
  ibeg = 1
  iend = 3
  jbeg = 1
  jend = 3

  kwf_rmin = 0
  kwf_rmax = 7d6
  kwf_thetamin = 0
  kwf_thetamax = pi

  dump_energy = .false.

  deflate_level = 5
  snap_dt = 20.
  dump_vtk = .false.
  dump_xdmf = dump_snaps_glob
  use_netcdf = .false.
  checkpointing = .false.
  nc_dumpbuffersize = 128
  nc_chunk_time_traces = .false.

  ! xdmf stuff
  i_n_xdmf = -1
  j_n_xdmf = -1
  allocate(i_arr_xdmf(1:npol_max+1))
  allocate(j_arr_xdmf(1:npol_max+1))
  i_arr_xdmf = -1
  j_arr_xdmf = -1
  xdmf_rmin = 0
  xdmf_rmax = 7d6
  xdmf_thetamin = 0
  xdmf_thetamax = pi


  keyword = ' '
  keyvalue = ' '

  if (mynum == 0) then
     if (verbose > 1) write(*, '(A)', advance='no') '    Reading inparam_advanced...'
     open(unit=iinparam_advanced, file='inparam_advanced', status='old', action='read', iostat=ioerr)
     if (ioerr /= 0) stop 'Check input file ''inparam_advanced''! Is it still there?'

     do
         read(iinparam_advanced, fmt='(a256)', iostat=ioerr) line
         if (ioerr < 0) exit
         if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
         read(line,*) keyword, keyvalue

         parameter_to_read : select case(trim(keyword))

         case('SAMPLING_PERIOD')
             read(keyvalue, *) seis_dt

         case('SOURCE_PERIOD')
             read(keyvalue, *) enforced_period

         case('SOURCE_function')
             read(keyvalue, *) stf_type
             stf_type = to_lower(stf_type)

         case('TIME_STEP')
             read(keyvalue, *) enforced_dt

         case('TIME_SCHEME')
             read(keyvalue,*) time_scheme
             time_scheme = to_lower(time_scheme)

         case('DATA_DIR')
             !datapath = keyvalue

         case('INFO_DIR')
             infopath = keyvalue

         case('RECIP_KH_INTEGRAL')
             read(keyvalue, *) storage_for_recip_KH_integral

         case('DIAGNOSTIC_FILE_OUTPUT')
             read(keyvalue,*) diagfiles

         case('MESH_TEST')
             read(keyvalue,*) do_mesh_tests

         case('KERNEL_WAVEFIELDS')
             read(keyvalue,*) dump_wavefields

         case('KERNEL_DUMPTYPE')
             read(keyvalue,*) dump_type
             dump_type = to_lower(dump_type)
             if (trim(dump_type) /= 'fullfields' .and. &
                 trim(dump_type) /= 'displ_only' .and. &
                 trim(dump_type) /= 'strain_only' .and. &
                 trim(dump_type) /= 'coupling' .and. &                  !! SB
                 trim(dump_type) /= 'coupling_box' .and. &        !! SB
                 trim(dump_type) /= 'displ_velo') then
                   write(*,*) dump_type
                   stop 'ERROR: invalid value for KERNEL_DUMPTYPE!'
             endif

         case('KERNEL_SPP')
             read(keyvalue,*) strain_samp

         case('KERNEL_SOURCE')
             read(keyvalue,*) src_dump_type
             src_dump_type = to_lower(src_dump_type)

         case('KERNEL_IBEG')
             read(keyvalue,*) ibeg

         case('KERNEL_IEND')
             read(keyvalue,*) iend

         case('KERNEL_JBEG')
             read(keyvalue,*) jbeg

         case('KERNEL_JEND')
             read(keyvalue,*) jend

         case('KERNEL_RMIN')
             read(keyvalue, *) kwf_rmin
             kwf_rmin = kwf_rmin * 1000

         case('KERNEL_RMAX')
             read(keyvalue, *) kwf_rmax
             kwf_rmax = kwf_rmax * 1000

         case('KERNEL_COLAT_MIN')
             read(keyvalue, *) kwf_thetamin
             kwf_thetamin = kwf_thetamin * pi / 180.

         case('KERNEL_COLAT_MAX')
             read(keyvalue, *) kwf_thetamax
             kwf_thetamax = kwf_thetamax * pi / 180.

         case('SAVE_ENERGY')
             read(keyvalue,*) dump_energy

         case('USE_NETCDF')
             read(keyvalue, *) use_netcdf

         case('CHECKPOINTING')
             read(keyvalue, *) checkpointing

         case('NETCDF_DUMP_BUFFER')
             read(keyvalue, *) nc_dumpbuffersize

         case('CHUNK_TIME_TRACES')
             read(keyvalue, *) nc_chunk_time_traces

         case('DEFLATE_LEVEL')
             read(keyvalue,*) deflate_level

         case('SNAPSHOT_DT')
             read(keyvalue, *) snap_dt

         case('SNAPSHOTS_FORMAT')
             if (dump_snaps_glob) then
               select case (to_lower(trim(keyvalue)))
                 case('xdmf')
                     dump_xdmf = .true.
                     dump_vtk = .false.
                 case('vtk')
                     dump_xdmf = .false.
                     dump_vtk = .true.
                 case('both')
                     dump_xdmf = .true.
                     dump_vtk = .true.
                 case default
                     stop 'invalid value for snapshots format!'
                 end select
             endif

         case('XDMF_RMIN')
             read(keyvalue, *) xdmf_rmin
             xdmf_rmin = xdmf_rmin * 1000

         case('XDMF_RMAX')
             read(keyvalue, *) xdmf_rmax
             xdmf_rmax = xdmf_rmax * 1000

         case('XDMF_COLAT_MIN')
             read(keyvalue, *) xdmf_thetamin
             xdmf_thetamin = xdmf_thetamin * pi / 180.

         case('XDMF_COLAT_MAX')
             read(keyvalue, *) xdmf_thetamax
             xdmf_thetamax = xdmf_thetamax * pi / 180.

         case('XDMF_GLL_I')
             read(keyvalue, *) i_n_xdmf
             read(keyvalue, *) i_n_xdmf, i_arr_xdmf(1:i_n_xdmf)

         case('XDMF_GLL_J')
             read(keyvalue, *) j_n_xdmf
             read(keyvalue, *) j_n_xdmf, j_arr_xdmf(1:j_n_xdmf)

         end select parameter_to_read
     enddo
  endif


  call broadcast_dble(seis_dt, 0)
  call broadcast_dble(enforced_period, 0)
  call broadcast_char(stf_type, 0)
  call broadcast_dble(enforced_dt, 0)

  call broadcast_char(time_scheme, 0)
  call broadcast_char(datapath, 0)
  call broadcast_char(infopath, 0)

  call broadcast_log(diagfiles, 0)
  call broadcast_log(do_mesh_tests, 0)
  call broadcast_log(dump_wavefields, 0)
  call broadcast_char(dump_type, 0)

  call broadcast_dble(strain_samp, 0)
  call broadcast_char(src_dump_type, 0)

  call broadcast_dble(kwf_rmin, 0)
  call broadcast_dble(kwf_rmax, 0)
  call broadcast_dble(kwf_thetamin, 0)
  call broadcast_dble(kwf_thetamax, 0)

  call broadcast_int(ibeg, 0)
  call broadcast_int(iend, 0)
  call broadcast_int(jbeg, 0)
  call broadcast_int(jend, 0)

  call broadcast_log(dump_energy, 0)

  call broadcast_int(deflate_level, 0)
  call broadcast_dble(snap_dt, 0)

  call broadcast_log(dump_energy, 0)
  call broadcast_log(dump_snaps_glob, 0)
  call broadcast_log(dump_xdmf, 0)
  call broadcast_log(dump_vtk, 0)
  call broadcast_log(use_netcdf, 0)
  call broadcast_log(checkpointing, 0)
  call broadcast_log(storage_for_recip_KH_integral, 0)
  call broadcast_int(nc_dumpbuffersize, 0)
  call broadcast_log(nc_chunk_time_traces, 0)

  call broadcast_dble(xdmf_rmin, 0)
  call broadcast_dble(xdmf_rmax, 0)
  call broadcast_dble(xdmf_thetamin, 0)
  call broadcast_dble(xdmf_thetamax, 0)

  call broadcast_int(i_n_xdmf, 0)
  call broadcast_int(j_n_xdmf, 0)

  call broadcast_int_arr(i_arr_xdmf, 0)
  call broadcast_int_arr(j_arr_xdmf, 0)

  lfdata = index(datapath,' ') - 1
  lfinfo = index(infopath,' ') - 1
  if (lpr .and. verbose > 1) print *, 'done'

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Getting information like code revision, username and hostname
subroutine get_runinfo
  integer  :: iget_runinfo = 500, ioerr

  hostname = 'UNKNOWN'
  username = 'UNKNOWN'
  git_hash = 'UNKNOWN'

  if (lpr .and. verbose > 1) write(*, '(A)', advance='no') '    Reading runinfo... '
  open(unit=iget_runinfo, file='runinfo', status='old', action='read',  iostat=ioerr)
  if (ioerr /= 0) then
     if (lpr .and. verbose > 1) &
        write(*,*) 'No file ''runinfo'' found, continuing without.'
  else
     read(iget_runinfo,*) git_hash
     read(iget_runinfo,*) username
     read(iget_runinfo,*) hostname
     read(iget_runinfo,'(A)') fflags
     read(iget_runinfo,'(A)') cflags
     read(iget_runinfo,'(A)') ldflags
     if (lpr .and. verbose > 1) print *, 'done'
  endif
  close(iget_runinfo)
  openmp = 'no'
  !$ openmp = 'yes'
  !>@TODO: Include support for more compilers. Also, Fortran2008 has the
  !!       intrinsic procedures: COMPILER_OPTIONS and COMPILER_VERSION.
  !!       Alas, they are only supported for gcc>4.6 and some ifort>13.
  !!       At a later stage, this whole module here could be simplified
  !!       using them.
  !!
  !!       some more compiler variables:
  !!       https://github.com/adobe-flash/crossbridge/blob/master/cmake-2.8.10.1/Modules/CMakeFortranCompilerId.F.in

#if defined(__gfortran__)
   compiler = 'gfortran'
#define gfortranversion __VERSION__
   compilerversion = gfortranversion
#undef gfortranversion
#endif

#if defined(__INTEL_COMPILER)
   compiler = 'ifort'
#define ifortversion __INTEL_COMPILER
   write(compilerversion, *) ifortversion
#undef ifortversion
#endif

#if defined(_CRAYFTN)
   compiler = 'crayFortran'
#define crayfortversion _CRAYFTN
   write(compilerversion, *) crayfortversion
#undef crayfortversion
#endif

#if defined(__PGI)
   compiler = 'pg Fortran'
#define pgfortversion __PGI
   write(compilerversion, *) pgfortversion
#undef pgfortversion
#endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Checking the consistency of some of the input parameters
subroutine check_basic_parameters
  use data_mesh
  character(len=1024) :: errmsg

  errmsg = 'SIMULATION_TYPE is not defined in input file inparam_basic'
  call pcheck(trim(simtype) == 'undefined', errmsg)

  errmsg = 'MESHNAME is not defined in input file inparam_basic'
  call pcheck(trim(meshname) == 'undefined', errmsg)

#ifndef enable_netcdf
  errmsg = 'trying to use netcdf IO but axisem was compiled without netcdf'
  call pcheck(use_netcdf, errmsg)
#endif

  errmsg = 'Analytical source wavefield dump not implemented YET!'
  call pcheck(src_dump_type == 'anal', errmsg)

  write(errmsg ,'(a,a,i4,a,i4,a,a,i4,a)') &
        'PROBLEM with REAL data kind!\n', &
        '... can only handle real kinds', sp, ' or ', dp, '.\n', &
        'real kind here:', realkind, &
        '\nchange parameter realkind in global_parameters.f90'
  call pcheck((realkind /= sp .and. realkind /= dp), errmsg)

!!!!!!!!!!!!!!!!!!!  OLD VERSION OF AXISEM SB
! errmsg = "!!!!!! NOT GOING ANY FURTHER !!!!!!\n" &
!        // "  It's just too much to save 10 frames of strain & velocity\n" &
!        // "  per source period! Choose something reasonable."
!  call pcheck(strain_samp > 30, errmsg) !! VM VM change 15 to 30
!!!!!!!!!!!!!!!!!!! OLD VERSION OF AXISEM

  errmsg = "Need indices for GLL points to dump xdmf. Set XDMF_GLL_* in inparam_advanced"
  call pcheck(dump_xdmf .and. (i_n_xdmf == -1 .or. j_n_xdmf == -1 ), errmsg)

  if (enforced_dt > zero) then
     if (lpr) then
        write(*,*)
        write(*,14) 'maximal time step', enforced_dt
     endif
  endif

  if (enforced_period > zero) then
     if (lpr .and. verbose > 1) then
        write(*,*)
        write(*,14) 'min. source period', enforced_period
     endif
  endif

14 format('  WARNING: Overriding',a19,' with:',f8.3,' seconds')

7 format(04x,a62)

  if (verbose > 1) then
      if (realkind == 4) then
          if (lpr) then
              write(*,7)
              write(*,7)'44444444444444444444444444444444444444444444444444444444444444'
              write(*,7)'444   Running the solver time loop with SINGLE PRECISION   444'
              write(*,7)'44444444444444444444444444444444444444444444444444444444444444'
              write(*,7)
          endif
      else if (realkind == 8) then
          if (lpr) then
              write(*,7)
              write(*,7)'88888888888888888888888888888888888888888888888888888888888888'
              write(*,7)'888   Running the solver time loop with DOUBLE PRECISION   888'
              write(*,7)'88888888888888888888888888888888888888888888888888888888888888'
              write(*,7)
          endif
      endif
  endif

end subroutine check_basic_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Compute numerical parameters, like time step, snapshot frequency
subroutine compute_numerical_parameters
  use attenuation, only: dump_memory_vars
  use data_mesh

  real(kind=dp)         :: s,z,r,theta,s_max,dshift
  real(kind=dp)         :: dsaxis(0:npol-1,0:npol), dzaxis(0:npol-1)
  real(kind=dp)         :: minds(nelem),maxds(nelem),mindz(nelem),maxdz(nelem)
  integer               :: ielem,ipol,jpol,i
  logical               :: found_shift
  character(len=1024)   :: errmsg

  if (lpr .and. verbose > 1) write(*,'(/,a)') '  Computing numerical parameters...'

  ! Overwrite time step or source period if demanded by input
  if (enforced_dt > zero) then
     if (time_scheme == 'newmark2' .and. enforced_dt > deltat .or. &
         time_scheme == 'symplec4' .and. enforced_dt > 1.5 * deltat .or. &
         time_scheme == 'SS_35o10' .and. enforced_dt > 3.0 * deltat  ) then
        write(errmsg,*) &
              'PROBLEM: Time step larger than allowed by mesh!\n', &
              'Chosen value (in inparam file) [s] :', enforced_dt, '\n', &
              'Maximal time step for this mesh [s]:', deltat, '\n', &
              'Change time step in input file to lower than this max\n', &
              'or to zero to use precalculated (recommended)'
        call pcheck(.true., errmsg)
     else
        if (lpr .and. verbose > 1) then
           write(*,'(/,a)')'    WARNING: Time step smaller than necessary by mesh!'
           write(*,20) enforced_dt, deltat
           write(*,19) 100. - enforced_dt / deltat * 100.
        endif
        deltat = enforced_dt
     endif
  else
     ! mesher suggestion is based on newmark scheme, the higher order schemes are stable
     ! for longer time steps
     select case (time_scheme)
     case ('newmark')
        deltat = deltat
     case ('symplec4')
        deltat = deltat * 1.5
     case ('SS_35o10')
        deltat = deltat * 3.0
     end select
     if (lpr .and. verbose > 1) then
        write(*,'(/,a)')'    Using time step precalculated by the mesher:',deltat
     endif
  endif
20 format('     Chosen/maximal time step [s]:',2(f7.3))
19 format('     ...lengthens this simulation by',f6.2,' percent!')

  ! source period
  if (enforced_period > zero &
       .and. (trim(stf_type) /= 'dirac_0' .or. trim(stf_type) /= 'quheavi') ) then
     !@TODO: does this catch a period shorter then mesh period but with e.g.
     !       Gaussian stf? (MvD)
     if (enforced_period < period) then
        if (lpr) then
           write(*,*)
           write(*,*) '    ERROR: Period smaller than necessary by mesh!'
           write(*,*) '    A pulse of this (short) chosen half width will produce numerical '
           write(*,*) '    noise on this (coarse) mesh'
           write(*,21)'   Chosen value (in inparam file):',enforced_period
           write(*,21)'   Minimal period for this mesh  :',period
           write(*,*) '    Change period in input file to larger than this min.'
           write(*,*) '    or to zero to use precalculated (recommended)'
           stop
        endif
     else
        if (lpr) then
           write(*,*)
           write(*,*) '    WARNING: Using larger period than necessary by mesh!'
           write(*,23) enforced_period, period
        endif
        t_0 = enforced_period
     endif
  else
     if (trim(stf_type) /= 'dirac_0' .or. trim(stf_type) /= 'quheavi') then
        !@TODO: so t_0 = period anyway, but some useless info to the output file in case
        !       of a dirac? (MvD)
        if (lpr) then
           write(*,*)
           write(*,*)'    Using period of the mesh:', period
        endif
        t_0 = period
     else
        t_0 = period ! Just for consistency
     endif
  endif

21 format(a36,f8.3,' s')
23 format('     Chosen/minimal period   [s]:',2(f7.3))


  ! Compute number of iterations in time loop
  niter = ceiling((seislength_t + smallval_dble) / deltat)
  if (lpr) then
     write(*,*)
     write(*,22) '    desired simulation length  :', seislength_t, ' seconds'
     write(*,22) '    offered simulation length  :', niter * deltat, ' seconds'
     write(*,11) '    number time loop iterations:', niter
  endif

  ! Compute seismogram sampling rate in time steps
  if (lpr) then
     write(*,*)
     write(*,22) '    desired seismogram sampling:', seis_dt, ' seconds'
  endif

  if (seis_dt > 0.0 .and. seis_dt >= deltat ) then
     seis_it = floor((seis_dt + smallval_dble) / deltat)
     seis_dt = deltat * seis_it
  else if (seis_dt > 0.0 .and. seis_dt < deltat) then
     write(errmsg,*) 'seismogram sampling cannot be smaller than time step... \n', &
                     ' seismogram sampling:  ', seis_dt, &
                     '\ nsimulation time step: ', deltat
     call pcheck(.true.,errmsg)
  else
     seis_it = 1
     seis_dt = deltat
  endif
  deltat_coarse = seis_dt

  nseismo = floor(real(niter) / real(seis_it)) + 1

  ! Frequency of checkpointing. Hardcoded to every 5% of runtime
  check_it = niter / 20

  if (lpr) then
     write(*,22)'    offered seismogram sampling:', deltat * seis_it, ' seconds'
     write(*,13)'    ...that is, every          :', seis_it, ' timesteps'
     write(*,11)'    number of samples          :', nseismo
  endif
22 format(a33,f9.2,a10)

  ! snapshot output, convert from interval given in seconds to
  ! incremental time steps
  if (dump_vtk .or. dump_xdmf .or. dump_memory_vars) then
     snap_it = floor(snap_dt / deltat)
     open(unit=2900+mynum, file=datapath(1:lfdata)//'/snap_info.dat'//appmynum)
     nsnap = floor(real(niter) / real(snap_it)) + 1

     write(2900+mynum,*) nsnap
     do ielem=1, nsnap
        write(2900+mynum,*) real(ielem) * snap_dt, ielem * snap_it
     enddo
     close(2900+mynum)

     if (lpr) then
        write(*,*)
        write(*,11)'    Number of snapshots        :',nsnap
        write(*,12)'    ...approximately every     :',snap_dt,'seconds'
        write(*,13)'    ...that is, every          :',snap_it,'timesteps'
     endif
11   format(a33,i8)
12   format(a33,f8.2,a10)
13   format(a33,i8,a10)

  endif

  ! Source time function
  if (lpr) then
     print *
     write(*,*)'  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-'
     write(*,*)'  SOURCE TIME function: ',trim(stf_type)
     ! TODO: discrete_dirac and discrete choice are not set when arriving here (MvD)
     if (discrete_dirac) write(*,*)'    discrete Dirac type: ', discrete_choice
  endif

  period_vs_discrete_halfwidth = 8.
  sampling_per_a = 1
  decay = 3.5d0

  if (.not. dump_wavefields &
      .and. ( trim(stf_type) == 'dirac_0' .or. trim(stf_type) == 'queavi') &
      .and. deltat_coarse > 1.9 * deltat ) then
     period_vs_discrete_halfwidth = period / (2. * deltat_coarse)
     if (period_vs_discrete_halfwidth < 15.) period_vs_discrete_halfwidth=15.
     if (lpr) then
        write(*,*)'    No wavefield dump, but discrete Dirac due to seismogram downsampling'
        write(*,*)'    Set discrete Dirac half width to [s]:', period / period_vs_discrete_halfwidth
        write(*,*)'    ... i.e. this part in the mesh period:', period_vs_discrete_halfwidth
     endif
     !@TODO: this is not very transparent to the user...
  endif

  if (trim(stf_type) == 'dirac_0' .or. trim(stf_type) == 'quheavi') then
     discrete_dirac = .true.

     if (dump_wavefields .or. deltat_coarse > 1.9 * deltat) then
        discrete_choice = 'gaussi'
     else
        discrete_choice = '1dirac'
     endif

     if (20.*seis_dt > period ) then
        if (lpr) then
           write(*,*)'   +++++++++++++++++++ W A R N I N G +++++++++++++++++++++ '
           write(*,*)'   The sampling period of seismograms is quite coarse given the'
           write(*,*)'   Dirac delta source time function. We suggest to use at least'
           write(*,*)'   20 points per period to ensure accurate results, i.e. seis_dt <= ',period/20.
           write(*,*)'   +++++++++++++++ E N D  o f  W A R N I N G +++++++++++++ '
        endif
     endif
  else
     discrete_dirac = .false.
  endif

  if (dump_wavefields) then
     ! define coarse time step for strain wavefield dumps
     !@TODO: why read strain_samp from inparam advanced and then overwrite
     !       here?? (MvD)
     if (discrete_dirac) strain_samp = ceiling(period_vs_discrete_halfwidth * sampling_per_a)
     strain_it     = floor(period / strain_samp / deltat)
     deltat_coarse = max(deltat_coarse, deltat * strain_it)
     strain_samp   = period / deltat_coarse
     strain_it     = floor(t_0 / strain_samp / deltat)
     deltat_coarse = deltat * strain_it
     !@TODO: This is a mess, but this way it is at least consistent between output and the actual value
     !       WTF? was there some endif lost or so? What's happening here??? (MvD)
     if (lpr) then
       write(*,*)'   dumping wavefields at sampling rate and deltat:', strain_samp, deltat_coarse
     endif
  else
     strain_it = seis_it
  endif

  if (lpr) write(*,*)'   coarsest dump every', strain_it, 'th time step, dt:', deltat_coarse

  if (discrete_dirac) then
     discrete_dirac_halfwidth = period / period_vs_discrete_halfwidth
     t_0 = discrete_dirac_halfwidth
     if (lpr) then
        write(*,*)
        write(*,*)'    DISCRETE DIRAC DEFINITIONS:'
        write(*,*)'    Period discrete Dirac, mesh, simul. [s]:', &
                       real(discrete_dirac_halfwidth),real(period),real(t_0)
        write(*,*)'    period mesh/period discrete Dirac:',real(period_vs_discrete_halfwidth)
        write(*,*)"    deltat SEM, seis, coarse [s]:",real(deltat),real(seis_dt),real(deltat_coarse)
        write(*,*)'    # seismogram points per mesh,Dirac period:', &
                       real(period/seis_dt),real(discrete_dirac_halfwidth/seis_dt)
        if (dump_wavefields) then
           write(*,*)'    # coarse points per mesh, Dirac period (int) :',strain_samp,sampling_per_a
           write(*,*)'    # coarse points per mesh, Dirac period (calc):', &
                real(period/deltat_coarse),real(discrete_dirac_halfwidth/deltat_coarse)
        endif
        write(*,*)'    # SEM points per mesh, Dirac period:', &
                       real(period/deltat),real(discrete_dirac_halfwidth/deltat)
     endif

     ! parameters for source-time function time shift
     found_shift = .false.
     do i=1, ceiling(4. * discrete_dirac_halfwidth / deltat)
        dshift = deltat * ceiling(4. * discrete_dirac_halfwidth / deltat) + real(i) * deltat
        if (.not. found_shift &
             .and. abs(nint(dshift / deltat_coarse) - dshift / deltat_coarse) < 0.01 * deltat &
             .and. abs(nint(dshift / deltat) - dshift / deltat) < 0.01 * deltat &
             .and. abs(nint(dshift / seis_dt) - dshift / seis_dt) < 0.01 * deltat) then
           shift_fact_discrete_dirac = deltat_coarse * ceiling(dshift / deltat_coarse)
           found_shift = .true.
        endif
     enddo
     shift_fact = shift_fact_discrete_dirac

  else ! smooth source time function, e.g. Gauss
     shift_fact = deltat_coarse * ceiling(1.5 * t_0 / deltat_coarse)
  endif

  if (lpr) then
     print *
     write(*,*)'  SHIFT FACTOR of source time function [s]:', shift_fact
     write(*,*)'   # SEM, seis, coarse points per shift factor:', &
        real(shift_fact / deltat), real(shift_fact / seis_dt), real(shift_fact / deltat_coarse)
     write(*,*)'   # simul. half widths per shift factor:', real(shift_fact / t_0)
     if (discrete_dirac) &
        write(*,*)'   # mesh halfwidths per shift fact', real(shift_fact / period)
     write(*,*)'  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-'
     print *
  endif

  ! strain tensor output, convert from num of dumps per period into
  ! incremental time steps
  if (dump_wavefields) then
     nstrain = floor(real(niter)/real(strain_it)) + 1

     ! This causes problems for massively parallel jobs on parallel file systems
     ! GPFS allows only a few hundreds of files in one directory
     if (diagfiles) then
        open(unit=2900+mynum,file=datapath(1:lfdata)//'/strain_info.dat'//appmynum)
        write(2900+mynum,*) nstrain
        do ielem = 1, nstrain
           write(2900+mynum,*)real(ielem)*t_0/real(strain_samp),ielem*strain_it
        enddo
        close(2900+mynum)
     endif

     if (lpr) then
        write(*,*)
        write(*,11)'    Number of wavefield dumps  :', nstrain

        write(*,12)'    ...approximately every     :', &
                   t_0/real(strain_samp), &
                   'seconds'
        write(*,13)'    ...that is, every          :',strain_it,'timestep'
     endif

     ndumppts_el = (iend - ibeg + 1) * (jend - jbeg + 1)
     if (lpr) then
        write(*,*)'    Define limitation of GLL points in the dumped fields:'
        write(*,*)'      ibeg=', ibeg, 'iend=', iend
        write(*,*)'      jbeg=', jbeg, 'jend=', jend
        write(*,*)'      # points saved within an element:', ndumppts_el
     endif
  endif

  s_max = zero

  ! Set some parameters for faster access in time loop
  half_dt = half * deltat
  half_dt_sq = half * deltat**2

  ! mesh info: coordinates of elements and collocation points
  if (diagfiles) then
      open(2222+mynum,file=infopath(1:lfinfo)//'/axial_points.dat'//appmynum)
      open(3333+mynum,file=infopath(1:lfinfo)//'/axial_ds_dz.dat'//appmynum)

      do ielem = 1, nelem

         ! write out axial points
         if (axis(ielem)) then
            call compute_coordinates(s, z, r, theta, ielem, 0, npol)
            do jpol=0, npol
               call compute_coordinates(s, z, r, theta, ielem, 0, jpol)
               write(2222+mynum,122) ielem, jpol, s, z, r, theta / pi * 180.
            enddo

            ! write out profile of grid spacing along Northern axis
            if (north(ielem)) then

               do jpol=0, npol
                  do ipol=0, npol-1
                     dsaxis(ipol,jpol) = dsqrt((scoord(ipol,jpol,ielem)-&
                          scoord(ipol+1,jpol,ielem))**2+&
                          (zcoord(ipol,jpol,ielem)-&
                          zcoord(ipol+1,jpol,ielem))**2 )
                  enddo
               enddo

               do jpol=0,npol-1
                  dzaxis(jpol) = dsqrt( (scoord(0,jpol,ielem) - &
                       scoord(0,jpol+1,ielem) )**2 + &
                       (zcoord(0,jpol,ielem) - &
                       zcoord(0,jpol+1,ielem))**2 )
               enddo
               minds(naxel) = minval(dsaxis)
               maxds(naxel) = maxval(dsaxis)
               mindz(naxel) = minval(dzaxis)
               maxdz(naxel) = maxval(dzaxis)
               call compute_coordinates(s,z,r,theta,ielem,0,npol/2)
               write(3333+mynum,123)ielem,r,minds(naxel),maxds(naxel), &
                    mindz(naxel),maxdz(naxel),maxds(naxel)/minds(naxel)
            endif
         endif
      enddo
122   format(i8,i3,4(1pe15.5))
123   format(i7,1pe14.3,5(1pe14.3))
      close(2222+mynum)
      close(3333+mynum)
  endif

  if (lpr) write(*,*)

end subroutine compute_numerical_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Writes out relevant simulation parameters, to simulation.info and the NetCDF
!! Output file.
subroutine write_parameters

    use data_comm
    use nc_helpers
    use data_mesh

    integer          :: iel,curvel,linel,seminoel,semisoel
    integer          :: curvel_solid,linel_solid,seminoel_solid,semisoel_solid
    integer          :: curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid
    integer          :: ipol,jpol,hmaxloc1(3),hminloc1(3)
    integer          :: maxprocssend_solid,maxprocsrecv_solid
    integer          :: maxprocssend_fluid,maxprocsrecv_fluid
    real(kind=dp)    :: dis1(0:npol-1,0:npol-1,nelem),dis2(0:npol-1,0:npol-1,nelem)
    real(kind=dp)    :: s,z,r,theta,rminglob,thetaminglob,rmaxglob,thetamaxglob
    real(kind=dp)    :: mysmin,myzmin,mysmax,myzmax
    real(kind=dp)    :: myrmin,mythetamin,myrmax,mythetamax
    real(kind=dp)    :: hmax,hmaxglob,hmin,hminglob
    character(len=7) :: clogic

    character(len=8)  :: mydate
    character(len=10) :: mytime
    character(len=5)  :: myzone
    character(len=24) :: mydatetime

    call date_and_time(mydate, mytime, myzone)
    write(mydatetime,1212) mydate(1:4), mydate(5:6), mydate(7:8), mytime(1:2), mytime(3:4), mytime(5:6), myzone

1212 format(A4,'-',A2,'-',A2,'T', A2,':',A2,':',A2, A5)

    if (verbose > 1) then
       write(69,*)'  writing out all relevant simulation parameters...'
       write(69,*)'  number of respective element types...'
    endif

    curvel=0; linel=0; seminoel=0; semisoel=0
    do iel=1,nelem
       if (eltype(iel) == 'curved') curvel=curvel+1
       if (eltype(iel) == 'linear') linel=linel+1
       if (eltype(iel) == 'semino') seminoel=seminoel+1
       if (eltype(iel) == 'semiso') semisoel=semisoel+1
    enddo

    curvel_solid=0; linel_solid=0; seminoel_solid=0; semisoel_solid=0
    do iel=1,nel_solid
       if (eltype(ielsolid(iel)) == 'curved') curvel_solid=curvel_solid+1
       if (eltype(ielsolid(iel)) == 'linear') linel_solid=linel_solid+1
       if (eltype(ielsolid(iel)) == 'semino') seminoel_solid=seminoel_solid+1
       if (eltype(ielsolid(iel)) == 'semiso') semisoel_solid=semisoel_solid+1
    enddo

    curvel_fluid=0; linel_fluid=0; seminoel_fluid=0; semisoel_fluid=0
    do iel=1,nel_fluid
       if (eltype(ielfluid(iel)) == 'curved') curvel_fluid=curvel_fluid+1
       if (eltype(ielfluid(iel)) == 'linear') linel_fluid=linel_fluid+1
       if (eltype(ielfluid(iel)) == 'semino') seminoel_fluid=seminoel_fluid+1
       if (eltype(ielfluid(iel)) == 'semiso') semisoel_fluid=semisoel_fluid+1
    enddo

    if (verbose > 1) write(69,*)'  grid spacing min/max...'
    do iel=1,nelem
       do ipol=0,npol-1
          do jpol=0,npol-1
             dis1(ipol,jpol,iel) = dsqrt(&
                  (scoord(ipol,jpol,iel)-scoord(ipol+1,jpol,iel))**2&
                  +(zcoord(ipol,jpol,iel)-zcoord(ipol+1,jpol,iel))**2)
             dis2(ipol,jpol,iel) = dsqrt(&
                  (scoord(ipol,jpol,iel)-scoord(ipol,jpol+1,iel))**2&
                  +(zcoord(ipol,jpol,iel)-zcoord(ipol,jpol+1,iel))**2)
          enddo
       enddo
    enddo

    if (verbose > 1) write(69,*)'  calculating hmax...'
    hmax = max(maxval(dis1),maxval(dis2))
    if (verbose > 1) write(69,*)'  hmaxstuff:',minval(dis1),minval(dis2),hmax
    hmaxglob = pmax(hmax)
    if (verbose > 1) write(69,*)'  hmaxglob:',hmaxglob
    hmaxloc1 = maxloc(dis1)
    if (maxval(dis2) < maxval(dis1)) hmaxloc1=maxloc(dis2)
    call compute_coordinates(s,z,rmaxglob,thetamaxglob,hmaxloc1(3), &
         hmaxloc1(1)-1,hmaxloc1(2)-1)
    if (verbose > 1) write(69,*)' rmax,thetamax:',rmaxglob,thetamaxglob

    if (verbose > 1) write(69,*)'  calculating hmin...'
    hmin = min(minval(dis1),minval(dis2))
    if (verbose > 1) write(69,*)'  hminstuff:',minval(dis1),minval(dis2),hmin
    hminglob = pmin(hmin)
    if (verbose > 1) write(69,*)'  hminglob:',hminglob

    hminloc1=minloc(dis1)
    if (minval(dis2) < minval(dis1)) hminloc1=minloc(dis2)
    call compute_coordinates(s,z,rminglob,thetaminglob,hminloc1(3), &
            hminloc1(1)-1,hminloc1(2)-1)

    ! Checking potential issues with input parameter consistency
    if (lpr) write(*,*)
    if (lpr) write(*,*)'  checking input parameters for consistency...'
    call check_parameters(hmaxglob,hminglob,curvel,linel,seminoel,semisoel, &
                         curvel_solid,linel_solid,seminoel_solid,semisoel_solid, &
                         curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid)

    maxprocssend_solid = pmax_int(sizesend_solid)
    maxprocsrecv_solid = pmax_int(sizerecv_solid)
    maxprocssend_fluid = pmax_int(sizesend_fluid)
    maxprocsrecv_fluid = pmax_int(sizerecv_fluid)


    ! output to stdout, only by proc nproc-1
    if (lpr) then

        write(*,*)
        write(*,*)':::::::::::::::: SIMULATION PARAMETERS::::::::::::::::::::::::'

        write(*,*)'  Code information_____________________________________'
        write(*,12)'     svn revision      :', git_hash
        write(*,12)'     username          :', username
        write(*,12)'     hostname          :', hostname
        write(*,12)'     compiler          :', compiler
        write(*,12)'     compilerversion   :', compilerversion
        write(*,20)'     FFLAGS            :', fflags
        write(*,20)'     CFLAGS            :', cflags
        write(*,20)'     LDFLAGS           :', ldflags
        write(*,12)'     OpenMP            :', openmp
        write(*,*)'  Global mesh information______________________________'
        write(*,12)'     Background model  :',bkgrdmodel
        write(*,10)'     # discontinuities :',ndisc
        write(*,13)'     Have fluid region ?',have_fluid
        write(*,11)'     Outer rad.     [m]:',router
        write(*,11)'     Inner rad.     [m]:',rmin
        write(*,10)'     Polynomial order  :',npol
        write(*,10)'     # control nodes   :',npoin
        write(*,10)'     Total elements    :',nelem
        write(*,10)'     Total # points    :',npoint
        write(*,10)'     # global numbers  :',nglob
        write(*,10)'     # axial elements  :',naxel
        write(*,10)'     # curved elements :',curvel
        write(*,10)'     # linear elements :',linel
        write(*,10)'     # mixed elements  :',seminoel+semisoel
        write(*,11)'     Min. distance  [m]:',min_distance_dim
        write(*,11)'     Min. distance/r0  :',min_distance_nondim

        write(*,*)'  Grid spacing, velocities etc.________________________'
        write(*,17)'     Min. (pre,comp)[m]:',hmin_glob,hminglob
        write(*,17)'     Max. (pre,comp)[m]:',hmax_glob,hmaxglob
        write(*,17)'     Min. vp[m/s], r[m]:',vpmin,vpminr
        write(*,17)'     Min. vs[m/s], r[m]:',vsmin,vsminr
        write(*,17)'     Max. vp[m/s], r[m]:',vpmax,vpmaxr
        write(*,17)'     Max. vs[m/s], r[m]:',vsmax,vsmaxr
        write(*,11)'     Max. lead time [s]:',char_time_max
        write(*,17)'     r [m], theta [deg]:',char_time_max_rad*router, &
                                              char_time_max_theta
        write(*,11)'     Min. lead time [s]:',char_time_min
        write(*,17)'     r [m], theta [deg]:',char_time_min_rad*router, &
                                              char_time_min_theta

        write(*,*)'  Solid-Fluid configuration____________________________'
        write(*,15)'     S/F elements      :',nel_solid,nel_fluid
        write(*,15)'     S/F # points      :',npoint_solid,npoint_fluid
        write(*,15)'     S/F global numbers:',nglob_solid,nglob_fluid
        write(*,15)'     S/F # axial elems :',naxel_solid,naxel_fluid
        write(*,10)'     # S/F boundary els:',nel_bdry
        write(*,15)'     S/F curved elems  :',curvel_solid,curvel_fluid
        write(*,15)'     S/F linear elems  :',linel_solid,linel_fluid
        write(*,15)'     S/F mixed elements:',seminoel_solid+semisoel_solid, &
                                              seminoel_fluid+semisoel_fluid

        write(*,*)'  Solid message passing_________________________________'
        write(*,10)'     # processors      :',nproc
        write(*,10)'     max. sent messages:',maxprocssend_solid
        write(*,10)'     max. sent size    :',sizemsgsendmax_solid
        write(*,10)'     nax. recv messages:',maxprocsrecv_solid
        write(*,10)'     max. recv size    :',sizemsgrecvmax_solid

        if (have_fluid) then
            write(*,*)'  Fluid message passing_________________________________'
            write(*,10)'     max. sent messages:',maxprocssend_fluid
            write(*,10)'     max. sent size    :',sizemsgsendmax_fluid
            write(*,10)'     nax. recv messages:',maxprocsrecv_fluid
            write(*,10)'     max. recv size    :',sizemsgrecvmax_fluid
        endif

        write(*,*)'  Source information___________________________________'
        write(*,16)'     Source type       :',src_type(1),src_type(2)
        write(*,11)'     Source depth   [m]:',zsrc
        write(*,11)'     Source colat [deg]:',srccolat*180./pi
        write(*,11)'     Source long  [deg]:',srclon*180./pi
        write(*,11)'     Magnitude    [N/m]:',magnitude
        write(*,12)'     Source time fct   :',trim(stf_type)
        write(*,11)'     Dom. period    [s]:',t_0
        write(*,*)'  Receiver information___________________________________'
        write(*,12)'     Receiver file type',rec_file_type
        write(*,19)'     Sum seismograms  :',sum_seis
        write(*,*)'  General numerical parameters_________________________'
        write(*,11)'     # elems/wavelength:',pts_wavelngth
        write(*,11)'     Courant number    :',courant
        write(*,11)'     Time step [s]     :',deltat
        write(*,10)'     # iterations      :',niter
        write(*,11)'     seismo length [s] :',niter*deltat
        write(*,12)'     time extrapolation:',time_scheme
        write(*,*)'  Input/Output information_____________________________'
        write(*,12)'     Output data path  :',trim(datapath)
        write(*,12)'     Output info path  :',trim(infopath)
        write(*,19)'     Sum wavefields:', sum_fields
        write(*,19)'     Dump energy       :',dump_energy
        write(*,18)'     XDMF VTK          :', dump_xdmf
        if (dump_vtk .or. dump_xdmf) then
            write(*,11)'     snap interval [s] :',snap_dt
            write(*,10)'     # snaps           :',snap_it
        endif
        write(*,19)'     Dump wavefields   :',dump_wavefields
        if (dump_wavefields) then
            write(*,12)'     Dumping type      :',dump_type
            write(*,11)'     dump interval [s] :',deltat_coarse !period/real(strain_samp)
            write(*,10)'     # wavefield dumps :',strain_it
        endif
        write(*,19)'     Need fluid displ. :',need_fluid_displ
        write(*,*)
        write(*,*)':::::::::::::::: END SIMULATION PARAMETERS::::::::::::::::::::'
        write(*,*)
        call flush(6)

        ! additionally write a header for the kernel software
        ! not used anymore
        !if (dump_wavefields) call create_kernel_header


        ! write generic simulation info file
        open(unit=55,file='simulation.info')
        write(55,23)trim(bkgrdmodel),'background model'
        write(55,21)deltat,'time step [s]'
        write(55,22)niter,'number of time steps'
        write(55,23)trim(src_type(1)),'source type'
        write(55,23)trim(src_type(2)),'source type'
        write(55,23)trim(stf_type),'source time function'
        write(55,23)trim(simtype),'simtype'
        write(55,21)period,'dominant source period'
        write(55,21)src_depth/1000.,'source depth [km]'

        write(55,21) srccolat, 'Source colatitude'
        write(55,21) srclon, 'Source longitude'

        write(55,25)magnitude,'scalar source magnitude'
        write(55,22)num_rec_tot,'number of receivers'
        write(55,22)nseismo,'length of seismogram [time samples]'
        write(55,21)real(deltat)*real(seis_it),'seismogram sampling [s]'
        if (dump_wavefields) then
            write(55,22) nstrain,'number of strain dumps'
            write(55,21) deltat_coarse, 'strain dump sampling rate [s]'
        else
            write(55,22)0,'number of strain dumps'
            write(55,21)0.,'strain dump sampling rate [s]'
        endif
        if (dump_vtk .or. dump_xdmf) then
            write(55,22) nsnap,'number of snapshot dumps'
            write(55,21)deltat*real(snap_it),'snapshot dump sampling rate [s]'
        else
            write(55,22)0,'number of snapshot dumps'
            write(55,21)0.,'snapshot dump sampling rate [s]'
        endif
        ! just not to cause trouble when reading simulation.info linewise:
        write(55,23)'cyl','receiver components '
        write(55,22)ibeg,'  ibeg: beginning GLL index for wavefield dumps'
        write(55,22)iend,'iend: end GLL index for wavefield dumps'
        write(55,21)shift_fact,'source shift factor [s]'
        write(55,22)nint(shift_fact/deltat),'source shift factor for deltat'
        write(55,22)nint(shift_fact/seis_dt),'source shift factor for seis_dt'
        write(55,22)nint(shift_fact/deltat_coarse),'source shift factor for deltat_coarse'
        write(55,23)trim(rec_file_type),'receiver file type'
        write(55,21)dtheta_rec,'receiver spacing (0 if not even)'
        write(55,24)use_netcdf,'use netcdf for wavefield output?'

        write(55,22) nelem, 'nelem'
        write(55,22) nel_fluid, 'nel_fluid'
        write(55,22) nproc, 'nproc'

        close(55) ! simulation.info

        write(*,*)
        write(*,*)'  wrote general simulation info into "simulation.info"'

21      format(f22.7,a45)
22      format(i20,a45)
23      format(a20,a45)
24      format(l20,a45)
25      format(1pe15.5,a45)

    endif ! lpr

! in case of parallel IO, all ranks write attributes
#ifdef enable_parallel_netcdf
    if ((use_netcdf)) then
#else
    if ((mynum == 0) .and. (use_netcdf)) then !Only proc0 has the netcdf file open at that point
#endif
        ! write generic simulation info file
        if (mynum == 0) write(*,*) ' Writing simulation info to netcdf file attributes'
        call nc_write_att_int(  7,                     'file version')
        call nc_write_att_char( trim(bkgrdmodel),      'background model')
        call nc_write_att_int(  merge(1, 0, do_anel),  'attenuation') ! merge: hacky conversion of logical to int
        call nc_write_att_dble( router / 1000,         'planet radius')
        call nc_write_att_char( trim(mydatetime),      'datetime')
        call nc_write_att_char( trim(git_hash),        'git commit hash')
        call nc_write_att_char( trim(username),        'user name')
        call nc_write_att_char( trim(hostname),        'host name')
        call nc_write_att_char( trim(compiler),        'compiler brand')
        call nc_write_att_char( trim(compilerversion), 'compiler version')
        call nc_write_att_char( trim(fflags),          'FFLAGS')
        call nc_write_att_char( trim(cflags),          'CFLAGS')
        call nc_write_att_char( trim(ldflags),         'LDFLAGS')
        call nc_write_att_char( trim(openmp),          'OpenMP')
        call nc_write_att_char( trim(time_scheme),     'time scheme')
        call nc_write_att_real( real(deltat),          'time step in sec')
        call nc_write_att_int(  niter,                 'number of time steps')
        call nc_write_att_int(  npol,                  'npol')
        call nc_write_att_char( trim(src_type(1)),     'excitation type')
        call nc_write_att_char( trim(src_type(2)),     'source type')
        call nc_write_att_char( trim(stf_type),        'source time function')
        call nc_write_att_char( trim(simtype),         'simulation type')
        call nc_write_att_real( real(period),          'dominant source period')
        call nc_write_att_real( real(src_depth/1000.), 'source depth in km')

        call nc_write_att_real( real(srccolat),        'Source colatitude')
        call nc_write_att_real( real(srclon),          'Source longitude' )

        call nc_write_att_real( real(magnitude),       'scalar source magnitude')
        call nc_write_att_int(  num_rec_tot,           'number of receivers')
        call nc_write_att_int(  nseismo,               'length of seismogram  in time samples')
        call nc_write_att_real( real(deltat)*real(seis_it), 'seismogram sampling in sec')
        if (dump_wavefields) then
           call nc_write_att_int(nstrain,              'number of strain dumps')
           call nc_write_att_dble(deltat_coarse,       'strain dump sampling rate in sec')
           call nc_write_att_char(trim(dump_type),     'dump type (displ_only, displ_velo, fullfields)')
           if (dump_type == 'displ_only' .or. dump_type == 'strain_only') then
              call nc_write_att_dble(kwf_rmin / 1000,  'kernel wavefield rmin')
              call nc_write_att_dble(min(kwf_rmax, router) / 1000,  'kernel wavefield rmax')
              call nc_write_att_dble(kwf_thetamin / pi * 180.,  'kernel wavefield colatmin')
              call nc_write_att_dble(kwf_thetamax / pi * 180.,  'kernel wavefield colatmax')
           else
              call nc_write_att_dble(0d0,               'kernel wavefield rmin')
              call nc_write_att_dble(router / 1000,     'kernel wavefield rmax')
              call nc_write_att_dble(0d0,               'kernel wavefield colatmin')
              call nc_write_att_dble(0d0,               'kernel wavefield colatmax')
           endif
        else
           call nc_write_att_int(0,                    'number of strain dumps')
           call nc_write_att_dble(0.d0,                'strain dump sampling rate in sec' )
        endif
        if (dump_vtk) then
           call nc_write_att_int(nsnap,                'number of snapshot dumps')
           call nc_write_att_real(real(deltat)*real(snap_it), 'snapshot dump sampling rate in sec')
        else
           call nc_write_att_int(0,                    'number of snapshot dumps')
           call nc_write_att_real(0.,                  'snapshot dump sampling rate in sec')
        endif
        call nc_write_att_char( 'cyl',                 'receiver components ')
        call nc_write_att_int(  ibeg,                  'ibeg')
        call nc_write_att_int(  iend,                  'iend')
        call nc_write_att_int(  jbeg,                  'jbeg')
        call nc_write_att_int(  jend,                  'jend')
        call nc_write_att_real( shift_fact,            'source shift factor in sec')
        call nc_write_att_int(  nint(shift_fact/deltat),  'source shift factor for deltat')
        call nc_write_att_int(  nint(shift_fact/seis_dt), 'source shift factor for seis_dt')
        call nc_write_att_int(  nint(shift_fact/deltat_coarse), 'source shift factor for deltat_coarse')
        call nc_write_att_char( trim(rec_file_type),   'receiver file type')
        call nc_write_att_real( dtheta_rec,            'receiver spacing (0 if not even)')
        write(clogic,*) use_netcdf
        call nc_write_att_char( clogic,                'use netcdf for wavefield output?')

        call nc_write_att_int( 0,                      'percent completed')
        call nc_write_att_int( 0,                      'finalized')
    endif


    ! output for each processor==============================================

    ! extract processor location
    mysmin = router
    myzmin = router
    mysmax = zero
    myzmax = zero
    myrmin = router
    mythetamin = 10.*pi
    myrmax = zero
    mythetamax = zero

    do iel=1,nelem
        do ipol=0,npol
            do jpol=0,npol
                call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
                if (s < mysmin) mysmin = s
                if (s > mysmax) mysmax = s
                if (r < myrmin) myrmin = r
                if (z < myzmin) myzmin = z
                if (z > myzmax) myzmax = z
                if (r > myrmax) myrmax = r
                if (theta < mythetamin) mythetamin = theta
                if (theta > mythetamax) mythetamax = theta
            enddo
        enddo
    enddo

    if (verbose > 1) then
       write(69,*)
       write(69,15)'My rank, total procs    :', mynum,nproc
       write(69,17)'Min./max. s [m]         :', mysmin,mysmax
       write(69,17)'Min./max. z [m]         :', myzmin,myzmax
       write(69,17)'Min./max. r [m]         :', myrmin,myrmax
       write(69,17)'Min./max. theta [deg]   :', mythetamin*180./pi, mythetamax*180./pi

       write(69,10)'Axial total elems       :', naxel
       write(69,10)'Axial solid elems       :', naxel_solid
       write(69,10)'Axial fluid elems       :', naxel_fluid

       write(69,13)'Have source             ?', have_src
       if (have_src) then
           write(69,11)'Depth asked for      [m]:', src_depth
           write(69,11)'Computed depth       [m]:', router - &
                                         zcoord(ipol_src,jpol_src,ielsolid(iel_src))
       endif
       write(69,13)'Have boundary els       ?', have_bdry_elem
       if (have_bdry_elem) then
           write(69,10)'# boundary elements     :', nel_bdry
       endif

       write(69,*)
       write(69,*)'Solid message passing_____________________________'
       write(69,10)' # recv messages        :', sizerecv_solid
       write(69,10)' Max. size recv messages:', sizemsgrecvmax_solid
       write(69,10)' # sent messages        :', sizesend_solid
       write(69,10)' Max. size sent messages:', sizemsgsendmax_solid

       if (have_fluid) then
           write(69,*)'Fluid message passing_____________________________'
           write(69,10)' # recv messages        :', sizerecv_fluid
           write(69,10)' Max. size recv messages:', sizemsgrecvmax_fluid
           write(69,10)' # sent messages        :', sizesend_fluid
           write(69,10)' Max. size sent messages:', sizemsgsendmax_fluid
       endif !have_fluid

       call flush(69)
    endif

10  format(a25,i14)
11  format(a25,1pe14.5)
12  format(a25,'   ',a18)
13  format(a25,L14)
15  format(a25,i14,i9)
16  format(a25,' ',a12,a10)
17  format(a25,2(1pe13.3))
18  format(a25,2(L14))
19  format(a25,L14)
20  format(a25,'   ', a51)

    ! write post processing file==============================================

    if (lpr) then

        if (trim(simtype) == 'moment' .and. src_type(2) == 'mrr') then
           open(unit=9, file="../param_post_processing")
        else if (trim(simtype) == 'single') then
           open(unit=9, file="param_post_processing")
        endif

        if ((trim(simtype) == 'moment' .and. src_type(2) == 'mrr') &
            .or. (trim(simtype) == 'single')) then

           write(*,*)'  Writing post processing input file: param_post_processing'
           write(*,*)'  ... mainly based on guessing from the current simulation, make sure to edit!'
           write(9,'(a,/,a,/,a,/)') &
                    '# receiver coordinate system', &
                    '# one of: enz, sph, cyl, xyz, src', &
                    'REC_COMP_SYS    enz'

           write(9,'(a,/,a,/,a)') &
                    '# period of source time funtion to be convolved', &
                    '# should be larger then the mesh period', &
                    '# 0. to not convolve'
           if (stf_type == 'dirac_0' .or. stf_type == 'quheavi' ) then
              write(9,'(a,f8.4,/)') &
                    'CONV_PERIOD     ', period
           else
              write(9,'(a,/)') &
                    'CONV_PERIOD     0.'
           endif

           write(9,'(6(a,/))') &
                    '# source time function (moment function) to convolve with', &
                    '# ''gauss_0'': Gaussian with dominant period SOURCE_PERIOD', &
                    '# ''gauss_1'': 1st derivative of gauss_0', &
                    'CONV_STF        gauss_0'

           write(9,'(a,/,a,l1/)') &
                    '# make 3D plots of the wavefield', &
                    'LOAD_SNAPS      ', dump_vtk

           write(9,'(a,/,a,/)') &
                    '# OUTPUT PATH', &
                    'DATA_DIR        "./Data_Postprocessing"'

           write(9,'(a,/,a,/)') &
                    '# Write out intermediate seismograms (processed, but not summed)', &
                    'DETAILED_OUTPUT false'

           write(9,'(a,/,a,/,a,/,/)') &
                    '# output seismograms at negative time', &
                    '# (to correct for finite width of the source time function', &
                    'NEGATIVE_TIME   T'

           write(9,'(a,/,a,/,a,/,a,/,a,/)') &
                    '#############################################', &
                    '# options the 3D wavefield plots', &
                    '# crossection location, starting and ending phi', &
                    '3D_PHI_START     0.', &
                    '3D_PHI_END      85.'

           write(9,'(a,/,a,f6.0/,a,/)') &
                    '# radius of top and bottom layer in km', &
                    '3D_RTOP        ', router/1000., &
                    '3D_RBOT         3190.'

           write(9,'(a,/,a,/,a,/)') &
                    '# switches for bottom, top and meridonial surface', &
                    '3D_PLOT_TOP     T', &
                    '3D_PLOT_BOT     T'

           write(9,'(a,/,a,/,a,i5,/,a)') &
                    '# time snapshot selection:', &
                    '3D_SNAP_BEG      1', &
                    '3D_SNAP_END  ', nsnap, &
                    '3D_SNAP_STRIDE   1'
           close(9)
        endif

        write(*,*)'    ... wrote file param_post_processing'
    endif

end subroutine write_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! < Checking some mesh parameters and message parsing
subroutine check_parameters(hmaxglob, hminglob, curvel, linel, seminoel, semisoel, &
                       curvel_solid, linel_solid, seminoel_solid, semisoel_solid, &
                       curvel_fluid, linel_fluid, seminoel_fluid, semisoel_fluid)

    use data_comm
    use data_mesh


    real(kind=dp)   , intent(in) :: hmaxglob,hminglob
    integer, intent(in) :: curvel,linel,seminoel,semisoel
    integer, intent(in) :: curvel_solid,linel_solid,seminoel_solid,semisoel_solid
    integer, intent(in) :: curvel_fluid,linel_fluid,seminoel_fluid,semisoel_fluid

    if (verbose > 1) write(*,*) procstrg, 'Checking solid message-passing...'
    if (nproc == 1 .and. psum_int(sizesend_solid) > 0 ) then
       write(*,*)'Problem: Have only one proc but want to send messages..'
       stop
    endif

    if (nproc == 1 .and. psum_int(sizerecv_solid) > 0 ) then
       write(*,*)'Problem: Have only one proc but want to receive messages...'
       stop
    endif

    if (nproc > 1 .and. psum_int(sizesend_solid) == 0 ) then
       write(*,*)'Problem: No proc is willing to send anything....'
       stop
    endif

    if (nproc > 1 .and. psum_int(sizesend_solid) == 0 ) then
       write(*,*)'Problem: No proc is willing to receive anything....'
       stop
    endif

    if (psum_int(sizesend_solid) < nproc-1 ) then
       write(*,*)'Problem: Some proc(s) not willing to send anything...'
       stop
    endif

    if (psum_int(sizerecv_solid) < nproc-1 ) then
       write(*,*)'Problem: Some proc(s) not willing to receive anything...'
       stop
    endif

    if (have_fluid) then
       if (verbose > 1) write(*,*)procstrg,'Checking fluid message-passing...'
       if (nproc == 1 .and. psum_int(sizesend_fluid) > 0 ) then
          write(*,*)'Problem: Have only one proc but want to send messages..'
          stop
       endif

       if (nproc == 1 .and. psum_int(sizerecv_fluid) > 0 ) then
          write(*,*)'Problem: Have only one proc but want to receive messages...'
          stop
       endif

       if (nproc > 1 .and. psum_int(sizesend_fluid) == 0 ) then
          write(*,*)'Problem: No proc is willing to send anything....'
          stop
       endif

       if (nproc > 1 .and. psum_int(sizesend_fluid) == 0 ) then
          write(*,*)'Problem: No proc is willing to receive anything....'
          stop
       endif

       if (psum_int(sizesend_fluid) < nproc-1 ) then
          write(*,*)'Problem: Some proc(s) not willing to send anything...'
          stop
       endif

       if (psum_int(sizerecv_fluid) < nproc-1 ) then
          write(*,*)'Problem: Some proc(s) not willing to receive anything...'
          stop
       endif

    endif !have_fluid

  ! Even more tests.............
  ! stop if difference between loaded and on-the-fly mesh larger than 1 METER....
  if ( (hmin_glob-hminglob) > 1.) then
      write(*,*)
      write(*,*)mynum,'Problem with minimal global grid spacing!'
      write(*,*)mynum,'Value from mesher        :',hmin_glob
      write(*,*)mynum,'Value computed on-the-fly:',hminglob
      stop
  endif

  ! stop if difference between loaded and on-the-fly mesh larger than 1 METER....
  if ((hmax_glob-hmaxglob) > 1.) then
      write(*,*)
      write(*,*)mynum,'Problem with maximal global grid spacing!'
      write(*,*)mynum,'Value from mesher        :',hmax_glob
      write(*,*)mynum,'Value computed on-the-fly:',hmaxglob
      stop
  endif

  ! stop if sum of respective element types do not sum up to nelem
  if (curvel+linel+seminoel+semisoel /= nelem) then
      write(*,*)
      write(*,*)mynum,'Problem with number of assigned global element types!'
      write(*,*)mynum,'curved,lin,semi(N),semi(S):',curvel,linel,seminoel, &
                                                    semisoel
      write(*,*)mynum,'Total # elements          :',nelem
      stop
  endif

  if (curvel_solid+linel_solid+seminoel_solid+semisoel_solid /= nel_solid) then
      write(*,*)
      write(*,*)mynum,'Problem with number of assigned solid element types!'
      write(*,*)mynum,'curved,lin,semi(N),semi(S):',curvel_solid,linel_solid, &
                                                   seminoel_solid,semisoel_solid
      write(*,*)mynum,'Total # elements          :',nel_solid
      stop
  endif

  if (curvel_fluid+linel_fluid+seminoel_fluid+semisoel_fluid /= nel_fluid) then
      write(*,*)
      write(*,*)mynum,'Problem with number of assigned fluid element types!'
      write(*,*)mynum,'curved,lin,semi(N),semi(S):',curvel_fluid,linel_fluid, &
                                                   seminoel_fluid,semisoel_fluid
      write(*,*)mynum,'Total # elements          :',nel_fluid
      stop
  endif

end subroutine check_parameters
!-----------------------------------------------------------------------------------------

end module parameters
!=========================================================================================
