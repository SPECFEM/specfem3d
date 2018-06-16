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
!
!=========================================================================================
!
!> Contains all functions for the wave propagation. prepare_waves has to be
!! called beforehand and then time_loop is the only allowed entry point to start
!! wave propagation.
module time_evol_wave

  use global_parameters
  use data_proc
  use data_mesh
  use data_source
  use data_time
  use seismograms
  use rotations
  use data_io, only: verbose
  use coupling_mod, only: dump_field_1d_cp !, dump_wavefields_mesh_1d_cp VM

  implicit none
  public :: prepare_waves, time_loop
  private

contains

!-----------------------------------------------------------------------------------------
!> Contains all the preliminaries to propagate waves; such as the
!! background model, the stiffness and mass terms, the source and receiver
!! parameters, and preparations for I/O (dumping meshes, opening files).
subroutine prepare_waves

  use data_io
  use parameters
  use def_precomp_terms
  use source
  use clocks_mod
  use meshes_io
  use attenuation, only: dump_memory_vars
  use commun

  character(len=120) :: fname

  if (lpr) then
     print *
     write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     if (have_fluid) then
        write(*,*)'++++++++    SEISMIC WAVE PROPAGATION: SOLID-FLUID CASE  ++++++++'
     else
        write(*,*)'+++++++++++  SEISMIC WAVE PROPAGATION: SOLID CASE  +++++++++++++'
     endif
     write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif

  ! read source parameters from sourceparams.dat
  call read_sourceparams

  ! rotation if source is not on the axis
  if (rot_src ) call def_rot_matrix

  ! build mapping to avoid duplicate points at element boundaries
  if (use_netcdf .and. (trim(dump_type) == 'displ_only' &
                        .or. trim(dump_type) == 'strain_only')) &
     call build_kwf_grid()

  ! compute/output some more parameters
  call compute_numerical_parameters

  ! Define velocity/density model (velocities in m/s, density in kg/m^3 ) AND
  ! compute all global matrices (Jacobian, mapping, mass matrix, S/F boundary)
  ! for solid and fluid domains respectively
  call read_model_compute_terms

  ! compute source time function
  call compute_stf

  ! compute source location within mesh and moment tensor/single force field
  call compute_src

  ! Prepare output
  if ( dump_energy .and. lpr) then ! only one proc dumps
     write(*,*)'  opening files for kinetic/potential energy...'
     if (have_fluid) then
        open(unit=4444,file=datapath(1:lfdata)//'/energy_sol.dat')
        open(unit=4445,file=datapath(1:lfdata)//'/energy_flu.dat')
     endif
     open(unit=4446,file=datapath(1:lfdata)//'/energy_glob.dat')
  endif

  if (dump_wavefields) then
     if (lpr) write(*,*)'  dumping strain mesh and associated fields...'
     !!!!! SB VM has changed with a case fullfields and coupling
     !!!! OLD VERSION :
     !if (dump_wavefields) then
     !   if (lpr) write(*,*)'  dumping strain mesh and associated fields...'
     !   select case (dump_type)
     !   case ('fullfields')
     !      call dump_wavefields_mesh_1d
     !      !case ('coupling') ! VM VM dump local mesh only
     !      !   call dump_wavefields_mesh_1d_cp
     !   end select
     !endif
     !!!!! END OLD VERSION

     call dump_wavefields_mesh_1d
  endif

  if (dump_vtk) then
     if (lpr) write(*,*)'  dumping global grids for snapshots...'
     call dump_glob_grid_midpoint(ibeg,iend,jbeg,jend)
  endif

  if (dump_xdmf) then
     if (lpr) write(*,*)'  dumping mesh for xdmf snapshots...'

     if (.not. use_netcdf) then
         fname = datapath(1:lfdata)//'/xdmf_snap_s_' //appmynum//'.dat'

#if defined(_CRAYFTN)
         open(13100, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', position='append')
#else
         open(13100, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', convert='big_endian', position='append')
#endif

         if (.not. src_type(1) == 'monopole') then
             fname = datapath(1:lfdata)//'/xdmf_snap_p_' //appmynum//'.dat'
#if defined(_CRAYFTN)
             open(13101, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', position='append')
#else
             open(13101, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', convert='big_endian', position='append')
#endif
         endif

         fname = datapath(1:lfdata)//'/xdmf_snap_z_' //appmynum//'.dat'
#if defined(_CRAYFTN)
         open(13102, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', position='append')
#else
         open(13102, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', convert='big_endian', position='append')
#endif

         fname = datapath(1:lfdata)//'/xdmf_snap_trace_' //appmynum//'.dat'
#if defined(_CRAYFTN)
         open(13103, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', position='append')
#else
         open(13103, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', convert='big_endian', position='append')
#endif

         fname = datapath(1:lfdata)//'/xdmf_snap_curlip_' //appmynum//'.dat'
#if defined(_CRAYFTN)
         open(13104, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', position='append')
#else
         open(13104, file=trim(fname), access='stream', status='unknown', &
             form='unformatted', convert='big_endian', position='append')
#endif
     endif
     call dump_xdmf_grid()
  endif

  if (anel_true .and. dump_memory_vars) &
     call prepare_mesh_memoryvar_vtk()

  ! Various seismogram output preparations...
  call prepare_seismograms
  call open_hyp_epi_equ_anti

  ! allow for different types of receiver files
  call prepare_from_recfile_seis

  if (lpr) then ! This has to be called by just one processor. Since 0 will have to
                ! do more stuff further below, let's assign lpr to this task
     ! dump meshes for displ_only kwf output
     if (dump_wavefields .and. dump_type == "displ_only") then
        call dump_kwf_midpoint_xdmf(datapath(1:lfdata)//'/axisem_output.nc4', &
                                    npoint_kwf_global, nelem_kwf_global)
        call dump_kwf_fem_xdmf(datapath(1:lfdata)//'/axisem_output.nc4', &
                                    npoint_kwf_global, nelem_kwf_global)
        call dump_kwf_sem_xdmf(datapath(1:lfdata)//'/axisem_output.nc4', &
                                    npoint_kwf_global, nelem_kwf_global)

     else if (dump_wavefields .and. dump_type == "strain_only") then
        call dump_kwf_gll_xdmf(datapath(1:lfdata)//'/axisem_output.nc4', &
                                    npoint_kwf_global)
     endif
  endif !lpr

  ! Need to reload old seismograms and add results
  if (isim > 1 .and. sum_seis ) then
     if (lpr) write(*,*)' Running multiple simulations and summing seismograms'
     if (lpr) write(*,*)' ...implementation of multiple simulations not finished'
     stop
  endif

  ! Need to reload old seismograms and add results
  if (isim > 1 .and. sum_fields) then
     if (lpr) write(*,*)' Running multiple simulations and summing wavefields'
     if (lpr) write(*,*)' ...implementation of multiple simulations not finished'
     stop
  endif

  ! write out seismic & numerical information on the simulation
  ! and run some tests on consistency of mesh/spacing/element types/messaging
  call write_parameters

  call barrier
  if (lpr) write(*,*) 'done preparing waves.'

end subroutine prepare_waves
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Entry point into the module time_evol_wave. Calls specific time loop
! !functions, either for newmark or symplectic time scheme.
subroutine time_loop

  use clocks_mod, only: tick

  iclockold = tick()

  if (time_scheme == 'newmark2') then
     call sf_time_loop_newmark
  else
     call symplectic_time_loop
  endif

  iclockold = tick(id=idold, since=iclockold)

end subroutine time_loop
!-----------------------------------------------------------------------------------------

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           T I M E   E X T R A P O L A T I O N   R O U T I N E S
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!-----------------------------------------------------------------------------------------
!> The conventional explicit, acceleration-driven Newmark scheme of 2nd order.
!! (e.g. Chaljub & Valette, 2004). The mass matrix is diagonal; we only store
!! its pre-assembled inverse at the stage of the time loop.
!! Explicit axial masking follows Nissen-Meyer et al. 2007, GJI,
!! "Spherical-earth Frechet sensitivity kernels" eqs. (80)-(82).
!! Note that the ordering (starting inside the fluid) is crucial such that no
!! iterations for the boundary terms are necessary.
!! Also note that our definition of the fluid potential is different from
!! Chaljub & Valette and the code SPECFEM by an inverse density factor.
!! This is the correct choice for our case of non-gravitating Earth models,
!! but shall be altered once gravity is taken into account.
subroutine sf_time_loop_newmark

  use commun
  use global_parameters
  use apply_masks
  use stiffness_mono
  use stiffness_di
  use stiffness_quad
  use stiffness_fluid
  use clocks_mod
  use data_matr, only: inv_mass_rho, inv_mass_fluid
  use attenuation, only: time_step_memvars
  use attenuation, only: n_sls_attenuation
  use attenuation, only: att_coarse_grained
  use data_mesh

  ! Solid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: disp, velo
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: acc0, acc1

  ! solid memory variables
  real(kind=realkind), allocatable :: memory_var(:,:,:,:,:)
  real(kind=realkind), allocatable :: memory_var_cg4(:,:,:,:)

  ! Fluid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: chi, dchi
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: ddchi0, ddchi1

  integer           :: iseismo = 0 ! < current seismogram sample
  integer           :: istrain = 0 ! < current kernel wavefield sample
  integer           :: isnap   = 0 ! < current wavefield sample for movies
  integer           :: iter

  if (lpr) then
     write(*,*)
     write(*,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
     write(*,*)'TTTT  2nd-order, acceleration-driven Newmark time scheme TTTTT'
     write(*,*)'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
     write(*,*)
  endif

  if (anel_true) then
     if (att_coarse_grained) then
        allocate(memory_var_cg4(1:4,6,n_sls_attenuation,nel_solid))
        memory_var_cg4 = 0
     else
        allocate(memory_var(0:npol,0:npol,6,n_sls_attenuation,nel_solid))
        memory_var = 0
     endif
  endif

  ! INITIAL CONDITIONS
  ! initializiation with a small number prevents performance loss (~factor 3) due to
  ! denormal floats in the onset of the p-wave (going from zero to some finite value)
  !
  ! alternatively, compilerflags -ffast-math (gfortran) or -ftz (ifort) might
  ! give the same speedup, but seem to be unstable on some systems
  !
  ! another alternative (IMPLEMENTED NOW, see ftz.c, called in main.f90):
  ! http://software.intel.com/en-us/articles/how-to-avoid-performance-penalties-for-gradual-underflow-behavior
  !
  ! some more reading:
  ! http://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x

  disp = zero
  velo = zero
  acc0 = zero
  acc1 = zero

  chi = zero
  dchi = zero
  ddchi0 = zero
  ddchi1 = zero

  t = zero

  if (lpr) write(*,*) '************ S T A R T I N G   T I M E   L O O P *************'
  if (verbose > 1) write(69,*) &
        '************ S T A R T I N G   T I M E   L O O P *************'

  iclockdump = tick()
  call dump_stuff(0, iseismo, istrain, isnap, disp, velo, chi, dchi, ddchi0, t)
  iclockdump = tick(id=iddump, since=iclockdump)

  do iter = 1, niter

     t = t + deltat
     call runtime_info(iter,disp,chi)

     chi = chi +  deltat * dchi + half_dt_sq * ddchi0

     disp(:,:,:,1) = disp(:,:,:,1) + deltat * velo(:,:,:,1) + half_dt_sq * acc0(:,:,:,1)
     if (src_type(1) /= 'monopole') &
        disp(:,:,:,2) = disp(:,:,:,2) + deltat * velo(:,:,:,2) + half_dt_sq * acc0(:,:,:,2)
     disp(:,:,:,3) = disp(:,:,:,3) + deltat * velo(:,:,:,3) + half_dt_sq * acc0(:,:,:,3)

     if (src_type(1) /= 'monopole') &
        call apply_axis_mask_scal(chi, nel_fluid, ax_el_fluid, naxel_fluid)

     iclockstiff = tick()
     call glob_fluid_stiffness(ddchi1, chi)
     iclockstiff = tick(id=idstiff, since=iclockstiff)

     call bdry_copy2fluid(ddchi1, disp)

     if (src_type(1) /= 'monopole') &
        call apply_axis_mask_scal(ddchi1, nel_fluid, ax_el_fluid, naxel_fluid)

     iclockcomm = tick()
     call pdistsum_fluid(ddchi1, phase=1)
     iclockcomm = tick(id=idcomm, since=iclockcomm)

     iclockstiff = tick()
     select case (src_type(1))
     case ('monopole')
        call apply_axis_mask_onecomp(disp, nel_solid, ax_el_solid, naxel_solid)
        call glob_stiffness_mono(acc1, disp)

        if (anel_true) then
           iclockanelst = tick()
           call glob_anel_stiffness_mono(acc1, memory_var, memory_var_cg4, &
                                         att_coarse_grained)
           iclockanelst = tick(id=idanelst, since=iclockanelst)
        endif

     case ('dipole')
        call apply_axis_mask_twocomp(disp, nel_solid, ax_el_solid, naxel_solid)
        call glob_stiffness_di(acc1, disp)

        if (anel_true) then
           iclockanelst = tick()
           call glob_anel_stiffness_di(acc1, memory_var, memory_var_cg4, &
                                         att_coarse_grained)
           iclockanelst = tick(id=idanelst, since=iclockanelst)
        endif

     case ('quadpole')
        call apply_axis_mask_threecomp(disp, nel_solid, ax_el_solid, naxel_solid)
        call glob_stiffness_quad(acc1, disp)

        if (anel_true) then
           iclockanelst = tick()
           call glob_anel_stiffness_quad(acc1, memory_var, memory_var_cg4, &
                                         att_coarse_grained)
           iclockanelst = tick(id=idanelst, since=iclockanelst)
        endif
     end select
     iclockstiff = tick(id=idstiff, since=iclockstiff)

     iclockcomm = tick()
     call pdistsum_fluid(ddchi1, phase=2)
     iclockcomm = tick(id=idcomm, since=iclockcomm)

     ddchi1 = - inv_mass_fluid * ddchi1
     call bdry_copy2solid(acc1, ddchi1)

     select case (src_type(1))
     case ('monopole')
        call apply_axis_mask_onecomp(acc1, nel_solid, ax_el_solid, naxel_solid)

     case ('dipole')
        call apply_axis_mask_twocomp(acc1, nel_solid, ax_el_solid, naxel_solid)

     case ('quadpole')
        call apply_axis_mask_threecomp(acc1, nel_solid, ax_el_solid, naxel_solid)
     end select

     iclockcomm = tick()
     call pdistsum_solid(acc1, phase=1)
     iclockcomm = tick(id=idcomm, since=iclockcomm)

     if (anel_true) then
        iclockanelts = tick()
        call time_step_memvars(memory_var, memory_var_cg4, disp, att_coarse_grained)
        iclockanelts = tick(id=idanelts, since=iclockanelts)
     endif

     dchi = dchi + half_dt * (ddchi0 + ddchi1)
     ddchi0 = ddchi1

     iclockcomm = tick()
     call pdistsum_solid(acc1, phase=2)
     iclockcomm = tick(id=idcomm, since=iclockcomm)

     call add_source(acc1, stf(iter))

     acc1(:,:,:,1) = - inv_mass_rho * acc1(:,:,:,1)

     if (src_type(1) /= 'monopole') &
        acc1(:,:,:,2) = - inv_mass_rho * acc1(:,:,:,2)

     if (src_type(1) == 'dipole') then
        ! for the factor 2 compare eq 32 in TNM (2006)
        acc1(:,:,:,3) = - two * inv_mass_rho * acc1(:,:,:,3)
     else
        acc1(:,:,:,3) = - inv_mass_rho * acc1(:,:,:,3)
     endif

     velo(:,:,:,1) = velo(:,:,:,1) + half_dt * (acc0(:,:,:,1) + acc1(:,:,:,1))
     if (src_type(1) /= 'monopole') &
        velo(:,:,:,2) = velo(:,:,:,2) + half_dt * (acc0(:,:,:,2) + acc1(:,:,:,2))
     velo(:,:,:,3) = velo(:,:,:,3) + half_dt * (acc0(:,:,:,3) + acc1(:,:,:,3))

     acc0(:,:,:,1) = acc1(:,:,:,1)
     if (src_type(1) /= 'monopole') &
        acc0(:,:,:,2) = acc1(:,:,:,2)
     acc0(:,:,:,3) = acc1(:,:,:,3)

     !!!!! CHANGES SB
     !!!! OLD
     !iclockdump = tick()
     !call dump_stuff(iter, disp, velo, chi, dchi, ddchi0, t) !! VM VM ecriture des champs sur disque
     !iclockdump = tick(id=iddump, since=iclockdump)
     !!!!! END OLD
     iclockdump = tick()
     call dump_stuff(iter, iseismo, istrain, isnap, disp, velo, chi, dchi, ddchi0, t)
     iclockdump = tick(id=iddump, since=iclockdump)
     ! SB  something tochange ?

  enddo ! time loop

end subroutine sf_time_loop_newmark
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> SOLVE coupled solid-fluid system of temporal ODEs:
!!   M*\dot{u}    = -K*u - B*\ddot{\chi} + F (solid)
!!   M*\ddot{chi} = -K*\chi - B*u (fluid)
!! using symplectic time integration schemes of 4th, 6th, 8th, 10th order
!!
!! The time step can be chosen 1.5 times larger than in Newmark, resulting in CPU
!! times about 2.5 times longer than Newmark, but considerably more accurate.
!! Consult Ampuero & Nissen-Meyer (to be submitted 2007) for examples of when this
!! choice might be more appropriate. Generally, for long propagation distances
!! (say, > 100 wavelengths), it is worthwhile considering this.
subroutine symplectic_time_loop

  use global_parameters
  use commun
  use apply_masks
  use stiffness_mono
  use stiffness_di
  use stiffness_quad
  use stiffness_fluid
  use clocks_mod
  use source, only: compute_stf_t
  use data_matr, only: inv_mass_rho,inv_mass_fluid
  use attenuation, only: time_step_memvars
  use attenuation, only: n_sls_attenuation
  use attenuation, only: att_coarse_grained
  use data_mesh

  ! solid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: disp, velo, acc

  ! solid memory variables
  real(kind=realkind), allocatable :: memory_var(:,:,:,:,:)
  real(kind=realkind), allocatable :: memory_var_cg4(:,:,:,:)

  ! fluid fields
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: chi, dchi, ddchi

  integer   :: iseismo = 0 ! < current seismogram sample
  integer   :: istrain = 0 ! < current kernel wavefield sample
  integer   :: isnap   = 0 ! < current wavefield sample for movies
  integer   :: iter, i

  ! symplectic stuff
  real(kind=dp), allocatable, dimension(:) :: coefd, coeff, coefv, subdt
  real(kind=dp), allocatable, dimension(:) :: stf_symp

  ! choose symplectic scheme (4,6,8,10th order) and compute coefficients
  call symplectic_coefficients(coefd, coeff, coefv)
  allocate(subdt(nstages), stf_symp(nstages))

  if (anel_true) then
     if (att_coarse_grained) then
        allocate(memory_var_cg4(1:4,6,n_sls_attenuation,nel_solid))
        memory_var_cg4 = 0
     else
        allocate(memory_var(0:npol,0:npol,6,n_sls_attenuation,nel_solid))
        memory_var = 0
     endif
  endif

  ! INITIAL CONDITIONS
  disp = zero
  velo = zero
  acc = zero

  chi = zero
  dchi = zero
  ddchi = zero

  t = zero
  if (lpr) write(*,*)'*********** S T A R T I N G   T I M E   L O O P ************'
  if (verbose > 1) write(69,*) &
        '*********** S T A R T I N G   T I M E   L O O P ************'

  iclockdump = tick()
  call dump_stuff(0, iseismo, istrain, isnap, disp,velo,chi,dchi,ddchi,t)
  iclockdump = tick(id=iddump, since=iclockdump)

  do iter=1, niter

     t = t + deltat
     call runtime_info(iter,disp,chi)

     ! ::::::::::::::::::::::::: ACTUAL SYMPLECTIC SOLVER :::::::::::::::::::::::::

     ! compute external force/source time function at 4 time intervals coeff(1:4)
     subdt = t - deltat + coeff
     call compute_stf_t(nstages, subdt, stf_symp)

     do i = 1, nstages  ! substages

        chi = chi + dchi * coefd(i)

        disp(:,:,:,1) = disp(:,:,:,1) + velo(:,:,:,1) * coefd(i)
        if (src_type(1) /= 'monopole') &
           disp(:,:,:,2) = disp(:,:,:,2) + velo(:,:,:,2) * coefd(i)
        disp(:,:,:,3) = disp(:,:,:,3) + velo(:,:,:,3) * coefd(i)

        if (src_type(1) /= 'monopole') &
             call apply_axis_mask_scal(chi, nel_fluid, ax_el_fluid,naxel_fluid)

        iclockstiff = tick()
        call glob_fluid_stiffness(ddchi, chi)
        iclockstiff = tick(id=idstiff, since=iclockstiff)

        call bdry_copy2fluid(ddchi, disp)

        if (src_type(1) /= 'monopole') &
             call apply_axis_mask_scal(ddchi, nel_fluid, ax_el_fluid, naxel_fluid)

        iclockcomm = tick()
        call pdistsum_fluid(ddchi, phase=1)
        iclockcomm = tick(id=idcomm, since=iclockcomm)

        iclockstiff = tick()
        select case (src_type(1))
           case ('monopole')
              call apply_axis_mask_onecomp(disp,nel_solid, ax_el_solid,naxel_solid)
              call glob_stiffness_mono(acc,disp)

              if (anel_true) then
                 iclockanelst = tick()
                 call glob_anel_stiffness_mono(acc, memory_var, memory_var_cg4, &
                                               att_coarse_grained)
                 iclockanelst = tick(id=idanelst, since=iclockanelst)
              endif

           case ('dipole')
              call apply_axis_mask_twocomp(disp,nel_solid, ax_el_solid,naxel_solid)
              call glob_stiffness_di(acc,disp)

              if (anel_true) then
                 iclockanelst = tick()
                 call glob_anel_stiffness_di(acc, memory_var, memory_var_cg4, &
                                               att_coarse_grained)
                 iclockanelst = tick(id=idanelst, since=iclockanelst)
              endif

           case ('quadpole')
              call apply_axis_mask_threecomp(disp,nel_solid, ax_el_solid,naxel_solid)
              call glob_stiffness_quad(acc,disp)

              if (anel_true) then
                 iclockanelst = tick()
                 call glob_anel_stiffness_quad(acc, memory_var, memory_var_cg4, &
                                               att_coarse_grained)
                 iclockanelst = tick(id=idanelst, since=iclockanelst)
              endif
        end select
        iclockstiff = tick(id=idstiff, since=iclockstiff)

        iclockcomm = tick()
        call pdistsum_fluid(ddchi, phase=2)
        iclockcomm = tick(id=idcomm, since=iclockcomm)

        ddchi = - ddchi * inv_mass_fluid

        call bdry_copy2solid(acc, ddchi)
        select case (src_type(1))
           case ('monopole')
              call apply_axis_mask_onecomp(acc,nel_solid, ax_el_solid,naxel_solid)

           case ('dipole')
              call apply_axis_mask_twocomp(acc,nel_solid, ax_el_solid,naxel_solid)

           case ('quadpole')
              call apply_axis_mask_threecomp(acc,nel_solid, ax_el_solid,naxel_solid)
        end select

        iclockcomm = tick()
        call pdistsum_solid(acc, phase=1)
        iclockcomm = tick(id=idcomm, since=iclockcomm)

        dchi = dchi + coefv(i) * ddchi

        iclockcomm = tick()
        call pdistsum_solid(acc, phase=2)
        iclockcomm = tick(id=idcomm, since=iclockcomm)

        call add_source(acc, real(stf_symp(i), kind=realkind))

        velo(:,:,:,1) = velo(:,:,:,1) - acc(:,:,:,1) * coefv(i) * inv_mass_rho

        if (src_type(1) /= 'monopole') &
           velo(:,:,:,2) = velo(:,:,:,2) - acc(:,:,:,2) * coefv(i) * inv_mass_rho

        if (src_type(1) == 'dipole') then !factor 2 b/c inv_rho has 1/2 embedded
           velo(:,:,:,3) = velo(:,:,:,3) - two * acc(:,:,:,3) * coefv(i) * inv_mass_rho
        else
           velo(:,:,:,3) = velo(:,:,:,3) - acc(:,:,:,3) * coefv(i) * inv_mass_rho
        endif

     enddo ! ... nstages substages


     chi = chi + dchi * coefd(nstages+1)

     disp(:,:,:,1) = disp(:,:,:,1) + velo(:,:,:,1) * coefd(nstages+1)
     if (src_type(1) /= 'monopole') &
        disp(:,:,:,2) = disp(:,:,:,2) + velo(:,:,:,2) * coefd(nstages+1)
     disp(:,:,:,3) = disp(:,:,:,3) + velo(:,:,:,3) * coefd(nstages+1)

     if (anel_true) then
        iclockanelts = tick()
        call time_step_memvars(memory_var, memory_var_cg4, disp, att_coarse_grained)
        iclockanelts = tick(id=idanelts, since=iclockanelts)
     endif

     ! ::::::::::::::::::::::::: END SYMPLECTIC SOLVER ::::::::::::::::::::::::::

     iclockdump = tick()
     call dump_stuff(iter, iseismo, istrain, isnap, disp,velo,chi,dchi,ddchi,t)
     iclockdump = tick(id=iddump, since=iclockdump)

  enddo ! time loop

end subroutine symplectic_time_loop
!-----------------------------------------------------------------------------------------

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     E N D   O F   T I M E   E X T R A P O L A T I O N   R O U T I N E S
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!-----------------------------------------------------------------------------------------
subroutine symplectic_coefficients(coefd,coeff,coefv)

  use commun, only: barrier,pend

  real(kind=dp), allocatable, dimension(:), intent(out) :: coefd, coeff, coefv
  real(kind=dp), allocatable, dimension(:) :: g
  real(kind=dp) :: zeta_symp, iota_symp, kappa_symp
  real(kind=dp) :: rho, theta, nu, lambda
  integer       :: Q, n, i
  real          :: B, C

  select case (time_scheme)

  !444444444444444444444444444444444444444444444444444444444444444444444444444
  case('symplec4') ! position extended Forest-Ruth like
                   ! (Omelyan, Mryglod and Folk, 2002)
     nstages = 4
     allocate(coefd(nstages+1), coeff(nstages), coefv(nstages))

     if (lpr) then
        write(*,*)
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*) 'TTTT  4th-order symplectic PEFRL scheme  TTTTTTTTTTTTTTT'
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*)
     endif

     ! symplectic parameters
     zeta_symp  = +0.1786178958448091
     iota_symp  = -0.2123418310626054
     kappa_symp = -0.06626458266981849

     ! extrapolation coefficients
     coefd(1:nstages+1) = (/zeta_symp, &
                            kappa_symp, &
                            dble(1.d0-2.d0*(zeta_symp+kappa_symp)), &
                            kappa_symp, &
                            zeta_symp/)

     coefv(1:nstages) = (/dble(half - iota_symp), &
                          iota_symp, &
                          iota_symp, &
                          dble(half - iota_symp)/)

     Q = 4
     B = 1 / 12500.d0
     C = 2.97633; ! empirical


  case('ML_SO4m5')  ! Order 4, S, m=5 in Table 2 of McLachlan (1995), preferred

     nstages=5
     allocate(coefd(nstages+1), coeff(nstages), coefv(nstages))

     if (lpr) then
        write(*,*)
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*) 'TTTT  4th-order symplectic McLachlan scheme  TTTTTTTTTTT'
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*)
     endif

     rho = (14.d0 - dsqrt(19.d0)) / 108.d0;
     theta = (20.d0 - 7.d0 * dsqrt(19.d0)) / 108.d0;
     nu = 2.d0 / 5.d0;
     lambda  = -1.d0 / 10.d0;

     coefd(1) = rho
     coefd(2) = theta
     coefd(3) = 1.d0 / 2.d0 - rho - theta
     coefd(4) = 1.d0 / 2.d0 - rho - theta
     coefd(5) = theta
     coefd(6) = rho

     coefv(1) = nu
     coefv(2) = lambda
     coefv(3) = 1.d0 - 2.d0 * (nu + lambda)
     coefv(4) = lambda
     coefv(5) = nu

     Q = 4
     B = 1.49e-05
     C = 3.035

  !6666666666666666666666666666666666666666666666666666666666666666666666666
  case('ML_SO6m7') ! best order 6 so far!
                   ! other order 6 are not better than the best order 4

     nstages = 7
     allocate(coefd(nstages+1), coeff(nstages), coefv(nstages))

     if (lpr) then
        write(*,*)
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*) 'TTTT  6th-order symplectic McLachlan scheme  TTTTTTTTTTT'
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*)
     endif

     coefd(1) = -1.01308797891717472981
     coefd(2) =  1.18742957373254270702
     coefd(3) = -0.01833585209646059034
     coefd(4) =  0.34399425728109261313
     do i=5, 8
        coefd(i) = coefd(nstages+2-i)
     enddo

     coefv(1) =  0.00016600692650009894
     coefv(2) = -0.37962421426377360608
     coefv(3) =  0.68913741185181063674
     coefv(4) =  0.38064159097092574080
     do i=5, 7
        coefv(i) = coefv(nstages+1-i)
     enddo

     Q = 6
     B = 1.3e-06
     C = 3.067

  !888888888888888888888888888888888888888888888888888888888888888888888888
  case('KL_O8m17') ! Kahan and Li (1997), improved on McLachlan (1995)

     n = 8
     nstages = 2 * n + 1
     allocate(g(n), coefd(nstages+1), coeff(nstages), coefv(nstages))

     if (lpr) then
        write(*,*)
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*) 'TTTT  8th-order symplectic Kahan/Li scheme  TTTTTTTTTTTT'
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*)
     endif

     g(1) =  0.13020248308889008088
     g(2) =  0.56116298177510838456
     g(3) = -0.38947496264484728641
     g(4) =  0.15884190655515560090
     g(5) = -0.39590389413323757734
     g(6) =  0.18453964097831570709
     g(7) =  0.25837438768632204729
     g(8) =  0.29501172360931029887

     call SS_scheme(n,coefd,coefv,g)

     Q = 8
     B = -100000.
     C = 3.

  !101010101010101010101010101010101010101010101010101010101010101010101010
  case('SS_35o10') ! Sofroniou and Spaletta (2004), best order 10

     n = 17
     nstages = 2 * n + 1
     allocate(g(n), coefd(nstages+1), coeff(nstages), coefv(nstages))

     if (lpr) then
        write(*,*)
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*) 'TTTT  10th-order symplectic Sofroniou/Spaletta scheme TT'
        write(*,*) 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        write(*,*)
     endif

    g(1)  =  0.078795722521686419263907679337684
    g(2)  =  0.31309610341510852776481247192647
    g(3)  =  0.027918383235078066109520273275299
    g(4)  = -0.22959284159390709415121339679655
    g(5)  =  0.13096206107716486317465685927961
    g(6)  = -0.26973340565451071434460973222411
    g(7)  =  0.074973343155891435666137105641410
    g(8)  =  0.11199342399981020488957508073640
    g(9)  =  0.36613344954622675119314812353150
    g(10) = -0.39910563013603589787862981058340
    g(11) =  0.10308739852747107731580277001372
    g(12) =  0.41143087395589023782070411897608
    g(13) = -0.0048663605831352617621956593099771
    g(14) = -0.39203335370863990644808193642610
    g(15) =  0.051942502962449647037182904015976
    g(16) =  0.050665090759924496335874344156866
    g(17) =  0.049674370639729879054568800279461

    call SS_scheme(n,coefd,coefv,g)

    Q = 10
    B = 4.58e-10
    C = 5.973

  case default
     write(*,*) procstrg, 'reporting ERROR ::::::::::'
     write(*,*) procstrg, time_scheme,'Time scheme unknown'; call pend; stop

  end select

  coefd = coefd * deltat
  coefv = coefv * deltat

  do i=1, nstages
     coeff(i) = sum(coefd(1:i))
  enddo

  call barrier
  if (mynum == 0) then
     write(*,*)'  :::::::::::::::: Symplectic coefficients :::::::::::::::::'
     write(*,*)'   order,stages:', Q, nstages
     write(*,*)'   dispersion error coeff:', B
     write(*,*)'   CFL factor (wrt Newmark deltat):', C/2.
     do i=1, nstages
        write(*,12) i, coefd(i), coefv(i), coeff(i)
     enddo
     write(*,12)nstages+1,coefd(nstages+1)
     write(*,*)'  :::::::::::::: End Symplectic coefficients :::::::::::::::'
     write(*,*)
  endif
  call barrier

12 format('   ',i3,' coeffd,coeffv,sub_dt:',3(1pe10.2))

end subroutine symplectic_coefficients
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> coefficients for symmetric compositions of symmetric methods
subroutine SS_scheme(n,a,b,g)

  integer, intent(in) :: n
  real(kind=dp), intent(in)  :: g(n)
  real(kind=dp), intent(out) :: a(nstages+1), b(nstages)
  integer :: i

  a(1) = g(1) / 2.d0
  a(2:n) = (g(1:n-1) + g(2:n)) / 2.d0
  a(n+1)= 1 / 2 - sum(a(1:n))
  do i=n+2, 2*n+2
     a(i) = a(2*n+3-i)
  enddo

  b(1:n) = g(1:n)
  b(n+1)= 1. - 2.d0 * sum(g)
  do i=n+2,2*n+1
     b(i) = b(2*n+2-i)
  enddo

end subroutine SS_scheme
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Print time step, time, min/max displacement and potential values
!! and stop the simulation if displacements blow up beyond acceptable...
subroutine runtime_info(iter, disp, chi)

  use commun, only: pend
  !use data_mesh

  integer, intent(in)             :: iter
  real(kind=realkind), intent(in) :: disp(0:,0:,:,:)
  real(kind=realkind), intent(in) :: chi(0:,0:,:)
  integer                         :: iblow(4), check_disp, check_iter, time_stamp
  character(len=4)                :: appistamp

  check_iter = 100 ! printing time/percentage done every check_iter-th time step.
  check_disp = 200 ! printing min/max displ. every check_disp-th time step.
  time_stamp = floor(real(niter)/100.) ! a file every time_stamp-th timestep.

  ! Stdout time step/percentage announcements
  if (lpr .and. mod(iter,check_iter) == 0 ) then
     write(*,13) iter, t, real(iter) / real(niter) * 100.
     call flush(6)
13   format('  time step:',i6,'; t=',f8.2,' s (',f5.1,'%)')
  endif

  ! Time stamp to file
  if (lpr .and. mod(iter,time_stamp) == 0) then
     call define_io_appendix(appistamp, floor(real(iter)/real(time_stamp)))
     open(unit=110, file='timestamp'//appistamp//'.txt')
     write(110,13) iter, t, real(iter) / real(niter) * 100.
     close(110)
  endif

  ! Check on min/max. displacement/potential values globally
  if (verbose > 1) then
     if ( mod(iter,check_disp) == 0 ) then
        if (iter == check_disp) then
           write(69,14) 'time', 'absmax(us)', 'absmax(up)', 'absmax(uz)', 'absmax(chi)'
        endif
        write(69,15) t, maxval(abs(disp(:,:,:,1))), maxval(abs(disp(:,:,:,2))), &
                     maxval(abs(disp(:,:,:,3))), maxval(abs(chi))
     endif
  endif
14 format(a7,4(a13))
15 format(f7.1,4(1pe12.3))

  ! Stop simulation if displacement exceeds source magnitude
  if ( maxval(abs(disp(1,1,:,:))) > 10*abs(magnitude) ) then
     write(*,*) procstrg,'!!!!!!!!!!!!!!! DISPLACEMENTS BLEW UP !!!!!!!!!!!!!!!'
     write(*,*) procstrg,'  Time step & time:', iter, t
     write(*,*) procstrg,'  Proc. num, displ value',mynum,maxval(abs(disp))
     iblow(1:4) = maxloc(abs(disp))
     write(*,*) procstrg,'iel,comp   :',iblow(3),iblow(4)
     write(*,*) procstrg,'elem r, th :',mean_rad_colat_solid(iblow(3),1), &
                                       mean_rad_colat_solid(iblow(3),2)
     write(*,*) procstrg,'ipol,jpol  :',iblow(1)-1,iblow(2)-1
     write(*,*) procstrg,'axis       :',axis_solid(iblow(3))
     write(*,*) procstrg,''
     stop
  endif

end subroutine runtime_info
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Add source term inside source elements only if source time function non-zero
!! and I have the source.
pure subroutine add_source(acc1, stf1)

  real(kind=realkind), intent(in)    :: stf1
  real(kind=realkind), intent(inout) :: acc1(0:,0:,:,:)
  integer                            :: iel,i

  i = 0
  if ( stf1 /= zero) then
     do iel=1, nelsrc
        i = i + 1
        acc1(:,:,ielsrc(iel),:) = acc1(:,:,ielsrc(iel),:) &
             - source_term_el(:,:,i,:) * stf1
     enddo
  endif

end subroutine add_source
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Includes all output action done during the time loop such as
!! various receiver definitions, wavefield snapshots, velocity field & strain
!! tensor for 3-D kernels
subroutine dump_stuff(iter, iseismo, istrain, isnap, &
                      disp, velo, chi, dchi, ddchi, time)

  use data_io
  use data_mesh
  use wavefields_io
  use attenuation, only: n_sls_attenuation, dump_memory_vars
  use nc_routines, only: nc_rec_checkpoint, nc_dump_strain

  integer, intent(in)            :: iter
  integer, intent(inout)         :: iseismo, istrain, isnap
  real(kind=dp),intent(in)       :: time

  real(kind=realkind),intent(in) :: disp(0:, 0:, :, :)
  real(kind=realkind),intent(in) :: velo(0:, 0:, :, :)
  real(kind=realkind),intent(in) :: chi(0:,  0:, :)
  real(kind=realkind),intent(in) :: dchi(0:, 0:, :)
  real(kind=realkind),intent(in) :: ddchi(0:,0:, :)

  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  !^-^-^-^-^-^- Time series^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^

  if (mod(iter,seis_it) == 0) then

     iseismo = iseismo + 1
     if (use_netcdf) then
        call nc_compute_recfile_seis_bare(disp, iseismo)
     else
        call compute_recfile_seis_bare(disp)
     endif

     ! Generic synthetics at hypo-/epicenter, equator, antipode (including time)
     call compute_hyp_epi_equ_anti(real(time, kind=dp), disp)

  endif

  if ((mod(iter, check_it) == 0) .and. (iter > 0)) then
     if (checkpointing .and. use_netcdf) then
        call nc_rec_checkpoint()
     endif
  endif

  ! Compute kinetic and potential energy globally
  if (dump_energy) call energy(disp, velo, dchi, ddchi)

  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^
  !^-^-^-^-^-^ Wavefield snapshots-^-^^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^
  if (dump_vtk) then
    if (mod(iter, snap_it) == 0) then
       isnap = isnap + 1
       if (lpr) then
          write(*,*)
          write(*,*) 'Writing global snap to file: ', isnap
          write(*,*)
       endif
       call glob_snapshot_midpoint(disp, chi, ibeg, iend, jbeg, jend, isnap)
     endif
  endif

  if (dump_xdmf) then
    if (mod(iter, snap_it) == 0) then
        if (.not. (dump_vtk)) isnap=isnap+1
        if (lpr) then
           write(*,*)
           write(*,*)'Writing global xdmf snap to file:',isnap
           write(*,*)
        endif
        call glob_snapshot_xdmf(disp, chi, time, isnap)
     endif
  endif

  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  ! Velocity field and strain tensor wavefields for 3-D kernels^-^-^-^-^
  !^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^^-^-^-^-^-^-^-^-^-^-^-^
  !
  ! At this point, we offer two end-member versions of dumping the fields
  ! to eventually calculate waveform kernels.
  !
  ! The FIRST one ('displ_only') dumps a minimal amount and requires extensive
  ! post-processing when calculating the kernels, but optimizes the SEM
  ! simulation in terms of memory, storage amount and CPU time.
  ! Note that this method IS only POSSIBLE IF ENTIRE SEM MESH IS DUMPED!!!
  !
  ! The SECOND one ('fullfields') computes the entire strain tensor and velocity
  ! field on-the-fly, resulting in more output (9 rather than 6 fields),
  ! more memory and CPU time during the SEM, but no post-processing necessary.
  ! Any kind of spatial distribution can be dumped, meaning in the long run
  ! this should be the more effective choice.
  !
  ! Possible cases in between these dumpsters are considered but not yet
  ! implemented (e.g. dumping 1/s, inverse fluid density, but compute derivatives
  ! on-the-fly to warrant grid flexibility).

  if (dump_wavefields) then

    if (mod(iter,strain_it) == 0) then

      ! dump displacement and velocity in each surface element
      ! for netcdf people set .true. in inparam to use it instead of the standard
      ! the update of the strain has to preceed the call to the function.
      ! It starts from
      istrain = istrain + 1

      call compute_surfelem(disp, velo)

      select case (trim(dump_type))

        case ('displ_only')
          ! Only dump the 3-comp displacement in solid and fluid.
          ! Minimal permanent storage, minimal run-time memory, minimal CPU time,
          ! but extensive post-processing (need to compute strain tensor and
          ! time derivatives, if needed).
             call dump_disp_global(disp, chi, istrain)       ! displacement globally

        case ('strain_only')
          ! Compute strain tensor on-the-fly here and dump the 6 components.
          ! Also compute corresponding fields in the fluid.
            call compute_strain(disp, chi, istrain)    ! strain globally

        case ('displ_velo')
          ! Only dump the 3-comp displacement and velocity fields in solid
          ! and potential & its derivative in fluid.
          ! Minimal permanent storage, minimal run-time memory, minimal CPU time,
          ! but extensive post-processing (need to compute strain tensor).
             call dump_disp(disp, chi, istrain)       ! displacement in solid, chi in fluid
             call dump_velo_dchi(velo, dchi, istrain) ! velocity in solid, dchi in fluid

        case ('fullfields')
          ! Compute strain tensor on-the-fly here and dump the 6 components.
          ! Also compute corresponding fields in the fluid.
          ! Maximal permanent storage, maximal run-time memory, maximal CPU time,
          ! but no post-processeing necessary as these are the fields that
          ! constitute density and elastic kernels.
            call compute_strain(disp, chi, istrain)    ! strain globally
            call dump_velo_global(velo, dchi, istrain) ! velocity globally

!!!!!!  SB
        case ('coupling') ! VM VM dump veloc and stress
           !call compute_strain_cp(disp,velo,chi)  !! VM VM
           !write(*,*) 'DUMP COUPLING FILES'
           call compute_stress_cp(disp,velo,chi,istrain) !! VM VM
           !!     call dump_velo_global_cp(velo,dchi)  !!  _cp means coupling
        case ('coupling_box')
           call compute_stress_cp(disp,velo,chi,istrain) !! VM VM

!!!!! END SB
        end select

        !> Check, whether it is time to dump the buffer variables to disk and if so,
        !! do so.
        if (use_netcdf) call nc_dump_strain(istrain)

    endif ! dumping interval strain_it

endif   ! dump_wavefields?

end subroutine dump_stuff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Compute the full, global strain tensor on-the-fly. Each of 6 (monopole: 4)
!! components is stored separately for solid and fluid domains respectively.
!! The dipole case is transfered to the (s,phi,z) system here.
!!
!! Dumping Ekk, E11, E22, E13, E23, and E12, this has the advantage
!! that if only lambda/bulk sound speed are of interest, then only a
!! scalar Ekk needs to be loaded. Be aware that diuj in the variable names does
!! NOT stand for partial derivatives, but rather the ij component of the
!! strain.
subroutine compute_strain(u, chi, istrain)

  use data_pointwise, only: inv_rho_fluid
  use data_source, only: src_type
  use pointwise_derivatives, only: axisym_gradient_fluid_add
  use pointwise_derivatives, only: axisym_gradient_fluid
  use pointwise_derivatives, only: axisym_gradient_solid_add
  use pointwise_derivatives, only: axisym_gradient_solid
  use pointwise_derivatives, only: f_over_s_solid
  use pointwise_derivatives, only: f_over_s_fluid
  use wavefields_io, only: dump_field_1d

  use data_mesh

  real(kind=realkind), intent(in) :: u(0:,0:,:,:)
  real(kind=realkind), intent(in) :: chi(0:,0:,:)
  integer,             intent(in) :: istrain

  real(kind=realkind)             :: grad_sol(0:npol,0:npol,nel_solid,2)
  real(kind=realkind)             :: buff_solid(0:npol,0:npol,nel_solid)
  real(kind=realkind)             :: usz_fluid(0:npol,0:npol,nel_fluid,2)
  real(kind=realkind)             :: grad_flu(0:npol,0:npol,nel_fluid,2)
  character(len=5)                :: appisnap
  real(kind=realkind), parameter  :: two_rk = real(2, kind=realkind)

  call define_io_appendix(appisnap, istrain)

  ! SSSSSSS Solid region SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  ! s,z components, identical for all source types..........................
  if (src_type(1) == 'dipole') then
     call axisym_gradient_solid(u(:,:,:,1) + u(:,:,:,2), grad_sol)
  else
     call axisym_gradient_solid(u(:,:,:,1), grad_sol) ! 1: dsus, 2: dzus
  endif

  call dump_field_1d(grad_sol(:,:,:,1), '/strain_dsus_sol', appisnap, nel_solid) !E33

  call axisym_gradient_solid_add(u(:,:,:,3), grad_sol) ! 1:dsuz+dzus,2:dzuz+dsus

  ! calculate entire E31 term: (dsuz+dzus)/2
  grad_sol(:,:,:,1) = grad_sol(:,:,:,1) / two_rk
  call dump_field_1d(grad_sol(:,:,:,1), '/strain_dsuz_sol', appisnap, nel_solid) !E31

  ! Components involving phi....................................................
  if (src_type(1) == 'monopole') then
     buff_solid = f_over_s_solid(u(:,:,:,1)) ! us/s
     call dump_field_1d(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !E22
     call dump_field_1d(buff_solid + grad_sol(:,:,:,2), '/straintrace_sol', appisnap, &
                        nel_solid) !Ekk = dzuz + dsus + us/s


  else if (src_type(1) == 'dipole') then
     buff_solid = two_rk * f_over_s_solid(u(:,:,:,2))   ! 2 u-/s
     call dump_field_1d(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !E22
     call dump_field_1d(buff_solid + grad_sol(:,:,:,2), '/straintrace_sol', appisnap, &
                        nel_solid) !Ekk = dzuz + dsus + 2 u-/s

     call axisym_gradient_solid(u(:,:,:,1) - u(:,:,:,2), grad_sol) !1:dsup,2:dzup

     call dump_field_1d(- f_over_s_solid(u(:,:,:,2)) - grad_sol(:,:,:,1) / two_rk, &
                        '/strain_dsup_sol', appisnap, nel_solid) !E12 = - 1/2 (dsup + u-/s)

     call dump_field_1d(- (f_over_s_solid(u(:,:,:,3)) +  grad_sol(:,:,:,2)) / two_rk, &
                        '/strain_dzup_sol', appisnap, nel_solid) !E23 = - 1/2 ( uz/s + dzup )

  else if (src_type(1) == 'quadpole') then
     buff_solid = f_over_s_solid(u(:,:,:,1) - two_rk * u(:,:,:,2)) ! us/s - 2 up/s
     call dump_field_1d(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !E22
     call dump_field_1d(buff_solid + grad_sol(:,:,:,2), & !Ekk = us/s - 2 up/s + dzuz + dsus
                        '/straintrace_sol', appisnap, nel_solid)

     call axisym_gradient_solid(u(:,:,:,2), grad_sol) ! 1: dsup, 2: dzup

     call dump_field_1d(- f_over_s_solid(u(:,:,:,1) + u(:,:,:,2) / two_rk) &
                            - grad_sol(:,:,:,1) / two_rk, &
                        '/strain_dsup_sol', appisnap, nel_solid) !E12 = -us/s + up/2 - dsup/2

     call dump_field_1d(- f_over_s_solid(u(:,:,:,3)) - grad_sol(:,:,:,2) / two_rk, &
                        '/strain_dzup_sol',appisnap, nel_solid) !E23 = -uz/s - dzup / 2
  endif

  ! FFFFFF Fluid region FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  !
  ! Fluid-region strain tensor is computed just like in the solid but for
  ! displacement components ds(chi), dz(chi).

  if (have_fluid) then

     ! construct displacements in the fluid
     call axisym_gradient_fluid(chi, usz_fluid)
     usz_fluid(:,:,:,1) = usz_fluid(:,:,:,1) * inv_rho_fluid
     usz_fluid(:,:,:,2) = usz_fluid(:,:,:,2) * inv_rho_fluid

     ! gradient of s component
     call axisym_gradient_fluid(usz_fluid(:,:,:,1), grad_flu)   ! 1:dsus, 2:dzus

     call dump_field_1d(grad_flu(:,:,:,1), '/strain_dsus_flu', appisnap, nel_fluid) ! E11

     ! gradient of z component added to s-comp gradient for strain trace and E13
     call axisym_gradient_fluid_add(usz_fluid(:,:,:,2), grad_flu)   !1:dsuz+dzus
                                                                    !2:dzuz+dsus

     ! calculate entire E31 term: (dsuz+dzus)/2
     grad_flu(:,:,:,1) = grad_flu(:,:,:,1) / two_rk
     call dump_field_1d(grad_flu(:,:,:,1), '/strain_dsuz_flu', appisnap, nel_fluid) ! E31

     ! Components involving phi................................................

     if (src_type(1) == 'monopole') then
        ! Calculate us/s and straintrace
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)), '/strain_dpup_flu', &
                           appisnap, nel_fluid) ! E22
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)) + grad_flu(:,:,:,2), &
                           '/straintrace_flu', appisnap, nel_fluid) ! Ekk

     else if (src_type(1) == 'dipole') then
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)), &
                           '/strain_dpup_flu', appisnap, nel_fluid)  !E22
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)) &
                            + grad_flu(:,:,:,2), &
                            '/straintrace_flu', appisnap, nel_fluid)  !Ekk

        call dump_field_1d((- f_over_s_fluid(usz_fluid(:,:,:,1))) / two_rk, &
                           '/strain_dsup_flu', appisnap, nel_fluid)   ! E12

        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,2)) / two_rk, &
                           '/strain_dzup_flu', appisnap, nel_fluid)  ! E23

     else if (src_type(1) == 'quadpole') then
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)), &  !E22
                           '/strain_dpup_flu', appisnap, nel_fluid)
        call dump_field_1d(f_over_s_fluid(usz_fluid(:,:,:,1)) &  !Ekk
                            + grad_flu(:,:,:,2), &
                           '/straintrace_flu', appisnap, nel_fluid)

        call dump_field_1d((- f_over_s_fluid(usz_fluid(:,:,:,1))), &
                           '/strain_dsup_flu', appisnap, nel_fluid)   !E12


        call dump_field_1d( - f_over_s_fluid(usz_fluid(:,:,:,2)), &
                           '/strain_dzup_flu', appisnap, nel_fluid)   !E23
     endif   !src_type

  endif   ! have_fluid

end subroutine compute_strain
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Computes the kinetic and potential/elastic/stored energy in the solid and
!! fluid subdomains separately. This involves one additional evaluation of
!! the stiffness system (for the velocity vector rather than displacement)
!! in both domains, but additional cost is rather low if only performed every
!! 5th time step such as the default for the time being.
!! Although non-applicable for accuracy measures of the numerical
!! approximation, the preserved total energy over time
!! (only sufficiently after the source time function is active!) is a useful
!! check on the consistency of the used time scheme and also reveals whether
!! the system leaks, i.e. whether all boundary conditions are correct.
subroutine energy(disp1,vel,dchi1,ddchi)

  use data_source, only: src_type
  use data_matr, only: unassem_mass_rho_solid, unassem_mass_lam_fluid
  use stiffness_mono
  use stiffness_di
  use stiffness_quad
  use stiffness_fluid
  use apply_masks
  use commun
  use data_mesh

  real(kind=realkind), dimension(0:,0:,:,:),intent(in) :: disp1
  real(kind=realkind), dimension(0:,0:,:,:),intent(in) :: vel
  real(kind=realkind), dimension(0:,0:,:),intent(in)   :: dchi1
  real(kind=realkind), dimension(0:,0:,:),intent(in)   :: ddchi

  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: disp
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: dchi
  real(kind=realkind), dimension(0:npol,0:npol,nel_solid,3) :: stiff
  real(kind=realkind), dimension(0:npol,0:npol,nel_fluid)   :: stiff_flu
  real(kind=realkind) :: ekin_sol, epot_sol, ekin_flu, epot_flu

  ekin_sol = zero
  epot_sol = zero
  ekin_flu = zero
  epot_flu = zero
  stiff = zero
  stiff_flu=zero

  disp = disp1
  dchi = dchi1

  ! ssssss SOLID REGION sssssssssssssssssssssssssssssssssssssssssssssssssssss

  ! Compute potential/stored/elastic energy in solid region.
  select case (src_type(1))
     case ('monopole')
           call apply_axis_mask_onecomp(disp,nel_solid, ax_el_solid,naxel_solid)
           call glob_stiffness_mono(stiff,disp)
           call apply_axis_mask_onecomp(stiff,nel_solid, ax_el_solid,naxel_solid)
     case ('dipole')
           call apply_axis_mask_twocomp(disp,nel_solid, ax_el_solid,naxel_solid)
           call glob_stiffness_di(stiff,disp)
           call apply_axis_mask_twocomp(stiff,nel_solid, ax_el_solid,naxel_solid)
     case ('quadpole')
           call apply_axis_mask_threecomp(disp,nel_solid, ax_el_solid,naxel_solid)
           call glob_stiffness_quad(stiff,disp)
           call apply_axis_mask_threecomp(stiff,nel_solid, ax_el_solid,naxel_solid)
  end select

  stiff = stiff * disp
  epot_sol = sum(stiff)
  epot_sol = two*pi*psum(epot_sol)

  ! Compute kinetic energy in solid region
  stiff(:,:,:,1) = vel(:,:,:,1)**2 * unassem_mass_rho_solid
  stiff(:,:,:,2) = vel(:,:,:,2)**2 * unassem_mass_rho_solid
  stiff(:,:,:,3) = vel(:,:,:,3)**2 * unassem_mass_rho_solid

  if (src_type(1) == 'dipole') stiff(:,:,:,3)=two*stiff(:,:,:,3)

  ekin_sol = sum(stiff)
  ekin_sol = two * pi * psum(ekin_sol)

  ! fffff FLUID REGION ffffffffffffffffffffffffffffffffffffffffffffffffffffff
  if (have_fluid) then

     ! Compute potential/stored/elastic energy in fluid region.
     ! Using the wave equation in the fluid here to avoid spatial derivatives.

     stiff_flu=ddchi**2 * unassem_mass_lam_fluid
     epot_flu=sum(stiff_flu)
     epot_flu=two*pi*psum(epot_flu)

     ! Fluid stiffness acting on derivative of potential
     stiff_flu = zero
     if (src_type(1) /= 'monopole') &
          call apply_axis_mask_scal(dchi,nel_fluid,ax_el_fluid,naxel_fluid)
     call glob_fluid_stiffness(stiff_flu,dchi)
     if (src_type(1) /= 'monopole') &
         call apply_axis_mask_scal(stiff_flu,nel_fluid,ax_el_fluid,naxel_fluid)

     ! Compute kinetic energy in fluid region
     stiff_flu=stiff_flu*dchi
     ekin_flu=sum(stiff_flu)
     ekin_flu=two*pi*psum(ekin_flu)

  endif ! have_fluid

  ! Save into time series, only one processor needs to do this since global.
  if (lpr) then
     if (have_fluid) then
        write(4444,8) t, ekin_sol, epot_sol
        write(4445,8) t, epot_flu, ekin_flu
     endif
     write(4446,9) t, epot_sol + epot_flu, ekin_sol + ekin_flu, &
                   half * (epot_sol + epot_flu + ekin_sol + ekin_flu)
  endif
8 format(3(1pe16.6))
9 format(4(1pe16.6))

end subroutine energy
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> FLUID: solid-fluid boundary term (solid displ.) added to fluid stiffness
!!        minus sign in accordance with the definition of bdry_matr
pure subroutine bdry_copy2fluid(uflu,usol)

  use data_matr, only: bdry_matr, solflubdry_radius
  use data_source, only: src_type


  real(kind=realkind), intent(inout) :: uflu(0:,0:,:)
  real(kind=realkind), intent(in)    :: usol(0:,0:,:,:)
  integer                            :: iel,jpols,jpolf,iels,ielf

  if (src_type(1) == 'dipole') then

     do iel=1,nel_bdry
        jpols = bdry_jpol_solid(iel)
        jpolf = bdry_jpol_fluid(iel)
        iels = bdry_solid_el(iel)
        ielf = bdry_fluid_el(iel)

        uflu(:,jpolf,ielf) = uflu(:,jpolf,ielf) - &
                             bdry_matr(:,iel,1) * &
                             (usol(:,jpols,iels,1) + usol(:,jpols,iels,2) ) - &
                             bdry_matr(:,iel,2) * usol(:,jpols,iels,3)
     enddo

  else

     do iel=1,nel_bdry
        jpols = bdry_jpol_solid(iel)
        jpolf = bdry_jpol_fluid(iel)
        iels = bdry_solid_el(iel)
        ielf = bdry_fluid_el(iel)

        uflu(:,jpolf,ielf) = uflu(:,jpolf,ielf) - &
                             bdry_matr(:,iel,1) * usol(:,jpols,iels,1) - &
                             bdry_matr(:,iel,2) * usol(:,jpols,iels,3)
     enddo

  endif

end subroutine bdry_copy2fluid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> SOLID: solid-fluid boundary term (fluid potential) added to solid stiffness
!!        plus sign in accordance with definition of bdry_matr
pure subroutine bdry_copy2solid(usol,uflu)

  use data_matr, only: bdry_matr, solflubdry_radius
  use data_source, only: src_type


  real(kind=realkind), intent(inout) :: usol(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: uflu(0:,0:,:)
  integer                            :: iel, jpols, jpolf, iels, ielf

  do iel = 1,nel_bdry
     jpols = bdry_jpol_solid(iel)
     jpolf = bdry_jpol_fluid(iel)
     iels = bdry_solid_el(iel)
     ielf = bdry_fluid_el(iel)

     if (src_type(1) == 'dipole') then

        usol(:,jpols,iels,1) = usol(:,jpols,iels,1) &
                                + bdry_matr(:,iel,1) * uflu(:,jpolf,ielf)
        usol(:,jpols,iels,2) = usol(:,jpols,iels,2) &
                                + bdry_matr(:,iel,1) * uflu(:,jpolf,ielf)

     else
        usol(:,jpols,iels,1) = usol(:,jpols,iels,1) &
                                + bdry_matr(:,iel,1) * uflu(:,jpolf,ielf)
     endif

     usol(:,jpols,iels,3) = usol(:,jpols,iels,3) &
                            + bdry_matr(:,iel,2) * uflu(:,jpolf,ielf)

  enddo


end subroutine bdry_copy2solid

!####################################################################################################!
!!
!!------------------------------------------------------------------------------------------!!
!!
!!--------- VM & CD, subroutine(s) for coupling/KH integral methods with specfem -----------!!
!!
!!------------------------------------------------------------------------------------------!!
!!
!####################################################################################################!


subroutine compute_stress_cp(u, v, chi, istrain)


  use data_source, only: src_type
  use pointwise_derivatives, only: axisym_gradient_solid_add
  use pointwise_derivatives, only: axisym_gradient_solid
  use pointwise_derivatives, only: f_over_s_solid
  use wavefields_io, only: dump_field_1d

  use data_mesh
  use coupling_mod

!  use data_pointwise, only: inv_rho_fluid, prefac_inv_s_rho_fluid
!  use pointwise_derivatives, only: axisym_gradient_fluid_add
!  use pointwise_derivatives, only: axisym_gradient_fluid
!  use pointwise_derivatives, only: f_over_s_fluid


  real(kind=realkind), intent(in) :: u(0:,0:,:,:),v(0:,0:,:,:)
  real(kind=realkind), intent(in) :: chi(0:,0:,:)

  real(kind=realkind)             :: strain(0:npol,0:npol,nel_solid,6)
  real(kind=realkind)             :: stress(0:npol,0:npol,nel_solid,6)
  real(kind=realkind)             :: grad_sol(0:npol,0:npol,nel_solid,2)
  real(kind=realkind)             :: grad_sol_save(0:npol,0:npol,nel_solid)
  real(kind=realkind)             :: buff_solid(0:npol,0:npol,nel_solid)
  real(kind=realkind)             :: usz_fluid(0:npol,0:npol,nel_fluid,2)
  real(kind=realkind)             :: up_fluid(0:npol,0:npol,nel_fluid)
  real(kind=realkind)             :: grad_flu(0:npol,0:npol,nel_fluid,2)

  !! CD CD add this = => for the calculation fo derivatives one by one
  real(kind=realkind)             :: D_us(0:npol,0:npol,nel_solid,3) ! = dus/ds, dus/dp, dus/dz
  real(kind=realkind)             :: D_up(0:npol,0:npol,nel_solid,3) ! = dup/ds, dup/dp, dup/dz
  real(kind=realkind)             :: D_uz(0:npol,0:npol,nel_solid,3) ! = duz/ds, duz/dp, duz/dz
  real(kind=realkind)             :: grad_sol_tmp(0:npol,0:npol,nel_solid,2)
  real(kind=realkind)             :: s_tmp(0:npol,0:npol,nel_solid)

  integer                         :: istrain
  character(len=5)                :: appisnap
  real(kind=realkind), parameter  :: two_rk = real(2, kind=realkind)
  logical, parameter              :: have_fluid_cp = .false. !! FOR NOW THE BOX IS ASSUMED TO BE IN SOLID REGION !!

  integer iel0, iel, i, j
  real(kind=realkind) l, m


!!--- Initialization

  strain = 0.
  stress = 0.

  D_us   = 0.
  D_up   = 0.
  D_uz   = 0.

!!--- WARNING !! = => Voigt convention for the strain/stress  = => 1:ss, 2:pp, 3:zz, 4:sp, 5:sz, 6:pz !!


  call define_io_appendix(appisnap,istrain)

! SSSSSSSSSSSSSSSSSSSS Solid region SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  select case (trim(src_type(1)))

!
!---
!

  case ('monopole')

!
!-------------------------------- calculation of strain --------------------------------!
!

    call axisym_gradient_solid(u(:,:,:,1), grad_sol) ! 1: dsus, 2: dzus

    strain(:,:,:,1)=grad_sol(:,:,:,1)

    grad_sol_save=grad_sol(:,:,:,1)
    call axisym_gradient_solid_add(u(:,:,:,3), grad_sol) ! 1:dsuz+dzus,2:dzuz+dsus

    strain(:,:,:,5)=grad_sol(:,:,:,1) / two_rk

    buff_solid = f_over_s_solid(u(:,:,:,1))

    strain(:,:,:,2)= buff_solid

    strain(:,:,:,3)=grad_sol(:,:,:,2) - grad_sol_save

!
!-------------------------------- computation of stress --------------------------------!
!

    do iel=1,nel_solid
      iel0=ielsolid(iel)
      do j=0,npol
        do i=0,npol

          l=lambda_cp(i,j,iel0)
          m=mu_cp(i,j,iel0)

          stress(i,j,iel,1) = (l+2.*m)*(strain(i,j,iel,1) + strain(i,j,iel,2)) -2.*m*strain(i,j,iel,2) + l*strain(i,j,iel,3)
          stress(i,j,iel,2) = (l+2.*m)*(strain(i,j,iel,1) + strain(i,j,iel,2)) -2.*m*strain(i,j,iel,1) + l*strain(i,j,iel,3)
          stress(i,j,iel,3) = l*(strain(i,j,iel,1) + strain(i,j,iel,2)) + (l+2*m)*strain(i,j,iel,3)
          stress(i,j,iel,5) = 2.*m*strain(i,j,iel,5)

        enddo
      enddo
    enddo

!
!-------------- compute derivatives one by one, to store it in case of KH --------------!
!

    if (storage_for_recip_KH_integral) then

      call axisym_gradient_solid(u(:,:,:,1), grad_sol_tmp)

      D_us(:,:,:,1) = grad_sol_tmp(:,:,:,1) !! dus/ds
      D_us(:,:,:,3) = grad_sol_tmp(:,:,:,2) !! dus/dz

      call axisym_gradient_solid(u(:,:,:,3), grad_sol_tmp)

      D_uz(:,:,:,1) = grad_sol_tmp(:,:,:,1) !! duz/ds
      D_uz(:,:,:,3) = grad_sol_tmp(:,:,:,2) !! duz/dz

!---

      D_up(:,:,:,2) = strain(:,:,:,2) !! (dup/dp)/s + us/s

      D_us(:,:,:,2) = 0. !! (dus/dp)/s
      D_up(:,:,:,3) = 0. !! dup/dz
      D_uz(:,:,:,2) = 0. !! (duz/dp)/s
      D_up(:,:,:,1) = 0. !! (dup/ds) - up/s

      call dump_field_1d_cp(D_us(:,:,:,1),'/pgrad_cyl_Dus_1', appisnap, nel_solid)
      call dump_field_1d_cp(D_us(:,:,:,3),'/pgrad_cyl_Dus_3', appisnap, nel_solid)
      call dump_field_1d_cp(D_uz(:,:,:,1),'/pgrad_cyl_Duz_1', appisnap, nel_solid)
      call dump_field_1d_cp(D_uz(:,:,:,3),'/pgrad_cyl_Duz_3', appisnap, nel_solid)
      call dump_field_1d_cp(D_up(:,:,:,2),'/pgrad_cyl_Dup_2', appisnap, nel_solid)

    endif

!
!---------------------------------------------------------------------------------------!
!

    call dump_field_1d_cp(u(:,:,:,1),'/displacement_us', appisnap, nel_solid)
    call dump_field_1d_cp(u(:,:,:,3),'/displacement_uz', appisnap, nel_solid)
    call dump_field_1d_cp(v(:,:,:,1),'/velocityfiel_us', appisnap, nel_solid)
    call dump_field_1d_cp(v(:,:,:,3),'/velocityfiel_uz', appisnap, nel_solid)

    call dump_field_1d_cp(stress(:,:,:,1),'/stress_Sg11_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,2),'/stress_Sg22_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,3),'/stress_Sg33_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,5),'/stress_Sg13_sol', appisnap, nel_solid)
!!     !call dump_field_1d_cp(stress(:,:,:,4),'/stress_Sg12_sol', appisnap, nel_solid)
!!     !call dump_field_1d_cp(stress(:,:,:,6),'/stress_Sg23_sol', appisnap, nel_solid)

!
!---
!

  case ('dipole')

!
!-------------------------------- calculation of strain --------------------------------!
!

    call axisym_gradient_solid(u(:,:,:,1) + u(:,:,:,2), grad_sol)

    strain(:,:,:,1)=grad_sol(:,:,:,1)

    grad_sol_save=grad_sol(:,:,:,1)
    call axisym_gradient_solid_add(u(:,:,:,3), grad_sol) ! 1:dsuz+dzus,2:dzuz+dsus

    strain(:,:,:,5)=grad_sol(:,:,:,1) / two_rk

    buff_solid = two_rk * f_over_s_solid(u(:,:,:,2))

    strain(:,:,:,2)= buff_solid

    strain(:,:,:,3)=grad_sol(:,:,:,2) - grad_sol_save

    call axisym_gradient_solid(u(:,:,:,1) - u(:,:,:,2), grad_sol) !1:dsup,2:dzup

    strain(:,:,:,4) =  f_over_s_solid(u(:,:,:,2)) + grad_sol(:,:,:,1) / two_rk !! VM VM

    strain(:,:,:,6) = ( f_over_s_solid(u(:,:,:,3)) +  grad_sol(:,:,:,2) ) / two_rk !! VM VM corrected a sign

!
!-------------------------------- computation of stress --------------------------------!
!

    do iel=1,nel_solid
      iel0=ielsolid(iel)
      do j=0,npol
        do i=0,npol

          l=lambda_cp(i,j,iel0)
          m=mu_cp(i,j,iel0)

          stress(i,j,iel,1) = (l+2.*m)*(strain(i,j,iel,1) + strain(i,j,iel,2)) -2.*m*strain(i,j,iel,2) + l*strain(i,j,iel,3)
          stress(i,j,iel,2) = (l+2.*m)*(strain(i,j,iel,1) + strain(i,j,iel,2)) -2.*m*strain(i,j,iel,1) + l*strain(i,j,iel,3)
          stress(i,j,iel,3) = l*(strain(i,j,iel,1) + strain(i,j,iel,2)) + (l+2*m)*strain(i,j,iel,3)
          stress(i,j,iel,4) = 2.*m*strain(i,j,iel,4)
          stress(i,j,iel,5) = 2.*m*strain(i,j,iel,5)
          stress(i,j,iel,6) = 2.*m*strain(i,j,iel,6)

        enddo
      enddo
    enddo

!
!-------------- compute derivatives one by one, to store it in case of KH --------------!
!

    if (storage_for_recip_KH_integral) then

      call axisym_gradient_solid(u(:,:,:,1)+u(:,:,:,2), grad_sol_tmp)

      D_us(:,:,:,1) = grad_sol_tmp(:,:,:,1) !! dus/ds
      D_us(:,:,:,3) = grad_sol_tmp(:,:,:,2) !! dus/dz

      call axisym_gradient_solid(u(:,:,:,1)-u(:,:,:,2), grad_sol_tmp)

      D_up(:,:,:,1) = grad_sol_tmp(:,:,:,1) - f_over_s_solid(u(:,:,:,1)-u(:,:,:,2)) !! (dup/ds) - up/s
      D_up(:,:,:,3) = grad_sol_tmp(:,:,:,2) !! dup/dz

      call axisym_gradient_solid(u(:,:,:,3), grad_sol_tmp)

      D_uz(:,:,:,1) = grad_sol_tmp(:,:,:,1) !! duz/ds
      D_uz(:,:,:,3) = grad_sol_tmp(:,:,:,2) !! duz/dz


      do iel=1,nel_solid
        iel0=ielsolid(iel)
        do j=0,npol
          do i=0,npol
            s_tmp(i,j,iel) = scoord(i,j,iel0)
          enddo
        enddo
      enddo

      D_us(:,:,:,2) = f_over_s_solid(u(:,:,:,1) + u(:,:,:,2))
      D_up(:,:,:,2) = strain(:,:,:,2)
      D_uz(:,:,:,2) = f_over_s_solid(u(:,:,:,3))

      call dump_field_1d_cp(D_us(:,:,:,1),'/pgrad_cyl_Dus_1', appisnap, nel_solid)
      call dump_field_1d_cp(D_us(:,:,:,2),'/pgrad_cyl_Dus_2', appisnap, nel_solid)
      call dump_field_1d_cp(D_us(:,:,:,3),'/pgrad_cyl_Dus_3', appisnap, nel_solid)

      call dump_field_1d_cp(D_up(:,:,:,1),'/pgrad_cyl_Dup_1', appisnap, nel_solid)
      call dump_field_1d_cp(D_up(:,:,:,2),'/pgrad_cyl_Dup_2', appisnap, nel_solid)
      call dump_field_1d_cp(D_up(:,:,:,3),'/pgrad_cyl_Dup_3', appisnap, nel_solid)

      call dump_field_1d_cp(D_uz(:,:,:,1),'/pgrad_cyl_Duz_1', appisnap, nel_solid)
      call dump_field_1d_cp(D_uz(:,:,:,2),'/pgrad_cyl_Duz_2', appisnap, nel_solid)
      call dump_field_1d_cp(D_uz(:,:,:,3),'/pgrad_cyl_Duz_3', appisnap, nel_solid)

    endif

!
!---------------------------------------------------------------------------------------!
!

    call dump_field_1d_cp(u(:,:,:,1)+u(:,:,:,2),'/displacement_us', appisnap, nel_solid)
    call dump_field_1d_cp(u(:,:,:,1)-u(:,:,:,2),'/displacement_up', appisnap, nel_solid)
    call dump_field_1d_cp(u(:,:,:,3),'/displacement_uz', appisnap, nel_solid)

    call dump_field_1d_cp(v(:,:,:,1)+v(:,:,:,2),'/velocityfiel_us', appisnap, nel_solid)
    call dump_field_1d_cp(v(:,:,:,1)-v(:,:,:,2),'/velocityfiel_up', appisnap, nel_solid)
    call dump_field_1d_cp(v(:,:,:,3),'/velocityfiel_uz', appisnap, nel_solid)

    call dump_field_1d_cp(stress(:,:,:,1),'/stress_Sg11_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,2),'/stress_Sg22_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,3),'/stress_Sg33_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,4),'/stress_Sg12_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,5),'/stress_Sg13_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,6),'/stress_Sg23_sol', appisnap, nel_solid)

!
!---
!

  case('quadpole')

!
!-------------------------------- calculation of strain --------------------------------!
!

    call axisym_gradient_solid(u(:,:,:,1), grad_sol) ! 1: dsus, 2: dzus

    strain(:,:,:,1)=grad_sol(:,:,:,1)

    grad_sol_save=grad_sol(:,:,:,1)
    call axisym_gradient_solid_add(u(:,:,:,3), grad_sol)

    strain(:,:,:,5)= grad_sol(:,:,:,1) / two_rk

    buff_solid = f_over_s_solid(u(:,:,:,1) - two_rk * u(:,:,:,2))

    strain(:,:,:,2)= buff_solid

    strain(:,:,:,3)=grad_sol(:,:,:,2) - grad_sol_save

    call axisym_gradient_solid(u(:,:,:,2), grad_sol) ! 1: dsup, 2: dzup

    strain(:,:,:,4) = (f_over_s_solid(- u(:,:,:,2) + 2 * u(:,:,:,1)) + grad_sol(:,:,:,1) ) /two_rk
    !f_over_s_solid(u(:,:,:,1) + u(:,:,:,2) / two_rk) - grad_sol(:,:,:,1) / two_rk

    strain(:,:,:,6) =  f_over_s_solid(u(:,:,:,3)) + grad_sol(:,:,:,2) / two_rk

!
!-------------------------------- computation of stress --------------------------------!
!

    do iel=1,nel_solid
      iel0=ielsolid(iel)
      do j=0,npol
        do i=0,npol

          l=lambda_cp(i,j,iel0)
          m=mu_cp(i,j,iel0)

          stress(i,j,iel,1) = (l+2.*m)*(strain(i,j,iel,1) + strain(i,j,iel,2)) -2.*m*strain(i,j,iel,2) + l*strain(i,j,iel,3)
          stress(i,j,iel,2) = (l+2.*m)*(strain(i,j,iel,1) + strain(i,j,iel,2)) -2.*m*strain(i,j,iel,1) + l*strain(i,j,iel,3)
          stress(i,j,iel,3) = l*(strain(i,j,iel,1) + strain(i,j,iel,2)) + (l+2*m)*strain(i,j,iel,3)
          stress(i,j,iel,4) = 2.*m*strain(i,j,iel,4)
          stress(i,j,iel,5) = 2.*m*strain(i,j,iel,5)
          stress(i,j,iel,6) = 2.*m*strain(i,j,iel,6)
        enddo
      enddo
    enddo

!
!---------------------------------------------------------------------------------------!
!

    call dump_field_1d_cp(u(:,:,:,1),'/displacement_us', appisnap, nel_solid)
    call dump_field_1d_cp(u(:,:,:,2),'/displacement_up', appisnap, nel_solid)
    call dump_field_1d_cp(u(:,:,:,3),'/displacement_uz', appisnap, nel_solid)

    call dump_field_1d_cp(V(:,:,:,1),'/velocityfiel_us', appisnap, nel_solid)
    call dump_field_1d_cp(V(:,:,:,2),'/velocityfiel_up', appisnap, nel_solid)
    call dump_field_1d_cp(V(:,:,:,3),'/velocityfiel_uz', appisnap, nel_solid)

    call dump_field_1d_cp(stress(:,:,:,1),'/stress_Sg11_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,2),'/stress_Sg22_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,3),'/stress_Sg33_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,4),'/stress_Sg12_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,5),'/stress_Sg13_sol', appisnap, nel_solid)
    call dump_field_1d_cp(stress(:,:,:,6),'/stress_Sg23_sol', appisnap, nel_solid)

  end select

end subroutine compute_stress_cp

!------------------------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!--- compute_strain_cp not used for the moment --- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------!
!!-!!
!!-!!subroutine compute_strain_cp(u, v, chi, istrain)
!!-!!
!!-!!  use coupling_mod, only: lambda_cp,mu_cp
!!-!!  !use data_pointwise, only: inv_rho_fluid, prefac_inv_s_rho_fluid
!!-!!  use data_source, only: src_type
!!-!!  !use pointwise_derivatives, only: axisym_gradient_fluid_add
!!-!!  !use pointwise_derivatives, only: axisym_gradient_fluid
!!-!!  use pointwise_derivatives, only: axisym_gradient_solid_add
!!-!!  use pointwise_derivatives, only: axisym_gradient_solid
!!-!!  use pointwise_derivatives, only: f_over_s_solid
!!-!!  !use pointwise_derivatives, only: f_over_s_fluid
!!-!!  use wavefields_io, only: dump_field_1d
!!-!!
!!-!!  use data_mesh
!!-!!
!!-!!
!!-!!  real(kind=realkind), intent(in) :: u(0:,0:,:,:),v(0:,0:,:,:)
!!-!!  real(kind=realkind), intent(in) :: chi(0:,0:,:)
!!-!!
!!-!!  real(kind=realkind)             :: grad_sol(0:npol,0:npol,nel_solid,2)
!!-!!  real(kind=realkind)             :: grad_sol_save(0:npol,0:npol,nel_solid)
!!-!!  real(kind=realkind)             :: buff_solid(0:npol,0:npol,nel_solid)
!!-!!  real(kind=realkind)             :: usz_fluid(0:npol,0:npol,nel_fluid,2)
!!-!!  real(kind=realkind)             :: up_fluid(0:npol,0:npol,nel_fluid)
!!-!!  real(kind=realkind)             :: grad_flu(0:npol,0:npol,nel_fluid,2)
!!-!!  character(len=5)                :: appisnap
!!-!!  real(kind=realkind), parameter  :: two_rk = real(2, kind=realkind)
!!-!!  logical, parameter :: have_fluid_cp=.false. !! FOR NOW THE BOX IS ASSUMED TO BE IN SOLID REGION
!!-!!
!!-!!  integer :: istrain
!!-!!
!!-!!  call define_io_appendix(appisnap,istrain)
!!-!!
!!-!!  ! SSSSSSS Solid region SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!!-!!
!!-!!  ! s,z components, identical for all source types..........................
!!-!!  if (src_type(1)=='dipole') then
!!-!!     call axisym_gradient_solid(u(:,:,:,1) + u(:,:,:,2), grad_sol)
!!-!!  else
!!-!!     call axisym_gradient_solid(u(:,:,:,1), grad_sol) ! 1: dsus, 2: dzus
!!-!!  endif
!!-!!  !                                        '/stress_Sg11_sol'
!!-!!  call dump_field_1d_cp(grad_sol(:,:,:,1), '/strain_dsus_sol', appisnap, nel_solid) !E11
!!-!!  grad_sol_save=grad_sol(:,:,:,1)
!!-!!
!!-!!  call axisym_gradient_solid_add(u(:,:,:,3), grad_sol) ! 1:dsuz+dzus,2:dzuz+dsus
!!-!!
!!-!!  ! calculate entire E31 term: (dsuz+dzus)/2
!!-!!  grad_sol(:,:,:,1) = grad_sol(:,:,:,1) / two_rk
!!-!!  call dump_field_1d_cp(grad_sol(:,:,:,1), '/strain_dsuz_sol', appisnap, nel_solid) !E31
!!-!!
!!-!!  ! Components involving phi....................................................
!!-!!  if (src_type(1) == 'monopole') then
!!-!!     buff_solid = f_over_s_solid(u(:,:,:,1))
!!-!!     call dump_field_1d_cp(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !E22
!!-!!     !call dump_field_1d_cp(buff_solid - grad_sol_save, '/strain_dzuz_sol', appisnap, nel_solid) !E33
!!-!!     call dump_field_1d_cp(buff_solid + grad_sol(:,:,:,2), '/straintrace_sol', appisnap, &
!!-!!                        nel_solid) !Ekk
!!-!!
!!-!!     !! VM VM add displ dump
!!-!!     call dump_field_1d_cp(u(:,:,:,1),'/displacement_us', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(u(:,:,:,3),'/displacement_uz', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(v(:,:,:,1),'/velocityfiel_us', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(v(:,:,:,3),'/velocityfiel_uz', appisnap, nel_solid)
!!-!!
!!-!!  else if (src_type(1) == 'dipole') then
!!-!!     buff_solid = two_rk * f_over_s_solid(u(:,:,:,2))
!!-!!     call dump_field_1d_cp(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !E22
!!-!!     !call dump_field_1d_cp(buff_solid - grad_sol_save, '/strain_dzuz_sol', appisnap, nel_solid) !E33
!!-!!     call dump_field_1d_cp(buff_solid + grad_sol(:,:,:,2), '/straintrace_sol', appisnap, &
!!-!!                        nel_solid) ! Ekk
!!-!!
!!-!!     call axisym_gradient_solid(u(:,:,:,1) - u(:,:,:,2), grad_sol) !1:dsup,2:dzup
!!-!!
!!-!!     call dump_field_1d_cp(- f_over_s_solid(u(:,:,:,2)) - grad_sol(:,:,:,1) / two_rk, &
!!-!!                        '/strain_dsup_sol', appisnap, nel_solid) !E12
!!-!!
!!-!!     call dump_field_1d_cp(- (f_over_s_solid(u(:,:,:,3)) +  grad_sol(:,:,:,2)) / two_rk, &
!!-!!                        '/strain_dzup_sol', appisnap, nel_solid) !E23
!!-!!
!!-!!     !! VM VM add displ dump
!!-!!     call dump_field_1d_cp(u(:,:,:,1)+u(:,:,:,2),'/displacement_us', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(u(:,:,:,1)-u(:,:,:,2),'/displacement_up', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(u(:,:,:,3),'/displacement_uz', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(v(:,:,:,1)+v(:,:,:,2),'/velocityfiel_us', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(v(:,:,:,1)-v(:,:,:,2),'/velocityfiel_up', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(v(:,:,:,3),'/velocityfiel_uz', appisnap, nel_solid)
!!-!!
!!-!!  else if (src_type(1) == 'quadpole') then
!!-!!     buff_solid = f_over_s_solid(u(:,:,:,1) - two_rk * u(:,:,:,2))
!!-!!     call dump_field_1d_cp(buff_solid, '/strain_dpup_sol', appisnap, nel_solid) !E22
!!-!!     !call dump_field_1d_cp(buff_solid - grad_sol_save, '/strain_dzuz_sol', appisnap, nel_solid) !E33
!!-!!     call dump_field_1d_cp(buff_solid + grad_sol(:,:,:,2), & !Ekk
!!-!!                        '/straintrace_sol', appisnap, nel_solid)
!!-!!
!!-!!     call axisym_gradient_solid(u(:,:,:,2), grad_sol) ! 1: dsup, 2: dzup
!!-!!
!!-!!     call dump_field_1d_cp(- f_over_s_solid(u(:,:,:,1) + u(:,:,:,2) / two_rk) &
!!-!!                            - grad_sol(:,:,:,1) / two_rk, &
!!-!!                        '/strain_dsup_sol', appisnap, nel_solid) !E12
!!-!!
!!-!!     call dump_field_1d_cp(- f_over_s_solid(u(:,:,:,3)) - grad_sol(:,:,:,2) / two_rk, &
!!-!!                        '/strain_dzup_sol',appisnap, nel_solid) !E23
!!-!!
!!-!!
!!-!!
!!-!!     !! VM VM add displ dump
!!-!!     call dump_field_1d_cp(u(:,:,:,1),'/displacement_us', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(u(:,:,:,2),'/displacement_up', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(u(:,:,:,3),'/displacement_uz', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(V(:,:,:,1),'/velocityfiel_us', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(V(:,:,:,2),'/velocityfiel_up', appisnap, nel_solid)
!!-!!     call dump_field_1d_cp(V(:,:,:,3),'/velocityfiel_uz', appisnap, nel_solid)
!!-!!  endif
!!-!!!!$
!!-!!!!$  ! FFFFFF Fluid region FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
!!-!!!!$  !
!!-!!!!$  ! Fluid-region strain tensor is computed just like in the solid but for
!!-!!!!$  ! displacement components ds(chi), dz(chi).
!!-!!!!$
!!-!!!!$  if (have_fluid_cp) then
!!-!!!!$
!!-!!!!$     ! construct displacements in the fluid
!!-!!!!$     call axisym_gradient_fluid(chi, usz_fluid)
!!-!!!!$     usz_fluid(:,:,:,1) = usz_fluid(:,:,:,1) * inv_rho_fluid
!!-!!!!$     usz_fluid(:,:,:,2) = usz_fluid(:,:,:,2) * inv_rho_fluid
!!-!!!!$
!!-!!!!$     ! gradient of s component
!!-!!!!$     call axisym_gradient_fluid(usz_fluid(:,:,:,1), grad_flu)   ! 1:dsus, 2:dzus
!!-!!!!$
!!-!!!!$     call dump_field_1d_cp(grad_flu(:,:,:,1), '/strain_dsus_flu', appisnap, nel_fluid) ! E11
!!-!!!!$
!!-!!!!$     ! gradient of z component added to s-comp gradient for strain trace and E13
!!-!!!!$     call axisym_gradient_fluid_add(usz_fluid(:,:,:,2), grad_flu)   !1:dsuz+dzus
!!-!!!!$                                                                    !2:dzuz+dsus
!!-!!!!$
!!-!!!!$     ! calculate entire E31 term: (dsuz+dzus)/2
!!-!!!!$     grad_flu(:,:,:,1) = grad_flu(:,:,:,1) / two_rk
!!-!!!!$     call dump_field_1d_cp(grad_flu(:,:,:,1), '/strain_dsuz_flu', appisnap, nel_fluid) ! E31
!!-!!!!$
!!-!!!!$     ! Components involving phi................................................
!!-!!!!$
!!-!!!!$     if (src_type(1) == 'monopole') then
!!-!!!!$        ! Calculate us/s and straintrace
!!-!!!!$        call dump_field_1d_cp(f_over_s_fluid(usz_fluid(:,:,:,1)), '/strain_dpup_flu', &
!!-!!!!$                           appisnap, nel_fluid) ! E22
!!-!!!!$        call dump_field_1d_cp(f_over_s_fluid(usz_fluid(:,:,:,1)) + grad_flu(:,:,:,2), &
!!-!!!!$                           '/straintrace_flu', appisnap, nel_fluid) ! Ekk
!!-!!!!$
!!-!!!!$     else if (src_type(1) == 'dipole') then
!!-!!!!$        up_fluid = prefac_inv_s_rho_fluid * chi
!!-!!!!$        call dump_field_1d_cp(f_over_s_fluid(usz_fluid(:,:,:,1) - up_fluid), &
!!-!!!!$                           '/strain_dpup_flu', appisnap, nel_fluid)  !E22
!!-!!!!$        call dump_field_1d_cp(f_over_s_fluid(usz_fluid(:,:,:,1) - up_fluid) &
!!-!!!!$                            + grad_flu(:,:,:,2), &
!!-!!!!$                            '/straintrace_flu', appisnap, nel_fluid)  !Ekk
!!-!!!!$
!!-!!!!$        ! gradient of phi component
!!-!!!!$        call axisym_gradient_fluid(up_fluid, grad_flu)   ! 1:dsup, 2:dzup
!!-!!!!$
!!-!!!!$        call dump_field_1d_cp((- grad_flu(:,:,:,1) &
!!-!!!!$                            - f_over_s_fluid(usz_fluid(:,:,:,1) &
!!-!!!!$                                - up_fluid)) / two_rk, &
!!-!!!!$                           '/strain_dsup_flu', appisnap, nel_fluid)   ! E12
!!-!!!!$
!!-!!!!$        call dump_field_1d_cp(- (grad_flu(:,:,:,2) - f_over_s_fluid(usz_fluid(:,:,:,2))) &
!!-!!!!$                            / two_rk, '/strain_dzup_flu', appisnap, nel_fluid)  ! E23
!!-!!!!$
!!-!!!!$     else if (src_type(1) == 'quadpole') then
!!-!!!!$        up_fluid = prefac_inv_s_rho_fluid * chi
!!-!!!!$        call dump_field_1d_cp(f_over_s_fluid(usz_fluid(:,:,:,1) &
!!-!!!!$                                           - two_rk * up_fluid), &  !E22
!!-!!!!$                           '/strain_dpup_flu', appisnap, nel_fluid)
!!-!!!!$        call dump_field_1d_cp(f_over_s_fluid(usz_fluid(:,:,:,1)&
!!-!!!!$                                           - two_rk * up_fluid) &  !Ekk
!!-!!!!$                            + grad_flu(:,:,:,2), &
!!-!!!!$                           '/straintrace_flu', appisnap, nel_fluid)
!!-!!!!$
!!-!!!!$        ! gradient of phi component
!!-!!!!$        call axisym_gradient_fluid(up_fluid, grad_flu)   ! 1:dsup, 2:dzup
!!-!!!!$
!!-!!!!$        call dump_field_1d_cp((- grad_flu(:,:,:,1) &
!!-!!!!$                             - f_over_s_fluid(usz_fluid(:,:,:,1) - up_fluid)), &
!!-!!!!$                           '/strain_dsup_flu', appisnap, nel_fluid)   !E12
!!-!!!!$
!!-!!!!$
!!-!!!!$        call dump_field_1d_cp(- grad_flu(:,:,:,2) / two_rk &
!!-!!!!$                            - f_over_s_fluid(usz_fluid(:,:,:,2)), &
!!-!!!!$                           '/strain_dzup_flu', appisnap, nel_fluid)   !E23
!!-!!!!$     endif   !src_type
!!-!!!!$
!!-!!!!$  endif   ! have_fluid
!!-!!
!!-!!end subroutine compute_strain_cp


end module time_evol_wave

!=========================================================================================
