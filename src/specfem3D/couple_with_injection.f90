!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


! *********************************************************************************
!
! coupling with an injection boundary (FK, DSM, AxiSEM, ..)
!
! some routines here were added by Ping Tong (TP / Tong Ping) for the FK3D calculation.
!
! when using this technique, please also reference the original work:
!
! "Three-dimensional full waveform inversion of short-period teleseismic wavefields based upon the SEM-DSM hybrid method"
! Vadim Monteiller, Sebastien Chevrot, Dimitri Komatitsch, Yi Wang
! Geophysical Journal International, Volume 202, Issue 2, 1 August 2015, Pages 811-827,
! https://doi.org/10.1093/gji/ggv189
!
!
! "A 3-D spectral-element and frequency-wave number hybrid method for high-resolution seismic array imaging"
! Tong, P; Komatitsch, D; Tseng, TL; Hung, SH; Chen, CW; Basini, P; Liu, QY
! GEOPHYSICAL RESEARCH LETTERS, Volume: 41  Issue: 20,  Pages: 7025-7034, OCT 28 2014, DOI: 10.1002/2014GL061644
! http://onlinelibrary.wiley.com/doi/10.1002/2014GL061644/abstract
!
!
! "High-resolution seismic array imaging based on an SEM-FK hybrid method"
! Ping Tong, Chin-wu Chen, Dimitri Komatitsch, Piero Basini and Qinya Liu
! Geophysical Journal International, Volume 197, Issue 1, 1 April 2014, Pages 369-395,
! https://doi.org/10.1093/gji/ggt508
!
! *********************************************************************************


  subroutine couple_with_injection_setup()

  use specfem_par
  use specfem_par_coupling

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname_trac

  ! checks if anything to do
  if (.not. IO_compute_task) return

  ! for coupling with EXTERNAL CODE
  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '**********************************************'
      write(IMAIN,*) '**** USING HYBRID METHOD                  ****'
      write(IMAIN,*) '**********************************************'
      write(IMAIN,*)
      write(IMAIN,*) 'using coupling with injection technique'
      select case(INJECTION_TECHNIQUE_TYPE)
      case (INJECTION_TECHNIQUE_IS_DSM)
        write(IMAIN,*) 'type of injection technique is DSM'
      case (INJECTION_TECHNIQUE_IS_AXISEM)
        write(IMAIN,*) 'type of injection technique is AXISEM'
      case (INJECTION_TECHNIQUE_IS_FK)
        write(IMAIN,*) 'type of injection technique is FK'
      case (INJECTION_TECHNIQUE_IS_SPECFEM)
        write(IMAIN,*) 'type of injection technique is SPECFEM'
      case default
        stop 'Invalid INJECTION_TECHNIQUE_TYPE chosen, must be 1 == DSM, 2 == AXISEM, 3 == FK or 4 == SPECFEM'
      end select
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! creates name for each rank process
    ! format: {TRACTION_PATH} / proc***_
    call create_name_database(prname_trac,myrank,TRACTION_PATH)
  endif

  ! allocates arrays for coupling
  ! note: num_abs_boundary_faces needs to be set
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
    ! coupling
    select case (INJECTION_TECHNIQUE_TYPE)
    case (INJECTION_TECHNIQUE_IS_DSM)
      ! DSM coupling
      allocate(Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2190')
      Veloc_dsm_boundary(:,:,:,:) = 0._CUSTOM_REAL

      allocate(Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2191')
      Tract_dsm_boundary(:,:,:,:) = 0._CUSTOM_REAL

      if (old_DSM_coupling_from_Vadim) then
        open(unit=IIN_veloc_dsm,file=trim(prname_trac) // 'vel.bin',status='old', &
             action='read',form='unformatted',iostat=ier)
        open(unit=IIN_tract_dsm,file=trim(prname_trac) // 'tract.bin',status='old', &
             action='read',form='unformatted',iostat=ier)
      else
        !! To verify for NOBU version (normally, remains empty)
      endif

    case (INJECTION_TECHNIQUE_IS_AXISEM)
      ! AxiSEM coupling
      allocate(Veloc_axisem(3,NGLLSQUARE*num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2192')
      Veloc_axisem(:,:) = 0._CUSTOM_REAL

      allocate(Tract_axisem(3,NGLLSQUARE*num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2193')
      Tract_axisem(:,:) = 0._CUSTOM_REAL

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  tractions: opening files ', trim(prname_trac) // 'sol_axisem'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
      ! debug
      !write(*,*) 'OPENING ', trim(prname_trac) // 'sol_axisem'

      open(unit=IIN_veloc_dsm,file=trim(prname_trac) // 'sol_axisem',status='old', &
           action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error: could not open file ',trim(prname_trac) // 'sol_axisem'
        print *,'Please check if traction file exists for coupling with AxiSEM...'
        stop 'Error opening tractions file proc****_sol_axisem'
      endif

      !! CD CD added this
      if (RECIPROCITY_AND_KH_INTEGRAL) then
        allocate(Displ_axisem_time(3,NGLLSQUARE*num_abs_boundary_faces,NSTEP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2194')
        allocate(Tract_axisem_time(3,NGLLSQUARE*num_abs_boundary_faces,NSTEP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2195')
        allocate(Tract_specfem_time(3,NGLLSQUARE*num_abs_boundary_faces,NSTEP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2196')
        allocate(Displ_specfem_time(3,NGLLSQUARE*num_abs_boundary_faces,NSTEP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2197')

        if (.not. SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
          !! We only read Specfem Tract and Displ, and Axisem Displ (Axisem Tract is read in compute_stacey_visco...)
          !! This is only for KH integral
          !! The unit numbers are here temporary
          ! user output
          if (myrank == 0) then
            write(IMAIN,*) '  KH integral: opening files ', trim(prname_trac) // 'axisem_displ_for_int_KH'
            write(IMAIN,*)
            call flush_IMAIN()
          endif
          ! debug
          !write(*,*) 'OPENING ', trim(prname_trac) // 'axisem_displ_for_int_KH, and the specfem disp and tract'

          open(unit=IIN_displ_axisem,file=trim(prname_trac) // 'axisem_displ_for_int_KH', &
            status='old',action='read',form='unformatted',iostat=ier)

          open(unit=237,file=trim(prname_trac) // 'specfem_displ_for_int_KH', &
            status='old',action='read',form='unformatted',iostat=ier)

          open(unit=238,file=trim(prname_trac) // 'specfem_tract_for_int_KH', &
            status='old',action='read',form='unformatted',iostat=ier)
        endif
      endif

    case (INJECTION_TECHNIQUE_IS_SPECFEM)
      ! SPECFEM coupling
      ! will do main setup in prepare stage when calling routine couple_with_injection_prepare_specfem_files()
    end select

  else
    ! no coupling
    ! dummy arrays (needed as function call arguments)
    allocate(Veloc_dsm_boundary(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2198')
    allocate(Tract_dsm_boundary(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2199')
    allocate(Veloc_axisem(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2200')
    allocate(Tract_axisem(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2201')
    allocate(Veloc_specfem(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2202')
    allocate(Tract_specfem(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2203')
  endif

  !! CD CD add this :
  !! We perform a first run of Specfem to save displacement and tractions of Specfem for the computation of KH integral
  !! The displ, tract, and veloc of Axisem have also to be stored
  if (SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  KH integral: opening files ', trim(prname_trac) // 'specfem_displ_for_int_KH'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! debug
    !write(*,*) 'OPENING ', trim(prname_trac) // 'specfem_displ_for_int_KH, and the specfem tract to SAVE IT'

    open(unit=237,file=trim(prname_trac) // 'specfem_displ_for_int_KH',form='unformatted')
    open(unit=238,file=trim(prname_trac) // 'specfem_tract_for_int_KH',form='unformatted')
  endif

  end subroutine couple_with_injection_setup

!
!-------------------------------------------------------------------------------------------------
!

  subroutine couple_with_injection_prepare_boundary()

  use specfem_par
  use specfem_par_coupling

  implicit none

  ! local parameters
  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  integer :: ier
  integer,dimension(:), allocatable :: npt_local_per_proc

  !! for FK point for intialization injected wavefield
  real(kind=CUSTOM_REAL) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
  real(kind=CUSTOM_REAL) :: ray_p,Tg,DF_FK

  real(kind=CUSTOM_REAL), parameter :: TOL_ZERO_TAKEOFF = 1.e-14

  ! checks if anything to do
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

  ! for forward simulation only
  if (SIMULATION_TYPE /= 1) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing injection boundary"
    call flush_IMAIN()
  endif

  select case (INJECTION_TECHNIQUE_TYPE)
  case (INJECTION_TECHNIQUE_IS_FK)
    ! FK boundary
    ! initial setup for future FK3D calculations
    ! get MPI starting time for FK
    tstart = wtime()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  using FK injection technique"
      call flush_IMAIN()
    endif

    call FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)
    call ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

    ! send FK parameters to others MPI slices
    call bcast_all_singlei(type_kpsv_fk)
    call bcast_all_singlei(nlayer)

    if (myrank > 0) then
      allocate(alpha_FK(nlayer), &
               beta_FK(nlayer), &
               rho_FK(nlayer), &
               mu_FK(nlayer), &
               h_FK(nlayer),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating arrays 2206')
      alpha_FK(:) = 0._CUSTOM_REAL; beta_FK(:) = 0._CUSTOM_REAL; rho_FK(:) = 0._CUSTOM_REAL
      mu_FK(:) = 0._CUSTOM_REAL; h_FK(:) = 0._CUSTOM_REAL
    endif

    call bcast_all_cr(alpha_FK, nlayer)
    call bcast_all_cr(beta_FK, nlayer)
    call bcast_all_cr(rho_FK, nlayer)
    call bcast_all_cr(mu_FK, nlayer)
    call bcast_all_cr(h_FK, nlayer)

    call bcast_all_singlecr(phi_FK)
    call bcast_all_singlecr(theta_FK)

    call bcast_all_singlecr(ff0)
    call bcast_all_singlecr(freq_sampling_fk)
    call bcast_all_singlecr(amplitude_fk)

    call bcast_all_singlecr(xx0)
    call bcast_all_singlecr(yy0)
    call bcast_all_singlecr(zz0)
    call bcast_all_singlecr(Z_REF_for_FK)

    call bcast_all_singlecr(tt0)
    call bcast_all_singlecr(tmax_fk)

    ! converts origin point Z to reference framework depth for FK,
    ! where top of lower half-space has to be at z==0
    zz0 = zz0 - Z_REF_for_FK

    ! converts to rad
    phi_FK   = phi_FK * PI/180.d0    ! azimuth
    theta_FK = theta_FK * PI/180.d0  ! take-off

    ! ray parameter p (according to Snell's law: sin(theta1)/v1 == sin(theta2)/v2)
    if (type_kpsv_fk == 1) then
      ! P-wave
      ray_p = sin(theta_FK)/alpha_FK(nlayer)    ! for vp (i.e., alpha)
    else if (type_kpsv_fk == 2) then
      ! SV-wave
      ray_p = sin(theta_FK)/beta_FK(nlayer)     ! for vs (i.e., beta)
    endif

    ! note: vertical incident (theta==0 -> p==0) is not handled.
    !       here, it limits ray parameter p to a very small value to handle the calculations
    if (abs(ray_p) < TOL_ZERO_TAKEOFF) ray_p = sign(TOL_ZERO_TAKEOFF,ray_p)

    ! maximum period
    Tg  = 1.d0 / ff0

    ! counts total number of (local) GLL points on absorbing boundary
    call count_num_boundary_points(num_abs_boundary_faces,abs_boundary_ispec,npt)

    ! get coupling points from each process
    allocate(npt_local_per_proc(0:NPROC-1),stat=ier)
    if (ier /= 0) stop 'Error allocating npt_local_per_proc array'
    npt_local_per_proc(:) = 0

    call gather_all_all_singlei(npt,npt_local_per_proc,NPROC)

    ! assumed total number of boundary points
    npoints_total = sum(npt_local_per_proc(:))

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  total number of coupling points         : ',npoints_total
      write(IMAIN,*) '  maximum number of local coupling points : ',maxval(npt_local_per_proc)
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    deallocate(npt_local_per_proc)

    ! safety check
    if (npoints_total == 0) then
      print *,'Error: no coupling boundary points defined for this simulation!'
      print *,'       Please check if mesh has absorbing boundaries setup properly to use COUPLE_WITH_INJECTION_TECHNIQUE'
      call exit_MPI(myrank,'No coupling boundary points found')
    endif

    !! compute the bottom middle point of the domain

    ! deallocate in case of severals runs occurs in inverse_problem program
    if (allocated(ipt_table)) deallocate(ipt_table)
    if (allocated(Veloc_FK))  deallocate(Veloc_FK)
    if (allocated(Tract_FK))  deallocate(Tract_FK)

    ! allocate memory for FK solution
    if (npt > 0) then
      allocate(ipt_table(NGLLSQUARE,num_abs_boundary_faces), stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2202')
    else
      ! dummy
      allocate(ipt_table(1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2204')
    endif
    ipt_table(:,:) = 0

    call find_size_of_working_arrays(deltat, freq_sampling_fk, tmax_fk, NF_FOR_STORING, &
                                     NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, DF_FK)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  computed FK parameters:'
      write(IMAIN,*) '    frequency sampling rate        = ', freq_sampling_fk,"(Hz)"
      write(IMAIN,*) '    number of frequencies to store = ', NF_FOR_STORING
      write(IMAIN,*) '    number of frequencies for FFT  = ', NF_FOR_FFT
      write(IMAIN,*) '    power of 2 for FFT             = ', NPOW_FOR_INTERP
      write(IMAIN,*)
      write(IMAIN,*) '    simulation time step           = ', deltat,"(s)"
      write(IMAIN,*) '    total simulation length        = ', NSTEP*deltat,"(s)"
      write(IMAIN,*)
      write(IMAIN,*) '    FK time resampling rate        = ', NP_RESAMP
      write(IMAIN,*) '    new time step for F-K          = ', NP_RESAMP * deltat,"(s)"
      write(IMAIN,*) '    new time window length         = ', tmax_fk,"(s)"
      write(IMAIN,*)
      write(IMAIN,*) '    frequency step for F-K         = ', DF_FK,"(Hz)"
      write(IMAIN,*)
      write(IMAIN,*) '  total number of points on boundary = ',npt
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! safety check with number of simulation time steps
    if (NSTEP/NP_RESAMP > NF_FOR_STORING + NP_RESAMP) then
      if (myrank == 0) then
        print *,'Error: FK time window length ',tmax_fk,' and NF_for_storing ',NF_FOR_STORING
        print *,'       are too small for chosen simulation length with NSTEP = ',NSTEP
        print *
        print *,'       you could use a smaller NSTEP <= ',NF_FOR_STORING*NP_RESAMP
        print *,'       or'
        print *,'       increase FK window length larger than ',(NSTEP/NP_RESAMP - NP_RESAMP) * NP_RESAMP * deltat
        print *,'       to have a NF for storing  larger than ',(NSTEP/NP_RESAMP - NP_RESAMP)
      endif
      stop 'Invalid FK setting'
    endif

    ! safety check
    if (NP_RESAMP == 0) then
      if (myrank == 0) then
        print *,'Error: FK resampling rate ',NP_RESAMP,' is invalid for frequency sampling rate ',freq_sampling_fk
        print *,'       and the chosen simulation DT = ',deltat
        print *
        print *,'       you could use a higher frequency sampling rate>',1./(deltat)
        print *,'       (or increase the time stepping size DT if possible)'
      endif
      stop 'Invalid FK setting'
    endif

    ! limits resampling sizes
    if (NP_RESAMP > 10000) then
      if (myrank == 0) then
        print *,'Error: FK resampling rate ',NP_RESAMP,' is too high for frequency sampling rate ',freq_sampling_fk
        print *,'       and the chosen simulation DT = ',deltat
        print *
        print *,'       you could use a higher frequency sampling rate > ',1./(10000*deltat)
        print *,'       (or increase the time stepping size DT if possible)'
      endif
      stop 'Invalid FK setting'
    endif

    if (npt > 0) then
      !! arrays for storing FK solution --------------------------------------------
      allocate(Veloc_FK(NDIM, npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2210')
      if (ier /= 0) stop 'error while allocating Veloc_FK'
      Veloc_FK(:,:,:) = 0._CUSTOM_REAL

      allocate(Tract_FK(NDIM, npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2210')
      if (ier /= 0) stop 'error while allocating Veloc_FK'
      Tract_FK(:,:,:) = 0._CUSTOM_REAL

      call FK3D(type_kpsv_fk, nlayer, NSTEP, npt, &
                ray_p, phi_FK, xx0, yy0, zz0, Tg, &
                tt0, alpha_FK, beta_FK, rho_FK, h_FK, &
                NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)
    else
      ! dummy
      allocate(Veloc_FK(1,1,1),Tract_FK(1,1,1))
    endif

    call synchronize_all()

    ! get MPI ending time for FK
    tCPU = wtime() - tstart

    ! user output
    if (myrank == 0) then
        write(IMAIN,'(a35,1x, f20.2, a7)')  " Elapsed time for FK computation : ",  tCPU, " sec. "
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    deallocate(alpha_FK, beta_FK, rho_FK, mu_FK, h_FK)

    ! * end of initial setup for future FK3D calculations *

    case (INJECTION_TECHNIQUE_IS_SPECFEM)
      ! SPECFEM coupling
      call couple_with_injection_prepare_specfem_files()

  end select

  end subroutine couple_with_injection_prepare_boundary

!
!-------------------------------------------------------------------------------------------------
!

!! count the number of point in the mesh partition boundary : npt

  subroutine count_num_boundary_points(num_abs_boundary_faces, abs_boundary_ispec, npt)

  use constants, only: NGLLSQUARE,myrank

  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  implicit none

  integer,                                    intent(inout)     :: npt
  integer,                                    intent(in)        :: num_abs_boundary_faces
  ! absorbing boundary surface
  integer, dimension(num_abs_boundary_faces), intent(in)        :: abs_boundary_ispec

  ! local parameters
  integer :: ispec, iface

  ! total number of injection points
  npt = 0

  do iface = 1, num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      npt = npt + NGLLSQUARE
    endif

    if (ispec_is_acoustic(ispec)) then
      ! reference GLL points on boundary face
      npt = npt + NGLLSQUARE
    endif

    if (ispec_is_poroelastic(ispec)) then
      ! poroelastic domains not supported yet
      print *,'Error: rank ',myrank,' has injection point in poroelastic domain - not supported yet'
      stop 'Wavefield injection for poroelastic domains not supported yet'
    endif
  enddo

  end subroutine count_num_boundary_points

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FK3D(kpsv, nlayer, NSTEP, npt, &
                  ray_p, phi, xx0, yy0, zz0, Tg, &
                  tt0, alpha_FK, beta_FK, rho_FK, h_FK, &
                  NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

  use constants

  use specfem_par, only: xstore, ystore, zstore, kappastore, mustore, rhostore

  use specfem_par, only: ibool, deltat, &
                         abs_boundary_ijk, abs_boundary_normal, &
                         abs_boundary_ispec, num_abs_boundary_faces

  use specfem_par_coupling, only: xx, yy, zz, xi1, xim, bdlambdamu, &
                                  nmx, nmy, nmz, Z_REF_for_FK, &
                                  ipt_table

  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  ! for plotting/debugging
  use specfem_par, only: DT,t0
  use specfem_par_coupling, only: Veloc_FK,Tract_FK

  implicit none

  integer,intent(in)   :: kpsv,nlayer,NSTEP,npt
  integer,intent(in)   :: NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP

  ! source
  real(kind=CUSTOM_REAL),intent(in) :: ray_p,phi,xx0,yy0,zz0,Tg,tt0
  real(kind=CUSTOM_REAL),intent(in) :: DF_FK
  ! model
  real(kind=CUSTOM_REAL),dimension(nlayer),intent(in) :: alpha_FK,beta_FK,rho_FK,h_FK

  ! local parameters
  real(kind=CUSTOM_REAL) :: rho_tmp,kappa_tmp,mu_tmp,xi
  integer :: ispec,iglob,i,j,k,iface,igll,ier,ipt

  integer :: ii, kk, iim1, iip1, iip2, it_tmp
  real(kind=CUSTOM_REAL) :: cs1,cs2,cs3,cs4,w,time_t
  real(kind=CUSTOM_REAL) :: vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
  real(kind=CUSTOM_REAL) :: x_loc,y_loc,z_loc
  character(len=128) :: filename1,filename2

  ! checks if anything to do
  if (npt == 0) return

  ! allocates temporary arrays
  allocate(xx(npt),yy(npt),zz(npt),xi1(npt),xim(npt),bdlambdamu(npt),nmx(npt),nmy(npt),nmz(npt),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2216')
  xx(:) = 0.0; yy(:) = 0.0; zz(:) = 0.0
  xi1(:) = 0.0; xim(:) = 0.0; bdlambdamu(:) = 0.0
  nmx(:) = 0.0; nmy(:) = 0.0; nmz(:) = 0.0

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  ! loops over boundary points
  ipt = 0
  do iface = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)
    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        iglob = ibool(i,j,k,ispec)

        ipt = ipt + 1
        ipt_table(igll,iface) = ipt

        xx(ipt) = xstore(iglob)
        yy(ipt) = ystore(iglob)
        zz(ipt) = zstore(iglob) - Z_REF_for_FK  !! VM VM put z in FK system of coordinate (z == 0 at top of lower half-space)

        nmx(ipt) = abs_boundary_normal(1,igll,iface)
        nmy(ipt) = abs_boundary_normal(2,igll,iface)
        nmz(ipt) = abs_boundary_normal(3,igll,iface)

        rho_tmp   = rhostore(i,j,k,ispec)
        kappa_tmp = kappastore(i,j,k,ispec)
        mu_tmp    = mustore(i,j,k,ispec)

        xi       = mu_tmp/(kappa_tmp + 4.0/3.0 * mu_tmp)
        xi1(ipt) = 1.0 - 2.0 * xi
        xim(ipt) = (1.0 - xi) * mu_tmp
        bdlambdamu(ipt) = (3.0 * kappa_tmp - 2.0 * mu_tmp) / (6.0 * kappa_tmp + 2.0 * mu_tmp)  ! Poisson's ratio 3K-2G/[2(3K+G)]
      enddo
    endif ! elastic

    if (ispec_is_acoustic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        iglob = ibool(i,j,k,ispec)

        ipt = ipt + 1
        ipt_table(igll,iface) = ipt

        xx(ipt) = xstore(iglob)
        yy(ipt) = ystore(iglob)
        zz(ipt) = zstore(iglob) - Z_REF_for_FK  !! VM VM put z in FK system of coordinate (z == 0 at top of lower half-space)

        nmx(ipt) = abs_boundary_normal(1,igll,iface)
        nmy(ipt) = abs_boundary_normal(2,igll,iface)
        nmz(ipt) = abs_boundary_normal(3,igll,iface)

        rho_tmp   = rhostore(i,j,k,ispec)
        kappa_tmp = kappastore(i,j,k,ispec)
        mu_tmp    = 0.0

        xi       = 0.0              ! xi = mu_tmp/(kappa_tmp + 4.0/3.0 * mu_tmp)
        xi1(ipt) = 1.0              ! xi1 = 1.0 - 2.0 * xi
        xim(ipt) = 0.0              ! xim = (1.0 - xi) * mu_tmp
        bdlambdamu(ipt) = 0.5       ! bdlambdamu = (3.0 * kappa_tmp - 2.0 * mu_tmp) / (6.0 * kappa_tmp + 2.0 * mu_tmp)
      enddo
    endif ! acoustic

    if (ispec_is_poroelastic(ispec)) then
      ! poroelastic domains not supported yet
      print *,'Error: rank ',myrank,' has injection point in poroelastic domain - not supported yet'
      stop 'Wavefield injection for poroelastic domains not supported yet'
    endif
  enddo

  ! saftey check
  if (ipt /= npt) then
    print *,'Error: rank ',myrank,' has invalid number of injection points ',ipt,' should be ',npt
    stop 'Error invalid number of injection points'
  endif

  ! FK wavefield
  call FK(alpha_FK, beta_FK, rho_FK, h_FK, nlayer, &
          Tg, ray_p, phi, xx0, yy0, zz0, &
          tt0, deltat, NSTEP, npt, &
          kpsv, NF_FOR_STORING, NPOW_FOR_FFT,  NP_RESAMP, DF_FK)

  ! file output for plotting
  if (myrank == 0) then
    ! output
    write(IMAIN,*) "  creating sample files:"
    do i = 1,2
      ! filename
      write(filename1,'(a,i1,a)') "plot_FK_Veloc.",i,".dat"
      write(filename2,'(a,i1,a)') "plot_FK_Tract.",i,".dat"
      ! user output
      write(IMAIN,*) "    ",trim(filename1)," and ",trim(filename2)

      ! files
      open(unit=8888,file="OUTPUT_FILES/"//trim(filename1),status='unknown',iostat=ier)
      if (ier /= 0) stop 'Error opening FK Veloc plot file'
      open(unit=8889,file="OUTPUT_FILES/"//trim(filename2),status='unknown',iostat=ier)
      if (ier /= 0) stop 'Error opening FK Tract plot file'

      ! boundary point index
      ! first and last boundary point
      if (i == 1) then
        ipt = ipt_table(1,1)
      else
        ipt = ipt_table(1,num_abs_boundary_faces)
      endif

      ! point locations
      x_loc = xx(ipt)
      y_loc = yy(ipt)
      z_loc = zz(ipt) + Z_REF_for_FK  ! original location

      ! header
      write(8888,*) "# FK Velocity - data point"
      write(8888,*) "# point id      : ",ipt
      write(8888,*) "# point location: x/y/z = ",x_loc,y_loc,z_loc
      write(8888,*) "# point param   : xi1/xim/bdlambdamu = ",xi1(ipt),xim(ipt),bdlambdamu(ipt)
      write(8888,*) "# line format   : #time #Vx #Vy #Vz"

      write(8889,*) "# FK Traction - data point"
      write(8889,*) "# point id      : ",ipt
      write(8889,*) "# point location: x/y/z = ",x_loc,y_loc,z_loc
      write(8889,*) "# point param   : xi1/xim/bdlambdamu = ",xi1(ipt),xim(ipt),bdlambdamu(ipt)
      write(8889,*) "# line format   : #time #Tx #Ty #Tz"

      ! data
      do it_tmp = 1,NSTEP
        ! FK coupling
        !! find indices
        ! example:
        !   np_resamp = 1 and it = 1,2,3,4,5,6, ..
        !   --> ii = 1,2,3,4,5,6,..
        !   np_resamp = 2 and it = 1,2,3,4,5,6, ..
        !   --> ii = 1,1,2,2,3,3,..
        ii = floor( real(it_tmp + NP_RESAMP - 1) / real( NP_RESAMP))
        ! example:
        !       kk = 1,2,1,2,1,2,,..
        kk = it_tmp - (ii-1) * NP_RESAMP
        ! example:
        !       w = 0,1/2,0,1/2,..
        w = dble(kk-1) / dble(NP_RESAMP)

        ! Cubic spline values
        cs4 = w*w*w/6.d0
        cs1 = 1.d0/6.d0 + w*(w-1.d0)/2.d0 - cs4
        cs3 = w + cs1 - 2.d0*cs4
        cs2 = 1.d0 - cs1 - cs3 - cs4

        ! interpolation indices
        iim1 = ii-1        ! 0,..
        iip1 = ii+1        ! 2,..
        iip2 = ii+2        ! 3,..

        ! interpolates velocity/stress
        ! velocity
        vx_FK = cs1 * Veloc_FK(1,ipt,iim1) + cs2 * Veloc_FK(1,ipt,ii) + cs3 * Veloc_FK(1,ipt,iip1) + cs4 * Veloc_FK(1,ipt,iip2)
        vy_FK = cs1 * Veloc_FK(2,ipt,iim1) + cs2 * Veloc_FK(2,ipt,ii) + cs3 * Veloc_FK(2,ipt,iip1) + cs4 * Veloc_FK(2,ipt,iip2)
        vz_FK = cs1 * Veloc_FK(3,ipt,iim1) + cs2 * Veloc_FK(3,ipt,ii) + cs3 * Veloc_FK(3,ipt,iip1) + cs4 * Veloc_FK(3,ipt,iip2)

        ! stress
        tx_FK = cs1 * Tract_FK(1,ipt,iim1) + cs2 * Tract_FK(1,ipt,ii) + cs3 * Tract_FK(1,ipt,iip1) + cs4 * Tract_FK(1,ipt,iip2)
        ty_FK = cs1 * Tract_FK(2,ipt,iim1) + cs2 * Tract_FK(2,ipt,ii) + cs3 * Tract_FK(2,ipt,iip1) + cs4 * Tract_FK(2,ipt,iip2)
        tz_FK = cs1 * Tract_FK(3,ipt,iim1) + cs2 * Tract_FK(3,ipt,ii) + cs3 * Tract_FK(3,ipt,iip1) + cs4 * Tract_FK(3,ipt,iip2)

        ! time
        time_t = (it_tmp-1) * DT - t0

        ! file output
        write(8888,*) time_t, vx_FK, vy_FK, vz_FK
        write(8889,*) time_t, tx_FK, ty_FK, tz_FK
      enddo
      ! closes files
      close(8888)
      close(8889)
    enddo
    write(IMAIN,*)
  endif

  ! free temporary arrays
  deallocate(xx, yy, zz, xi1, xim, bdlambdamu, nmx, nmy, nmz)

  end subroutine FK3D

!
!-------------------------------------------------------------------------------------------------
!

  ! FK elastic + acoustic
  subroutine FK(vp, vs, rho, H, nlayer, &
                Tg, ray_p, phi, x0, y0, z0, &
                t0, dt, npts, npt, &
                kpsv, NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

  use constants, only: myrank,CUSTOM_REAL,IMAIN,PI,TINYVAL,NGLLSQUARE

  use specfem_par_coupling, only: Veloc_FK, Tract_FK, &
                                  xx, yy, zz, xi1, xim, bdlambdamu, &
                                  nmx, nmy, nmz, NPTS_STORED, NPTS_INTERP, &
                                  amplitude_fk, ipt_table,Z_REF_for_FK

  use specfem_par, only: num_abs_boundary_faces,abs_boundary_ispec
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par, only: ACOUSTIC_SIMULATION

  implicit none

  integer,                parameter                         :: CUSTOM_CMPLX = 8
  real(kind=CUSTOM_REAL), parameter                         :: zign_neg = -1.0

  ! input and output
  integer,                                     intent(in)   :: nlayer, npts,npt, kpsv
  integer                                                   :: NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP
  ! model
  real(kind=CUSTOM_REAL),  dimension(nlayer),  intent(in)   :: vp(nlayer),vs(nlayer),rho(nlayer),H(nlayer)

  ! source
  real(kind=CUSTOM_REAL),                      intent(in)    :: dt, ray_p, phi, x0, y0, z0, Tg, t0, DF_FK


  ! local parameters
  real(kind=CUSTOM_REAL),     dimension(:),   allocatable    :: fvec, dtmp
  real(kind=CUSTOM_REAL),     dimension(:,:), allocatable    :: field
  real(kind=CUSTOM_REAL),     dimension(:),   allocatable    :: tmp_t1, tmp_t2, tmp_t3, tmp_it1

  complex(kind=CUSTOM_CMPLX), dimension(:,:), allocatable    :: coeff, field_f
  complex(kind=CUSTOM_CMPLX), dimension(:),   allocatable    :: tmp_f1, tmp_f2, tmp_f3
  complex(kind=CUSTOM_CMPLX)                                 :: C_3,stf_coeff,a,b,c,d,delta_mat
  complex(kind=CUSTOM_CMPLX)                                 :: dx_f,dz_f,txz_f,tzz_f
  complex(kind=CUSTOM_CMPLX)                                 :: N_mat(4,4),N1_mat(2,2)

  real(kind=CUSTOM_REAL)                                     :: sigma_rr,sigma_rt,sigma_rz,sigma_tt,sigma_tz,sigma_zz
  real(kind=CUSTOM_REAL)                                     :: Txx_tmp, Txy_tmp, Txz_tmp, Tyy_tmp, Tyz_tmp, Tzz_tmp
  real(kind=CUSTOM_REAL)                                     :: dt_fk,df,om,Tdelay,f_Nyquist,C_1

  integer                                                    :: npow,npts2,nf,nf2,nn,ii,ipt,i,j,lpts
  integer                                                    :: npoints2
  integer                                                    :: ier,iface,igll,ispec,ilayer

  ! spline work array
  double precision, dimension(:), allocatable                :: tmp_c

  ! matrices for computation
  complex(kind=CUSTOM_CMPLX)                                :: eta_alpha(nlayer),eta_beta(nlayer),gamma0(nlayer),gamma1(nlayer)
  complex(kind=CUSTOM_CMPLX)                                :: E_mat(4,4),Qmat(2,2),Pmat(4,4),Qmat_I(2,2)
  logical                                                   :: have_fluid_layer
  real(kind=CUSTOM_REAL)                                    :: height,two_mul
  complex(kind=CUSTOM_CMPLX)                                :: eta_p,eta_s
  complex(kind=CUSTOM_CMPLX)                                :: MM,bot_vec(4),G_mat(4,4)
  real(kind=CUSTOM_REAL),parameter                          :: THRESHOLD_VS = 1.0e-6
  integer                                                   :: ilayer_ac

!! DK DK here is the hardwired maximum size of the array
!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
  real(kind=CUSTOM_REAL) :: mpow(30)

  ! taper
  ! idea: tapers the onset of the injection traces to diminish numerical noise.
  !       the inverse FFT can lead to non-zero onsets for coarse frequency sampling.
  logical, parameter :: USE_TAPERED_BEGINNING = .true.
  integer, parameter :: taper_nlength = 20    ! tapers first 20 steps
  real(kind=CUSTOM_REAL) :: taper

  ! fixed parameters
  integer, parameter     :: nvar = 5
  logical, parameter     :: comp_stress = .true.

  ! initializations
  !! new way to do time domain resampling
  df    = DF_FK
  nf2   = NF_FOR_STORING+1   ! number of positive frequency sample points
  nf    = 2*NF_FOR_STORING   ! number of total frequencies after symmetrisation
  npts2 = nf                 ! number of samples in time serie

  !! VM VM recompute new values for new way to do
  npow = ceiling(log(npts2*1.0)/log(2.0))
  npts2 = 2**npow
  NPOW_FOR_FFT = npow

  !! number of points for resampled vector
  npoints2 = NP_RESAMP*(npts2-1)+1

  dt_fk = 1.0_CUSTOM_REAL/(df*(npts2-1))
  f_Nyquist = 1.0_CUSTOM_REAL/(2.0_CUSTOM_REAL * dt)    ! Nyquist frequency of specfem time serie

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  Entering the FK synthetics program:'
    write(IMAIN,*) '    Number of samples stored for FK solution = ', NF_FOR_STORING
    write(IMAIN,*) '    Number of points used for FFT            = ', npts2
    write(IMAIN,*) '    Total time length used for FK            = ', t0+(npts2-1)*dt_fk,'(s)'
    write(IMAIN,*)
    write(IMAIN,*) '    simulation Nyquist frequency             = ', f_Nyquist,'(Hz)'
    write(IMAIN,*) '    FK time step                             = ', dt_fk
    write(IMAIN,*) '    FK frequency step                        = ', df
    write(IMAIN,*) '    power of 2 for FFT                       = ', npow
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  !! check if dt_fk is compatible with dt_specfem
  allocate(fvec(nf2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2218')
  if (ier /= 0) stop 'error while allocating'
  fvec(:) = 0.0_CUSTOM_REAL
  do ii = 1, nf2
    fvec(ii) = (ii-1)*df
  enddo

  allocate(coeff(2,nf2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2219')
  if (ier /= 0) stop 'error while allocating'
  coeff(:,:) = (0.d0,0.d0)

  allocate(field_f(nf,nvar),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2220')
  if (ier /= 0) stop 'error while allocating'
  field_f(:,:) = (0.d0,0.d0)

  allocate(field(npts2,nvar),dtmp(npts),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2221')
  if (ier /= 0) stop 'error while allocating'
  field(:,:) = 0._CUSTOM_REAL; dtmp(:) = 0._CUSTOM_REAL

  !! allocate debug vectors
  allocate(tmp_f1(npts2), tmp_f2(npts2), tmp_f3(npts2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2222')
  if (ier /= 0) stop 'error while allocating'
  tmp_f1(:) = (0.d0,0.d0)
  tmp_f2(:) = (0.d0,0.d0)
  tmp_f3(:) = (0.d0,0.d0)

  if (ier /= 0) stop 'error while allocating'
  allocate(tmp_t1(npts2), tmp_t2(npts2), tmp_t3(npts2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2223')
  if (ier /= 0) stop 'error while allocating'
  tmp_t1(:) = 0.0_CUSTOM_REAL
  tmp_t2(:) = 0.0_CUSTOM_REAL
  tmp_t3(:) = 0.0_CUSTOM_REAL

  allocate(tmp_it1(npoints2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2224')
  if (ier /= 0) stop 'error while allocating'
  tmp_it1(:) = 0.0_CUSTOM_REAL

  ! temporary work array for splines
  allocate(tmp_c(npts2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2225')
  tmp_c(:) = 0.d0

  NPTS_STORED = npts2
  NPTS_INTERP = npoints2

  nn = int(-t0/dt) ! what if this is not an integer number?

  !! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
  !! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
  do i = 1,npow
    mpow(i) = 2**(npow-i)
  enddo


  ! check if top layers are in fluid material
  ! note: indexing assumes that last layer (nlayer) being the bottom, lower halfspace,
  !       and the first layer (1) being at the top surface
  have_fluid_layer = .false.
  do j = nlayer,1,-1
    if (vs(j) < THRESHOLD_VS) then
      ilayer_ac = j      ! bottom layer of acoustic layer(s)
      have_fluid_layer = .true.
      exit
    endif
  enddo

  ! check no elastic layers for ilayer < ilayer_ac
  if (have_fluid_layer .and. (ilayer_ac > 1)) then
    if (sum(vs(1:ilayer_ac)) / ilayer_ac > THRESHOLD_VS) then
      if (myrank == 0) print *,'Also check the FK model, make sure fluid layers are on the top of elastic layers'
      stop
    endif
  endif

  ! make sure ACOUSTIC_SIMULATION is enabled
  if (have_fluid_layer .and. (.not. ACOUSTIC_SIMULATION)) then
    if (myrank == 0) then
      print *,'FK model contains fluid layers, but ACOUSTIC_SIMULATION is not enabled'
    endif
    stop
  endif

  if (have_fluid_layer) then
    if (myrank == 0) write(IMAIN,*) 'FK simulation: acoustic + elastic'
  else
    if (myrank == 0) write(IMAIN,*) 'FK simulation: elastic'
  endif

  ! compute temporary variables
  do i = 1,nlayer
    ! P
    ! see (A5): E_23 = -i nu_p / k = -i omega sqrt(1/alpha^2 - p^2) / k
    !           factor eta_alpha = -i sqrt(1/alpha^2 - p^2)
    eta_alpha(i) = -cmplx(0,1) * sqrt( 1.0 / vp(i)**2 - ray_p**2 )

    ! SV
    ! see (A5): E_11 = -i nu_s / k = -i omega sqrt(1/beta^2 - p^2) / k
    !           factor eta_beta = -i sqrt(1/beta^2 - p^2)
    if (vs(i) < THRESHOLD_VS) then
      eta_beta(i) = 0.
    else
      eta_beta(i) = -cmplx(0,1) * sqrt( 1.0 / vs(i)**2 - ray_p**2 )
    endif

    ! auxiliary variables
    gamma0(i) = 2.0 * vs(i)**2 * ray_p**2
    if (vs(i) < THRESHOLD_VS) then
      gamma1(i) = 0.
    else
      gamma1(i) = 1.0 - 1.0/gamma0(i)
    endif
  enddo

  if (myrank == 0) then
    write(IMAIN,*) '    starting from ',nn,' points before time 0'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! amplitude in half space
  if (kpsv == 1) then
    ! P-SV
    C_3 = amplitude_fk * cmplx(0,1.) * ray_p * vp(nlayer)      ! amp. of incoming P in the bot. layer
    eta_p = sqrt(1.0/vp(nlayer)**2 - ray_p**2)                 ! vertical slowness for lower layer
    if (myrank == 0) write(IMAIN,*) '  Incoming P : C_3,  ray_p, eta = ', C_3, ray_p, eta_p
  else
    ! SV-wave
    ! for C_2 = sin(inc) (u=[cos(inc), sin(inc)])
    C_1 = amplitude_fk * ray_p * vs(nlayer)                   ! amp. of incoming S in the bot. layer
    eta_s = sqrt(1.0/vs(nlayer)**2 - ray_p**2)                ! vertical slowness for lower layer

    if (myrank == 0 ) write(IMAIN,*) '  Incoming S :  C_1,  ray_p, eta = ', C_1, ray_p, eta_s
  endif

  ! pre-computed factor for half-space (layer with index nlayer)
  two_mul = 2.0 * rho(nlayer) * vs(nlayer) * vs(nlayer)

  ! note: vertical incident (p=0) is not handled in Tong et al. explanations.
  if (abs(ray_p) < 1.e-15) stop 'Invalid ray parameter p (cannot be zero) in FK() routine'

  ! E matrix
  ! initializes matrix
  E_mat(:,:) = (1.0,0.0)

  ! Tong et al. (2014), appendix (A10) E_0:
  ! note: E_mat is not omega dependent
  E_mat(1,1) =  eta_beta(nlayer) / ray_p
  E_mat(1,2) = -E_mat(1,1)
  E_mat(2,3) =  eta_alpha(nlayer) / ray_p
  E_mat(2,4) = -E_mat(2,3)
  E_mat(3,1) = two_mul * gamma1(nlayer)
  E_mat(3,2) = E_mat(3,1)
  E_mat(3,3) = two_mul * eta_alpha(nlayer) / ray_p
  E_mat(3,4) = -E_mat(3,3)
  E_mat(4,1) = two_mul * eta_beta(nlayer) / ray_p
  E_mat(4,2) = -E_mat(4,1)
  E_mat(4,3) = E_mat(3,1)
  E_mat(4,4) = E_mat(3,1)

  ! now loop every frequency to determine coefs in half space
  do ii = 1,nf2
    om = 2.0 * PI * fvec(ii)

    ! apply propagation matrix in elastic layers
    N_mat(:,:) = E_mat(:,:)

    ! note: indexing assumes that last layer (nlayer) being the bottom, lower halfspace,
    !       and the first layer (1) being at the top surface
    ilayer = 1
    if (have_fluid_layer) ilayer = ilayer_ac + 1

    do i = nlayer-1,ilayer,-1
      call fk_propagator_psv(om,eta_alpha(i),eta_beta(i),rho(i), &
                             vs(i),H(i),ray_p,gamma1(i),Pmat)
      ! resulting matrix
      N_mat = matmul(Pmat,N_mat) * gamma0(i)
    enddo

    ! apply propagation matrix in acoustic layers
    if (have_fluid_layer) then
      Qmat_I(:,:) = 0.0_CUSTOM_CMPLX
      Qmat_I(1,1) = cmplx(1.0,0.0,kind=CUSTOM_CMPLX)
      Qmat_I(2,2) = cmplx(1.0,0.0,kind=CUSTOM_CMPLX)
      ! note: assumes that ilayer_ac is the bottom layer of acoustic layers
      do j = ilayer_ac,1,-1
        call fk_propagator_ac(om,eta_alpha(j),rho(j),H(j),ray_p,Qmat)
        ! resulting matrix
        Qmat_I = matmul(Qmat,Qmat_I)
      enddo
      Qmat = Qmat_I
    endif

    ! determine coefs in half space
    if (.not. have_fluid_layer) then
      ! inverse matrix
      a = N_mat(3,2); b = N_mat(3,4); c = N_mat(4,2); d = N_mat(4,4)
      delta_mat = a*d - b*c
      if (abs(delta_mat) > TINYVAL) then
        if (kpsv == 1) then
          coeff(1,ii) = -(d*N_mat(3,3) - b*N_mat(4,3)) / delta_mat * C_3
          coeff(2,ii) = -(-c*N_mat(3,3) + a*N_mat(4,3)) / delta_mat * C_3
        else
          coeff(1,ii) = -(d*N_mat(3,1) - b*N_mat(4,1)) / delta_mat * C_1
          coeff(2,ii) = -(-c*N_mat(3,1) + a*N_mat(4,1)) / delta_mat * C_1
        endif
      else
        coeff(1,ii) = (0.d0,0.d0)
        coeff(2,ii) = (0.d0,0.d0)
      endif
    else
      N1_mat(1,1) = N_mat(3,2)
      N1_mat(1,2) = N_mat(3,4)
      if (kpsv == 1) then
        N1_mat(2,1) = Qmat(2,1) * N_mat(2,2) - Qmat(2,2) * N_mat(4,2)
        N1_mat(2,2) = Qmat(2,1) * N_mat(2,4) - Qmat(2,2) * N_mat(4,4)
        MM = Qmat(2,1) * N_mat(2,3) - Qmat(2,2) * N_mat(4,3)
      else
        N1_mat(2,1) = Qmat(2,1) * N_mat(2,2) - Qmat(2,2) * N_mat(4,2)
        N1_mat(2,2) = Qmat(2,1) * N_mat(2,4) - Qmat(2,2) * N_mat(4,4)
        MM = Qmat(2,1) * N_mat(2,1) - Qmat(2,2) * N_mat(4,1)
      endif
      a = N1_mat(1,1); b = N1_mat(1,2); c = N1_mat(2,1); d = N1_mat(2,2)
      delta_mat = a*d - b*c
      if (abs(delta_mat) > TINYVAL) then
        if (kpsv == 1) then
          coeff(1,ii) = -( d * N_mat(3,3) - b * MM) / delta_mat * C_3
          coeff(2,ii) = -(-c * N_mat(3,3) + a * MM) / delta_mat * C_3
        else
          coeff(1,ii) = -( d * N_mat(3,1) - b * MM) / delta_mat * C_1
          coeff(2,ii) = -(-c * N_mat(3,1) + a * MM) / delta_mat * C_1
        endif
      else
        coeff(1,ii) = (0.d0,0.d0)
        coeff(2,ii) = (0.d0,0.d0)
      endif
    endif

      ! if (ii == 10 .and. myrank == 0) then
      !   print *,'debug: Rayleigh coeff ',coeff(1,ii),coeff(2,ii),delta_mat
      !   print *,'N_mat'
      !   print *,N_mat
      ! endif
  enddo

  ! loop every point to calculate stress/velocity
  do iface = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)

    ! GLL points on boundary face
    do igll = 1,NGLLSQUARE
      ! point index using table lookup
      ipt = ipt_table(igll,iface)

      ! initializes
      field_f(:,:) = (0.d0,0.d0)

      ! time delay with respect to top of lower half-space (set to be at z==0)
      Tdelay = ray_p * (xx(ipt)-x0) * cos(phi) + ray_p * (yy(ipt)-y0) * sin(phi) + eta_p * (0.0-z0)

      do ii = 1, nf2
        om = 2.0 * PI * fvec(ii)                                 !! pulsation

        stf_coeff = exp(-(om * Tg/2)**2)                         !! apodization window
        stf_coeff = stf_coeff * exp(cmplx(0,-1)*om*Tdelay)

        ! bottom vector
        bot_vec(:) = (0.,0.)
        if (kpsv == 1) then
          bot_vec(2) = coeff(1,ii)
          bot_vec(3) = C_3
          bot_vec(4) = coeff(2,ii)
        else
          bot_vec(1) = C_1
          bot_vec(3) = coeff(1,ii)
          bot_vec(4) = coeff(2,ii)
        endif

        ! find which layer this point is in
        if (ispec_is_elastic(ispec)) then
          ! elastic element
          if (zz(ipt) <= 0.0) then
            ! in lower half space
            G_mat(:,:) = (0.0,0.0)
            ! incident, up-going S-wave
            ! (A5): Gamma_11 = e^(-i nu_s z) = e^( (-i nu_s) * z ) = e^(nu_be * z)
            G_mat(1,1) = exp(om * eta_beta(nlayer) * zz(ipt))
            ! reflected, down-going S-wave
            ! (A5): Gamma_22 = e^(i nu_s z)  = e^( -(-i nu_s) * z ) = e^(-nu_be * z)
            G_mat(2,2) = exp(-om * eta_beta(nlayer) * zz(ipt))
            ! incident, up-going P-wave
            ! (A5): Gamma_33 = e^(-i nu_p z) = e^( (-i nu_p) * z ) = e^(nu_al * z)
            G_mat(3,3) = exp(om * eta_alpha(nlayer) * zz(ipt))
            ! reflected, down-going P-wave
            ! (A5): Gamma_44 = e^(i nu_p z) = e^( -(-i nu_p) * z ) = e^(-nu_al * z)
            G_mat(4,4) = exp(-om * eta_alpha(nlayer) * zz(ipt))
            N_mat = matmul(E_mat,G_mat)
          else
            ! in layers
            ! determines layer in which the point lies
            ! note: indexing assumes that last layer (nlayer) being the bottom, lower halfspace,
            !       and the first layer (1) being at the top surface
            ilayer = nlayer
            do j = nlayer-1 , 1 , -1
              if (zz(ipt) <= sum(H(j:nlayer-1))) then
                ilayer = j; exit
              endif
            enddo

            if (have_fluid_layer .and. ilayer <= ilayer_ac) then
              print *,'Error: point cannot be in acoustic domain'
              print *,'  z = ',zz(ipt),'z_ref = ',zz(ipt) + Z_REF_for_FK,'H_layers = ', sum(H(ilayer_ac+1:nlayer-1))
              print *,'  layer = ',ilayer,'ilayer_ac',ilayer_ac,'nlayer = ',nlayer
              print *,'  H = ',H(:)
              print *,'exiting...'
              stop
            endif

            ! new thickness
            height = zz(ipt) - sum(H(ilayer+1:nlayer-1))

            ! propagation to this point
            N_mat(:,:) = E_mat(:,:)
            do j = nlayer-1,ilayer,-1
              if ( j > ilayer) then
                call fk_propagator_psv(om,eta_alpha(j),eta_beta(j),rho(j), &
                                       vs(j),H(j),ray_p,gamma1(j),Pmat)
              else
                call fk_propagator_psv(om,eta_alpha(j),eta_beta(j),rho(j), &
                                       vs(j),height,ray_p,gamma1(j),Pmat)
              endif
              ! resulting matrix
              N_mat = gamma0(j) * matmul(Pmat,N_mat)
            enddo
          endif ! endif (zz(ipt) <= 0.0)

          !! zz(ipt) is the height of point with respect to the lower layer
          bot_vec = matmul(N_mat,bot_vec)
          dx_f = bot_vec(1) ! y_1
          dz_f = bot_vec(2) ! y_3

          ! for the Stacey boundary contribution, we need velocity = (i om) displacement (in frequency domain)
          field_f(ii,1) = stf_coeff * dx_f * cmplx(0,-1) * cmplx(0,om)             ! (i om)u_x
          field_f(ii,2) = stf_coeff * dz_f * cmplx(0,om)                           ! (i om)u_z

          ! stress
          if (comp_stress) then
            txz_f = bot_vec(3)      ! tilde{y}_4
            tzz_f = bot_vec(4)      ! tilde{y}_6

            field_f(ii,3) = stf_coeff * om * ray_p * (xi1(ipt)*tzz_f - 4.0*xim(ipt)*dx_f) ! T_xx
            field_f(ii,4) = stf_coeff * om * ray_p * txz_f * cmplx(0,-1)                  ! T_xz
            field_f(ii,5) = stf_coeff * om * ray_p * tzz_f                                ! T_zz
          endif

        else if (ispec_is_acoustic(ispec)) then
          ! acoustic element
          if (nlayer == 1 .and. ilayer_ac == 1) then
            ! single acoustic layer (nlayer=1 and ilayer_ac=1)
            ! in this case, all points are within the single acoustic half-space.
            ilayer = 1
            height = zz(ipt)  ! Point's coordinate relative to the top of the acoustic half-space (z=0 in FK system)
          else
            ! multi-layered acoustic models
            ! in this case, points should be in fluid layers
            ! determines layer in which the point lies
            ! note: indexing assumes that last layer (ilayer_ac) being the bottom acoustic layer,
            !       and the first layer (1) being at the top surface
            ilayer = ilayer_ac + 1
            do j = ilayer_ac + 1 , 1 , -1
              if (zz(ipt) <= sum(H(j:nlayer-1))) then
                ilayer = j; exit
              endif
            enddo
            height = zz(ipt) - sum(H(ilayer+1:nlayer-1))

            ! checks position
            if (height < 0.) then
              print *,'Error: please check, point is in the air'
              print *,'  z = ',zz(ipt),'z_ref = ',zz(ipt) + Z_REF_for_FK,'H_layers = ',sum(H(ilayer+1:nlayer-1))
              print *,'  layer = ',ilayer,'ilayer_ac',ilayer_ac,'nlayer = ',nlayer
              print *,'  H = ',H(:)
              print *,'exiting...'
              stop
            endif
          endif

          ! propagate to this point - acoustic layers
          Qmat_I(:,:) = 0.0_CUSTOM_CMPLX
          Qmat_I(1,1) = cmplx(1.0,0.0,kind=CUSTOM_CMPLX)
          Qmat_I(2,2) = cmplx(1.0,0.0,kind=CUSTOM_CMPLX)
          do j = ilayer_ac,ilayer,-1
            if (j > ilayer) then
              call fk_propagator_ac(om,eta_alpha(j),rho(j),H(j),ray_p,Qmat(:,:))
            else
              call fk_propagator_ac(om,eta_alpha(j),rho(j),height,ray_p,Qmat(:,:))
            endif
            ! resulting matrix
            Qmat_I = matmul(Qmat,Qmat_I)
          enddo
          Qmat = Qmat_I

          ! propagation to this point - elastic layers
          N_mat(:,:) = E_mat(:,:)
          do j = nlayer-1,ilayer_ac+1,-1
            call fk_propagator_psv(om,eta_alpha(j),eta_beta(j),rho(j), &
                                   vs(j),H(j),ray_p,gamma1(j),Pmat)
            ! resulting matrix
            N_mat = gamma0(j) * matmul(Pmat,N_mat)
          enddo

          !! zz(ipt) is the height of point with respect to the lower layer
          bot_vec = matmul(N_mat,bot_vec)

          bot_vec(4) = -bot_vec(4) ! P = - \sigma_zz
          bot_vec(1) = bot_vec(2) ! uz
          bot_vec(2) = bot_vec(4) ! P / k
          bot_vec(1:2) = matmul(Qmat,bot_vec(1:2))

          ! make sure free boundary condition is satisfied
          if (abs(height - H(1)) < 1.0e-6 * H(1)) bot_vec(2) = 0.

          ! for acoustic wave we cache dchi/displ
          ! ux/uz
          dz_f = bot_vec(1) ! uz
          dx_f = cmplx(0,-1) * ray_p**2 / rho(ilayer) * bot_vec(2)  ! ux = -ik P / (rho om^2)

          ! for the Stacey boundary contribution, we need displ
          field_f(ii,1) = stf_coeff * dx_f              ! u_x
          field_f(ii,2) = stf_coeff * dz_f              ! u_z

          ! dchi
          ! we set 3:5 the same value to match the api: FFTINV
          field_f(ii,3:5) = cmplx(0.,1.) * bot_vec(2) * stf_coeff * ray_p ! dchi = - P / (i om) = i vec(2) * rayp
        endif ! end ispec_is_elastic(ispec)

      enddo ! end ii for freq

      ! pad negative f, and convert to time series
      do ii = 2, nf2-1
        field_f(nf+2-ii,:) = conjg(field_f(ii,:))
      enddo

      ! inverse fast fourier transform
      field(:,:) = 0.0
      do j = 1, nvar
        ! inverse FFT
        call FFTinv(npow,field_f(:,j),zign_neg,dt,field(:,j),mpow)

        ! wrap around to start from t0: here one has to be careful if t0/dt is not
        ! exactly an integer, assume nn > 0
        if (nn > 0) then
          dtmp(1:nn) = field(npts2-nn+1:npts2,j)
          field(nn+1:npts2,j) = field(1:npts2-nn,j)
          field(1:nn,j) = dtmp(1:nn)
        else if (nn < 0) then
          dtmp(1:nn) = field(1:nn,j)
          field(1:npts-nn,j) = field(nn+1:npts,j)
          field(npts-nn+1:npts,j) = dtmp(1:nn)
        endif
      enddo

      !! store undersampled version of velocity  FK solution
      if (ispec_is_elastic(ispec)) then
        ! elastic element
        tmp_t1(:) = field(:,1) * cos(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        Veloc_FK(1,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,1) * sin(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        Veloc_FK(2,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,2)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        Veloc_FK(3,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

      else if (ispec_is_acoustic(ispec)) then
        ! acoustic element
        ! nqdu comment: veloc_FK now saves displ instead of velocity
        !               and tract_FK saves chi_dot instead of traction
        tmp_t1(:) = field(:,1) * cos(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        veloc_FK(1,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,1) * sin(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        veloc_FK(2,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,2)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        veloc_FK(3,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        ! traction -> chi_dot
        tmp_t1(:) = field(:,3)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        do j = 1,3
          Tract_FK(j,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)
        enddo
      endif

      !! compute traction
      if (comp_stress .and. ispec_is_elastic(ispec)) then
        do lpts = 1, NF_FOR_STORING
          sigma_rr = field(lpts,3)
          sigma_rt = 0.0
          sigma_rz = field(lpts,4)
          sigma_zz = field(lpts,5)
          sigma_tt = bdlambdamu(ipt) * (sigma_rr + sigma_zz)
          sigma_tz = 0.0

          Txx_tmp = sigma_rr * cos(phi) * cos(phi) + sigma_tt * sin(phi) * sin(phi)
          Txy_tmp = cos(phi) * sin(phi) * (sigma_rr - sigma_tt)
          Txz_tmp = sigma_rz * cos(phi)
          Tyy_tmp = sigma_rr * sin(phi) * sin(phi) + sigma_tt * cos(phi) * cos(phi)
          Tyz_tmp = sigma_rz * sin(phi)
          Tzz_tmp = sigma_zz

          !! store directly the traction
          Tract_FK(1,ipt,lpts) = Txx_tmp * nmx(ipt) +  Txy_tmp * nmy(ipt) +  Txz_tmp * nmz(ipt)
          Tract_FK(2,ipt,lpts) = Txy_tmp * nmx(ipt) +  Tyy_tmp * nmy(ipt) +  Tyz_tmp * nmz(ipt)
          Tract_FK(3,ipt,lpts) = Txz_tmp * nmx(ipt) +  Tyz_tmp * nmy(ipt) +  Tzz_tmp * nmz(ipt)
        enddo

        !! store undersamped version of tractions FK solution
        tmp_t1(1:NF_FOR_STORING) = Tract_FK(1,ipt,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        Tract_FK(1,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING) = Tract_FK(2,ipt,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        Tract_FK(2,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING) = Tract_FK(3,ipt,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2, tmp_c)
        Tract_FK(3,ipt,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)
      endif

      ! user output
      if (myrank == 0 .and. npt > 1000) then
        if (mod(ipt,(npt/10)) == 0) then
          write(IMAIN,*) '  done ',ipt/(npt/10)*10,'% points out of ',npt
          call flush_IMAIN()
        endif
      endif

    enddo
  enddo

  ! taper
  if (USE_TAPERED_BEGINNING) then
    ! check if taper length is reasonable compared to wavefield storage traces
    if (NF_FOR_STORING > 2 * taper_nlength) then
      do i = 1,taper_nlength
        ! cosine taper, otherwise using a constant (1.0) instead
        taper = (1.0 - cos(PI*(i-1)/taper_nlength)) * 0.5        ! between [0,1[

        ! tapers traces
        Veloc_FK(:,:,i) = taper * Veloc_FK(:,:,i)
        if (comp_stress) then
          Tract_FK(:,:,i) = taper * Tract_FK(:,:,i)
        endif
      enddo
    endif
  endif

  ! free temporary arrays
  deallocate(fvec,coeff, field_f, field, dtmp)
  deallocate(tmp_f1, tmp_f2, tmp_f3, tmp_t1, tmp_t2, tmp_t3)
  deallocate(tmp_c)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "  FK computing passed "
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine FK

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fk_propagator_psv(om,eta_alpha,eta_beta,rho,vs,H,ray_p,gamma1,Pmat)

  use specfem_par, only: CUSTOM_REAL
  implicit none

  integer, parameter                       :: CUSTOM_CMPLX = 8
  real(kind=CUSTOM_REAL),intent(in)        :: om,rho,vs,ray_p,H
  complex(kind=CUSTOM_CMPLX),intent(in)    :: eta_beta,eta_alpha,gamma1
  complex(kind=CUSTOM_CMPLX),intent(out)   :: Pmat(4,4)
  complex(kind=CUSTOM_CMPLX)               :: c1,ca,sa,xa,ya,c2,cb,sb,xb,yb
  complex(kind=CUSTOM_CMPLX)               :: g1,g1_sq,mul,nu_al,nu_be,two_mul

  ! compute propagation matrix

  ! factor nu_al = -i nu_p = -i omega sqrt(1/alpha^2 - p^2)
  !                        = omega (-i sqrt(1/alpha^2 - p^2)
  !                        = omega (eta_alpha)
  nu_al = om * eta_alpha

  ! factor nu_be = -i nu_s = -i omega sqrt(1/beta^2 - p^2)
  !                        = omega (-i sqrt(1/beta^2 - p^2)
  !                        = omega (eta_beta)
  nu_be = om * eta_beta

  ! matrix variables
  ! C_alpha = cos(nu_p h) = cos( omega sqrt(1/alpha^2 - p^2) h )
  ! S_alpha = -sin(nu_p h) = -sin( omega sqrt(1/alpha^2 - p^2) h )
  ! with nu_p h = omega sqrt(1/alpha^2 - p^2) h
  !
  ! note: there might be some sign conflict and confusion between sin and sinh in the definition after (A9) of Tong et al.
  !       instead of S_alpha = -sin(nu_p h) as stated in the paper, here sa becomes [i sin(nu_p h)] as imaginary number.
  !
  c1 = nu_al * H
  ca = (exp(c1) + exp(-c1))/2.0
  sa = (exp(c1) - exp(-c1))/2.0     ! imaginary part
  xa = eta_alpha * sa / ray_p
  ya = ray_p * sa / eta_alpha

  ! C_beta = cos(nu_s h) = cos( omega sqrt(1/beta^2 - p^2) h )
  ! S_beta = -sin(nu_s h) = -sin( omega sqrt(1/beta^2 - p^2) h )
  ! with nu_s h = omega sqrt(1/beta^2 - p^2) h
  c2 = nu_be * H
  cb = (exp(c2) + exp(-c2))/2.0  ! cos(nu_s h)
  sb = (exp(c2) - exp(-c2))/2.0  ! sin(nu_s h) imaginary part
  xb = eta_beta * sb / ray_p
  yb = ray_p * sb / eta_beta

  ! layer factors
  g1 = gamma1
  mul = cmplx(rho * vs * vs ,0.0)
  g1_sq = g1 * g1
  two_mul = 2. * mul

  ! Tong et al. (2014), appendix (A8)
  ! propagation matrix P_n
  Pmat(1,1) = ca - g1*cb
  Pmat(1,2) = xb - g1*ya
  !org: P_layer(1,3) = (ya - xb)/(2*mul)          ! misses factor 1/k
  !     P_layer(1,4) = (cb - ca)/(2*mul)          ! misses factor 1/k
  Pmat(1,3) = (ya - xb) / two_mul
  Pmat(1,4) = (cb - ca) / two_mul

  Pmat(2,1) = xa - g1*yb
  Pmat(2,2) = cb - g1*ca
  !org: P_layer(2,3) = (ca - cb)/(2*mul)          ! misses factor 1/k
  !     P_layer(2,4) = (yb - xa)/(2*mul)          ! misses factor 1/k
  Pmat(2,3) = (ca - cb) / two_mul
  Pmat(2,4) = (yb - xa) / two_mul

  !org: P_layer(3,1) = 2*mul * (xa - g1**2 * yb)  ! misses factor k
  !     P_layer(3,2) = 2*mul * g1 * (cb - ca)     ! misses factor k
  Pmat(3,1) = two_mul * (xa - g1_sq * yb)
  Pmat(3,2) = two_mul * g1 * (cb - ca)
  Pmat(3,3) = ca - g1*cb
  Pmat(3,4) = g1*yb - xa

  !org: P_layer(4,1) = 2*mul * g1 * (ca - cb)     ! misses factor k
  !     P_layer(4,2) = 2*mul * (xb - g1**2 * ya)  ! misses factor k
  Pmat(4,1) = two_mul * g1 * (ca - cb)
  Pmat(4,2) = two_mul * (xb - g1_sq * ya)
  Pmat(4,3) = g1*ya - xb
  Pmat(4,4) = cb - g1*ca

  end subroutine fk_propagator_psv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fk_propagator_ac(om,eta_alpha,rho,H,ray_p,Qmat)

  use specfem_par, only: CUSTOM_REAL
  implicit none

  integer, parameter                       :: CUSTOM_CMPLX = 8
  real(kind=CUSTOM_REAL),intent(in)        :: om , rho, ray_p,H
  complex(kind=CUSTOM_CMPLX),intent(in)    :: eta_alpha
  complex(kind=CUSTOM_CMPLX),intent(out)   :: Qmat(2,2)
  complex(kind=CUSTOM_CMPLX)               :: c1,ca,sa
  complex(kind=CUSTOM_CMPLX)               :: nu_al

  ! factor nu_al = -i nu_p = -i omega sqrt(1/alpha^2 - p^2)
  !                        = omega (-i sqrt(1/alpha^2 - p^2)
  !                        = omega (eta_alpha)
  nu_al = om * eta_alpha

  ! matrix variables
  ! C_alpha = cos(nu_p h) = cos( omega sqrt(1/alpha^2 - p^2) h )
  ! S_alpha = -sin(nu_p h) = -sin( omega sqrt(1/alpha^2 - p^2) h )
  ! with nu_p h = omega sqrt(1/alpha^2 - p^2) h
  c1 = nu_al * H
  ca = (exp(c1) + exp(-c1))/2.0
  sa = (exp(c1) - exp(-c1))/2.0     ! imaginary part

  Qmat(1,1) = ca
  Qmat(2,2) = ca
  Qmat(1,2) = -sa * eta_alpha * ray_p / rho! missing k
  Qmat(2,1) = sa * rho / (ray_p  * eta_alpha)

  end subroutine fk_propagator_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_size_of_working_arrays(deltat, freq_sampling_fk, tmax_fk, NF_FOR_STORING, &
                                         NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, DF_FK)

  use constants, only: CUSTOM_REAL
  use specfem_par, only: NSTEP

  implicit none
  real(kind=CUSTOM_REAL),intent(in)    :: deltat,freq_sampling_fk
  real(kind=CUSTOM_REAL),intent(inout) :: tmax_fk
  real(kind=CUSTOM_REAL),intent(inout) :: DF_FK
  integer,               intent(inout) :: NF_FOR_STORING, NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP
  ! local parameter
  real(kind=CUSTOM_REAL)               :: df, dt_min_fk, Frq_ech_Fk
  real(kind=CUSTOM_REAL)               :: tmax_samp, tmax_use

  ! sampling frequency to store fk solution
  Frq_ech_Fk = freq_sampling_fk

  ! sampling time step to store fk solution
  dt_min_fk = 1.0 /  Frq_ech_Fk

  ! compute resampling rate
  NP_RESAMP  = floor(dt_min_fk / deltat)

  ! checks sampling rate
  if (NP_RESAMP == 0) NP_RESAMP = 1

  ! update dt for fk with respect to integer np_resampling
  dt_min_fk = NP_RESAMP * deltat  !! this is the time step sampling for FK storage

  !! time window lengths
  ! maximum simulation time length would be:
  !   tmax_simulation = NSTEP * deltat
  ! FK signal trace is symmetric after inverse FFT, thus NF_FOR_STORING must match at least match NSTEP * 2 with resampling rate.
  ! maximum time length with respect to simulation NSTEP and resampling rate:
  tmax_samp = (2*NSTEP/NP_RESAMP + 1) * dt_min_fk

  ! take user defined window length as initial default
  tmax_use = tmax_fk

  ! limits window length to actual simulation length
  if (tmax_use > tmax_samp) tmax_use = tmax_samp

  ! compute number of time steps to store
  NF_FOR_STORING  = ceiling( tmax_use / dt_min_fk)

  ! in power of two
  NF_FOR_STORING  =   ceiling(log(real(NF_FOR_STORING))/log(2.))

  ! multiply by 2 in order to do an inverse FFT
  NF_FOR_FFT      =   2**(NF_FOR_STORING+1)

  NPOW_FOR_INTERP =   NF_FOR_STORING+1
  NF_FOR_STORING  =   2**NF_FOR_STORING

  ! now we have this new time window
  tmax_fk = dt_min_fk * (NF_FOR_FFT - 1)

  ! step in frequency for fk
  df = 1.0 / tmax_fk
  DF_FK = df

  !debug
  !print *,'debug: tmax ',tmax_fk,tmax_samp,tmax_use,'np_resamp',NP_RESAMP,'NSTEP',NSTEP,2*NSTEP/NP_RESAMP, &
  !                      'NF',NF_FOR_STORING,ceiling( tmax_use / dt_min_fk)

  end subroutine find_size_of_working_arrays


!
!-------------------------------------------------------------------------------------------------
!

!! #################  INTERPOLATION ROUTINES IN TIME DOMAIN ######################################

!! compute and store spline coefficients

  subroutine compute_spline_coef_to_store(Sig, npts, spline_coeff, tmp_c)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,                                             intent(in)     :: npts
  real(kind=CUSTOM_REAL), dimension(npts),             intent(in)     :: Sig
  real(kind=CUSTOM_REAL), dimension(npts),             intent(inout)  :: spline_coeff
  double precision, dimension(npts),intent(out)                       :: tmp_c  ! temporary work array

  !! computation in double precision
  double precision,parameter                                          :: error = 1.d-24
  double precision,parameter                                          :: z1 = dsqrt(3.d0) - 2.d0
  double precision                                                    :: zn, sumc
  integer                                                             :: i, n_init

  ! Compute pole value
  tmp_c(:) = dble(Sig(:)) * (1.d0 - z1) * ( 1.d0 - 1.d0/z1)

  ! Initialisation causal filter
  n_init = ceiling(log(error)/log(abs(z1)))

  ! check limits: by default is n_init==42, for very short test simulations this might be larger than npts
  if (n_init > npts) n_init = npts

  sumc = tmp_c(1)
  zn = z1
  do i = 1,n_init
     sumc = sumc + zn * tmp_c(i)
     zn = zn * z1
  enddo
  tmp_c(1) = sumc

  ! Causal filter
  do i = 2,npts
     tmp_c(i) = tmp_c(i) + z1 * tmp_c(i-1)
  enddo

  ! Initialisation anti-causal filter
  tmp_c(npts) = ( z1 / (z1 - 1.d0) ) * tmp_c(npts)
  do i = npts-1,1,-1
     tmp_c(i) = z1 * (tmp_c(i+1) - tmp_c(i))
  enddo

  !! store spline coeff in CUSTOM_REAL precision
  spline_coeff(:) = tmp_c(:)

  end subroutine compute_spline_coef_to_store


!
!-------------------------------------------------------------------------------------------------
!

!! VM VM READING INPUT FILE FOR FK MODEL

  subroutine ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

  use specfem_par
  use specfem_par_coupling

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
  integer                :: ioerr
  character(len=100)     :: keyword, keyvalue, line
  character(len=100)     :: keyword_tmp, incident_wave

  real(kind=CUSTOM_REAL) :: rho_layer, vp_layer, vs_layer, ztop_layer
  real(kind=CUSTOM_REAL) :: Radius_box, wave_length_at_bottom
  real(kind=CUSTOM_REAL), dimension(:), allocatable  :: rho_fk_input, vp_fk_input, vs_fk_input, ztop_fk_input
  integer,  dimension(:), allocatable  :: ilayer_fk_input
  integer  :: ilayer,ier
  logical  :: position_of_wavefront_not_read

  !!--------------------------------------------------------------
  ! # model description :
  ! NLAYER   n # number of layers
  ! LAYER 1  rho, vp ,vs, ztop
  ! LAYER 2  rho, vp, vs, ztop
  ! ...
  ! LAYER n  rho, vp, vs, ztop  # homogenoeus half-space
  ! #
  ! # incident wave description:
  ! INCIDENT_WAVE  "p" or "sv"
  ! BACK_AZITUTH    bazi
  ! INCIDENCE       inc
  ! ORIGIN_WAVEFRONT xx0, yy0, zz0
  ! ORIGIN_TIME      tt0
  ! FREQUENCY_MAX    ff0
  ! FREQUENCY_SAMPLING freq_sampling_fk
  ! TIME_WINDOW      tmax_fk
  ! AMPLITUDE        amplitude_fk
  !!----------------------------------------------------------------

  ! only main process reads
  if (myrank /= 0) return

  !! set default values
  tt0 = 0.0                 ! origin time
  tmax_fk = 128.0           ! time window length

  ff0 = 0.1                 ! frequency maximum in Hz
  freq_sampling_fk = 10.0   ! frequency sampling rate in Hz
  amplitude_fk = 1.0        ! amplitude factor

  type_kpsv_fk = 1  ! 1 == P-wave / 2 == SV-wave

  position_of_wavefront_not_read = .true.

  !! READING input file
  open(85,file=trim(FKMODEL_FILE))
  do
     read(85, fmt='(a100)',iostat=ioerr) line
     !call remove_begin_blanks(line)
     if (ioerr < 0) exit
     if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
     read(line,*) keyword, keyvalue

     select case(trim(keyword))
     case('NLAYER')
        read(line, *) keyword_tmp, nlayer

        allocate(alpha_FK(nlayer), beta_FK(nlayer), rho_FK(nlayer), mu_FK(nlayer), h_FK(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating arrays 2226')
        alpha_FK(:) = 0._CUSTOM_REAL; beta_FK(:) = 0._CUSTOM_REAL; rho_FK(:) = 0._CUSTOM_REAL
        mu_FK(:) = 0._CUSTOM_REAL; h_FK(:) = 0._CUSTOM_REAL

        allocate(rho_fk_input(nlayer), &
                 vp_fk_input(nlayer), &
                 vs_fk_input(nlayer), &
                 ztop_fk_input(nlayer+1), &
                 ilayer_fk_input(nlayer+1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating arrays 2227')
        rho_fk_input(:) = 0._CUSTOM_REAL
        vp_fk_input(:) = 0._CUSTOM_REAL
        vs_fk_input(:) = 0._CUSTOM_REAL
        ztop_fk_input(:) = 0._CUSTOM_REAL
        ilayer_fk_input(:) = -1

     case('LAYER')
        read(line, *) keyword_tmp, ilayer, rho_layer, vp_layer, vs_layer, ztop_layer

        ilayer_fk_input(ilayer) = ilayer
        rho_fk_input(ilayer) = rho_layer
        vp_fk_input(ilayer) = vp_layer
        vs_fk_input(ilayer) = vs_layer
        ztop_fk_input(ilayer) = ztop_layer

        ! acoustic case: FK for zero-shear velocity not handled, here we set it to a (very) small value
        !if (abs(vs_fk_input(ilayer)) < TINYVAL) vs_fk_input(ilayer) = 1.e-5

     case('INCIDENT_WAVE')
         read(line,*)  keyword_tmp, incident_wave

         select case(trim(incident_wave))
            case ('p', 'P')
              type_kpsv_fk = 1
            case('sv','SV')
              type_kpsv_fk = 2
            case default
              type_kpsv_fk = 1
            end select

     case('BACK_AZIMUTH')
        read(line,*)  keyword_tmp, phi_FK
        phi_FK = - phi_FK - 90.0

     case('AZIMUTH')
        read(line,*)  keyword_tmp, phi_FK
        phi_FK = 90.0 - phi_FK

     case('TAKE_OFF')
        read(line,*)  keyword_tmp, theta_FK

     case('ORIGIN_WAVEFRONT')
         read(line,*)  keyword_tmp, xx0, yy0, zz0
         position_of_wavefront_not_read = .false.

      case('ORIGIN_TIME')
         read(line,*)  keyword_tmp, tt0

     case('FREQUENCY_MAX')
        read(line,*)  keyword_tmp, ff0

     case('FREQUENCY_SAMPLING')
        read(line,*)  keyword_tmp, freq_sampling_fk

     case('TIME_WINDOW')
        read(line,*)  keyword_tmp, tmax_fk

      case('AMPLITUDE')
        read(line,*)  keyword_tmp, amplitude_fk

     end select
  !!------------------------------------------------------------------------------------------------------
  enddo
  close(85)

  if (allocated(ilayer_fk_input)) then

     ilayer_fk_input(nlayer+1) = ilayer_fk_input(nlayer)
     ztop_fk_input(nlayer+1) = ztop_fk_input(nlayer)

     ! Z coordinate reference for top of lower-halfspace (last layer input)
     Z_REF_for_FK = ztop_fk_input(nlayer)

     do ilayer = 1, nlayer
        alpha_FK(ilayer) = vp_fk_input(ilayer)
        beta_FK(ilayer) = vs_fk_input(ilayer)
        rho_FK(ilayer) = rho_fk_input(ilayer)

        mu_FK(ilayer) = rho_fk_input(ilayer) * vs_fk_input(ilayer)**2
        h_FK(ilayer) =  ztop_fk_input(ilayer) - ztop_fk_input(ilayer+1)

        ! checks if layer has values
        if (ilayer_fk_input(ilayer) == -1) then
           write(*,*) " ERROR READING FK INPUT FILE "
           write(*,*) " MISSING LAYER ", ilayer
           stop 'Missing layer in FK model'
        endif

        ! checks vp,rho are strictly positive
        ! if (alpha_FK(ilayer) <= 0.0_CUSTOM_REAL .or. rho_FK(ilayer) <= 0.0_CUSTOM_REAL) then
        !   print *,'Error: invalid elastic material for FK coupling in layer ',ilayer
        !   print *,'  vp = ',alpha_FK(ilayer),' vs = ',beta_FK(ilayer),' rho = ',rho_FK(ilayer)
        !   print *,'vp & rho must be strictly positive'
        !   stop 'Invalid material for FK coupling found in ReadFKModelInput() routine'
        ! endif

        ! checks if vs is strictly positive for SV coupling
        ! if (type_kpsv_fk == 2) then
        !   if (beta_FK(ilayer) <= 0.0_CUSTOM_REAL) then
        !     print *,'Error: invalid elastic material for FK coupling in layer ',ilayer
        !     print *,'  vp = ',alpha_FK(ilayer),' vs = ',beta_FK(ilayer),' rho = ',rho_FK(ilayer)
        !     print *,'vp, vs & rho must be strictly positive (in particular vs must be > 0) for SV incident wave'
        !     stop 'Invalid elastic material for FK coupling found in ReadFKModelInput() routine'
        !   endif
        ! endif

     enddo

     deallocate(ilayer_fk_input, rho_fk_input, vp_fk_input, vs_fk_input, ztop_fk_input)

  else

     write(*,*) " ERROR READING FK INPUT FILE "
     write(*,*) " NOT BE ABLE TO READ MODEL PROPERTIES "
     stop 'Error reading FK input file'

  endif

  !! compute position of wave front
  if (position_of_wavefront_not_read) then
    ! sets center point of box
    xx0 = 0.5*(Xmin_box + Xmax_box)
    yy0 = 0.5*(Ymin_box + Ymax_box)

    Radius_box = sqrt( (Xmin_box - xx0)**2 + (Ymin_box - yy0)**2)

    if (type_kpsv_fk == 1) then
      ! P-wave
      wave_length_at_bottom = alpha_FK(nlayer) / ff0    ! vp
    else if (type_kpsv_fk == 2) then
      ! SV-wave
      wave_length_at_bottom = beta_FK(nlayer) / ff0     ! vs
    endif

    ! depth position below bottom
    zz0 = Zmin_box - Radius_box * sin ( abs(theta_FK) * (PI / 180.d0) )  &
           - 3.0 * wave_length_at_bottom * cos ( abs(theta_FK) * (PI / 180.d0) )
  endif

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) "**********************************************"
  write(IMAIN,*) "         USING FK INJECTION TECHNIQUE         "
  write(IMAIN,*) "**********************************************"
  write(IMAIN,*)
  write(IMAIN,*) "  Model : " , nlayer , " layers "
  write(IMAIN,*)
  do ilayer = 1, nlayer
     write(IMAIN,'(a9, i3, 3(a6,2x,f8.3), 3x, a9, f18.5 )') &
          "   layer ", ilayer, &
          " rho =",   rho_FK(ilayer), &
          " vp  =",   alpha_FK(ilayer), &
          " vs  =",   beta_FK(ilayer), &
          " Height =", h_FK(ilayer)
  enddo
  write(IMAIN,*)
  write(IMAIN,*) "  FK (azimuth) angle phi   = ",   phi_FK,'(deg)'
  write(IMAIN,*) "  FK  take-off angle theta = ", theta_FK,'(deg)'
  write(IMAIN,*)
  write(IMAIN,*) "  Origin wavefront point FK           : ", xx0, yy0, zz0
  write(IMAIN,*) "  Time shift  FK                      : ", tt0
  write(IMAIN,*) "  Z reference for FK routine          : ", Z_REF_for_FK
  write(IMAIN,*)
  write(IMAIN,*) "  Time window for FK computing        : ", tmax_fk
  write(IMAIN,*) "  Frequency max                       : ", ff0
  write(IMAIN,*) "  Frequency sampling rate             : ", freq_sampling_fk
  write(IMAIN,*) "  Amplitude factor                    : ", amplitude_fk
  write(IMAIN,*)
  write(IMAIN,*) "  Model boundary box                  : Xmin/Xmax = ",Xmin_box,Xmax_box
  write(IMAIN,*) "                                      : Ymin/Ymax = ",Ymin_box,Ymax_box
  write(IMAIN,*) "                                      : Zmin/Zmax = ",Zmin_box,Zmax_box
  write(IMAIN,*)
  write(IMAIN,*) "  Type of incoming wave (1=P), (2=SV) : ",type_kpsv_fk
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine ReadFKModelInput

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

  use specfem_par

  implicit none

  real(kind=CUSTOM_REAL), intent(in out) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
  real(kind=CUSTOM_REAL)                 :: Xmin_loc, Xmax_loc, Ymin_loc, Ymax_loc, Zmin_loc, Zmax_loc

  Xmin_loc = minval (xstore(:))
  Xmax_loc = maxval (xstore(:))
  Ymin_loc = minval (ystore(:))
  Ymax_loc = maxval (ystore(:))
  Zmin_loc = minval (zstore(:))
  Zmax_loc = maxval (zstore(:))

  call min_all_all_cr(Xmin_loc, Xmin_box)
  call max_all_all_cr(Xmax_loc, Xmax_box)
  call min_all_all_cr(Ymin_loc, Ymin_box)
  call max_all_all_cr(Ymax_loc, Ymax_box)
  call min_all_all_cr(Zmin_loc, Zmin_box)
  call max_all_all_cr(Zmax_loc, Zmax_box)

  end subroutine FindBoundaryBox


!-------------------------------------------------------------------------------------------------
!
! routines for SPECFEM injection technique
!
!-------------------------------------------------------------------------------------------------

! routines for storing wavefield solution on coupling boundary points

  subroutine locate_coupling_points(num_coupling_points_total,coupling_points)

! locates all coupling points (given in specfem_coupling_points.bin file) in current mesh
!
! similar as in routine locate_receivers()

  use constants, only: HUGEVAL,CUSTOM_REAL,IMAIN,NDIM

  use specfem_par, only: ibool,myrank,NSPEC_AB,NGLOB_AB,NPROC, &
    xstore,ystore,zstore

  use specfem_par_coupling

  implicit none

  integer, intent(in) :: num_coupling_points_total
  real(kind=CUSTOM_REAL),dimension(6,num_coupling_points_total),intent(in) :: coupling_points

  ! local parameters
  integer :: ipoin,islice,ispec,ier

  ! temporary arrays
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, dimension(:), allocatable :: final_distance
  double precision :: final_distance_max
  integer, dimension(:),allocatable :: idomain

  double precision, dimension(:), allocatable :: xi_point,eta_point,gamma_point
  double precision, dimension(:,:,:), allocatable :: nu_point

  ! location search
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  double precision :: x,y,z,x_new,y_new,z_new
  double precision :: xi,eta,gamma,final_distance_squared
  double precision, dimension(NDIM,NDIM) :: nu_found

  integer :: ispec_found,idomain_found

  ! subset size
  integer, parameter :: NPOIN_SUBSET_MAX = 500

  ! subset arrays
  double precision, dimension(NPOIN_SUBSET_MAX) :: xi_poin_subset,eta_poin_subset,gamma_poin_subset
  double precision, dimension(NPOIN_SUBSET_MAX) :: x_found_subset,y_found_subset,z_found_subset
  double precision, dimension(NPOIN_SUBSET_MAX) :: final_distance_subset
  double precision, dimension(NDIM,NDIM,NPOIN_SUBSET_MAX) :: nu_subset
  integer, dimension(NPOIN_SUBSET_MAX) :: ispec_selected_poin_subset,idomain_subset
  integer :: npoin_subset_current_size,ipoin_in_this_subset,ipoin_already_done
  integer :: num_output_info

  ! timer MPI
  double precision, external :: wtime
  double precision :: tstart,tCPU

  ! get MPI starting time
  tstart = wtime()

  ! initializes
  npoints_total = num_coupling_points_total

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  locating coupling points:'
    call flush_IMAIN()

    ! output frequency for large number of points
    ! number to output in steps of 1000/10000/.. depending on how large the number of entries is
    if (npoints_total > 500000) then
        num_output_info = min(100000,int(10**floor(log10(dble(npoints_total)))))
    else if (npoints_total > 50000) then
        num_output_info = min(10000,int(10**floor(log10(dble(npoints_total)))))
    else
        num_output_info = min(1000,int(10**floor(log10(dble(npoints_total)))))
    endif
    ! number to output about ~50 steps, rounds to the next multiple of 1000
    !num_output_info = max(1000,int(ceiling(ceiling(npoints_total/50.0)/1000.0)*1000))
  endif

  ! compute typical size of elements
  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! allocate memory for point arrays, i.e. points given in specfem_coupling_points.bin file
  allocate(islice_selected_point(npoints_total), &
           ispec_selected_point(npoints_total),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for coupling points'
  ! initializes arrays
  islice_selected_point(:) = -1
  ispec_selected_point(:) = 0

  ! allocate memory for additional arrays using number of stations
  allocate(x_found(npoints_total), &
           y_found(npoints_total), &
           z_found(npoints_total), &
           xi_point(npoints_total), &
           eta_point(npoints_total), &
           gamma_point(npoints_total), &
           nu_point(NDIM,NDIM,npoints_total), &
           final_distance(npoints_total), &
           idomain(npoints_total),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for coupling points x_found,..'
  x_found(:) = 0.d0; y_found(:) = 0.d0; z_found(:) = 0.d0
  xi_point(:) = 0.d0; eta_point(:) = 0.d0; gamma_point(:) = 0.d0; nu_point(:,:,:) = 0.0d0
  final_distance(:) = HUGEVAL
  idomain(:) = 0

  ! note: we loop over subsets of points to reduce the MPI communication for each point
  !
  ! loop on all the stations to locate the stations
  do ipoin_already_done = 0, npoints_total, NPOIN_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    npoin_subset_current_size = min(NPOIN_SUBSET_MAX, npoints_total - ipoin_already_done)

    ! initializes search results
    final_distance_subset(:) = HUGEVAL

    ! loop over the stations within this subset
    do ipoin_in_this_subset = 1,npoin_subset_current_size

      ! mapping from station number in current subset to real station number in all the subsets
      ipoin = ipoin_in_this_subset + ipoin_already_done

      x = coupling_points(1,ipoin)
      y = coupling_points(2,ipoin)
      z = coupling_points(3,ipoin)

      ! locates point in mesh
      call locate_point_in_mesh(x, y, z, &
                                .true., elemsize_max_glob, &
                                ispec_found, xi, eta, gamma, &
                                x_new, y_new, z_new, &
                                idomain_found, nu_found, final_distance_squared)

      ispec_selected_poin_subset(ipoin_in_this_subset) = ispec_found

      x_found_subset(ipoin_in_this_subset) = x_new
      y_found_subset(ipoin_in_this_subset) = y_new
      z_found_subset(ipoin_in_this_subset) = z_new

      xi_poin_subset(ipoin_in_this_subset) = xi
      eta_poin_subset(ipoin_in_this_subset) = eta
      gamma_poin_subset(ipoin_in_this_subset) = gamma

      idomain_subset(ipoin_in_this_subset) = idomain_found
      nu_subset(:,:,ipoin_in_this_subset) = nu_found(:,:)
      final_distance_subset(ipoin_in_this_subset) = sqrt(final_distance_squared)

      ! user output progress
      if (myrank == 0) then
        if (mod(ipoin,num_output_info) == 0) then
          tCPU = wtime() - tstart
          write(IMAIN,*) '    located point ',ipoin,'out of',npoints_total,' - elapsed time: ',sngl(tCPU),'s'
          call flush_IMAIN()
        endif
      endif
    enddo ! loop over subset

    ! main process locates best location in all slices
    call locate_MPI_slice(npoin_subset_current_size,ipoin_already_done, &
                          ispec_selected_poin_subset, &
                          x_found_subset, y_found_subset, z_found_subset, &
                          xi_poin_subset,eta_poin_subset,gamma_poin_subset, &
                          idomain_subset,nu_subset,final_distance_subset, &
                          npoints_total,ispec_selected_point, islice_selected_point, &
                          x_found,y_found,z_found, &
                          xi_point, eta_point, gamma_point, &
                          idomain,nu_point,final_distance)

  enddo ! loop over stations

  ! bcast from main process
  call bcast_all_i(islice_selected_point,npoints_total)

  ! note: in principle, only islice must be updated on all secondary processes, the ones containing the best location
  !       could have valid entries in all other arrays set before already.
  !       nevertheless, for convenience we broadcast all needed arrays back to the secondarys
  call bcast_all_i(idomain,npoints_total)
  call bcast_all_i(ispec_selected_point,npoints_total)

  call bcast_all_dp(xi_point,npoints_total)
  call bcast_all_dp(eta_point,npoints_total)
  call bcast_all_dp(gamma_point,npoints_total)

  call bcast_all_dp(x_found,npoints_total)
  call bcast_all_dp(y_found,npoints_total)
  call bcast_all_dp(z_found,npoints_total)

  call bcast_all_dp(nu_point,NDIM*NDIM*npoints_total)
  call bcast_all_dp(final_distance,npoints_total)

  ! checks if we got valid elements
  ! locate point might return a zero ispec if point outside/above mesh
  do ipoin = 1,npoints_total
    ! slice the point is in
    islice = islice_selected_point(ipoin)

    ! checks slice
    if (islice < 0 .or. islice > NPROC-1) then
      print *,'Error locating point # ',ipoin
      print *,'  found in an invalid slice: ',islice
      call exit_MPI(myrank,'Error something is wrong with the slice number of coupling point')
    endif

    ! checks found element
    if (myrank == islice_selected_point(ipoin)) then
      ! element index
      ispec = ispec_selected_point(ipoin)
      ! checks if valid
      if (ispec < 1 .or. ispec > NSPEC_AB) then
        ! invalid element
        print *,'Error locating coupling point # ',ipoin
        print *,'  original x: ',coupling_points(1,ipoin)
        print *,'  original y: ',coupling_points(2,ipoin)
        print *,'  original z: ',coupling_points(3,ipoin)
        print *,'  found in an invalid element: slice         :',islice_selected_point(ipoin)
        print *,'                               ispec         :',ispec_selected_point(ipoin),'out of ',NSPEC_AB
        print *,'                               domain        :',idomain(ipoin)
        print *,'                               xi/eta/gamma  :',xi_point(ipoin),eta_point(ipoin),gamma_point(ipoin)
        print *,'                               final_distance:',final_distance(ipoin)
        print *
        print *,'Please check your coupling point positions, and move them inside the (coarse) mesh geometry - exiting...'
        print *
        call exit_MPI(myrank,'Error locating coupling points')
      endif
    endif
  enddo
  call synchronize_all()

  ! this is executed by main process only
  if (myrank == 0) then
    ! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))

    ! checks stations location
    do ipoin = 1,npoints_total
      if (final_distance(ipoin) == HUGEVAL) then
        write(IMAIN,*) 'Error locating coupling point # ',ipoin,' has huge distance'
        call exit_MPI(myrank,'Error locating receiver')
      endif
    enddo

    ! output info for larger distances
    if (final_distance_max > 0.1 * elemsize_min_glob) then
      do ipoin = 1,npoints_total
        if (final_distance(ipoin) > 0.99 * final_distance_max) then
          write(IMAIN,*)
          write(IMAIN,*) '  poor location for coupling point # ',ipoin,' with final_distance: ',final_distance(ipoin)
          write(IMAIN,*) '    original x/y/z = ',coupling_points(1:3,ipoin)
          write(IMAIN,*) '    found at x/y/z = ',sngl(x_found(ipoin)),sngl(y_found(ipoin)),sngl(z_found(ipoin))
          write(IMAIN,*) '    in slice :',islice_selected_point(ipoin),' ispec :',ispec_selected_point(ipoin)
          write(IMAIN,*) '       xi/eta/gamma :',sngl(xi_point(ipoin)),sngl(eta_point(ipoin)),sngl(gamma_point(ipoin))
          call flush_IMAIN()
        endif
      enddo
    endif

    ! display maximum error for all the receivers
    write(IMAIN,*)
    write(IMAIN,*) '  maximum error in location of all the coupling points: ',sngl(final_distance_max),' m'
    write(IMAIN,*)

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
    if (final_distance_max > elemsize_max_glob) then
      write(IMAIN,*) '  ************************************************************'
      write(IMAIN,*) '  ************************************************************'
      write(IMAIN,*) '  ***** WARNING: at least one point is very poorly located ***'
      write(IMAIN,*) '  ************************************************************'
      write(IMAIN,*) '  ************************************************************'
      write(IMAIN,*)
    endif
  endif    ! end of section executed by main process only

  ! sets up point arrays for storing wavefield solution
  call determine_local_coupling_points(num_coupling_points_total,coupling_points, &
                                       xi_point,eta_point,gamma_point, &
                                       x_found,y_found,z_found)

  ! deallocate temporary arrays
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(idomain)
  deallocate(xi_point,eta_point,gamma_point,nu_point)

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - tstart
    write(IMAIN,*) '  Elapsed time for coupling points detection in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine locate_coupling_points

!
!-------------------------------------------------------------------------------------------------

  subroutine determine_local_coupling_points(num_coupling_points_total,coupling_points, &
                                             xi_point,eta_point,gamma_point, &
                                             x_found,y_found,z_found)

! pre-computes arrays for coupling points in this slice

  use constants, only: myrank,CUSTOM_REAL,MAX_STRING_LEN,NDIM,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: xigll,yigll,zigll

  use shared_parameters, only: TRACTION_PATH,NPROC,SAVE_MESH_FILES

  use specfem_par_coupling

  implicit none

  integer, intent(in) :: num_coupling_points_total
  real(kind=CUSTOM_REAL),dimension(6,num_coupling_points_total),intent(in) :: coupling_points

  double precision, dimension(num_coupling_points_total), intent(in) :: xi_point,eta_point,gamma_point
  double precision, dimension(num_coupling_points_total), intent(in) :: x_found,y_found,z_found

  ! local parameters
  integer :: ipoin,ipoin_local,ier

  ! Lagrange interpolators
  real(kind=CUSTOM_REAL) :: hxi(NGLLX),heta(NGLLY),hgamma(NGLLZ)

  ! vtk output
  character(len=MAX_STRING_LEN) :: prname_trac,filename
  integer, dimension(:), allocatable :: iglob_tmp
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_tmp,ystore_tmp,zstore_tmp

  ! initializes
  npoints_local = 0

  ! count local coupling points in this slice
  ipoin_local = 0
  do ipoin = 1,npoints_total
    if (islice_selected_point(ipoin) == myrank) then
      ! counter
      ipoin_local = ipoin_local + 1
    endif
  enddo

  ! sets number of local points
  npoints_local = ipoin_local

  ! for collecting wavefields on main process from each process
  allocate(nb_points_local_per_proc(0:NPROC-1),stat=ier)
  if (ier /= 0) stop 'Error allocating nb_points_local_per_proc array'
  nb_points_local_per_proc(:) = 0

  call gather_all_singlei(npoints_local,nb_points_local_per_proc,NPROC)

  ! for main process only to collect all point wavefield
  if (myrank == 0) then
    ! array for all points
    allocate(veloc_total(NDIM,npoints_total), &
             traction_total(NDIM,npoints_total),stat=ier)
    if (ier /= 0) stop 'Error allocating veloc_total & traction_total array'
    veloc_total(:,:) = 0.0_CUSTOM_REAL
    traction_total(:,:) = 0.0_CUSTOM_REAL

    ! buffer array to collect local wavefields
    allocate(buffer_points_recv(NDIM,maxval(nb_points_local_per_proc(:))),stat=ier)
    if (ier /= 0) stop 'Error allocating buffer_points_recv array'
    buffer_points_recv(:,:) = 0.0_CUSTOM_REAL
  endif

  ! setup local points arrays
  if (npoints_local > 0) then
    ! velocity & traction
    allocate(veloc_points(NDIM,npoints_local), &
             traction_points(NDIM,npoints_local), &
             normal_points(NDIM,npoints_local),stat=ier)
    if (ier /= 0) stop 'Error allocating local veloc & traction array'
    veloc_points(:,:) = 0.0_CUSTOM_REAL
    traction_points(:,:) = 0.0_CUSTOM_REAL
    normal_points(:,:) = 0.0_CUSTOM_REAL

    ! allocate interpolator arrays
    allocate(hxi_point_store(NGLLX,npoints_local), &
             heta_point_store(NGLLY,npoints_local), &
             hgamma_point_store(NGLLZ,npoints_local),stat=ier)
    if (ier /= 0) stop 'Error allocating local interpolator arrays'
    hxi_point_store(:,:) = 0.0_CUSTOM_REAL
    heta_point_store(:,:) = 0.0_CUSTOM_REAL
    hgamma_point_store(:,:) = 0.0_CUSTOM_REAL

    ! get interpolators
    ipoin_local = 0
    do ipoin = 1,npoints_total
      ! pre-computes values for local points
      if (islice_selected_point(ipoin) == myrank) then
        ! calculates interpolators
        ! (1-D Lagrange interpolators)
        call lagrange_only_interpol(NGLLX,xi_point(ipoin),xigll,hxi)
        call lagrange_only_interpol(NGLLY,eta_point(ipoin),yigll,heta)
        call lagrange_only_interpol(NGLLZ,gamma_point(ipoin),zigll,hgamma)

        ! local point counter
        ipoin_local = ipoin_local + 1

        ! stores local interpolators
        hxi_point_store(:,ipoin_local) = hxi(:)
        heta_point_store(:,ipoin_local) = heta(:)
        hgamma_point_store(:,ipoin_local) = hgamma(:)

        ! stores normal
        normal_points(1,ipoin_local) = coupling_points(4,ipoin)  ! nx
        normal_points(2,ipoin_local) = coupling_points(5,ipoin)  ! ny
        normal_points(3,ipoin_local) = coupling_points(6,ipoin)  ! nz
      endif
    enddo

    ! vtk output of points found
    if (SAVE_MESH_FILES) then
      ! saves points coupling points
      allocate(iglob_tmp(npoints_local), &
               xstore_tmp(npoints_local), &
               ystore_tmp(npoints_local), &
               zstore_tmp(npoints_local),stat=ier)
      if (ier /= 0) stop 'error allocating array iglob_tmp arrays'
      iglob_tmp(:) = 0
      xstore_tmp(:) = 0.0_CUSTOM_REAL; ystore_tmp(:) = 0.0_CUSTOM_REAL; zstore_tmp(:) = 0.0_CUSTOM_REAL

      ipoin_local = 0
      do ipoin = 1,npoints_total
        if (islice_selected_point(ipoin) == myrank) then
          ! local point counter
          ipoin_local = ipoin_local + 1

          iglob_tmp(ipoin_local) = ipoin_local
          xstore_tmp(ipoin_local) = x_found(ipoin)
          ystore_tmp(ipoin_local) = y_found(ipoin)
          zstore_tmp(ipoin_local) = z_found(ipoin)
        endif
      enddo

      ! creates name for each rank process
      ! format: {TRACTION_PATH} / proc***_
      call create_name_database(prname_trac,myrank,TRACTION_PATH)

      filename = trim(prname_trac) // 'specfem_coupling_points_found'
      call write_VTK_data_points(npoints_local,xstore_tmp,ystore_tmp,zstore_tmp, &
                                 iglob_tmp,npoints_local,filename)
      deallocate(iglob_tmp,xstore_tmp,ystore_tmp,zstore_tmp)
    endif
  endif

contains

    subroutine lagrange_only_interpol(NGLL,xi,xigll,h)

    ! subroutine to compute only the Lagrange interpolants (h) at any point xi in [-1,1]
    ! (mixed CUSTOM_REAL and dp inputs)

    use constants, only: CUSTOM_REAL

    implicit none

    integer,intent(in) :: NGLL
    double precision,intent(in) :: xi
    double precision,dimension(NGLL),intent(in) :: xigll

    real(kind=CUSTOM_REAL),dimension(NGLL),intent(out) :: h

    ! local parameters
    integer :: dgr,i
    double precision :: prod1,prod2
    double precision :: x0,x

    do dgr = 1,NGLL

      prod1 = 1.d0
      prod2 = 1.d0

      ! lagrangian interpolants
      x0 = xigll(dgr)

      do i = 1,NGLL
        if (i /= dgr) then
          x = xigll(i)
          prod1 = prod1 * (xi - x)
          prod2 = prod2 * (x0 - x)
        endif
      enddo

      ! convert to custom real
      h(dgr) = real(prod1 / prod2,kind=CUSTOM_REAL)
    enddo

    end subroutine lagrange_only_interpol

  end subroutine determine_local_coupling_points

!
!-------------------------------------------------------------------------------------------------
!

  subroutine store_coupling_points_wavefield()

! stores velocity and traction on coupling points

  use specfem_par
  use specfem_par_coupling

  use specfem_par_acoustic, only: ispec_is_acoustic,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic
  use specfem_par_elastic, only: ispec_is_elastic,veloc
  use specfem_par_poroelastic, only: ispec_is_poroelastic,velocs_poroelastic

  implicit none

  ! local parameters
  integer :: ipoin,ipoin_local
  integer :: ispec

  real(kind=CUSTOM_REAL) :: vxd,vyd,vzd
  real(kind=CUSTOM_REAL) :: txd,tyd,tzd
  real(kind=CUSTOM_REAL) :: nx,ny,nz

  ! Lagrange interpolators
  real(kind=CUSTOM_REAL) :: hxi(NGLLX),heta(NGLLY),hgamma(NGLLZ)

  ! for collecting local wavefields on main
  integer :: npoints_local_proc,iproc

  ! checks subsampling recurrence
  if (mod(it-1,NTSTEP_BETWEEN_OUTPUT_SAMPLE) /= 0) return

  ! gets wavefield
  if (npoints_local > 0) then
    ! GPU simulations
    if (GPU_MODE) then
      ! not implemented yet
      call exit_MPI(myrank,'computing tractions is not supported yet for GPU simulations')

      !TODO: need to transfer velocity & stresses
      ! we need to transfer fields for computing velocity & traction in this routine
      ! poroelastic not supported yet on GPU
      if (ELASTIC_SIMULATION) then
        call transfer_veloc_from_device(NDIM*NGLOB_AB,veloc,Mesh_pointer)
      endif
      if (ACOUSTIC_SIMULATION) then
        call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                            potential_dot_acoustic, potential_dot_dot_acoustic, &
                                            Mesh_pointer)
      endif
    endif

    ! computes wavefield at point locations
    ipoin_local = 0

    do ipoin = 1,npoints_total
      ! slice with point stores current velocity and traction
      if (islice_selected_point(ipoin) == myrank) then
        ! local point counter
        ipoin_local = ipoin_local + 1

        ! gets local interpolators
        ! (1-D Lagrange interpolators)
        hxi(:) = hxi_point_store(:,ipoin_local)
        heta(:) = heta_point_store(:,ipoin_local)
        hgamma(:) = hgamma_point_store(:,ipoin_local)

        ! element
        ispec = ispec_selected_point(ipoin)

        ! velocity
        call compute_point_velocity(ispec,hxi,heta,hgamma,vxd,vyd,vzd)

        ! store velocity for local points
        veloc_points(1,ipoin_local) = vxd
        veloc_points(2,ipoin_local) = vyd
        veloc_points(3,ipoin_local) = vzd

        ! get normal for point
        nx = normal_points(1,ipoin_local)
        ny = normal_points(2,ipoin_local)
        nz = normal_points(3,ipoin_local)

        ! traction
        call compute_point_traction(ispec,nx,ny,nz,hxi,heta,hgamma,txd,tyd,tzd)

        ! store traction for local points
        traction_points(1,ipoin_local) = txd
        traction_points(2,ipoin_local) = tyd
        traction_points(3,ipoin_local) = tzd
      endif  ! islice
    enddo
  endif   ! npoints_local

  ! main process collect all points
  if (myrank == 0) then
    ! fill local point wavefield into total array for main process
    if (npoints_local > 0) then
      ipoin_local = 0
      do ipoin = 1,npoints_total
        if (islice_selected_point(ipoin) == 0) then
          ipoin_local = ipoin_local + 1
          ! store in full arrays
          veloc_total(:,ipoin) = veloc_points(:,ipoin_local)
          traction_total(:,ipoin) = traction_points(:,ipoin_local)
        endif
      enddo
    endif

    ! collect wavefield from all other processes
    do iproc = 1,NPROC-1
      ! receive local wavefield
      npoints_local_proc = nb_points_local_per_proc(iproc)
      if (npoints_local_proc > 0) then
        ! velocity
        ! receive array
        call recvv_cr(buffer_points_recv,NDIM*npoints_local_proc,iproc,itag)
        ! fill into total array
        ipoin_local = 0
        do ipoin = 1,npoints_total
          if (islice_selected_point(ipoin) == iproc) then
            ipoin_local = ipoin_local + 1
            ! store in full array
            veloc_total(:,ipoin) = buffer_points_recv(:,ipoin_local)
          endif
        enddo

        ! traction
        ! receive array
        call recvv_cr(buffer_points_recv,NDIM*npoints_local_proc,iproc,itag)
        ! fill into total array
        ipoin_local = 0
        do ipoin = 1,npoints_total
          if (islice_selected_point(ipoin) == iproc) then
            ipoin_local = ipoin_local + 1
            ! store in full array
            traction_total(:,ipoin) = buffer_points_recv(:,ipoin_local)
          endif
        enddo
      endif
    enddo
  else
    ! secondary processes send arrays
    if (npoints_local > 0) then
      call sendv_cr(veloc_points,NDIM*npoints_local,0,itag)
      call sendv_cr(traction_points,NDIM*npoints_local,0,itag)
    endif
  endif

  ! main process writes out velocity and traction
  if (myrank == 0) then
    write(IOUT_COUP) veloc_total,traction_total
  endif

contains

    subroutine compute_point_velocity(ispec,hxi,heta,hgamma,vxd,vyd,vzd)

    implicit none
    integer, intent(in) :: ispec
    real(kind=CUSTOM_REAL), intent(in) :: hxi(NGLLX),heta(NGLLY),hgamma(NGLLZ)
    real(kind=CUSTOM_REAL), intent(out) :: vxd,vyd,vzd

    ! local parameters
    real(kind=CUSTOM_REAL) :: hlagrange
    real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: veloc_element

    ! local parameters
    integer :: i,j,k,iglob

    ! initializes
    vxd = 0.0_CUSTOM_REAL
    vyd = 0.0_CUSTOM_REAL
    vzd = 0.0_CUSTOM_REAL

    ! elastic domain
    if (ispec_is_elastic(ispec)) then
      ! interpolates wavefield at exact location
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            hlagrange = hxi(i) * heta(j) * hgamma(k)
            ! velocity
            vxd = vxd + veloc(1,iglob) * hlagrange
            vyd = vyd + veloc(2,iglob) * hlagrange
            vzd = vzd + veloc(3,iglob) * hlagrange
          enddo
        enddo
      enddo
    endif

    ! acoustic domain
    if (ispec_is_acoustic(ispec)) then
      ! velocity vector
      call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,veloc_element)
      ! interpolates seismograms at exact receiver locations
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            hlagrange = hxi(i) * heta(j) * hgamma(k)
            ! velocity
            vxd = vxd + hlagrange * veloc_element(1,i,j,k)
            vyd = vyd + hlagrange * veloc_element(2,i,j,k)
            vzd = vzd + hlagrange * veloc_element(3,i,j,k)
          enddo
        enddo
      enddo
    endif

    ! poroelastic domain
    if (ispec_is_poroelastic(ispec)) then
      ! interpolates wavefield at exact location
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            hlagrange = hxi(i) * heta(j) * hgamma(k)
            ! velocity
            vxd = vxd + velocs_poroelastic(1,iglob) * hlagrange
            vyd = vyd + velocs_poroelastic(2,iglob) * hlagrange
            vzd = vzd + velocs_poroelastic(3,iglob) * hlagrange
          enddo
        enddo
      enddo
    endif

    end subroutine compute_point_velocity

    !-----------------------------------------------------------------------------------------------------

    subroutine compute_point_traction(ispec,nx,ny,nz,hxi,heta,hgamma,txd,tyd,tzd)

    use specfem_par_movie, only: stress_xx,stress_yy,stress_zz,stress_xy,stress_xz,stress_yz

    implicit none
    integer, intent(in) :: ispec
    real(kind=CUSTOM_REAL), intent(in) :: nx,ny,nz
    real(kind=CUSTOM_REAL), intent(in) :: hxi(NGLLX),heta(NGLLY),hgamma(NGLLZ)
    real(kind=CUSTOM_REAL), intent(out) :: txd,tyd,tzd

    ! local parameters
    real(kind=CUSTOM_REAL) :: hlagrange
    !real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: stress_element
    real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

    ! local parameters
    integer :: i,j,k

    ! initializes
    txd = 0.0_CUSTOM_REAL
    tyd = 0.0_CUSTOM_REAL
    tzd = 0.0_CUSTOM_REAL

    ! elastic domain
    if (ispec_is_elastic(ispec)) then
      ! interpolates wavefield at exact location
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            hlagrange = hxi(i) * heta(j) * hgamma(k)

            ! stress
            sigma_xx = stress_xx(i,j,k,ispec)
            sigma_yy = stress_yy(i,j,k,ispec)
            sigma_zz = stress_zz(i,j,k,ispec)
            sigma_xy = stress_xy(i,j,k,ispec)
            sigma_xz = stress_xz(i,j,k,ispec)
            sigma_yz = stress_yz(i,j,k,ispec)

            ! assumes symmetric stress tensor
            sigma_yx = sigma_xy
            sigma_zx = sigma_xz
            sigma_zy = sigma_yz

            ! traction
            txd = txd + (sigma_xx * nx + sigma_xy * ny + sigma_xz * nz) * hlagrange
            tyd = tyd + (sigma_yx * nx + sigma_yy * ny + sigma_yz * nz) * hlagrange
            tzd = tzd + (sigma_zx * nx + sigma_zy * ny + sigma_zz * nz) * hlagrange
          enddo
        enddo
      enddo
    endif

    ! acoustic domain
    if (ispec_is_acoustic(ispec)) then
      ! not implemented yet
      call exit_MPI(myrank,'computing traction is not supported yet for acoustic elements')
    endif

    ! poroelastic domain
    if (ispec_is_poroelastic(ispec)) then
      ! not implemented yet
      call exit_MPI(myrank,'computing traction is not supported yet for poroelastic elements')
    endif

    end subroutine compute_point_traction

  end subroutine store_coupling_points_wavefield

!
!-------------------------------------------------------------------------------------------------
!

! routines for running solver simulation with wavefield injection

  subroutine couple_with_injection_prepare_specfem_files()

  use constants, only: myrank,NDIM,MAX_STRING_LEN,NGLLSQUARE,IMAIN,IOUT,IIN,IIN_veloc_dsm

  use shared_parameters, only: DT,NSTEP,NPROC,TRACTION_PATH

  use specfem_par, only: t0,num_abs_boundary_faces

  use specfem_par_coupling

  implicit none

  ! local parameters
  integer :: ipoin,ipoin_local,i,ier
  character(len=MAX_STRING_LEN) :: prname_trac,filename

  ! file reading
  double precision :: start_time,dt_incr
  integer :: num_coupling_points_total,ntimesteps
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tmp_veloc,tmp_traction
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: tmp_veloc_points_timeseries,tmp_traction_points_timeseries
  integer :: offset_proc,reclen

  double precision :: tmin_file,tmax_file,tmin_curr,tmax_curr
  double precision :: sizeval

  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  ! synchronize first
  call synchronize_all()

  ! get MPI starting time
  tstart = wtime()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  using SPECFEM coupling'
    call flush_IMAIN()
  endif

  ! SPECFEM coupling
  !
  ! note: this setup routine will be called for the small-scale "local" simulation in the solver.
  !       at this point, we should have the wavefield solution from the coarse simulation saved in TRACTION_PATH folder,
  !       as a single file `specfem_coupling_solution.bin`.
  !
  !       the solution must still be interpolated for the time step size and duration of this small-scale simulation,
  !       similar as it is done for AxiSEM coupling with using the external tool `xreformat` in the modified AxiSEM folder.
  !
  !       for this coupling with injection implementation, the coupling boundary points are defined by the Stacey absorbing
  !       boundary faces, i.e., by num_abs_boundary_faces and corresponding surface arrays.
  !       we will assume that these boundary points are the ones used as coupling boundary points for which the solution
  !       was computed. That is, we won't re-locate the points in this solver simulation in case the small-scale simulation
  !       setup would have changed.
  !
  !       to inject the wavefield, we first need to interpolate in time the wavefield solution and store them
  !       for each process as local arrays into new files. this needs only to be done for slices that have boundary points.

  ! assumed number of local boundary points
  npoints_local = num_abs_boundary_faces * NGLLSQUARE

  ! get coupling points from each process
  allocate(nb_points_local_per_proc(0:NPROC-1),stat=ier)
  if (ier /= 0) stop 'Error allocating nb_points_local_per_proc array'
  nb_points_local_per_proc(:) = 0

  call gather_all_all_singlei(npoints_local,nb_points_local_per_proc,NPROC)

  ! assumed total number of boundary points
  npoints_total = sum(nb_points_local_per_proc(:))

  ! offset in total points array for this slice
  if (myrank == 0) then
    offset_proc = 0
  else
    offset_proc = sum(nb_points_local_per_proc(0:myrank-1))
  endif

  ! read full wavefield solution
  filename = trim(TRACTION_PATH) // '/' // 'specfem_coupling_solution.bin'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading solution file : ',trim(filename)
    call flush_IMAIN()
  endif

  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open file ',trim(filename)
    print *,'       Please check if file exists and was computed with a previous (coarse) regional forward simulation...'
    stop 'Error opening specfem_coupling_solution.bin'
  endif

  ! first line: #num_points #step_size #start_time #ntimesteps
  read(IIN,iostat=ier) num_coupling_points_total, dt_incr, start_time, ntimesteps

  ! check
  if (ier /= 0) then
    print *,'Error: rank ',myrank,' failed to read file ',trim(filename)
    print *,'       Please check that the file contains data...'
    stop 'Error reading specfem_coupling_solution.bin'
  endif

  ! solution file time range
  tmin_file = start_time  ! start_time is negative
  tmax_file = dble(ntimesteps-1) * dt_incr + start_time

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    total number of coupling points         : ',num_coupling_points_total
    write(IMAIN,*) '    number of time steps                    : ',ntimesteps
    write(IMAIN,*) '    time step size                          : ',dt_incr,'s'
    write(IMAIN,*) '    start time                              : ',tmin_file,'s'
    write(IMAIN,*) '    end time                                : ',tmax_file,'s'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check number of points (in case setup would have changed)
  if (num_coupling_points_total /= npoints_total) then
    print *,'Error: rank ',myrank,' found an invalid number of total boundary points ',num_coupling_points_total
    print *,'       total number of boundary points for this simulation should be ',npoints_total
    print *,'       Please check if solution file was computed for this small-scale simulation.'
    call exit_MPI(myrank,'Invalid number of boundary points found')
  endif

  ! safety check
  if (npoints_total == 0) then
    print *,'Error: no coupling boundary points defined for this simulation!'
    print *,'       Please check if mesh has absorbing boundaries setup properly to use COUPLE_WITH_INJECTION_TECHNIQUE'
    call exit_MPI(myrank,'No coupling boundary points found')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  interpolating time series:'
    write(IMAIN,*) '    maximum number of local coupling points : ',maxval(nb_points_local_per_proc)
    write(IMAIN,*)
    ! estimates array memory tmp_veloc_points_timeseries & tmp_traction_points_timeseries
    sizeval = 2.d0 * dble(NDIM) * dble(maxval(nb_points_local_per_proc)) * dble(ntimesteps) * dble(CUSTOM_REAL)
    ! adds memory for read arrays tmp_veloc & tmp_traction
    sizeval = sizeval + 2.d0 * dble(NDIM) * dble(npoints_total) * dble(CUSTOM_REAL)
    ! info output
    write(IMAIN,*) '    required maximum array size per slice   = ',sngl(sizeval / 1024.d0 / 1024.d0),'MB'
    write(IMAIN,*) '                                            = ',sngl(sizeval / 1024.d0 / 1024.d0 / 1024.d0),'GB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate temporary read arrays
  ! all points
  allocate(tmp_veloc(NDIM,npoints_total), &
           tmp_traction(NDIM,npoints_total),stat=ier)
  if (ier /= 0) stop 'Error allocating tmp_veloc,.. arrays'
  tmp_veloc(:,:) = 0.0_CUSTOM_REAL
  tmp_traction(:,:) = 0.0_CUSTOM_REAL

  ! local points, all times
  allocate(tmp_veloc_points_timeseries(NDIM,ntimesteps,npoints_local), &
           tmp_traction_points_timeseries(NDIM,ntimesteps,npoints_local),stat=ier)
  if (ier /= 0) stop 'Error allocating tmp_veloc_points_timeseries,.. arrays'
  tmp_veloc_points_timeseries(:,:,:) = 0.0_CUSTOM_REAL
  tmp_traction_points_timeseries(:,:,:) = 0.0_CUSTOM_REAL

  ! read in points and save local time series
  do i = 1,ntimesteps
    ! reads current time step values for all points
    read(IIN,iostat=ier) tmp_veloc,tmp_traction

    ! check
    if (ier /= 0) then
      print *,'Error: rank ',myrank,' failed to read wavefields in file ',trim(filename),' at timestep ',i
      print *,'       Please check that the file contains all wavefield data...'
      stop 'Error reading specfem_coupling_solution.bin'
    endif

    ! assigns values on local points
    if (npoints_local > 0) then
      ipoin_local = 0
      do ipoin = 1,npoints_total
        ! takes range between [offset,offset + npoints_local]
        if (ipoin > offset_proc .and. ipoin <= offset_proc + npoints_local) then
          ipoin_local = ipoin_local + 1
          tmp_veloc_points_timeseries(:,i,ipoin_local) = tmp_veloc(:,ipoin)
          tmp_traction_points_timeseries(:,i,ipoin_local) = tmp_traction(:,ipoin)
        endif
      enddo
      ! check
      if (ipoin_local /= npoints_local) call exit_MPI(myrank,'Invalid number of local points found')
    endif
  enddo

  close(IIN)

  ! free temporary arrays
  deallocate(tmp_veloc,tmp_traction)

  ! current simulation time range
  tmin_curr = - t0     ! t0 is positive for negative start times
  tmax_curr = dble((NSTEP-1)) * DT - t0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    current simulation time steps           : ',NSTEP
    write(IMAIN,*) '                       time step size       : ',DT,'s'
    write(IMAIN,*) '                       start time           : ',tmin_curr,'s'
    write(IMAIN,*) '                       end time             : ',tmax_curr,'s'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check time range of current simulation versus time span covered by solution file
  if (tmin_curr >= tmax_file .or. tmax_curr <= tmin_file) then
    print *,'Error: Invalid time range of current simulation: time is outside of time range of solution file'
    print *,'       time range of solution file (tmin/tmax)   = ',tmin_file,'/',tmax_file
    print *,'       current simulation time range (tmin/tmax) = ',tmin_curr,'/',tmax_curr
    print *,'       Please check simulation setup and/or use INJECTION_START_TIME in Par_file to specify, now exiting...'
    call exit_MPI(myrank,'Invalid time range for coupled injection simulation')
  endif

  if (tmin_curr < tmin_file) then
    ! user warning
    if (myrank == 0) then
      write(IMAIN,*) 'Warning: current simulation starts before start time of solution file!'
      write(IMAIN,*) '         time range of solution file (tmin/tmax)   = ',tmin_file,'/',tmax_file
      write(IMAIN,*) '         current simulation time range (tmin/tmax) = ',tmin_curr,'/',tmax_curr
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! synchronizes
  call synchronize_all()

  ! creates name for each rank process
  ! format: {TRACTION_PATH} / proc***_
  call create_name_database(prname_trac,myrank,TRACTION_PATH)

  filename = trim(prname_trac) // 'sol_specfem.bin'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    preparing tractions files...'
    write(IMAIN,*) '    for slice 0 file is ',trim(filename)
    write(IMAIN,*)
    ! filesize estimation
    sizeval = 2.d0 * dble(NDIM) * dble(maxval(nb_points_local_per_proc)) * dble(NSTEP) * dble(CUSTOM_REAL)
    write(IMAIN,*) '    estimated maximum file size per proc    = ',sngl(sizeval / 1024.d0 / 1024.d0),'MB'
    write(IMAIN,*) '                                            = ',sngl(sizeval / 1024.d0 / 1024.d0 / 1024.d0),'GB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! boundary arrays for reading/writing
  allocate(Veloc_specfem(NDIM,npoints_local),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'error allocating array 2192')
  Veloc_specfem(:,:) = 0.0_CUSTOM_REAL

  allocate(Tract_specfem(NDIM,npoints_local),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'error allocating array 2193')
  Tract_specfem(:,:) = 0.0_CUSTOM_REAL

  ! synchronizes to be sure all processes had enough memory
  call synchronize_all()

  ! file with direct access
  ! gets size of single record
  inquire(iolength=reclen) Veloc_specfem
  ! check integer size limit: size of reclen must fit onto an 4-byte integer
  if (reclen > int(2147483646.0 / 2)) then
    print *,'reclen needed exceeds integer 4-byte limit: ',reclen
    print *,'  ',reclen,' CUSTOM_REAL/ndim/npoints_local = ',CUSTOM_REAL, NDIM, npoints_local
    print *,'bit size Fortran: ',bit_size(reclen)
    call exit_MPI(myrank,"Error reclen integer limit")
  endif
  ! record size for velocity & traction fields
  reclen = 2 * reclen

  ! open files to store wavefield solution for this slice and simulation time stepping
  !open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  ! w/ direct access
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted', &
       access='direct',recl=reclen,iostat=ier)

  ! check
  if (ier /= 0) then
    print *,'Error: could not open file ',trim(filename)
    print *,'       Please check if path exists...'
    stop 'Error opening database proc***_sol_specfem.bin'
  endif

  ! interpolate wavefield solution in time
  !
  ! Sinc (Whittaker-Shannon) interpolation:
  !   continuous over whole range, most accurate, also very accurate in frequency domain
  !
  ! Cubic splines interpolation:
  !   continuous derivatives over whole range
  !
  ! Catmull-Rom (cubic):
  !   cubic interpolation, local on neighbors, smooth on equidistant samples and exact at nodal points
  !
  ! cubic local:
  !   cubic interpolation, local on neighbors
  !
  ! time series interpolation - timing benchmark:
  !   given input: total number of coupling points (npoints_total)         = 68800
  !                number of time steps (ntimesteps)                       = 2500
  !
  !   interpolating time series for file output:
  !                maximum number of local coupling points (npoints_local) = 17200
  !                current simulation time steps (NSTEP)                   = 10000
  !
  !
  ! -> interpolation with sinc (Whittaker–Shannon) takes around ~ 601.24 s
  !                       cubic splines                         ~ 586.76 s
  !                       Catmull-Rom                           ~  22.68 s
  !                       cubic local                           ~  17.71 s
  !
  ! as a conclusion, unless we need a very high accurary in frequency-domain (sinc interpolation),
  ! the Catmull-Rom interpolation provides the best trade-off between accuracy and speed for time series interpolation.
  !
  ! used for testing...
  ! sinc (Whittaker–Shannon) interpolation in time
  !call interpolate_with_sinc(ntimesteps,dt_incr,start_time, &
  !                           tmp_veloc_points_timeseries,tmp_traction_points_timeseries)
  !
  ! cubic splines
  !call interpolate_with_cubic_splines(ntimesteps,dt_incr,start_time, &
  !                                    tmp_veloc_points_timeseries,tmp_traction_points_timeseries)
  !
  ! cubic
  !call interpolate_with_cubic_local(ntimesteps,dt_incr,start_time, &
  !                                  tmp_veloc_points_timeseries,tmp_traction_points_timeseries)
  !

  ! Catmull-Rom interpolation (preferred)
  call interpolate_with_Catmull_Rom(ntimesteps,dt_incr,start_time, &
                                    tmp_veloc_points_timeseries,tmp_traction_points_timeseries)

  ! closes file
  close(IOUT)

  ! free temporary arrays
  deallocate(tmp_veloc_points_timeseries,tmp_traction_points_timeseries)

  ! synchronizes MPI processes
  call synchronize_all()

  ! clear arrays for reading
  Veloc_specfem(:,:) = 0.0_CUSTOM_REAL
  Tract_specfem(:,:) = 0.0_CUSTOM_REAL

  ! open the file created before
  ! to read in corresponding time slices in routine compute_Stacey_elastic()
  !open(unit=IIN_veloc_dsm,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  !
  ! with direct access
  open(unit=IIN_veloc_dsm,file=trim(filename),status='old',action='read',form='unformatted', &
       access='direct',recl=reclen,iostat=ier)

  ! check
  if (ier /= 0) then
    print *,'Error: could not open file ',trim(filename)
    print *,'Please check if traction file has been created for coupling with SPECFEM solution...'
    stop 'Error opening tractions file proc****_sol_specfem.bin'
  endif

  ! synchronizes MPI processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - tstart
    write(IMAIN,*)
    write(IMAIN,*) '  Elapsed time for preparing coupling files in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine couple_with_injection_prepare_specfem_files

!
!-------------------------------------------------------------------------------------------------
!

  subroutine interpolate_with_Catmull_Rom(ntimesteps,dt_incr,start_time, &
                                            tmp_veloc_points_timeseries,tmp_traction_points_timeseries)

! interpolates velocity & traction for current simulation times using Catmull-Rom (cubic) interpolation

  use constants, only: myrank,CUSTOM_REAL,NDIM,IMAIN,IOUT,OUTPUT_FILES

  use shared_parameters, only: DT,NSTEP

  use specfem_par, only: t0

  use specfem_par_coupling, only: npoints_local,Veloc_specfem,Tract_specfem

  implicit none

  integer, intent(in) :: ntimesteps
  double precision, intent(in) :: dt_incr,start_time
  real(kind=CUSTOM_REAL), dimension(NDIM,ntimesteps,npoints_local), intent(in) :: tmp_veloc_points_timeseries, &
                                                                                  tmp_traction_points_timeseries

  ! local parameters
  integer :: it_tmp,ipoin_local,i,ier
  double precision :: current_time
  double precision :: tmin_file,tmax_file

  ! interpolation
  real(kind=CUSTOM_REAL),dimension(NDIM) :: vel_p0,vel_p1,vel_p2,vel_p3
  real(kind=CUSTOM_REAL),dimension(NDIM) :: vel_a,vel_b,vel_c,vel_d
  real(kind=CUSTOM_REAL),dimension(NDIM) :: tract_p0,tract_p1,tract_p2,tract_p3
  real(kind=CUSTOM_REAL),dimension(NDIM) :: tract_a,tract_b,tract_c,tract_d

  real(kind=CUSTOM_REAL) :: xi
  integer :: idx,idx0,idx1,idx2,idx3

  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  ! debugging - turn on to check accuracy of interpolated time serie
  logical, parameter :: DEBUG = .false.     ! only works if myrank == 0 has local points to interpolate

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    using Catmull-Rom (cubic) interpolation'
    call flush_IMAIN()
  endif

  ! checks if anything to do in this slice
  if (npoints_local == 0) return

  ! debug file output
  if (DEBUG .and. myrank == 0) then
    ! time series from solution file
    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel.txt',status='unknown',action='write')
    do i = 1,ntimesteps
      ! format: # time # vel_Z
      write(87,*) dble((i-1)) * dt_incr + start_time, tmp_veloc_points_timeseries(3,i,1000)  ! at local point 1000
    enddo
    close(87)
    ! interpolated time series
    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel_interpol.txt',status='unknown',action='write')
  endif

  ! get MPI starting time
  tstart = wtime()

  ! solution file time range
  tmin_file = start_time  ! start_time is negative
  tmax_file = dble(ntimesteps-1) * dt_incr + start_time

  ! loop over time steps of current "small-scale" simulation
  do it_tmp = 1,NSTEP
    ! current simulation time
    current_time = dble((it_tmp-1)) * DT - t0

    ! checks if current_time is out of bounds
    if (current_time < tmin_file .or. current_time > tmax_file) then
      ! time outside of provided time range
      ! set to zero
      Veloc_specfem(:,:) = 0.0_CUSTOM_REAL
      Tract_specfem(:,:) = 0.0_CUSTOM_REAL
    else
      ! Catmull-Rom (cubic) interpolation:
      ! needs 4 points p0,p1,p2,p3 and interpolates point p(t) between point p1 and p2
      ! see: https://www.iquilezles.org/www/articles/minispline/minispline.htm

      ! determine time interval index
      idx = int( (current_time - tmin_file) / dt_incr) + 1

      ! checks bounds
      if (idx < 1) idx = 1
      if (idx > ntimesteps) idx = ntimesteps

      ! relative offset between [0,1] in interval [idx-1,idx] of time series
      xi = real( ((current_time - tmin_file) - (idx-1) * dt_incr) / dt_incr, kind=CUSTOM_REAL)

      ! interpolation points
      ! at the boundaries idx == 1, idx < ntimesteps-1, we take the limited values and mimick interpolation in these cases.
      if (idx == 1) then
        idx0 = idx    ! no previous one, use same as idx1
        idx1 = idx
        idx2 = idx+1
        idx3 = idx+2
      else if (idx == ntimesteps - 1) then
        idx0 = idx-1
        idx1 = idx
        idx2 = idx+1
        idx3 = idx+1  ! no later one, use same as idx2
      else if (idx == ntimesteps) then
        idx0 = idx-1
        idx1 = idx
        idx2 = idx    ! no later ones
        idx3 = idx
      else
        ! it > 1 .and. it < ntimesteps-1
        idx0 = idx-1
        idx1 = idx
        idx2 = idx+1
        idx3 = idx+2
      endif

      do ipoin_local = 1,npoints_local
        ! velocity
        ! interpolation points
        vel_p0(:) = tmp_veloc_points_timeseries(:,idx0,ipoin_local)
        vel_p1(:) = tmp_veloc_points_timeseries(:,idx1,ipoin_local)
        vel_p2(:) = tmp_veloc_points_timeseries(:,idx2,ipoin_local)
        vel_p3(:) = tmp_veloc_points_timeseries(:,idx3,ipoin_local)

        ! Catmull-Rom interpolation
        !   a = 2 * p1
        !   b = p2 - p0
        !   c = 2 * p0 - 5 * p1 + 4 * p2 - p3
        !   d = -p0 + 3 * p1 - 3 * p2 + p3
        !   and interpolate value at xi
        !   val = 1/2 * ( a + b * xi + c * xi^2 + d * xi^3 )
        !
        vel_a(:) = 2.0_CUSTOM_REAL * vel_p1(:)
        vel_b(:) = vel_p2(:) - vel_p0(:)
        vel_c(:) = 2.0_CUSTOM_REAL * vel_p0(:) - 5.0_CUSTOM_REAL * vel_p1(:) + 4.0_CUSTOM_REAL * vel_p2(:) - vel_p3(:)
        vel_d(:) = -vel_p0(:) + 3.0_CUSTOM_REAL * vel_p1(:) - 3.0_CUSTOM_REAL * vel_p2(:) + vel_p3(:)

        ! cubic polynomial: a + b * t + c * t^2 + d * t^3 = a + t * (b + t * (c + t * d))
        Veloc_specfem(:,ipoin_local) = 0.5_CUSTOM_REAL * ( vel_a(:) + xi * (vel_b(:) + xi * (vel_c(:) + xi * vel_d(:))))

        ! traction
        ! interpolation points
        tract_p0(:) = tmp_traction_points_timeseries(:,idx0,ipoin_local)
        tract_p1(:) = tmp_traction_points_timeseries(:,idx1,ipoin_local)
        tract_p2(:) = tmp_traction_points_timeseries(:,idx2,ipoin_local)
        tract_p3(:) = tmp_traction_points_timeseries(:,idx3,ipoin_local)

        ! Catmull-Rom interpolation
        tract_a(:) = 2.0_CUSTOM_REAL * tract_p1(:)
        tract_b(:) = tract_p2(:) - tract_p0(:)
        tract_c(:) = 2.0_CUSTOM_REAL * tract_p0(:) - 5.0_CUSTOM_REAL * tract_p1(:) + 4.0_CUSTOM_REAL * tract_p2(:) - tract_p3(:)
        tract_d(:) = -tract_p0(:) + 3.0_CUSTOM_REAL * tract_p1(:) - 3.0_CUSTOM_REAL * tract_p2(:) + tract_p3(:)

        ! cubic polynomial: a + b * t + c * t^2 + d * t^3 = a + t * (b + t * (c + t * d))
        Tract_specfem(:,ipoin_local) = 0.5_CUSTOM_REAL * ( tract_a(:) + xi * (tract_b(:) + xi * (tract_c(:) + xi * tract_d(:))))
      enddo
    endif

    ! store in file
    !write(IOUT,iostat=ier) Veloc_specfem,Tract_specfem
    ! w/ direct access
    write(IOUT,rec=it_tmp,iostat=ier) Veloc_specfem,Tract_specfem

    ! check
    if (ier /= 0) then
      print *,'Error: rank ',myrank,' failed to write out interpolated wavefields'
      print *,'       Please check if enough file space is available...'
      stop 'Error output interpolated wavefields'
    endif

    ! user output
    if (myrank == 0) then
      if (mod(it_tmp,NSTEP/10) == 0) then
        tCPU = wtime() - tstart
        write(IMAIN,*) '    interpolated time step ',it_tmp,'out of',NSTEP,' - elapsed time:',sngl(tCPU),'s'
        call flush_IMAIN()
      endif
    endif

    ! debug file output
    if (DEBUG .and. myrank == 0) then
      ! format: # time # vel_Z
      write(87,*) current_time, Veloc_specfem(3,1000)  ! at local point 1000
    endif
  enddo

  ! debug file output
  if (DEBUG .and. myrank == 0) close(87)

  end subroutine interpolate_with_Catmull_Rom

!
!-------------------------------------------------------------------------------------------------
!

! used for testing...

!  subroutine interpolate_with_sinc(ntimesteps,dt_incr,start_time, &
!                                   tmp_veloc_points_timeseries,tmp_traction_points_timeseries)
!
!! interpolates velocity & traction for current simulation times using Sinc (Whittaker–Shannon) interpolation
!
!  use constants, only: myrank,CUSTOM_REAL,NDIM,IMAIN,IOUT,OUTPUT_FILES
!
!  use shared_parameters, only: DT,NSTEP
!
!  use specfem_par, only: t0
!
!  use specfem_par_coupling, only: npoints_local,Veloc_specfem,Tract_specfem
!
!  implicit none
!
!  integer, intent(in) :: ntimesteps
!  double precision, intent(in) :: dt_incr,start_time
!  real(kind=CUSTOM_REAL), dimension(NDIM,ntimesteps,npoints_local), intent(in) :: tmp_veloc_points_timeseries, &
!                                                                                  tmp_traction_points_timeseries
!
!  ! local parameters
!  integer :: it_tmp,ipoin_local,i
!  double precision :: current_time
!  double precision :: tmin_file,tmax_file
!
!  real(kind=CUSTOM_REAL), dimension(ntimesteps) :: sinc_table
!  real(kind=CUSTOM_REAL) :: vel(NDIM),trac(NDIM)
!
!  ! timing
!  double precision :: tstart,tCPU
!  double precision, external :: wtime
!
!  ! debugging
!  logical, parameter :: DEBUG = .true.     ! only works if myrank == 0 has local points to interpolate
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) '    using sinc interpolation'
!    call flush_IMAIN()
!  endif
!
!  ! checks if anything to do in this slice
!  if (npoints_local == 0) return
!
!  ! debug file output
!  if (DEBUG .and. myrank == 0) then
!    ! time series from solution file
!    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel.txt',status='unknown',action='write')
!    do i = 1,ntimesteps
!      ! format: # time # vel_Z
!      write(87,*) dble((i-1)) * dt_incr + start_time, tmp_veloc_points_timeseries(3,i,1000)  ! at local point 1000
!    enddo
!    close(87)
!    ! interpolated time series
!    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel_interpol.txt',status='unknown',action='write')
!  endif
!
!  ! get MPI starting time
!  tstart = wtime()
!
!  ! solution file time range
!  tmin_file = start_time  ! start_time is negative
!  tmax_file = dble(ntimesteps-1) * dt_incr + start_time
!
!  ! loop over time steps of current "small-scale" simulation
!  do it_tmp = 1,NSTEP
!    ! current simulation time
!    current_time = dble((it_tmp-1)) * DT - t0
!
!    ! checks if current_time is out of bounds
!    if (current_time < tmin_file .or. current_time > tmax_file) then
!      ! time outside of provided time range
!      ! set to zero
!      Veloc_specfem(:,:) = 0.0_CUSTOM_REAL
!      Tract_specfem(:,:) = 0.0_CUSTOM_REAL
!    else
!      ! sinc interpolation
!      ! setup sinc table to determine values at current_time
!      call compute_sinc_table(sinc_table,current_time,ntimesteps,dt_incr,start_time)
!
!      ! loop over points in this slice
!      do ipoin_local = 1,npoints_local
!        ! interpolate over solution time span
!        vel(:) = 0.0_CUSTOM_REAL
!        trac(:) = 0.0_CUSTOM_REAL
!        do i = 1,ntimesteps
!          ! velocity
!          vel(:) = vel(:) + sinc_table(i) * tmp_veloc_points_timeseries(:,i,ipoin_local)
!          ! traction
!          trac(:) = trac(:) + sinc_table(i) * tmp_traction_points_timeseries(:,i,ipoin_local)
!        enddo
!
!        ! store interpolated arrays for the current time
!        Veloc_specfem(:,ipoin_local) = vel(:)
!        Tract_specfem(:,ipoin_local) = trac(:)
!      enddo
!    endif
!
!    ! store in file
!    write(IOUT) Veloc_specfem,Tract_specfem
!
!    ! user output
!    if (myrank == 0) then
!      if (mod(it_tmp,NSTEP/10) == 0) then
!        tCPU = wtime() - tstart
!        write(IMAIN,*) '    interpolated time step ',it_tmp,'out of',NSTEP,' - elapsed time:',sngl(tCPU),'s'
!        call flush_IMAIN()
!      endif
!    endif
!
!    ! debug file output
!    if (DEBUG .and. myrank == 0) then
!      ! format: # time # vel_Z
!      write(87,*) current_time, Veloc_specfem(3,1000)  ! at local point 1000
!    endif
!  enddo
!
!  ! debug file output
!  if (DEBUG .and. myrank == 0) close(87)
!
!  end subroutine interpolate_with_sinc
!
!!
!!-------------------------------------------------------------------------------------------------
!!
!
!  subroutine compute_sinc_table(sinc_tab,tcurrent,n,dt,t0)
!
!  ! computes table for sinc (Whittaker–Shannon) interpolation
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer, intent(in) :: n
!  real(kind=CUSTOM_REAL), intent(out) :: sinc_tab(n)
!  double precision, intent(in) :: tcurrent,dt,t0
!
!  ! local parameters
!  integer :: i
!  real(kind=CUSTOM_REAL) :: t_sinc, dt_inv
!
!  ! inverted time step (multiplication will be faster than division)
!  dt_inv = real(1.d0 / dt,kind=CUSTOM_REAL)
!
!  do i = 1,n
!    ! time diff argument: sinc( (t - nT)/T ) with T the time step of the given time series
!    t_sinc = real((tcurrent - (dble(i-1) * dt + t0) ),kind=CUSTOM_REAL) * dt_inv
!
!    ! store interpolation table
!    sinc_tab(i) = func_mysinc(t_sinc)
!  enddo
!
!contains
!
!    ! Sinc function
!    function func_mysinc(x) result (sincval)
!
!    use constants, only: CUSTOM_REAL
!
!    implicit none
!    real(kind=CUSTOM_REAL), intent(in) :: x
!    real(kind=CUSTOM_REAL) :: sincval
!
!    ! local parameters
!    real(kind=CUSTOM_REAL) :: xn
!    real(kind=CUSTOM_REAL), parameter :: PI = 3.141592653589793_CUSTOM_REAL
!
!    ! normalized sinc function
!    if (abs(x) > 1.e-6) then
!      ! using PI * x instead of just x as in the unnormalized sinc function
!      ! note: for x = 1.d-6 sin( PI*x ) / ( PI*x) ~ 0.999999999998
!      xn = PI * x
!      sincval = sin(xn) / xn
!    else
!      ! sinc function value at x == 0
!      sincval = 1.0_CUSTOM_REAL
!    endif
!
!    end function func_mysinc
!
!  end subroutine compute_sinc_table

!
!-------------------------------------------------------------------------------------------------
!

! used for testing...

!  subroutine interpolate_with_cubic_splines(ntimesteps,dt_incr,start_time, &
!                                            tmp_veloc_points_timeseries,tmp_traction_points_timeseries)
!
!! interpolates velocity & traction for current simulation times using cubic splines
!
!  use constants, only: myrank,CUSTOM_REAL,NDIM,IMAIN,IOUT,OUTPUT_FILES
!
!  use shared_parameters, only: DT,NSTEP
!
!  use specfem_par, only: t0
!
!  use specfem_par_coupling, only: npoints_local,Veloc_specfem,Tract_specfem
!
!  implicit none
!
!  integer, intent(in) :: ntimesteps
!  double precision, intent(in) :: dt_incr,start_time
!  real(kind=CUSTOM_REAL), dimension(NDIM,ntimesteps,npoints_local), intent(in) :: tmp_veloc_points_timeseries, &
!                                                                                  tmp_traction_points_timeseries
!
!  ! local parameters
!  integer :: it_tmp,ipoin_local,i,icomp,ier
!  double precision :: current_time
!  double precision :: tmin_file,tmax_file
!
!  ! cubic splines
!  ! note: coefficient a is equal to initial values from timeseries
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: vel_coeffs_b, vel_coeffs_c, vel_coeffs_d
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: tract_coeffs_b, tract_coeffs_c, tract_coeffs_d
!
!  real(kind=CUSTOM_REAL), dimension(ntimesteps) :: b, c, d, y
!  real(kind=CUSTOM_REAL) :: xi
!  integer :: idx
!
!  ! timing
!  double precision :: tstart,tCPU
!  double precision, external :: wtime
!
!  ! debugging
!  logical, parameter :: DEBUG = .true.     ! only works if myrank == 0 has local points to interpolate
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) '    using cubic spline interpolation'
!    call flush_IMAIN()
!  endif
!
!  ! checks if anything to do in this slice
!  if (npoints_local == 0) return
!
!  ! debug file output
!  if (DEBUG .and. myrank == 0) then
!    ! time series from solution file
!    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel.txt',status='unknown',action='write')
!    do i = 1,ntimesteps
!      ! format: # time # vel_Z
!      write(87,*) dble((i-1)) * dt_incr + start_time, tmp_veloc_points_timeseries(3,i,1000)  ! at local point 1000
!    enddo
!    close(87)
!    ! interpolated time series
!    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel_interpol.txt',status='unknown',action='write')
!  endif
!
!  ! get MPI starting time
!  tstart = wtime()
!
!  ! pre-compute spline coefficients
!  ! allocate arrays
!  allocate(vel_coeffs_b(NDIM,ntimesteps,npoints_local), &
!           vel_coeffs_c(NDIM,ntimesteps,npoints_local), &
!           vel_coeffs_d(NDIM,ntimesteps,npoints_local),stat=ier)
!  if (ier /= 0) stop 'Error allocating spline coefficients for velocity'
!  vel_coeffs_b(:,:,:) = 0.0_CUSTOM_REAL
!  vel_coeffs_c(:,:,:) = 0.0_CUSTOM_REAL
!  vel_coeffs_d(:,:,:) = 0.0_CUSTOM_REAL
!
!  allocate(tract_coeffs_b(NDIM,ntimesteps,npoints_local), &
!           tract_coeffs_c(NDIM,ntimesteps,npoints_local), &
!           tract_coeffs_d(NDIM,ntimesteps,npoints_local),stat=ier)
!  if (ier /= 0) stop 'Error allocating spline coefficients for traction'
!  tract_coeffs_b(:,:,:) = 0.0_CUSTOM_REAL
!  tract_coeffs_c(:,:,:) = 0.0_CUSTOM_REAL
!  tract_coeffs_d(:,:,:) = 0.0_CUSTOM_REAL
!
!  ! loop over points in this slice
!  do ipoin_local = 1,npoints_local
!    ! compute spline coefficients
!    ! velocity
!    do icomp = 1,NDIM
!      do i = 1,ntimesteps
!        y(i) = tmp_veloc_points_timeseries(icomp,i,ipoin_local)
!      enddo
!      ! spline coefficients
!      call compute_spline_timeseries(y,ntimesteps,dt_incr,b,c,d)
!      ! store coefficients for evaluation
!      do i = 1,ntimesteps
!        vel_coeffs_b(icomp,i,ipoin_local) = b(i)
!        vel_coeffs_c(icomp,i,ipoin_local) = c(i)
!        vel_coeffs_d(icomp,i,ipoin_local) = d(i)
!      enddo
!    enddo
!    ! traction
!    do icomp = 1,NDIM
!      do i = 1,ntimesteps
!        y(i) = tmp_traction_points_timeseries(icomp,i,ipoin_local)
!      enddo
!      ! spline coefficients
!      call compute_spline_timeseries(y,ntimesteps,dt_incr,b,c,d)
!      ! store coefficients for evaluation
!      do i = 1,ntimesteps
!        tract_coeffs_b(icomp,i,ipoin_local) = b(i)
!        tract_coeffs_c(icomp,i,ipoin_local) = c(i)
!        tract_coeffs_d(icomp,i,ipoin_local) = d(i)
!      enddo
!    enddo
!  enddo
!
!  ! user output
!  if (myrank == 0) then
!    tCPU = wtime() - tstart
!    write(IMAIN,*) '    setting up spline coefficients took ',sngl(tCPU),'s'
!    call flush_IMAIN()
!  endif
!
!  ! solution file time range
!  tmin_file = start_time  ! start_time is negative
!  tmax_file = dble(ntimesteps-1) * dt_incr + start_time
!
!  ! loop over time steps of current "small-scale" simulation
!  do it_tmp = 1,NSTEP
!    ! current simulation time
!    current_time = dble((it_tmp-1)) * DT - t0
!
!    ! checks if current_time is out of bounds
!    if (current_time < tmin_file .or. current_time > tmax_file) then
!      ! time outside of provided time range
!      ! set to zero
!      Veloc_specfem(:,:) = 0.0_CUSTOM_REAL
!      Tract_specfem(:,:) = 0.0_CUSTOM_REAL
!    else
!      ! cubic spline interpolation
!
!      ! determine time interval index
!      idx = int( (current_time - tmin_file) / dt_incr) + 1
!
!      ! offset in interval interval [tmin,tmax] of time series
!      xi = real( (current_time - tmin_file) - (idx-1) * dt_incr, kind=CUSTOM_REAL)
!
!      ! interpolate
!      ! velocity
!      do ipoin_local = 1,npoints_local
!        do icomp = 1,NDIM
!          Veloc_specfem(icomp,ipoin_local) = tmp_veloc_points_timeseries(icomp,idx,ipoin_local) + &
!                          xi * (vel_coeffs_b(icomp,idx,ipoin_local) + &
!                          xi * (vel_coeffs_c(icomp,idx,ipoin_local) + xi * vel_coeffs_d(icomp,idx,ipoin_local)))
!        enddo
!      enddo
!
!      ! traction
!      do ipoin_local = 1,npoints_local
!        do icomp = 1,NDIM
!          Tract_specfem(icomp,ipoin_local) = tmp_traction_points_timeseries(icomp,idx,ipoin_local) + &
!                          xi * (tract_coeffs_b(icomp,idx,ipoin_local) + &
!                          xi * (tract_coeffs_c(icomp,idx,ipoin_local) + xi * tract_coeffs_d(icomp,idx,ipoin_local)))
!        enddo
!      enddo
!    endif
!
!    ! store in file
!    write(IOUT) Veloc_specfem,Tract_specfem
!
!    ! user output
!    if (myrank == 0) then
!      if (mod(it_tmp,NSTEP/10) == 0) then
!        tCPU = wtime() - tstart
!        write(IMAIN,*) '    interpolated time step ',it_tmp,'out of',NSTEP,' - elapsed time:',sngl(tCPU),'s'
!        call flush_IMAIN()
!      endif
!    endif
!
!    ! debug file output
!    if (DEBUG .and. myrank == 0) then
!      ! format: # time # vel_Z
!      write(87,*) current_time, Veloc_specfem(3,1000)  ! at local point 1000
!    endif
!  enddo
!
!  ! debug file output
!  if (DEBUG .and. myrank == 0) close(87)
!
!contains
!
!    subroutine compute_spline_timeseries(y,n,dt,b,c,d)
!
!    use constants, only: CUSTOM_REAL
!
!    implicit none
!
!    integer, intent(in) :: n
!    real(kind=CUSTOM_REAL), dimension(n),intent(in) :: y
!    double precision, intent(in) :: dt
!    real(kind=CUSTOM_REAL), dimension(n),intent(out) :: b,c,d
!
!    ! local parameters
!    real(kind=CUSTOM_REAL), dimension(n) :: z
!    ! Precomputed constants
!    real(kind=CUSTOM_REAL) :: three_dt_inv,dt_inv
!    integer :: i
!
!    ! Precompute constants
!    dt_inv = 1.0_CUSTOM_REAL / real(dt_incr,kind=CUSTOM_REAL)
!    three_dt_inv = 3.0_CUSTOM_REAL * dt_inv
!
!    ! Step 1: Compute z values directly
!    z(1) = 0.0_CUSTOM_REAL  ! initial z array to zero
!    z(n) = 0.0_CUSTOM_REAL
!    do i = 2, n-1
!      z(i) = three_dt_inv * (y(i+1) - 2.0d0 * y(i) + y(i-1))
!    enddo
!
!    ! Step 2: Compute coefficients b, c, d directly
!
!    ! note: coefficient a(:) = y(:)
!
!    do i = 1, n-1
!      b(i) = dt_inv * (y(i+1) - y(i)) - dt * (2.0d0 * z(i) + z(i+1)) / 3.0d0
!      c(i) = z(i)
!      d(i) = (z(i+1) - z(i)) / (3.0d0 * dt)
!    enddo
!
!    ! final point not needed in evaluation as last point value will be == a(n) == y(n)
!    b(n) = 0.0_CUSTOM_REAL
!    c(n) = 0.0_CUSTOM_REAL
!    d(n) = 0.0_CUSTOM_REAL
!
!    end subroutine compute_spline_timeseries
!
!  end subroutine interpolate_with_cubic_splines

!
!-------------------------------------------------------------------------------------------------
!

! used for testing...
!
!  subroutine interpolate_with_cubic_local(ntimesteps,dt_incr,start_time, &
!                                            tmp_veloc_points_timeseries,tmp_traction_points_timeseries)
!
!! interpolates velocity & traction for current simulation times using cubic interpolation locally with 4 points
!
!  use constants, only: myrank,CUSTOM_REAL,NDIM,IMAIN,IOUT,OUTPUT_FILES
!
!  use shared_parameters, only: DT,NSTEP
!
!  use specfem_par, only: t0
!
!  use specfem_par_coupling, only: npoints_local,Veloc_specfem,Tract_specfem
!
!  implicit none
!
!  integer, intent(in) :: ntimesteps
!  double precision, intent(in) :: dt_incr,start_time
!  real(kind=CUSTOM_REAL), dimension(NDIM,ntimesteps,npoints_local), intent(in) :: tmp_veloc_points_timeseries, &
!                                                                                  tmp_traction_points_timeseries
!
!  ! local parameters
!  integer :: it_tmp,ipoin_local,i
!  double precision :: current_time
!  double precision :: tmin_file,tmax_file
!
!  ! interpolation
!  real(kind=CUSTOM_REAL),dimension(NDIM) :: p0,p1,p2,p3
!  real(kind=CUSTOM_REAL) :: cs1,cs2,cs3,cs4
!
!  real(kind=CUSTOM_REAL) :: xi
!  integer :: idx,idx0,idx1,idx2,idx3
!
!  ! timing
!  double precision :: tstart,tCPU
!  double precision, external :: wtime
!
!  ! debugging
!  logical, parameter :: DEBUG = .true.     ! only works if myrank == 0 has local points to interpolate
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) '    using local (cubic) interpolation'
!    call flush_IMAIN()
!  endif
!
!  ! checks if anything to do in this slice
!  if (npoints_local == 0) return
!
!  ! debug file output
!  if (DEBUG .and. myrank == 0) then
!    ! time series from solution file
!    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel.txt',status='unknown',action='write')
!    do i = 1,ntimesteps
!      ! format: # time # vel_Z
!      write(87,*) dble((i-1)) * dt_incr + start_time, tmp_veloc_points_timeseries(3,i,1000)  ! at local point 1000
!    enddo
!    close(87)
!    ! interpolated time series
!    open(87,file=trim(OUTPUT_FILES)//'/debug_timeseries_vel_interpol.txt',status='unknown',action='write')
!  endif
!
!  ! get MPI starting time
!  tstart = wtime()
!
!  ! solution file time range
!  tmin_file = start_time  ! start_time is negative
!  tmax_file = dble(ntimesteps-1) * dt_incr + start_time
!
!  ! loop over time steps of current "small-scale" simulation
!  do it_tmp = 1,NSTEP
!    ! current simulation time
!    current_time = dble((it_tmp-1)) * DT - t0
!
!    ! checks if current_time is out of bounds
!    if (current_time < tmin_file .or. current_time > tmax_file) then
!      ! time outside of provided time range
!      ! set to zero
!      Veloc_specfem(:,:) = 0.0_CUSTOM_REAL
!      Tract_specfem(:,:) = 0.0_CUSTOM_REAL
!    else
!      ! cubic interpolation
!
!      ! determine time interval index
!      idx = int( (current_time - tmin_file) / dt_incr) + 1
!
!      ! checks bounds
!      if (idx < 1) idx = 1
!      if (idx > ntimesteps) idx = ntimesteps
!
!      ! relative offset between [0,1] in interval [idx-1,idx] of time series
!      xi = real( ((current_time - tmin_file) - (idx-1) * dt_incr) / dt_incr, kind=CUSTOM_REAL)
!
!      ! Cubic (spline?) values
!      !
!      ! w == 0 -> c1 = 1/6, c2 = 2/3, c3 = 1/6, c4 = 0
!      ! w == 1 -> c1 = 0  , c2 = 1/6, c3 = 2/3, c4 = 1/6
!      cs4 = xi * xi * xi / 6.d0
!      cs1 = 1.d0/6.d0 + xi * (xi-1.d0) / 2.d0 - cs4
!      cs3 = xi + cs1 - 2.d0 * cs4
!      cs2 = 1.d0 - cs1 - cs3 - cs4
!
!      ! interpolation points
!      ! at the boundaries idx == 1, idx < ntimesteps-1, we take the limited values and mimick interpolation in these cases.
!      if (idx == 1) then
!        idx0 = idx    ! no previous one, use same as idx1
!        idx1 = idx
!        idx2 = idx+1
!        idx3 = idx+2
!      else if (idx == ntimesteps - 1) then
!        idx0 = idx-1
!        idx1 = idx
!        idx2 = idx+1
!        idx3 = idx+1  ! no later one, use same as idx2
!      else if (idx == ntimesteps) then
!        idx0 = idx-1
!        idx1 = idx
!        idx2 = idx    ! no later ones
!        idx3 = idx
!      else
!        ! it > 1 .and. it < ntimesteps-1
!        idx0 = idx-1
!        idx1 = idx
!        idx2 = idx+1
!        idx3 = idx+2
!      endif
!
!      ! interpolation indices
!      !idx0 = idx-1        ! 0,..
!      !idx1 = idx
!      !idx2 = idx+1        ! 2,..
!      !idx3 = idx+2        ! 3,..
!      !if (idx0 < 1) iim1 = 1
!      !if (idx2 > ntimesteps) idx2 = ntimesteps
!      !if (idx3 > ntimesteps) idx3 = ntimesteps
!
!      ! interpolates velocity/traction on all local coupling points
!      do ipoin_local = 1,npoints_local
!        ! velocity
!        ! interpolation points
!        p0(:) = tmp_veloc_points_timeseries(:,idx0,ipoin_local)
!        p1(:) = tmp_veloc_points_timeseries(:,idx1,ipoin_local)
!        p2(:) = tmp_veloc_points_timeseries(:,idx2,ipoin_local)
!        p3(:) = tmp_veloc_points_timeseries(:,idx3,ipoin_local)
!
!        ! interpolated velocity
!        Veloc_specfem(:,ipoin_local) = cs1 * p0(:) + cs2 * p1(:) + cs3 * p2(:) + cs4 * p3(:)
!
!        ! traction
!        ! interpolation points
!        p0(:) = tmp_traction_points_timeseries(:,idx0,ipoin_local)
!        p1(:) = tmp_traction_points_timeseries(:,idx1,ipoin_local)
!        p2(:) = tmp_traction_points_timeseries(:,idx2,ipoin_local)
!        p3(:) = tmp_traction_points_timeseries(:,idx3,ipoin_local)
!
!        ! interpolated traction
!        Tract_specfem(:,ipoin_local) = cs1 * p0(:) + cs2 * p1(:) + cs3 * p2(:) + cs4 * p3(:)
!      enddo
!    endif
!
!    ! store in file
!    write(IOUT) Veloc_specfem,Tract_specfem
!
!    ! user output
!    if (myrank == 0) then
!      if (mod(it_tmp,NSTEP/10) == 0) then
!        tCPU = wtime() - tstart
!        write(IMAIN,*) '    interpolated time step ',it_tmp,'out of',NSTEP,' - elapsed time:',sngl(tCPU),'s'
!        call flush_IMAIN()
!      endif
!    endif
!
!    ! debug file output
!    if (DEBUG .and. myrank == 0) then
!      ! format: # time # vel_Z
!      write(87,*) current_time, Veloc_specfem(3,1000)  ! at local point 1000
!    endif
!  enddo
!
!  ! debug file output
!  if (DEBUG .and. myrank == 0) close(87)
!
!  end subroutine interpolate_with_cubic_local


!
!-------------------------------------------------------------------------------------------------
!

  subroutine couple_with_injection_cleanup()

  use constants
  use shared_parameters, only: INJECTION_TECHNIQUE_TYPE,RECIPROCITY_AND_KH_INTEGRAL
  use specfem_par_coupling

  implicit none
  ! free coupling arrays and close files
  select case(INJECTION_TECHNIQUE_TYPE)

  case (INJECTION_TECHNIQUE_IS_DSM)
    deallocate(Veloc_dsm_boundary,Tract_dsm_boundary)
    if (old_DSM_coupling_from_Vadim) then
      close(IIN_veloc_dsm)
      close(IIN_tract_dsm)
    endif

  case (INJECTION_TECHNIQUE_IS_AXISEM)
    deallocate(Veloc_axisem,Tract_axisem)
    close(IIN_veloc_dsm)

    if (RECIPROCITY_AND_KH_INTEGRAL) then
      deallocate(Displ_axisem_time,Tract_axisem_time)
      deallocate(Displ_specfem_time,Tract_specfem_time)
      if (.not. SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
        close(IIN_displ_axisem)
      endif
    endif

  case (INJECTION_TECHNIQUE_IS_FK)
    deallocate(Veloc_FK,Tract_FK)

  case (INJECTION_TECHNIQUE_IS_SPECFEM)
    deallocate(Veloc_specfem,Tract_specfem)
    close(IIN_veloc_dsm)

  end select

  end subroutine couple_with_injection_cleanup

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fetch_injection_wavefield()

! read injection file

  use constants
  use specfem_par, only: SIMULATION_TYPE,GPU_MODE,Mesh_pointer
  use specfem_par, only: num_abs_boundary_faces,it

  ! boundary coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE,RECIPROCITY_AND_KH_INTEGRAL
  use specfem_par_coupling, only: it_dsm, &
    Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem, Tract_axisem_time, &
    Veloc_specfem, Tract_specfem
  use specfem_par_coupling, only: NP_RESAMP, Veloc_FK, Tract_FK

  implicit none
  integer :: ii, kk, iim1, iip1, iip2,npts
  real(kind=CUSTOM_REAL) :: cs1,cs2,cs3,cs4,w
  ! temporary arrays coupling
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: veloc_inj,tract_inj

  ! safety checks
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE ) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! only for forward wavefield
  if (SIMULATION_TYPE /= 1) return

  ! gets velocity & stress for boundary points
  select case(INJECTION_TECHNIQUE_TYPE)
  case (INJECTION_TECHNIQUE_IS_DSM)
    ! DSM coupling
    if (old_DSM_coupling_from_Vadim) then
      if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
        call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
      endif
    else
      !! MODIFS DE NOBU 2D
    endif

  case (INJECTION_TECHNIQUE_IS_AXISEM)
    ! AxiSEM coupling
    call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)

    !! CD CD add this
    if (RECIPROCITY_AND_KH_INTEGRAL) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)

  case (INJECTION_TECHNIQUE_IS_SPECFEM)
    ! SPECFEM coupling
    call read_specfem_file(Veloc_specfem,Tract_specfem,num_abs_boundary_faces*NGLLSQUARE,it)
  end select

  ! return for cpu mode
  if (.not. GPU_MODE) return

  ! copy injection field to GPU
  ! allcoate space
  allocate(veloc_inj(NDIM,num_abs_boundary_faces*NGLLSQUARE), &
           tract_inj(NDIM,num_abs_boundary_faces*NGLLSQUARE))

  select case(INJECTION_TECHNIQUE_TYPE)
  case (INJECTION_TECHNIQUE_IS_DSM)
    ! DSM
    veloc_inj(:,:) = reshape(Veloc_dsm_boundary(:,it_dsm,:,:),(/NDIM,num_abs_boundary_faces*NGLLSQUARE/))
    tract_inj(:,:) = reshape(Tract_dsm_boundary(:,it_dsm,:,:),(/NDIM,num_abs_boundary_faces*NGLLSQUARE/))
    it_dsm = it_dsm + 1

  case (INJECTION_TECHNIQUE_IS_AXISEM)
    ! AxiSEM coupling
    veloc_inj(:,:) = Veloc_axisem(:,:)
    tract_inj(:,:) = Tract_axisem(:,:)

    !! CD CD add this
    if (RECIPROCITY_AND_KH_INTEGRAL) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)

  case (INJECTION_TECHNIQUE_IS_SPECFEM)
    ! SPECFEM coupling
    veloc_inj(:,:) = Veloc_specfem(:,:)
    tract_inj(:,:) = Tract_specfem(:,:)

  case (INJECTION_TECHNIQUE_IS_FK)
    ! FK coupling
    !! find indices
    ! example:
    !   np_resamp = 1 and it = 1,2,3,4,5,6, ..
    !   --> ii = 1,2,3,4,5,6,..,NSTEP
    !   np_resamp = 2 and it = 1,2,3,4,5,6, ..
    !   --> ii = 1,1,2,2,3,3,..,NSTEP/2
    ii = floor( real(it + NP_RESAMP - 1) / real( NP_RESAMP))
    ! example:
    !       kk = 1,2,1,2,1,2,,..
    kk = it - (ii-1) * NP_RESAMP
    ! example:
    !       w = 0,1/2,0,1/2,..
    w = real(dble(kk-1) / dble(NP_RESAMP),kind=CUSTOM_REAL)

    ! Cubic spline values
    cs4 = w*w*w/6.d0
    cs1 = 1.d0/6.d0 + w*(w-1.d0)/2.d0 - cs4
    cs3 = w + cs1 - 2.d0*cs4
    cs2 = 1.d0 - cs1 - cs3 - cs4

    ! interpolation indices
    iim1 = ii-1        ! 0,..
    iip1 = ii+1        ! 2,..
    iip2 = ii+2        ! 3,..

    ! B-spline approximation
    veloc_inj(:,:) = cs1 * Veloc_FK(:,:,iim1) + cs2 * Veloc_FK(:,:,ii) + &
                     cs3 * Veloc_FK(:,:,iip1) + cs4 * Veloc_FK(:,:,iip2)
    tract_inj(:,:) = cs1 * Tract_FK(:,:,iim1) + cs2 * Tract_FK(:,:,ii) + &
                     cs3 * Tract_FK(:,:,iip1) + cs4 * Tract_FK(:,:,iip2)
  end select

  ! copy injection_fields to device
  veloc_inj = -veloc_inj
  tract_inj = -tract_inj
  npts = num_abs_boundary_faces * NGLLSQUARE * NDIM

  call transfer_injection_field_to_device(npts,veloc_inj,tract_inj,Mesh_pointer)

  ! free memory
  deallocate(veloc_inj,tract_inj)

  end subroutine fetch_injection_wavefield

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)

  use constants

  implicit none

  integer igll,it_dsm
  integer iface,num_abs_boundary_faces,i,j
  real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  real(kind=CUSTOM_REAL) :: dsm_boundary_tmp(3,Ntime_step_dsm,NGLLX,NGLLY)

  it_dsm = 1

  do iface = 1,num_abs_boundary_faces
    igll = 0
    do j = 1,NGLLY
      do i = 1,NGLLX
        igll = igll + 1
        read(IIN_veloc_dsm) dsm_boundary_tmp(:,:,i,j)
        Veloc_dsm_boundary(:,:,igll,iface) = dsm_boundary_tmp(:,:,i,j)
        read(IIN_tract_dsm) dsm_boundary_tmp(:,:,i,j)
        Tract_dsm_boundary(:,:,igll,iface) = dsm_boundary_tmp(:,:,i,j)
      enddo
    enddo
  enddo

  end subroutine read_dsm_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_axisem_file(Veloc_axisem,Tract_axisem,nb)

  use constants

  implicit none

  integer :: nb
  real(kind=CUSTOM_REAL) :: Veloc_axisem(3,nb)
  real(kind=CUSTOM_REAL) :: Tract_axisem(3,nb)

  read(IIN_veloc_dsm) Veloc_axisem, Tract_axisem

  end subroutine read_axisem_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_specfem_file(Veloc_specfem,Tract_specfem,nb,it)

  use constants, only: myrank,CUSTOM_REAL,NDIM,IIN_veloc_dsm

  implicit none

  integer, intent(in) :: nb
  real(kind=CUSTOM_REAL), intent(out) :: Veloc_specfem(NDIM,nb)
  real(kind=CUSTOM_REAL), intent(out) :: Tract_specfem(NDIM,nb)
  integer,intent(in) :: it

  ! local parameters
  integer :: ier

  ! reads velocity & traction
  !read(IIN_veloc_dsm,iostat=ier) Veloc_specfem, Tract_specfem
  ! w/ direct access
  read(IIN_veloc_dsm,rec=it,iostat=ier) Veloc_specfem, Tract_specfem

  ! check
  if (ier /= 0) then
    print *,'Error: rank ',myrank,'failed to read in wavefield (veloc & traction) at time step it = ',it
    print *,'       Please check that the wavefield injection files proc***_sol_specfem.bin exists.'
    call exit_MPI(myrank,'Error reading injection wavefield file')
  endif

  end subroutine read_specfem_file

