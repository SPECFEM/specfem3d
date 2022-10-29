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
  character(len=MAX_STRING_LEN) :: dsmname

! for coupling with EXTERNAL CODE !! CD CD modify here
  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '**********************************************'
      write(IMAIN,*) '      **** USING HYBRID METHOD  ****'
      write(IMAIN,*) '**********************************************'
      write(IMAIN,*)
      write(IMAIN,*) ' using coupling with injection technique:'
      select case(INJECTION_TECHNIQUE_TYPE)
      case (INJECTION_TECHNIQUE_IS_DSM)
        write(IMAIN,*) ' type of injection technique is DSM'
      case (INJECTION_TECHNIQUE_IS_AXISEM)
        write(IMAIN,*) ' type of injection technique is AXISEM'
      case (INJECTION_TECHNIQUE_IS_FK)
        write(IMAIN,*) ' type of injection technique is FK'
      case default
        stop 'Invalid INJECTION_TECHNIQUE_TYPE chosen, must be 1 == DSM, 2 == AXISEM or 3 == FK'
      end select
      write(IMAIN,*)
      write(IMAIN,*) '**********************************************'
      write(IMAIN,*)
    endif

    call create_name_database(dsmname,myrank,TRACTION_PATH)
  endif

  ! allocates arrays for coupling
  ! note: num_abs_boundary_faces needs to be set
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then

    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then

      allocate(Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2190')
      allocate(Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2191')

      if (old_DSM_coupling_from_Vadim) then
        open(unit=IIN_veloc_dsm,file=dsmname(1:len_trim(dsmname))//'vel.bin',status='old', &
             action='read',form='unformatted',iostat=ier)
        open(unit=IIN_tract_dsm,file=dsmname(1:len_trim(dsmname))//'tract.bin',status='old', &
             action='read',form='unformatted',iostat=ier)
      else
        !! To verify for NOBU version (normally, remains empty)
      endif

    else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then

      allocate(Veloc_axisem(3,NGLLSQUARE*num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2192')
      allocate(Tract_axisem(3,NGLLSQUARE*num_abs_boundary_faces),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2193')

      open(unit=IIN_veloc_dsm,file=dsmname(1:len_trim(dsmname))//'sol_axisem',status='old', &
           action='read',form='unformatted',iostat=ier)
      write(*,*) 'OPENING ', dsmname(1:len_trim(dsmname))//'sol_axisem'

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
          open(unit=IIN_displ_axisem,file=dsmname(1:len_trim(dsmname))//'axisem_displ_for_int_KH', &
            status='old',action='read',form='unformatted',iostat=ier)

          open(unit=237,file=dsmname(1:len_trim(dsmname))//'specfem_displ_for_int_KH', &
            status='old',action='read',form='unformatted',iostat=ier)

          open(unit=238,file=dsmname(1:len_trim(dsmname))//'specfem_tract_for_int_KH', &
            status='old',action='read',form='unformatted',iostat=ier)

          write(*,*) 'OPENING ', dsmname(1:len_trim(dsmname))//'axisem_displ_for_int_KH, and the specfem disp and tract'
        endif

      endif

    endif

  else
    ! dummy arrays
    allocate(Veloc_dsm_boundary(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2198')
    allocate(Tract_dsm_boundary(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2199')
    allocate(Veloc_axisem(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2200')
    allocate(Tract_axisem(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2201')
  endif

  !! CD CD add this :
  !! We perform a first run of Specfem to save displacement and tractions of Specfem for the computation of KH integral
  !! The displ, tract, and veloc of Axisem have also to be stored
  if (SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
    open(unit=237,file=dsmname(1:len_trim(dsmname))//'specfem_displ_for_int_KH',form='unformatted')
    open(unit=238,file=dsmname(1:len_trim(dsmname))//'specfem_tract_for_int_KH',form='unformatted')
    write(*,*) 'OPENING ', dsmname(1:len_trim(dsmname))//'specfem_displ_for_int_KH, and the specfem tract to SAVE IT'
  endif

  end subroutine couple_with_injection_setup

!
!-------------------------------------------------------------------------------------------------
!

  subroutine couple_with_injection_prepare_boundary()

  use specfem_par
  use specfem_par_coupling
  use specfem_par_elastic, only: ispec_is_elastic

  implicit none

  ! local parameters
  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  integer :: ier

  !! for FK point for intialization injected wavefield
  real(kind=CUSTOM_REAL) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
  real(kind=CUSTOM_REAL) :: ray_p,Tg,DF_FK

  real(kind=CUSTOM_REAL), parameter :: TOL_ZERO_TAKEOFF = 1.e-14

!  initial setup for future FK3D calculations

  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) then

    if (myrank == 0) then
      write(IMAIN,*) "preparing injection boundary"
      call flush_IMAIN()
    endif

    ! FK boundary
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then

      ! get MPI starting time for FK
      tstart = wtime()

      ! counts total number of (local) GLL points on absorbing boundary
      call nbound(NSPEC_AB,num_abs_boundary_faces,abs_boundary_ispec,ispec_is_elastic,npt)

      !! compute the bottom midle point of the domain

      !! VM VM dealocate in case of severals runs occurs in inverse_problem program
      if (allocated(nbdglb)) deallocate(nbdglb)
      if (allocated(vx_FK))  deallocate(vx_FK)
      if (allocated(vy_FK))  deallocate(vy_FK)
      if (allocated(vz_FK))  deallocate(vz_FK)
      if (allocated(tx_FK))  deallocate(tx_FK)
      if (allocated(ty_FK))  deallocate(ty_FK)
      if (allocated(tz_FK))  deallocate(tz_FK)
      if (allocated(VX_t))   deallocate(VX_t)
      if (allocated(VY_t))   deallocate(VY_t)
      if (allocated(VZ_t))   deallocate(VZ_t)
      if (allocated(TX_t))   deallocate(TX_t)
      if (allocated(TY_t))   deallocate(TY_t)
      if (allocated(TZ_t))   deallocate(TZ_t)

      !! allocate memory for FK solution
      if (npt > 0) then
        allocate(nbdglb(npt),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2202')
        allocate(vx_FK(npt),vy_FK(npt),vz_FK(npt),tx_FK(npt),ty_FK(npt),tz_FK(npt),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2203')
      else
        allocate(nbdglb(1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2204')
        allocate(vx_FK(1),vy_FK(1),vz_FK(1),tx_FK(1),ty_FK(1),tz_FK(1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2205')
      endif

      call FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

      call ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box)

      ! send FK parameters to others MPI slices
      call bcast_all_singlei(kpsv)
      call bcast_all_singlei(nlayer)

      if (myrank > 0) then
        allocate(alpha_FK(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2206')
        allocate(beta_FK(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2207')
        allocate(mu_FK(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2208')
        allocate(h_FK(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2209')
      endif

      call bcast_all_cr(alpha_FK, nlayer)
      call bcast_all_cr(beta_FK, nlayer)
      call bcast_all_cr(mu_FK, nlayer)
      call bcast_all_cr(h_FK, nlayer)

      call bcast_all_singlecr(phi_FK)
      call bcast_all_singlecr(theta_FK)
      call bcast_all_singlecr(ff0)
      call bcast_all_singlecr(xx0)
      call bcast_all_singlecr(yy0)
      call bcast_all_singlecr(zz0)
      call bcast_all_singlecr(tt0)
      call bcast_all_singlecr(Z_REF_for_FK)
      call bcast_all_singlecr(tmax_fk)
      call bcast_all_singlel(stag)

      ! converts to rad
      phi_FK   = phi_FK * PI/180.d0
      theta_FK = theta_FK * PI/180.d0

      ! ray parameter p (according to Snell's law: sin(theta1)/v1 == sin(theta2)/v2)
      if (kpsv == 1) then
        ! P-wave
        ray_p = sin(theta_FK)/alpha_FK(nlayer)    ! for vp (i.e., alpha)
      else if (kpsv == 2) then
        ! SV-wave
        ray_p = sin(theta_FK)/beta_FK(nlayer)     ! for vs (i.e., beta)
      endif

      ! note: vertical incident (theta==0 -> p==0) is not handled.
      !       here, it limits ray parameter p to a very small value to handle the calculations
      if (abs(ray_p) < TOL_ZERO_TAKEOFF) ray_p = sign(TOL_ZERO_TAKEOFF,ray_p)

      Tg  = 1.d0 / ff0

      if (npt > 0) then

        call find_size_of_working_arrays(deltat, tmax_fk, NF_FOR_STORING, NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, DF_FK)

        if (myrank == 0) then
          write(IMAIN,*) '  number of frequencies to store = ', NF_FOR_STORING
          write(IMAIN,*) '  number of frequencies for FFT  = ', NF_FOR_FFT
          write(IMAIN,*) '  power of 2 for FFT             = ', NPOW_FOR_INTERP
          write(IMAIN,*) '  resampling rate                = ', NP_RESAMP
          write(IMAIN,*)
          write(IMAIN,*) '  new time step for F-K          = ', NP_RESAMP * deltat
          write(IMAIN,*) '  new time window length         = ', tmax_fk
          write(IMAIN,*) '  frequency step for F-K         = ', DF_FK
          write(IMAIN,*)
          write(IMAIN,*) '  total number of points on boundary = ',npt
          call flush_IMAIN()
        endif

        ! safety check with number of simulation time steps
        if (NSTEP / 2 > NF_FOR_STORING + NP_RESAMP) then
          print *,'Error: FK time window length ',tmax_fk,' and NF_for_storing ',NF_FOR_STORING
          print *,'       are too small for chosen simulation length with NSTEP = ',NSTEP
          print *
          print *,'       you could use a smaller NSTEP <= ',NF_FOR_STORING*2
          print *,'       or'
          print *,'       increase FK window length larger than ',(NSTEP/2 - NP_RESAMP) * NP_RESAMP * deltat
          print *,'       to have a NF for storing  larger than ',(NSTEP/2 - NP_RESAMP)
          stop 'Invalid FK setting'
        endif

        !! arrays for storing FK solution --------------------------------------------

        allocate(VX_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2210')
        if (ier /= 0) stop 'error while allocating VX_t'
        VX_t(:,:) = 0._CUSTOM_REAL

        allocate(VY_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2211')
        if (ier /= 0) stop 'error while allocating VY_t'
        VY_t(:,:) = 0._CUSTOM_REAL

        allocate(VZ_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2212')
        if (ier /= 0) stop 'error while allocating VZ_t'
        VZ_t(:,:) = 0._CUSTOM_REAL

        allocate(TX_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2213')
        if (ier /= 0) stop 'error while allocating TX_t'
        TX_t(:,:) = 0._CUSTOM_REAL

        allocate(TY_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2214')
        if (ier /= 0) stop 'error while allocating TY_t'
        TY_t(:,:) = 0._CUSTOM_REAL

        allocate(TZ_t(npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2215')
        if (ier /= 0) stop 'error while allocating TZ_t'
        TZ_t(:,:) = 0._CUSTOM_REAL

        call FK3D(NSPEC_AB, ibool, abs_boundary_ijk, abs_boundary_normal, &
                  abs_boundary_ispec, num_abs_boundary_faces, ispec_is_elastic, &
                  kpsv, nlayer, nstep, npt, nbdglb, &
                  ray_p, phi_FK, xx0, yy0, zz0, Tg, &
                  tt0, alpha_FK, beta_FK, mu_FK, h_FK, deltat, &
                  NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

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

      deallocate(alpha_FK, beta_FK, mu_FK, h_FK)
   endif
  endif

! * end of initial setup for future FK3D calculations *

  end subroutine couple_with_injection_prepare_boundary

!
!-------------------------------------------------------------------------------------------------
!

!! count the number of point in the mesh partition boundary : npt

  subroutine nbound(NSPEC_AB, num_abs_boundary_faces, abs_boundary_ispec, ispec_is_elastic, npt)

  use constants

  implicit none

  integer,                                    intent(inout)     :: npt
  integer,                                    intent(in)        :: NSPEC_AB, num_abs_boundary_faces
  ! elastic domain flag
  logical, dimension(NSPEC_AB),               intent(in)        :: ispec_is_elastic
  ! absorbing boundary surface
  integer, dimension(num_abs_boundary_faces), intent(in)        :: abs_boundary_ispec
  ! local parameters
  integer                                                       :: ispec, iface

  npt = 0
  do iface = 1, num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)
    if ( ispec_is_elastic(ispec) ) then
      ! reference GLL points on boundary face
      npt = npt + NGLLSQUARE
    endif
  enddo

  end subroutine nbound

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FK3D(NSPEC_AB, ibool, abs_boundary_ijk, abs_boundary_normal, &
                  abs_boundary_ispec, num_abs_boundary_faces, ispec_is_elastic, &
                  kpsv, nlayer, nstep, npt, nbdglb, &
                  ray_p, phi, xx0, yy0, zz0, Tg, &
                  tt0, alpha_FK, beta_FK, mu_FK, h_FK, deltat, &
                  NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

  use constants

  use specfem_par, only: xstore, ystore, zstore, kappastore, mustore, rhostore
  use specfem_par_coupling, only: xx, yy, zz, xi1, xim, bdlambdamu, &
                                   nmx, nmy, nmz,  Z_REF_for_FK

  implicit none

  integer              :: NSPEC_AB,kpsv,nlayer,npt,nstep,ipt
  integer              :: NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP
  ! global index
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  integer, dimension(npt) :: nbdglb

  ! source
  real(kind=CUSTOM_REAL) :: ray_p,phi,xx0,yy0,zz0,tt0
  real(kind=CUSTOM_REAL) :: DF_FK

  ! model
  real(kind=CUSTOM_REAL),dimension(nlayer) :: alpha_FK,beta_FK,mu_FK,h_FK

  real(kind=CUSTOM_REAL) :: rhotmp,kappatmp,mutmp,xi,deltat,Tg
  logical, dimension(NSPEC_AB) :: ispec_is_elastic

  ! absorbing boundary surface
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),dimension(3,NGLLSQUARE,num_abs_boundary_faces) :: abs_boundary_normal

  ! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll,ier

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  if (npt > 0) then
     allocate(xx(npt),yy(npt),zz(npt),xi1(npt),xim(npt),bdlambdamu(npt),nmx(npt),nmy(npt),nmz(npt),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2216')
  else
     allocate(xx(1),yy(1),zz(1),xi1(1),xim(1),bdlambdamu(1),nmx(1),nmy(1),nmz(1),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 2217')
  endif

  nbdglb(:) = 0
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
         nbdglb(ipt) = iglob

         xx(ipt) = xstore(iglob)
         yy(ipt) = ystore(iglob)
         zz(ipt) = zstore(iglob) -  Z_REF_for_FK  !! VM VM put z in FK system of coordinate

         nmx(ipt) = abs_boundary_normal(1,igll,iface)
         nmy(ipt) = abs_boundary_normal(2,igll,iface)
         nmz(ipt) = abs_boundary_normal(3,igll,iface)

         rhotmp   = rhostore(i,j,k,ispec)
         kappatmp = kappastore(i,j,k,ispec)
         mutmp    = mustore(i,j,k,ispec)

         xi       = mutmp/(kappatmp + 4.0/3.0 * mutmp)
         xi1(ipt) = 1.0 - 2.0 * xi
         xim(ipt) = (1.0 - xi) * mutmp
         bdlambdamu(ipt) = (3.0 * kappatmp - 2.0 * mutmp) / (6.0 * kappatmp + 2.0 * mutmp)

      enddo

    endif ! ispec_is_elastic
  enddo

  call FK(alpha_FK, beta_FK, mu_FK, h_FK, nlayer, &
          Tg, ray_p, phi, xx0, yy0, zz0, &
          tt0, deltat, nstep, npt, &
          kpsv, NF_FOR_STORING, NPOW_FOR_FFT,  NP_RESAMP, DF_FK)

  deallocate(xx, yy, zz, xi1, xim, bdlambdamu, nmx, nmy, nmz)

  end subroutine FK3D

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FK(al, be, mu, h, nlayer, &
                Tg, ray_p, phi, x0, y0, z0, &
                t0, dt, npts, np, &
                kpsv, NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

  use constants, only: myrank,CUSTOM_REAL,IMAIN,PI

  use specfem_par_coupling, only: VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t, &
                                   xx, yy, zz, xi1, xim, bdlambdamu, &
                                   nmx, nmy, nmz, NPTS_STORED, NPTS_INTERP

  implicit none

  integer,                parameter                         :: CUSTOM_CMPLX = 8
  real(kind=CUSTOM_REAL), parameter                         :: zign_neg = -1.0

  ! input and output
  integer,                                     intent(in)   :: nlayer, np, npts, kpsv
  integer                                                   :: NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP

  ! model
  real(kind=CUSTOM_REAL),  dimension(nlayer),  intent(in)   :: al(nlayer),be(nlayer),mu(nlayer),H(nlayer)

  ! source
  real(kind=CUSTOM_REAL),                      intent(in)    :: dt, ray_p, phi, x0, y0, z0, Tg, t0, DF_FK

  real(kind=CUSTOM_REAL),     dimension(:),   allocatable    :: fvec, dtmp
  complex(kind=CUSTOM_CMPLX), dimension(:,:), allocatable    :: coeff, field_f
  real(kind=CUSTOM_REAL),     dimension(:,:), allocatable    :: field
  complex(kind=CUSTOM_CMPLX), dimension(:),   allocatable    :: tmp_f1, tmp_f2, tmp_f3
  real(kind=CUSTOM_REAL),     dimension(:),   allocatable    :: tmp_t1, tmp_t2, tmp_t3, tmp_it1

  complex(kind=CUSTOM_CMPLX)                                 :: C_3,stf_coeff,a,b,c,d,delta_mat,N_mat(4,4),dx_f,dz_f,txz_f,tzz_f
  real(kind=CUSTOM_REAL)                                     :: epsil,dt_fk
  real(kind=CUSTOM_REAL)                                     :: sigma_rr,sigma_rt,sigma_rz,sigma_tt,sigma_tz,sigma_zz
  real(kind=CUSTOM_REAL)                                     :: Txx_tmp, Txy_tmp, Txz_tmp, Tyy_tmp, Tyz_tmp, Tzz_tmp
  real(kind=CUSTOM_REAL)                                     :: df,om,tdelay,eta_p,eta_s,fmax,C_1
  integer                                                    :: npow,npts2,nf,nf2,nn,ii,ip,i,j,nvar,lpts
  integer                                                    :: npoints2
  integer                                                    :: ier
  logical                                                    :: comp_stress, pout

!! DK DK here is the hardwired maximum size of the array
!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
  real(kind=CUSTOM_REAL) :: mpow(30)

  epsil = 1.0e-7
  comp_stress = .true.
  nvar = 5
  pout = .false.

  fmax = 1.0_CUSTOM_REAL/(2.0_CUSTOM_REAL * dt)    ! Nyquist frequency of specfem time serie

  !! new way to do time domain resampling
  df    = DF_FK
  nf2   = NF_FOR_STORING+1   ! number of positive frequency sample points
  nf    = 2*NF_FOR_STORING   ! number of total frequencies after symetrisation
  npts2 = nf                 ! number of samples in time serie

  !! VM VM recompute new values for new way to do
  npow = ceiling(log(npts2*1.0)/log(2.0))
  npts2 = 2**npow
  NPOW_FOR_FFT = npow

  dt_fk = 1.0_CUSTOM_REAL/(df*(npts2-1))

  !! number of points for resmpled vector
  npoints2 = NP_RESAMP*(npts2-1)+1

  ! user output
  if (myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'Entering the FK synthetics program:'
     write(IMAIN,*) '  Number of points used for FFT            = ', npts2
     write(IMAIN,*) '  Number of samples stored for FK solution = ', NF_FOR_STORING
     write(IMAIN,*) '  Total time length used for FK            = ', t0+(npts2-1)*dt_fk
     write(IMAIN,*) '  FK time step       = ', dt_fk
     write(IMAIN,*) '  FK frequency step  = ', df
     write(IMAIN,*) '  power of 2 for FFT = ', npow
     call flush_IMAIN()
  endif

  !! check if dt_fk is compatible with dt_specfem
  !!
  !!
  !!

  allocate(fvec(nf2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2218')
  fvec = 0.0_CUSTOM_REAL
  do ii = 1, nf2
    fvec(ii) = (ii-1)*df
  enddo

  allocate(coeff(2,nf2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2219')
  if (ier /= 0) stop 'error while allocating'

  allocate(field_f(nf,nvar),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2220')
  if (ier /= 0) stop 'error while allocating'

  allocate(field(npts2,nvar),dtmp(npts),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2221')
  if (ier /= 0) stop 'error while allocating'

  !! allocate debug vectors
  allocate(tmp_f1(npts2), tmp_f2(npts2), tmp_f3(npts2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2222')

  if (ier /= 0) stop 'error while allocating'
  allocate(tmp_t1(npts2), tmp_t2(npts2), tmp_t3(npts2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2223')

  if (ier /= 0) stop 'error while allocating'
  allocate(tmp_it1(npoints2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2224')

  NPTS_STORED = npts2
  NPTS_INTERP = npoints2

  tmp_t1(:) = 0.0_CUSTOM_REAL
  tmp_t2(:) = 0.0_CUSTOM_REAL
  tmp_t3(:) = 0.0_CUSTOM_REAL

  tmp_f1(:) = (0.,0.)
  tmp_f2(:) = (0.,0.)
  tmp_f3(:) = (0.,0.)

  nn = int(-t0/dt) ! what if this is not an integer number?

!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
  do i = 1,npow
    mpow(i) = 2**(npow-i)
  enddo

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  starting from ',nn,' points before time 0'
    call flush_IMAIN()
  endif

  if (kpsv == 1) then
     ! P-wave

     ! for C_3=i sin(inc) (u=[sin(inc), cos(inc)])
     C_3 = cmplx(0,1.) * ray_p * al(nlayer)      ! amp. of incoming P in the bot. layer
     eta_p = sqrt(1.0/al(nlayer)**2 - ray_p**2)  ! vertical slowness for lower layer

     if (myrank == 0) write(IMAIN,*) '  Incoming P : C_3,  ray_p, eta = ', C_3, ray_p, eta_p

     N_mat(:,:) = (0.0,0.0)

     ! find out the wave coefficients in the bottom layer for all freqs -------------------------------
     do ii = 1, nf2
        om = 2.0 * PI * fvec(ii)
        ! propagation matrix
        call compute_N_rayleigh(al,be,mu,H,nlayer,om,ray_p,sum(H(1:nlayer-1)),N_mat) !total-thickness=sum(h)

        a = N_mat(3,2); b = N_mat(3,4); c = N_mat(4,2); d = N_mat(4,4)
        delta_mat = a*d - b*c
        if (delta_mat /= 0.0) then
          coeff(1,ii) = -(d*N_mat(3,3) - b*N_mat(4,3)) / delta_mat * C_3
          coeff(2,ii) = -(-c*N_mat(3,3) + a*N_mat(4,3)) / delta_mat * C_3
        else
          coeff(1,ii) = 0.0
          coeff(2,ii) = 0.0
        endif

        !debug
        !if (ii == 1 .and. myrank == 0) then
        !  print *,'debug: rayleigh coeff ',coeff(1,ii),coeff(2,ii),delta_mat
        !  print *,'Nmat'
        !  print *,N_mat
        !endif
     enddo

     ! loop over all data points -------------------------------------------------
     do ip = 1, np  ! maybe this can be run faster by shifting t for diff. x of fixed z

        field_f = 0.0
        tdelay = ray_p * (xx(ip)-x0) * cos(phi) + ray_p * (yy(ip)-y0) * sin(phi) + eta_p * (0-z0)

        do ii = 1, nf2
           om = 2.0 * PI * fvec(ii)               !! pulsation
           stf_coeff = exp(-(om*tg/2)**2)   !! apodization window
           stf_coeff = stf_coeff * exp(cmplx(0,-1)*om*tdelay)

           !! zz(ip) is the height of point with respect to the lower layer
           call compute_N_rayleigh(al,be,mu,H,nlayer,om,ray_p,zz(ip),N_mat)

           dx_f = N_mat(1,2)*coeff(1,ii) + N_mat(1,4)*coeff(2,ii) + N_mat(1,3)*C_3  ! y_1
           dz_f = N_mat(2,2)*coeff(1,ii) + N_mat(2,4)*coeff(2,ii) + N_mat(2,3)*C_3  ! y_3

           field_f(ii,1) = stf_coeff * dx_f * cmplx(0,-1) * cmplx(0,om)             ! (i om)u_x
           field_f(ii,2) = stf_coeff * dz_f * cmplx(0,om)                           ! (i om)u_z

           if (comp_stress) then
              txz_f = N_mat(3,2)*coeff(1,ii) + N_mat(3,4)*coeff(2,ii) + N_mat(3,3)*C_3 ! tilde{y}_4
              tzz_f = N_mat(4,2)*coeff(1,ii) + N_mat(4,4)*coeff(2,ii) + N_mat(4,3)*C_3 ! tilde{y}_6
              field_f(ii,3) = stf_coeff * om * ray_p * (xi1(ip)*tzz_f - 4.0*xim(ip)*dx_f) ! T_xx
              field_f(ii,4) = stf_coeff * om * ray_p * txz_f * cmplx(0,-1)                ! T_xz
              field_f(ii,5) = stf_coeff * om * ray_p * tzz_f                              ! T_zz
           endif

        enddo

        ! pad negative f, and convert to time series
        do ii = 2, nf2-1
           field_f(nf+2-ii,:) = conjg(field_f(ii,:))
        enddo

        !! inverse fast fourier transform
        field = 0.0
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
        tmp_t1(:) = field(:,1) * cos(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vx_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,1) * sin(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vy_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,2)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vz_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        !! compute traction
        do lpts = 1, NF_FOR_STORING

           sigma_rr = field(lpts,3)
           sigma_rt = 0.0
           sigma_rz = field(lpts,4)
           sigma_zz = field(lpts,5)
           sigma_tt = bdlambdamu(ip)*(sigma_rr+sigma_zz)
           sigma_tz = 0.0

           Txx_tmp = sigma_rr * cos(phi) * cos(phi) + sigma_tt * sin(phi) * sin(phi)
           Txy_tmp = cos(phi) * sin(phi) * (sigma_rr - sigma_tt)
           Txz_tmp = sigma_rz * cos(phi)
           Tyy_tmp = sigma_rr * sin(phi) * sin(phi) + sigma_tt * cos(phi) * cos(phi)
           Tyz_tmp = sigma_rz * sin(phi)
           Tzz_tmp = sigma_zz

           !! store directly the traction
           Tx_t(ip,lpts) = Txx_tmp*nmx(ip) +  Txy_tmp*nmy(ip) +  Txz_tmp*nmz(ip)
           Ty_t(ip,lpts) = Txy_tmp*nmx(ip) +  Tyy_tmp*nmy(ip) +  Tyz_tmp*nmz(ip)
           Tz_t(ip,lpts) = Txz_tmp*nmx(ip) +  Tyz_tmp*nmy(ip) +  Tzz_tmp*nmz(ip)

        enddo

        !! store undersamped version of tractions FK solution
        tmp_t1(1:NF_FOR_STORING) = Tx_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tx_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING) = Ty_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Ty_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING) = Tz_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tz_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        ! user output
        if (myrank == 0 .and. np > 1000) then
          if (mod(ip,(np/10)) == 0) then
            write(IMAIN,*) '  done ',ip/(np/10)*10,'% points out of ',np
            call flush_IMAIN()
          endif
        endif
     enddo

  else if (kpsv == 2) then
     ! SV-wave

     ! for C_2= sin(inc) (u=[cos(inc), sin(inc)])
     C_1 = ray_p * be(nlayer)                   ! amp. of incoming S in the bot. layer
     eta_s = sqrt(1.0/be(nlayer)**2 - ray_p**2) ! vertical slowness for lower layer

     if (myrank == 0) write(IMAIN,*) '  Incoming S :  C_1,  ray_p, eta = ', C_1, ray_p, eta_s

     N_mat(:,:) = (0.0,0.0)

     ! find out the wave coefficients in the bottom layer for all freqs
     do ii = 1, nf2
        om = 2.0 * PI * fvec(ii)
        ! propagation matrix
        !if (ii == nf2) pout = .true.
        call compute_N_rayleigh(al,be,mu,H,nlayer,om,ray_p,sum(H(1:nlayer-1)),N_mat) !total-thickness=sum(h)

        a = N_mat(3,2); b = N_mat(3,4); c = N_mat(4,2); d = N_mat(4,4)
        delta_mat = a*d-b*c
        if (delta_mat /= 0.0) then
          coeff(1,ii) = -(d*N_mat(3,1) - b*N_mat(4,1)) / delta_mat * C_1
          coeff(2,ii) = -(-c*N_mat(3,1) + a*N_mat(4,1)) / delta_mat * C_1
        else
          coeff(1,ii) = 0.0
          coeff(2,ii) = 0.0
        endif
     enddo

     ! loop over all data points
     do ip = 1, np  ! maybe this can be run faster by shifting t for diff. x of fixed z

        field_f = 0.
        tdelay = ray_p * (xx(ip)-x0) * cos(phi) + ray_p * (yy(ip)-y0) * sin(phi) + eta_s * (0-z0)

        do ii = 1, nf2
           om = 2.0 * PI * fvec(ii)
           stf_coeff = exp(-(om*tg/2)**2) * exp(cmplx(0,-1)*om*tdelay)

           ! z is the height of position with respect to the lowest layer interface.
           call compute_N_rayleigh(al,be,mu,H,nlayer,om,ray_p,zz(ip),N_mat)

           dx_f = N_mat(1,2)*coeff(1,ii) + N_mat(1,4)*coeff(2,ii) + N_mat(1,1)*C_1  ! y_1
           dz_f = N_mat(2,2)*coeff(1,ii) + N_mat(2,4)*coeff(2,ii) + N_mat(2,1)*C_1  ! y_3
           field_f(ii,1) = stf_coeff * dx_f * cmplx(0,-1) * cmplx(0,om)  ! (i om)u_x(1.20)
           field_f(ii,2) = stf_coeff * dz_f * cmplx(0,om)                ! (i om)u_z

           if (comp_stress) then
              txz_f = N_mat(3,2)*coeff(1,ii) + N_mat(3,4)*coeff(2,ii) + N_mat(3,1)*C_1 ! tilde{y}_4
              tzz_f = N_mat(4,2)*coeff(1,ii) + N_mat(4,4)*coeff(2,ii) + N_mat(4,1)*C_1 ! tilde{y}_6

              field_f(ii,3) = stf_coeff * om * ray_p * (xi1(ip)*tzz_f - 4.0*xim(ip)*dx_f) ! T_xx
              field_f(ii,4) = stf_coeff * om * ray_p * txz_f * cmplx(0,-1)                ! T_xz
              field_f(ii,5) = stf_coeff * om * ray_p * tzz_f                              ! T_zz
           endif

        enddo

        ! pad negative f, and convert to time series
        do ii = 2, nf2-1
           field_f(nf+2-ii,:) = conjg(field_f(ii,:))
        enddo

        field = 0.
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
        tmp_t1(:) = field(:,1)*cos(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vx_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,1)*sin(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vy_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:) = field(:,2)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vz_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        !! compute traction
        do lpts = 1, NF_FOR_STORING

           sigma_rr = field(lpts,3)
           sigma_rt = 0.0
           sigma_rz = field(lpts,4)
           sigma_zz = field(lpts,5)
           sigma_tt = bdlambdamu(ip)*(sigma_rr+sigma_zz)
           sigma_tz = 0.0

           Txx_tmp = sigma_rr * cos(phi) * cos(phi) + sigma_tt * sin(phi) * sin(phi)
           Txy_tmp = cos(phi) * sin(phi) * (sigma_rr - sigma_tt)
           Txz_tmp = sigma_rz * cos(phi)
           Tyy_tmp = sigma_rr * sin(phi) * sin(phi) + sigma_tt * cos(phi) * cos(phi)
           Tyz_tmp = sigma_rz * sin(phi)
           Tzz_tmp = sigma_zz

           !! store directly the traction
           Tx_t(ip,lpts) = Txx_tmp*nmx(ip) +  Txy_tmp*nmy(ip) +  Txz_tmp*nmz(ip)
           Ty_t(ip,lpts) = Txy_tmp*nmx(ip) +  Tyy_tmp*nmy(ip) +  Tyz_tmp*nmz(ip)
           Tz_t(ip,lpts) = Txz_tmp*nmx(ip) +  Tyz_tmp*nmy(ip) +  Tzz_tmp*nmz(ip)
        enddo

        !! store undersamped version of tractions FK solution
        tmp_t1(1:NF_FOR_STORING) = Tx_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tx_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING) = Ty_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Ty_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING) = Tz_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tz_t(ip,1:NF_FOR_STORING) = tmp_t2(1:NF_FOR_STORING)

     enddo
  endif

  deallocate(fvec,coeff, field_f, field, dtmp)
  deallocate(tmp_f1, tmp_f2, tmp_f3, tmp_t1, tmp_t2, tmp_t3)

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

  subroutine compute_N_rayleigh(alpha, beta, mu, H, nlayer, om, ray_p, ht, Nmat)

  ! assumes that ht = 0 is the bottom interface

  use constants, only: CUSTOM_REAL

  implicit none

  ! precision for complex
  integer, parameter :: CUSTOM_CMPLX = 8

  ! input
  integer,                                       intent(in)   :: nlayer
  real(kind=CUSTOM_REAL),     dimension(nlayer), intent(in)   :: alpha, beta, mu, H
  real(kind=CUSTOM_REAL),                        intent(in)   :: om,ray_p,ht

  ! output
  complex(kind=CUSTOM_CMPLX), dimension(4,4),    intent(inout) :: Nmat(4,4)

  ! local vars
  integer                                         :: i, j, ilayer
  complex(kind=CUSTOM_CMPLX), dimension(nlayer)   :: eta_alpha, eta_beta, nu_al, nu_be
  complex(kind=CUSTOM_CMPLX), dimension(4,4)      :: Emat, Gmat, Player
  complex(kind=CUSTOM_CMPLX)                      :: ca, sa, xa, ya, cb, sb, xb, yb, g1, mul, c1, c2
  real(kind=CUSTOM_REAL),     dimension(nlayer)   :: gamma0, gamma1
  real(kind=CUSTOM_REAL),     dimension(nlayer)   :: hh

  if (nlayer < 1) stop 'nlayer has to be larger than or equal to 1'

! see: Tang et al. (2014),
!      High-resolution seismic array imaging based on an SEM-FK hybrid method,
!      GJI, 197, 369-395.
! details in the appendix

  ! note: vertical incident (p=0) is not handled in Tang et al. explanations.
  !       here, it limits p to a very small value to handle the calculations
  if (abs(ray_p) < 1.e-15) stop 'Invalid ray parameter p (cannot be zero) in compute_N_rayleigh() routine'

  do i = 1,nlayer
    eta_alpha(i) = -cmplx(0,1) * sqrt(1.0/alpha(i) + ray_p) * sqrt(1.0/alpha(i) - ray_p) ! i*vertical slowness, purely imaginary
    eta_beta(i) = -cmplx(0,1) * sqrt(1.0/beta(i) + ray_p) * sqrt(1.0/beta(i) - ray_p)
    nu_al(i) = om * eta_alpha(i)
    nu_be(i) = om * eta_beta(i) ! i * vertical wavenumber
    gamma0(i) = 2.0 * beta(i)**2 * ray_p**2
    if (gamma0(i) /= 0.0) then
      gamma1(i) = 1.0 - 1.0/gamma0(i)
    else
      gamma1(i) = - huge(1.0_CUSTOM_REAL) ! very large negative gamma
    endif
  enddo

  ! Tang et al. (2014), appendix (A10) E_0:
  ! note Emat is not omega dependent
  Emat(1,1) =  eta_beta(nlayer) / ray_p
  Emat(1,2) = -Emat(1,1)
  Emat(1,3) = 1
  Emat(1,4) = 1
  Emat(2,1) = 1
  Emat(2,2) = 1
  Emat(2,3) =  eta_alpha(nlayer) / ray_p
  Emat(2,4) = -Emat(2,3)

  Emat(3,1) = 2.0 * mu(nlayer) * gamma1(nlayer)
  Emat(3,2) = Emat(3,1)
  Emat(3,3) = 2.0 * mu(nlayer) * eta_alpha(nlayer) / ray_p
  Emat(3,4) = -Emat(3,3)
  Emat(4,1) = 2.0 * mu(nlayer) * eta_beta(nlayer) / ray_p
  Emat(4,2) = -Emat(4,1)
  Emat(4,3) = Emat(3,1)
  Emat(4,4) = Emat(3,1)

  if (ht > sum(h(1:nlayer-1))) then
    write(*,*) ' FK error '
    write(*,*) ' Z point is located in the air above the surface rather than in the solid!'
    write(*,*) ' current z :', ht, ' max z allowed : ',  sum(h(1:nlayer-1))
    stop
  endif

  ! figure out the location z with respect to layer stack
  if (ht <= 0) then ! in lower half space
    Gmat = 0.0
    Gmat(1,1) = exp(nu_be(nlayer) * ht)
    Gmat(2,2) = exp(-nu_be(nlayer) * ht)
    Gmat(3,3) = exp(nu_al(nlayer) * ht)
    Gmat(4,4) = exp(-nu_al(nlayer) * ht)
    Nmat = matmul(Emat,Gmat)
  else ! in layers
    hh = H
    ilayer = nlayer
    do j = nlayer-1 , 1 , -1
      if (ht <= sum(H(j:nlayer-1))) then
        ilayer = j; exit
      endif
    enddo
    hh(ilayer+1:nlayer-1) = H(ilayer+1:nlayer-1)
    hh(ilayer) = ht - sum(H(ilayer+1:nlayer-1))
    if (hh(ilayer) < 0) stop 'Error setting layer thickness'

    ! compute propagation matrices
    Nmat(:,:) = Emat(:,:)

    do j = nlayer-1, ilayer, -1
      c1 = nu_al(j) * hh(j)
      ca = (exp(c1) + exp(-c1))/2.0  ! cos(vp h)
      sa = (exp(c1) - exp(-c1))/2.0  ! sin(vp h)

      xa = eta_alpha(j) * sa / ray_p
      ya = ray_p * sa / eta_alpha(j)

      c2 = nu_be(j) * hh(j)
      cb = (exp(c2) + exp(-c2))/2.0  ! cos(vs h)
      sb = (exp(c2) - exp(-c2))/2.0  ! sin(vs h)

      xb = eta_beta(j) * sb / ray_p
      yb = ray_p * sb / eta_beta(j)
      g1 = gamma1(j)
      mul = mu(j)

      ! Tang et al. (2014), appendix (A8)
      Player(1,1) = ca - g1*cb
      Player(1,2) = xb - g1*ya
      Player(1,3) = (ya - xb)/(2*mul)
      Player(1,4) = (cb - ca)/(2*mul)

      Player(2,1) = xa - g1*yb
      Player(2,2) = cb - g1*ca
      Player(2,3) = (ca - cb)/(2*mul)
      Player(2,4) = (yb - xa)/(2*mul)

      Player(3,1) = 2*mul * (xa - g1**2 * yb)
      Player(3,2) = 2*mul * g1 * (cb - ca)
      Player(3,3) = ca - g1*cb
      Player(3,4) = g1*yb - xa

      Player(4,1) = 2*mul * g1 * (ca - cb)
      Player(4,2) = 2*mul * (xb - g1**2 * ya)
      Player(4,3) = g1*ya - xb
      Player(4,4) = cb - g1*ca

      !debug
      !print *,'debug: j,g1,xa,xb,ya,yb,ca,cb',j,g1,xa,xb,ya,yb,ca,cb
      !print *,'debug: j,player',j,Player

      Nmat = gamma0(j) * matmul(Player,Nmat)
    enddo
  endif

  ! debug
  !print *,'debug: Nmat '
  !do j = 1,4
  !  print *,Nmat(:,j)
  !enddo

  end subroutine compute_N_rayleigh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FFT(npow,xi,zign,dtt,mpow)

! Fourier transform
! This inputs AND outputs a complex function.
! The convention is FFT --> e^(-iwt)
! numerical factor for Plancherel theorem: planch_fac = dble(NPT * dt * dt)

  use constants, only: CUSTOM_REAL

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer :: npow
  integer :: lblock,k,FK,jh,ii,istart
  integer :: l,iblock,nblock,i,lbhalf,j,lx

  complex(kind=CUSTOM_CMPLX),dimension(*) :: xi
  complex(kind=CUSTOM_CMPLX) :: wk, hold, q

  real(kind=CUSTOM_REAL) :: zign,flx,inv_of_flx,v,dtt

!! DK DK here is the hardwired maximum size of the array
!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
  real(kind=CUSTOM_REAL) :: mpow(30)

  real(kind=CUSTOM_REAL), parameter :: PI = acos(-1.d0)

!! DK DK added this sanity check
  if (npow > 30) stop 'error: the FK FTT routine has an hardwired maximum of 30 levels'
!! DK DK in any case the line below imposes a maximum of 31, otherwise the integer 2**n will overflow

  lx = 2**npow

  do l = 1,npow

    nblock = 2**(l-1)
    lblock = lx/nblock
    lbhalf = lblock/2

    k = 0

    do iblock = 1,nblock

      FK = k
      flx = lx

      v = zign * 2.0_CUSTOM_REAL * PI * FK/flx         ! Fourier convention
      ! - sign: MATLAB convention: forward e^{-i om t}
      ! + sign: engineering convention: forward e^{i om t}
      wk = cmplx(cos(v),-sin(v))   ! sign change to -sin(v) or sin(v)
      istart = lblock*(iblock-1)

      do i = 1,lbhalf
        j  = istart+i
        jh = j+lbhalf
        q = xi(jh)*wk
        xi(jh) = xi(j)-q
        xi(j)  = xi(j)+q
      enddo

      do i = 2,npow
        ii = i
        if (k < mpow(i)) goto 4
        k = k-mpow(i)
      enddo

    4 k = k+mpow(ii)

    enddo
  enddo

  k = 0

  do j = 1,lx
    if (k < j) goto 5

    hold = xi(j)
    xi(j) = xi(k+1)
    xi(k+1) = hold

5   do i = 1,npow
      ii = i
      if (k < mpow(i)) goto 7
      k = k-mpow(i)
    enddo

7   k = k+mpow(ii)
  enddo

  ! final steps deal with dt factors
  if (zign > 0.) then      ! FORWARD FFT

    xi(1:lx) = xi(1:lx) * dtt    ! multiplication by dt

  else                     ! REVERSE FFT

    flx = flx*dtt
    inv_of_flx = 1._CUSTOM_REAL / flx

!! DK DK Aug 2016: changed to multiplication by the precomputed inverse to make the routine faster
!       xi(1:lx) = xi(1:lx) / flx         ! division by dt
    xi(1:lx) = xi(1:lx) * inv_of_flx  ! division by dt

  endif

  end subroutine FFT

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FFTinv(npow,s,zign,dtt,r,mpow)

! inverse Fourier transform -- calls FFT

  use constants, only: CUSTOM_REAL

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer :: npow, nsmp, nhalf
  real(kind=CUSTOM_REAL)  :: dtt,zign

  complex(kind=CUSTOM_CMPLX), intent(in) :: s(*)
  real(kind=CUSTOM_REAL), intent(out) :: r(*)   ! note that this is real, not double precision

!! DK DK here is the hardwired maximum size of the array
!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
  real(kind=CUSTOM_REAL) :: mpow(30)

  nsmp = 2**npow
  nhalf = nsmp/2

  call rspec(s,nhalf)   ! restructuring
  call FFT(npow,s,zign,dtt,mpow)    ! Fourier transform

  r(1:nsmp) = real(s(1:nsmp))     ! take the real part

  end subroutine FFTinv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rspec(s,np2)

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer :: np2,n,n1,i

  complex(kind=CUSTOM_CMPLX) :: s(*)

  n = 2*np2
  n1 = np2+1

  s(n1) = 0.0
  s(1)  = cmplx(real(s(1)),0.0)

  do i = 1,np2
    s(np2+i) = conjg(s(np2+2-i))
  enddo

  end subroutine rspec

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_size_of_working_arrays(deltat, tmax_fk, NF_FOR_STORING, &
                                         NF_FOR_FFT, NPOW_FOR_INTERP, np_resampling, DF_FK)

  use constants, only: CUSTOM_REAL

  implicit none
  real(kind=CUSTOM_REAL),intent(inout) :: tmax_fk
  real(kind=CUSTOM_REAL),intent(inout) :: DF_FK, deltat
  integer,               intent(inout) :: NF_FOR_STORING, NF_FOR_FFT, NPOW_FOR_INTERP, np_resampling

  real(kind=CUSTOM_REAL)               :: df, dt_min_fk, Frq_ech_Fk

  !! sampling frequency to store fk solution
  Frq_ech_Fk = 10._CUSTOM_REAL  !! WM WM TO DO PUT THIS IN PARAMETER

  !! sampling time step to store fk solution
  dt_min_fk = 1. /  Frq_ech_Fk

  !!  compute resampling rate
  np_resampling  = floor(dt_min_fk / deltat)

  !! update dt for fk with respect to integer np_resampling
  dt_min_fk = np_resampling * deltat  !! this is the time step sampling for FK storage

  !! compute number of time steps to store
  NF_FOR_STORING  = ceiling( tmax_fk / dt_min_fk)

  !! in power of two
  NF_FOR_STORING  =   ceiling(log(real(NF_FOR_STORING))/log(2.))

  !! multiply by 2 in order to do an inverse FFT
  NF_FOR_FFT      =   2** (NF_FOR_STORING+1)

  NPOW_FOR_INTERP =   NF_FOR_STORING+1
  NF_FOR_STORING  =   2** NF_FOR_STORING

  !! now we have this new time window
  tmax_fk = dt_min_fk * (NF_FOR_FFT - 1)

  !! step in frequency for fk
  df = 1. / tmax_fk

  DF_FK = df

  end subroutine find_size_of_working_arrays


!
!-------------------------------------------------------------------------------------------------
!

!! #################  INTERPOLATION ROUTINES IN TIME DOMAIN ######################################

!! compute and store spline coefficients

  subroutine compute_spline_coef_to_store(Sig, npts, spline_coeff)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,                                             intent(in)     :: npts
  real(kind=CUSTOM_REAL), dimension(npts),             intent(in)     :: Sig
  real(kind=CUSTOM_REAL), dimension(npts),             intent(inout)  :: spline_coeff


  !! computation in double precision
  double precision                                                    :: error=1.d-24
  double precision                                                    :: z1, zn, sumc
  double precision, dimension(:), allocatable                         :: c
  integer                                                             :: i, n_init, ier

  allocate(c(npts),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2225')

  ! Compute pole value
  z1 = dsqrt(3.d0)-2.d0
  c(:) = dble(Sig(:)) * (1.d0-z1) *( 1.d0 - 1.d0/ z1)

  ! Initialisation causal filter
  n_init = ceiling(log(error)/log(abs(z1)))
  sumc = c(1)
  zn = z1
  do i = 1,n_init
     sumc = sumc + zn*c(i)
     zn = zn*z1
  enddo
  c(1) = sumc

  ! Causal filter
  do i = 2,npts
     c(i) = c(i) + z1* c(i-1)
  enddo

  ! Initialisation anti-causal filter
  c(npts) = ( z1 / (z1-1.d0) ) *c(npts)
  do i = npts-1,1,-1
     c(i) = z1*(c(i+1)-c(i))
  enddo

  !! store spline coeff in CUSTOM_REAL precision
  spline_coeff(:)=c(:)

  deallocate(c)

  end subroutine compute_spline_coef_to_store


!
!-------------------------------------------------------------------------------------------------
!

!! VM VM READING INPUT FILE FOR FK MODEL

  subroutine ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box)

  use specfem_par
  use specfem_par_coupling

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box
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
  ! TIME_WINDOW      tmax_fk
  !!----------------------------------------------------------------

  ! only main process reads
  if (myrank /= 0) return

  !! set default values
  tt0 = 0.
  tmax_fk = 128.
  ff0 = 0.1
  kpsv = 1  ! 1 == P-wave / 2 == SV-wave
  position_of_wavefront_not_read = .true.
  stag = .false.

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

        allocate(alpha_FK(nlayer), beta_FK(nlayer), mu_FK(nlayer), h_FK(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2226')
        allocate(rho_fk_input(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2227')
        allocate(vp_fk_input(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2228')
        allocate(vs_fk_input(nlayer),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2229')
        allocate(ztop_fk_input(nlayer+1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2230')
        allocate(ilayer_fk_input(nlayer+1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2231')
        ilayer_fk_input(:) = -1

     case('LAYER')
        read(line, *) keyword_tmp, ilayer, rho_layer, vp_layer, vs_layer, ztop_layer

        ilayer_fk_input(ilayer) = ilayer
        rho_fk_input(ilayer) = rho_layer
        vp_fk_input(ilayer) = vp_layer
        vs_fk_input(ilayer) = vs_layer
        ztop_fk_input(ilayer) = ztop_layer

        ! checks vp,vs,rho are strictly positive
        if (vp_layer <= 0.0_CUSTOM_REAL .or. vs_layer <= 0.0_CUSTOM_REAL .or. rho_layer <= 0.0_CUSTOM_REAL) then
          print *,'Error: invalid elastic material for FK coupling in layer ',ilayer
          print *,'  vp = ',vp_layer,' vs = ',vs_layer,' rho = ',rho_layer
          print *,'must be strictly positive (in particular, vs must be > 0)'
          stop 'Invalid elastic material for FK coupling found in ReadFKModelInput() routine'
        endif

     case('INCIDENT_WAVE')
         read(line,*)  keyword_tmp, incident_wave

         select case(trim(incident_wave))
            case ('p', 'P')
               kpsv = 1
            case('sv','SV')
               kpsv = 2
            case default
               kpsv = 1
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

     case('TIME_WINDOW')
        read(line,*)  keyword_tmp, tmax_fk

     end select
  !!------------------------------------------------------------------------------------------------------
  enddo

  if (allocated(ilayer_fk_input)) then

     ilayer_fk_input(nlayer+1) = ilayer_fk_input(nlayer)
     ztop_fk_input(nlayer+1) = ztop_fk_input(nlayer)
     Z_ref_for_FK = ztop_fk_input(nlayer)

     do ilayer = 1, nlayer
        alpha_FK(ilayer) = vp_fk_input(ilayer)
        beta_FK(ilayer) = vs_fk_input(ilayer)

        mu_FK(ilayer) = rho_fk_input(ilayer) * vs_fk_input(ilayer)**2
        h_FK(ilayer) =  ztop_fk_input(ilayer) - ztop_fk_input(ilayer+1)

        if (ilayer_fk_input(ilayer) == -1) then
           write(*,*) " ERROR READING FK INPUT FILE "
           write(*,*) " MISSING LAYER ", ilayer
           stop 'Missing layer in FK model'
        endif
     enddo

     deallocate(ilayer_fk_input, rho_fk_input, vp_fk_input, vs_fk_input, ztop_fk_input)

  else

     write(*,*) " ERROR READING FK INPUT FILE "
     write(*,*) " NOT BE ABLE TO READ MODEL PROPERTIES "
     stop 'Error reading FK input file'

  endif

  !! compute position of wave front
  if (position_of_wavefront_not_read) then
    xx0 = 0.5*(Xmin_box + Xmax_box)
    yy0 = 0.5*(Ymin_box + Ymax_box)
    Radius_box = sqrt( (Xmin_box - xx0)**2 + (Ymin_box - yy0)**2)

    if (kpsv == 1) then
      ! P-wave
      wave_length_at_bottom = alpha_FK(nlayer) / ff0    ! vp
    else if (kpsv == 2) then
      ! SV-wave
      wave_length_at_bottom = beta_FK(nlayer) / ff0     ! vs
    endif

    zz0 = Zmin_box - Radius_box * sin ( abs (theta_FK) * (PI / 180.d0)  ) -  &
           3.0 * wave_length_at_bottom * cos ( abs (theta_FK) * (PI / 180.d0)  ) -  &
           Z_ref_for_FK

  endif

  write(IMAIN,*) " ********************************************** "
  write(IMAIN,*) "         USING FK INJECTION TECHNIQUE           "
  write(IMAIN,*) " ********************************************** "

  write(IMAIN,*)
  write(IMAIN,*) "         Model : " , nlayer , " layers "
  write(IMAIN,*)

  do ilayer = 1, nlayer
     write(IMAIN,'(a7, i3, 3(a6,2x,f8.3), 3x, a9, f18.5 )') &
          'layer ' , ilayer, &
          " rho  =",   mu_FK(ilayer) /  beta_FK(ilayer)**2, &
          " vp   =",   alpha_FK(ilayer), &
          " vs   =",   beta_FK(ilayer), &
          " Height =", h_FK(ilayer)
  enddo

  write(IMAIN,*)
  write(IMAIN,*)
  write(IMAIN,*) " FK phi   = ",   phi_FK,'(deg)'
  write(IMAIN,*) " FK theta = ", theta_FK,'(deg)'
  write(IMAIN,*)
  write(IMAIN,*) " Origin wavefront point FK  : ", xx0, yy0, zz0
  write(IMAIN,*) " time shift  FK             : ", tt0
  write(IMAIN,*) " Z reference for FK routine : ", Z_ref_for_FK
  write(IMAIN,*) " Window for FK computing    : ", tmax_fk
  write(IMAIN,*) " frequency max              : ", ff0
  write(IMAIN,*) " type of incoming wave (1=p), (2=sv) :",kpsv
  write(IMAIN,*)
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
