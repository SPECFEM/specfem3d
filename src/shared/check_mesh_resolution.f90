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

  subroutine check_mesh_resolution(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                   ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                                   kappastore,mustore,rhostore, &
                                   phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                   DT, model_speed_max,min_resolved_period)

! check the mesh, stability and resolved period for acoustic / elastic domains
!
! returns: maximum velocity in model ( model_speed_max ), minimum_period_resolved

  use constants
  use shared_parameters, only: SAVE_MESH_FILES,LOCAL_PATH, &
   ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: kappastore,mustore,rhostore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: phistore,tortstore
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vpI,rho_vpII,rho_vsI

  double precision,intent(in) :: DT
  real(kind=CUSTOM_REAL),intent(inout) :: model_speed_max,min_resolved_period

  ! local parameters
  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax,vpmin_glob,vpmax_glob,vsmin_glob,vsmax_glob
  real(kind=CUSTOM_REAL) :: vp2min,vp2max,vp2min_glob,vp2max_glob

  real(kind=CUSTOM_REAL) :: poissonmin,poissonmax,poissonmin_glob,poissonmax_glob
  real(kind=CUSTOM_REAL) :: distance_min,distance_max,distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min,x_max,x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min,y_max,y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min,z_max,z_min_glob,z_max_glob
  real(kind=CUSTOM_REAL) :: cmax,cmax_glob,pmax,pmax_glob
  real(kind=CUSTOM_REAL) :: dt_suggested,dt_suggested_glob,avg_distance
  real(kind=CUSTOM_REAL) :: vel_min,vel_max

  logical:: DT_PRESENT

  integer :: NSPEC_AB_global_min,NSPEC_AB_global_max,NSPEC_AB_global_sum
  integer :: NGLOB_AB_global_min,NGLOB_AB_global_max,NGLOB_AB_global_sum
  integer :: ispec

  logical :: has_vs_zero
  logical :: has_vp2_zero
  real(kind=CUSTOM_REAL),dimension(1) :: tmp_val

  ! debug: for vtk output
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp1,tmp2
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename,prname

  logical :: IMAIN_opened

  ! timing
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! checks if IMAIN is opened for writing out
  inquire(unit=IMAIN,opened=IMAIN_opened)

! note: this routine can take a long time for large meshes...
  if (myrank == 0 .and. IMAIN_opened) then
    write(IMAIN,*) "Mesh resolution:"
    call flush_IMAIN()
  endif

  ! timing gets MPI starting time
  time_start = wtime()

  ! initializations
  if (DT <= 0.0d0) then
    DT_PRESENT = .false.
  else
    DT_PRESENT = .true.
  endif

  vpmin_glob = HUGEVAL
  vpmax_glob = -HUGEVAL

  vp2min_glob = HUGEVAL
  vp2max_glob = -HUGEVAL

  vsmin_glob = HUGEVAL
  vsmax_glob = -HUGEVAL

  poissonmin_glob = HUGEVAL
  poissonmax_glob = -HUGEVAL

  distance_min_glob = HUGEVAL
  distance_max_glob = -HUGEVAL

  x_min_glob = HUGEVAL
  x_max_glob = -HUGEVAL

  y_min_glob = HUGEVAL
  y_max_glob = -HUGEVAL

  z_min_glob = HUGEVAL
  z_max_glob = -HUGEVAL

  elemsize_min_glob = HUGEVAL
  elemsize_max_glob = -HUGEVAL

  cmax_glob = -HUGEVAL
  pmax_glob = -HUGEVAL

  dt_suggested_glob = HUGEVAL

  has_vs_zero = .false.
  has_vp2_zero = .false.

  ! debug: for vtk output
  if (SAVE_MESH_FILES) then
    allocate(tmp1(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1229')
    allocate(tmp2(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1230')
    if (ier /= 0) stop 'error allocating array tmp'
    tmp1(:) = 0.0
    tmp2(:) = 0.0
  endif

!! DK DK May 2009: added this to print the minimum and maximum number of elements
!! DK DK May 2009: and points in the CUBIT + SCOTCH mesh
  call min_all_i(NSPEC_AB,NSPEC_AB_global_min)
  call max_all_i(NSPEC_AB,NSPEC_AB_global_max)
  call sum_all_i(NSPEC_AB,NSPEC_AB_global_sum)

  call min_all_i(NGLOB_AB,NGLOB_AB_global_min)
  call max_all_i(NGLOB_AB,NGLOB_AB_global_max)
  call sum_all_i(NGLOB_AB,NGLOB_AB_global_sum)

! outputs infos
  if (myrank == 0 .and. IMAIN_opened) then
    write(IMAIN,*)
    write(IMAIN,*) '********'
    write(IMAIN,*) 'minimum and maximum number of elements'
    write(IMAIN,*) 'and points in the CUBIT + SCOTCH mesh:'
    write(IMAIN,*)
    write(IMAIN,*) 'NSPEC_global_min = ',NSPEC_AB_global_min
    write(IMAIN,*) 'NSPEC_global_max = ',NSPEC_AB_global_max
    write(IMAIN,*) 'NSPEC_global_max / NSPEC_global_min imbalance = ',sngl(dble(NSPEC_AB_global_max) / dble(NSPEC_AB_global_min)), &
                      ' = ',sngl((dble(NSPEC_AB_global_max) / dble(NSPEC_AB_global_min) - 1.d0) * 100.d0),' %'
    write(IMAIN,*) 'NSPEC_global_sum = ',NSPEC_AB_global_sum
    write(IMAIN,*)
    write(IMAIN,*) 'NGLOB_global_min = ',NGLOB_AB_global_min
    write(IMAIN,*) 'NGLOB_global_max = ',NGLOB_AB_global_max
    write(IMAIN,*) 'NGLOB_global_max / NGLOB_global_min imbalance = ',sngl(dble(NGLOB_AB_global_max) / dble(NGLOB_AB_global_min)), &
                      ' = ',sngl((dble(NGLOB_AB_global_max) / dble(NGLOB_AB_global_min) - 1.d0) * 100.d0),' %'
    write(IMAIN,*) 'NGLOB_global_sum = ',NGLOB_AB_global_sum
    write(IMAIN,*)
    write(IMAIN,*) 'If you have elements of a single type (all acoustic, all elastic, all poroelastic, and without CPML)'
    write(IMAIN,*) 'in the whole mesh, then there should be no significant imbalance in the above numbers.'
    write(IMAIN,*) 'Otherwise, it is normal to have imbalance in elements and points because the domain decomposer'
    write(IMAIN,*) 'compensates for the different cost of different elements by partitioning them unevenly among processes.'
    write(IMAIN,*) '********'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks Courant number and minimum resolved period for each grid cell
  do ispec = 1,NSPEC_AB

    ! determines minimum/maximum velocities within this element
    if (ispec_is_elastic(ispec) .or. ispec_is_acoustic(ispec)) then
      ! elastic/acoustic
      call get_vpvs_minmax(vpmin,vpmax,vsmin,vsmax,poissonmin,poissonmax, &
                           ispec,has_vs_zero, &
                           NSPEC_AB,kappastore,mustore,rhostore)
      ! minimum/maximum velocity
      vel_min = min( vpmin,vsmin )
      vel_max = max( vpmax,vsmax )

      ! Poisson ratio
      poissonmin_glob = min(poissonmin_glob,poissonmin)
      poissonmax_glob = max(poissonmax_glob,poissonmax)

    else if (ispec_is_poroelastic(ispec)) then
      ! poroelastic
      call get_vpvs_minmax_poro(vpmin,vpmax,vp2min,vp2max,vsmin,vsmax, &
                                ispec,has_vs_zero,has_vp2_zero, &
                                NSPEC_AB,phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI)
      ! minimum/maximum velocity
      if (has_vp2_zero) then
        vel_min = min( vpmin,vsmin )
      else
        vel_min = min( vpmin,vp2min,vsmin )
      endif
      vel_max = max( vpmax,vp2max,vsmax )

      ! vp2
      vp2min_glob = min ( vp2min_glob, vp2min)
      vp2max_glob = max ( vp2max_glob, vp2max)

    else
      stop 'Error element type not recognized in check mesh resolution'
    endif

    ! min/max for whole cpu partition
    vpmin_glob = min(vpmin_glob, vpmin)
    vpmax_glob = max(vpmax_glob, vpmax)

    vsmin_glob = min(vsmin_glob, vsmin)
    vsmax_glob = max(vsmax_glob, vsmax)

    ! computes minimum and maximum size of this grid cell
    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                             NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    elemsize_min_glob = min(elemsize_min_glob, elemsize_min)
    elemsize_max_glob = max(elemsize_max_glob, elemsize_max)

    ! estimation of minimum period resolved
    ! based on average GLL distance within element and minimum velocity
    !
    ! rule of thumb (Komatitsch et al. 2005):
    ! "average number of points per minimum wavelength in an element should be around 5."

    ! average distance between GLL points within this element
    avg_distance = elemsize_max / ( NGLLX - 1)  ! since NGLLX = NGLLY = NGLLZ

    ! largest possible minimum period such that number of points per minimum wavelength
    ! npts = ( min(vpmin,vsmin)  * pmax ) / avg_distance  is about ~ NPTS_PER_WAVELENGTH
    !
    ! note: obviously, this estimation depends on the choice of points per wavelength
    !          which is empirical at the moment.
    !          also, keep in mind that the minimum period is just an estimation and
    !          there is no such sharp cut-off period for valid synthetics.
    !          seismograms become just more and more inaccurate for periods shorter than this estimate.
    pmax = avg_distance / vel_min * NPTS_PER_WAVELENGTH
    pmax_glob = max(pmax_glob,pmax)

    ! old: based on GLL distance, i.e. on maximum ratio ( gridspacing / velocity )
    !pmax = distance_max / min( vpmin,vsmin ) * NELEM_PER_WAVELENGTH
    !pmax_glob = max(pmax_glob,pmax)

    ! computes minimum and maximum distance of neighbor GLL points in this grid cell
    call get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
                                NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    distance_min_glob = min(distance_min_glob, distance_min)
    distance_max_glob = max(distance_max_glob, distance_max)

    ! Courant number
    ! based on minimum GLL point distance and maximum velocity
    ! i.e. on the maximum ratio of ( velocity / gridsize )
    if (DT_PRESENT) then
      cmax = vel_max * DT / distance_min
      cmax_glob = max(cmax_glob,cmax)

      ! debug: for vtk output
      if (SAVE_MESH_FILES) tmp1(ispec) = cmax
    endif

    ! suggested timestep
    dt_suggested = COURANT_SUGGESTED * distance_min / vel_max
    dt_suggested_glob = min( dt_suggested_glob, dt_suggested)

    ! debug: for vtk output
    if (SAVE_MESH_FILES) tmp2(ispec) = pmax
  enddo

  ! Vp velocity
  vpmin = vpmin_glob
  vpmax = vpmax_glob
  call min_all_cr(vpmin,vpmin_glob)
  call max_all_cr(vpmax,vpmax_glob)

  ! Vp2 velocity (relevant for poroelastic cases)
  if (POROELASTIC_SIMULATION) then
    vp2min = vp2min_glob
    if (has_vp2_zero) vp2min = 0.0
    vp2max = vp2max_glob
    call min_all_cr(vp2min,vp2min_glob)
    call max_all_cr(vp2max,vp2max_glob)
  endif

  ! Vs velocity
  vsmin = vsmin_glob
  if (has_vs_zero) vsmin = 0.0
  vsmax = vsmax_glob
  call min_all_cr(vsmin,vsmin_glob)
  call max_all_cr(vsmax,vsmax_glob)

  ! Poisson's ratio
  if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION) then
    poissonmin = poissonmin_glob
    poissonmax = poissonmax_glob
    call min_all_cr(poissonmin,poissonmin_glob)
    call max_all_cr(poissonmax,poissonmax_glob)
  endif

  ! outputs infos
  if (myrank == 0 .and. IMAIN_opened) then
    write(IMAIN,*)
    write(IMAIN,*) '********'
    write(IMAIN,*) 'Model: P   velocity min,max = ',vpmin_glob,vpmax_glob
    if (POROELASTIC_SIMULATION) then
      write(IMAIN,*) 'Model: PII velocity min,max = ',vp2min_glob,vp2max_glob
    endif
    write(IMAIN,*) 'Model: S   velocity min,max = ',vsmin_glob,vsmax_glob
    write(IMAIN,*)
    if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION) then
      write(IMAIN,*) 'Model: Poisson''s ratio min,max = ',poissonmin_glob,poissonmax_glob
      ! Poisson's ratio must be between -1 and +1/2
      if (poissonmin_glob < -1.0000001d0 .or. poissonmax_glob > 0.50000001d0) then
        write(IMAIN,*)
        write(IMAIN,*) '       Error: Poisson''s ratio is out of range (should be between -1 and +0.5).'
        write(IMAIN,*)
        stop 'Poisson''s ratio out of range'
      endif
    endif
    write(IMAIN,*) '********'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! GLL point distance
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  call min_all_cr(distance_min,distance_min_glob)
  call max_all_cr(distance_max,distance_max_glob)

  ! element size
  elemsize_min = elemsize_min_glob
  elemsize_max = elemsize_max_glob
  call min_all_cr(elemsize_min,elemsize_min_glob)
  call max_all_cr(elemsize_max,elemsize_max_glob)

  ! model dimensions
  x_min_glob = minval(xstore)
  x_max_glob = maxval(xstore)

  y_min_glob = minval(ystore)
  y_max_glob = maxval(ystore)

  z_min_glob = minval(zstore)
  z_max_glob = maxval(zstore)

  ! min and max dimensions of the model
  x_min = x_min_glob
  x_max = x_max_glob
  call min_all_cr(x_min,x_min_glob)
  call max_all_cr(x_max,x_max_glob)

  y_min = y_min_glob
  y_max = y_max_glob
  call min_all_cr(y_min,y_min_glob)
  call max_all_cr(y_max,y_max_glob)

  z_min = z_min_glob
  z_max = z_max_glob
  call min_all_cr(z_min,z_min_glob)
  call max_all_cr(z_max,z_max_glob)

  ! minimum period
  pmax = pmax_glob
  call max_all_cr(pmax,pmax_glob)

  ! time step
  dt_suggested = dt_suggested_glob
  call min_all_cr(dt_suggested,dt_suggested_glob)

  ! determines global min/max values from all cpu partitions
  if (DT_PRESENT) then
    ! Courant number
    cmax = cmax_glob
    call max_all_cr(cmax,cmax_glob)
  endif

  ! outputs infos
  if (myrank == 0 .and. IMAIN_opened) then
    write(IMAIN,*) '*********************************************'
    write(IMAIN,*) '*** Verification of simulation parameters ***'
    write(IMAIN,*) '*********************************************'
    write(IMAIN,*)
    write(IMAIN,*) '*** Xmin and Xmax of the model = ',x_min_glob,x_max_glob
    write(IMAIN,*) '*** Ymin and Ymax of the model = ',y_min_glob,y_max_glob
    write(IMAIN,*) '*** Zmin and Zmax of the model = ',z_min_glob,z_max_glob
    write(IMAIN,*)
    write(IMAIN,*) '*** Max GLL point distance = ',distance_max_glob
    write(IMAIN,*) '*** Min GLL point distance = ',distance_min_glob
    write(IMAIN,*) '*** Max/min ratio = ',distance_max_glob/distance_min_glob
    write(IMAIN,*)
    write(IMAIN,*) '*** Max element size = ',elemsize_max_glob
    write(IMAIN,*) '*** Min element size = ',elemsize_min_glob
    write(IMAIN,*) '*** Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
    write(IMAIN,*)
    write(IMAIN,*) '*** Minimum period resolved = ',pmax_glob
    write(IMAIN,*) '*** Maximum suggested time step = ',dt_suggested_glob
    write(IMAIN,*)
    if (DT_PRESENT) then
      write(IMAIN,*) '*** for DT : ',DT
      write(IMAIN,*) '*** Max stability for wave velocities = ',cmax_glob
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! checks velocities
  if (myrank == 0) then
    if (vpmin_glob <= 0.0_CUSTOM_REAL) then
      call exit_mpi(myrank,"error: vp minimum velocity")
    endif
    if (vpmax_glob >= HUGEVAL) then
      call exit_mpi(myrank,"error: vp maximum velocity")
    endif
    if (POROELASTIC_SIMULATION) then
      if (vp2min_glob < 0.0_CUSTOM_REAL) then
        call exit_mpi(myrank,"error: vp2 minimum velocity")
      endif
      if (vp2max_glob >= HUGEVAL) then
        call exit_mpi(myrank,"error: vp2 maximum velocity")
      endif
    endif
    if (vsmin_glob < 0.0_CUSTOM_REAL) then
      call exit_mpi(myrank,"error: vs minimum velocity")
    endif
    if (vsmax_glob >= HUGEVAL) then
      call exit_mpi(myrank,"error: vs maximum velocity")
    endif
  endif

  ! checks mesh
  if (myrank == 0) then
    if (distance_min_glob <= 0.0_CUSTOM_REAL) then
      call exit_mpi(myrank,"error: GLL points minimum distance")
    endif
    if (distance_max_glob >= HUGEVAL) then
      call exit_mpi(myrank,"error: GLL points maximum distance")
    endif
    if (elemsize_min_glob <= 0.0_CUSTOM_REAL) then
      call exit_mpi(myrank,"error: element minimum size")
    endif
    if (elemsize_max_glob >= HUGEVAL) then
      call exit_mpi(myrank,"error: element maximum size")
    endif
  endif

  ! returns the maximum velocity
  if (myrank == 0) then
    if (vpmax_glob > vsmax_glob) then
      model_speed_max = vpmax_glob
    else
      model_speed_max = vsmax_glob
    endif
  endif
  tmp_val(1) = model_speed_max
  call bcast_all_cr(tmp_val,1)
  model_speed_max = tmp_val(1)

  ! returns minimum period
  if (myrank == 0) min_resolved_period = pmax_glob
  tmp_val(1) = min_resolved_period
  call bcast_all_cr(tmp_val,1)
  min_resolved_period = tmp_val(1)

  ! timing
  tCPU = wtime() - time_start
  if (myrank == 0) then
    write(IMAIN,*) "Elapsed time for checking mesh resolution in seconds = ", tCPU
    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()
  endif

  ! debug: for vtk output
  if (SAVE_MESH_FILES) then
    call create_name_database(prname,myrank,LOCAL_PATH)

    ! user output
    if (myrank == 0 .and. IMAIN_opened) then
      write(IMAIN,*) 'saving VTK files for Courant number and minimum period'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! Courant number
    if (DT_PRESENT) then
      filename = trim(prname)//'res_Courant_number'
      call write_VTU_data_elem_cr_binary(NSPEC_AB,NGLOB_AB, &
                                         xstore,ystore,zstore,ibool, &
                                         tmp1,filename)
    endif
    ! minimum period estimate
    filename = trim(prname)//'res_minimum_period'
    call write_VTU_data_elem_cr_binary(NSPEC_AB,NGLOB_AB, &
                                       xstore,ystore,zstore,ibool, &
                                       tmp2,filename)

    deallocate(tmp1,tmp2)
  endif

  end subroutine check_mesh_resolution

!
!-------------------------------------------------------------------------------------------------
!

! unused, left for reference...
!
!  subroutine check_mesh_resolution_poro(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
!                                        DT, model_speed_max,min_resolved_period, &
!                                        phistore,tortstore,rhoarraystore, &
!                                        rho_vpI,rho_vpII,rho_vsI, &
!                                        LOCAL_PATH,SAVE_MESH_FILES)
!
!! check the mesh, stability and resolved period for poroelastic domains
!!
!! returns: maximum velocity in model ( model_speed_max ), minimum_period_resolved
!
!  use constants
!
!  implicit none
!
!  integer :: NSPEC_AB,NGLOB_AB
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vpI,rho_vpII,rho_vsI
!  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: phistore,tortstore
!  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
!  double precision :: DT
!  real(kind=CUSTOM_REAL) :: model_speed_max,min_resolved_period
!
!  character(len=MAX_STRING_LEN) :: LOCAL_PATH
!  logical :: SAVE_MESH_FILES
!
!  ! local parameters
!  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax,vpmin_glob,vpmax_glob,vsmin_glob,vsmax_glob
!  real(kind=CUSTOM_REAL) :: vp2min,vp2max,vp2min_glob,vp2max_glob
!  real(kind=CUSTOM_REAL) :: distance_min,distance_max,distance_min_glob,distance_max_glob
!  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_glob,elemsize_max_glob
!  real(kind=CUSTOM_REAL) :: cmax,cmax_glob,pmax,pmax_glob
!  real(kind=CUSTOM_REAL) :: dt_suggested,dt_suggested_glob,avg_distance
!
!  logical:: DT_PRESENT
!
!  integer :: NSPEC_AB_global_min,NSPEC_AB_global_max,NSPEC_AB_global_sum
!  integer :: NGLOB_AB_global_min,NGLOB_AB_global_max,NGLOB_AB_global_sum
!  integer :: ispec
!
!  logical :: has_vs_zero,has_vp2_zero
!  real(kind=CUSTOM_REAL),dimension(1) :: tmp_val
!
!  ! debug: for vtk output
!  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp1,tmp2
!
!  integer:: ier
!  character(len=MAX_STRING_LEN) :: filename,prname
!
!  ! timing
!  double precision, external :: wtime
!  double precision :: time_start,tCPU
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) "Mesh resolution: (w/ poroelastic elements)"
!    call flush_IMAIN()
!  endif
!
!  ! timing gets MPI starting time
!  time_start = wtime()
!
!  ! initializations
!  if (DT <= 0.0d0) then
!    DT_PRESENT = .false.
!  else
!    DT_PRESENT = .true.
!  endif
!
!  vpmin_glob = HUGEVAL
!  vpmax_glob = -HUGEVAL
!
!  vp2min_glob = HUGEVAL
!  vp2max_glob = -HUGEVAL
!
!  vsmin_glob = HUGEVAL
!  vsmax_glob = -HUGEVAL
!
!  distance_min_glob = HUGEVAL
!  distance_max_glob = -HUGEVAL
!
!  elemsize_min_glob = HUGEVAL
!  elemsize_max_glob = -HUGEVAL
!
!  cmax_glob = -HUGEVAL
!  pmax_glob = -HUGEVAL
!
!  dt_suggested_glob = HUGEVAL
!
!  has_vs_zero = .false.
!  has_vp2_zero = .false.
!
!  ! debug: for vtk output
!  if (SAVE_MESH_FILES) then
!    allocate(tmp1(NSPEC_AB),stat=ier)
!    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1231')
!    allocate(tmp2(NSPEC_AB),stat=ier)
!    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1232')
!    if (ier /= 0) stop 'error allocating array tmp'
!    tmp1(:) = 0.0
!    tmp2(:) = 0.0
!  endif
!
!  ! checks Courant number and minimum resolved period for each grid cell
!  do ispec = 1,NSPEC_AB
!
!    ! determines minimum/maximum velocities within this element
!    call get_vpvs_minmax_poro(vpmin,vpmax,vp2min,vp2max,vsmin,vsmax,ispec,has_vs_zero, &
!                              has_vp2_zero,NSPEC_AB, &
!                              phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI)
!
!    ! min/max for whole cpu partition
!    vpmin_glob = min ( vpmin_glob, vpmin)
!    vpmax_glob = max ( vpmax_glob, vpmax)
!
!    vp2min_glob = min ( vp2min_glob, vp2min)
!    vp2max_glob = max ( vp2max_glob, vp2max)
!
!    vsmin_glob = min ( vsmin_glob, vsmin)
!    vsmax_glob = max ( vsmax_glob, vsmax)
!
!    ! computes minimum and maximum size of this grid cell
!    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
!                             NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)
!
!    elemsize_min_glob = min( elemsize_min_glob, elemsize_min)
!    elemsize_max_glob = max( elemsize_max_glob, elemsize_max)
!
!    ! estimation of minimum period resolved
!    ! based on average GLL distance within element and minimum velocity
!    !
!    ! rule of thumb (Komatitsch et al. 2005):
!    ! "average number of points per minimum wavelength in an element should be around 5."
!
!    ! average distance between GLL points within this element
!    avg_distance = elemsize_max / ( NGLLX - 1)  ! since NGLLX = NGLLY = NGLLZ
!
!    ! largest possible minimum period such that number of points per minimum wavelength
!    ! npts = ( min(vpmin,vp2min,vsmin)  * pmax ) / avg_distance  is about ~ NPTS_PER_WAVELENGTH
!    !
!    ! note: obviously, this estimation depends on the choice of points per wavelength
!    !          which is empirical at the moment.
!    !          also, keep in mind that the minimum period is just an estimation and
!    !          there is no such sharp cut-off period for valid synthetics.
!    !          seismograms become just more and more inaccurate for periods shorter than this estimate.
!    if (has_vp2_zero) then
!      pmax = avg_distance / min( vpmin,vsmin ) * NPTS_PER_WAVELENGTH
!    else
!      pmax = avg_distance / min( vpmin,vp2min,vsmin ) * NPTS_PER_WAVELENGTH
!    endif
!    pmax_glob = max(pmax_glob,pmax)
!
!    ! old: based on GLL distance, i.e. on maximum ratio ( gridspacing / velocity )
!    !pmax = distance_max / min( vpmin,vsmin ) * NELEM_PER_WAVELENGTH
!    !pmax_glob = max(pmax_glob,pmax)
!
!    ! computes minimum and maximum distance of neighbor GLL points in this grid cell
!    call get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
!                                NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)
!
!    distance_min_glob = min( distance_min_glob, distance_min)
!    distance_max_glob = max( distance_max_glob, distance_max)
!
!    ! Courant number
!    ! based on minimum GLL point distance and maximum velocity
!    ! i.e. on the maximum ratio of ( velocity / gridsize )
!    if (DT_PRESENT) then
!      cmax = max( vpmax,vp2max,vsmax ) * DT / distance_min
!      cmax_glob = max(cmax_glob,cmax)
!
!      ! debug: for vtk output
!      if (SAVE_MESH_FILES) tmp1(ispec) = cmax
!    endif
!
!    ! suggested timestep
!    dt_suggested = COURANT_SUGGESTED * distance_min / max( vpmax,vp2max,vsmax )
!    dt_suggested_glob = min( dt_suggested_glob, dt_suggested)
!
!    ! debug: for vtk output
!    if (SAVE_MESH_FILES) tmp2(ispec) = pmax
!  enddo
!
!  ! determines global min/max values from all cpu partitions
!  if (DT_PRESENT) then
!    ! Courant number
!    cmax = cmax_glob
!    call max_all_cr(cmax,cmax_glob)
!  endif
!
!  ! minimum period
!  pmax = pmax_glob
!  call max_all_cr(pmax,pmax_glob)
!
!  ! time step
!  dt_suggested = dt_suggested_glob
!  call min_all_cr(dt_suggested,dt_suggested_glob)
!
!  ! Vp velocity
!  vpmin = vpmin_glob
!  vpmax = vpmax_glob
!  call min_all_cr(vpmin,vpmin_glob)
!  call max_all_cr(vpmax,vpmax_glob)
!
!  ! Vp2 velocity (relevant for poroelastic cases)
!  vp2min = vp2min_glob
!  if (has_vp2_zero) vp2min = 0.0
!
!  vp2max = vp2max_glob
!  call min_all_cr(vp2min,vp2min_glob)
!  call max_all_cr(vp2max,vp2max_glob)
!
!  ! Vs velocity
!  vsmin = vsmin_glob
!  if (has_vs_zero) vsmin = 0.0
!
!  vsmax = vsmax_glob
!  call min_all_cr(vsmin,vsmin_glob)
!  call max_all_cr(vsmax,vsmax_glob)
!
!  ! GLL point distance
!  distance_min = distance_min_glob
!  distance_max = distance_max_glob
!  call min_all_cr(distance_min,distance_min_glob)
!  call max_all_cr(distance_max,distance_max_glob)
!
!  ! element size
!  elemsize_min = elemsize_min_glob
!  elemsize_max = elemsize_max_glob
!  call min_all_cr(elemsize_min,elemsize_min_glob)
!  call max_all_cr(elemsize_max,elemsize_max_glob)
!
!!! DK DK May 2009: added this to print the minimum and maximum number of elements
!!! DK DK May 2009: and points in the CUBIT + SCOTCH mesh
!  call min_all_i(NSPEC_AB,NSPEC_AB_global_min)
!  call max_all_i(NSPEC_AB,NSPEC_AB_global_max)
!  call sum_all_i(NSPEC_AB,NSPEC_AB_global_sum)
!
!  call min_all_i(NGLOB_AB,NGLOB_AB_global_min)
!  call max_all_i(NGLOB_AB,NGLOB_AB_global_max)
!  call sum_all_i(NGLOB_AB,NGLOB_AB_global_sum)
!
!! outputs infos
!  if (myrank == 0) then
!    write(IMAIN,*)
!    write(IMAIN,*) '********'
!    write(IMAIN,*) 'minimum and maximum number of elements'
!    write(IMAIN,*) 'and points in the CUBIT + SCOTCH mesh:'
!    write(IMAIN,*)
!    write(IMAIN,*) 'NSPEC_AB_global_min = ',NSPEC_AB_global_min
!    write(IMAIN,*) 'NSPEC_AB_global_max = ',NSPEC_AB_global_max
!    write(IMAIN,*) 'NSPEC_AB_global_sum = ',NSPEC_AB_global_sum
!    write(IMAIN,*)
!    write(IMAIN,*) 'NGLOB_AB_global_min = ',NGLOB_AB_global_min
!    write(IMAIN,*) 'NGLOB_AB_global_max = ',NGLOB_AB_global_max
!    write(IMAIN,*) 'NGLOB_AB_global_sum = ',NGLOB_AB_global_sum
!    write(IMAIN,*)
!    write(IMAIN,*) '********'
!    write(IMAIN,*) 'Model: P velocity min,max = ',vpmin_glob,vpmax_glob
!    write(IMAIN,*) 'Model: PII velocity min,max = ',vp2min_glob,vp2max_glob
!    write(IMAIN,*) 'Model: S velocity min,max = ',vsmin_glob,vsmax_glob
!    write(IMAIN,*) '********'
!    write(IMAIN,*)
!    write(IMAIN,*) '*********************************************'
!    write(IMAIN,*) '*** Verification of simulation parameters ***'
!    write(IMAIN,*) '*********************************************'
!    write(IMAIN,*)
!    write(IMAIN,*) '*** Max GLL point distance = ',distance_max_glob
!    write(IMAIN,*) '*** Min GLL point distance = ',distance_min_glob
!    write(IMAIN,*) '*** Max/min ratio = ',distance_max_glob/distance_min_glob
!    write(IMAIN,*) '*** Max element size = ',elemsize_max_glob
!    write(IMAIN,*) '*** Min element size = ',elemsize_min_glob
!    write(IMAIN,*) '*** Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
!    write(IMAIN,*)
!    write(IMAIN,*) '*** Minimum period resolved = ',pmax_glob
!    write(IMAIN,*) '*** Maximum suggested time step = ',dt_suggested_glob
!    write(IMAIN,*)
!    if (DT_PRESENT) then
!      write(IMAIN,*) '*** for DT : ',DT
!      write(IMAIN,*) '*** Max stability for wave velocities = ',cmax_glob
!      write(IMAIN,*)
!    endif
!    call flush_IMAIN()
!  endif
!
!  ! checks velocities
!  if (myrank == 0) then
!    if (vpmin_glob <= 0.0_CUSTOM_REAL) then
!      call exit_mpi(myrank,"error: vp minimum velocity")
!    endif
!    if (vpmax_glob >= HUGEVAL) then
!      call exit_mpi(myrank,"error: vp maximum velocity")
!    endif
!    if (vp2min_glob < 0.0_CUSTOM_REAL) then
!      call exit_mpi(myrank,"error: vp2 minimum velocity")
!    endif
!    if (vp2max_glob >= HUGEVAL) then
!      call exit_mpi(myrank,"error: vp2 maximum velocity")
!    endif
!    if (vsmin_glob < 0.0_CUSTOM_REAL) then
!      call exit_mpi(myrank,"error: vs minimum velocity")
!    endif
!    if (vsmax_glob >= HUGEVAL) then
!      call exit_mpi(myrank,"error: vs maximum velocity")
!    endif
!  endif
!
!  ! checks mesh
!  if (myrank == 0) then
!    if (distance_min_glob <= 0.0_CUSTOM_REAL) then
!      call exit_mpi(myrank,"error: GLL points minimum distance")
!    endif
!    if (distance_max_glob >= HUGEVAL) then
!      call exit_mpi(myrank,"error: GLL points maximum distance")
!    endif
!    if (elemsize_min_glob <= 0.0_CUSTOM_REAL) then
!      call exit_mpi(myrank,"error: element minimum size")
!    endif
!    if (elemsize_max_glob >= HUGEVAL) then
!      call exit_mpi(myrank,"error: element maximum size")
!    endif
!  endif
!
!  ! returns the maximum velocity
!  if (myrank == 0) then
!    if (vpmax_glob > vsmax_glob) then
!      model_speed_max = vpmax_glob
!    else
!      model_speed_max = vsmax_glob
!    endif
!  endif
!  tmp_val(1) = model_speed_max
!  call bcast_all_cr(tmp_val,1)
!  model_speed_max = tmp_val(1)
!
!  ! returns minimum period
!  if (myrank == 0) min_resolved_period = pmax_glob
!  tmp_val(1) = min_resolved_period
!  call bcast_all_cr(tmp_val,1)
!  min_resolved_period = tmp_val(1)
!
!  ! timing
!  tCPU = wtime() - time_start
!  if (myrank == 0) then
!    write(IMAIN,*) "Elapsed time for checking mesh resolution in seconds = ", tCPU
!    ! flushes file buffer for main output file (IMAIN)
!    call flush_IMAIN()
!  endif
!
!  ! debug: for vtk output
!  if (SAVE_MESH_FILES) then
!    call create_name_database(prname,myrank,LOCAL_PATH)
!    ! Courant number
!    if (DT_PRESENT) then
!      filename = trim(prname)//'res_Courant_number'
!      call write_VTU_data_elem_cr_binary(NSPEC_AB,NGLOB_AB, &
!                                  xstore,ystore,zstore,ibool, &
!                                  tmp1,filename)
!    endif
!    ! minimum period estimate
!    filename = trim(prname)//'res_minimum_period'
!    call write_VTU_data_elem_cr_binary(NSPEC_AB,NGLOB_AB, &
!                                xstore,ystore,zstore,ibool, &
!                                tmp2,filename)
!    deallocate(tmp1,tmp2)
!  endif
!
!  end subroutine check_mesh_resolution_poro

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_vpvs_minmax(vpmin,vpmax,vsmin,vsmax,poissonmin,poissonmax, &
                             ispec,has_vs_zero, &
                             NSPEC_AB,kappastore,mustore,rhostore)

! calculates the min/max size of the specified element (ispec) for acoustic / elastic domains

  use constants

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: vpmin,vpmax,vsmin,vsmax,poissonmin,poissonmax

  integer :: ispec
  logical :: has_vs_zero

  integer :: NSPEC_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore,mustore,rhostore

  ! local parameters
  real(kind=CUSTOM_REAL) :: vp_sq,vs_sq,poisson
  integer :: i,j,k
  integer :: incrx,incry,incrz

  ! element looping
  ! note: the mesh resolution check is only approximative, it is thus sufficient to check only a few element points
  logical,parameter :: MIDPOINT_CHECK_ONLY = .true.

  ! initializes
  vpmin = HUGEVAL
  vpmax = -HUGEVAL
  vsmin = HUGEVAL
  vsmax = -HUGEVAL
  poissonmin = HUGEVAL
  poissonmax = -HUGEVAL

  ! looping increments
  if (MIDPOINT_CHECK_ONLY) then
    ! only corners and midpoints
    incrx = int(NGLLX/2.0)
    incry = int(NGLLY/2.0)
    incrz = int(NGLLZ/2.0)
  else
    ! all GLL points
    incrx = 1
    incry = 1
    incrz = 1
  endif

  ! loops over GLL points
  ! (avoids temporary 3d arrays)
  do k = 1,NGLLZ,incrz
    do j = 1,NGLLY,incry
      do i = 1,NGLLX,incrx

        ! vp squared
        if (rhostore(i,j,k,ispec) > TINYVAL) then
          vp_sq = (FOUR_THIRDS * mustore(i,j,k,ispec) + kappastore(i,j,k,ispec)) / rhostore(i,j,k,ispec)
        else
          vp_sq = 0.0_CUSTOM_REAL
        endif
        ! min/max
        if (vp_sq < vpmin) vpmin = vp_sq
        if (vp_sq > vpmax) vpmax = vp_sq

        ! vs squared
        if (rhostore(i,j,k,ispec) > TINYVAL) then
          vs_sq = mustore(i,j,k,ispec) / rhostore(i,j,k,ispec)
        else
          vs_sq = 0.0_CUSTOM_REAL
        endif
        ! ignore fluid regions with Vs = 0
        if (vs_sq > TINYVAL) then
          if (vs_sq < vsmin) vsmin = vs_sq
        else
          has_vs_zero = .true.
        endif
        if (vs_sq > vsmax) vsmax = vs_sq

        ! Poisson solid: for Poisson solid, the Lame parameters lambda == mu,
        !                and vp/vs = sqrt(3) => vp = sqrt(3) * vs ~ 1.73 * vs and Poisson's ratio == 0.25 (1/4)
        !
        !                note if vs = 0, then Poisson's ratio = 0.5. this value indicates a fluid or a material
        !                that maintains constant volume regardless of stress, also known as ideal incompressible solid.
        !
        !                typical values:  0.0 cork
        !                                 0.06 - 0.27 concrete
        !                                 0.27 - 0.30 steel
        !                                 0.2 sandstone
        !                                 0.3 carbonate rocks
        !                                 0.3+ shale
        !                                 0.4 coal
        !                                 0.42 - 0.44 gold
        !                                 0.5 rubber
        !
        !                for more crustal rocks, see e.g.:
        !                N. Christensen, 1996, Poisson's ratio and crustal seismology, JGR Solid Earth, vol. 101, B2.
        !                http://onlinelibrary.wiley.com/doi/10.1029/95JB03446/abstract
        !
        !                                 0.265 average continental crust
        !                                 0.253 average upper crust
        !                                 0.279 average lower crust
        !                                 0.30  average oceanic crust
        !                                 0.24 - 0.29 average crustal rocks
        !
        !                                 0.1  quartzite
        !                                 0.24 - 0.33 olivine (forsterite to fayalite)
        !                                 0.24 granite
        !                                 0.29 basalt
        !                                 0.29 gabbro
        !                                 0.34 serpentinite
        !
        !                theoretical limiting values of Poisson's ratio are -1 < sigma < 1/2
        if (has_vs_zero) then
          poisson = 0.5_CUSTOM_REAL
        else
          ! Poisson's ratio for vp and vs: \nu = 1/2 \frac{(vp/vs)^2 - 2}{(vp/vs)^2 - 1} = \frac{vp^2 - 2 vs^2}{2 vp^2 - 2 vs^2}
          if (vp_sq > TINYVAL) then
            poisson = 0.5_CUSTOM_REAL * (vp_sq - 2.0_CUSTOM_REAL * vs_sq) / (vp_sq - vs_sq)
            ! Poisson's ratio for kappa and mu: \nu = 1/2 (3 kappa - 2 mu)/(3 kappa + mu)
            ! poisson = 0.5d0 * (3.d0*kappa - 2.d0*mu)/(3.d0*kappa + mu)
          else
            ! vp not defined
            poisson = 1.0_CUSTOM_REAL
          endif
        endif
        ! min/max
        if (poisson < poissonmin) poissonmin = poisson
        if (poisson > poissonmax) poissonmax = poisson
      enddo
    enddo
  enddo

  ! takes square roots
  vpmin = sqrt(vpmin)
  vpmax = sqrt(vpmax)

  if (.not. has_vs_zero) vsmin = sqrt(vsmin)
  vsmax = sqrt(vsmax)

  end subroutine get_vpvs_minmax


!
!-------------------------------------------------------------------------------------------------
!

! unused, left for reference...

!  subroutine get_vpvs_minmax_alt(vpmin,vpmax,vsmin,vsmax,ispec,has_vs_zero, &
!                                 NSPEC_AB,kappastore,mustore,rho_vp,rho_vs)
!
!! calculates the min/max size of the specified element (ispec) for acoustic / elastic domains
!
!  use constants
!
!  implicit none
!
!  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax
!
!  integer :: ispec
!  logical :: has_vs_zero
!
!  integer :: NSPEC_AB
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
!    kappastore,mustore,rho_vp,rho_vs
!
!  ! local parameters
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ):: vp_elem,vs_elem
!  real(kind=CUSTOM_REAL), dimension(1) :: val_min,val_max
!
!  ! initializes
!  vpmin = HUGEVAL
!  vpmax = -HUGEVAL
!  vsmin = HUGEVAL
!  vsmax = -HUGEVAL
!
!  ! vp
!  where( rho_vp(:,:,:,ispec) > TINYVAL )
!    vp_elem(:,:,:) = (FOUR_THIRDS * mustore(:,:,:,ispec) &
!                    + kappastore(:,:,:,ispec)) / rho_vp(:,:,:,ispec)
!  elsewhere
!    vp_elem(:,:,:) = 0.0_CUSTOM_REAL
!  endwhere
!
!  val_min = minval(vp_elem(:,:,:))
!  val_max = maxval(vp_elem(:,:,:))
!
!  vpmin = min(vpmin,val_min(1))
!  vpmax = max(vpmax,val_max(1))
!
!  ! vs
!  where( rho_vs(:,:,:,ispec) > TINYVAL )
!    vs_elem(:,:,:) = mustore(:,:,:,ispec) / rho_vs(:,:,:,ispec)
!  elsewhere
!    vs_elem(:,:,:) = 0.0_CUSTOM_REAL
!  endwhere
!
!  val_min = minval(vs_elem(:,:,:))
!  val_max = maxval(vs_elem(:,:,:))
!
!  ! ignore fluid regions with Vs = 0
!  if (val_min(1) > 0.0001) then
!    vsmin = min(vsmin,val_min(1))
!  else
!    has_vs_zero = .true.
!  endif
!  vsmax = max(vsmax,val_max(1))
!
!  end subroutine get_vpvs_minmax_alt

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_vpvs_minmax_poro(vpmin,vpmax,vp2min,vp2max,vsmin,vsmax,ispec,has_vs_zero, &
                                  has_vp2_zero,NSPEC_AB, &
                                  phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI)

! calculates the min/max size of the specified  element (ispec) for poroelastic domains
! [CM] Note: in case of coupled acoustic-poro-elastic, rho_vpI,rho_vpII,rho_vsI, etc have
! been appropriately defined in /src/generate_databases/get_model.f90

  use constants

  implicit none

  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vp2min,vp2max,vsmin,vsmax

  integer :: ispec
  logical :: has_vs_zero
  logical :: has_vp2_zero

  integer :: NSPEC_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vpI,rho_vpII,rho_vsI
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: phistore,tortstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: vp,vp2,vs,tort,rho_s,rho_f,phi
  integer :: i,j,k
  integer :: incrx,incry,incrz

  ! element looping
  ! note: the mesh resolution check is only approximative, it is thus sufficient to check only a few element points
  logical,parameter :: MIDPOINT_CHECK_ONLY = .true.

  ! initializes
  vpmin = HUGEVAL
  vpmax = -HUGEVAL
  vp2min = HUGEVAL
  vp2max = -HUGEVAL
  vsmin = HUGEVAL
  vsmax = -HUGEVAL

  ! looping increments
  if (MIDPOINT_CHECK_ONLY) then
    ! only corners and midpoints
    incrx = int(NGLLX/2.0)
    incry = int(NGLLY/2.0)
    incrz = int(NGLLZ/2.0)
  else
    ! all GLL points
    incrx = 1
    incry = 1
    incrz = 1
  endif

  ! loops over GLL points
  ! (avoids temporary 3d arrays)
  do k = 1,NGLLZ,incrz
    do j = 1,NGLLY,incry
      do i = 1,NGLLX,incrx

        rho_s = rhoarraystore(1,i,j,k,ispec)
        rho_f = rhoarraystore(2,i,j,k,ispec)
        tort = tortstore(i,j,k,ispec)
        phi = phistore(i,j,k,ispec)

        ! vs, vpI, vpII
        if (abs(tort) > 1.e-20) then
          vs   = rho_vsI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f - phi/tort * rho_f )
          vp  = rho_vpI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f - phi/tort * rho_f )
          vp2 = rho_vpII(i,j,k,ispec)/ ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f - phi/tort * rho_f )
        else
          vs   = rho_vsI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f )
          vp  = rho_vpI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f )
          vp2 = rho_vpII(i,j,k,ispec)/ ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f )
        endif

        ! statistics
        ! vp
        vpmin = min(vpmin,vp)
        vpmax = max(vpmax,vp)

        ! vp2
        ! ignore non porous region with vp2 = 0
        vp2max = max(vp2max,vp2)
        if (vp2 > TINYVAL) then
          vp2min = min(vp2min,vp2)
        else
          has_vp2_zero = .true.
        endif

        ! vs
        ! ignore fluid regions with Vs = 0
        vsmax = max(vsmax,vs)
        if (vs > TINYVAL) then
          vsmin = min(vsmin,vs)
        else
          has_vs_zero = .true.
        endif
      enddo
    enddo
  enddo

  end subroutine get_vpvs_minmax_poro

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                  x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                                  elemsize_min_glob,elemsize_max_glob, &
                                  distance_min_glob,distance_max_glob)

! checks the mesh sizes
!
! returns: global dimensions, element size and GLL point distances

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL

  implicit none

  integer,intent(in) :: myrank
  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  real(kind=CUSTOM_REAL),intent(out) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL),intent(out) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL),intent(out) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL),intent(out) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL),intent(out) :: z_min_glob,z_max_glob

  ! local parameters
  real(kind=CUSTOM_REAL) :: distance_min,distance_max
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max
  real(kind=CUSTOM_REAL) :: x_min,x_max
  real(kind=CUSTOM_REAL) :: y_min,y_max
  real(kind=CUSTOM_REAL) :: z_min,z_max
  integer :: ispec

  ! initializes
  x_min_glob = HUGEVAL
  x_max_glob = -HUGEVAL

  y_min_glob = HUGEVAL
  y_max_glob = -HUGEVAL

  z_min_glob = HUGEVAL
  z_max_glob = -HUGEVAL

  elemsize_min_glob = HUGEVAL
  elemsize_max_glob = -HUGEVAL

  distance_min_glob = HUGEVAL
  distance_max_glob = -HUGEVAL

  ! model dimensions
  x_min_glob = minval(xstore)
  x_max_glob = maxval(xstore)

  y_min_glob = minval(ystore)
  y_max_glob = maxval(ystore)

  z_min_glob = minval(zstore)
  z_max_glob = maxval(zstore)

  ! min and max dimensions of the model
  x_min = x_min_glob
  x_max = x_max_glob
  call min_all_cr(x_min,x_min_glob)
  call max_all_cr(x_max,x_max_glob)

  y_min = y_min_glob
  y_max = y_max_glob
  call min_all_cr(y_min,y_min_glob)
  call max_all_cr(y_max,y_max_glob)

  z_min = z_min_glob
  z_max = z_max_glob
  call min_all_cr(z_min,z_min_glob)
  call max_all_cr(z_max,z_max_glob)

  ! gets distances for each grid cell
  do ispec=1,NSPEC_AB

    ! computes minimum and maximum size of this grid cell
    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                             NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    elemsize_min_glob = min(elemsize_min_glob, elemsize_min)
    elemsize_max_glob = max(elemsize_max_glob, elemsize_max)

    ! computes minimum and maximum distance of neighbor GLL points in this grid cell
    call get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
                                NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    distance_min_glob = min(distance_min_glob, distance_min)
    distance_max_glob = max(distance_max_glob, distance_max)

  enddo

  ! element size
  elemsize_min = elemsize_min_glob
  elemsize_max = elemsize_max_glob
  call min_all_cr(elemsize_min,elemsize_min_glob)
  call max_all_cr(elemsize_max,elemsize_max_glob)

  ! GLL point distance
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  call min_all_cr(distance_min,distance_min_glob)
  call max_all_cr(distance_max,distance_max_glob)

  ! checks mesh
  if (myrank == 0) then
    if (distance_min_glob <= 0.0_CUSTOM_REAL) then
      call exit_mpi(myrank,"Error GLL points minimum distance invalid")
    endif
    if (distance_max_glob >= HUGEVAL) then
      call exit_mpi(myrank,"Error GLL points maximum distance invalid")
    endif
    if (elemsize_min_glob <= 0.0_CUSTOM_REAL) then
      call exit_mpi(myrank,"Error element minimum size invalid")
    endif
    if (elemsize_max_glob >= HUGEVAL) then
      call exit_mpi(myrank,"Error element maximum size invalid")
    endif
  endif

  end subroutine check_mesh_distances

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
                          NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

! calculates the min/max distances between neighboring GLL points within the specified element (ispec);
! we purposely do not include the distance along the diagonals of the element, only along its three coordinate axes.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL

  implicit none

  real(kind=CUSTOM_REAL) :: distance_min,distance_max

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: dist
  real(kind=CUSTOM_REAL) :: x1,x2,y1,y2,z1,z2
  integer :: i,j,k,iglob1,iglob2

  ! initializes
  distance_min = HUGEVAL
  distance_max = -HUGEVAL

  ! loops over all GLL points
  ! (combines directions to speed up calculations)
  do k=1,NGLLZ-1
    do j=1,NGLLY-1
      do i=1,NGLLX-1
        ! reference point
        iglob1 = ibool(i,j,k,ispec)
        x1 = xstore(iglob1)
        y1 = ystore(iglob1)
        z1 = zstore(iglob1)

        ! along X
        iglob2 = ibool(i+1,j,k,ispec)
        x2 = xstore(iglob2)
        y2 = ystore(iglob2)
        z2 = zstore(iglob2)

        dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

        if (dist < distance_min) distance_min = dist
        if (dist > distance_max) distance_max = dist

        ! along Y
        iglob2 = ibool(i,j+1,k,ispec)
        x2 = xstore(iglob2)
        y2 = ystore(iglob2)
        z2 = zstore(iglob2)

        dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

        if (dist < distance_min) distance_min = dist
        if (dist > distance_max) distance_max = dist

        ! along Z
        iglob2 = ibool(i,j,k+1,ispec)
        x2 = xstore(iglob2)
        y2 = ystore(iglob2)
        z2 = zstore(iglob2)

        dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

        if (dist < distance_min) distance_min = dist
        if (dist > distance_max) distance_max = dist

      enddo
    enddo
  enddo

  distance_min = sqrt( distance_min )
  distance_max = sqrt( distance_max )

  end subroutine get_GLL_minmaxdistance

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                          NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

! calculates the min/max size of an edge of the specified element (ispec);
! we purposely do not include the distance along the diagonals of the element, only the size of its edges.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL

  implicit none

  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: dist
  real(kind=CUSTOM_REAL) :: x1,x2,y1,y2,z1,z2
  integer :: i,j,k
  integer :: i1,i2,j1,j2,k1,k2
  integer :: iglob1,iglob2

  ! initializes
  elemsize_min = HUGEVAL
  elemsize_max = -HUGEVAL

  ! loops over the four edges that are along X
  i1 = 1
  i2 = NGLLX
  do k = 1, NGLLZ, NGLLZ-1
    do j = 1, NGLLY, NGLLY-1
      iglob1 = ibool(i1,j,k,ispec)
      iglob2 = ibool(i2,j,k,ispec)

      x1 = xstore(iglob1)
      y1 = ystore(iglob1)
      z1 = zstore(iglob1)

      x2 = xstore(iglob2)
      y2 = ystore(iglob2)
      z2 = zstore(iglob2)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

      if (dist < elemsize_min) elemsize_min = dist
      if (dist > elemsize_max) elemsize_max = dist
    enddo
  enddo

  ! loops over the four edges that are along Y
  j1 = 1
  j2 = NGLLY
  do k = 1, NGLLZ, NGLLZ-1
    do i = 1, NGLLX, NGLLX-1
      iglob1 = ibool(i,j1,k,ispec)
      iglob2 = ibool(i,j2,k,ispec)

      x1 = xstore(iglob1)
      y1 = ystore(iglob1)
      z1 = zstore(iglob1)

      x2 = xstore(iglob2)
      y2 = ystore(iglob2)
      z2 = zstore(iglob2)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

      if (dist < elemsize_min) elemsize_min = dist
      if (dist > elemsize_max) elemsize_max = dist
    enddo
  enddo

  ! loops over the four edges that are along Z
  k1 = 1
  k2 = NGLLZ
  do j = 1, NGLLY, NGLLY-1
    do i = 1, NGLLX, NGLLX-1
      iglob1 = ibool(i,j,k1,ispec)
      iglob2 = ibool(i,j,k2,ispec)

      x1 = xstore(iglob1)
      y1 = ystore(iglob1)
      z1 = zstore(iglob1)

      x2 = xstore(iglob2)
      y2 = ystore(iglob2)
      z2 = zstore(iglob2)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

      if (dist < elemsize_min) elemsize_min = dist
      if (dist > elemsize_max) elemsize_max = dist
    enddo
  enddo

  elemsize_min = sqrt( elemsize_min )
  elemsize_max = sqrt( elemsize_max )

  end subroutine get_elem_minmaxsize

