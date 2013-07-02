!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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


! sum_preconditioned_kernels
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! this version uses the approximate Hessian stored as
!   - proc***_reg1_hess_kernel.bin
!          to precondition the transverse isotropic kernel files
!
! input file: kernels_list.txt
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernels_list.txt")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory

module sum_par

  include 'constants.h'

  ! USER PARAMETERS

  ! by default, this algorithm uses transverse isotropic (bulk,bulk_betav,bulk_betah,eta) kernels to sum up
  ! if you prefer using isotropic kernels, set flags below accordingly

  ! if you prefer using isotropic kernels (bulk,bulk_beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ISO_KERNELS = .false.

  ! if you prefer isotropic  (alpha,beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ALPHA_BETA_RHO = .false.

  ! 1 permille of maximum for inverting hessian
  real(kind=CUSTOM_REAL),parameter :: THRESHOLD_HESS = 1.e-3

  ! sums all hessians before inverting and preconditioning
  logical, parameter :: USE_HESS_SUM = .true.

  ! uses source mask to blend out source elements
  logical, parameter :: USE_SOURCE_MASK = .false.

  ! maximum number of kernels listed
  integer, parameter :: MAX_NUM_NODES = 1000

  ! default list name
  character(len=150),parameter :: kernel_file_list = '../kernels_list.txt'

  ! mesh size
  integer :: NSPEC_AB, NGLOB_AB

end module sum_par

!
!-------------------------------------------------------------------------------------------------
!

program sum_kernels

  use :: mpi
  use sum_par
  implicit none

  include 'precision.h'


  character(len=150) :: kernel_list(MAX_NUM_NODES)
  character(len=150) :: sline, kernel_name,prname_lp
  integer :: nker, myrank, sizeprocs,  ier
  integer :: ios

  double precision :: DT
  double precision :: HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML
  integer :: NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP, &
            UTM_PROJECTION_ZONE,SIMULATION_TYPE,NGNOD,NGNOD2D
  integer :: NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  integer :: NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO
  integer :: MOVIE_TYPE,IMODEL
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
            USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  logical :: ATTENUATION,USE_OLSEN_ATTENUATION, &
            APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,USE_FORCE_POINT_SOURCE
  logical :: STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD,STACEY_INSTEAD_OF_FREE_SURFACE
  logical :: ANISOTROPY,SAVE_MESH_FILES,USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION
  logical :: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE,FULL_ATTENUATION_SOLID
  character(len=256) LOCAL_PATH,TOMOGRAPHY_PATH,TRAC_PATH

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if(myrank==0) then
    write(*,*) 'sum_preconditioned_kernels:'
    write(*,*)
    write(*,*) 'reading kernel list: '
  endif
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! reads in event list
  nker=0
  open(unit = 20, file = trim(kernel_file_list), status = 'old',iostat = ios)
  if (ios /= 0) then
     print *,'Error opening ',trim(kernel_file_list),myrank
     stop 1
  endif
  do while (1 == 1)
     read(20,'(a)',iostat=ios) sline
     if (ios /= 0) exit
     nker = nker+1
     kernel_list(nker) = sline
  enddo
  close(20)
  if( myrank == 0 ) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! needs local_path for mesh files
  call read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                        USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
                        PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,FULL_ATTENUATION_SOLID,TRAC_PATH)

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROC) then
    if( myrank == 0 ) then
      print*,''
      print*,'error: run xsum_kernels with the same number of MPI processes '
      print*,'       as specified in Par_file by NPROC when slices were created'
      print*,''
      print*,'for example: mpirun -np ',NPROC,' ./xsum_kernels ...'
      print*,''
    endif
    call exit_mpi(myrank,'Error total number of slices')
  endif
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! reads mesh file
  !
  ! needs to get array dimensions

  ! opens external mesh file
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
  open(unit=27,file=trim(prname_lp),&
          status='old',action='read',form='unformatted',iostat=ios)
  if( ier /= 0 ) then
    print*,'error: could not open database '
    print*,'path: ',trim(prname_lp)
    call exit_mpi(myrank, 'error reading external mesh file')
  endif

  ! gets number of elements and global points for this partition
  read(27) NSPEC_AB
  read(27) NGLOB_AB

  close(27)

  ! user output
  if(myrank == 0) then
    print*
    print*,'  rank:',myrank,'  summing kernels'
    print*,kernel_list(1:nker)
    print*
  endif

  ! synchronizes
  call mpi_barrier(MPI_COMM_WORLD,ier)

  ! sums up kernels
  if( USE_ISO_KERNELS ) then

    if( myrank == 0 ) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    !  isotropic kernels
    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_beta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  else if( USE_ALPHA_BETA_RHO ) then

    ! isotropic kernels
    if( myrank == 0 ) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'reg1_alpha_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_beta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  else

    ! transverse isotropic kernels
    if( myrank == 0 ) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah,eta'

    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betav_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betah_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_eta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  endif

  if(myrank==0) write(*,*) 'done writing all kernels'

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program sum_kernels

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel_pre(kernel_name,kernel_list,nker,myrank)

  use :: mpi
  use sum_par
  implicit none

  include 'precision.h'

  real(kind=CUSTOM_REAL) :: norm,norm_sum
  character(len=150) :: kernel_name,kernel_list(MAX_NUM_NODES)
  integer :: nker,myrank,ier

  ! local parameters
  character(len=150) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    kernel,hess,total_kernel
  integer :: iker,ios
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: total_hess,mask_source

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          hess(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB))


  if( USE_HESS_SUM ) then
    allocate( total_hess(NGLLX,NGLLY,NGLLZ,NSPEC_AB) )
    total_hess(:,:,:,:) = 0.0_CUSTOM_REAL
  endif

  if( USE_SOURCE_MASK ) then
    allocate( mask_source(NGLLX,NGLLY,NGLLZ,NSPEC_AB) )
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! loops over all event kernels
  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if(myrank==0) then
      write(*,*) 'reading in event kernel for: ',trim(kernel_name)
      write(*,*) '  and preconditioner: ','reg1_hess_kernel'
      write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
    kernel = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_'//trim(kernel_name)//'.bin'

    open(12,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
      write(*,*) '  kernel not found:',trim(k_file)
      cycle
    endif
    read(12) kernel
    close(12)

    ! outputs norm of kernel
    norm = sum( kernel * kernel )
    call mpi_reduce(norm,norm_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
    if( myrank == 0 ) then
      print*,'  norm kernel: ',sqrt(norm_sum)
    endif


    ! approximate Hessian
    hess = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_reg1_hess_kernel.bin'

    open(12,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
      write(*,*) '  hess not found:',trim(k_file)
    else
      read(12) hess
      close(12)
    endif

    ! outputs norm of preconditioner
    norm = sum( hess * hess )
    call mpi_reduce(norm,norm_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
    if( myrank == 0 ) then
      print*,'  norm preconditioner: ',sqrt(norm_sum)
      print*
    endif

    ! note: we take absolute values for hessian (as proposed by Yang)
    hess = abs(hess)

    ! source mask
    if( USE_SOURCE_MASK ) then
      ! reads in mask
      write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                            //'/proc',myrank,'_reg1_mask_source.bin'
      open(12,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
      if( ios /= 0 ) then
        write(*,*) '  file not found:',trim(k_file)
        cycle
      endif
      read(12) mask_source
      close(12)

      ! masks source elements
      kernel = kernel * mask_source

    endif

    ! precondition
    if( USE_HESS_SUM ) then

      ! sums up hessians first
      total_hess = total_hess + hess

    else

      ! inverts hessian
      call invert_hess( myrank,hess )

      ! preconditions each event kernel with its hessian
      kernel = kernel * hess

    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel

  enddo

  ! preconditions summed kernels with summed hessians
  if( USE_HESS_SUM ) then

      ! inverts hessian matrix
      call invert_hess( myrank,total_hess )

      ! preconditions kernel
      total_kernel = total_kernel * total_hess

  endif

  ! stores summed kernels
  if(myrank==0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)

  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,'_'//trim(kernel_name)//'.bin'

  open(12,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ios)
  if( ios /= 0 ) then
    write(*,*) 'ERROR kernel not written:',trim(k_file)
    stop 'error kernel write'
  endif
  write(12) total_kernel
  close(12)

  if(myrank==0) write(*,*)

  ! frees memory
  deallocate(kernel,hess,total_kernel)
  if( USE_HESS_SUM ) deallocate(total_hess)
  if( USE_SOURCE_MASK ) deallocate(mask_source)

end subroutine sum_kernel_pre

!
!-------------------------------------------------------------------------------------------------
!

subroutine invert_hess( myrank,hess_matrix )

! inverts the hessian matrix
! the approximate hessian is only defined for diagonal elements: like
! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
! on all GLL points, which are indexed (i,j,k,ispec)

  use :: mpi
  use sum_par
  implicit none

  include 'precision.h'

  integer :: myrank
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
    hess_matrix

  ! local parameters
  real(kind=CUSTOM_REAL) :: maxh,maxh_all
  integer :: ier

  ! maximum value of hessian
  maxh = maxval( abs(hess_matrix) )

  ! determines maximum from all slices on master
  call mpi_allreduce(maxh,maxh_all,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,ier)

  ! user output
  if( myrank == 0 ) then
    print*
    print*,'hessian maximum: ',maxh_all
    print*
  endif

  ! normalizes hessian
  if( maxh_all < 1.e-18 ) then
    ! hessian is zero, re-initializes
    hess_matrix = 1.0_CUSTOM_REAL
    !call exit_mpi(myrank,'error hessian too small')
  else
    ! since hessian has absolute values, this scales between [0,1]
    hess_matrix = hess_matrix / maxh_all
  endif


  ! inverts hessian values
  where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
    hess_matrix = 1.0_CUSTOM_REAL / hess_matrix
  elsewhere
    hess_matrix = 1.0_CUSTOM_REAL / THRESHOLD_HESS
  endwhere

  ! rescales hessian
  !hess_matrix = hess_matrix * maxh_all

end subroutine invert_hess
