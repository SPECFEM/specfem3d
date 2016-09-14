program add_model

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

  ! ======================================================

  integer, parameter :: NSPEC=NSPEC_AB
  logical, parameter :: MINMAX_THRESHOLD_OLD = .false.   ! threshold the old model ("current model")
  logical, parameter :: MINMAX_THRESHOLD_NEW = .true.    ! threshold the new model

  character(len=150) :: sline, m_file, fname
  character(len=150) :: ftag_file_list, ftag_list(10)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model_vp, model_vs, model_vp_new, model_vs_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model_dC, model_dB, total_model
  real(kind=CUSTOM_REAL) :: alpha0, beta0, bulk0, step_fac

  real(kind=CUSTOM_REAL) :: VS_MIN, VS_MAX, VP_MIN, VP_MAX
  real(kind=CUSTOM_REAL) :: vsmin_before, vsmax_before, vpmin_before, vpmax_before
  real(kind=CUSTOM_REAL) :: vsmin_after, vsmax_after, vpmin_after, vpmax_after
  real(kind=CUSTOM_REAL) :: vsmin_new_before, vsmax_new_before, vpmin_new_before, vpmax_new_before
  real(kind=CUSTOM_REAL) :: vsmin_new_after, vsmax_new_after, vpmin_new_after, vpmax_new_after

  integer :: nfile, myrank, sizeprocs,  ier
  integer :: i,j,k,ispec,iglob, ishell, n, it, j1, ib, npts_sem, ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  ! KEY: subjective multiplication factor for the model update
  ! SOCAL: dm00 1.0, dm01 1.0, dm02 0.125, dm03 0.08, dm04 0.05, dm05 0.08
  !        dm06 0.09, dm07 0.05, dm08 0.035, dm09 0.04, dm10 0.04
  !        dm11 0.03, dm12  0.05, dm13 0.03, dm14 0.03, dm15 0.05
  step_fac = 0.05

  open(20,file='step_fac',status='unknown')
  write(20,'(1e24.12)') step_fac
  close(20)

  if (MINMAX_THRESHOLD_OLD .or. MINMAX_THRESHOLD_NEW) then
     ! minmax wavespeed values for southern california simulations
     VS_MIN = 600.0
     VS_MAX = 4700.0
     VP_MIN = 1500.0
     VP_MAX = 8200.0

     if (myrank == 0) then
        open(19,file='VS_VP_MINMAX',status='unknown')
        write(19,'(4e24.12)') VS_MIN, VS_MAX, VP_MIN, VP_MAX
        close(19)
     endif
  endif

!!$  alpha0 = 5000.0
!!$  beta0  = 3000.0
!!$  bulk0  = sqrt( alpha0**2 - (4.0/3.0)*beta0**2 )
!!$
!!$  open(19,file='wavespeed_reference',status='unknown')
!!$  write(19,'(3e24.12)') alpha0, beta0, bulk0
!!$  close(19)

  !-----------------------------------------------------

  ! read in list of file tags showing what model iteration you are on
  nfile=0
  open(unit=20, file='INPUT/ftags', status='old',iostat=ios)
  if (ios /= 0) then
    print *,'Error opening ',trim(ftag_file_list)
    stop
  endif
  do while (1 == 1)
    read(20,'(a)',iostat=ios) sline
    if (ios /= 0) exit
    nfile=nfile+1
    ftag_list(nfile) = sline
  enddo
  close(20)

  !-----------------------------------------------------
  ! read in current model and model updates: vp, vs, dC, dB

  fname = trim(ftag_list(1))
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted')
  read(12) model_vp(:,:,:,1:nspec)
  close(12)

  fname = trim(ftag_list(2))
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted')
  read(12) model_vs(:,:,:,1:nspec)
  close(12)

  fname = trim(ftag_list(3))
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted')
  read(12) model_dC(:,:,:,1:nspec)
  close(12)

  fname = trim(ftag_list(4))
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted')
  read(12) model_dB(:,:,:,1:nspec)
  close(12)

  ! multiply model updates by a subjective factor that will change the step
  model_dC = step_fac * model_dC
  model_dB = step_fac * model_dB

  !-----------------------------------------------------

  ! compute minmax values of current model
  ! NOTE: mpi_reduce operates on the values from all procs,
  !       but the reduced value only exists on the root proc.
  call mpi_reduce(minval(model_vs(:,:,:,1:nspec)), vsmin_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs(:,:,:,1:nspec)), vsmax_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp(:,:,:,1:nspec)), vpmin_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp(:,:,:,1:nspec)), vpmax_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  if (myrank == 0) then
     open(19,file='vs_vp_minmax_before',status='unknown')
     write(19,'(4e24.12)') vsmin_before, vsmax_before, vpmin_before, vpmax_before
     close(19)
  endif

  ! threshold current model and write out the modified version
  if (MINMAX_THRESHOLD_OLD) then
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 if (model_vs(i,j,k,ispec) < VS_MIN) model_vs(i,j,k,ispec) = VS_MIN
                 if (model_vs(i,j,k,ispec) > VS_MAX) model_vs(i,j,k,ispec) = VS_MAX
                 if (model_vp(i,j,k,ispec) < VP_MIN) model_vp(i,j,k,ispec) = VP_MIN
                 if (model_vp(i,j,k,ispec) > VP_MAX) model_vp(i,j,k,ispec) = VP_MAX
              enddo
           enddo
        enddo
     enddo

     ! compute minmax values of the thresholded current model
     call mpi_reduce(minval(model_vs(:,:,:,1:nspec)), vsmin_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
     call mpi_reduce(maxval(model_vs(:,:,:,1:nspec)), vsmax_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
     call mpi_reduce(minval(model_vp(:,:,:,1:nspec)), vpmin_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
     call mpi_reduce(maxval(model_vp(:,:,:,1:nspec)), vpmax_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)

     if (myrank == 0) then
        open(19,file='vs_vp_minmax_after',status='unknown')
        write(19,'(4e24.12)') vsmin_after, vsmax_after, vpmin_after, vpmax_after
        close(19)
     endif

     fname = 'vs'
     write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
     open(12,file=trim(m_file),form='unformatted')
     write(12) model_vs(:,:,:,1:nspec)
     close(12)

     fname = 'vp'
     write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
     open(12,file=trim(m_file),form='unformatted')
     write(12) model_vp(:,:,:,1:nspec)
     close(12)

  endif

!!$  !-----------------------------------------------------
!!$  ! compute current model in terms of B = ln(beta/beta0) and C = log(c/c0)
!!$
!!$  total_model = 0.
!!$  total_model = log( model_vs(:,:,:,1:nspec) / beta0 )
!!$  fname = 'Beta'
!!$  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
!!$  open(12,file=trim(m_file),form='unformatted')
!!$  write(12) total_model(:,:,:,1:nspec)
!!$  close(12)
!!$
!!$  total_model = 0.
!!$  total_model = 0.5 * log( (model_vp**2 - (4.0/3.0)*model_vs**2) / (alpha0**2 - (4.0/3.0)*beta0**2) )
!!$  fname = 'Bulk'
!!$  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
!!$  open(12,file=trim(m_file),form='unformatted')
!!$  write(12) total_model(:,:,:,1:nspec)
!!$  close(12)

  !-----------------------------------------------------

  ! S wavespeed model
  model_vs_new = 0.
  model_vs_new = model_vs * exp( model_dB )

  ! P wavespeed model
  model_vp_new = 0.
  model_vp_new = sqrt( (4.0/3.0)* model_vs**2 * exp( 2.0*model_dB ) + &
                      (model_vp**2 - (4.0/3.0)* model_vs**2) * exp( 2.0*model_dC ) )

  ! compute minmax values of new model
  call mpi_reduce(minval(model_vs_new(:,:,:,1:nspec)), vsmin_new_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs_new(:,:,:,1:nspec)), vsmax_new_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp_new(:,:,:,1:nspec)), vpmin_new_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp_new(:,:,:,1:nspec)), vpmax_new_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)

  if (myrank == 0) then
     open(19,file='vs_vp_new_minmax_before',status='unknown')
     write(19,'(4e24.12)') vsmin_new_before, vsmax_new_before, vpmin_new_before, vpmax_new_before
     close(19)
  endif

  !-----------------------------------------------------
  ! threshold model according to minmax values specified above

  if (MINMAX_THRESHOLD_NEW) then
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 if (model_vs_new(i,j,k,ispec) < VS_MIN) model_vs_new(i,j,k,ispec) = VS_MIN
                 if (model_vs_new(i,j,k,ispec) > VS_MAX) model_vs_new(i,j,k,ispec) = VS_MAX
                 if (model_vp_new(i,j,k,ispec) < VP_MIN) model_vp_new(i,j,k,ispec) = VP_MIN
                 if (model_vp_new(i,j,k,ispec) > VP_MAX) model_vp_new(i,j,k,ispec) = VP_MAX
              enddo
           enddo
        enddo
     enddo
  endif

  !-----------------------------------------------------
  ! write out new models and their global minmax values

  ! compute minmax values of new model
  call mpi_reduce(minval(model_vs_new(:,:,:,1:nspec)), vsmin_new_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs_new(:,:,:,1:nspec)), vsmax_new_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp_new(:,:,:,1:nspec)), vpmin_new_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp_new(:,:,:,1:nspec)), vpmax_new_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)

  ! this should only be different if using MINMAX_THRESHOLD_NEW
  if (myrank == 0) then
     open(19,file='vs_vp_new_minmax_after',status='unknown')
     write(19,'(4e24.12)') vsmin_new_after, vsmax_new_after, vpmin_new_after, vpmax_new_after
     close(19)
  endif

  fname = 'vs_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) model_vs_new(:,:,:,1:nspec)
  close(12)

  fname = 'vp_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) model_vp_new(:,:,:,1:nspec)
  close(12)

  !-----------------------------------------------------
  ! compute Poisson's ratio and bulk wavespeed of the current model and new model

  ! Poisson's ratio of current model
  total_model = 0.
  total_model = ( model_vp**2 - 2.0*model_vs**2 ) / ( 2.0*model_vp**2 - 2.0*model_vs**2 )
  fname = 'poisson'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! Poisson's ratio of new model
  total_model = 0.
  total_model = ( model_vp_new**2 - 2.0*model_vs_new**2 ) / &
                ( 2.0*model_vp_new**2 - 2.0*model_vs_new**2 )
  fname = 'poisson_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! bulk wavespeed of current model
  total_model = 0.
  total_model = sqrt( model_vp**2 - (4.0/3.0)*model_vs**2 )
  fname = 'vb'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! bulk wavespeed of new model
  total_model = 0.
  total_model = sqrt( model_vp_new**2 - (4.0/3.0)*model_vs_new**2 )
  fname = 'vb_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  !-----------------------------------------------------

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program add_model


