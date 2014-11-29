program model_pert

  implicit none
  include 'mpif.h'
  include 'constants.h'
  include 'precision.h'
  include 'values_from_mesher.h'

  ! ======================================================

  integer, parameter :: NSPEC=NSPEC_AB

  character(len=150) :: sline, m_file, fname
  character(len=150) :: ftag_file_list, ftag_list(10)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model_A, model_B, total_model
  integer :: nfile, myrank, sizeprocs, ier, n, ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  ! read in list of file tags
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
  ! read in model A and model B

  fname = trim(ftag_list(1))
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted')
  read(12) model_A(:,:,:,1:nspec)
  close(12)

  fname = trim(ftag_list(2))
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted')
  read(12) model_B(:,:,:,1:nspec)
  close(12)

  !-----------------------------------------------------
  ! compute ln(A/B)

  total_model = 0.
  total_model = log( model_A(:,:,:,1:nspec) / model_B(:,:,:,1:nspec) )
  fname = trim(ftag_list(3))
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  !-----------------------------------------------------

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program model_pert


