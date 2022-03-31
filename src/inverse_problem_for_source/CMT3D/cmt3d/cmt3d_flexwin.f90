program cmt3d_flexwin

  use cmt3d_sub
  use cmt3d_sub2
  use cmt3d_sub3
  use cmt3d_sub4

  implicit none

  character(len=150) :: par_file
  integer :: ier, i
  real*8, dimension(:,:), allocatable :: A
  real*8, dimension(:), allocatable :: b,dm

  integer :: yr,mo,jda,ho,mi
  real*8:: sec, t_cmt, hdur, elat, elon, depth, moment_tensor(NM)


  ! read and print parameters
  par_file = 'INVERSION.PAR'
  print *
  print *, 'Reading inversion parameters ...'
  call set_parameters(par_file)

  ! read cmt file
  print *, 'Reading cmtsolution parameters ...'
  call get_cmt(cmt_file,yr,mo,jda,ho,mi,sec,t_cmt,hdur,elat,elon,depth, &
               moment_tensor)
  ! assign pars to an array
  do i = 1, NPARMAX
     if (i <= NM) then
        cmt_par(i) = moment_tensor(i)
     else if (i == 7) then
        cmt_par(i) = depth
     else if (i == 8) then
        cmt_par(i) = elon
     else if (i == 9) then
        cmt_par(i) = elat
     else if (i == 10) then
        cmt_par(i) = t_cmt
     else
        cmt_par(i) = hdur
     endif
  enddo
  if (global_coord) call rotate_cmt(cmt_par,npar,elon,elat,1)

  ! data weights
  print *, 'Compute weights associated with each window ...'
  call setup_data_weights

! allocate arrays
  allocate(A(npar,npar),b(npar),dm(npar),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1093')
  if (ier /= 0) stop 'Error allocating '
  print *, 'Set up inversion matrix ...'
  call setup_matrix(A,b,npar)

! invert for the scaled parameters
  print *
  print *, 'Invert cmt parameters ...'
  call invert_cmt(A,b,dm,npar)

! calcuate misfit reduction based upon the new solution
  print *
  print *, 'Calculate formal misfit reduction ...'
  call variance_reduction(dm,npar)

! write the new cmtsolution file
  print *
  print *, 'Write new cmt solution file ...'
  if (global_coord) call rotate_cmt(new_cmt_par,npar,elon,elat,-1)

  call write_new_cmtsolution(cmt_file,new_cmt_file,new_cmt_par)

! deallocate arrays and close files
  deallocate(A,b,dm,stat=ier)
  if (ier /= 0) stop 'Error deallocating '

  print *, ' '
  print *, 'Done with cmt3d inversion'

end program cmt3d_flexwin
