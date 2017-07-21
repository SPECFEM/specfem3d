! This program is used to compute lbfgs update direction
! Author: Hejun Zhu, hejunzhu@princeton.edu
! Princeton University, New Jersey, USA
! Last modified: Tue Aug 21 17:48:28 EDT 2012

module globe_parameter
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE
  integer,parameter:: NGLOB=NGLOB_CRUST_MANTLE
  integer,parameter:: NKERNEL=4
  integer,parameter:: m_store=5   ! stored model step 3 <= m_store <= 7

  integer:: myrank, sizeprocs,ier
  integer::iker,ispec,i,j,k
  character(len=512)::filename,dirname

  character(len=256):: kernel_name(NKERNEL)
  character(len=256):: model_name(NKERNEL)
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: ibool

end module globe_parameter

program xcompute_direction_lbfgs
  use globe_parameter
  implicit none

  include 'mpif.h'
  include '../../SHARE_FILES/HEADER_FILES/precision.h'

  integer:: iter_start,iter_current,iter_store,istore
  character(len=128):: s_iter_start,s_iter_current

  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB):: q_vector,r_vector
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB):: gradient1,gradient0,model1,model0
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB):: gradient_diff,model_diff
  real(kind=CUSTOM_REAL),dimension(128):: p,a
  real(kind=CUSTOM_REAL):: p_tmp,p_sum,a_tmp,a_sum,b_tmp,b_sum
  real(kind=CUSTOM_REAL):: b,p_k_up,p_k_down,p_k_up_sum,p_k_down_sum,p_k


  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  ! read in parameters
  call get_command_argument(1,s_iter_start)
  call get_command_argument(2,s_iter_current)
  read(s_iter_start,*) iter_start
  read(s_iter_current,*) iter_current
  if (myrank == 0) print *, 'starting iteration for this period band:',iter_start
  if (myrank == 0) print *, 'current iteration:',iter_current

  iter_store = iter_current-m_store
  if ( iter_store <= iter_start ) then
        iter_store = iter_start
  endif
  if (myrank == 0) print *, 'stored iteration:',iter_store


  kernel_name=(/"reg1_bulk_betah_kernel_precond_smooth","reg1_bulk_betav_kernel_precond_smooth","reg1_eta_kernel_precond_smooth","reg1_bulk_c_kernel_precond_smooth"/)
  model_name=(/"reg1_vsh","reg1_vsv","reg1_eta","reg1_bulk"/)

  ! initialize arrays
  a(:)=0.0
  p(:)=0.0
  gradient1(:)=0.0
  gradient0(:)=0.0
  model1(:)=0.0
  model0(:)=0.0
  gradient_diff(:)=0.0
  model_diff(:)=0.0
  q_vector(:)=0.0
  r_vector(:)=0.0

  call get_ibool
  call get_gradient(iter_current,q_vector)

  if (myrank == 0) then
     print *,'************************************************'
     print *,'*******starting backward store *****************'
     print *,'************************************************'
  endif

  do istore=iter_current-1,iter_store,-1
     call get_gradient(istore+1,gradient1)
     call get_gradient(istore,gradient0)
     call get_model(istore+1,model1)
     call get_model(istore,model0)
     gradient_diff=gradient1-gradient0
     model_diff=model1-model0

     p_tmp=sum(gradient_diff*model_diff)
     call mpi_allreduce(p_tmp,p_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
     p(istore)=1.0/p_sum

     a_tmp=sum(model_diff*q_vector)
     call mpi_allreduce(a_tmp,a_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
     a(istore)=p(istore)*a_sum

     if (myrank == 0) print *,'a,p:',a(istore),p(istore)
     q_vector=q_vector-a(istore)*gradient_diff
  enddo

  istore=iter_current-1
  call get_gradient(istore+1,gradient1)
  call get_gradient(istore,gradient0)
  call get_model(istore+1,model1)
  call get_model(istore,model0)
  gradient_diff=gradient1-gradient0
  model_diff=model1-model0

! this implements Algorithm 7.4 and equation (7.20) on page 178 of the book of
! Jorge Nocedal and Stephen Wright, "Numerical Optimization", Springer, second edition (2006)
  p_k_up=sum(gradient_diff*model_diff)
  p_k_down=sum(gradient_diff*gradient_diff)
  call mpi_allreduce(p_k_up,p_k_up_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
  call mpi_allreduce(p_k_down,p_k_down_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
  p_k=p_k_up_sum/p_k_down_sum

  if ( myrank == 0) print *,'p_k:',p_k
  r_vector=p_k*q_vector
  !r_vector=1.0*q_vector

  if (myrank == 0) then
     print *,'******************************************'
     print *,'********starting forward store ***********'
     print *,'******************************************'
  endif

  do istore=iter_store,iter_current-1,1
     call get_gradient(istore+1,gradient1)
     call get_gradient(istore,gradient0)
     call get_model(istore+1,model1)
     call get_model(istore,model0)

     gradient_diff=gradient1-gradient0
     model_diff=model1-model0

     b_tmp=sum(gradient_diff*r_vector)
     call mpi_allreduce(b_tmp,b_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
     b=p(istore)*b_sum

     if (myrank == 0) print *,'a,b:',a(istore),b

     r_vector=r_vector+model_diff*(a(istore)-b)

  enddo
  r_vector=-1.0*r_vector

  call write_gradient(iter_current,r_vector)

  call MPI_FINALIZE(ier)
end program xcompute_direction_lbfgs



subroutine get_ibool
  use globe_parameter
  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLOB)::tmp

  write(dirname,'(a)') '/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE'
  write(filename,'(a,i6.6,a)') trim(dirname)//'/proc',myrank,'_reg1_solver_data_2.bin'
  open(1001,file=trim(filename),status='old',form='unformatted',iostat=ier)
  if ( ier /= 0 ) call exit_mpi(myrank,'error opening solver2 file')
  read(1001) tmp(1:NGLOB)
  read(1001) tmp(1:NGLOB)
  read(1001) tmp(1:NGLOB)
  read(1001) ibool(:,:,:,1:NSPEC)
  close(1001)
end subroutine get_ibool


subroutine get_gradient(iter,gradient)
  use globe_parameter
  implicit none
  integer::iter,iglob
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB)::gradient
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC)::vector
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLOB)::vector_gll

  do iker=1,NKERNEL
     write(dirname,'(a,i2.2)') '../SUMMED_KERNEL_M',iter
     write(filename,'(a,i6.6,a)') trim(dirname)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
     open(1001,file=trim(filename),status='old',form='unformatted',iostat=ier)
     if ( myrank == 0) print *,'reading gradient:',trim(filename)
     if (ier /= 0 ) then
        print *,'error reading:',trim(filename)
        call exit_mpi(myrank,'file not found')
     endif
     read(1001) vector(:,:,:,1:NSPEC)
     close(1001)
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob=ibool(i,j,k,ispec)
                 vector_gll(iker,iglob)=vector(i,j,k,ispec)
              enddo
            enddo
         enddo
      enddo
  enddo
  gradient(1:NGLOB)=vector_gll(1,1:NGLOB)
  gradient(NGLOB+1:2*NGLOB)=vector_gll(2,1:NGLOB)
  gradient(2*NGLOB+1:3*NGLOB)=vector_gll(3,1:NGLOB)
  gradient(3*NGLOB+1:4*NGLOB)=vector_gll(4,1:NGLOB)
end subroutine get_gradient


subroutine get_model(iter,model)
  use globe_parameter
  implicit none
  integer::iter,iglob
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB):: model
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLOB):: vector_gll
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: vector

  do iker=1,NKERNEL
     write(dirname,'(a,i2.2)') '../MODEL_M',iter
     write(filename,'(a,i6.6,a)') trim(dirname)//'/proc',myrank,'_'//trim(model_name(iker))//'.bin'
     open(1001,file=trim(filename),status='old',form='unformatted',iostat=ier)
     if ( myrank == 0) print *,'reading model:',trim(filename)
     if ( ier /= 0) then
        print *,'error reading:',trim(filename)
        call exit_mpi(myrank,'file not found')
     endif
     read(1001) vector(:,:,:,1:NSPEC)
     close(1001)
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob=ibool(i,j,k,ispec)
                 vector_gll(iker,iglob)=vector(i,j,k,ispec)
              enddo
           enddo
        enddo
     enddo
  enddo
  model(1:NGLOB)=log(vector_gll(1,1:NGLOB))
  model(NGLOB+1:2*NGLOB)=log(vector_gll(2,1:NGLOB))
  model(2*NGLOB+1:3*NGLOB)=log(vector_gll(3,1:NGLOB))
  model(3*NGLOB+1:4*NGLOB)=log(vector_gll(4,1:NGLOB))
end subroutine get_model


subroutine write_gradient(iter,gradient)
  use globe_parameter
  implicit none

  integer::iter,iglob
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB)::gradient
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLLX,NGLLY,NGLLZ,NSPEC)::vector
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLOB)::vector_gll

  vector_gll(1,1:NGLOB)=gradient(1:NGLOB)
  vector_gll(2,1:NGLOB)=gradient(NGLOB+1:2*NGLOB)
  vector_gll(3,1:NGLOB)=gradient(2*NGLOB+1:3*NGLOB)
  vector_gll(4,1:NGLOB)=gradient(3*NGLOB+1:4*NGLOB)

  do iker=1,NKERNEL
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob=ibool(i,j,k,ispec)
                 vector(iker,i,j,k,ispec)=vector_gll(iker,iglob)
              enddo
            enddo
         enddo
      enddo

      write(dirname,'(a,i2.2)') '../DIRECTION_LBFGS_M',iter
      write(filename,'(a,i6.6,a)') trim(dirname)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
      open(1001,file=trim(filename),form='unformatted',action='write')
      if ( myrank == 0) print *,'writing direct:',filename
      write(1001) vector(iker,:,:,:,1:NSPEC)
      close(1001)
  enddo
end subroutine write_gradient
