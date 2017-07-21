program xcompute_direction_cg
  implicit none

  include 'mpif.h'
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'
  include '../../SHARE_FILES/HEADER_FILES/precision.h'

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE
  integer,parameter:: NKERNEL=4
  integer:: myrank, sizeprocs,ier
  integer:: iker,ispec,i,j,k

  character(len=512):: direction_0_dir, direction_1_dir, gradient_0_dir, gradient_1_dir
  character(len=512):: direction_0_file, direction_1_file, gradient_0_file, gradient_1_file
  character(len=256):: kernel_name(NKERNEL)

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: direction_0, direction_1,gradient_0,gradient_1
  real(kind=CUSTOM_REAL)::beta,beta_upper,beta_down,beta_upper_all_tmp,beta_down_all_tmp
  real(kind=CUSTOM_REAL),dimension(NKERNEL)::beta_upper_all,beta_down_all

  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  call get_command_argument(1,direction_0_dir)
  call get_command_argument(2,direction_1_dir)
  call get_command_argument(3,gradient_0_dir)
  call get_command_argument(4,gradient_1_dir)

  if (trim(direction_0_dir) == '' .or. trim(direction_1_dir) == '' &
        .or. trim(gradient_0_dir) == '' .or. trim(gradient_1_dir) == '') then
        call exit_MPI(myrank,'USAGE: xcompute_direction_cg direction_0_dir direction_1_dir gradient_0_dir gradient_1_dir')
  endif

  kernel_name=(/"reg1_bulk_betah_kernel_precond_smooth","reg1_bulk_betav_kernel_precond_smooth","reg1_bulk_c_kernel_precond_smooth","reg1_eta_kernel_precond_smooth"/)

  do iker = 1,NKERNEL
        write(gradient_0_file,'(a,i6.6,a)') trim(gradient_0_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
        write(gradient_1_file,'(a,i6.6,a)') trim(gradient_1_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'


        open(1001,file=trim(gradient_0_file),status='old',form='unformatted',iostat=ier)
        if (myrank == 0) print *, 'reading gradient0:',trim(gradient_0_file)
        if ( ier /= 0) then
                print *, 'error opening:',trim(gradient_0_file)
                call exit_mpi(myrank,'file not found')
        endif
        read(1001) gradient_0(:,:,:,1:NSPEC)
        close(1001)

        open(1001,file=trim(gradient_1_file),status='old',form='unformatted',iostat=ier)
        if (myrank == 0) print *, 'reading gradient1:',trim(gradient_1_file)
        if (ier /= 0) then
                print *, 'error opening:',trim(gradient_1_file)
                call exit_mpi(myrank,'file not found')
        endif
        read(1001) gradient_1(:,:,:,1:NSPEC)
        close(1001)

        beta_upper=sum(gradient_1*(gradient_1-gradient_0))
        beta_down=sum(gradient_0*gradient_0)

        call mpi_barrier(MPI_COMM_WORLD,ier)
        call mpi_allreduce(beta_upper,beta_upper_all_tmp,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)
        call mpi_allreduce(beta_down,beta_down_all_tmp,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ier)

        beta_upper_all(iker)=beta_upper_all_tmp
        beta_down_all(iker)=beta_down_all_tmp
  enddo

  beta=sum(beta_upper_all)/sum(beta_down_all)
  if (myrank == 0 ) then
        print *,'before zero',myrank,beta
  endif
  if ( beta < 0.0 ) then
        beta=0.0
  endif


  if (myrank == 0 ) then
        print *,myrank,beta
  endif


  do iker = 1,NKERNEL
        direction_1=0._CUSTOM_REAL

        write(direction_0_file,'(a,i6.6,a)') trim(direction_0_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
        write(direction_1_file,'(a,i6.6,a)') trim(direction_1_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
        write(gradient_0_file,'(a,i6.6,a)') trim(gradient_0_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
        write(gradient_1_file,'(a,i6.6,a)') trim(gradient_1_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'

        open(1001,file=trim(direction_0_file),status='old',form='unformatted',iostat=ier)
        if ( myrank == 0) print *,'reading direction0:',trim(direction_0_file)
        if (ier /= 0 ) then
                print *, 'error opening:',trim(direction_0_file)
                call exit_mpi(myrank,'file not found')
        endif
        read(1001) direction_0(:,:,:,1:NSPEC)
        close(1001)

        open(1001,file=trim(gradient_0_file),status='old',form='unformatted',iostat=ier)
        if (myrank == 0) print *, 'reading gradient0:',trim(gradient_0_file)
        if ( ier /= 0) then
                print *, 'error opening:',trim(gradient_0_file)
                call exit_mpi(myrank,'file not found')
        endif
        read(1001) gradient_0(:,:,:,1:NSPEC)
        close(1001)

        open(1001,file=trim(gradient_1_file),status='old',form='unformatted',iostat=ier)
        if (myrank == 0) print *, 'reading gradient1:',trim(gradient_1_file)
        if (ier /= 0) then
                print *, 'error opening:',trim(gradient_1_file)
                call exit_mpi(myrank,'file not found')
        endif
        read(1001) gradient_1(:,:,:,1:NSPEC)
        close(1001)




        do ispec=1,NSPEC
           do k = 1,NGLLZ
              do j = 1,NGLLY
                 do i = 1,NGLLX

                    direction_1(i,j,k,ispec)=-gradient_1(i,j,k,ispec)+beta*direction_0(i,j,k,ispec)
                 enddo ! i
              enddo  ! j
            enddo ! k
        enddo ! ispec
        open(1001,file=trim(direction_1_file),form='unformatted',action='write')
        if (myrank == 0) print *, 'writing direction1:',direction_1_file
        write(1001) direction_1
        close(1001)

  enddo ! kernel type

  call MPI_FINALIZE(ier)

end program xcompute_direction_cg
