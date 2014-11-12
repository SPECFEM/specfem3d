! This subroutine is used to sum all event kernels to get misfit kernels 
! Last modified: Fri Sep 14 08:50:36 EDT 2012

program sum_kernels
  implicit none 
  include 'mpif.h' 
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h' 

  integer,parameter:: NSPEC=NSPEC_CRUST_MANTLE 
  integer,parameter:: NKERNEL=6    !bulk_betah, bulk_betav, bulk_c, eta, rho, hess 

  integer:: myrank, sizeprocs,ier 
  integer:: ios,nevent,ievent,iker 
  character(len=350):: eventid,line,kernel_file,input_dir,output_dir 
  character(len=150):: event_list(1000), kernel_name(NKERNEL) 
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: kernel,total_kernel 


  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier) 

  call getarg(1,input_dir) 
  call getarg(2,output_dir)
  call getarg(3,eventid) 


  if (trim(input_dir) == '' & 
     .or. trim(output_dir) == '' & 
     .or. trim(eventid) == '' ) then 
     call exit_MPI(myrank,'USAGE: xsum_kernels input_dir output_dir eventid')
  end if 

  if (myrank == 0) then 
     write(*,*) 'SUM EVENT KERNELS TO GET MISFIT KENRELS'
     write(*,*) 'INPUT DIRECTION:',input_dir 
     write(*,*) 'OUTPUT DIRECTION:',output_dir 
     write(*,*) 'INPUT EVENTFILE:',eventid 
  end if 


  kernel_name=(/"reg1_bulk_betah_kernel","reg1_bulk_betav_kernel","reg1_bulk_c_kernel","reg1_eta_kernel","reg1_rho_kernel","reg1_hess_kernel"/) 

  nevent=0 
  open(unit=1001,file=trim(eventid),status='old',iostat=ios) 
  if ( ios /= 0 ) then 
     print*, 'ERROR OPENING', trim(eventid)
     stop 
  end if 
  do while ( 1 == 1) 
     read(1001,'(a)',iostat=ios) line 
     if ( ios /=0) exit
     nevent=nevent+1 
     event_list(nevent)=line 
  end do 
  close(1001) 


  do iker=1,NKERNEL 
     total_kernel=0.

     do ievent=1,nevent 

        if (myrank==0) write(*,*) 'READING IN EVENT KERNEL:',trim(kernel_name(iker)),' FOR ',trim(event_list(ievent)) 

        write(kernel_file,'(a,i6.6,a)') trim(input_dir)//'/'//trim(event_list(ievent))//'/KERNEL/'//'proc',myrank,'_'//trim(kernel_name(iker))//'.bin' 
                
        open(unit=1002,file=trim(kernel_file),status='old',form='unformatted')
        read(1002) kernel(:,:,:,1:NSPEC) 
        close(1002)

        ! sum ther kernel
        if (iker == 6 ) then ! for hessian , sum the absolute value 
           total_kernel(:,:,:,1:NSPEC)=total_kernel(:,:,:,1:NSPEC) + abs(kernel(:,:,:,1:NSPEC))
        else 
           total_kernel(:,:,:,1:NSPEC)=total_kernel(:,:,:,1:NSPEC) + kernel(:,:,:,1:NSPEC)
        end if 

     end do 
     if (myrank==0) write(*,*) 'WRITING MISFIT KERNELS:',trim(kernel_name(iker))
     write(kernel_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name(iker))//'.bin'
        
     open(1002,file=trim(kernel_file),form='unformatted')
     write(1002) total_kernel(:,:,:,1:NSPEC)
     close(1002) 
  end do

  if (myrank==0) write(*,*) 'done summing all the kernels' 

  call MPI_FINALIZE(ier) 

end program sum_kernels
