module dsm_coupling
  use specfem_par
  implicit none
  include 'mpif.h'
  include 'precision.h'

  integer  nb_DSM_sources
  character(len=250), allocatable :: dir_dsm_name(:)
  real(kind=CUSTOM_REAL) , allocatable :: sources_position(:,:),incident_long_wave(:),&
       incident_direction(:,:),incident_position(:,:),incident_params(:,:)
  contains

    subroutine read_DSM_coupling_parameters()
      implicit none
      integer i,ier

      if (myrank==0) then
         open(27,file='inversion.par')
         read(27,*) nb_DSM_sources
         close(27)
         if (COUPLING_WITH_DSM) then
            open(27,file='coupling_parameters.par',status='old')
            !read(27,*) nb_DSM_sources
            allocate(dir_dsm_name(nb_DSM_sources))
            do i=1,nb_DSM_sources
               read(27,'(a)') dir_dsm_name(i)
               write(*,*) dir_dsm_name(i)
            enddo
            close(27)
         else if (COUPLING_WITH_PLAN) then
            open(27,file='plan_parameters.par',status='old')
            !read(27,*) nb_DSM_sources
            allocate(incident_long_wave(nb_DSM_sources))
            allocate(incident_direction(3,nb_DSM_sources))
            allocate(incident_position(3,nb_DSM_sources))
            allocate(incident_params(3,nb_DSM_sources))
            do i=1,nb_DSM_sources
               read(27,*) incident_long_wave(i)
               read(27,*) incident_direction(1,i),incident_direction(2,i),incident_direction(3,i)
               read(27,*) incident_position(1,i),incident_position(2,i),incident_position(3,i)
               read(27,*) incident_params(1,i),incident_params(2,i),incident_params(3,i)
            enddo
            close(27)
         else
            open(27,file='source_parameters.par',status='old')
            !read(27,*) nb_DSM_sources
            allocate(sources_position(3,nb_DSM_sources))
            do i=1,nb_DSM_sources
               read(27,*) sources_position(1,i),sources_position(2,i),sources_position(3,i)
            enddo
         endif
      endif

      if (COUPLING_WITH_DSM) then
         call mpi_bcast(nb_DSM_sources,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
         if (myrank > 0) allocate(dir_dsm_name(nb_DSM_sources))
         if (myrank==0) write(*,*) 'Bcst DSM dir',dir_dsm_name
         call mpi_bcast( dir_dsm_name, 250*nb_DSM_sources,MPI_CHARACTER ,0,MPI_COMM_WORLD,ier)
      else if (COUPLING_WITH_PLAN) then
         call mpi_bcast(nb_DSM_sources,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
         if (myrank > 0) then
            allocate(incident_long_wave(nb_DSM_sources))
            allocate(incident_direction(3,nb_DSM_sources))
            allocate(incident_position(3,nb_DSM_sources))
            allocate(incident_params(3,nb_DSM_sources))
         endif

         call mpi_bcast(incident_long_wave,nb_DSM_sources,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
         call mpi_bcast(incident_direction,3*nb_DSM_sources,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
         call mpi_bcast(incident_position,3*nb_DSM_sources,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
         call mpi_bcast(incident_params,3*nb_DSM_sources,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      else
         call mpi_bcast(nb_DSM_sources,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
         if (myrank > 0) allocate(sources_position(3,nb_DSM_sources))
         call mpi_bcast(sources_position,3*nb_DSM_sources,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
      endif
    end subroutine read_DSM_coupling_parameters

    subroutine setup_DSM_tractions()
      implicit none


      call read_DSM_coupling_parameters()
      allocate(Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces))
      allocate(Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces))
      if (myrank==0) write(*,*) 'DSM sources :',nb_DSM_sources
    end subroutine setup_DSM_tractions


    subroutine open_DSM_tractions(i)

      integer i
      integer ier

      !write(*,*) 'open :' ,trim(dir_dsm_name(i))//'/'//trim(dsmname)//'vel.bin'
      open(unit=IIN_veloc_dsm,file=trim(dir_dsm_name(i))//'/'//trim(dsmname)//'vel.bin',status='old',&
           action='read',form='unformatted',iostat=ier)

      open(unit=IIN_tract_dsm,file=trim(dir_dsm_name(i))//'/'//trim(dsmname)//'tract.bin',status='old',&
           action='read',form='unformatted',iostat=ier)

    end subroutine open_DSM_tractions

    subroutine close_DSM_tractions(i)
      integer i
      close(IIN_veloc_dsm)
      close(IIN_tract_dsm)

    end subroutine close_DSM_tractions

end module dsm_coupling
