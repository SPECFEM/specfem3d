program random_model

  implicit none

  include "mpif.h"

  integer,parameter :: NGLLX=5,NGLLY=5,NGLLZ=5,IOUT=20
  character(len=512),parameter :: LOCAL_PATH='../in_out_files/DATABASES_MPI/'

  integer :: myrank,ier,nspec,nglob,NPROC,j,ios
  double precision :: percent
  character(len=256) prname,arg
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read,random


  call MPI_Init(ier)
  call MPI_Comm_Rank(MPI_COMM_WORLD,myrank,ier)
  call MPI_Comm_Size(MPI_COMM_WORLD,NPROC,ier)

  !! input parameters
  if( iargc() .ne. 1 ) stop 'Usage: ./xrandom_model percent [must be small enough (~1d-5) for F*dm=S(m+dm)-S(m) to be valid]'
  j=1;  call getarg(j, arg); read(arg,*,iostat=ios) percent;   if (ios /= 0) stop 'Error reading percent'

  ! processors name
  write(prname,'(a,i6.6,a)') trim(LOCAL_PATH)//'proc',myrank,'_'
  ! nspec & nglob
  open(unit=IOUT,file=trim(adjustl(prname))//'external_mesh.bin',status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening database proc######_external_mesh.bin'
  read(IOUT) nspec
  read(IOUT) nglob
  close(IOUT)
  ! rho
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array rho_read'
  open(unit=IOUT,file=trim(adjustl(prname))//'rho.bin',status='old',action='read',form='unformatted')
  read(IOUT) rho_read
  close(IOUT)
  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vp_read'
  open(unit=IOUT,file=trim(adjustl(prname))//'vp.bin',status='old',action='read',form='unformatted')
  read(IOUT) vp_read
  close(IOUT)
  ! vs
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vs_read'
  open(unit=IOUT,file=trim(adjustl(prname))//'vs.bin',status='old',action='read',form='unformatted')
  read(IOUT) vs_read
  close(IOUT)

  ! perturb model randomly
  allocate( random(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array random'
  CALL RANDOM_SEED()
  
  CALL RANDOM_NUMBER(random)
  random=random/maxval(abs(random))*2.0-1.0
  rho_read=rho_read*(1.0+percent*random)

  CALL RANDOM_NUMBER(random)
  random=random/maxval(abs(random))*2.0-1.0
  vp_read= vp_read*(1.0+percent*random)

  CALL RANDOM_NUMBER(random)
  random=random/maxval(abs(random))*2.0-1.0
  vs_read= vs_read*(1.0+percent*random)

  ! store perturbed model
  ! rho
  open(unit=IOUT,file=trim(adjustl(prname))//'rho.bin',status='old',action='write',form='unformatted')
  write(IOUT) rho_read
  close(IOUT)
  ! vp
  open(unit=IOUT,file=trim(adjustl(prname))//'vp.bin',status='old',action='write',form='unformatted')
  write(IOUT) vp_read
  close(IOUT)
  ! vs
  open(unit=IOUT,file=trim(adjustl(prname))//'vs.bin',status='old',action='write',form='unformatted')
  write(IOUT) vs_read
  close(IOUT)

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  call MPI_Finalize(ier)

end program random_model
