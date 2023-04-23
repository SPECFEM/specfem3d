program read_absorbing_interfaces

  implicit none
  include 'mpif.h'
  include 'constants.h'
  integer NSPEC_AB,NGLOB_AB
  integer NSPEC,NGLOB
  integer ispec,iglob,nbproc,myrank,inodes,nbnode,ier,ierr,i,j,k
  integer iface,num_abs_boundary_faces,igll
  character(len=27) filename


! name of the database file
  character(len=256) prname,dsmname,trname,LOCAL_PATH,Modele1D,datafile,tablefile,tracfile,TRAC_PATH
  character(len=256) dsmdirectory,meshdirectory
  integer, dimension(:,:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore
  real(kind=SIZE_DOUBLE), dimension (:,:,:,:), allocatable :: xstore_ref,ystore_ref,zstore_ref
  real(kind=CUSTOM_REAL) x0,y0,z0,x1,y1,z1
  !
  real(kind=CUSTOM_REAL), allocatable :: abs_boundary_normal(:,:,:),abs_boundary_jacobian2Dw(:,:)
  integer, allocatable  :: abs_boundary_ijk(:,:,:),abs_boundary_ispec(:)
  logical, allocatable  :: ispec_is_inner(:), ispec_is_elastic(:)
  real(kind=CUSTOM_REAL) x,y,z,nx,ny,nz

!
  double precision r,theta,phi,lon0,lat0,azi0,xd,yd,zd,zmin,xx(3),xlat,ylon
  double precision, allocatable :: vp(:,:),vs(:,:),rho(:,:),zlayer(:)
  integer ndeg,nlayer,ilayer
  integer code_face,iref,jref,kref
!
  integer, allocatable :: flag_boundary_xmin(:),flag_boundary_xmax(:),flag_boundary_ymin(:)&
  ,flag_boundary_ymax(:),flag_boundary_zmin(:)
  integer n2d_xmin,n2d_xmax,n2d_ymin,n2d_ymax,n2d_zmin,i_code_face,ispec2D
  integer max_spec,ispec_glob,max_spec_local,iproc,ispec_global
  integer, allocatable :: IndLoc2Glob(:,:),IndSpec_face(:,:)
!
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NTIMESTEP
  integer ibloc,nbloc,it_local,i1,i2,i3,i4,iii,iblocmin
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable  :: TracXmin,TracXmax,TracYmin,TracYmax,TracZmin
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable  :: VelXmin,VelXmax,VelYmin,VelYmax,VelZmin
  double precision, dimension(:,:,:,:), allocatable :: TracsgXmin,Tracsgxmax,Tracsgymin,Tracsgymax,Tracsgzmin
  double precision, dimension(:,:,:,:), allocatable :: VelsgXmin,Velsgxmax,Velsgymin,Velsgymax,Velsgzmin
  integer, dimension(:,:,:,:), allocatable :: igxmin,igxmax,igymin,igymax,igzmin
  integer nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin
  integer nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax
  integer nbrcymin,Nblocymin,Nbrecymin,nseisymin
  integer nbrcymax,Nblocymax,Nbrecymax,nseisymax
  integer kkkk

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)

!
!------------------- PARAMETRES ----------------------------
  NTIMESTEP = 100
  open(10,file='ParFileInterface')
  read(10,'(a)') dsmdirectory
  read(10,'(a)') meshdirectory
  read(10,'(a)') LOCAL_PATH
  read(10,'(a)') TRAC_PATH
  read(10,*) iblocmin
  close(10)
  Modele1D=trim(meshdirectory)//"/model_1D.in"
  dsmname=dsmdirectory
  1000 format(6f30.10,i5)
  1001 format(3f30.10,3i15,i5)

  write(*,*) myrank,trim(dsmdirectory)
  write(*,*) myrank,trim(meshdirectory)
  write(*,*) myrank,trim(LOCAL_PATH)
  write(*,*) myrank,trim(TRAC_PATH)

! --- lecture du modele 1D
  open(10,file=Modele1D)
  read(10,*) nlayer,ndeg

  allocate(vp(nlayer,ndeg))
  allocate(vs(nlayer,ndeg))
  allocate(rho(nlayer,ndeg))
  allocate(zlayer(nlayer))

  do ilayer = 1, nlayer

     read(10,*) zlayer(ilayer)
     read(10,*) vp(ilayer,:)
     read(10,*) vs(ilayer,:)
     read(10,*) rho(ilayer,:)

  enddo
  read(10,*) zmin
  read(10,*) lon0,lat0,azi0
  close(10)

!----------------------- lecture des tables de correspondances ---------

  open(90,file=trim(meshdirectory)//'/flags_boundary.txt')  ! table de correspondance ispec2D < - > ispec

  open(91,file=trim(meshdirectory)//'/Nb_ielm_faces.txt')   ! taille tableaux pour allocation memoire
  read(91,*) n2d_xmin
  read(91,*) n2d_xmax
  read(91,*) n2d_ymin
  read(91,*) n2d_ymax
  read(91,*) n2d_zmin
  close(91)
  nspec2d_xmin=n2d_xmin
  nspec2d_xmax=n2d_xmax
  nspec2d_ymin=n2d_ymin
  nspec2d_ymax=n2d_ymax
  nspec2d_bottom=n2d_zmin

  allocate(flag_boundary_xmin(n2d_xmin))
  allocate(flag_boundary_xmax(n2d_xmax))
  allocate(flag_boundary_ymin(n2d_ymin))
  allocate(flag_boundary_ymax(n2d_ymax))
  allocate(flag_boundary_zmin(n2d_zmin))

  max_spec = 0
  do while (1 == 1)
    read(90,*,end=100) ispec,ispec2D,code_face
    if (code_face == 1) flag_boundary_xmin(ispec2D)=ispec
    if (code_face == 2) flag_boundary_xmax(ispec2D)=ispec
    if (code_face == 3) flag_boundary_ymin(ispec2D)=ispec
    if (code_face == 4) flag_boundary_ymax(ispec2D)=ispec
    if (code_face == 5) flag_boundary_zmin(ispec2D)=ispec
    max_spec = max(max_spec,ispec)
  enddo
  100 close(90)
  allocate(IndSpec_face(max_spec,5))
  IndSpec_face(:,:) = 0
  open(90,file=trim(meshdirectory)//'/flags_boundary.txt')  ! table de correspondance ispec2D < - > ispec
  do while (1 == 1)
    read(90,*,end=101) ispec,ispec2D,code_face
    IndSpec_face(ispec,code_face)=ispec2D
  enddo
  101 close(90)

  open(90,file=trim(meshdirectory)//'/Numglob2loc_elmn.txt')
  max_spec_local=0
  do while(1 == 1)
     read(90,*,end=102) ispec_glob,ispec,iproc
     max_spec_local=max(max_spec_local,ispec)
  enddo
  102  close(90)
  allocate(IndLoc2Glob(max_spec_local,nbproc))
  IndLoc2Glob(:,:) = 0
  open(90,file=trim(meshdirectory)//'/Numglob2loc_elmn.txt')
  max_spec_local=0
  do while(1 == 1)
     read(90,*,end=103) ispec_glob,ispec,iproc
     IndLoc2Glob(ispec,iproc+1) = ispec_glob
  enddo
  103  close(90)

  ! tracctions et vitesse  VM
  write(*,*) 'NtimeStep', NtimeStep
  allocate(TracXmin(3,NtimeStep,NGLLY,NGLLZ,nspec2D_xmin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(TracXmax(3,NtimeStep,NGLLY,NGLLZ,nspec2D_xmax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(TracYmin(3,NtimeStep,NGLLX,NGLLZ,nspec2D_ymin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(TracYmax(3,NtimeStep,NGLLX,NGLLZ,nspec2D_ymax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(TracZmin(3,NtimeStep,NGLLX,NGLLY,nspec2D_bottom),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  ! velocity
  allocate(VelXmin(3,NtimeStep,NGLLY,NGLLZ,nspec2D_xmin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(VelXmax(3,NtimeStep,NGLLY,NGLLZ,nspec2D_xmax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(VelYmin(3,NtimeStep,NGLLX,NGLLZ,nspec2D_ymin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(VelYmax(3,NtimeStep,NGLLX,NGLLZ,nspec2D_ymax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(VelZmin(3,NtimeStep,NGLLX,NGLLY,nspec2D_bottom),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'

! premiere ouverutre des fichiers tractions et vitesses et table de correspondance -----------

!--------------
  if (myrank == 0) then
     open(unit=41,file=dsmname(1:len_trim(dsmname))//'velxmin.bin',form='unformatted',action='read')
     read(41) nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin
     open(unit=31,file=dsmname(1:len_trim(dsmname))//'tractxmin.bin',status='old',form='unformatted',action='read')
     read(31) nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin
  endif
  call MPI_Bcast(nbrcxmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nblocxmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nbrecxmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nseisxmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !write(41) nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin
  !write(*,*) nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin
  allocate(igxmin(3,NGLLY,NGLLZ,nspec2D_xmin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Tracsgxmin(NtimeStep,3,nseisxmin,Nbrecxmin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
   allocate(Velsgxmin(NtimeStep,3,nseisxmin,Nbrecxmin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  if (myrank == 0) then
     open(27,file=meshdirectory(1:len_trim(meshdirectory))//'IgXmin')
     call ReadIg(igxmin,NGLLY,NGLLZ,nspec2D_xmin,Nbrecxmin)
     close(27)
  endif
  call MPI_Bcast(igxmin,3*NGLLY*NGLLZ*nspec2D_xmin,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (myrank == 0) then
     open(unit=42,file=dsmname(1:len_trim(dsmname))//'velxmax.bin',form='unformatted',action='read')
     read(42) nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax
     open(unit=32,file=dsmname(1:len_trim(dsmname))//'tractxmax.bin',status='old',form='unformatted',action='read')
     read(32) nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax
  endif
  call MPI_Bcast(nbrcxmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nblocxmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nbrecxmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nseisxmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !write(42) nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax
  !write(*,*) 'Tailles :',nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax
  allocate(igxmax(3,NGLLY,NGLLZ,nspec2D_xmax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Tracsgxmax(NtimeStep,3,nseisxmax,Nbrecxmax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Velsgxmax(NtimeStep,3,nseisxmax,Nbrecxmax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  if (myrank == 0) then
     open(27,file=meshdirectory(1:len_trim(meshdirectory))//'IgXmax')
     call ReadIg(igxmax,NGLLY,NGLLZ,nspec2D_xmax,Nbrecxmax)
     close(27)
  endif
  call MPI_Bcast(igxmax,3*NGLLY*NGLLZ*nspec2D_xmax,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (myrank == 0) then
     open(unit=43,file=dsmname(1:len_trim(dsmname))//'velymin.bin',form='unformatted')
     read(43) nbrcymin,Nblocymin,Nbrecymin,nseisymin
     open(unit=33,file=dsmname(1:len_trim(dsmname))//'tractymin.bin',status='old',form='unformatted')
     read(33) nbrcymin,Nblocymin,Nbrecymin,nseisymin
  endif
  call MPI_Bcast(nbrcymin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nblocymin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nbrecymin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nseisymin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !write(43) nbrcymin,Nblocymin,Nbrecymin,nseisymin
  allocate(igymin(3,NGLLX,NGLLZ,nspec2D_ymin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Tracsgymin(NtimeStep,3,nseisymin,Nbrecymin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Velsgymin(NtimeStep,3,nseisymin,Nbrecymin),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  if (myrank == 0) then
     open(27,file=meshdirectory(1:len_trim(meshdirectory))//'IgYmin')
     call ReadIg(igymin,NGLLX,NGLLZ,nspec2D_ymin,Nbrecymin)
     close(27)
  endif
  call MPI_Bcast(igymin,3*NGLLX*NGLLZ*nspec2D_ymin,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (myrank == 0) then
     open(unit=44,file=dsmname(1:len_trim(dsmname))//'velymax.bin',form='unformatted',action='read')
     read(44) nbrcymax,Nblocymax,Nbrecymax,nseisymax
     open(unit=34,file=dsmname(1:len_trim(dsmname))//'tractymax.bin',status='old',form='unformatted',action='read')
     read(34) nbrcymax,Nblocymax,Nbrecymax,nseisymax
  endif
  call MPI_Bcast(nbrcymax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nblocymax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nbrecymax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nseisymax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !write(44) nbrcymax,Nblocymax,Nbrecymax,nseisymax
  allocate(igymax(3,NGLLX,NGLLZ,nspec2D_ymax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Tracsgymax(NtimeStep,3,nseisymax,Nbrecymax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Velsgymax(NtimeStep,3,nseisymax,Nbrecymax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  if (myrank == 0) then
     open(27,file=meshdirectory(1:len_trim(meshdirectory))//'IgYmax')
     call ReadIg(igymax,NGLLX,NGLLZ,nspec2D_ymax,Nbrecymax)
     close(27)
  endif
  call MPI_Bcast(igymax,3*NGLLX*NGLLZ*nspec2D_ymax,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (myrank == 0) then
     open(unit=45,file=dsmname(1:len_trim(dsmname))//'velzmin.bin',form='unformatted')
     read(45) i1,i2,i3,i4
     open(unit=35,file=dsmname(1:len_trim(dsmname))//'tractzmin.bin',status='old',form='unformatted')
     read(35) i1,i2,i3,i4
  endif
  !write(45) i1,i2,i3,i4
  allocate(igzmin(3,NGLLX,NGLLY,nspec2D_bottom))
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Tracsgzmin(NtimeStep,3,nseisxmin,nseisymin))
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  allocate(Velsgzmin(NtimeStep,3,nseisxmin,nseisymin))
  if (ier /= 0) stop 'error: not enough memory to allocate array'
  if (myrank == 0) then
     open(27,file=meshdirectory(1:len_trim(meshdirectory))//'IgZmin')
     call ReadIgZ(igzmin,NGLLX,NGLLY,nspec2D_bottom)
     close(27)
  endif
  nbloc = Nblocymax
  call MPI_Bcast(igzmin,3*NGLLX*NGLLY*nspec2D_bottom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
!---------------------- lecture du maillage
!
!  ouverture d'un fichier traction par thread MPI
!
!
     ! lecture des facettes pour un processeur
     write(*,*) LOCAL_PATH
     call create_name_database(prname,myrank,LOCAL_PATH)
     call create_name_database(trname,myrank,TRAC_PATH)
     write(*,*) LOCAL_PATH
     write(*,*) prname(1:len_trim(prname))//'external_mesh.bin'

     open(27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old', &
          action='read',form='unformatted',iostat=ier)
     read(27) NSPEC_AB
     read(27) NGLOB_AB

     allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
     allocate(xstore(NGLOB_AB), &
          ystore(NGLOB_AB), &
          zstore(NGLOB_AB),stat=ier)
     allocate(ispec_is_inner(nspec_ab))
     allocate(ispec_is_elastic(nspec_ab))

     read(27) ibool
     read(27) xstore
     read(27) ystore
     read(27) zstore
     close(27)

     write(*,*) 'STEP 1'

     write(*,*) prname(1:len_trim(prname))//'inner'
     open(27,file=prname(1:len_trim(prname))//'inner',status='old', &
          action='read',form='unformatted',iostat=ier)
     read(27) ispec_is_inner
     read(27) ispec_is_elastic
     close(27)

     write(*,*) 'STEP 2'

     write(*,*) prname(1:len_trim(prname))//'absorb_dsm'
     open(27,file=prname(1:len_trim(prname))//'absorb_dsm',status='old', &
          action='read',form='unformatted',iostat=ier)
     read(27) num_abs_boundary_faces
     write(*,*) num_abs_boundary_faces
     allocate(abs_boundary_ispec(num_abs_boundary_faces))
     allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces))
     allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces))
     allocate(abs_boundary_normal(3,NGLLSQUARE,num_abs_boundary_faces))

     read(27) abs_boundary_ispec
     read(27) abs_boundary_ijk
     read(27) abs_boundary_jacobian2Dw
     read(27) abs_boundary_normal
     close(27)
     ! fin de lecture

     write(*,*) 'STEP 3'

     ! ouverture du fichier traction
     !open(27,file=trname(1:len_trim(trname))//'tract.indx')
     open(28,file=trname(1:len_trim(trname))//'tract.bin',form='unformatted',action='write')
     open(29,file=trname(1:len_trim(trname))//'vel.bin',form='unformatted',action='write')


     do ibloc = 1 , nbloc
       if (myrank == 0) then

          write(*,*) 'read bloc', ibloc,nbloc

       !endif
          ! read traction
          call ReadTract(31,Tracsgxmin,TracXmin,NtimeStep,NGLLY,NGLLZ,nspec2D_xmin,nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin,igxmin)
          call ReadTract(32,Tracsgxmax,TracXmax,NtimeStep,NGLLY,NGLLZ,nspec2D_xmax,nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax,igxmax)
          call ReadTract(33,Tracsgymin,TracYmin,NtimeStep,NGLLX,NGLLZ,nspec2D_ymin,nbrcymin,Nblocymin,Nbrecymin,nseisymin,igymin)
          call ReadTract(34,Tracsgymax,TracYmax,NtimeStep,NGLLX,NGLLZ,nspec2D_ymax,nbrcymax,Nblocymax,Nbrecymax,nseisymax,igymax)
          call ReadTract(35,Tracsgzmin,TracZmin,NtimeStep,NGLLX,NGLLY,nspec2D_bottom,nbrcxmin,Nblocxmin,nseisxmin,nseisymin,igzmin)
          ! pour zmin : on lit d'abord nseisymin points et ensuite nseisxmin (ceci est l'ordre
          ! de la boucle de lecture)
!
          ! read velocity
          call ReadVel(41,Velsgxmin,VelXmin,NtimeStep,NGLLY,NGLLZ,nspec2D_xmin,nbrcxmin,Nblocxmin,Nbrecxmin,nseisxmin,igxmin)
          call ReadVel(42,Velsgxmax,VelXmax,NtimeStep,NGLLY,NGLLZ,nspec2D_xmax,nbrcxmax,Nblocxmax,Nbrecxmax,nseisxmax,igxmax)
          call ReadVel(43,Velsgymin,VelYmin,NtimeStep,NGLLX,NGLLZ,nspec2D_ymin,nbrcymin,Nblocymin,Nbrecymin,nseisymin,igymin)
          call ReadVel(44,Velsgymax,VelYmax,NtimeStep,NGLLX,NGLLZ,nspec2D_ymax,nbrcymax,Nblocymax,Nbrecymax,nseisymax,igymax)
          call ReadVel(45,Velsgzmin,VelZmin,NtimeStep,NGLLX,NGLLY,nspec2D_bottom,nbrcxmin,Nblocxmin,nseisxmin,nseisymin,igzmin)
          ! pour zmin : on lit d'abord nseisymin points et ensuite nseisxmin (ceci est l'ordre
          ! de la boucle de lecture)

       endif

       call MPI_Bcast(TracXmin,3*NtimeStep*NGLLY*NGLLZ*nspec2D_xmin   ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(TracXmax,3*NtimeStep*NGLLY*NGLLZ*nspec2D_xmax   ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(TracYmin,3*NtimeStep*NGLLX*NGLLZ*nspec2D_ymin   ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(TracYmax,3*NtimeStep*NGLLX*NGLLZ*nspec2D_ymax   ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(TracZmin,3*NtimeStep*NGLLX*NGLLY*nspec2D_bottom ,MPI_REAL,0,MPI_COMM_WORLD,ierr)

       call MPI_Bcast(VelXmin,3*NtimeStep*NGLLY*NGLLZ*nspec2D_xmin  ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(VelXmax,3*NtimeStep*NGLLY*NGLLZ*nspec2D_xmax  ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(VelYmin,3*NtimeStep*NGLLX*NGLLZ*nspec2D_ymin  ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(VelYmax,3*NtimeStep*NGLLX*NGLLZ*nspec2D_ymax  ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(VelZmin,3*NtimeStep*NGLLX*NGLLY*nspec2D_bottom,MPI_REAL,0,MPI_COMM_WORLD,ierr)

       do iface = 1,num_abs_boundary_faces ! boucle sur les facettes

          ispec = abs_boundary_ispec(iface)  ! numerotation locale
          ispec_global = IndLoc2Glob(ispec,myrank+1)
          if (ispec_global == 0) then
            write(*,*) 'element',ispec,myrank
            write(*,*) 'pas defini'
            stop
          endif
          ! test pour connaitre la face sur laquelle on se trouve
          igll = 13
          iref = abs_boundary_ijk(1,igll,iface)
          jref = abs_boundary_ijk(2,igll,iface)
          kref = abs_boundary_ijk(3,igll,iface)

          if (iref == 1) code_face = 1
          if (iref == NGLLX) code_face = 2
          if (jref == 1) code_face = 3
          if (jref == NGLLY) code_face = 4
          if (kref == 1) code_face=5

          ispec2d = IndSpec_face(ispec_global,code_face)

          ! xmin
          if (code_face == 1) then
             do igll=1,NGLLSQUARE
                i = abs_boundary_ijk(1,igll,iface)
                j = abs_boundary_ijk(2,igll,iface)
                k = abs_boundary_ijk(3,igll,iface)

                iglob = ibool(i,j,k,ispec)
                !x0=xstore_ref(i,j,k,ispec_global)
                !y0=ystore_ref(i,j,k,ispec_global)
                !z0=zstore_ref(i,j,k,ispec_global)

                x1=xstore(iglob)
                y1=ystore(iglob)
                z1=zstore(iglob)

                if (i /= 1) then
                  write(*,*) 'le point',i,j,k,iface
                  write(*,*) 'n est pas sur la face xmin'
                  stop
               endif
               !write(27,'(7i10)') i,j,k,iface,ispec2d,Igxmin(1,j,k,ispec2d), Igxmin(2,j,k,ispec2d) ! ispec2d,code_face
               !write(27,'(3f20.5,10x,3f20.5)') x0,y0,z0,x1,y1,z1
               if (ibloc >= iblocmin) write(28) TracXmin(:,:,j,k,ispec2d)
               if (ibloc >= iblocmin) write(29) VelXmin(:,:,j,k,ispec2d)
             enddo
          endif


          ! xmax
          if (code_face == 2) then
             do igll=1,NGLLSQUARE
               i = abs_boundary_ijk(1,igll,iface)
               j = abs_boundary_ijk(2,igll,iface)
               k = abs_boundary_ijk(3,igll,iface)

               iglob = ibool(i,j,k,ispec)
               !x0=xstore_ref(i,j,k,ispec_global)
               !y0=ystore_ref(i,j,k,ispec_global)
               !z0=zstore_ref(i,j,k,ispec_global)

               x1=xstore(iglob)
               y1=ystore(iglob)
               z1=zstore(iglob)



               if (i /= NGLLX) then
                 write(*,*) 'le point',i,j,k,iface
                 write(*,*) 'n est pas sur la face xmax'
                 stop
               endif
               !write(27,'(7i10)') i,j,k,iface,ispec2d,Igxmax(1,j,k,ispec2d), Igxmax(2,j,k,ispec2d) ! ispec2d,code_face
               if (ibloc >= iblocmin) write(28) TracXmax(:,:,j,k,ispec2d)
               if (ibloc >= iblocmin) write(29) VelXmax(:,:,j,k,ispec2d)
            enddo
         endif


         ! ymin
         if (code_face == 3) then
            do igll=1,NGLLSQUARE
               i = abs_boundary_ijk(1,igll,iface)
               j = abs_boundary_ijk(2,igll,iface)
               k = abs_boundary_ijk(3,igll,iface)

               iglob = ibool(i,j,k,ispec)
               !x0=xstore_ref(i,j,k,ispec_global)
               !y0=ystore_ref(i,j,k,ispec_global)
               !z0=zstore_ref(i,j,k,ispec_global)

               x1=xstore(iglob)
               y1=ystore(iglob)
               z1=zstore(iglob)


               if (j /= 1) then
                 write(*,*) 'le point',i,j,k,iface
                 write(*,*) 'n est pas sur la face ymin'
                 stop
               endif
               !write(27,'(7i10)') i,j,k,iface,ispec2d,Igymin(1,i,k,ispec2d), Igymin(2,i,k,ispec2d) ! ispec2d,code_face
               !write(27,'(3f20.5,10x,3f20.5)') x0,y0,z0,x1,y1,z1
               if (ibloc >= iblocmin) write(28) TracYmin(:,:,i,k,ispec2d)
               if (ibloc >= iblocmin) write(29) VelYmin(:,:,i,k,ispec2d)
               !write(27,*) 'YMIN:',ispec2d,i,j,k
               !write(27,*) 'PT :',x1,y1,z1
               !do kkkk=1,50
               !write(27,*) TracYmin(1,kkkk,i,k,ispec2d),TracYmin(2,kkkk,i,k,ispec2d), TracYmin(3,kkkk,i,k,ispec2d)
               !enddo
               !write(27,*) '-------------------------------'
            enddo
         endif

         ! ymax
         if (code_face == 4) then
           do igll=1,NGLLSQUARE
              i = abs_boundary_ijk(1,igll,iface)
              j = abs_boundary_ijk(2,igll,iface)
              k = abs_boundary_ijk(3,igll,iface)

              iglob = ibool(i,j,k,ispec)
              !x0=xstore_ref(i,j,k,ispec_global)
              !y0=ystore_ref(i,j,k,ispec_global)
              !z0=zstore_ref(i,j,k,ispec_global)

              x1=xstore(iglob)
              y1=ystore(iglob)
              z1=zstore(iglob)


              if (j /= NGLLY) then
                write(*,*) 'le point',i,j,k,iface
                write(*,*) 'n est pas sur la face ymax'
                stop
              endif
              !write(27,'(7i10)') i,j,k,iface,ispec2d,Igymax(1,i,k,ispec2d), Igymax(2,i,k,ispec2d) ! ispec2d,code_face
              !write(27,'(3f20.5,10x,3f20.5)') x0,y0,z0,x1,y1,z1
              if (ibloc >= iblocmin) write(28) TracYmax(:,:,i,k,ispec2d)
              if (ibloc >= iblocmin) write(29) VelYmax(:,:,i,k,ispec2d)
              !write(27,*) 'YMAX:',ispec2d,i,j,k
              !write(27,*) 'PT :',x1,y1,z1
              !do kkkk=1,50
              !write(27,*)TracYmax(1,kkkk,i,k,ispec2d),TracYmax(2,kkkk,i,k,ispec2d),TracYmax(3,kkkk,i,k,ispec2d)
              !enddo
              !write(27,*) '-------------------------------'
           enddo
         endif

         ! zmin
         if (code_face == 5) then
            do igll=1,NGLLSQUARE
              i = abs_boundary_ijk(1,igll,iface)
              j = abs_boundary_ijk(2,igll,iface)
              k = abs_boundary_ijk(3,igll,iface)

              iglob = ibool(i,j,k,ispec)
              !x0=xstore_ref(i,j,k,ispec_global)
              !y0=ystore_ref(i,j,k,ispec_global)
              !z0=zstore_ref(i,j,k,ispec_global)

              x1=xstore(iglob)
              y1=ystore(iglob)
              z1=zstore(iglob)



              if (k /= 1) then
                write(*,*) 'le point',i,j,k,iface
                write(*,*) 'n est pas sur la face zmin'
                stop
              endif
              !write(27,'(7i10)') i,j,k,iface,ispec2d,Igzmin(1,i,j,ispec2d), Igzmin(2,i,j,ispec2d) ! ispec2d,code_face
              !write(27,'(3f20.5,10x,3f20.5)') x0,y0,z0,x1,y1,z1
              if (ibloc >= iblocmin) write(28) TracZmin(:,:,i,j,ispec2d)
              if (ibloc >= iblocmin) write(29) VelZmin(:,:,i,j,ispec2d)
            enddo
         endif

     enddo
   enddo
     !close(27)  ! fermeture du fichier traction
     close(28)
     close(29)
     close(31)
     close(32)
     close(33)
     close(34)
     close(35)
     close(41)
     close(42)
     close(43)
     close(44)
     close(45)

  call MPI_FINALIZE(ierr)

!!!  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!----------- STEP read_absorbing fin ----------!!!!!!!!!!!'

end program read_absorbing_interfaces



!=====================================================================

subroutine create_name_database(prname,iproc,LOCAL_PATH)

! create the name of the database for the mesher and the solver

  implicit none

  integer iproc

! name of the database file
  character(len=256) prname,procname,LOCAL_PATH,clean_LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_')") iproc

! suppress white spaces if any
  clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full name with path
  prname = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // procname

end subroutine create_name_database


subroutine Cart2Sph(lat,lon,X)
  implicit none
  double precision lat,lon,X(3)
  double precision deg2rad

  deg2rad = 3.141592653589793d0/180.d0

  lat = 90.d0 - dacos(X(3)) / deg2rad
  lon = datan2(X(2),X(1)) / deg2rad

  !write(*,*) lat,lon

end subroutine Cart2Sph

subroutine Sph2Cart(lat,lon,X)
  implicit none
  double precision lat1,lon1,lat,lon,X(3)
  double precision deg2rad

  deg2rad = 3.141592653589793d0/180.d0


  lat1 = deg2rad * (90.d0-lat)
  lon1 = deg2rad * lon

  X(1)=dcos(lon1)*dsin(lat1)
  X(2)=dsin(lon1)*dsin(lat1)
  X(3)=dcos(lat1)
  !write(*,* ) 'll ', lat,lon
  !write(*,*) X

end subroutine Sph2Cart

subroutine Rotate(lat,lon,Olat,Olon)
  implicit none
  double precision Olon, Olat, lat, lon, deg2rad
  double precision Rot(3,3),X(3),X1(3)
  integer i,j,k

  deg2rad  =  3.141592653589793d0/180.d0

  Rot(1,1)=dcos(deg2rad*Olat)
  Rot(1,2)=0.d0
  Rot(1,3)=-dsin(deg2rad*Olat)

  Rot(2,1)=0.d0
  Rot(2,2)=1.d0
  Rot(2,3)=0.d0

  Rot(3,1)=dsin(deg2rad*Olat)
  Rot(3,2)=0.d0
  Rot(3,3)=dcos(deg2rad*Olat)

  call Sph2Cart(lat,lon,X)

  do i=1,3
     X1(i)=0.d0
  enddo

  do i=1,3
     do k=1,3
        X1(i)=X1(i) + Rot(i,k) * X(k)
     enddo
  enddo
  call Cart2Sph(lat,lon,X1)

  lon = lon + olon


end subroutine Rotate




subroutine ReadIg(ig,NGLL1,NGLL2,nspec2D,Nb)
    implicit none
    integer NGLL1,NGLL2,nspec2D,Nb
    integer i,i1,i2,i3,i4,i5,i6
    integer ig(3,NGLL1,NGLL2,nspec2D)
    do i = 1, NGLL1*NGLL2*nspec2D
       read(27,*) i1,i2,i3,i4,i5,i6
       ig(1,i1,i2,i3)=i4
       ig(2,i1,i2,i3)=1 + Nb - i5 ! il faut inverser car les profs sont rangees
                                  ! dans le sens oppose dans dsm / specfem3D
       ig(3,i1,i2,i3)=i6
       !write(*,*) ig(1,i1,i2,i3),ig(2,i1,i2,i3)
    enddo
end subroutine ReadIg

subroutine ReadIgZ(ig,NGLL1,NGLL2,nspec2D)
    implicit none
    integer NGLL1,NGLL2,nspec2D
    integer i,i1,i2,i3,i4,i5,i6
    integer ig(3,NGLL1,NGLL2,nspec2D)
    do i = 1, NGLL1*NGLL2*nspec2D
       read(27,*) i1,i2,i3,i4,i5,i6
       ig(1,i1,i2,i3)=i4   ! i4
       ig(2,i1,i2,i3)=i5   ! i5
       ig(3,i1,i2,i3)=i6
       !write(*,*) ig(1,i1,i2,i3),ig(2,i1,i2,i3)
    enddo
end subroutine ReadIgZ




  subroutine ReadVel(iunit,Trcsg,Trac,NtimeStep,NGLL1,NGLL2,nspec2D,nbrc,Nbloc,Nbrec,nseis,Ig)
    implicit none
    include 'constants.h'
    integer iunit,NtimeStep,NGLL1,NGLL2,nspec2D,ibloc,iii
    integer nbrc,Nbloc,Nbrec,nseis,irec,ns,komp,i,j,k,isp
    real(kind=CUSTOM_REAL) Trac(3,NtimeStep,NGLL1,NGLL2,nspec2D)
    double precision Trcsg(NtimeStep,3,nseis,Nbrec), corrfac
    integer Ig(3,NGLL1,NGLL2,nspec2D)

    !corrfac=dble(1e-4)  ! pour gemini
    corrfac=1.d0

    corrfac = corrfac * dble(1e6)  ! multiplication pour ne pas avoir des
                                   ! residus trop
                                   ! petit numeriquement
    !corrfac=dble(1e6) ! pour DSM
    !if (ibloc==1) then
    !   read(iunit) nbrc,Nbloc,Nbrec,nseis
    !endif
    !write(*,*) 'Read Vel'
    do irec = 1,Nbrec
       do ns = 1, nseis
          do komp = 1, 3
             read(iunit) (Trcsg(i,komp,ns,irec),i=1,nbrc)
             !write(*,*) (Trcsg(i,komp,ns,irec),i=1,nbrc)
          enddo
       enddo
    enddo
!!$    do i=1,nbrc
!!$       write(*,*) Trcsg(i,1,37,83),Trcsg(i,2,37,83),Trcsg(i,3,37,83)
!!$    enddo
    !read(*,*) iii
    do isp = 1, nspec2D
       do k=1,NGLL2
          do j=1,NGLL1
             do i=1,nbrc
                do komp=1,3
                !write(*,*) Ig(1,j,k,isp),Ig(2,j,k,isp)
                   Trac(komp,i,j,k,isp) =  Trcsg(i,komp,Ig(1,j,k,isp),Ig(2,j,k,isp))  * corrfac
                enddo
             enddo
          enddo
       enddo
    enddo
    return
!!$    do i=1,nbrc
!!$       write(*,*) Trac(1,i,1,1,10),Trac(2,i,1,1,10),Trac(3,i,1,1,10)
!!$    enddo
!!$    read(*,*) iii
  end subroutine ReadVel



 subroutine ReadTract(iunit,Trcsg,Trac,NtimeStep,NGLL1,NGLL2,nspec2D,nbrc,Nbloc,Nbrec,nseis,Ig)
    implicit none
    include 'constants.h'
    integer iunit,NtimeStep,NGLL1,NGLL2,nspec2D,ibloc,iii
    integer nbrc,Nbloc,Nbrec,nseis,irec,ns,komp,i,j,k,isp
    real(kind=CUSTOM_REAL) Trac(3,NtimeStep,NGLL1,NGLL2,nspec2D)
    double precision Trcsg(NtimeStep,3,nseis,Nbrec), corrfac
    integer Ig(3,NGLL1,NGLL2,nspec2D)

    !corrfac=dble(1e-4)  ! pour gemini
    corrfac=dble(1e6) ! pour DSM
    corrfac = corrfac * dble(1e6)  ! multiplication pour ne pas avoir des residus trop
                        ! petit numeriquement
    !if (ibloc==1) then
    !   read(iunit) nbrc,Nbloc,Nbrec,nseis
    !endif
    !write(*,*) 'Read Tract',iunit,Nbrec
    do irec = 1,Nbrec
       do ns = 1, nseis
          do komp = 1, 3
             read(iunit) (Trcsg(i,komp,ns,irec),i=1,nbrc)
             !write(*,*) (Trcsg(i,komp,ns,irec),i=1,nbrc)
          enddo
       enddo
    enddo
!!$    do i=1,nbrc
!!$       write(*,*) Trcsg(i,1,37,83),Trcsg(i,2,37,83),Trcsg(i,3,37,83)
!!$    enddo
    !read(*,*) iii
    do isp = 1, nspec2D
       do k=1,NGLL2
          do j=1,NGLL1
             do i=1,nbrc
                do komp=1,3
                !write(*,'(7i12)') komp,i,j,k,isp,Ig(1,j,k,isp),Ig(2,j,k,isp)
                   Trac(komp,i,j,k,isp) =  Trcsg(i,komp,Ig(1,j,k,isp),Ig(2,j,k,isp))  * corrfac
                enddo
             enddo
          enddo
       enddo
    enddo
    return
!!$    do i=1,nbrc
!!$       write(*,*) Trac(1,i,1,1,10),Trac(2,i,1,1,10),Trac(3,i,1,1,10)
!!$    enddo
!!$    read(*,*) iii
  end subroutine ReadTract


