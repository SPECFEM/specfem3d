program create_input_files
  implicit none
  character(len=250) input_file,path
  double precision dp(6)
  integer id(3)
 !
  integer nb_proc_by_colors,reste_by_colors, nb_frq_by_colors,reste_by_frq
  integer  reste_ifrq,ncpus, sub_comm, ifrqmax, ifrq_by_proc, icpu,i,j,k,imin,imax,ifq

  integer ifmax_for_dsm,iamax_for_dsm
  double precision tlen_dsm,fhigh_dsm,acc_level
!
  double precision X,Y,Z,long,lat,rotation_matrix(3,3),ZREF
  double precision ANGULAR_WIDTH_XI,ANGULAR_WIDTH_ETA,DEPTH_CHUNK
  double precision lon_center_chunk,lat_center_chunk,chunk_azi

  ! input file name
  read(*,'(a)') input_file
  call  open_parameter_file(input_file)


  ! parameters --------
  call read_value_double_precision(tlen_dsm,'TLEN')
  call read_value_double_precision(fhigh_dsm,'FHIGH')
  call read_value_double_precision(acc_level,'ACCURACY_LEVEL')
  tlen_dsm = tlen_dsm / 10.
  ifmax_for_dsm = 1.4d0 * fhigh_dsm * tlen_dsm
  iamax_for_dsm = acc_level * ifmax_for_dsm

  ! param.in -------------------------------------------------------------------
  open(10,file='params.in')
  call read_value_string(path, 'DSM_BINARY_PATH')
  write(10,'(a)') 'BIN='//trim(path)
  call read_value_string(path, 'SPECFEM3D_BINARY_PATH')
  write(10,'(a)') 'BINSEM='//trim(path)
  call read_value_string(path, 'SHELL_SCRIPT_PATH')
  write(10,'(a)') 'SCRIPTS='//trim(path)
  write(10,'(a)') 'DSM_tractions=$(pwd)'//'/DSM_tractions/'
  write(10,'(a)') 'OUT=$DSM_tractions'
  write(10,'(a)') 'REP=Tract/'
  call read_value_string(path, 'DSM_INPUT_DIR')
  write(10,'(a)') 'IN_DSM=$(pwd)/'//trim(path)
  call read_value_string(path, 'MESH_DIR')
  write(10,'(a)') 'MESH=$(pwd)/'//trim(path)
  call read_value_string(path,'FILE_OUT_DSM')
  write(10,'(a)') 'MODELE_1D='//trim(path)
  close(10)

  ! ParFileInterface ---------------------------------------------------------
  open(10,file='ParFileInterface')
  write(10,*) './Tract/'
  write(10,*) './MESH/'
  write(10,*) './OUTPUT_FILES/DATABASES_MPI/'
  write(10,*) './DATA/DSM_tractions_for_specfem3D/'
  write(10,*) '1'
  close(10)

 ! input file for DSM write coef --------------------------------------------------------
 open(10,file='input_dsm_for_write_coef')
 write(10,'(a)') './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,'(a)') trim(path)
 write(10,'(a)') 'st'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 !call read_value_integer(id(1),'IMAX')
 id(1)=iamax_for_dsm
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6(f10.4,1x))') dp(1:6)
 write(10,'(a3)') 'end'

 ! input file xmin ------------------------------------------
 open(10,file='input_dsm_for_read_xmin')
 write(10,'(a)') './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,'(a)') trim(path)
 write(10,'(a6)') 'stxmin'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 !call read_value_integer(id(1),'IMAX')
 id(1)=iamax_for_dsm
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6(f10.4,1x))') dp(1:6)
 write(10,'(a)') '../'
 !call  read_value_integer(id(1),'IFRQMIN')
 !call  read_value_integer(id(2),'IFRQMAX')
 id(1)=0
 id(2)=ifmax_for_dsm
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 !call read_value_double_precision(dp(1),'FGAUSS')
 !write(10,*) dp(1)
 write(10,'(a3)') 'end'

 ! input file xmax ------------------------------------------
 open(10,file='input_dsm_for_read_xmax')
 write(10,'(a)') './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,'(a)') trim(path)
 write(10,'(a6)') 'stxmax'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 !call read_value_integer(id(1),'IMAX')
 id(1)=iamax_for_dsm
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6(f10.4,1x))') dp(1:6)
 write(10,'(a)') '../'
 !call  read_value_integer(id(1),'IFRQMIN')
 !call  read_value_integer(id(2),'IFRQMAX')
 id(1)=0
 id(2)=ifmax_for_dsm
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 !call read_value_double_precision(dp(1),'FGAUSS')
 !write(10,*) dp(1)
 write(10,'(a3)') 'end'

! input file ymin ------------------------------------------
 open(10,file='input_dsm_for_read_ymin')
 write(10,'(a)') './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,'(a)') trim(path)
 write(10,'(a6)') 'stymin'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 !call read_value_integer(id(1),'IMAX')
 id(1)=iamax_for_dsm
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6(f10.4,1x))') dp(1:6)
 write(10,'(a)') '../'
 !call  read_value_integer(id(1),'IFRQMIN')
 !call  read_value_integer(id(2),'IFRQMAX')
 id(1)=0
 id(2)=ifmax_for_dsm
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 !call read_value_double_precision(dp(1),'FGAUSS')
 !write(10,*) dp(1)
 write(10,'(a3)') 'end'

! input file ymax ------------------------------------------
 open(10,file='input_dsm_for_read_ymax')
 write(10,'(a)') './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,'(a)') trim(path)
 write(10,'(a6)') 'stymax'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 !call read_value_integer(id(1),'IMAX')
 id(1)=iamax_for_dsm
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6(f10.4,1x))') dp(1:6)
 write(10,'(a)') '../'
 !call  read_value_integer(id(1),'IFRQMIN')
 !call  read_value_integer(id(2),'IFRQMAX')
 id(1)=0
 id(2)=ifmax_for_dsm
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 !call read_value_double_precision(dp(1),'FGAUSS')
 !write(10,*) dp(1)
 write(10,'(a3)') 'end'

 ! input file zmin ------------------------------------------
 open(10,file='input_dsm_for_read_zmin')
 write(10,'(a)') './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,'(a)') trim(path)
 write(10,'(a6)') 'stzmin'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 !call read_value_integer(id(1),'IMAX')
 id(1)=iamax_for_dsm
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6(f10.4,1x))') dp(1:6)
 write(10,'(a)') '../'
 !call  read_value_integer(id(1),'IFRQMIN')
 !call  read_value_integer(id(2),'IFRQMAX')
 id(1)=0
 id(2)=ifmax_for_dsm
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 !call read_value_double_precision(dp(1),'FGAUSS')
 !write(10,*) dp(1)
 write(10,'(a3)') 'end'

 ! INPUT FILE FOR MESHER
 open(10,file='ParFileMeshChunk')
 write(10,'(a1)') '#'
 call read_value_double_precision(dp(1),'ANGULAR_WIDTH_XI_RAD')
 call read_value_double_precision(dp(2),'ANGULAR_WIDTH_ETA_RAD')
 write(10,*) dp(1),dp(2)
 write(10,'(a1)') '#'
 call read_value_double_precision(dp(1),'LON_CENTER')
 call read_value_double_precision(dp(2),'LAT_CENTER')
 call read_value_double_precision(dp(3),'AZI_CHUNK')
 write(10,*) dp(1),dp(2),dp(3)
 write(10,'(a1)') '#'
 call read_value_double_precision(dp(1),'DEPTH_CHUNK')
 write(10,*) dp(1)
 write(10,'(a1)') '#'
 call  read_value_integer(id(1),'NEL_LON')
 call  read_value_integer(id(2),'NEL_LAT')
 call  read_value_integer(id(3),'NEL_DEPTH')
 write(10,*) id(1),id(2),id(3)
 write(10,'(a1)') '#'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,'(a)') trim(path)
 close(10)

 !----- parallelisation file -------
 !call read_value_integer(ifrqmax,'IFRQMAX')
 ifrqmax=ifmax_for_dsm
 call read_value_integer(ncpus,'MPI_CPUS')
 call read_value_integer(sub_comm,'SUB_COMM')

 ifrq_by_proc = int((ifrqmax+1)/ncpus) ! quotient
 reste_ifrq   = (ifrqmax + 1) - ifrq_by_proc *  ncpus

 ! step 1
 open(10,file='./input_dsm/FrqsMpi.txt')
 write(10,*) '2'


 ifq=-1;imin=1;imax=ifrq_by_proc+1
 do icpu = 1,reste_ifrq
    write(10,*) icpu-1,ifq+imin,ifq+imax
    ifq=ifq+imax
 enddo

 imax=ifrq_by_proc
 do icpu = reste_ifrq+1,ncpus
    write(10,*) icpu-1,ifq+imin,ifq+imax
    ifq=ifq+imax
 enddo
 close(10)

 ! step 2
 open(10,file='./input_dsm/Double_para.txt')

 nb_proc_by_colors = ncpus / sub_comm
 reste_by_colors =  ncpus - sub_comm * nb_proc_by_colors

 nb_frq_by_colors =  (ifrqmax+1) / sub_comm
 reste_by_frq = (ifrqmax+1) - sub_comm * nb_frq_by_colors

 write(10,*) sub_comm,  ncpus
 imin=0;k=0
 do i=1, sub_comm

    if (i <= reste_by_frq) then
       imax = imin+nb_frq_by_colors-1 + 1
    else
       imax = imin+nb_frq_by_colors-1
    endif

    write(10,*) i-1,imin,imax,nb_proc_by_colors
    imin=imax+1
    if (j <= reste_by_colors) then
       do  j=1,nb_proc_by_colors+1
          write(10,*) k
          k=k+1
       enddo
    else
       do  j=1,nb_proc_by_colors
          write(10,*) k
          k=k+1
       enddo
    endif
 enddo
 close(10)

 call read_value_double_precision(ANGULAR_WIDTH_XI,'ANGULAR_WIDTH_XI_RAD')
 call read_value_double_precision(ANGULAR_WIDTH_ETA,'ANGULAR_WIDTH_ETA_RAD')
 call read_value_double_precision(DEPTH_CHUNK,'DEPTH_CHUNK')
 call read_value_double_precision(lon_center_chunk,'LON_CENTER')
 call read_value_double_precision(lat_center_chunk,'LAT_CENTER')
 call read_value_double_precision(chunk_azi,'AZI_CHUNK')
 call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)
!!$ call compute_ZREF(ZREF,ANGULAR_WIDTH_XI,ANGULAR_WIDTH_ETA,DEPTH_CHUNK)
 ZREF=0.d0
!-------- station file (lat long)----
 open(10,file='station_for_simulation.txt')
 open(20,file='DATA/STATIONS_test')

 do
    read(10,*,end=99) lat,long
    call  Geogr2Cart(X,Y,Z,long,lat,rotation_matrix,ZREF)
    write(20,*) x,y,z
 enddo
99 close(10)

end program create_input_files



! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(value_to_read, name)

  implicit none

  integer value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(value_to_read, name)

  implicit none

  double precision value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(value_to_read, name)

  implicit none

  logical value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(value_to_read, name)

  implicit none

  character(len=*) value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  value_to_read = string_read

  end subroutine read_value_string

 subroutine open_parameter_file(filename)
   character(len=250) filename
   integer ier

   ier = 0
   call param_open(filename, len(filename), ier);
 end subroutine open_parameter_file


!----------------------------------------------------------------------------------------------------------
!! routine for convertion coordiantes



  subroutine  compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)
  implicit none
  double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk,chunk_azi
  double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)

    ! je met le chunk en 0,0
    axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
    call rotation_matrix_axe(R00,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)

   ! rotation de l'azimuth du chunk
    axe_rotation(1)=-1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
    call rotation_matrix_axe(R0,axe_rotation,90.-chunk_azi)

    ! on met le chunk a la bonne latitude
    axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
    call rotation_matrix_axe(R1,axe_rotation,lat_center_chunk)

    ! on met le chunk a la bonne longitude
    axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=-1.d0
    call rotation_matrix_axe(R2,axe_rotation, lon_center_chunk)



    ! rotation resultante
    call compose4matrix(rotation_matrix,R2,R1,R0,R00)


  end subroutine compute_rotation_matrix

!
!   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
!
!       Vadim Monteiller Mars 2013
!
!---------------------------------------------------------------------------------------------------------------------------------------------------------
! matrice de rotation 3D d'axe "axe" et d'angle theta (degres)
!
subroutine rotation_matrix_axe(R,axe,theta)
  implicit none
  double precision axe(3),theta,pi,deg2rad
  double precision R(3,3)
  double precision c,s,ux,uy,uz,norme_axe
  integer i,j

  pi=3.1415926535897932d0
  deg2rad = pi / 180.d0
  ! on normalise l'axe
  norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)

  ! composantes de l'axe
  ux=axe(1)/norme_axe
  uy=axe(2)/norme_axe
  uz=axe(3)/norme_axe

  ! on calcule le cos et sin
  c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)

  ! matrice de rotation complexe
  R(1,1)=(ux**2 + (1.d0-ux**2)*c)
  R(1,2)=(ux*uy*(1.d0-c)-uz*s)
  R(1,3)=(ux*uy*(1.d0-c)+uy*s)

  R(2,1)=(ux*uy*(1.d0-c)+uz*s)
  R(2,2)=(uy**2+(1.d0-uy**2)*c)
  R(2,3)=(uy*uz*(1.d0-c)-ux*s)

  R(3,1)=(ux*uz*(1.d0-c)-uy*s)
  R(3,2)=(uy*uz*(1.d0-c)+ux*s)
  R(3,3)=(uz**2+(1.d0-uz**2)*c)

!!$  write(49,*) ' MATRICE ROTATION '
!!$  write(49,*) R(1,:)
!!$  write(49,*) R(2,:)
!!$  write(49,*) R(3,:)
!!$  write(49,*)
end subroutine rotation_matrix_axe

!-------------------------------------------------------------------------------------------------------------------------------------------------------
! R=R2*R1*R0
subroutine compose4matrix(R,R00,R0,R1,R2)
  implicit none
  double precision R(3,3),R0(3,3),R1(3,3),R2(3,3),R00(3,3),Rtmp(3,3)
  integer i,j,k


 R(:,:)=0.d0
  ! multiplication R=R0*R00
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R0(i,k)*R00(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R1*R
 Rtmp=R
 R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R1(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R2*R
 Rtmp=R
 R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R2(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

!!$  write(49,*) ' MATRICE ROTATION COMPLETE '
!!$  write(49,*) R(1,:)
!!$  write(49,*) R(2,:)
!!$  write(49,*) R(3,:)
!!$  write(49,*)
end subroutine compose4matrix

!---------------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine Geogr2Cart(X,Y,Z,long,lat,rotation_matrix,ZREF)
  implicit none
  double precision vector_ori(3),vector_rotated(3),rotation_matrix(3,3)
  double precision X,Y,Z,long,lat
  double precision  radius_earth, deg2rad,ZREF
  integer i,j,NDIM


  NDIM=3
  deg2rad = 3.141592653589793238d0/ 180.d0
  radius_earth = 6371000.d0


  vector_ori(1) = radius_earth * cos(long*deg2rad) * cos(lat*deg2rad)
  vector_ori(2) = radius_earth * sin(long*deg2rad) * cos(lat*deg2rad)
  vector_ori(3) = radius_earth *                     sin(lat*deg2rad)

  do i = 1,NDIM
     vector_rotated(i) = 0.d0
     do j = 1,NDIM
        vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
     enddo
  enddo

  x = vector_rotated(1)
  y = vector_rotated(2)
  z = vector_rotated(3) - ZREF


end subroutine Geogr2Cart

!!$
!!$subroutine compute_ZREF(ZREF,ANGULAR_WIDTH_XI,ANGULAR_WIDTH_ETA,DEPTH_CHUNK)
!!$  implicit none
!!$  deg2rad=3.141592653589793d0/180.d0
!!$  ANGULAR_WIDTH_XI = deg2rad * ANGULAR_WIDTH_XI
!!$  ANGULAR_WIDTH_ETA = deg2rad * ANGULAR_WIDTH_ETA
!!$  DEPTH_CHUNK  = DEPTH_CHUNK * 1000.d0
!!$
!!$
!!$
!!$end subroutine compute_ZREF
