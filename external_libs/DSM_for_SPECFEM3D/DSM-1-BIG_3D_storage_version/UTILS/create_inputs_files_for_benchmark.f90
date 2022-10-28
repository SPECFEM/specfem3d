program create_input_files_for_benchmark
  implicit none
  character(len=250) input_file,path
  double precision dp(6)
  integer id(3)

  ! input file name
  read(*,'(a)') input_file
  call  open_parameter_file(input_file)

  ! param.in -------------------------------------------------------------------
  open(10,file='params.in')
  call read_value_string(path, 'DSM_BINARY_PATH')
  write(10,*) 'BIN='//trim(path)
  call read_value_string(path, 'SPECFEM3D_BINARY_PATH')
  write(10,*) 'BINSEM='//trim(path)
  call read_value_string(path, 'SHELL_SCRIPT_PATH')
  write(10,*) 'SCRIPTS='//trim(path)
  write(10,*) 'DSM_tractions=$(pwd)'//'/DSM_tractions/'
  write(10,*) 'OUT=$DSM_tractions'
  write(10,*) 'REP=Tract/'
  call read_value_string(path, 'DSM_INPUT_DIR')
  write(10,*) 'IN_DSM=$(pwd)/'//trim(path)
  call read_value_string(path, 'MESH_DIR')
  write(10,*) 'MESH=$(pwd)/'//trim(path)
  close(10)

  ! ParFileInterface ---------------------------------------------------------
  open(10,file='ParFileInterface')
  write(10,*) '../Tract/'
  write(10,*) '../MESH/'
  write(10,*) '../OUTPUT_FILES/DATABASES_MPI/'
  write(10,*) '../DATA/DSM_tractions_for_specfem3D/'
  close(10)

 ! input file for DSM write coef --------------------------------------------------------
 open(10,file='input_dsm_for_write_coef')
 write(10,*) './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,*) trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,*) trim(path)
 write(10,*) 'st'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_integer(id(1),'IMAX')
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6f4.1)') dp(1:6)
 write(10,'(a3)') 'end'



 ! input file xmin ------------------------------------------
 open(10,file='input_dsm_for_read_xmin')
 write(10,*) './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,*) trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,*) trim(path)
 write(10,'(a6)') 'stxmin'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_integer(id(1),'IMAX')
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6f4.1)') dp(1:6)
 write(10,*) '../'
 call  read_value_integer(id(1),'IFRQMIN')
 call  read_value_integer(id(2),'IFRQMAX')
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 write(10,'(a3)') 'end'

 ! input file xmax ------------------------------------------
 open(10,file='input_dsm_for_read_xmax')
 write(10,*) './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,*) trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,*) trim(path)
 write(10,'(a6)') 'stxmax'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_integer(id(1),'IMAX')
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6f4.1)') dp(1:6)
 write(10,*) '../'
 call  read_value_integer(id(1),'IFRQMIN')
 call  read_value_integer(id(2),'IFRQMAX')
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 write(10,'(a3)') 'end'

! input file ymin ------------------------------------------
 open(10,file='input_dsm_for_read_ymin')
 write(10,*) './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,*) trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,*) trim(path)
 write(10,'(a6)') 'stymin'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_integer(id(1),'IMAX')
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6f4.1)') dp(1:6)
 write(10,*) '../'
 call  read_value_integer(id(1),'IFRQMIN')
 call  read_value_integer(id(2),'IFRQMAX')
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 write(10,'(a3)') 'end'

! input file ymax ------------------------------------------
 open(10,file='input_dsm_for_read_ymax')
 write(10,*) './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,*) trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,*) trim(path)
 write(10,'(a6)') 'stymax'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_integer(id(1),'IMAX')
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6f4.1)') dp(1:6)
 write(10,*) '../'
 call  read_value_integer(id(1),'IFRQMIN')
 call  read_value_integer(id(2),'IFRQMAX')
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
 write(10,'(a3)') 'end'

 ! input file zmin ------------------------------------------
 open(10,file='input_dsm_for_read_zmin')
 write(10,*) './'
 call read_value_string(path, 'FILE_MODEL_1D')
 write(10,*) trim(path)
 call read_value_string(path, 'FILE_OUT_DSM')
 write(10,*) trim(path)
 write(10,'(a6)') 'stzmin'
 call read_value_double_precision(dp(1),'TLEN')
 write(10,*) dp(1)/10.d0
 call read_value_double_precision(dp(1),'SRC_DEPTH')
 call read_value_double_precision(dp(2),'SRC_LAT')
 call read_value_double_precision(dp(3),'SRC_LON')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_integer(id(1),'IMAX')
 write(10,*) '0',id(1)
 write(10,*) '0'
 call read_value_double_precision(dp(1),'MRR')
 call read_value_double_precision(dp(2),'MTT')
 call read_value_double_precision(dp(3),'MPP')
 call read_value_double_precision(dp(4),'MRT')
 call read_value_double_precision(dp(5),'MRP')
 call read_value_double_precision(dp(6),'MTP')
 write(10,'(6f4.1)') dp(1:6)
 write(10,*) '../'
 call  read_value_integer(id(1),'IFRQMIN')
 call  read_value_integer(id(2),'IFRQMAX')
 write(10,*) id(1),id(2)
 call read_value_double_precision(dp(1),'SAMPLING')
 call read_value_double_precision(dp(2),'TSTART')
 call read_value_double_precision(dp(3),'TEND')
 write(10,*) dp(1),dp(2),dp(3)
 call read_value_double_precision(dp(1),'FLOW')
 call read_value_double_precision(dp(2),'FHIGH')
 write(10,*) dp(1),dp(2)
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
 close(10)
end program create_input_files_for_benchmark



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
   call param_open(filename, len(filename), ierr);
 end subroutine open_parameter_file
