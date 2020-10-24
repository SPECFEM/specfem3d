!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

subroutine read_model_iso()

! reads in current isotropic model: vp & vs & rho

  use tomography_model_iso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading isotropic model...'

  ! allocate arrays for storing the databases
  allocate(model_vp(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1053')
  allocate(model_vs(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1054')
  allocate(model_rho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1055')
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vp = 0.0_CUSTOM_REAL
  model_vs = 0.0_CUSTOM_REAL
  model_rho = 0.0_CUSTOM_REAL

  ! reads in current vp & vs & rho model files
  ! vp model
  fname = 'vp'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vp(:,:,:,1:nspec)
  close(IIN)

  ! vs model
  fname = 'vs'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vs(:,:,:,1:nspec)
  close(IIN)

  ! rho model
  fname = 'rho'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_rho(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  call min_all_cr(minval(model_vp),min_vp)
  call max_all_cr(maxval(model_vp),max_vp)

  call min_all_cr(minval(model_vs),min_vs)
  call max_all_cr(maxval(model_vs),max_vs)

  call min_all_cr(minval(model_rho),min_rho)
  call max_all_cr(maxval(model_rho),max_rho)

  if (myrank == 0) then
    print *
    print *,'initial models:'
    print *,'  vp min/max : ',min_vp,max_vp
    print *,'  vs min/max : ',min_vs,max_vs
    print *,'  rho min/max: ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax',status='unknown')
    write(IOUT,*) '#min_vs #max_vs #min_vp #max_vp #min_rho #max_rho'
    write(IOUT,'(6e24.12)') min_vs,max_vs,min_vp,max_vp,min_rho,max_rho
    close(IOUT)
  endif

  ! global addressing
  call read_model_database()

end subroutine read_model_iso

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model_tiso()

! reads in current transverse isotropic model: vpv.. & vsv.. & eta & rho

  use tomography_model_tiso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,min_vsv,min_vsh, &
    max_vpv,max_vph,max_vsv,max_vsh,min_eta,max_eta,min_rho,max_rho
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading model...'

  ! allocate arrays for storing the databases
  allocate(model_vpv(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1056')
  allocate(model_vph(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1057')
  allocate(model_vsv(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1058')
  allocate(model_vsh(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1059')
  allocate(model_eta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1060')
  allocate(model_rho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1061')
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vpv = 0.0_CUSTOM_REAL
  model_vph = 0.0_CUSTOM_REAL
  model_vsv = 0.0_CUSTOM_REAL
  model_vsh = 0.0_CUSTOM_REAL
  model_eta = 0.0_CUSTOM_REAL
  model_rho = 0.0_CUSTOM_REAL

  ! reads in model files
  ! vpv model
  fname = 'vpv'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vpv(:,:,:,1:nspec)
  close(IIN)

  ! vph model
  fname = 'vph'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vph(:,:,:,1:nspec)
  close(IIN)

  ! vsv model
  fname = 'vsv'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vsv(:,:,:,1:nspec)
  close(IIN)

  ! vsh model
  fname = 'vsh'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vsh(:,:,:,1:nspec)
  close(IIN)

  ! eta model
  fname = 'eta'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_eta(:,:,:,1:nspec)
  close(IIN)

  ! rho model
  fname = 'rho'
  write(m_file,'(a,i6.6,a)') trim(INPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_rho(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  call min_all_cr(minval(model_vpv),min_vpv)
  call max_all_cr(maxval(model_vpv),max_vpv)

  call min_all_cr(minval(model_vph),min_vph)
  call max_all_cr(maxval(model_vph),max_vph)

  call min_all_cr(minval(model_vsv),min_vsv)
  call max_all_cr(maxval(model_vsv),max_vsv)

  call min_all_cr(minval(model_vsh),min_vsh)
  call max_all_cr(maxval(model_vsh),max_vsh)

  call min_all_cr(minval(model_eta),min_eta)
  call max_all_cr(maxval(model_eta),max_eta)

  call min_all_cr(minval(model_rho),min_rho)
  call max_all_cr(maxval(model_rho),max_rho)

  if (myrank == 0) then
    print *
    print *,'initial models:'
    print *,'  vpv min/max: ',min_vpv,max_vpv
    print *,'  vph min/max: ',min_vph,max_vph
    print *,'  vsv min/max: ',min_vsv,max_vsv
    print *,'  vsh min/max: ',min_vsh,max_vsh
    print *,'  eta min/max: ',min_eta,max_eta
    print *,'  rho min/max: ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax',status='unknown')
    write(IOUT,*) '#min_vsv #max_vsv #min_vsh #max_vsh #min_vpv #max_vpv #min_vph #max_vph ' &
               // '#min_eta #max_eta #min_rho #max_rho'
    write(IOUT,'(12e24.12)') min_vsv,max_vsv,min_vsh,max_vsh,min_vpv,max_vpv,min_vph,max_vph, &
                            min_eta,max_eta,min_rho,max_rho
    close(IOUT)
  endif

  ! global addressing
  call read_model_database()

end subroutine read_model_tiso


!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model_database()

! reads in current transverse isotropic model: vpv.. & vsv.. & eta & rho

  use tomography_par

  implicit none
  integer :: ival,ier
  character(len=MAX_STRING_LEN) :: m_file

  ! global addressing
  write(m_file,'(a,i6.6,a)') trim(INPUT_DATABASES_DIR)//'proc',myrank,trim(REG)//'external_mesh.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif

  read(IIN) ival !nspec
  if (ival /= nspec) then
    print *,'Error: invalid nspec ',ival,' found, should be ',nspec
    call exit_mpi(myrank,'Error invalid nspec value in external_mesh.bin')
  endif

  read(IIN) ival !nglob
  if (ival /= nglob) then
    print *,'Error: invalid nglob ',ival,' found, should be ',nglob
    call exit_mpi(myrank,'Error invalid nglob value in external_mesh.bin')
  endif

  read(IIN) ival ! skip nspec_irregular

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1062')
  if (ier /= 0) stop 'Error allocating ibool array for databases'

  ! mesh node locations
  allocate(x(NGLOB),y(NGLOB),z(NGLOB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1063')
  if (ier /= 0) stop 'Error allocating x/y/z arrays for mesh nodes'

  read(IIN) ibool(:,:,:,1:nspec)

  read(IIN) x(1:nglob)
  read(IIN) y(1:nglob)
  read(IIN) z(1:nglob)
  close(IIN)

end subroutine read_model_database

