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

!--------------------------------------------------------------------------------------------------
!
! IPATI
!
! based on given rho and vp structure for GLL files
!
!--------------------------------------------------------------------------------------------------

  subroutine model_ipati(myrank,nspec,LOCAL_PATH)

  use constants, only: IMAIN,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IIN
  use create_regions_mesh_ext_par

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname_lp,filename

  ! -----------------------------------------------------------------------------

  ! note: vp not vs structure is available (as is often the case in exploration seismology),
  ! scaling factor
  real, parameter :: SCALING_FACTOR = 1.0/1.8

  ! -----------------------------------------------------------------------------

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using IPATI model from: ',trim(LOCAL_PATH)
    write(IMAIN,*) '     scaling factor: ',SCALING_FACTOR
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'proc',myrank,'_'

  ! density
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 609')
  if (ier /= 0) stop 'error allocating array rho_read'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(IIN) rho_read
  close(IIN)

  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 610')
  if (ier /= 0) stop 'error allocating array vp_read'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(IIN) vp_read
  close(IIN)

  ! vs scaled from vp
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 611')
  if (ier /= 0) stop 'error allocating array vs_read'

  ! scaling
  vs_read = vp_read * SCALING_FACTOR

  ! isotropic model parameters
  rhostore    = rho_read
  kappastore  = rhostore * ( vp_read * vp_read - FOUR_THIRDS * vs_read * vs_read )
  mustore     = rhostore * vs_read * vs_read
  rho_vp = rhostore * vp_read
  rho_vs = rhostore * vs_read

  ! free memory
  deallocate( rho_read,vp_read,vs_read)

  end subroutine model_ipati

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_ipati_water(myrank,nspec,LOCAL_PATH)

  use constants, only: IMAIN,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IIN
  use create_regions_mesh_ext_par

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ispec,ier
  character(len=MAX_STRING_LEN) :: prname_lp,filename

  ! -----------------------------------------------------------------------------

  ! note: vp not vs structure is available (as is often the case in exploration seismology),
  ! scaling factor
  real, parameter :: SCALING_FACTOR = 1.0/1.8

  ! -----------------------------------------------------------------------------

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using IPATI_WATER model from: ',trim(LOCAL_PATH)
    write(IMAIN,*) '     scaling factor: ',SCALING_FACTOR
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'proc',myrank,'_'

  ! density
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 612')
  if (ier /= 0) stop 'error allocating array rho_read'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(IIN) rho_read
  close(IIN)

  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 613')
  if (ier /= 0) stop 'error allocating array vp_read'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(IIN) vp_read
  close(IIN)

  ! vs scaled from vp
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 614')
  if (ier /= 0) stop 'error allocating array vs_read'

  ! scaling
  vs_read = vp_read * SCALING_FACTOR

  ! overwrites only elastic elements
  do ispec=1,nspec
    ! assumes water layer with acoustic elements are set properly
    ! only overwrites elastic elements
    if (ispec_is_elastic(ispec)) then
      ! isotropic model parameters
      rhostore(:,:,:,ispec) = rho_read(:,:,:,ispec)
      kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * ( vp_read(:,:,:,ispec) * vp_read(:,:,:,ispec) &
                                    - FOUR_THIRDS * vs_read(:,:,:,ispec) * vs_read(:,:,:,ispec) )
      mustore(:,:,:,ispec) = rhostore(:,:,:,ispec) * vs_read(:,:,:,ispec) * vs_read(:,:,:,ispec)
      rho_vp(:,:,:,ispec) = rhostore(:,:,:,ispec) * vp_read(:,:,:,ispec)
      rho_vs(:,:,:,ispec) = rhostore(:,:,:,ispec) * vs_read(:,:,:,ispec)
    endif
  enddo

  ! free memory
  deallocate( rho_read,vp_read,vs_read)

  end subroutine model_ipati_water



