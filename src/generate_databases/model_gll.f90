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
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!
!--------------------------------------------------------------------------------------------------

  subroutine model_gll(myrank,nspec,LOCAL_PATH)

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,MAX_STRING_LEN,ATTENUATION

  use create_regions_mesh_ext_par, only: rhostore,kappastore,mustore,rho_vp,rho_vs,qkappa_attenuation_store,qmu_attenuation_store

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname_lp,filename

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using GLL model from: ',trim(LOCAL_PATH)
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)// '/' //'proc',myrank,'_'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! if only vp structure is available (as is often the case in exploration seismology),
  !!! use lines for vp only

  ! density
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 647')
  if (ier /= 0) stop 'error allocating array rho_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: rho.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(28) rho_read
  close(28)

  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
  if (ier /= 0) stop 'error allocating array vp_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vp.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(28) vp_read
  close(28)

  ! vs
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 649')
  if (ier /= 0) stop 'error allocating array vs_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vs.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vs.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vs.bin file'
  endif

  read(28) vs_read
  close(28)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where density structure is not given
  !!! modify according to your desire

  !  rho_read = 1000.0
  !  where ( mustore > 100.0 )  &
  !           rho_read = (1.6612 * (vp_read / 1000.0)     &
  !                      -0.4720 * (vp_read / 1000.0)**2  &
  !                      +0.0671 * (vp_read / 1000.0)**3  &
  !                      -0.0043 * (vp_read / 1000.0)**4  &
  !                      +0.000106*(vp_read / 1000.0)**5)*1000.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where shear wavespeed structure is not given
  !!! modify according to your desire

  !   vs_read = 0.0
  !   where ( mustore > 100.0 )       vs_read = vp_read / sqrt(3.0)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! update arrays that will be saved and used in the solver xspecfem3D
  !!! the following part is neccessary if you uncommented something above

  rhostore(:,:,:,:) = rho_read(:,:,:,:)
  kappastore(:,:,:,:) = rhostore(:,:,:,:) * ( vp_read(:,:,:,:) * vp_read(:,:,:,:) &
                                              - FOUR_THIRDS * vs_read(:,:,:,:) * vs_read(:,:,:,:) )
  mustore(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:) * vs_read(:,:,:,:)

  rho_vp(:,:,:,:) = rhostore(:,:,:,:) * vp_read(:,:,:,:)
  rho_vs(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:)

  ! gets attenuation arrays from files
  if (ATTENUATION) then
    ! shear attenuation
    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: qmu.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'qmu.bin'
    open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading qmu.bin file'
    endif

    read(28) qmu_attenuation_store
    close(28)

    ! bulk attenuation
    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: qkappa.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'qkappa.bin'
    open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading qkappa.bin file'
    endif

    read(28) qkappa_attenuation_store
    close(28)
  endif

  ! free memory
  deallocate(rho_read,vp_read,vs_read)

  end subroutine model_gll

