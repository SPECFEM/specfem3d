!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------
!
!
! example model to couple with an injection technique (FK,DSM,AXISEM)
!

  module model_coupled_par

  ! VM VM my model for DSM coupling
  double precision, dimension (:,:), allocatable :: vpv_1D,vsv_1D,density_1D
  double precision, dimension (:), allocatable :: zlayer
  double precision, dimension (:), allocatable :: smooth_vp,smooth_vs
  integer :: ilayer,nlayer,ncoeff,ndeg_poly

  double precision :: ZREF,OLON,OLAT

  end module model_coupled_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_coupled_broadcast()

! standard routine to setup model

  use model_coupled_par

  use constants

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  ! safety check
  if (.not. (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH)) then
    print *,'Error: model coupling requires coupling with injection technique or mesh a chunk of the earth'
    stop 'Error model coupling'
  endif

  call read_model_for_coupling_or_chunk()

  end subroutine model_coupled_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_for_coupling_or_chunk()

  use constants, only: IMAIN,IN_DATA_FILES,myrank

  use model_coupled_par !! VM VM custom subroutine for coupling with DSM

  implicit none

  ! local parameters
  character(len=256) :: filename
  integer :: i,cc,ier
  double precision :: aa,bb

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading model for coupling or mesh a chunk of the earth...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'/coeff_poly_deg12'
  open(27,file=trim(filename),iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',filename
    stop 'Error opening coupled coeff_poly_deg12 file'
  endif

  read(27,*) ndeg_poly
  allocate(smooth_vp(0:ndeg_poly),smooth_vs(0:ndeg_poly),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 620')
  do i = ndeg_poly,0,-1
    read(27,*) aa,bb,cc
    smooth_vp(i) = aa
    smooth_vs(i) = bb
    ! write(*,*) a,b
  enddo
  close(27)

  !write(*,*) " Reading 1D model "

  filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'/model_1D.in'
  open(27,file=trim(filename),iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',filename
    stop 'Error opening coupled model_1D.in file'
  endif

  read(27,*) nlayer,ncoeff
  allocate(vpv_1D(nlayer,ncoeff),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 621')
  allocate(vsv_1D(nlayer,ncoeff),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 622')
  allocate(density_1D(nlayer,ncoeff),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 623')
  allocate(zlayer(nlayer),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 624')
  do i = 1,nlayer
    read(27,*) zlayer(i)
    read(27,*) vpv_1D(i,:)
    read(27,*) vsv_1D(i,:)
    read(27,*) density_1D(i,:)
  enddo
  read(27,*) ZREF
  read(27,*) OLON,OLAT
  close(27)

  end subroutine read_model_for_coupling_or_chunk

!----------------------------------------------------------------
!! !! ================= VM VM custom subroutine for DSM coupling
!----------------------------------------------------------------

  subroutine model_coupled_FindLayer(x,y,z)

  use model_coupled_par

  implicit none
  integer il
  double precision radius
  double precision :: x,y,z

  radius =  dsqrt(x**2 + y**2 + (z+ZREF)**2) / 1000.d0
  il = 1
  do while (radius > zlayer(il) .and. il < nlayer)
     il = il + 1
  enddo
  il = il - 1
  ilayer = il

  end subroutine model_coupled_FindLayer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_coupled_values(xmesh,ymesh,zmesh,rho,vp,vs)

! given a GLL point, returns super-imposed velocity model values

  use constants, only: CUSTOM_REAL

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL),intent(out) :: vp,vs,rho

  ! local parameters
  double precision :: x,y,z
  double precision :: radius

  ! safety check
  if (.not. (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH)) then
      print *,'Error: model coupling requires coupling with injection technique or mesh a chunk of the earth'
    stop 'Error model coupling'
  endif

  ! GLL point location converted to real
  x = xmesh
  y = ymesh
  z = zmesh

  call  model_1D_coupling(x,y,z,rho,vp,vs,radius)

  end subroutine model_coupled_values

!----------------------------------------------------------------

  subroutine model_1D_coupling(x_eval,y_eval,z_eval,rho_final,vp_final,vs_final,r1)

  use constants, only: CUSTOM_REAL

  use model_coupled_par

  implicit none

  double precision,intent(in) :: x_eval,y_eval,z_eval
  real(kind=CUSTOM_REAL),intent(out) :: rho_final,vp_final,vs_final

  ! local parameters
  double precision :: r1,radius
  double precision :: rho,vp,vs

  !double precision, parameter :: Xtol = 1d-2

  radius = dsqrt(x_eval**2 + y_eval**2 + (z_eval+ZREF)**2)
  radius = radius / 1000.d0
  r1 = radius

  ! get vp,vs and rho
  radius = radius / zlayer(nlayer)

  vp = Interpol(vpv_1D,ilayer,radius,nlayer)
  vs = Interpol(vsv_1D,ilayer,radius,nlayer)
  rho = Interpol(density_1D,ilayer,radius,nlayer)

  ! converts units to m/s
  vp_final = vp * 1000.d0
  vs_final = vs * 1000.d0
  rho_final = rho * 1000.d0

  contains

    function Interpol(v,i,x,nl)

    implicit none
    integer i,nl
    double precision Interpol,x,v(nl,4)

    Interpol = v(i,1)+x*(v(i,2)+x*(v(i,3)+x*v(i,4)))

    end function Interpol

  end subroutine model_1D_coupling


