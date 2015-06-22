!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine save_arrays_solver(prname,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                     kappav,muv,ibool,rmass,nspec,nglob,myrank,NPROCTOT,xstore,ystore,zstore)

  implicit none

  include "constants.h"

  integer :: nspec,nglob,myrank,NPROCTOT

! arrays with jacobian matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,kappav,muv

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! mass matrix
  real(kind=CUSTOM_REAL) rmass(nglob)

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

!!!!!!!!!!!!!!!!!!  integer :: i,j,k,ispec

  real(kind=CUSTOM_REAL) :: memory_size

! processor identification
  character(len=150) prname

! estimate total memory size (the size of a real number is 4 bytes)
! we perform the calculation in single precision rather than integer
! to avoid integer overflow in the case of very large meshes
  memory_size = 4. * ((3.*NDIM + 1.) * NGLOB + 12. * real(NGLLX*NGLLY*NGLLZ)*real(NSPEC))
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'approximate total memory size that will be used by the solver in each slice = ',memory_size/1024./1024.,' Mb'
    write(IMAIN,*) 'i.e. = ',memory_size/1024./1024./1024.,' Gb'
#ifdef USE_MPI
    write(IMAIN,*)
    write(IMAIN,*) 'approximate total memory size that will be used by the solver for all the slices = ', &
      NPROCTOT*memory_size/1024./1024.,' Mb'
    write(IMAIN,*) 'i.e. = ',NPROCTOT*memory_size/1024./1024./1024.,' Gb'
#else
!! DK DK this line is completely useless and meaningless but I put it to
!! DK DK avoid a compiler warning about an unused variable
    memory_size = NPROCTOT
#endif
    write(IMAIN,*)
  endif

! open(unit=IOUT,file=prname(1:len_trim(prname))//'database.dat',status='unknown')
  open(unit=IOUT,file=prname(1:len_trim(prname))//'database.dat',status='unknown',form='unformatted')

! write real numbers here
  write(IOUT) xix
  write(IOUT) xiy
  write(IOUT) xiz
  write(IOUT) etax
  write(IOUT) etay
  write(IOUT) etaz
  write(IOUT) gammax
  write(IOUT) gammay
  write(IOUT) gammaz
  write(IOUT) kappav
  write(IOUT) muv

! write an integer here
  write(IOUT) ibool

! store the mass matrix
  write(IOUT) rmass

! store the coordinates of the mesh
  write(IOUT) xstore
  write(IOUT) ystore
  write(IOUT) zstore

! do ispec = 1,NSPEC
!   do k=1,NGLLZ
!     do j=1,NGLLY
!       do i=1,NGLLX
! write real numbers here
! write(IOUT,*) xix(i,j,k,ispec)
! write(IOUT,*) xiy(i,j,k,ispec)
! write(IOUT,*) xiz(i,j,k,ispec)
! write(IOUT,*) etax(i,j,k,ispec)
! write(IOUT,*) etay(i,j,k,ispec)
! write(IOUT,*) etaz(i,j,k,ispec)
! write(IOUT,*) gammax(i,j,k,ispec)
! write(IOUT,*) gammay(i,j,k,ispec)
! write(IOUT,*) gammaz(i,j,k,ispec)
! write(IOUT,*) kappav(i,j,k,ispec)
! write(IOUT,*) muv(i,j,k,ispec)

! write an integer here
! write(IOUT,*) ibool(i,j,k,ispec)
!       enddo
!     enddo
!   enddo
! enddo

! store the mass matrix
! do i = 1,nglob
!   write(IOUT,*) rmass(i)
! enddo

  close(IOUT)

  end subroutine save_arrays_solver

