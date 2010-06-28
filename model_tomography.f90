!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

! generic tomography file
!
! note: the idea is to use an external, tomography velocity model 
!
! most of the routines here are place-holders, please add/implement your own routines
!

  module tomography

  include "constants.h"

  ! for external tomography....
  ! (regular spaced, xyz-block file in ascii)
  character (len=80) :: TOMO_FILENAME = 'DATA/veryfast_tomography_abruzzo_complete.xyz' 
  
  ! model dimensions
  double precision :: ORIG_X,ORIG_Y,ORIG_Z
  double precision :: END_X,END_Y,END_Z   
  double precision :: SPACING_X,SPACING_Y,SPACING_Z  

  ! model parameter records    
  real(kind=CUSTOM_REAL), dimension (:), allocatable :: vp_tomography,vs_tomography,rho_tomography,z_tomography 

  ! model entries
  integer :: NX,NY,NZ    
  integer :: nrecord

  ! min/max statistics
  double precision :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX      

  end module tomography

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_tomography_broadcast(myrank)

  implicit none

  ! include "constants.h"
  ! include "precision.h"
  ! include 'mpif.h'  
  integer :: myrank

  ! all processes read in same file
  ! note: for a high number of processes this might lead to a bottleneck
  call read_model_tomography(myrank)

  ! otherwise:
  
  ! only master reads in model file      
  !if(myrank == 0) call read_external_model()      
  ! broadcast the information read on the master to the nodes, e.g.
  !call MPI_BCAST(nrecord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  !if( myrank /= 0 ) allocate( vp_tomography(1:nrecord) )
  !call MPI_BCAST(vp_tomography,size(vp_tomography),CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)  

  end subroutine model_tomography_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_tomography(myrank)

! start magnoni 29/11/09
! read Vp Vs and rho from extracted text file

! assuming that only tomography undefined material is allowed.... 
! and all the tomographic regions are collect inside one file called TOMO_FILENAME with homogenous resolution
! this could be problematic for example if the tomographic regions have different resolution 
! leading to a waste of memory and cpu time in the partitioning process 

  use tomography
  
  implicit none

  integer :: myrank
  
  ! local parameters
  real(kind=CUSTOM_REAL) :: x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo      
  integer :: irecord,ier

  !TOMO_FILENAME='DATA/veryfast_tomography_abruzzo_complete.xyz'
  ! probably the simple position for the filename is the constat.h
  ! but it is also possible to include the name of the file in the material file (therefore in the undef_mat_prop)
  ! if we want more than one tomofile (Examples: 2 file with a differente resolution 
  ! as in los angeles case we need to loop over mat_ext_mesh(1,ispec)... 
  ! it is a possible solution )      
  !  magnoni 1/12/09
  open(unit=27,file=TOMO_FILENAME,status='old',iostat=ier) 
  if( ier /= 0 ) call exit_MPI(myrank,'error reading tomography file')
  
  ! reads in model dimensions
  read(27,*) ORIG_X, ORIG_Y, ORIG_Z, END_X, END_Y, END_Z  
  read(27,*) SPACING_X, SPACING_Y, SPACING_Z 
  read(27,*) NX, NY, NZ 
  read(27,*) VP_MIN, VP_MAX, VS_MIN, VS_MAX, RHO_MIN, RHO_MAX 

  nrecord = NX*NY*NZ   

  ! allocates model records
  allocate(vp_tomography(1:nrecord), &
          vs_tomography(1:nrecord), &
          rho_tomography(1:nrecord), &
          z_tomography(1:nrecord),stat=ier) 
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays') 

  ! reads in record sections
  do irecord = 1,nrecord   
    read(27,*) x_tomo,y_tomo,z_tomo,vp_tomo,vs_tomo,rho_tomo      
    
    ! stores record values
    vp_tomography(irecord) = vp_tomo
    vs_tomography(irecord) = vs_tomo
    rho_tomography(irecord) = rho_tomo
    z_tomography(irecord) = z_tomo
  enddo 
  
  close(27)   
                                                                
  end subroutine read_model_tomography


!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_tomography(x_eval,y_eval,z_eval, &                      
                             rho_final,vp_final,vs_final)

  use tomography

  implicit none

  !integer, intent(in) :: NX,NY,NZ
  !real(kind=CUSTOM_REAL), dimension(1:NX*NY*NZ), intent(in) :: vp_tomography,vs_tomography,rho_tomography,z_tomography
  !double precision, intent(in) :: ORIG_X,ORIG_Y,ORIG_Z,SPACING_X,SPACING_Y,SPACING_Z
  !double precision, intent(in) :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX  

  double precision, intent(in) :: x_eval,y_eval,z_eval
  real(kind=CUSTOM_REAL), intent(out) :: vp_final,vs_final,rho_final

  ! local parameters
  integer :: ix,iy,iz
  integer :: p0,p1,p2,p3,p4,p5,p6,p7

  double precision :: spac_x,spac_y,spac_z
  double precision :: gamma_interp_x,gamma_interp_y
  double precision :: gamma_interp_z1,gamma_interp_z2,gamma_interp_z3, &
    gamma_interp_z4,gamma_interp_z5,gamma_interp_z6,gamma_interp_z7,gamma_interp_z8
  real(kind=CUSTOM_REAL) :: vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8, &
    vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8,rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8

  ! determine spacing and cell for linear interpolation
  spac_x = (x_eval - ORIG_X) / SPACING_X
  spac_y = (y_eval - ORIG_Y) / SPACING_Y
  spac_z = (z_eval - ORIG_Z) / SPACING_Z

  ix = int(spac_x)
  iy = int(spac_y)
  iz = int(spac_z)

  gamma_interp_x = spac_x - dble(ix)
  gamma_interp_y = spac_y - dble(iy)

  ! suppress edge effects for points outside of the model SPOSTARE DOPO
  if(ix < 0) then
    ix = 0
    gamma_interp_x = 0.d0
  endif
  if(ix > NX-2) then
    ix = NX-2
    gamma_interp_x = 1.d0
  endif

  if(iy < 0) then
    iy = 0
    gamma_interp_y = 0.d0
  endif
  if(iy > NY-2) then
    iy = NY-2
    gamma_interp_y = 1.d0
  endif

  if(iz < 0) then
     iz = 0
  !   gamma_interp_z = 0.d0
  endif
  if(iz > NZ-2) then
     iz = NZ-2
  !  gamma_interp_z = 1.d0
  endif


  ! define 8 corners of interpolation element
  p0 = ix+iy*NX+iz*(NX*NY)
  p1 = (ix+1)+iy*NX+iz*(NX*NY)
  p2 = (ix+1)+(iy+1)*NX+iz*(NX*NY)
  p3 = ix+(iy+1)*NX+iz*(NX*NY)
  p4 = ix+iy*NX+(iz+1)*(NX*NY)
  p5 = (ix+1)+iy*NX+(iz+1)*(NX*NY)
  p6 = (ix+1)+(iy+1)*NX+(iz+1)*(NX*NY)
  p7 = ix+(iy+1)*NX+(iz+1)*(NX*NY)

  if(z_tomography(p4+1) == z_tomography(p0+1)) then
          gamma_interp_z1 = 1.d0
      else
          gamma_interp_z1 = (z_eval-z_tomography(p0+1))/(z_tomography(p4+1)-z_tomography(p0+1))   
  endif
  if(gamma_interp_z1 > 1.d0) then
          gamma_interp_z1 = 1.d0
  endif
  if(gamma_interp_z1 < 0.d0) then
          gamma_interp_z1 = 0.d0
  endif
      
     
  if(z_tomography(p5+1) == z_tomography(p1+1)) then
          gamma_interp_z2 = 1.d0
      else
          gamma_interp_z2 = (z_eval-z_tomography(p1+1))/(z_tomography(p5+1)-z_tomography(p1+1))
  endif
  if(gamma_interp_z2 > 1.d0) then
          gamma_interp_z2 = 1.d0
  endif
  if(gamma_interp_z2 < 0.d0) then
          gamma_interp_z2 = 0.d0
  endif
      
     
  if(z_tomography(p6+1) == z_tomography(p2+1)) then
          gamma_interp_z3 = 1.d0
      else
          gamma_interp_z3 = (z_eval-z_tomography(p2+1))/(z_tomography(p6+1)-z_tomography(p2+1))
  endif
  if(gamma_interp_z3 > 1.d0) then
          gamma_interp_z3 = 1.d0
  endif
  if(gamma_interp_z3 < 0.d0) then
          gamma_interp_z3 = 0.d0
  endif
      
     
  if(z_tomography(p7+1) == z_tomography(p3+1)) then
          gamma_interp_z4 = 1.d0
      else
          gamma_interp_z4 = (z_eval-z_tomography(p3+1))/(z_tomography(p7+1)-z_tomography(p3+1))
  endif
  if(gamma_interp_z4 > 1.d0) then
          gamma_interp_z4 = 1.d0
  endif
  if(gamma_interp_z4 < 0.d0) then
          gamma_interp_z4 = 0.d0
  endif
      
  gamma_interp_z5 = 1. - gamma_interp_z1
  gamma_interp_z6 = 1. - gamma_interp_z2
  gamma_interp_z7 = 1. - gamma_interp_z3
  gamma_interp_z8 = 1. - gamma_interp_z4

  vp1 = vp_tomography(p0+1)
  vp2 = vp_tomography(p1+1)
  vp3 = vp_tomography(p2+1)
  vp4 = vp_tomography(p3+1)
  vp5 = vp_tomography(p4+1)
  vp6 = vp_tomography(p5+1)
  vp7 = vp_tomography(p6+1)
  vp8 = vp_tomography(p7+1)

  vs1 = vs_tomography(p0+1)
  vs2 = vs_tomography(p1+1)
  vs3 = vs_tomography(p2+1)
  vs4 = vs_tomography(p3+1)
  vs5 = vs_tomography(p4+1)
  vs6 = vs_tomography(p5+1)
  vs7 = vs_tomography(p6+1)
  vs8 = vs_tomography(p7+1)

  rho1 = rho_tomography(p0+1)
  rho2 = rho_tomography(p1+1)
  rho3 = rho_tomography(p2+1)
  rho4 = rho_tomography(p3+1)
  rho5 = rho_tomography(p4+1)
  rho6 = rho_tomography(p5+1)
  rho7 = rho_tomography(p6+1)
  rho8 = rho_tomography(p7+1)

  ! use trilinear interpolation in cell to define Vp Vs and rho
  vp_final = &
     vp1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
     vp2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
     vp3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
     vp4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
     vp5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
     vp6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
     vp7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
     vp8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
    
  vs_final = &
     vs1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
     vs2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
     vs3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
     vs4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
     vs5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
     vs6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
     vs7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
     vs8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
         
  rho_final = &
     rho1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + &
     rho2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + &
     rho3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + &
     rho4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + &
     rho5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + &
     rho6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + &
     rho7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + &
     rho8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
              
  ! impose minimum and maximum velocity and density if needed
  if(vp_final < VP_MIN) vp_final = VP_MIN
  if(vs_final < VS_MIN) vs_final = VS_MIN
  if(rho_final < RHO_MIN) rho_final = RHO_MIN
  if(vp_final > VP_MAX) vp_final = VP_MAX
  if(vs_final > VS_MAX) vs_final = VS_MAX
  if(rho_final > RHO_MAX) rho_final = RHO_MAX

  end subroutine model_tomography
