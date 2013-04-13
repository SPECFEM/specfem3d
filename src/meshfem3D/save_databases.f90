!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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


  subroutine save_databases(prname,nspec,nglob,iproc_xi,iproc_eta, &
                            NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta,&
                            ibool,nodes_coords,true_material_num, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP,&
                            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top,&
                            NMATERIALS,material_properties)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  ! number of spectral elements in each block
  integer nspec

  ! number of vertices in each block
  integer nglob

  ! MPI cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8
  integer iproc_xi,iproc_eta
  integer NPROC_XI,NPROC_ETA
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)
  integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

  ! arrays with the mesh
  integer ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision :: nodes_coords(nglob,3)

  integer true_material_num(nspec)

  ! boundary parameters locator
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM)
  integer ibelm_top(NSPEC2D_TOP)

  ! material properties
  integer :: NMATERIALS
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
  double precision , dimension(NMATERIALS,6) ::  material_properties
  double precision , dimension(16) :: matpropl
  integer :: i,ispec,iglob,ier
  ! dummy_nspec_cpml is used here to match the read instructions in generate_databases/read_partition_files.f90
  integer :: dummy_nspec_cpml

  ! name of the database files
  character(len=256) prname

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max,idoubl
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  integer, parameter :: IIN_database = 15

  ! opens database file
  open(unit=IIN_database,file=prname(1:len_trim(prname))//'Database', &
        status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening Database file'

  write(IIN_database) nglob
  do iglob=1,nglob
     write(IIN_database) iglob,nodes_coords(iglob,1),nodes_coords(iglob,2),nodes_coords(iglob,3)
  enddo

  ! Materials properties
   write(IIN_database) NMATERIALS, 0
   do idoubl = 1,NMATERIALS
      !write(IIN_database,*) material_properties(idoubl,:)
      matpropl(:) = 0.d0
      matpropl(1:6) = material_properties(idoubl,1:6)
      ! pad dummy zeros to fill up 16 entries (poroelastic medium not allowed)
      write(IIN_database) matpropl
   enddo


  write(IIN_database) nspec
  do ispec=1,nspec
      !write(IIN_database,'(11i14)') ispec,true_material_num(ispec),1,ibool(1,1,1,ispec),ibool(2,1,1,ispec),&
      !     ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(1,1,2,ispec),&
      !     ibool(2,1,2,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)
      write(IIN_database) ispec,true_material_num(ispec),1,ibool(1,1,1,ispec),ibool(2,1,1,ispec),&
           ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(1,1,2,ispec),&
           ibool(2,1,2,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)
  enddo

  ! Boundaries
  write(IIN_database) 1,nspec2D_xmin
  write(IIN_database) 2,nspec2D_xmax
  write(IIN_database) 3,nspec2D_ymin
  write(IIN_database) 4,nspec2D_ymax
  write(IIN_database) 5,NSPEC2D_BOTTOM
  write(IIN_database) 6,NSPEC2D_TOP

  do i=1,nspec2D_xmin
     write(IIN_database) ibelm_xmin(i),ibool(1,1,1,ibelm_xmin(i)),ibool(1,NGLLY_M,1,ibelm_xmin(i)),&
          ibool(1,1,NGLLZ_M,ibelm_xmin(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
  enddo
  do i=1,nspec2D_xmax
     write(IIN_database) ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), &
          ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
  enddo
  do i=1,nspec2D_ymin
     write(IIN_database) ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)),&
          ibool(1,1,NGLLZ_M,ibelm_ymin(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
  enddo
  do i=1,nspec2D_ymax
     write(IIN_database) ibelm_ymax(i),ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i)),ibool(1,NGLLY_M,1,ibelm_ymax(i)), &
          ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
  enddo
  do i=1,NSPEC2D_BOTTOM
     write(IIN_database) ibelm_bottom(i),ibool(1,1,1,ibelm_bottom(i)),ibool(NGLLX_M,1,1,ibelm_bottom(i)), &
          ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i)),ibool(1,NGLLY_M,1,ibelm_bottom(i))
  enddo
  do i=1,NSPEC2D_TOP
     write(IIN_database) ibelm_top(i),ibool(1,1,NGLLZ_M,ibelm_top(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i)), &
          ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
  enddo

  ! JC JC todo: implement C-PML code in internal mesher
  ! dummy_nspec_cpml is used here to match the read instructions in generate_databases/read_partition_files.f90
  dummy_nspec_cpml = 0
  write(IIN_database) dummy_nspec_cpml

  ! MPI Interfaces

  if(NPROC_XI >= 2 .or. NPROC_ETA >= 2) then

  nb_interfaces = 4
  interfaces(W:N) = .true.
  interfaces(NW:SW) = .false.
  if(iproc_xi == 0) then
     nb_interfaces =  nb_interfaces -1
     interfaces(W) = .false.
  endif
  if(iproc_xi == NPROC_XI-1) then
     nb_interfaces =  nb_interfaces -1
     interfaces(E) = .false.
  endif
  if(iproc_eta == 0) then
     nb_interfaces =  nb_interfaces -1
     interfaces(S) = .false.
  endif
  if(iproc_eta == NPROC_ETA-1) then
     nb_interfaces =  nb_interfaces -1
     interfaces(N) = .false.
  endif

  if((interfaces(W) .eqv. .true.) .and. (interfaces(N) .eqv. .true.)) then
       interfaces(NW) = .true.
       nb_interfaces =  nb_interfaces +1
  endif
  if((interfaces(N) .eqv. .true.) .and. (interfaces(E) .eqv. .true.)) then
       interfaces(NE) = .true.
       nb_interfaces =  nb_interfaces +1
  endif
  if((interfaces(E) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
       interfaces(SE) = .true.
       nb_interfaces =  nb_interfaces +1
  endif
  if((interfaces(W) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
       interfaces(SW) = .true.
       nb_interfaces =  nb_interfaces +1
  endif

  nspec_interface(:) = 0
  if(interfaces(W))  nspec_interface(W) = count(iMPIcut_xi(1,:) .eqv. .true.)
  if(interfaces(E))  nspec_interface(E) = count(iMPIcut_xi(2,:) .eqv. .true.)
  if(interfaces(S))  nspec_interface(S) = count(iMPIcut_eta(1,:) .eqv. .true.)
  if(interfaces(N))  nspec_interface(N) = count(iMPIcut_eta(2,:) .eqv. .true.)
  if(interfaces(NW))  nspec_interface(NW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
  if(interfaces(NE))  nspec_interface(NE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
  if(interfaces(SE))  nspec_interface(SE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))
  if(interfaces(SW))  nspec_interface(SW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))

  nspec_interfaces_max = maxval(nspec_interface)

  write(IIN_database) nb_interfaces,nspec_interfaces_max

  if(interfaces(W)) then
     write(IIN_database) addressing(iproc_xi-1,iproc_eta),nspec_interface(W)
     do ispec = 1,nspec
        if(iMPIcut_xi(1,ispec))  write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(1,2,1,ispec), &
             ibool(1,1,2,ispec),ibool(1,2,2,ispec)
     enddo
  endif

  if(interfaces(E)) then
     write(IIN_database) addressing(iproc_xi+1,iproc_eta),nspec_interface(E)
     do ispec = 1,nspec
        if(iMPIcut_xi(2,ispec))  write(IIN_database) ispec,4,ibool(2,1,1,ispec),ibool(2,2,1,ispec), &
             ibool(2,1,2,ispec),ibool(2,2,2,ispec)
     enddo
  endif

   if(interfaces(S)) then
     write(IIN_database) addressing(iproc_xi,iproc_eta-1),nspec_interface(S)
     do ispec = 1,nspec
        if(iMPIcut_eta(1,ispec))  write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(2,1,1,ispec), &
             ibool(1,1,2,ispec),ibool(2,1,2,ispec)
     enddo
  endif

  if(interfaces(N)) then
     write(IIN_database) addressing(iproc_xi,iproc_eta+1),nspec_interface(N)
     do ispec = 1,nspec
        if(iMPIcut_eta(2,ispec))  write(IIN_database) ispec,4,ibool(2,2,1,ispec),ibool(1,2,1,ispec), &
             ibool(2,2,2,ispec),ibool(1,2,2,ispec)
     enddo
  endif

  if(interfaces(NW)) then
     write(IIN_database) addressing(iproc_xi-1,iproc_eta+1),nspec_interface(NW)
     do ispec = 1,nspec
        if((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.))  then
           write(IIN_database) ispec,2,ibool(1,2,1,ispec),ibool(1,2,2,ispec),-1,-1
        endif
     enddo
  endif

  if(interfaces(NE)) then
     write(IIN_database) addressing(iproc_xi+1,iproc_eta+1),nspec_interface(NE)
     do ispec = 1,nspec
        if((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.))  then
           write(IIN_database) ispec,2,ibool(2,2,1,ispec),ibool(2,2,2,ispec),-1,-1
        endif
     enddo
  endif

  if(interfaces(SE)) then
     write(IIN_database) addressing(iproc_xi+1,iproc_eta-1),nspec_interface(SE)
     do ispec = 1,nspec
        if((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.))  then
           write(IIN_database) ispec,2,ibool(2,1,1,ispec),ibool(2,1,2,ispec),-1,-1
        endif
     enddo
  endif

  if(interfaces(SW)) then
     write(IIN_database) addressing(iproc_xi-1,iproc_eta-1),nspec_interface(SW)
     do ispec = 1,nspec
        if((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.))  then
           write(IIN_database) ispec,2,ibool(1,1,1,ispec),ibool(1,1,2,ispec),-1,-1
        endif
     enddo
  endif

  else

     write(IIN_database) 0,0

  endif

  close(IIN_database)

  end subroutine save_databases
