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

  subroutine save_databases_hdf5(prname,nspec,nglob, &
                            iMPIcut_xi,iMPIcut_eta, &
                            ibool,nodes_coords,ispec_material_id, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            xstore, ystore, zstore)


  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC,SAVE_MESH_AS_CUBIT,NDIM,CUSTOM_REAL
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M,NUMBER_OF_MATERIAL_PROPERTIES

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,NGNOD,NGNOD2D,LOCAL_PATH

  use meshfem_par, only: NPROC, myrank, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NMATERIALS,material_properties,material_properties_undef,&
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions, &
    addressing, &
    iproc_xi_current,iproc_eta_current, &
    NPROC_XI,NPROC_ETA

  use phdf5_utils
  use my_mpi


  implicit none

  ! number of spectral elements in each block
  integer, intent(in) :: nspec

  ! number of vertices in each block
  integer, intent(in) :: nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8
  integer iproc
  !integer NPROC_XI,NPROC_ETA
  logical, intent(in) :: iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)
  !integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

  ! arrays with the mesh
  integer, intent(in) :: ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision, intent(in) :: nodes_coords(nglob,NDIM)

  !! VM VM add all GLL points for Axisem coupling
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xstore, ystore, zstore

  integer, intent(in) :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer, intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer, intent(in) :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer, intent(in) :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer, intent(in) :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer, intent(in) :: ibelm_top(NSPEC2D_TOP)

  ! material properties
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
  double precision , dimension(17,NMATERIALS) :: mat_prop


  ! CPML
  integer :: nspec_CPML_total,ispec_CPML

  integer :: i,ispec,iglob,ier

  ! name of the database files
  character(len=MAX_STRING_LEN), intent(in) :: prname

end subroutine save_databases_hdf5