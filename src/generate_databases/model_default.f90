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

!--------------------------------------------------------------------------------------------------
!
! model given by (CUBIT) mesh parameters
!
!--------------------------------------------------------------------------------------------------

  subroutine model_default(materials_ext_mesh,nmat_ext_mesh, &
                          undef_mat_prop,nundefMat_ext_mesh, &
                          imaterial_id,imaterial_def, &
                          xmesh,ymesh,zmesh, &
                          rho,vp,vs,iflag_aniso,qkappa_atten,qmu_atten,idomain_id, &
                          rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                          phi,tort,kxx,kxy,kxz,kyy,kyz,kzz)

! takes model values specified by mesh properties

  !use generate_databases_par, only: myrank
  use create_regions_mesh_ext_par

  implicit none

  integer, intent(in) :: nmat_ext_mesh
  double precision, dimension(16,nmat_ext_mesh),intent(in) :: materials_ext_mesh

  integer, intent(in) :: nundefMat_ext_mesh
  character (len=30), dimension(6,nundefMat_ext_mesh):: undef_mat_prop

  integer, intent(in) :: imaterial_id,imaterial_def

  double precision, intent(in) :: xmesh,ymesh,zmesh

  real(kind=CUSTOM_REAL) :: vp,vs,rho,qkappa_atten,qmu_atten

  integer :: iflag_aniso
  integer :: idomain_id

  real(kind=CUSTOM_REAL) :: kappa_s,kappa_f,kappa_fr,mu_fr,rho_s,rho_f,phi,tort,eta_f, &
                           kxx,kxy,kxz,kyy,kyz,kzz

  ! local parameters
  integer :: iflag,flag_below,flag_above
  integer :: iundef

  ! check if the material is known or unknown
  if( imaterial_id > 0 ) then
    ! gets velocity model as specified by (cubit) mesh files for elastic & acoustic
    ! or from nummaterial_poroelastic_file for poroelastic (too many arguments for cubit)

    ! material domain_id
    idomain_id = nint(materials_ext_mesh(6,imaterial_id))

    select case( idomain_id )

    case( IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC )
      ! (visco)elastic or acoustic

      ! density
      ! materials_ext_mesh format:
      ! #index1 = rho #index2 = vp #index3 = vs #index4 = Q_mu #index5 = iflag_aniso   #index7 = Q_kappa
      ! Q_kappa is not stored next to Q_mu for historical reasons, because it was added later
      rho = materials_ext_mesh(1,imaterial_id)

      ! isotropic values: vp, vs
      vp = materials_ext_mesh(2,imaterial_id)
      vs = materials_ext_mesh(3,imaterial_id)

      ! attenuation
      ! Q_kappa is not stored next to Q_mu for historical reasons, because it was added later
      qkappa_atten = materials_ext_mesh(7,imaterial_id)
      qmu_atten = materials_ext_mesh(4,imaterial_id)

      ! anisotropy
      iflag_aniso = nint(materials_ext_mesh(5,imaterial_id))

    case( IDOMAIN_POROELASTIC )
      ! poroelastic
      ! materials_ext_mesh format:
      ! rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz

      ! solid properties
      rho_s =  materials_ext_mesh(1,imaterial_id)
      kappa_s =  materials_ext_mesh(13,imaterial_id)
      ! fluid properties
      rho_f =  materials_ext_mesh(2,imaterial_id)
      kappa_f =  materials_ext_mesh(14,imaterial_id)
      eta_f =  materials_ext_mesh(5,imaterial_id)
      ! frame properties
      kappa_fr =  materials_ext_mesh(15,imaterial_id)
      mu_fr =  materials_ext_mesh(16,imaterial_id)
      phi =  materials_ext_mesh(3,imaterial_id)
      tort =  materials_ext_mesh(4,imaterial_id)
      kxx =  materials_ext_mesh(7,imaterial_id)
      kxy =  materials_ext_mesh(8,imaterial_id)
      kxz =  materials_ext_mesh(9,imaterial_id)
      kyy =  materials_ext_mesh(10,imaterial_id)
      kyz =  materials_ext_mesh(11,imaterial_id)
      kzz =  materials_ext_mesh(12,imaterial_id)

    case default
      print*,'error: domain id = ',idomain_id,'not recognized'
      stop 'error: domain not recognized'

    end select

 else if ( imaterial_def == 1 ) then

    stop 'material: interface not implemented yet'

    do iundef = 1,nundefMat_ext_mesh
       if(trim(undef_mat_prop(2,iundef)) == 'interface') then
          read(undef_mat_prop(3,iundef),'(1i3)') flag_below
          read(undef_mat_prop(4,iundef),'(1i3)') flag_above
       endif
    enddo

    ! see file model_interface_bedrock.f90: routine interface()
    !call interface(iflag,flag_below,flag_above,ispec,nspec,i,j,k,xstore,ystore,zstore,ibedrock)

    ! dummy: takes 1. defined material
    iflag = 1
    rho = materials_ext_mesh(1,iflag)
    vp = materials_ext_mesh(2,iflag)
    vs = materials_ext_mesh(3,iflag)
    ! Q_kappa is not stored next to Q_mu for historical reasons, because it was added later
    qkappa_atten = materials_ext_mesh(7,iflag)
    qmu_atten = materials_ext_mesh(4,iflag)

    iflag_aniso = nint(materials_ext_mesh(5,iflag))
    idomain_id = nint(materials_ext_mesh(6,iflag))

  else if ( imaterial_def == 2 ) then

    ! material definition undefined, uses definition from tomography model

    ! gets model values from tomography file
    call model_tomography(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten,imaterial_id)

    ! no anisotropy
    iflag_aniso = 0

    ! sets acoustic/elastic domain as given in materials properties
    iundef = - imaterial_id    ! iundef must be positive
    read(undef_mat_prop(6,iundef),*) idomain_id
    ! or
    !idomain_id = IDOMAIN_ELASTIC    ! forces to be elastic domain

  else

    stop 'material: not implemented yet'

  endif

  end subroutine model_default
