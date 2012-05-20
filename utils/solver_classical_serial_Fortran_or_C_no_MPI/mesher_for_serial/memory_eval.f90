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

! compute the approximate amount of static memory needed to run the solver

  subroutine memory_eval(OCEANS,ABSORBING_CONDITIONS,ATTENUATION,ANISOTROPIC_3D_MANTLE,&
                       TRANSVERSE_ISOTROPY,ANISOTROPIC_INNER_CORE,ROTATION,&
                       ONE_CRUST,doubling_index,this_region_has_a_doubling,&
                       ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_sampling_array,&
                       NSPEC,nglob,SIMULATION_TYPE,MOVIE_VOLUME,SAVE_FORWARD, &
         NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,static_memory_size)

  implicit none

  include "constants.h"

! input
  logical, intent(in) :: TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
             ROTATION,ATTENUATION,ONE_CRUST,OCEANS,ABSORBING_CONDITIONS,MOVIE_VOLUME,SAVE_FORWARD
  integer, dimension(MAX_NUM_REGIONS), intent(in) :: NSPEC, nglob
  integer, intent(in) :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA,SIMULATION_TYPE
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS), intent(in) :: doubling_index
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS), intent(in) :: this_region_has_a_doubling
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS), intent(in) :: ner,ratio_sampling_array

! output
  double precision, intent(out) :: static_memory_size

! variables
  integer :: ilayer,NUMBER_OF_MESH_LAYERS,ner_without_doubling,ispec_aniso

  integer, intent(out) :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

! generate the elements in all the regions of the mesh
  ispec_aniso = 0

  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif

! count anisotropic elements
  do ilayer = 1, NUMBER_OF_MESH_LAYERS
      if(doubling_index(ilayer) == IFLAG_220_80 .or. doubling_index(ilayer) == IFLAG_80_MOHO) then
          ner_without_doubling = ner(ilayer)
          if(this_region_has_a_doubling(ilayer)) then
              ner_without_doubling = ner_without_doubling - 2
              ispec_aniso = ispec_aniso + &
              (NSPEC_DOUBLING_SUPERBRICK*(NEX_PER_PROC_XI/ratio_sampling_array(ilayer)/2)* &
              (NEX_PER_PROC_ETA/ratio_sampling_array(ilayer)/2))
          endif
          ispec_aniso = ispec_aniso + &
          ((NEX_PER_PROC_XI/ratio_sampling_array(ilayer))*(NEX_PER_PROC_ETA/ratio_sampling_array(ilayer))*ner_without_doubling)
      endif
  enddo

! define static size of the arrays whose size depends on logical tests

  if(ANISOTROPIC_INNER_CORE) then
    NSPECMAX_ANISO_IC = NSPEC(IREGION_INNER_CORE)
  else
    NSPECMAX_ANISO_IC = 1
  endif

  if(ANISOTROPIC_3D_MANTLE) then
    NSPECMAX_ISO_MANTLE = 1
    NSPECMAX_TISO_MANTLE = 1
    NSPECMAX_ANISO_MANTLE = NSPEC(IREGION_CRUST_MANTLE)
  else

    NSPECMAX_ISO_MANTLE = NSPEC(IREGION_CRUST_MANTLE)
    if(TRANSVERSE_ISOTROPY) then
      NSPECMAX_TISO_MANTLE = ispec_aniso
    else
      NSPECMAX_TISO_MANTLE = 1
    endif

    NSPECMAX_ANISO_MANTLE = 1
  endif

! if attenuation is off, set dummy size of arrays to one
  if(ATTENUATION) then
    NSPEC_CRUST_MANTLE_ATTENUAT = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_ATTENUATION = NSPEC(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_ATTENUAT = 1
    NSPEC_INNER_CORE_ATTENUATION = 1
  endif

  if(ATTENUATION .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    NSPEC_CRUST_MANTLE_STR_OR_ATT = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STR_OR_ATT = NSPEC(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STR_OR_ATT = 1
    NSPEC_INNER_CORE_STR_OR_ATT = 1
  endif

  if(ATTENUATION .and. SIMULATION_TYPE == 3) then
    NSPEC_CRUST_MANTLE_STR_AND_ATT = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STR_AND_ATT = NSPEC(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STR_AND_ATT = 1
    NSPEC_INNER_CORE_STR_AND_ATT = 1
  endif


  if(SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STRAIN_ONLY = NSPEC(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = 1
    NSPEC_INNER_CORE_STRAIN_ONLY = 1
  endif

  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    NSPEC_CRUST_MANTLE_ADJOINT = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_OUTER_CORE_ADJOINT = NSPEC(IREGION_OUTER_CORE)
    NSPEC_INNER_CORE_ADJOINT = NSPEC(IREGION_INNER_CORE)

    NGLOB_CRUST_MANTLE_ADJOINT = NGLOB(IREGION_CRUST_MANTLE)
    NGLOB_OUTER_CORE_ADJOINT = NGLOB(IREGION_OUTER_CORE)
    NGLOB_INNER_CORE_ADJOINT = NGLOB(IREGION_INNER_CORE)

    if(ROTATION) then
      NSPEC_OUTER_CORE_ROT_ADJOINT = NSPEC(IREGION_OUTER_CORE)
    else
      NSPEC_OUTER_CORE_ROT_ADJOINT = 1
    endif
  else
    NSPEC_CRUST_MANTLE_ADJOINT = 1
    NSPEC_OUTER_CORE_ADJOINT = 1
    NSPEC_INNER_CORE_ADJOINT = 1

    NGLOB_CRUST_MANTLE_ADJOINT = 1
    NGLOB_OUTER_CORE_ADJOINT = 1
    NGLOB_INNER_CORE_ADJOINT = 1

    NSPEC_OUTER_CORE_ROT_ADJOINT = 1
   endif

! if absorbing conditions are off, set dummy size of arrays to one
  if(ABSORBING_CONDITIONS) then
    NSPEC_CRUST_MANTLE_STACEY = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_OUTER_CORE_STACEY = NSPEC(IREGION_OUTER_CORE)
  else
    NSPEC_CRUST_MANTLE_STACEY = 1
    NSPEC_OUTER_CORE_STACEY = 1
  endif

! if oceans are off, set dummy size of arrays to one
  if(OCEANS) then
    NGLOB_CRUST_MANTLE_OCEANS = NGLOB(IREGION_CRUST_MANTLE)
  else
    NGLOB_CRUST_MANTLE_OCEANS = 1
  endif

  if(ROTATION) then
    NSPEC_OUTER_CORE_ROTATION = NSPEC(IREGION_OUTER_CORE)
  else
    NSPEC_OUTER_CORE_ROTATION = 1
  endif

! add size of each set of static arrays multiplied by the number of such arrays

  static_memory_size = 0.d0

! R_memory_crust_mantle
! static_memory_size = static_memory_size + 5.d0*dble(N_SLS)*dble(NGLLX)* &
!   dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ATTENUAT*dble(CUSTOM_REAL)

! R_memory_inner_core
! static_memory_size = static_memory_size + 5.d0*dble(N_SLS)*dble(NGLLX)* &
!   dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_ATTENUATION*dble(CUSTOM_REAL)

! xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle
! etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,
! gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_CRUST_MANTLE)*9.d0*dble(CUSTOM_REAL)

! ibool_crust_mantle
  static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_CRUST_MANTLE)*dble(SIZE_INTEGER)

! xix_outer_core,xiy_outer_core,xiz_outer_core,
! etax_outer_core,etay_outer_core,etaz_outer_core,
! gammax_outer_core,gammay_outer_core,gammaz_outer_core
! rhostore_outer_core,kappavstore_outer_core
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_OUTER_CORE)*11.d0*dble(CUSTOM_REAL)

! ibool_outer_core
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_OUTER_CORE)*dble(SIZE_INTEGER)

! idoubling_crust_mantle
! static_memory_size = static_memory_size + NSPEC(IREGION_CRUST_MANTLE)*dble(SIZE_INTEGER)

!!!!!!!!!!!!! xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,rmass_crust_mantle
!!!!!!!!!!!!! static_memory_size = static_memory_size + NGLOB(IREGION_CRUST_MANTLE)*4.d0*dble(CUSTOM_REAL)
! rmass_crust_mantle
  static_memory_size = static_memory_size + NGLOB(IREGION_CRUST_MANTLE)*1.d0*dble(CUSTOM_REAL)

!!!!!!!!! rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle
!!!!!!!!  static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ISO_MANTLE*3.d0*dble(CUSTOM_REAL)
! kappavstore_crust_mantle,muvstore_crust_mantle
  static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ISO_MANTLE*2.d0*dble(CUSTOM_REAL)

! kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_TISO_MANTLE*3.d0*dble(CUSTOM_REAL)

! c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,
! c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle,
! c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle,
! c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle,
! c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle,
! c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle,
! c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ANISO_MANTLE*21.d0*dble(CUSTOM_REAL)

! displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
  static_memory_size = static_memory_size + dble(NDIM)*NGLOB(IREGION_CRUST_MANTLE)*3.d0*dble(CUSTOM_REAL)

! xstore_outer_core, ystore_outer_core, zstore_outer_core, rmass_outer_core, displ_outer_core, veloc_outer_core, accel_outer_core
! static_memory_size = static_memory_size + NGLOB(IREGION_OUTER_CORE)*7.d0*dble(CUSTOM_REAL)

! ibool_inner_core
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_INNER_CORE)*dble(SIZE_INTEGER)

! xix_inner_core,xiy_inner_core,xiz_inner_core,
! etax_inner_core,etay_inner_core,etaz_inner_core,
! gammax_inner_core,gammay_inner_core,gammaz_inner_core,
! rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_INNER_CORE)*12.d0*dble(CUSTOM_REAL)

! xstore_inner_core,ystore_inner_core,zstore_inner_core,rmass_inner_core
! static_memory_size = static_memory_size + NGLOB(IREGION_INNER_CORE)*4.d0*dble(CUSTOM_REAL)

! c11store_inner_core,c33store_inner_core,c12store_inner_core,c13store_inner_core,c44store_inner_core
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ANISO_IC*5.d0*dble(CUSTOM_REAL)

! displ_inner_core,veloc_inner_core,accel_inner_core
! static_memory_size = static_memory_size + dble(NDIM)*NGLOB(IREGION_INNER_CORE)*3.d0*dble(CUSTOM_REAL)

! A_array_rotation,B_array_rotation
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ROTATION*2.d0*dble(CUSTOM_REAL)

! if(ABSORBING_CONDITIONS) then

! rho_vp_crust_mantle,rho_vs_crust_mantle
!   static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_CRUST_MANTLE)*2.d0*dble(CUSTOM_REAL)

! vp_outer_core
!   static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

! endif

! if(OCEANS) then

! rmass_ocean_load
!   static_memory_size = static_memory_size + NGLOB(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

! updated_dof_ocean_load
!   static_memory_size = static_memory_size + NGLOB(IREGION_CRUST_MANTLE)*dble(SIZE_LOGICAL)

! endif

! add arrays used to save strain for attenuation or for adjoint runs

! epsilondev_crust_mantle
! static_memory_size = static_memory_size + 5.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_STR_OR_ATT*dble(CUSTOM_REAL)

! eps_trace_over_3_crust_mantle
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_STR_OR_ATT*dble(CUSTOM_REAL)

! epsilondev_inner_core
! static_memory_size = static_memory_size + 5.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_STR_OR_ATT*dble(CUSTOM_REAL)

! eps_trace_over_3_inner_core
! static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_STR_OR_ATT*dble(CUSTOM_REAL)

! add arrays used for adjoint runs only (LQY: not very accurate)

! b_R_memory_crust_mantle
! b_epsilondev_crust_mantle
! b_eps_trace_over_3_crust_mantle
! rho_kl_crust_mantle,beta_kl_crust_mantle, alpha_kl_crust_mantle
! static_memory_size = static_memory_size + (5.d0*dble(N_SLS) + 9.d0)* &
!     dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ADJOINT*dble(CUSTOM_REAL)

! b_div_displ_outer_core
! rho_kl_outer_core,alpha_kl_outer_core
! static_memory_size = static_memory_size + 3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ADJOINT*dble(CUSTOM_REAL)

! b_R_memory_inner_core
! b_epsilondev_inner_core
! b_eps_trace_over_3_inner_core
! rho_kl_inner_core,beta_kl_inner_core, alpha_kl_inner_core
! static_memory_size = static_memory_size + (5.d0*dble(N_SLS) + 9.d0)* &
!     dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_ADJOINT*dble(CUSTOM_REAL)

! b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle
! static_memory_size = static_memory_size + 3.d0*dble(NDIM)*NGLOB_CRUST_MANTLE_ADJOINT*dble(CUSTOM_REAL)

! b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core
! static_memory_size = static_memory_size + 3.d0*NGLOB_OUTER_CORE_ADJOINT*dble(CUSTOM_REAL)

! b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core
! static_memory_size = static_memory_size + 3.d0*dble(NDIM)*NGLOB_INNER_CORE_ADJOINT*dble(CUSTOM_REAL)

! b_A_array_rotation,b_B_array_rotation
! static_memory_size = static_memory_size + 2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ROT_ADJOINT*dble(CUSTOM_REAL)

  end subroutine memory_eval

