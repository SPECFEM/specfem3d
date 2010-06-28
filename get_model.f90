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

  subroutine get_model(myrank,nspec,ibool,mat_ext_mesh,nelmnts_ext_mesh, &
                        materials_ext_mesh,nmat_ext_mesh, &
                        undef_mat_prop,nundefMat_ext_mesh, &
                        ANISOTROPY)

  use create_regions_mesh_ext_par
  implicit none

  ! number of spectral elements in each block
  integer :: myrank,nspec

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! external mesh
  integer :: nelmnts_ext_mesh
  integer :: nmat_ext_mesh,nundefMat_ext_mesh

  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh
  double precision, dimension(6,nmat_ext_mesh) :: materials_ext_mesh
  character (len=30), dimension(6,nundefMat_ext_mesh):: undef_mat_prop

  ! anisotropy
  logical :: ANISOTROPY

  ! local parameters
  real(kind=CUSTOM_REAL) :: vp,vs,rho
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                        c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  integer :: ispec,i,j,k,iundef,iflag_atten
  integer :: iflag,flag_below,flag_above
  integer :: iflag_aniso,idomain_id,imaterial_id

  ! gll point location
  double precision :: xloc,yloc,zloc
  integer :: iglob

  ! initializes element domain flags
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! prepares tomography model if needed for elements with undefined material definitions
  if( nundefMat_ext_mesh > 0 ) then
    call model_tomography_broadcast(myrank)
  endif

  ! prepares external model values if needed
  if( USE_MODEL_EXTERNAL_VALUES ) then
    call model_external_broadcast(myrank)
  endif

! !  Piero, read bedrock file
! in case, see file model_interface_bedrock.f90:
!  call model_bedrock_broadcast(myrank)


  ! material properties on all GLL points: taken from material values defined for
  ! each spectral element in input mesh
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

           ! material index 1: associated material number
           imaterial_id = mat_ext_mesh(1,ispec)

           ! check if the material is known or unknown
           if( imaterial_id > 0) then
              ! gets velocity model as specified by (cubit) mesh files

              ! density
              ! materials_ext_mesh format:
              ! #index1 = rho #index2 = vp #index3 = vs #index4 = Q_flag #index5 = 0
              rho = materials_ext_mesh(1,imaterial_id)

              ! isotropic values: vp, vs
              vp = materials_ext_mesh(2,imaterial_id)
              vs = materials_ext_mesh(3,imaterial_id)

              ! attenuation
              iflag_atten = materials_ext_mesh(4,imaterial_id)
              !change for piero :
              !if(mat_ext_mesh(1,ispec) == 1) then
              !   iflag_attenuation_store(i,j,k,ispec) = 1
              !else
              !   iflag_attenuation_store(i,j,k,ispec) = 2
              !endif

              ! anisotropy
              iflag_aniso = materials_ext_mesh(5,imaterial_id)

              ! material domain_id
              idomain_id = materials_ext_mesh(6,imaterial_id)

           else if (mat_ext_mesh(2,ispec) == 1) then

              stop 'material: interface not implemented yet'

              do iundef = 1,nundefMat_ext_mesh
                 if(trim(undef_mat_prop(2,iundef)) == 'interface') then
                    read(undef_mat_prop(4,iundef),'(1i3)') flag_below
                    read(undef_mat_prop(5,iundef),'(1i3)') flag_above
                 endif
              enddo

              ! see file model_interface_bedrock.f90: routine interface()
              !call interface(iflag,flag_below,flag_above,ispec,nspec,i,j,k,xstore,ystore,zstore,ibedrock)

              iflag = 1
              rho = materials_ext_mesh(1,iflag)
              vp = materials_ext_mesh(2,iflag)
              vs = materials_ext_mesh(3,iflag)
              iflag_atten = materials_ext_mesh(4,iflag)
              !change for piero :
              !  if(iflag == 1) then
              !     iflag_attenuation_store(i,j,k,ispec) = 1
              !  else
              !     iflag_attenuation_store(i,j,k,ispec) = 2
              !  endif
              iflag_aniso = materials_ext_mesh(5,iflag)
              idomain_id = materials_ext_mesh(6,iflag)

           else if ( imaterial_id < 0 ) then
           
              ! material definition undefined, uses definition from tomography model
              ! GLL point location
              iglob = ibool(i,j,k,ispec)
              xloc = xstore_dummy(iglob)
              yloc = ystore_dummy(iglob)
              zloc = zstore_dummy(iglob)

              ! gets model values from tomography file
              call model_tomography(xloc,yloc,zloc, &
                                  rho,vp,vs)

              iflag_atten = 1   ! attenuation: would use IATTENUATION_SEDIMENTS_40
              iflag_aniso = 0   ! no anisotropy
              idomain_id = 2    ! elastic domain

           else

              stop 'material: not implemented yet'

           end if

           ! adds/gets velocity model as specified in model_external_values.f90
           if( USE_MODEL_EXTERNAL_VALUES ) then
             call model_external_values(i,j,k,ispec,idomain_id,imaterial_id, &
                            nspec,ibool, &
                            iflag_aniso,iflag_atten, &
                            rho,vp,vs, &
                            c11,c12,c13,c14,c15,c16, &
                            c22,c23,c24,c25,c26,c33, &
                            c34,c35,c36,c44,c45,c46, &
                            c55,c56,c66,ANISOTROPY)
           endif

           ! adds anisotropic default model
           if( ANISOTROPY .and. .not. USE_MODEL_EXTERNAL_VALUES ) then
             call model_aniso(iflag_aniso,rho,vp,vs,c11,c12,c13,c14,c15,c16, &
                     c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45, &
                     c46,c55,c56,c66)

           endif

           ! stores velocity model

           ! density
           rhostore(i,j,k,ispec) = rho

           ! kappa, mu
           kappastore(i,j,k,ispec) = rho*( vp*vp - FOUR_THIRDS*vs*vs )
           mustore(i,j,k,ispec) = rho*vs*vs

           ! attenuation
           iflag_attenuation_store(i,j,k,ispec) = iflag_atten

           ! Stacey, a completer par la suite
           rho_vp(i,j,k,ispec) = rho*vp
           rho_vs(i,j,k,ispec) = rho*vs
           !end pll

           ! adds anisotropic perturbation to vp, vs
           if( ANISOTROPY ) then
             !call model_aniso(iflag_aniso,rho,vp,vs,c11,c12,c13,c14,c15,c16, &
             !        c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
             c11store(i,j,k,ispec) = c11
             c12store(i,j,k,ispec) = c12
             c13store(i,j,k,ispec) = c13
             c14store(i,j,k,ispec) = c14
             c15store(i,j,k,ispec) = c15
             c16store(i,j,k,ispec) = c16
             c22store(i,j,k,ispec) = c22
             c23store(i,j,k,ispec) = c23
             c24store(i,j,k,ispec) = c24
             c25store(i,j,k,ispec) = c25
             c26store(i,j,k,ispec) = c26
             c33store(i,j,k,ispec) = c33
             c34store(i,j,k,ispec) = c34
             c35store(i,j,k,ispec) = c35
             c36store(i,j,k,ispec) = c36
             c44store(i,j,k,ispec) = c44
             c45store(i,j,k,ispec) = c45
             c46store(i,j,k,ispec) = c46
             c55store(i,j,k,ispec) = c55
             c56store(i,j,k,ispec) = c56
             c66store(i,j,k,ispec) = c66
           endif

           ! material domain
           !print*,'velocity model:',ispec,idomain_id
           if( idomain_id == IDOMAIN_ACOUSTIC ) then
             ispec_is_acoustic(ispec) = .true.
           else if( idomain_id == IDOMAIN_ELASTIC ) then
             ispec_is_elastic(ispec) = .true.
           else if( idomain_id == IDOMAIN_POROELASTIC ) then
             stop 'poroelastic material domain not implemented yet'
             ispec_is_poroelastic(ispec) = .true.
           else
             stop 'error material domain index'
           endif

        enddo
      enddo
    enddo
    !print*,myrank,'ispec:',ispec,'rho:',rhostore(1,1,1,ispec),'vp:',vpstore(1,1,1,ispec),'vs:',vsstore(1,1,1,ispec)
  enddo

  ! checks material domains
  do ispec=1,nspec
    if( (ispec_is_acoustic(ispec) .eqv. .false.) &
          .and. (ispec_is_elastic(ispec) .eqv. .false.) &
          .and. (ispec_is_poroelastic(ispec) .eqv. .false.) ) then
      print*,'error material domain not assigned to element:',ispec
      print*,'acoustic: ',ispec_is_acoustic(ispec)
      print*,'elastic: ',ispec_is_elastic(ispec)
      print*,'poroelastic: ',ispec_is_poroelastic(ispec)
      stop 'error material domain index element'
    endif
  enddo


! !! DK DK store the position of the six stations to be able to
! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
! in case, see file model_interface_bedrock.f90: routine model_bedrock_store()

  end subroutine get_model

