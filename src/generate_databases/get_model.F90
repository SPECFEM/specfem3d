!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine get_model(myrank)

  use generate_databases_par, only: IMODEL, &
    IMODEL_DEFAULT,IMODEL_GLL,IMODEL_1D_PREM,IMODEL_1D_CASCADIA,IMODEL_1D_SOCAL, &
    IMODEL_SALTON_TROUGH,IMODEL_TOMO,IMODEL_USER_EXTERNAL,IMODEL_IPATI,IMODEL_IPATI_WATER, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
    nspec => NSPEC_AB,ibool,mat_ext_mesh, &
    materials_ext_mesh,nmat_ext_mesh,undef_mat_prop,nundefMat_ext_mesh, &
    ANISOTROPY,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,TWO,IMAIN

  use create_regions_mesh_ext_par

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_4.F90"
#endif

  implicit none

  integer :: myrank

  ! local parameters
  real(kind=CUSTOM_REAL) :: vp,vs,rho,qkappa_atten,qmu_atten
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                        c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  real(kind=CUSTOM_REAL) :: kappa_s,kappa_f,kappa_fr,mu_fr,rho_s,rho_f,phi,tort,eta_f, &
                        kxx,kxy,kxz,kyy,kyz,kzz,rho_bar
  real(kind=CUSTOM_REAL) :: cpIsquare,cpIIsquare,cssquare,H_biot,M_biot,C_biot,D_biot, &
                        afactor,bfactor,cfactor

  integer :: ispec,i,j,k

  ! material domain
  integer :: idomain_id

  integer :: imaterial_id,imaterial_def

  ! GLL point location
  double precision :: xmesh,ymesh,zmesh
  integer :: iglob

  ! timing
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! initializes element domain flags
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! prepares tomographic models if needed for elements with undefined material definitions
#ifndef DEBUG_COUPLED
  if (nundefMat_ext_mesh > 0 .or. IMODEL == IMODEL_TOMO) call model_tomography_broadcast(myrank)
#endif

  ! prepares external model values if needed
  select case (IMODEL)
  case (IMODEL_USER_EXTERNAL )
    call model_external_broadcast(myrank)
  case (IMODEL_SALTON_TROUGH )
    call model_salton_trough_broadcast(myrank)
  end select

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_1.F90"
#endif

  ! get MPI starting time
  time_start = wtime()

  ! material properties on all GLL points: taken from material values defined for
  ! each spectral element in input mesh
  do ispec = 1, nspec

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_2.F90"
#endif

    ! loops over all GLL points in element
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initializes material
          vp = 0._CUSTOM_REAL
          vs = 0._CUSTOM_REAL
          rho = 0._CUSTOM_REAL

          rho_s = 0._CUSTOM_REAL
          kappa_s = 0._CUSTOM_REAL
          rho_f = 0._CUSTOM_REAL
          kappa_f = 0._CUSTOM_REAL
          eta_f = 0._CUSTOM_REAL
          kappa_fr = 0._CUSTOM_REAL
          mu_fr = 0._CUSTOM_REAL
          phi = 0._CUSTOM_REAL
          tort = 0._CUSTOM_REAL
          kxx = 0._CUSTOM_REAL
          kxy = 0._CUSTOM_REAL
          kxz = 0._CUSTOM_REAL
          kyy = 0._CUSTOM_REAL
          kyz = 0._CUSTOM_REAL
          kzz = 0._CUSTOM_REAL

          qkappa_atten = 9999._CUSTOM_REAL
          qmu_atten = 9999._CUSTOM_REAL

          c11 = 0._CUSTOM_REAL
          c12 = 0._CUSTOM_REAL
          c13 = 0._CUSTOM_REAL
          c14 = 0._CUSTOM_REAL
          c15 = 0._CUSTOM_REAL
          c16 = 0._CUSTOM_REAL
          c22 = 0._CUSTOM_REAL
          c23 = 0._CUSTOM_REAL
          c24 = 0._CUSTOM_REAL
          c25 = 0._CUSTOM_REAL
          c26 = 0._CUSTOM_REAL
          c33 = 0._CUSTOM_REAL
          c34 = 0._CUSTOM_REAL
          c35 = 0._CUSTOM_REAL
          c36 = 0._CUSTOM_REAL
          c44 = 0._CUSTOM_REAL
          c45 = 0._CUSTOM_REAL
          c46 = 0._CUSTOM_REAL
          c55 = 0._CUSTOM_REAL
          c56 = 0._CUSTOM_REAL
          c66 = 0._CUSTOM_REAL

          ! gets xyz coordinates of GLL point
          iglob = ibool(i,j,k,ispec)
          xmesh = xstore_dummy(iglob)
          ymesh = ystore_dummy(iglob)
          zmesh = zstore_dummy(iglob)

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_3.F90"
#endif

          ! material index 1: associated material number
          ! 1 = acoustic, 2 = elastic, 3 = poroelastic, -1 = undefined tomographic
          imaterial_id = mat_ext_mesh(1,ispec)

          ! material index 2: associated material definition
          ! 1 = interface, 2 = tomography material
          imaterial_def = mat_ext_mesh(2,ispec)

          ! assigns material properties
          call get_model_values(materials_ext_mesh,nmat_ext_mesh, &
                               undef_mat_prop,nundefMat_ext_mesh, &
                               imaterial_id,imaterial_def, &
                               xmesh,ymesh,zmesh, &
                               rho,vp,vs,qkappa_atten,qmu_atten,idomain_id, &
                               rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                               phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                               c11,c12,c13,c14,c15,c16, &
                               c22,c23,c24,c25,c26,c33, &
                               c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                               ANISOTROPY)

          ! stores velocity model
          if (idomain_id == IDOMAIN_ACOUSTIC .or. idomain_id == IDOMAIN_ELASTIC) then

            ! elastic or acoustic material

            ! density
            rhostore(i,j,k,ispec) = rho

            ! kappa, mu
            kappastore(i,j,k,ispec) = rho*( vp*vp - FOUR_THIRDS*vs*vs )
            mustore(i,j,k,ispec) = rho*vs*vs

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_5.F90"
#endif

            ! attenuation
            qkappa_attenuation_store(i,j,k,ispec) = qkappa_atten
            qmu_attenuation_store(i,j,k,ispec) = qmu_atten

            ! Stacey, a completer par la suite
            rho_vp(i,j,k,ispec) = rho*vp
            rho_vs(i,j,k,ispec) = rho*vs
            !
            rho_vpI(i,j,k,ispec) = rho*vp
            rho_vpII(i,j,k,ispec) = 0.d0
            rho_vsI(i,j,k,ispec) = rho*vs
            rhoarraystore(1,i,j,k,ispec) = rho
            rhoarraystore(2,i,j,k,ispec) = rho
            phistore(i,j,k,ispec) = 0.d0
            tortstore(i,j,k,ispec) = 1.d0
            !end pll

          else

            ! poroelastic material

            ! solid properties
            rhoarraystore(1,i,j,k,ispec) = rho_s
            kappaarraystore(1,i,j,k,ispec) = kappa_s
            ! fluid properties
            rhoarraystore(2,i,j,k,ispec) = rho_f
            kappaarraystore(2,i,j,k,ispec) = kappa_f
            etastore(i,j,k,ispec) = eta_f
            ! frame properties
            kappaarraystore(3,i,j,k,ispec) = kappa_fr
            mustore(i,j,k,ispec) = mu_fr
            phistore(i,j,k,ispec) = phi
            tortstore(i,j,k,ispec) = tort
            permstore(1,i,j,k,ispec) = kxx
            permstore(2,i,j,k,ispec) = kxy
            permstore(3,i,j,k,ispec) = kxz
            permstore(4,i,j,k,ispec) = kyy
            permstore(5,i,j,k,ispec) = kyz
            permstore(6,i,j,k,ispec) = kzz

            !Biot coefficients for the input phi
            D_biot = kappa_s*(1._CUSTOM_REAL + phi*(kappa_s/kappa_f - 1._CUSTOM_REAL))
            H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) &
                      + kappa_fr + 4._CUSTOM_REAL*mu_fr/3._CUSTOM_REAL
            C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
            M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)
            ! Approximated velocities (no viscous dissipation)
            rho_bar =  (1._CUSTOM_REAL - phi)*rho_s + phi*rho_f
            afactor = rho_bar - phi/tort*rho_f
            bfactor = H_biot + phi*rho_bar/(tort*rho_f)*M_biot - TWO*phi/tort*C_biot
            cfactor = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)
            cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
            cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
            cssquare = mu_fr/afactor

            ! AC based on cpI,cpII & cs
            rho_vpI(i,j,k,ispec) = (rho_bar - phi/tort*rho_f)*sqrt(cpIsquare)
            rho_vpII(i,j,k,ispec) = (rho_bar - phi/tort*rho_f)*sqrt(cpIIsquare)
            rho_vsI(i,j,k,ispec) = (rho_bar - phi/tort*rho_f)*sqrt(cssquare)

          endif !if (idomain_id == IDOMAIN_ACOUSTIC .or. idomain_id == IDOMAIN_ELASTIC)

          ! stores anisotropic parameters
          if (ANISOTROPY) then
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


          ! stores material domain
          select case (idomain_id )
          case (IDOMAIN_ACOUSTIC )
            ispec_is_acoustic(ispec) = .true.
          case (IDOMAIN_ELASTIC )
            ispec_is_elastic(ispec) = .true.
          case (IDOMAIN_POROELASTIC )
            ispec_is_poroelastic(ispec) = .true.
          case default
            stop 'Error material domain index'
          end select

        enddo
      enddo
    enddo

    ! user output
    if (myrank == 0) then
      if (mod(ispec,nspec/10) == 0) then
        tCPU = wtime() - time_start
        ! remaining
        tCPU = (10.0-ispec/(nspec/10.0))/ispec/(nspec/10.0)*tCPU
        write(IMAIN,*) "    ",ispec/(nspec/10) * 10," %", &
                      " time remaining:", tCPU,"s"

        ! flushes file buffer for main output file (IMAIN)
        call flush_IMAIN()

      endif
    endif
  enddo

  ! checks material domains
  do ispec=1,nspec
    ! checks if domain is set
    if ((ispec_is_acoustic(ispec) .eqv. .false.) .and. &
        (ispec_is_elastic(ispec) .eqv. .false.) .and. &
        (ispec_is_poroelastic(ispec) .eqv. .false.)) then
      print *,'Error material domain not assigned to element:',ispec
      print *,'acoustic: ',ispec_is_acoustic(ispec)
      print *,'elastic: ',ispec_is_elastic(ispec)
      print *,'poroelastic: ',ispec_is_poroelastic(ispec)
      stop 'Error material domain index element'
    endif
    ! checks if domain is unique
    if (((ispec_is_acoustic(ispec) .eqv. .true.) .and. (ispec_is_elastic(ispec) .eqv. .true.)) .or. &
       ((ispec_is_acoustic(ispec) .eqv. .true.) .and. (ispec_is_poroelastic(ispec) .eqv. .true.)) .or. &
       ((ispec_is_poroelastic(ispec) .eqv. .true.) .and. (ispec_is_elastic(ispec) .eqv. .true.)) .or. &
       ((ispec_is_acoustic(ispec) .eqv. .true.) .and. (ispec_is_elastic(ispec) .eqv. .true.) .and. &
       (ispec_is_poroelastic(ispec) .eqv. .true.))) then
      print *,'Error material domain assigned twice to element:',ispec
      print *,'acoustic: ',ispec_is_acoustic(ispec)
      print *,'elastic: ',ispec_is_elastic(ispec)
      print *,'poroelastic: ',ispec_is_poroelastic(ispec)
      stop 'Error material domain index element'
    endif
  enddo

  ! sets simulation domain flags
  ! all processes will have e.g. acoustic_simulation flag set if any flag is .true. somewhere
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )

  ! deallocates tomographic arrays
#ifndef DEBUG_COUPLED
  if ( (nundefMat_ext_mesh > 0 .or. IMODEL == IMODEL_TOMO) ) call deallocate_tomography_files()
#endif

  end subroutine get_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_model_values(materials_ext_mesh,nmat_ext_mesh, &
                              undef_mat_prop,nundefMat_ext_mesh, &
                              imaterial_id,imaterial_def, &
                              xmesh,ymesh,zmesh, &
                              rho,vp,vs,qkappa_atten,qmu_atten,idomain_id, &
                              rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                              phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                              c11,c12,c13,c14,c15,c16, &
                              c22,c23,c24,c25,c26,c33, &
                              c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              ANISOTROPY)

  use generate_databases_par, only: IMODEL,IMODEL_DEFAULT,IMODEL_GLL,IMODEL_1D_PREM,IMODEL_1D_CASCADIA,IMODEL_1D_SOCAL, &
    IMODEL_SALTON_TROUGH,IMODEL_TOMO,IMODEL_USER_EXTERNAL,IMODEL_IPATI,IMODEL_IPATI_WATER, &
    IMODEL_1D_PREM_PB,IMODEL_GLL, IMODEL_SEP, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,ATTENUATION_COMP_MAXIMUM

  use create_regions_mesh_ext_par

#ifdef DEBUG_COUPLED
  include "../../../add_to_get_model_9.F90"
#endif

  implicit none

  integer, intent(in) :: nmat_ext_mesh
  double precision, dimension(16,nmat_ext_mesh),intent(in) :: materials_ext_mesh

  integer, intent(in) :: nundefMat_ext_mesh
  character(len=MAX_STRING_LEN), dimension(6,nundefMat_ext_mesh) :: undef_mat_prop

  integer, intent(in) :: imaterial_id,imaterial_def

  double precision, intent(in) :: xmesh,ymesh,zmesh

  real(kind=CUSTOM_REAL) :: vp,vs,rho,qkappa_atten,qmu_atten

  integer :: idomain_id

  real(kind=CUSTOM_REAL) :: kappa_s,kappa_f,kappa_fr,mu_fr,rho_s,rho_f,phi,tort,eta_f, &
                           kxx,kxy,kxz,kyy,kyz,kzz

  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                        c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  logical :: ANISOTROPY

  ! local parameters
  integer :: iflag_aniso
  integer :: iundef,imaterial_PB
  logical :: has_tomo_value

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_7.F90"
#endif

  ! use acoustic domains for simulation
  logical,parameter :: USE_PURE_ACOUSTIC_MOD = .false.

  ! initializes with default values
  ! no anisotropy
  iflag_aniso = 0
  idomain_id = IDOMAIN_ELASTIC

  ! attenuation
  ! shear attenuation: arbitrary value, see maximum in constants.h
  qmu_atten = ATTENUATION_COMP_MAXIMUM
  ! bulk attenuation: arbitrary (high) default value
  qkappa_atten = 9999.0_CUSTOM_REAL

  ! selects chosen velocity model
  select case (IMODEL)

  case (IMODEL_DEFAULT,IMODEL_GLL,IMODEL_IPATI,IMODEL_IPATI_WATER, IMODEL_SEP)
    ! material values determined by mesh properties
    call model_default(materials_ext_mesh,nmat_ext_mesh, &
                       undef_mat_prop,nundefMat_ext_mesh, &
                       imaterial_id,imaterial_def, &
                       xmesh,ymesh,zmesh, &
                       rho,vp,vs, &
                       iflag_aniso,qkappa_atten,qmu_atten,idomain_id, &
                       rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                       phi,tort,kxx,kxy,kxz,kyy,kyz,kzz)

#ifdef DEBUG_COUPLED
    include "../../../add_to_get_model_8.F90"
#endif

  case (IMODEL_1D_PREM)
    ! 1D model profile from PREM
    call model_1D_prem_iso(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

  case (IMODEL_1D_PREM_PB)
    ! 1D model profile from PREM modified by Piero
    imaterial_PB = abs(imaterial_id)
    call model_1D_PREM_routine_PB(xmesh,ymesh,zmesh,rho,vp,vs,imaterial_PB)

    ! sets acoustic/elastic domain as given in materials properties
    iundef = - imaterial_id    ! iundef must be positive
    read(undef_mat_prop(6,iundef),*) idomain_id

  case (IMODEL_1D_CASCADIA)
    ! 1D model profile for Cascadia region
    call model_1D_cascadia(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

  case (IMODEL_1D_SOCAL)
    ! 1D model profile for Southern California
    call model_1D_socal(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

  case (IMODEL_SALTON_TROUGH)
    ! gets model values from tomography file
    call model_salton_trough(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

  case (IMODEL_TOMO)
    ! gets model values from tomography file
    call model_tomography(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten,imaterial_id,has_tomo_value)

    ! in case no tomography value defined for this region, fall back to defaults
    if (.not. has_tomo_value) then
      print *,'Error: tomography value not defined for model material id ',imaterial_id
      print *,'Please check if Par_file setting MODEL = tomo is applicable, or try using MODEL = default ...'
      stop 'Error tomo model not found for material'
    endif

  case (IMODEL_USER_EXTERNAL)

    ! Florian Schumacher, Germany, June 2015
    ! FS FS: added call to model_default here, before calling model_external_values in order to
    !        be able to superimpose a model onto the default one:

    ! material values determined by mesh properties
#ifndef DEBUG_COUPLED
    call model_default(materials_ext_mesh,nmat_ext_mesh, &
                       undef_mat_prop,nundefMat_ext_mesh, &
                       imaterial_id,imaterial_def, &
                       xmesh,ymesh,zmesh,rho,vp,vs, &
                       iflag_aniso,qkappa_atten,qmu_atten,idomain_id, &
                       rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                       phi,tort,kxx,kxy,kxz,kyy,kyz,kzz)
#endif

    ! user model from external routine
    ! adds/gets velocity model as specified in model_external_values.f90
    call model_external_values(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id)

  case default
    stop 'Error model not implemented yet'
  end select

  ! adds anisotropic default model
  if (ANISOTROPY) then
    call model_aniso(iflag_aniso,rho,vp,vs, &
                    c11,c12,c13,c14,c15,c16, &
                    c22,c23,c24,c25,c26,c33, &
                    c34,c35,c36,c44,c45,c46,c55,c56,c66)
  endif

  ! for pure acoustic simulations (a way of avoiding re-mesh, re-partition etc.)
  ! can be used to compare elastic & acoustic reflections in exploration seismology
  ! do NOT use it unless you are confident
  if (USE_PURE_ACOUSTIC_MOD) then
    idomain_id = IDOMAIN_ACOUSTIC
  endif

  ! checks if valid vp value
  if (idomain_id == IDOMAIN_ACOUSTIC .or. idomain_id == IDOMAIN_ELASTIC) then
    if (vp <= 0._CUSTOM_REAL) then
      print *,'Error: encountered zero Vp velocity in element! '
      print *,'domain id = ',idomain_id,' material id = ',imaterial_id, 'vp/vs/rho = ',vp,vs,rho
      stop 'Error zero Vp velocity found'
    endif
  endif

  end subroutine get_model_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_model_binaries(myrank,nspec,LOCAL_PATH)

! reads in material parameters from external binary files

  use generate_databases_par, only: IMAIN, IMODEL, IMODEL_GLL,IMODEL_IPATI,IMODEL_IPATI_WATER, IMODEL_SEP, ADIOS_FOR_MESH

  use create_regions_mesh_ext_par

  use model_ipati_adios_mod, only: model_ipati_adios,model_ipati_water_adios

  use model_sep_mod, only: model_sep

  implicit none

  ! number of spectral elements in each block
  integer :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! external GLL models
  ! variables for importing models from files in SPECFEM format, e.g.,  proc000000_vp.bin etc.
  ! can be used for importing updated model in iterative inversions

  ! note: we read in these binary files after mesh coloring, since mesh coloring is permuting arrays.
  !          here, the ordering in **_vp.bin etc. can be permuted as they are outputted when saving mesh files

  select case (IMODEL)
  case (IMODEL_GLL)
    ! note:
    ! import the model from files in SPECFEM format
    ! note that those those files should be saved in LOCAL_PATH
    if (ADIOS_FOR_MESH) then
      call model_gll_adios(myrank,nspec,LOCAL_PATH)
    else
      call model_gll(myrank,nspec,LOCAL_PATH)
    endif

  case (IMODEL_IPATI)
    ! import the model from modified files in SPECFEM format
    if (ADIOS_FOR_MESH) then
      call model_ipati_adios(myrank,nspec,LOCAL_PATH)
    else
      call model_ipati(myrank,nspec,LOCAL_PATH)
    endif

  case (IMODEL_IPATI_WATER)
    ! import the model from modified files in SPECFEM format
    if (ADIOS_FOR_MESH) then
      call model_ipati_water_adios(myrank,nspec,LOCAL_PATH)
    else
      call model_ipati_water(myrank,nspec,LOCAL_PATH)
    endif

  case (IMODEL_SEP)
    ! use values from SEP files
    call model_sep()

  case default
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '     no external binary model used '
    endif
  end select

  end subroutine get_model_binaries
