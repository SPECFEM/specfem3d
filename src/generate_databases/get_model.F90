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

  subroutine get_model()

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,TWO,IMAIN

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION

  use generate_databases_par, only: IMODEL, &
    IMODEL_SALTON_TROUGH,IMODEL_TOMO,IMODEL_USER_EXTERNAL, &
    IMODEL_COUPLED, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
    nspec => NSPEC_AB,ibool,mat_ext_mesh,nundefMat_ext_mesh, &
    ANISOTROPY

  use create_regions_mesh_ext_par

  ! injection technique
  use constants, only: INJECTION_TECHNIQUE_IS_FK,INJECTION_TECHNIQUE_IS_DSM,INJECTION_TECHNIQUE_IS_AXISEM
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH,INJECTION_TECHNIQUE_TYPE

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: vp,vs,rho,qkappa_atten,qmu_atten
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                            c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  real(kind=CUSTOM_REAL) :: kappa_s,kappa_f,kappa_fr,mu_fr,rho_s,rho_f,phi,tort,eta_f, &
                            kxx,kxy,kxz,kyy,kyz,kzz,rho_bar
  real(kind=CUSTOM_REAL) :: cpIsquare,cpIIsquare,cssquare,H_biot,M_biot,C_biot,D_biot, &
                            afactor,bfactor,cfactor
  real(kind=CUSTOM_REAL) :: vpI_poro,vpII_poro,vs_poro

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
  if (nundefMat_ext_mesh > 0 .or. IMODEL == IMODEL_TOMO) call model_tomography_broadcast(myrank)

  ! prepares external model values if needed
  select case (IMODEL)
  case (IMODEL_USER_EXTERNAL)
    call model_external_broadcast()
  case (IMODEL_SALTON_TROUGH)
    call model_salton_trough_broadcast()
  case (IMODEL_COUPLED)
    call model_coupled_broadcast()
  end select

! DP: not sure if this here should check if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) ..
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then

    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK .and. MESH_A_CHUNK_OF_THE_EARTH) &
      stop 'coupling with INJECTION_TECHNIQUE_IS_FK is incompatible with MESH_A_CHUNK_OF_THE_EARTH &
                      &because of Earth curvature not honored'

    !! VM VM for coupling with DSM
    !! find the layer in which the middle of the element is located
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*)
      write(IMAIN,*) '     USING A HYBRID METHOD (THE CODE IS COUPLED WITH AN INJECTION TECHNIQUE)'
      write(IMAIN,*)
      select case(INJECTION_TECHNIQUE_TYPE)
      case (INJECTION_TECHNIQUE_IS_DSM)
        write(IMAIN,*) '     INJECTION TECHNIQUE TYPE = ', INJECTION_TECHNIQUE_TYPE,' (DSM) '
      case (INJECTION_TECHNIQUE_IS_AXISEM)
        write(IMAIN,*) '     INJECTION TECHNIQUE TYPE = ', INJECTION_TECHNIQUE_TYPE,' (AXISEM) '
      case (INJECTION_TECHNIQUE_IS_FK)
        write(IMAIN,*) '     INJECTION TECHNIQUE TYPE = ', INJECTION_TECHNIQUE_TYPE,' (FK) '
      case default
        stop 'Invalid INJECTION_TECHNIQUE_TYPE chosen, must be 1 == DSM, 2 == AXISEM or 3 == FK'
      end select
      write(IMAIN,*)
      write(IMAIN,*)
    endif

    if (INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_FK) then
      ! MESH_A_CHUNK_OF_THE_EARTH, DSM or AxiSEM
      if ((NGLLX == 5) .and. (NGLLY == 5) .and. (NGLLZ == 5)) then
        ! gets xyz coordinates of GLL point
        iglob = ibool(3,3,3,1)
        xmesh = xstore_unique(iglob)
        ymesh = ystore_unique(iglob)
        zmesh = zstore_unique(iglob)
        call model_coupled_FindLayer(xmesh,ymesh,zmesh)
      else
        stop 'bad number of GLL points for coupling with an external code'
      endif
    endif
  endif

  ! get MPI starting time
  time_start = wtime()

  ! material properties on all GLL points: taken from material values defined for
  ! each spectral element in input mesh
  do ispec = 1, nspec

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      if (INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_FK) then
        iglob = ibool(3,3,3,ispec)
        xmesh = xstore_unique(iglob)
        ymesh = ystore_unique(iglob)
        zmesh = zstore_unique(iglob)
        call model_coupled_FindLayer(xmesh,ymesh,zmesh)
      endif
    endif

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
          xmesh = xstore_unique(iglob)
          ymesh = ystore_unique(iglob)
          zmesh = zstore_unique(iglob)

          !! VM VM for coupling with DSM
          !! find the layer in which the middle of the element is located
          if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
            if (INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_FK) then
              if (i == 3 .and. j == 3 .and. k == 3) call model_coupled_FindLayer(xmesh,ymesh,zmesh)
            endif
          endif

          ! material index 1: associated material number
          ! 1 = acoustic, 2 = elastic, 3 = poroelastic, -1 = undefined tomographic
          imaterial_id = mat_ext_mesh(1,ispec)

          ! material index 2: associated material definition
          ! 1 = interface, 2 = tomography material
          imaterial_def = mat_ext_mesh(2,ispec)

          ! assigns material properties
          call get_model_values(imaterial_id,imaterial_def, &
                                xmesh,ymesh,zmesh,ispec, &
                                rho,vp,vs,qkappa_atten,qmu_atten,idomain_id, &
                                rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                                phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                                c11,c12,c13,c14,c15,c16, &
                                c22,c23,c24,c25,c26,c33, &
                                c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                ANISOTROPY)

          ! stores velocity model

          ! converts model coefficients
          ! will only be used to store into array, but not used further
          select case(idomain_id)
          case (IDOMAIN_ACOUSTIC)
            ! parameters for poroelastic model arrays
            ! assuming vs = 0
            rho_s = 0.0     ! solid properties
            kappa_s = 0.0
            rho_f = rho     ! fluid properties
            kappa_f = rho*vp*vp
            eta_f = 0.0
            kappa_fr = rho*vp*vp  ! frame properties
            mu_fr = 0.0
            phi = 1.0       ! fluid
            tort = 0.0
            kxx = 0.0
            kxy = 0.0
            kxz = 0.0
            kyy = 0.0
            kyz = 0.0
            kzz = 0.0
          case (IDOMAIN_ELASTIC)
            ! solid properties
            rho_s = rho     ! solid properties
            kappa_s = rho*(vp*vp - FOUR_THIRDS*vs*vs)
            rho_f = rho     ! fluid properties
            kappa_f = rho*(vp*vp - FOUR_THIRDS*vs*vs)
            eta_f = 0.0
            kappa_fr = rho*(vp*vp - FOUR_THIRDS*vs*vs)  ! frame properties
            mu_fr = rho*vs*vs
            phi = 0.0       ! solid
            tort = 0.0
            kxx = 0.0
            kxy = 0.0
            kxz = 0.0
            kyy = 0.0
            kyz = 0.0
            kzz = 0.0
          case (IDOMAIN_POROELASTIC)
            ! parameters for acoustic/elastic model arrays
            ! selects solid parameters for elastic arrays
            rho = rho_s
            vp = sqrt((kappa_s + FOUR_THIRDS * mu_fr) / rho_s)
            vs = sqrt(mu_fr / rho_s)
          end select

          ! elastic or acoustic material
          ! density
          rhostore(i,j,k,ispec) = rho

          ! kappa, mu
          kappastore(i,j,k,ispec) = rho*( vp*vp - FOUR_THIRDS*vs*vs )
          mustore(i,j,k,ispec) = rho*vs*vs

          ! attenuation
          qkappa_attenuation_store(i,j,k,ispec) = qkappa_atten
          qmu_attenuation_store(i,j,k,ispec) = qmu_atten

          ! Stacey, a completer par la suite
          rho_vp(i,j,k,ispec) = rho*vp
          rho_vs(i,j,k,ispec) = rho*vs

          ! poroelastic material store
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
          if (abs(D_biot - kappa_fr) > 1.e-20) then
            H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS * mu_fr
            C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
            M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)
          else
            H_biot = kappa_s + FOUR_THIRDS * mu_fr
            C_biot = kappa_s
            M_biot = 0.0
          endif
          ! Approximated velocities (no viscous dissipation)
          rho_bar =  (1._CUSTOM_REAL - phi)*rho_s + phi*rho_f
          if (abs(tort) > 1.e-20) then
            afactor = rho_bar - phi/tort*rho_f
            bfactor = H_biot + phi*rho_bar/(tort*rho_f)*M_biot - TWO*phi/tort*C_biot
            cfactor = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)
          else
            afactor = 0.0
            bfactor = 0.0
            cfactor = 0.0
          endif

          if (abs(afactor) > 1.e-20) then
            cpIsquare  = (bfactor + sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
            cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
            cssquare = mu_fr/afactor
          else
            cpIsquare =  vp*vp
            cpIIsquare = vp*vp
            cssquare = mu_fr
          endif

          ! AC based on cpI,cpII & cs
          if (abs(tort) > 1.e-20) then
            rho_vpI(i,j,k,ispec)  = (rho_bar - phi/tort*rho_f) * sqrt(cpIsquare)
            rho_vpII(i,j,k,ispec) = (rho_bar - phi/tort*rho_f) * sqrt(cpIIsquare)
            rho_vsI(i,j,k,ispec)  = (rho_bar - phi/tort*rho_f) * sqrt(cssquare)

            vs_poro   = rho_vsI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f - phi/tort * rho_f )
            vpI_poro  = rho_vpI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f - phi/tort * rho_f )
            vpII_poro = rho_vpII(i,j,k,ispec)/ ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f - phi/tort * rho_f )
          else
            rho_vpI(i,j,k,ispec)  = rho_bar * sqrt(cpIsquare)
            rho_vpII(i,j,k,ispec) = rho_bar * sqrt(cpIIsquare)
            rho_vsI(i,j,k,ispec)  = rho_bar * sqrt(cssquare)

            vs_poro   = rho_vsI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f )
            vpI_poro  = rho_vpI(i,j,k,ispec) / ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f )
            vpII_poro = rho_vpII(i,j,k,ispec)/ ( (1.0_CUSTOM_REAL - phi) * rho_s + phi * rho_f )
          endif

          !debug
          if (.false.) then
            if (myrank == 0 .and. idomain_id == 1 .and. i == 1 .and. j == 1 .and. k == 1) &
              print *,'debug: model domain ',idomain_id,'ispec ',ispec,' vp/vs/rho = ',vp,vs,rho, 'rhovp/rhovs ',rho*vp,rho*vs, &
                      ' rhovpI/rhovpII/rhovsI = ',rho_vpI(i,j,k,ispec),rho_vpII(i,j,k,ispec),rho_vsI(i,j,k,ispec), &
                      ' vp_poro ',vpI_poro,' vpII_poro ',vpII_poro,' vs_poro ',vs_poro, &
                      ' D/H/C/M biot ',D_biot,H_biot,C_biot,M_biot,' rho_bar/a/b/c factor',rho_bar,afactor,bfactor,cfactor, &
                      ' cpIsq/cpIIsq/cssq ',cpIsquare,cpIIsquare,cssquare
          endif

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
          select case (idomain_id)
          case (IDOMAIN_ACOUSTIC)
            ispec_is_acoustic(ispec) = .true.
          case (IDOMAIN_ELASTIC)
            ispec_is_elastic(ispec) = .true.
          case (IDOMAIN_POROELASTIC)
            ispec_is_poroelastic(ispec) = .true.
          case default
            stop 'Error invalid material domain index, must be 1 == acoustic,2 == elastic or 3 == poroelastic'
          end select

        enddo
      enddo
    enddo

    ! user output
    if (myrank == 0) then
      if (mod(ispec,max(nspec/10,1)) == 0) then
        tCPU = wtime() - time_start
        ! remaining
        tCPU = (10.0-ispec/(nspec/10.0))/ispec/(nspec/10.0)*tCPU
        write(IMAIN,*) "    ",ispec/(max(nspec/10,1)) * 10," %", &
                      " time remaining:", tCPU,"s"

        ! flushes file buffer for main output file (IMAIN)
        call flush_IMAIN()

      endif
    endif
  enddo

  ! checks material domains
  do ispec = 1,nspec
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
  if ( (nundefMat_ext_mesh > 0 .or. IMODEL == IMODEL_TOMO) ) call deallocate_tomography_files()

  end subroutine get_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_model_values(imaterial_id,imaterial_def, &
                              xmesh,ymesh,zmesh,ispec, &
                              rho,vp,vs,qkappa_atten,qmu_atten,idomain_id, &
                              rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                              phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                              c11,c12,c13,c14,c15,c16, &
                              c22,c23,c24,c25,c26,c33, &
                              c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              ANISOTROPY)

  use generate_databases_par, only: IMODEL,IMODEL_DEFAULT,IMODEL_GLL,IMODEL_1D_PREM,IMODEL_1D_CASCADIA,IMODEL_1D_SOCAL, &
    IMODEL_SALTON_TROUGH,IMODEL_TOMO,IMODEL_USER_EXTERNAL,IMODEL_IPATI,IMODEL_IPATI_WATER, &
    IMODEL_1D_PREM_PB,IMODEL_GLL, IMODEL_SEP,IMODEL_COUPLED, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,ATTENUATION_COMP_MAXIMUM

  use generate_databases_par, only: undef_mat_prop

  use create_regions_mesh_ext_par

  use constants, only: INJECTION_TECHNIQUE_IS_FK
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE

  implicit none

  integer, intent(in) :: imaterial_id,imaterial_def

  double precision, intent(in) :: xmesh,ymesh,zmesh
  integer, intent(in) :: ispec

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

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation, local parameters
   real(kind=CUSTOM_REAL) :: tem1,tem2,ratio_x,ratio

!  flag indicating whether we smooth the model on the edges to make it 1D on the edges to match the FK calculation
   logical, parameter :: SMOOTH_THE_MODEL_EDGES_FOR_FK = .false.

! if the flag SMOOTH_THE_MODEL_EDGES_FOR_FK above is true, do we smooth the Moho only, or the rest of the model as well
! if the flag SMOOTH_THE_MODEL_EDGES_FOR_FK above is false, then this flag is ignored
   logical, parameter :: SMOOTH_THE_MOHO_ONLY_FOR_FK = .false.
! *********************************************************************************

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

  case (IMODEL_DEFAULT, IMODEL_GLL, IMODEL_IPATI, IMODEL_IPATI_WATER, IMODEL_SEP)
    ! material values determined by mesh properties
    call model_default(imaterial_id,imaterial_def, &
                       xmesh,ymesh,zmesh, &
                       rho,vp,vs, &
                       iflag_aniso,qkappa_atten,qmu_atten,idomain_id, &
                       rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                       phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                       c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  case (IMODEL_1D_PREM)
    ! 1D model profile from PREM
    call model_1D_prem_iso(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten,qkappa_atten)

  case (IMODEL_1D_PREM_PB)
    ! 1D model profile from PREM modified by Piero
    imaterial_PB = abs(imaterial_id)
    call model_1D_PREM_routine_PB(xmesh,ymesh,zmesh,rho,vp,vs,imaterial_PB)

    ! sets acoustic/elastic domain as given in materials properties
    iundef = - imaterial_id    ! iundef must be positive
    read(undef_mat_prop(6,iundef),*) idomain_id

  case (IMODEL_1D_CASCADIA)
    ! 1D model profile for Cascadia region
    call model_1D_cascadia(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten,qkappa_atten)

  case (IMODEL_1D_SOCAL)
    ! 1D model profile for Southern California
    call model_1D_socal(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten,qkappa_atten)

  case (IMODEL_SALTON_TROUGH)
    ! gets model values from tomography file
    call model_salton_trough(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten,qkappa_atten)

  case (IMODEL_TOMO)
    ! gets model values from tomography file
    call model_tomography(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten, &
                          c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                          imaterial_id,has_tomo_value)

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
    call model_default(imaterial_id,imaterial_def, &
                       xmesh,ymesh,zmesh,rho,vp,vs, &
                       iflag_aniso,qkappa_atten,qmu_atten,idomain_id, &
                       rho_s,kappa_s,rho_f,kappa_f,eta_f,kappa_fr,mu_fr, &
                       phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                       c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

    ! user model from external routine
    ! adds/gets velocity model as specified in model_external_values.f90
    call model_external_values(xmesh,ymesh,zmesh,ispec,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id)

  case (IMODEL_COUPLED)
    ! user model for coupling with injection method
    ! adds/gets velocity model as specified in model_coupled.f90
    call model_coupled_values(xmesh,ymesh,zmesh,rho,vp,vs)

  case default
    stop 'Error model not implemented yet'
  end select

  ! adds anisotropic default model
  if (IMODEL /= IMODEL_TOMO .and. ANISOTROPY) then
    call model_aniso(iflag_aniso,rho,vp,vs, &
                     c11,c12,c13,c14,c15,c16, &
                     c22,c23,c24,c25,c26,c33, &
                     c34,c35,c36,c44,c45,c46,c55,c56,c66)
  endif

  ! *********************************************************************************
  ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
  ! for smoothing the boundary portions of the model
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK .and. SMOOTH_THE_MODEL_EDGES_FOR_FK) then
      ! only applies to a default/gll model with specific layering
      ! includes user specifics, we might consider defining this as a separate model
      if (IMODEL == IMODEL_DEFAULT .or. IMODEL == IMODEL_GLL) then
        !! VM VM this seems to be an hardcoded model we should remove it
        if (.not. SMOOTH_THE_MOHO_ONLY_FOR_FK) then
          ratio = 1.0
          ! fast anomaly
          if (imaterial_id == 4) then
            tem1 = (-8.0/20.0)*xmesh + 100000.0
            tem2 = (-8.0/20.0)*xmesh + 140000.0
            if (zmesh > tem1 .and. zmesh < tem2) then
              if (xmesh < 30000.0) then
                ratio_x = 1.0 + 0.04*xmesh/30000.0
              else if (xmesh > 170000.0) then
                ratio_x = 1.0 + 0.04*(200000.0-xmesh)/30000.0
              else
                ratio_x = 1.04
              endif
              if (ymesh < 30000.0) then
                ratio = 1.0 + 0.04*ymesh/30000.0
              else if (ymesh > 170000.0) then
                ratio = 1.0 + 0.04*(200000.0-ymesh)/30000.0
              else
                ratio = 1.04
              endif
              if (ratio_x < ratio) then
                ratio = ratio_x
              endif
              rho = rho*ratio
              vp  = vp *ratio
              vs  = vs *ratio
            endif
          endif

          ! slow anomaly
          if (imaterial_id == 3) then
            tem1 = (-8.0/20.0)*xmesh + 140000.0
            tem2 = (-8.0/20.0)*xmesh + 155000.0
            if (zmesh > tem1 .and. zmesh < tem2) then
              if (xmesh < 30000.0) then
                ratio_x = 1.0 - 0.06*xmesh/30000.0
              else if (xmesh > 170000.0) then
                ratio_x = 1.0 - 0.06*(200000.0-xmesh)/30000.0
              else
                ratio_x = 0.94
              endif
              if (ymesh < 30000.0) then
                ratio = 1.0 - 0.06*ymesh/30000.0
              else if (ymesh > 170000.0) then
                ratio = 1.0 - 0.06*(200000.0-ymesh)/30000.0
              else
                ratio = 0.94
              endif
              if (ratio_x > ratio) then
                ratio = ratio_x
              endif
              rho = rho*ratio
              vp  = vp *ratio
              vs  = vs *ratio
            endif
          endif

        endif ! of if (.not. SMOOTH_THE_MOHO_ONLY_FOR_FK) then

        ! Moho section
        if (imaterial_id == 1) then
          tem1 = 160000.0
          tem2 = 170000.0
          if (zmesh > tem1 .and. zmesh < tem2) then
            if (ymesh < 30000.0) then
              vp  = ((30000.0-ymesh)/30000.0)*5800.0 + (ymesh/30000.0)*8080.0
              vs  = ((30000.0-ymesh)/30000.0)*3198.0 + (ymesh/30000.0)*4485.0
              rho = ((30000.0-ymesh)/30000.0)*2600.0 + (ymesh/30000.0)*3380.0
            else if (ymesh > 170000.0) then
              vp  = ((200000.0-ymesh)/30000.0)*5800.0+((ymesh-170000.0)/30000.0)*8080.0
              vs  = ((200000.0-ymesh)/30000.0)*3198.0+((ymesh-170000.0)/30000.0)*4485.0
              rho = ((200000.0-ymesh)/30000.0)*2600.0+((ymesh-170000.0)/30000.0)*3380.0
            endif
          endif
        endif
      endif
    endif
  endif
  ! *********************************************************************************

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

  subroutine get_model_binaries(nspec,LOCAL_PATH)

! reads in material parameters from external binary files

  use constants, only: myrank,IMAIN,IMODEL_GLL,IMODEL_IPATI,IMODEL_IPATI_WATER,IMODEL_SEP

  use generate_databases_par, only: IMODEL,ADIOS_FOR_MESH

  use create_regions_mesh_ext_par

  use model_sep_mod, only: model_sep

  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec
  character(len=MAX_STRING_LEN),intent(in) :: LOCAL_PATH

  ! external GLL models
  ! variables for importing models from files in SPECFEM format, e.g.,  proc000000_vp.bin etc.
  ! can be used for importing updated model in iterative inversions

  ! note: we read in these binary files after mesh coloring, since mesh coloring is permuting arrays.
  !          here, the ordering in **_vp.bin etc. can be permuted as they are outputted when saving mesh files

  select case (IMODEL)
  case (IMODEL_GLL)
    ! note:
    ! import the model from files in SPECFEM format
    ! note that those files should be saved in LOCAL_PATH
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
