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

  subroutine get_model(nspec,ibool,mat_ext_mesh,nelmnts_ext_mesh, &
                        materials_ext_mesh,nmat_ext_mesh, &
                        undef_mat_prop,nundefMat_ext_mesh, &
                        ANISOTROPY)

  use create_regions_mesh_ext_par 
  implicit none

  ! number of spectral elements in each block
  integer :: nspec

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
  
! !  Piero, read bedrock file
!  allocate(ibedrock(NX_TOPO_ANT,NY_TOPO_ANT))              
!  if(myrank == 0) then
!      call read_bedrock_file(ibedrock)
!  !    write(IMAIN,*)
!  !    write(IMAIN,*) 'regional bedrock file read ranges in m from ',minval(ibedrock),' to ',maxval(ibedrock)
!  !    write(IMAIN,*)
!   endif
!  ! broadcast the information read on the master to the nodes
!  ! call MPI_BCAST(ibedrock,NX_TOPO_ANT*NY_TOPO_ANT,MPI_REAL,0,MPI_COMM_WORLD,ier)
! call bcast_all_cr(ibedrock,NX_TOPO_ANT*NY_TOPO_ANT)

  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

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
           else

              stop 'material: tomography not implemented yet'
              ! call tomography()

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
!    utm_x_station(1) =  783500.6250000d0
!    utm_y_station(1) = -11828.7519531d0

!    utm_x_station(2) =  853644.5000000d0
!    utm_y_station(2) = -114.0138092d0

!    utm_x_station(3) = 863406.0000000d0
!    utm_y_station(3) = -53736.1640625d0

!    utm_x_station(4) =   823398.8125000d0
!    utm_y_station(4) = 29847.4511719d0

!    utm_x_station(5) = 863545.3750000d0
!    utm_y_station(5) = 19669.6621094d0

!    utm_x_station(6) =  817099.3750000d0
!    utm_y_station(6) = -24430.2871094d0

!  print*,myrank,'après store the position of the six stations'
!  call flush(6)

!  print*, myrank,minval(nodes_coords_ext_mesh(1,:))
!  call flush(6)


! print*, myrank,maxval(nodes_coords_ext_mesh(1,:))
!  call flush(6)


!  do ispec = 1, nspec

!     zmesh = zstore(2,2,2,ispec)

!    ! if(doubling_index == IFLAG_ONE_LAYER_TOPOGRAPHY) then
!     if(any(ibelm_top == ispec)) then
!     doubling_value_found_for_Piero = IFLAG_ONE_LAYER_TOPOGRAPHY
       
!     else if(zmesh < Z_23p4km) then
!        doubling_value_found_for_Piero = IFLAG_MANTLE_BELOW_23p4km
       
!     else if(zmesh < Z_14km) then
!        doubling_value_found_for_Piero = IFLAG_14km_to_23p4km
       
!     else
!        doubling_value_found_for_Piero = IFLAG_BEDROCK_down_to_14km
!     endif
!    idoubling(ispec) = doubling_value_found_for_Piero

!     do k = 1, NGLLZ
!       do j = 1, NGLLY
!         do i = 1, NGLLX

           
!            if(idoubling(ispec) == IFLAG_ONE_LAYER_TOPOGRAPHY .or. &
!               idoubling(ispec) == IFLAG_BEDROCK_down_to_14km) then
              
!               ! since we have suppressed UTM projection for Piero Basini, UTMx is the same as long
!               ! and UTMy is the same as lat
!               long = xstore(i,j,k,ispec)
!               lat = ystore(i,j,k,ispec)
              
!               ! get coordinate of corner in model
!               icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
!               icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1
              
!               ! avoid edge effects and extend with identical point if outside model
!               if(icornerlong < 1) icornerlong = 1
!               if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
!               if(icornerlat < 1) icornerlat = 1
!               if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1
              
!               ! compute coordinates of corner
!               long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
!               lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO
                   
!               ! compute ratio for interpolation
!               ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
!               ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO
                   
!               ! avoid edge effects
!               if(ratio_xi < 0.) ratio_xi = 0.
!               if(ratio_xi > 1.) ratio_xi = 1.
!               if(ratio_eta < 0.) ratio_eta = 0.
!               if(ratio_eta > 1.) ratio_eta = 1.
                   
!               ! interpolate elevation at current point
!               elevation_bedrock = &
!                    ibedrock(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!                    ibedrock(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!                    ibedrock(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!                    ibedrock(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta
                   
!               !! DK DK exclude circles around each station to make sure they are on the bedrock
!               !! DK DK and not in the ice
!               is_around_a_station = .false.
!               do istation = 1,NUMBER_OF_STATIONS
!                  if(sqrt((long - utm_x_station(istation))**2 + (lat - utm_y_station(istation))**2) < RADIUS_TO_EXCLUDE) then
!                     is_around_a_station = .true.
!                     exit
!                  endif
!               enddo
              
!               ! define elastic parameters in the model
              
!               ! we are above the bedrock interface i.e. in the ice, and not too close to a station
!               if(zmesh >= elevation_bedrock .and. .not. is_around_a_station) then
!                  vp = 3800.d0
!                  vs = 1900.d0
!                  rho = 900.d0
!                  iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_ICE
                 
!                  ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!               else
!                  vp = 5800.d0
!                  vs = 3200.d0
!                  rho = 2600.d0
!                  iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
!               endif
              
!            else if(idoubling(ispec) == IFLAG_14km_to_23p4km) then
!               vp = 6800.d0
!               vs = 3900.d0
!               rho = 2900.d0
!               iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
              
!            else if(idoubling(ispec) == IFLAG_MANTLE_BELOW_23p4km) then
!               vp = 8100.d0
!               vs = 4480.d0
!               rho = 3380.d0
!               iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
              
!            endif
           
!                 !pll  8/06
!                     if(CUSTOM_REAL == SIZE_REAL) then
!                        rhostore(i,j,k,ispec) = sngl(rho)
!                        vpstore(i,j,k,ispec) = sngl(vp)
!                        vsstore(i,j,k,ispec) = sngl(vs)
!                     else
!                        rhostore(i,j,k,ispec) = rho
!                        vpstore(i,j,k,ispec) = vp
!                        vsstore(i,j,k,ispec) = vs
!                     end if
                
!                 kappastore(i,j,k,ispec) = rhostore(i,j,k,ispec)*(vpstore(i,j,k,ispec)*vpstore(i,j,k,ispec) - &
!                      4.d0*vsstore(i,j,k,ispec)*vsstore(i,j,k,ispec)/3.d0)
!                 mustore(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)*&
!                      vsstore(i,j,k,ispec)
           
!                 ! Stacey, a completer par la suite  
!                 rho_vp(i,j,k,ispec) = rhostore(i,j,k,ispec)*vpstore(i,j,k,ispec)
!                 rho_vs(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)
!                 !end pll
                
!                 !      kappastore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
!                 !       (materials_ext_mesh(2,mat_ext_mesh(ispec))*materials_ext_mesh(2,mat_ext_mesh(ispec)) - &
!                 !        4.d0*materials_ext_mesh(3,mat_ext_mesh(ispec))*materials_ext_mesh(3,mat_ext_mesh(ispec))/3.d0)
!                 !      mustore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
!                                                         materials_ext_mesh(3,mat_ext_mesh(ispec))*&
!                 !  x    materials_ext_mesh(3,mat_ext_mesh(ispec))
!              enddo
!           enddo
!        enddo
!     enddo

  end subroutine get_model

