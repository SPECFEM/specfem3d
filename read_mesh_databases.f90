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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine read_mesh_databases()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none
  
  integer :: i,j,k,ispec,iglob
  integer :: iinterface,ier
  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  
! start reading the databasesa

! info about external mesh simulation
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',&
      action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error: could not open database '
    print*,'path: ',prname(1:len_trim(prname))//'external_mesh.bin'
    call exit_mpi(myrank,'error opening database')
  endif
  
  read(27) NSPEC_AB
  read(27) NGLOB_AB

  read(27) ibool
  
  read(27) xstore
  read(27) ystore
  read(27) zstore
  
  read(27) xix
  read(27) xiy
  read(27) xiz
  read(27) etax
  read(27) etay
  read(27) etaz
  read(27) gammax
  read(27) gammay
  read(27) gammaz
  read(27) jacobian

  read(27) kappastore
  read(27) mustore

  read(27) ispec_is_acoustic
  read(27) ispec_is_elastic
  read(27) ispec_is_poroelastic

  ! acoustic
  ! all processes will have acoustic_simulation set if any flag is .true.  
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  if( ACOUSTIC_SIMULATION ) then    
    ! potentials
    allocate(potential_acoustic(NGLOB_AB))
    allocate(potential_dot_acoustic(NGLOB_AB))
    allocate(potential_dot_dot_acoustic(NGLOB_AB))
    
    ! mass matrix, density
    allocate(rmass_acoustic(NGLOB_AB))
    allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    
    read(27) rmass_acoustic    
    read(27) rhostore            
  endif

  ! elastic
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  if( ELASTIC_SIMULATION ) then
    ! displacement,velocity,acceleration  
    allocate(displ(NDIM,NGLOB_AB))
    allocate(veloc(NDIM,NGLOB_AB))
    allocate(accel(NDIM,NGLOB_AB))

    allocate(rmass(NGLOB_AB))
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(iflag_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))
    allocate(c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO))

    read(27) rmass
    if( OCEANS ) then
      read(27) rmass_ocean_load
    endif
    !pll
    read(27) rho_vp
    read(27) rho_vs
    read(27) iflag_attenuation_store
    
  else    
    ! no elastic attenuation & anisotropy
    ATTENUATION = .false.
    ANISOTROPY = .false.
  endif
  
  ! poroelastic
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )  
  if( POROELASTIC_SIMULATION ) then
  
    stop 'not implemented yet: read rmass_solid_poroelastic .. '
    
    allocate(rmass_solid_poroelastic(NGLOB_AB))
    allocate(rmass_fluid_poroelastic(NGLOB_AB))

    read(27) rmass_solid_poroelastic
    read(27) rmass_fluid_poroelastic    
  endif

! checks simulation types are valid
  if( (.not. ACOUSTIC_SIMULATION ) .and. &
     (.not. ELASTIC_SIMULATION ) .and. &
     (.not. POROELASTIC_SIMULATION ) ) then
     close(27)
     call exit_mpi(myrank,'error no simulation type defined')
  endif
  
  ! checks attenuation flags: see integers defined in constants.h
  if( ATTENUATION ) then
    if( minval(iflag_attenuation_store(:,:,:,:)) < 1 ) then
      close(27)
      call exit_MPI(myrank,'error attenuation flag entry exceeds range')
    endif
    if( maxval(iflag_attenuation_store(:,:,:,:)) > NUM_REGIONS_ATTENUATION ) then
      close(27)
      call exit_MPI(myrank,'error attenuation flag entry exceeds range')
    endif
  endif        
  
! absorbing boundary surface
  read(27) num_abs_boundary_faces
  allocate(abs_boundary_ispec(num_abs_boundary_faces))
  allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces))
  allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces))
  allocate(abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces))
  read(27) abs_boundary_ispec
  read(27) abs_boundary_ijk
  read(27) abs_boundary_jacobian2Dw
  read(27) abs_boundary_normal

! free surface 
  read(27) num_free_surface_faces
  allocate(free_surface_ispec(num_free_surface_faces))
  allocate(free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces))
  allocate(free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces))
  allocate(free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces))
  read(27) free_surface_ispec
  read(27) free_surface_ijk
  read(27) free_surface_jacobian2Dw
  read(27) free_surface_normal

! acoustic-elastic coupling surface
  read(27) num_coupling_ac_el_faces
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_ispec(num_coupling_ac_el_faces))
  read(27) coupling_ac_el_ispec   
  read(27) coupling_ac_el_ijk
  read(27) coupling_ac_el_jacobian2Dw 
  read(27) coupling_ac_el_normal 
    
! MPI interfaces
  read(27) num_interfaces_ext_mesh
  read(27) max_nibool_interfaces_ext_mesh
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh))
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh))
  allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  read(27) my_neighbours_ext_mesh
  read(27) nibool_interfaces_ext_mesh
  read(27) ibool_interfaces_ext_mesh

  if( ANISOTROPY ) then
    read(27) c11store
    read(27) c12store
    read(27) c13store
    read(27) c14store
    read(27) c15store
    read(27) c16store
    read(27) c22store
    read(27) c23store
    read(27) c24store
    read(27) c25store
    read(27) c26store
    read(27) c33store
    read(27) c34store
    read(27) c35store
    read(27) c36store
    read(27) c44store
    read(27) c45store
    read(27) c46store
    read(27) c55store
    read(27) c56store
    read(27) c66store  
  endif
  
  close(27)

! MPI communications
  allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(request_send_vector_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_recv_vector_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_send_scalar_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_recv_scalar_ext_mesh(num_interfaces_ext_mesh))

! locate inner and outer elements
  allocate(ispec_is_inner(NSPEC_AB))
  allocate(iglob_is_inner(NGLOB_AB))
  ispec_is_inner(:) = .true.
  iglob_is_inner(:) = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    do i = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(i,iinterface)
      iglob_is_inner(iglob) = .false.
    enddo
  enddo
  do ispec = 1, NSPEC_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          ispec_is_inner(ispec) = iglob_is_inner(iglob) .and. ispec_is_inner(ispec)
        enddo
      enddo
    enddo
  enddo
  deallocate( iglob_is_inner )  

! sets up elements for loops in acoustic simulations
  if( ACOUSTIC_SIMULATION ) then
    ! counts inner and outer elements
    nspec_inner_acoustic = 0
    nspec_outer_acoustic = 0
    do ispec = 1, NSPEC_AB
      if( ispec_is_acoustic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          nspec_inner_acoustic = nspec_inner_acoustic + 1
        else
          nspec_outer_acoustic = nspec_outer_acoustic + 1
        endif
      endif
    enddo
        
    ! stores indices of inner and outer elements for faster(?) computation 
    num_phase_ispec_acoustic = max(nspec_inner_acoustic,nspec_outer_acoustic)
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2))
    nspec_inner_acoustic = 0
    nspec_outer_acoustic = 0
    do ispec = 1, NSPEC_AB
      if( ispec_is_acoustic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          nspec_inner_acoustic = nspec_inner_acoustic + 1
          phase_ispec_inner_acoustic(nspec_inner_acoustic,2) = ispec
        else
          nspec_outer_acoustic = nspec_outer_acoustic + 1
          phase_ispec_inner_acoustic(nspec_outer_acoustic,1) = ispec
        endif
      endif
    enddo
    !print *,'rank ',myrank,' acoustic inner spec: ',nspec_inner_acoustic
    !print *,'rank ',myrank,' acoustic outer spec: ',nspec_outer_acoustic
  endif

! sets up elements for loops in acoustic simulations
  if( ELASTIC_SIMULATION ) then
    ! counts inner and outer elements
    nspec_inner_elastic = 0
    nspec_outer_elastic = 0
    do ispec = 1, NSPEC_AB
      if( ispec_is_elastic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          nspec_inner_elastic = nspec_inner_elastic + 1
        else
          nspec_outer_elastic = nspec_outer_elastic + 1
        endif
      endif
    enddo
        
    ! stores indices of inner and outer elements for faster(?) computation 
    num_phase_ispec_elastic = max(nspec_inner_elastic,nspec_outer_elastic)
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2))
    nspec_inner_elastic = 0
    nspec_outer_elastic = 0
    do ispec = 1, NSPEC_AB
      if( ispec_is_elastic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          nspec_inner_elastic = nspec_inner_elastic + 1
          phase_ispec_inner_elastic(nspec_inner_elastic,2) = ispec
        else
          nspec_outer_elastic = nspec_outer_elastic + 1
          phase_ispec_inner_elastic(nspec_outer_elastic,1) = ispec
        endif
      endif
    enddo
    !print *,'rank ',myrank,' elastic inner spec: ',nspec_inner_elastic
    !print *,'rank ',myrank,' elastic outer spec: ',nspec_outer_elastic
  endif



! gets model dimensions  
  minl = minval( xstore )
  maxl = maxval( xstore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LONGITUDE_MIN = min_all
  LONGITUDE_MAX = max_all

  minl = minval( ystore )
  maxl = maxval( ystore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LATITUDE_MIN = min_all
  LATITUDE_MAX = max_all
  
! check courant criteria on mesh
  if( ELASTIC_SIMULATION ) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                        kappastore,mustore,rho_vp,rho_vs,DT,model_speed_max )
  else if( ACOUSTIC_SIMULATION ) then  
      allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      rho_vp = sqrt( kappastore / rhostore ) * rhostore
      rho_vs = 0.0_CUSTOM_REAL
      call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                        kappastore,mustore,rho_vp,rho_vs,DT,model_speed_max )
      deallocate(rho_vp,rho_vs)
  endif

! reads adjoint parameters
  call read_mesh_databases_adjoint()

  end subroutine read_mesh_databases
  
!
!-------------------------------------------------------------------------------------------------
!  

  subroutine read_mesh_databases_adjoint()

! reads in moho meshes

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none
  
  integer :: ier

! allocates adjoint arrays for elastic simulations
  if( ELASTIC_SIMULATION .and. SIMULATION_TYPE == 3 ) then
    ! backward displacement,velocity,acceleration fields  
    allocate(b_displ(NDIM,NGLOB_ADJOINT))
    allocate(b_veloc(NDIM,NGLOB_ADJOINT))
    allocate(b_accel(NDIM,NGLOB_ADJOINT))
  
    ! adjoint kernels

    ! primary, isotropic kernels
    ! density kernel
    allocate(rho_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT))
    ! shear modulus kernel
    allocate(mu_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT))
    ! compressional modulus kernel
    allocate(kappa_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT))

    ! derived kernels
    ! density prime kernel
    allocate(rhop_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT))
    ! vp kernel
    allocate(alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT))
    ! vs kernel
    allocate(beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT))
    
    ! MPI handling
    allocate(b_request_send_vector_ext_mesh(num_interfaces_ext_mesh))
    allocate(b_request_recv_vector_ext_mesh(num_interfaces_ext_mesh))    
    allocate(b_buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
    allocate(b_buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
    
  endif

! allocates adjoint arrays for acoustic simulations
  if( ACOUSTIC_SIMULATION .and. SIMULATION_TYPE == 3 ) then
    ! backward potentials  
    allocate(b_potential_acoustic(NGLOB_ADJOINT))
    allocate(b_potential_dot_acoustic(NGLOB_ADJOINT))
    allocate(b_potential_dot_dot_acoustic(NGLOB_ADJOINT))
    
    ! kernels
    allocate(rho_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT)) 
    allocate(rhop_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT)) 
    allocate(kappa_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT)) 
    allocate(alpha_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT)) 

    ! MPI handling
    allocate(b_request_send_scalar_ext_mesh(num_interfaces_ext_mesh))
    allocate(b_request_recv_scalar_ext_mesh(num_interfaces_ext_mesh))    
    allocate(b_buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
    allocate(b_buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
    
  endif

! allocates attenuation solids
  if( ATTENUATION .and. SIMULATION_TYPE == 3 ) then
    allocate(b_R_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS), &
            b_R_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS), &
            b_R_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS), &
            b_R_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS), &
            b_R_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS) )
            
    allocate(b_epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL), &
            b_epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL), &
            b_epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL), &
            b_epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL), &
            b_epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL) )    
  endif
  
! ADJOINT moho
! moho boundary
  if( ELASTIC_SIMULATION ) then
    allocate( is_moho_top(NSPEC_BOUN),is_moho_bot(NSPEC_BOUN) )

    if( SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3 ) then
    
      ! boundary elements
      !open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='unknown',form='unformatted')
      open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='old',&
            form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: could not open ibelm_moho '
        print*,'path: ',prname(1:len_trim(prname))//'ibelm_moho.bin'
        call exit_mpi(myrank,'error opening ibelm_moho')
      endif
      
      read(27) NSPEC2D_MOHO
      
      ! allocates arrays for moho mesh
      allocate(ibelm_moho_bot(NSPEC2D_MOHO))
      allocate(ibelm_moho_top(NSPEC2D_MOHO))
      allocate(normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO))
      allocate(normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO))
      allocate(ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO))
      allocate(ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO))

      read(27) ibelm_moho_top
      read(27) ibelm_moho_bot
      read(27) ijk_moho_top
      read(27) ijk_moho_bot
      
      close(27)

      ! normals
      open(unit=27,file=prname(1:len_trim(prname))//'normal_moho.bin',status='old',&
            form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: could not open normal_moho '
        print*,'path: ',prname(1:len_trim(prname))//'normal_moho.bin'
        call exit_mpi(myrank,'error opening normal_moho')
      endif
      
      read(27) normal_moho_top
      read(27) normal_moho_bot    
      close(27)

      ! flags    
      open(unit=27,file=prname(1:len_trim(prname))//'is_moho.bin',status='old',&
            form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: could not open is_moho '
        print*,'path: ',prname(1:len_trim(prname))//'is_moho.bin'
        call exit_mpi(myrank,'error opening is_moho')
      endif
      
      read(27) is_moho_top
      read(27) is_moho_bot    
      
      close(27)
      
      ! moho kernel
      allocate( moho_kl(NGLLSQUARE,NSPEC2D_MOHO) )      
      moho_kl = 0._CUSTOM_REAL
      
    else
      NSPEC2D_MOHO = 1
    endif
  
    allocate( dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
                                   dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
                                   b_dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
                                   b_dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO) )  
  endif
  
  end subroutine read_mesh_databases_adjoint
