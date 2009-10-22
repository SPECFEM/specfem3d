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
  
  implicit none

! start reading the databasesa

! info about external mesh simulation
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',action='read',form='unformatted')
  read(27) NSPEC_AB
  read(27) NGLOB_AB
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
  
  !pll
  read(27) rho_vp
  read(27) rho_vs
  read(27) iflag_attenuation_store

  ! checks attenuation flags: see integers defined in constants.h
  if( ATTENUATION ) then
    if( minval(iflag_attenuation_store(:,:,:,:)) < 1 ) then
      close(27)
      call exit_MPI(myrank,'something is wrong with the mesh attenuation: flag entry is invalid')
    endif
    if( maxval(iflag_attenuation_store(:,:,:,:)) > NUM_REGIONS_ATTENUATION ) then
      close(27)
      call exit_MPI(myrank,'something is wrong with the mesh attenuation: flag entry exceeds number of defined attenuation flags in constants.h')
    endif
  endif        
  
!  read(27) NSPEC2DMAX_XMIN_XMAX_ext 
!  read(27) NSPEC2DMAX_YMIN_YMAX_ext
!  allocate(nimin(2,NSPEC2DMAX_YMIN_YMAX_ext),nimax(2,NSPEC2DMAX_YMIN_YMAX_ext),nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX_ext))
!  allocate(njmin(2,NSPEC2DMAX_XMIN_XMAX_ext),njmax(2,NSPEC2DMAX_XMIN_XMAX_ext),nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX_ext))
!  read(27) nimin
!  read(27) nimax
!  read(27) njmin
!  read(27) njmax
!  read(27) nkmin_xi 
!  read(27) nkmin_eta
  !end pll

  read(27) kappastore
  read(27) mustore
  read(27) rmass
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore

! absorbing boundaries
  !pll
!  read(27) nspec2D_xmin
!  read(27) nspec2D_xmax
!  read(27) nspec2D_ymin
!  read(27) nspec2D_ymax
!  read(27) NSPEC2D_BOTTOM
!  read(27) NSPEC2D_TOP    
!  allocate(ibelm_xmin(nspec2D_xmin))
!  allocate(ibelm_xmax(nspec2D_xmax))
!  allocate(ibelm_ymin(nspec2D_ymin))
!  allocate(ibelm_ymax(nspec2D_ymax))
!  allocate(ibelm_bottom(NSPEC2D_BOTTOM))  
!  allocate(ibelm_top(NSPEC2D_TOP))
!
!  allocate(ibelm_gll_xmin(3,NGLLY,NGLLZ,nspec2D_xmin))
!  allocate(ibelm_gll_xmax(3,NGLLY,NGLLZ,nspec2D_xmax))
!  allocate(ibelm_gll_ymin(3,NGLLX,NGLLZ,nspec2D_ymin))
!  allocate(ibelm_gll_ymax(3,NGLLX,NGLLZ,nspec2D_ymax))
!  allocate(ibelm_gll_bottom(3,NGLLY,NGLLY,nspec2D_bottom))
!  allocate(ibelm_gll_top(3,NGLLY,NGLLY,nspec2D_top))
!              
!  allocate(jacobian2D_xmin(NGLLY,NGLLZ,nspec2D_xmin))
!  allocate(jacobian2D_xmax(NGLLY,NGLLZ,nspec2D_xmax))
!  allocate(jacobian2D_ymin(NGLLX,NGLLZ,nspec2D_ymin))
!  allocate(jacobian2D_ymax(NGLLX,NGLLZ,nspec2D_ymax))
!  allocate(jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM))
!  allocate(jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP))
!  
!  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,nspec2D_xmin))
!  allocate(normal_xmax(NDIM,NGLLY,NGLLZ,nspec2D_xmax))
!  allocate(normal_ymin(NDIM,NGLLX,NGLLZ,nspec2D_ymin))
!  allocate(normal_ymax(NDIM,NGLLX,NGLLZ,nspec2D_ymax))
!  allocate(normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM))
!  allocate(normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP))
!  read(27) ibelm_xmin
!  read(27) ibelm_xmax
!  read(27) ibelm_ymin
!  read(27) ibelm_ymax
!  read(27) ibelm_bottom
!  read(27) ibelm_top
!  
!  read(27) ibelm_gll_xmin
!  read(27) ibelm_gll_xmax
!  read(27) ibelm_gll_ymin
!  read(27) ibelm_gll_ymax
!  read(27) ibelm_gll_bottom
!  read(27) ibelm_gll_top
!  
!  read(27) normal_xmin
!  read(27) normal_xmax
!  read(27) normal_ymin
!  read(27) normal_ymax
!  read(27) normal_bottom
!  read(27) normal_top
!  read(27) jacobian2D_xmin
!  read(27) jacobian2D_xmax
!  read(27) jacobian2D_ymin
!  read(27) jacobian2D_ymax
!  read(27) jacobian2D_bottom
!  read(27) jacobian2D_top
!  !end pll

  read(27) num_absorbing_boundary_faces
  allocate(absorbing_boundary_ispec(num_absorbing_boundary_faces))
  allocate(absorbing_boundary_ijk(3,NGLLSQUARE,num_absorbing_boundary_faces))
  allocate(absorbing_boundary_jacobian2D(NGLLSQUARE,num_absorbing_boundary_faces))
  allocate(absorbing_boundary_normal(NDIM,NGLLSQUARE,num_absorbing_boundary_faces))
  read(27) absorbing_boundary_ispec
  read(27) absorbing_boundary_ijk
  read(27) absorbing_boundary_jacobian2D
  read(27) absorbing_boundary_normal

! free surface 
  read(27) NSPEC2D_TOP    
  allocate(ibelm_top(NSPEC2D_TOP))
  allocate(jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP))
  allocate(normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP))
  read(27) ibelm_top
  read(27) jacobian2D_top
  read(27) normal_top
  
! MPI interfaces
  read(27) num_interfaces_ext_mesh
  read(27) max_nibool_interfaces_ext_mesh
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh))
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh))
  allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  read(27) my_neighbours_ext_mesh
  read(27) nibool_interfaces_ext_mesh
  read(27) ibool_interfaces_ext_mesh

  allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  allocate(request_send_vector_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_recv_vector_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_send_scalar_ext_mesh(num_interfaces_ext_mesh))
  allocate(request_recv_scalar_ext_mesh(num_interfaces_ext_mesh))
  close(27)

! locate inner and outer elements
  allocate(ispec_is_inner_ext_mesh(NSPEC_AB))
  allocate(iglob_is_inner_ext_mesh(NGLOB_AB))
  ispec_is_inner_ext_mesh(:) = .true.
  iglob_is_inner_ext_mesh(:) = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    do i = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(i,iinterface)
      iglob_is_inner_ext_mesh(iglob) = .false.
    enddo
  enddo
  do ispec = 1, NSPEC_AB
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          ispec_is_inner_ext_mesh(ispec) = iglob_is_inner_ext_mesh(iglob) .and. ispec_is_inner_ext_mesh(ispec)
        enddo
      enddo
    enddo
  enddo

! counts inner and outer elements
!    nspec_inner = 0
!    nspec_outer = 0
!    do ispec = 1, NSPEC_AB
!      if( ispec_is_inner_ext_mesh(ispec) .eqv. .true. ) then
!        nspec_inner = nspec_inner + 1
!      else
!        nspec_outer = nspec_outer + 1
!      endif
!    enddo

! stores indices of inner and outer elements for faster(?) compute_forces_with_Deville routine
!    if( nspec_inner > 0 ) allocate( spec_inner(nspec_inner))
!    if( nspec_outer > 0 ) allocate( spec_outer(nspec_outer))
!    nspec_inner = 0
!    nspec_outer = 0
!    do ispec = 1, NSPEC_AB
!      if( ispec_is_inner_ext_mesh(ispec) .eqv. .true. ) then
!        nspec_inner = nspec_inner + 1
!        spec_inner(nspec_inner) = ispec
!      else
!        nspec_outer = nspec_outer + 1
!        spec_outer(nspec_outer) = ispec
!      endif
!    enddo
  !print *,'rank ',myrank,' inner spec: ',nspec_inner
  !print *,'rank ',myrank,' outer spec: ',nspec_outer

    
! check courant criteria on mesh
  call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    kappastore,mustore,rho_vp,rho_vs, &
                                    DT )

  end subroutine