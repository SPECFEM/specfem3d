!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  subroutine save_databases(nspec,nglob, &
                            iMPIcut_xi,iMPIcut_eta, &
                            nodes_coords,ispec_material_id, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
    SAVE_MESH_AS_CUBIT,NDIM,IMAIN,IIN_DB,myrank

  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,NGNOD,NGNOD2D

  use meshfem3D_par, only: ibool,xstore,ystore,zstore, &
    addressing,NPROC_XI,NPROC_ETA,iproc_xi_current,iproc_eta_current, &
    prname, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NMATERIALS,material_properties, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions

  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec

  ! number of vertices in each block
  integer,intent(in) :: nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  logical,intent(in) :: iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  ! arrays with the mesh
  double precision,intent(in) :: nodes_coords(nglob,NDIM)
  integer,intent(in) :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer,intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer,intent(in) :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer,intent(in) :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer,intent(in) :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer,intent(in) :: ibelm_top(NSPEC2D_TOP)

  ! local parameters
  integer :: nspec_CPML_total,ispec_CPML
  double precision , dimension(17) :: matpropl
  integer :: i,ispec,iglob,ier

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  integer, parameter :: IIN_database = IIN_DB

  integer :: ndef,nundef
  integer :: mat_id,domain_id
  ! there was a potential bug here if nspec is big
  !integer,dimension(2,nspec) :: material_index
  integer,dimension(:,:),allocatable :: material_index
  character(len=MAX_STRING_LEN), dimension(6,1) :: undef_mat_prop

  ! temporary array for local nodes (either NGNOD2D or NGNOD)
  integer, dimension(NGNOD) :: loc_node
  integer, dimension(NGNOD) :: anchor_iax,anchor_iay,anchor_iaz
  integer :: ia,inode

  ! assignes material index
  ! format: (1,ispec) = #material_id , (2,ispec) = #material_definition
  allocate(material_index(2,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1346')
  if (ier /= 0) stop 'Error allocating array material_index'
  material_index (:,:) = 0
  do ispec = 1, nspec
    ! material id
    material_index(1,ispec) = ispec_material_id(ispec)
    ! material definition: 1 = interface type / 2 = tomography type
    if (ispec_material_id(ispec) > 0) then
      ! dummy value, not used any further
      material_index(2,ispec) = 1
    else
      ! negative material ids
      ! by default, assumes tomography model material
      ! (note: interface type not implemented yet...)
      material_index(2,ispec) = 2
    endif
  enddo

  ! Materials properties
  ! counts defined/undefined materials
  ndef = 0
  nundef = 0
  do i = 1,NMATERIALS
    ! material_properties(:,:) array:
    !   first dimension  : material_id
    !   second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    !
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    mat_id = int(material_properties(i,8))
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1
  enddo
  !debug
  !print *,'materials def/undef: ',ndef,nundef

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'mesh files:'
    write(IMAIN,*) '  saving files: proc***_Database'
    call flush_IMAIN()
  endif

  ! opens database file
  open(unit=IIN_database,file=prname(1:len_trim(prname))//'Database', &
       status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening Database file: ',prname(1:len_trim(prname))//'Database'
    stop 'error opening Database file'
  endif

  ! global nodes
  write(IIN_database) nglob
  do iglob = 1,nglob
    ! format: #id #x #y #z
    write(IIN_database) iglob,nodes_coords(iglob,1),nodes_coords(iglob,2),nodes_coords(iglob,3)
  enddo

  ! materials
  ! format: #number of defined materials #number of undefined materials
  write(IIN_database) ndef, nundef

  ! writes out defined materials
  do i = 1,NMATERIALS
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    domain_id = int(material_properties(i,7))
    mat_id = int(material_properties(i,8))
    if (mat_id > 0) then
      ! pad dummy zeros to fill up 17 entries
      matpropl(:) = 0.d0
      select case(domain_id)
      case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC)
        ! material properties format:
        !#(1)rho  #(2)vp #(3)vs #(4)Q_Kappa #(5)Q_mu #(6)anisotropy_flag #(7)domain_id #(8)mat_id
        !
        ! output format for xgenerate_database (same as for cubit/trelis inputs):
        !   rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag,material_domain_id
        !
        ! skipping mat_id, not needed
        matpropl(1:7) = material_properties(i,1:7)
      case (IDOMAIN_POROELASTIC)
        ! material properties format:
        !#(1)rho_s #(2)rho_f #(3)phi #(4)tort #(5)eta #(6)0 #(7)domain_id #(8)mat_id
        !            .. #(9)kxx #(10)kxy #(11)kxz #(12)kyy #(13)kyz #(14)kzz #(15)kappa_s #(16)kappa_f #(17)kappa_fr #(18)mu_fr
        !
        ! output format for xgenerate_database (same as for cubit/trelis inputs):
        !   rhos,rhof,phi,tort,eta,0,material_domain_id,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,mufr
        matpropl(1:7) = material_properties(i,1:7)
        ! skipping mat_id, not needed
        matpropl(8:17) = material_properties(i,9:18)
      end select
      ! writes to database
      write(IIN_database) matpropl(:)
    endif
  enddo

  ! writes out undefined materials
  do i = 1,NMATERIALS
    domain_id = int(material_properties(i,7))
    mat_id = int(material_properties(i,8))
    if (mat_id < 0) then
      ! format:
      ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
      ! format example tomography: -1 tomography elastic tomography_model.xyz 0 2
      undef_mat_prop(:,:) = ''
      ! material id
      write(undef_mat_prop(1,1),*) mat_id
      ! name
      undef_mat_prop(2,1) = 'tomography'
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC)
        undef_mat_prop(3,1) = 'acoustic'
      case (IDOMAIN_ELASTIC)
        undef_mat_prop(3,1) = 'elastic'
      case (IDOMAIN_POROELASTIC)
        undef_mat_prop(3,1) = 'poroelastic'
      end select
      ! default name
      undef_mat_prop(4,1) = 'tomography_model.xyz'
      ! default tomo-id (unused)
      write(undef_mat_prop(5,1),*) 0
      ! domain-id
      write(undef_mat_prop(6,1),*) domain_id
      ! debug
      !print *,'undef mat: ',undef_mat_prop
      ! writes out properties
      write(IIN_database) undef_mat_prop
    endif
  enddo

  ! spectral elements
  !
  ! note: check with routine write_partition_database() to produce identical output
  write(IIN_database) nspec

  ! sets up node addressing
  call hex_nodes_anchor_ijk_NGLL(NGNOD,anchor_iax,anchor_iay,anchor_iaz,NGLLX_M,NGLLY_M,NGLLZ_M)

  do ispec = 1,nspec
    ! format: #ispec #material_id #dummy/tomo_id #iglob1 #iglob2 ..

    ! assumes NGLLX_M == NGLLY_M == NGLLZ_M == 2
    !write(IIN_database) ispec,material_index(1,ispec),material_index(2,ispec), &
    !       ibool(1,1,1,ispec),ibool(2,1,1,ispec),ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(1,1,2,ispec), &
    !       ibool(2,1,2,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)

    ! gets anchor nodes
    do ia = 1,NGNOD
      iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
      loc_node(ia) = iglob
    enddo
    ! output
    write(IIN_database) ispec,material_index(1,ispec),material_index(2,ispec),(loc_node(ia),ia = 1,NGNOD)
  enddo

  ! Boundaries
  !
  ! note: check with routine write_boundaries_database() to produce identical output
  write(IIN_database) 1,nspec2D_xmin
  write(IIN_database) 2,nspec2D_xmax
  write(IIN_database) 3,nspec2D_ymin
  write(IIN_database) 4,nspec2D_ymax
  write(IIN_database) 5,NSPEC2D_BOTTOM
  write(IIN_database) 6,NSPEC2D_TOP

  do i = 1,nspec2D_xmin
    ispec = ibelm_xmin(i)
    ! gets anchor nodes on xmin
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iax(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for xmin'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for xmin'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! assumes NGLLX_M == NGLLY_M == NGLLZ_M == 2
    !write(IIN_database) ispec,ibool(1,1,1,ispec),ibool(1,NGLLY_M,1,ispec), &
    !                          ibool(1,1,NGLLZ_M,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec)
  enddo

  do i = 1,nspec2D_xmax
    ispec = ibelm_xmax(i)
    ! gets anchor nodes on xmax
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iax(ia) == NGLLX_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for xmax'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for xmax'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), &
    !      ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
  enddo

  do i = 1,nspec2D_ymin
    ispec = ibelm_ymin(i)
    ! gets anchor nodes on ymin
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iay(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for ymin'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for ymin'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)), &
    !      ibool(1,1,NGLLZ_M,ibelm_ymin(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
  enddo

  do i = 1,nspec2D_ymax
    ispec = ibelm_ymax(i)
    ! gets anchor nodes on ymax
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iay(ia) == NGLLY_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for ymax'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for ymax'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_ymax(i),ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i)),ibool(1,NGLLY_M,1,ibelm_ymax(i)), &
    !      ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
  enddo

  do i = 1,NSPEC2D_BOTTOM
    ispec = ibelm_bottom(i)
    ! gets anchor nodes on bottom
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iaz(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for bottom'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for bottom'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_bottom(i),ibool(1,1,1,ibelm_bottom(i)),ibool(NGLLX_M,1,1,ibelm_bottom(i)), &
    !      ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i)),ibool(1,NGLLY_M,1,ibelm_bottom(i))
  enddo

  do i = 1,NSPEC2D_TOP
    ispec = ibelm_top(i)
    ! gets anchor nodes on top
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iaz(ia) == NGLLZ_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for top'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for top'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_top(i),ibool(1,1,NGLLZ_M,ibelm_top(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i)), &
    !      ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
  enddo

  ! CPML
  !
  ! note: check with routine write_cpml_database() to produce identical output
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call synchronize_all()
  call bcast_all_singlei(nspec_CPML_total)
  call synchronize_all()

  write(IIN_database) nspec_CPML_total

  if (nspec_CPML_total > 0) then
     write(IIN_database) nspec_CPML

     do ispec_CPML = 1,nspec_CPML
        write(IIN_database) CPML_to_spec(ispec_CPML), CPML_regions(ispec_CPML)
     enddo
     do ispec = 1,nspec
        write(IIN_database) is_CPML(ispec)
     enddo
  endif

  ! MPI Interfaces
  !
  ! note: check with routine write_interfaces_database() to produce identical output
  if (NPROC_XI > 1 .or. NPROC_ETA > 1) then
    ! determines number of MPI interfaces for each slice
    nb_interfaces = 4
    interfaces(W:N) = .true.
    interfaces(NW:SW) = .false.

    ! slices at model boundaries
    if (iproc_xi_current == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(W) = .false.
    endif
    if (iproc_xi_current == NPROC_XI-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(E) = .false.
    endif
    if (iproc_eta_current == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(S) = .false.
    endif
    if (iproc_eta_current == NPROC_ETA-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(N) = .false.
    endif

    ! slices in middle of model
    if ((interfaces(W) .eqv. .true.) .and. (interfaces(N) .eqv. .true.)) then
      interfaces(NW) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(N) .eqv. .true.) .and. (interfaces(E) .eqv. .true.)) then
      interfaces(NE) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(E) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
      interfaces(SE) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(W) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
      interfaces(SW) = .true.
      nb_interfaces =  nb_interfaces +1
    endif

    nspec_interface(:) = 0
    if (interfaces(W))  nspec_interface(W) = count(iMPIcut_xi(1,:) .eqv. .true.)
    if (interfaces(E))  nspec_interface(E) = count(iMPIcut_xi(2,:) .eqv. .true.)
    if (interfaces(S))  nspec_interface(S) = count(iMPIcut_eta(1,:) .eqv. .true.)
    if (interfaces(N))  nspec_interface(N) = count(iMPIcut_eta(2,:) .eqv. .true.)
    if (interfaces(NW))  nspec_interface(NW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(NE))  nspec_interface(NE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(SE))  nspec_interface(SE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))
    if (interfaces(SW))  nspec_interface(SW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))

    nspec_interfaces_max = maxval(nspec_interface)

    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    write(IIN_database) nb_interfaces,nspec_interfaces_max

    ! face elements
    if (interfaces(W)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current-1,iproc_eta_current),nspec_interface(W)
      do ispec = 1,nspec
        if (iMPIcut_xi(1,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          ! note: face outputs 4 corner points
          write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(1,NGLLY_M,1,ispec), &
                                      ibool(1,1,NGLLZ_M,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec)
        endif
      enddo
    endif

    if (interfaces(E)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current+1,iproc_eta_current),nspec_interface(E)
      do ispec = 1,nspec
        if (iMPIcut_xi(2,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,4,ibool(NGLLX_M,1,1,ispec),ibool(NGLLX_M,NGLLY_M,1,ispec), &
                                      ibool(NGLLX_M,1,NGLLZ_M,ispec),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
        endif
      enddo
    endif

    if (interfaces(S)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current,iproc_eta_current-1),nspec_interface(S)
      do ispec = 1,nspec
        if (iMPIcut_eta(1,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(NGLLX_M,1,1,ispec), &
                                      ibool(1,1,NGLLZ_M,ispec),ibool(NGLLX_M,1,NGLLZ_M,ispec)
        endif
      enddo
    endif

    if (interfaces(N)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current,iproc_eta_current+1),nspec_interface(N)
      do ispec = 1,nspec
        if (iMPIcut_eta(2,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,4,ibool(NGLLX_M,NGLLY_M,1,ispec),ibool(1,NGLLY_M,1,ispec), &
                                      ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec)
        endif
      enddo
    endif

    ! edge elements
    if (interfaces(NW)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current-1,iproc_eta_current+1),nspec_interface(NW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          ! note: edge elements output 2 corners and 2 dummy values
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(1,NGLLY_M,1,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(NE)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current+1,iproc_eta_current+1),nspec_interface(NE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(NGLLX_M,NGLLY_M,1,ispec),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SE)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current+1,iproc_eta_current-1),nspec_interface(SE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(NGLLX_M,1,1,ispec),ibool(NGLLX_M,1,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SW)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current-1,iproc_eta_current-1),nspec_interface(SW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(1,1,1,ispec),ibool(1,1,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif
  else
    ! single process execution, no MPI boundaries
    nb_interfaces = 0
    nspec_interfaces_max = 0
    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    write(IIN_database) nb_interfaces,nspec_interfaces_max
  endif

  close(IIN_database)

  ! CUBIT output
  if (SAVE_MESH_AS_CUBIT) then
    ! only for single process at the moment
    if (NPROC_XI == 1 .and. NPROC_ETA == 1) then
      !! VM VM add outputs as CUBIT
      call save_output_mesh_files_as_cubit(nspec,nglob, &
                                           nodes_coords, ispec_material_id, &
                                           nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                           ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                           nspec_CPML_total)
      ! output for AxiSEM coupling
      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        call save_output_mesh_files_for_coupled_model(nspec, &
                                                      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                                      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                                      xstore,ystore,zstore)
      endif
    endif
  endif

  deallocate(material_index)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done mesh files'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_databases

!---------------------------------------------------------------------------------------------------------------

  !! VM VM add subroutine to save meshes in case of a single MPI process
  subroutine save_output_mesh_files_as_cubit(nspec,nglob, &
                                             nodes_coords, ispec_material_id, &
                                             nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                             ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                             nspec_CPML_total)

  use constants, only: NDIM,IMAIN,myrank,IIN_DB
  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M
  use shared_parameters, only: NGNOD

  use meshfem3D_par, only: ibool, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NMATERIALS,material_properties, &
    nspec_CPML,CPML_to_spec,CPML_regions

  implicit none

  integer, parameter :: IIN_database = IIN_DB

  ! number of spectral elements in each block
  integer, intent(in):: nspec

  ! number of vertices in each block
  integer, intent(in) :: nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  !integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  ! arrays with the mesh
  double precision, intent(in) :: nodes_coords(nglob,NDIM)

  integer, intent(in) :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer, intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer, intent(in) :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer, intent(in) :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer, intent(in) :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer, intent(in) :: ibelm_top(NSPEC2D_TOP)

  integer, intent(in) :: nspec_CPML_total

  ! local parameters
  integer :: i,ispec,iglob,ier
  integer :: mat_id,domain_id
  double precision  :: z_bottom

  integer :: nnodes_mesh,inode
  integer, dimension(:),allocatable :: iglob_to_nodeid
  logical, dimension(:), allocatable :: mask_iglob

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  saving mesh files as cubit in directory: MESH/'
    if (NGNOD /= 8) then
      write(IMAIN,*) '    (Using these cubit files will require using NGNOD = 8 in Par_file for decomposer)'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  z_bottom = 0.d0 ! will shift coordinates in z-direction

  ! nodes needed by ibool for outputting anchor nodes
  allocate(mask_iglob(nglob),iglob_to_nodeid(nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating mask_iglob'
  mask_iglob(:) = .false.
  iglob_to_nodeid(:) = 0

  ! only output corner points as mesh, no internal points, otherwise decomposer will complain about unused nodes
  do ispec = 1,nspec
    ! corner points
    mask_iglob(ibool(1,1,1,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,1,1,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,NGLLZ_M,1,ispec)) = .true.
    mask_iglob(ibool(1,NGLLZ_M,1,ispec)) = .true.
    mask_iglob(ibool(1,1,NGLLZ_M,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,1,NGLLZ_M,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,NGLLZ_M,NGLLZ_M,ispec)) = .true.
    mask_iglob(ibool(1,NGLLZ_M,NGLLZ_M,ispec)) = .true.
  enddo

  ! corner points only
  nnodes_mesh = count(mask_iglob(:))
  if (nnodes_mesh < 1 .or. nnodes_mesh > nglob) then
    print *,'Error: nnodes_mesh ',nnodes_mesh,' is invalid for nglob = ',nglob
    stop 'Error nnodes_mesh invalid'
  endif

  ! maps iglob to new unique node IDs
  inode = 0
  do iglob = 1,nglob
    if (mask_iglob(iglob)) then
      inode = inode + 1
      iglob_to_nodeid(iglob) = inode
    endif
  enddo
  if (inode /= nnodes_mesh) then
    print *,'Error: inode count ',inode,' for nnodes_mesh ',nnodes_mesh,' is invalid; nglob = ',nglob
    stop 'Error inode count invalid'
  endif

  ! safety check just to re-visit this routine in case NGNOD changes
  if (NGNOD /= 8 .and. NGNOD /= 27) then
    stop 'Error invalid NGNOD for routine save_output_mesh_files_as_cubit()'
  endif

  ! outputs mesh as files in MESH/ directory
  ! (used for CUBIT models stored in a specfem-readable way; will need to run xdecompose_mesh with these files to continue)

  open(IIN_database, file = 'MESH/nummaterial_velocity_file',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ','MESH/nummaterial_velocity_file'
    print *,'Please check if directory MESH/ exists for saving mesh files as cubit...'
    stop 'Error opening mesh file'
  endif

  do i = 1,NMATERIALS
    domain_id = int(material_properties(i,7))
    mat_id =  int(material_properties(i,8))
    if (domain_id > 0) then
      ! format: #domain_id #material_id #rho #vp #vs #Qkappa #Qmu #anisotropy_flag
      write(IIN_database,'(2i6,5f15.5,i6)') domain_id,mat_id,material_properties(i,1:5),0
    else
      write(*,*) 'STOP: undefined mat not yet implemented'
      stop
    endif
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/materials_file')
  do ispec = 1, nspec
    ! format: #ispec #material_id
    write(IIN_database,*) ispec,ispec_material_id(ispec)
  enddo

  open(IIN_database,file='MESH/nodes_coords_file')
  write(IIN_database,*) nnodes_mesh
  do iglob = 1,nglob
    ! point coordinates
    if (mask_iglob(iglob)) then
      ! format: #iglob #x #y #z(corrected by z_bottom)
      write(IIN_database,'(i14,3x,3(f20.5,1x))') iglob_to_nodeid(iglob), &
                                                 nodes_coords(iglob,1), &
                                                 nodes_coords(iglob,2), &
                                                 nodes_coords(iglob,3)-z_bottom
    endif
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/mesh_file')
  write(IIN_database,*) nspec
  do ispec = 1,nspec
    ! corner points
    ! format: #ispec #iglob1 #iglob2 #iglob3 #iglob4 #iglob5 #iglob6 #iglob7 #iglob8
    write(IIN_database,'(9i15)')  ispec, &
                                  iglob_to_nodeid(ibool(1,1,1,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,1,1,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,NGLLZ_M,1,ispec)), &
                                  iglob_to_nodeid(ibool(1,NGLLZ_M,1,ispec)), &
                                  iglob_to_nodeid(ibool(1,1,NGLLZ_M,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,1,NGLLZ_M,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,NGLLZ_M,NGLLZ_M,ispec)), &
                                  iglob_to_nodeid(ibool(1,NGLLZ_M,NGLLZ_M,ispec))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_xmin')
  write(IIN_database,*) nspec2D_xmin
  do i = 1,nspec2D_xmin
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_xmin(i), &
                                      iglob_to_nodeid(ibool(1,1,1,ibelm_xmin(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,1,ibelm_xmin(i))), &
                                      iglob_to_nodeid(ibool(1,1,NGLLZ_M,ibelm_xmin(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_xmax')
  write(IIN_database,*) nspec2D_xmax
  do i = 1,nspec2D_xmax
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_xmax(i), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,1,ibelm_xmax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_ymin')
  write(IIN_database,*) nspec2D_ymin
  do i = 1,nspec2D_ymin
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_ymin(i), &
                                      iglob_to_nodeid(ibool(1,1,1,ibelm_ymin(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,1,ibelm_ymin(i))), &
                                      iglob_to_nodeid(ibool(1,1,NGLLZ_M,ibelm_ymin(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_ymax')
  write(IIN_database,*) nspec2D_ymax
  do i = 1,nspec2D_ymax
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_ymax(i), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,1,ibelm_ymax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i)))
  enddo


  open(IIN_database,file='MESH/absorbing_surface_file_bottom')
  write(IIN_database,*) NSPEC2D_BOTTOM
  do i = 1,NSPEC2D_BOTTOM
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_bottom(i), &
                                      iglob_to_nodeid(ibool(1,1,1,ibelm_bottom(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,1,ibelm_bottom(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,1,ibelm_bottom(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/free_or_absorbing_surface_file_zmax')
  write(IIN_database,*) NSPEC2D_TOP
  do i = 1,NSPEC2D_TOP
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_top(i), &
                                      iglob_to_nodeid(ibool(1,1,NGLLZ_M,ibelm_top(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i)))
  enddo
  close(IIN_database)

  if (nspec_CPML_total > 0) then
    open(IIN_database,file='MESH/absorbing_cpml_file')
    write(IIN_database,*) nspec_CPML
    do i = 1,nspec_CPML
      ! format: #cpml_ispec #cpml_region_id
      write(IIN_database,*) CPML_to_spec(i), CPML_regions(i)
    enddo
    close(IIN_database)
  endif

  deallocate(mask_iglob,iglob_to_nodeid)

  end subroutine save_output_mesh_files_as_cubit



!---------------------------------------------------------------------------------------------------------------

  !! VM VM add subroutine to save meshes in case of a single MPI process
  subroutine save_output_mesh_files_for_coupled_model(nspec, &
                                                      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                                      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                                      xgrid,ygrid,zgrid)

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM, ZERO, IMAIN, myrank, &
    INJECTION_TECHNIQUE_IS_AXISEM
  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use shared_parameters, only: NGNOD,COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE

  use meshfem3D_par, only: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  implicit none

  ! number of spectral elements in each block
  integer :: nspec

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  !integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  !! VM VM add all GLL points for Axisem coupling
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xgrid, ygrid, zgrid

  ! boundary parameters locator
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer :: ibelm_top(NSPEC2D_TOP)

  integer :: i,ispec

  double precision  :: z_bottom

  ! for axisem coupling case  ( only serial case for mesher use scotch after)
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)
  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0
  double precision   :: rotation_matrix(3,3)
  double precision   :: zlayer(nlayer), vpv(nlayer,4), vsv(nlayer,4), density(nlayer,4)
  integer            :: ilayer, updown(NGLLZ)

  !! GLL points
  double precision, dimension(:,:,:), allocatable  ::  longitud, latitud, radius
  double precision, dimension(:,:,:), allocatable  ::  xstore, ystore, zstore
   !! Element control points
  double precision, dimension(:), allocatable :: xelm, yelm, zelm
  double precision, dimension(NGLLZ) :: radius_Z ! for temporary copies

  !! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable    :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable  :: dershape3D
  !! GLL points and weights of integration
  double precision, dimension(:), allocatable  :: xigll, yigll, zigll, wxgll, wygll, wzgll

  double precision  :: ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD
  double precision  :: lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision  :: radius_of_box_top

  integer :: ielm, j,k, imin,imax,jmin,jmax,kmin,kmax,ier
  integer :: nel_lat, nel_lon, nel_depth
  logical :: buried_box

  character(len=10)  :: line
  character(len=250) :: model1D_file

  double precision,parameter  :: deg2rad = 3.141592653589793d0/180.d0

  ! safety check
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

!! VM VM add files in case of AxiSEM coupling
  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  saving mesh files for coupled model in directory: MESH/'
      call flush_IMAIN()
    endif

1000 format(3f30.10)

    z_bottom = 0.d0 ! will shift coordinates in z-direction

    allocate(longitud(NGLLX,NGLLY,NGLLZ), latitud(NGLLX,NGLLY,NGLLZ), radius(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1347')

    allocate(xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1348')

    allocate(xelm(NGNOD), yelm(NGNOD), zelm(NGNOD),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1349')

    allocate(xigll(NGLLX), yigll(NGLLY), zigll(NGLLZ), wxgll(NGLLX),wygll(NGLLY), wzgll(NGLLZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1350')

    !
    !--- set up coordinates of the Gauss-Lobatto-Legendre points
    !

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

    !
    !--- if number of points is odd, the middle abscissa is exactly zero
    !
    if (mod(NGLLX,2) /= 0) xigll((NGLLX - 1)/2 + 1) = ZERO
    if (mod(NGLLY,2) /= 0) yigll((NGLLY - 1)/2 + 1) = ZERO
    if (mod(NGLLZ,2) /= 0) zigll((NGLLZ - 1)/2 + 1) = ZERO

    !
    !--- get the 3-D shape functions
    !
    allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ),dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1351')

    call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ)

    !! reading parameters for coupling

    open(90, file='MESH/ParFileMeshChunk',action='read')
    read(90,'(a)') line
    read(90,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
    read(90,'(a)') line
    read(90,*) lon_center_chunk, lat_center_chunk, chunk_azi
    read(90,'(a)') line
    read(90,*) chunk_depth
    read(90,'(a)') line
    read(90,*) nel_lon,nel_lat, nel_depth
    read(90,'(a)') line
    read(90,'(a)') model1D_file
    read(90,'(a)') line
    read(90,*) buried_box
    if (buried_box) then
      read(90,'(a)') line
      read(90,*) radius_of_box_top
       radius_of_box_top =  radius_of_box_top * 1000.
    else
      radius_of_box_top = 6371000.
    endif
    model1D_file = 'MESH/'//trim(model1D_file)
    close(90)

    ! read 1D AxiSEM model
    call Read_dsm_model(model1D_file,vpv,vsv,density,zlayer,nlayer)

    ! modele 1D
    open(88,file='MESH/model_1D.in')
    write(88,*) nlayer,4
    do i=1,nlayer
      write(88,*) zlayer(i)
      write(88,'(4f20.10)') vpv(i,:)
      write(88,'(4f20.10)') vsv(i,:)
      write(88,'(4f20.10)') density(i,:)
    enddo
    z_bottom = minval(zgrid(:,:,:,:))
    write(88,*)  radius_of_box_top + z_bottom!6371000.+z_bottom
    write(88,*)  lon_center_chunk,  lat_center_chunk,  chunk_azi
    close(88)

    ! compute rotation matrix
    call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)

    open(91, file = 'MESH/list_ggl_boundary_spherical.txt')
    open(92, file = 'MESH/list_ggl_boundary_Cartesian.txt')
    open(89, file = 'MESH/flags_boundary.txt')

    open(90, file = 'MESH/Nb_ielm_faces.txt')
    write(90,*)  nspec2D_xmin
    write(90,*)  nspec2D_xmax
    write(90,*)  nspec2D_ymin
    write(90,*)  nspec2D_ymax
    write(90,*)  nspec2D_bottom
    write(90,*)  nspec2D_top
    close(90)

    ! xmin
    do ielm=1,nspec2D_xmin

      ispec=ibelm_xmin(ielm)

      write(89,*) ispec,ielm,1

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(NGLLX_M,1,1,ispec)
      xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
      xelm(4)=xgrid(1,NGLLY_M,1,ispec)
      xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
      xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
      xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(NGLLX_M,1,1,ispec)
      yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
      yelm(4)=ygrid(1,NGLLY_M,1,ispec)
      yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
      yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
      yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(NGLLX_M,1,1,ispec)
      zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
      zelm(4)=zgrid(1,NGLLY_M,1,ispec)
      zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
      zelm(6)=zgrid(NGLLX_M,1,NGLLY_M,ispec)
      zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top ! 6371000.
      radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
      call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

      imin = 1
      imax = 1
      jmin = 1
      jmax = NGLLY
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,1, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! xmax
    do ielm=1,nspec2D_xmax

      ispec=ibelm_xmax(ielm)

      write(89,*) ispec,ielm,2

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(NGLLX_M,1,1,ispec)
      xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
      xelm(4)=xgrid(1,NGLLY_M,1,ispec)
      xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
      xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
      xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(NGLLX_M,1,1,ispec)
      yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
      yelm(4)=ygrid(1,NGLLY_M,1,ispec)
      yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
      yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
      yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(NGLLX_M,1,1,ispec)
      zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
      zelm(4)=zgrid(1,NGLLY_M,1,ispec)
      zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
      zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
      zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top ! 6371000.
      radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
      call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

      imin = NGLLX
      imax = NGLLX
      jmin = 1
      jmax = NGLLY
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,2, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! ymin
    do ielm=1,nspec2D_ymin

      ispec=ibelm_ymin(ielm)

      write(89,*) ispec,ielm,3

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(NGLLX_M,1,1,ispec)
      xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
      xelm(4)=xgrid(1,NGLLY_M,1,ispec)
      xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
      xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
      xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(NGLLX_M,1,1,ispec)
      yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
      yelm(4)=ygrid(1,NGLLY_M,1,ispec)
      yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
      yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
      yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(NGLLX_M,1,1,ispec)
      zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
      zelm(4)=zgrid(1,NGLLY_M,1,ispec)
      zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
      zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
      zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top !6371000.
      radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
      call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = 1
      jmax = 1
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,3, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! ymax
    do ielm=1,nspec2D_ymax

      ispec=ibelm_ymax(ielm)

      write(89,*) ispec,ielm,4

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(NGLLX_M,1,1,ispec)
      xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
      xelm(4)=xgrid(1,NGLLY_M,1,ispec)
      xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
      xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
      xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(NGLLX_M,1,1,ispec)
      yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
      yelm(4)=ygrid(1,NGLLY_M,1,ispec)
      yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
      yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
      yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(NGLLX_M,1,1,ispec)
      zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
      zelm(4)=zgrid(1,NGLLY_M,1,ispec)
      zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
      zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
      zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top !6371000.
      radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
      call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = NGLLY
      jmax = NGLLY
      kmin = 1
      kmax = NGLLZ

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,4, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    ! bottom
    do ielm=1,nspec2D_BOTTOM

      ispec=ibelm_bottom(ielm)

      write(89,*) ispec,ielm,5

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(NGLLX_M,1,1,ispec)
      xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
      xelm(4)=xgrid(1,NGLLY_M,1,ispec)
      xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
      xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
      xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(NGLLX_M,1,1,ispec)
      yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
      yelm(4)=ygrid(1,NGLLY_M,1,ispec)
      yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
      yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
      yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(NGLLX_M,1,1,ispec)
      zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
      zelm(4)=zgrid(1,NGLLY_M,1,ispec)
      zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
      zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
      zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top ! 6371000.
      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top ! 6371000.
      radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
      call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = 1
      jmax = NGLLY
      kmin = 1
      kmax = 1

      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,5, &
                    ilayer,updown(k)
               write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            enddo
         enddo
      enddo
    enddo

    if (buried_box) then
      ! top
      do ielm=1,nspec2D_TOP

         ispec=ibelm_top(ielm)

         write(89,*) ispec,ielm,6

         xelm(1)=xgrid(1,1,1,ispec)
         xelm(2)=xgrid(NGLLX_M,1,1,ispec)
         xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
         xelm(4)=xgrid(1,NGLLY_M,1,ispec)
         xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
         xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
         xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
         xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

         yelm(1)=ygrid(1,1,1,ispec)
         yelm(2)=ygrid(NGLLX_M,1,1,ispec)
         yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
         yelm(4)=ygrid(1,NGLLY_M,1,ispec)
         yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
         yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
         yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
         yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

         zelm(1)=zgrid(1,1,1,ispec)
         zelm(2)=zgrid(NGLLX_M,1,1,ispec)
         zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
         zelm(4)=zgrid(1,NGLLY_M,1,ispec)
         zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
         zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
         zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
         zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

         call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
         zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top !6371000.
         call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
         zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top !6371000.
         radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
         call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

         imin = 1
         imax = NGLLX
         jmin = 1
         jmax = NGLLY
         kmin = NGLLZ
         kmax = NGLLZ

         do k=kmin,kmax
            do j=jmin,jmax
               do i=imin,imax
                  write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,6, &
                       ilayer,updown(k)
                  write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
               enddo
            enddo
         enddo
      enddo
    endif

    close(89)
    close(91)
    close(92)

    deallocate(shape3D,dershape3D)
  endif

  end subroutine save_output_mesh_files_for_coupled_model


