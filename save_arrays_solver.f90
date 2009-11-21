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


! for external mesh 

  subroutine save_arrays_solver_ext_mesh(nspec,nglob, &
                    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                    gammaxstore,gammaystore,gammazstore, &
                    jacobianstore, rho_vp,rho_vs,iflag_attenuation_store, &
                    rhostore,kappastore,mustore, &
                    rmass,rmass_acoustic,rmass_solid_poroelastic,rmass_fluid_poroelastic, &
                    OCEANS,rmass_ocean_load,NGLOB_OCEAN,&
                    ibool, &
                    xstore_dummy,ystore_dummy,zstore_dummy, &
                    abs_boundary_normal,abs_boundary_jacobian2Dw, &
                    abs_boundary_ijk,abs_boundary_ispec, &
                    num_abs_boundary_faces, &
                    free_surface_normal,free_surface_jacobian2Dw, &
                    free_surface_ijk,free_surface_ispec, &
                    num_free_surface_faces, &
                    coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                    coupling_ac_el_ijk,coupling_ac_el_ispec, &
                    num_coupling_ac_el_faces, &
                    num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    prname,SAVE_MESH_FILES, &
                    ANISOTROPY,NSPEC_ANISO, &
                    c11store,c12store,c13store,c14store,c15store,c16store, &
                    c22store,c23store,c24store,c25store,c26store,c33store, &
                    c34store,c35store,c36store,c44store,c45store,c46store, &
                    c55store,c56store,c66store, &
                    ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  implicit none

  include "constants.h"

  integer :: nspec,nglob

! jacobian  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xixstore,xiystore,xizstore, &
            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rho_vp,rho_vs

! attenuation
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: iflag_attenuation_store

! material
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore,kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass,rmass_acoustic, &
            rmass_solid_poroelastic,rmass_fluid_poroelastic
! ocean load
  logical :: OCEANS
  integer :: NGLOB_OCEAN
  real(kind=CUSTOM_REAL),dimension(NGLOB_OCEAN) :: rmass_ocean_load

! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  
! absorbing boundary surface  
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces) 
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces) 
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces) 
  
! free surface
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL) :: free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces)  
  real(kind=CUSTOM_REAL) :: free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

! acoustic-elastic coupling surface
  integer :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces) 
  real(kind=CUSTOM_REAL) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces) 
  integer :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ispec(num_coupling_ac_el_faces)   

! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

! file name
  character(len=256) prname
  logical :: SAVE_MESH_FILES

! anisotropy
  logical :: ANISOTROPY
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store

! material domain flags
  logical, dimension(nspec) :: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic
  
! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: v_tmp
  integer,dimension(:),allocatable :: v_tmp_i
  
  !real(kind=CUSTOM_REAL) :: minimum(1)
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: ier,i  
  logical :: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION
  character(len=256) :: filename

  integer, dimension(:), allocatable :: iglob_tmp
  integer :: j,inum
  
! saves mesh file proc***_external_mesh.bin
  filename = prname(1:len_trim(prname))//'external_mesh.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening database proc######_external_mesh.bin'
  
  write(IOUT) nspec
  write(IOUT) nglob

  write(IOUT) xixstore
  write(IOUT) xiystore
  write(IOUT) xizstore
  write(IOUT) etaxstore
  write(IOUT) etaystore
  write(IOUT) etazstore
  write(IOUT) gammaxstore
  write(IOUT) gammaystore
  write(IOUT) gammazstore
  write(IOUT) jacobianstore

  write(IOUT) ibool

  write(IOUT) xstore_dummy
  write(IOUT) ystore_dummy
  write(IOUT) zstore_dummy

  write(IOUT) kappastore
  write(IOUT) mustore

  write(IOUT) ispec_is_acoustic
  write(IOUT) ispec_is_elastic
  write(IOUT) ispec_is_poroelastic

! acoustic
! all processes will have acoustic_simulation set if any flag is .true. somewhere
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  if( ACOUSTIC_SIMULATION ) then    
    write(IOUT) rmass_acoustic
    write(IOUT) rhostore
  endif

! elastic
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  if( ELASTIC_SIMULATION ) then
    write(IOUT) rmass
    if( OCEANS) then
      write(IOUT) rmass_ocean_load
    endif
    !pll Stacey 
    write(IOUT) rho_vp
    write(IOUT) rho_vs
    write(IOUT) iflag_attenuation_store
  endif

! poroelastic  
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )  
  if( POROELASTIC_SIMULATION ) then
    write(IOUT) rmass_solid_poroelastic
    write(IOUT) rmass_fluid_poroelastic
  endif

! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  write(IOUT) abs_boundary_ispec
  write(IOUT) abs_boundary_ijk
  write(IOUT) abs_boundary_jacobian2Dw
  write(IOUT) abs_boundary_normal

! free surface 
  write(IOUT) num_free_surface_faces
  write(IOUT) free_surface_ispec
  write(IOUT) free_surface_ijk
  write(IOUT) free_surface_jacobian2Dw
  write(IOUT) free_surface_normal

! acoustic-elastic coupling surface
  write(IOUT) num_coupling_ac_el_faces
  write(IOUT) coupling_ac_el_ispec   
  write(IOUT) coupling_ac_el_ijk
  write(IOUT) coupling_ac_el_jacobian2Dw 
  write(IOUT) coupling_ac_el_normal 

!MPI interfaces
  write(IOUT) num_interfaces_ext_mesh
  write(IOUT) maxval(nibool_interfaces_ext_mesh(:))
  write(IOUT) my_neighbours_ext_mesh
  write(IOUT) nibool_interfaces_ext_mesh

  allocate(ibool_interfaces_ext_mesh_dummy(maxval(nibool_interfaces_ext_mesh(:)),num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'
  
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:maxval(nibool_interfaces_ext_mesh(:)),i)
  enddo
  write(IOUT) ibool_interfaces_ext_mesh_dummy

! anisotropy
  if( ANISOTROPY ) then
    write(IOUT) c11store
    write(IOUT) c12store
    write(IOUT) c13store
    write(IOUT) c14store
    write(IOUT) c15store
    write(IOUT) c16store
    write(IOUT) c22store
    write(IOUT) c23store
    write(IOUT) c24store
    write(IOUT) c25store
    write(IOUT) c26store
    write(IOUT) c33store
    write(IOUT) c34store
    write(IOUT) c35store
    write(IOUT) c36store
    write(IOUT) c44store
    write(IOUT) c45store
    write(IOUT) c46store
    write(IOUT) c55store
    write(IOUT) c56store
    write(IOUT) c66store
  endif

  close(IOUT)


! stores arrays in binary files
  if( SAVE_MESH_FILES ) then
    
    ! mesh arrays used for example in combine_vol_data.f90
    !--- x coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted')
    write(27) xstore_dummy
    close(27)

    !--- y coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted')
    write(27) ystore_dummy
    close(27)

    !--- z coordinate
    open(unit=27,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted')
    write(27) zstore_dummy
    close(27)

    ! ibool
    open(unit=27,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted')
    write(27) ibool
    close(27)

    allocate( v_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier); if( ier /= 0 ) stop 'error allocating array '

    ! vp (for checking the mesh and model)  
    !minimum = minval( abs(rho_vp) )
    !if( minimum(1) /= 0.0 ) then
    !  v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
    !else
    !  v_tmp = 0.0
    !endif  
    v_tmp = 0.0
    where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp    
    open(unit=27,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

    ! VTK file output    
    ! vp values
    filename = prname(1:len_trim(prname))//'vp'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp,filename)


    ! vs (for checking the mesh and model)
    !minimum = minval( abs(rho_vs) )
    !if( minimum(1) /= 0.0 ) then
    !  v_tmp = mustore / rho_vs
    !else  
    !  v_tmp = 0.0
    !endif
    v_tmp = 0.0
    where( rho_vs /= 0._CUSTOM_REAL )  v_tmp = mustore / rho_vs    
    open(unit=27,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

    ! VTK file output    
    ! vs values
    filename = prname(1:len_trim(prname))//'vs'
    call write_VTK_data_gll_cr(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp,filename)

    ! VTK file output
    ! saves attenuation flag assigned on each gll point into a vtk file 
    filename = prname(1:len_trim(prname))//'attenuation_flag'
    call write_VTK_data_gll_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        iflag_attenuation_store,&
                        filename)
    ! VTK file output  
    ! acoustic-elastic domains    
    if( ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION ) then
      ! saves points on acoustic-elastic coupling interface
      allocate( iglob_tmp(NGLLSQUARE*num_coupling_ac_el_faces))
      inum = 0
      iglob_tmp(:) = 0
      do i=1,num_coupling_ac_el_faces
        do j=1,NGLLSQUARE
          inum = inum+1
          iglob_tmp(inum) = ibool(coupling_ac_el_ijk(1,j,i), &
                                  coupling_ac_el_ijk(2,j,i), &
                                  coupling_ac_el_ijk(3,j,i), &
                                  coupling_ac_el_ispec(i) )
        enddo
      enddo
      filename = prname(1:len_trim(prname))//'coupling_acoustic_elastic'      
      call write_VTK_data_points(nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        iglob_tmp,NGLLSQUARE*num_coupling_ac_el_faces, &
                        filename)
      
      ! saves acoustic/elastic flag    
      allocate(v_tmp_i(nspec))                                  
      do i=1,nspec
        if( ispec_is_acoustic(i) ) then
          v_tmp_i(i) = 1
        else if( ispec_is_elastic(i) ) then
          v_tmp_i(i) = 2
        else
          v_tmp_i(i) = 0
        endif
      enddo
      filename = prname(1:len_trim(prname))//'acoustic_elastic_flag'
      call write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        v_tmp_i,filename)
    endif

    !! saves 1. MPI interface
    !    if( num_interfaces_ext_mesh >= 1 ) then
    !      filename = prname(1:len_trim(prname))//'MPI_1_points'
    !      call write_VTK_data_points(nglob, &
    !                        xstore_dummy,ystore_dummy,zstore_dummy, &
    !                        ibool_interfaces_ext_mesh_dummy(1:nibool_interfaces_ext_mesh(1),1), &
    !                        nibool_interfaces_ext_mesh(1), &
    !                        filename)
    !    endif
    !    

    deallocate(v_tmp)
    
  endif ! SAVE_MESH_FILES

! cleanup
  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier); if( ier /= 0 ) stop 'error deallocating array'


  end subroutine save_arrays_solver_ext_mesh
  
  

!=============================================================
!
!! old way
!! regular mesh
!
!  subroutine save_arrays_solver(flag_sediments,not_fully_in_bedrock,rho_vp,rho_vs,prname,xixstore,xiystore,xizstore, &
!            etaxstore,etaystore,etazstore, &
!            gammaxstore,gammaystore,gammazstore,jacobianstore, &
!            xstore,ystore,zstore,kappastore,mustore, &
!            ANISOTROPY, &
!            c11store,c12store,c13store,c14store,c15store,c16store, &
!            c22store,c23store,c24store,c25store,c26store,c33store,c34store,c35store,c36store, &
!            c44store,c45store,c46store,c55store,c56store,c66store, &
!            ibool,idoubling,rmass,rmass_ocean_load,npointot_oceans, &
!            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
!            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
!            normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
!            jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
!            jacobian2D_bottom,jacobian2D_top, &
!            iMPIcut_xi,iMPIcut_eta,nspec,nglob, &
!            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP,OCEANS)
!
!  implicit none
!
!  include "constants.h"
!
!  integer nspec,nglob
!  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP
!  integer npointot_oceans
!
!  logical OCEANS
!  logical ANISOTROPY
!
!! arrays with jacobian matrix
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
!    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
!    gammaxstore,gammaystore,gammazstore,jacobianstore
!
!! arrays with mesh parameters
!  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
!  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
!  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)
!
!  real(kind=CUSTOM_REAL) kappastore(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) mustore(NGLLX,NGLLY,NGLLZ,nspec)
!
!  real(kind=CUSTOM_REAL) c11store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c12store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c13store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c14store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c15store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c16store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c22store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c23store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c24store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c25store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c26store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c33store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c34store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c35store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c36store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c44store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c45store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c46store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c55store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c56store(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) c66store(NGLLX,NGLLY,NGLLZ,nspec)
!
!! Stacey
!  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec)
!  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,nspec)
!
!! flag indicating whether point is in the sediments
!  logical flag_sediments(NGLLX,NGLLY,NGLLZ,nspec)
!  logical not_fully_in_bedrock(nspec)
!
!  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)
!
!! doubling mesh flag
!  integer idoubling(nspec)
!
!! mass matrix
!  real(kind=CUSTOM_REAL) rmass(nglob)
!
!! additional ocean load mass matrix
!  real(kind=CUSTOM_REAL) rmass_ocean_load(npointot_oceans)
!
!! boundary parameters locator
!  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
!  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
!  integer ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP)
!
!! normals
!  real(kind=CUSTOM_REAL) normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
!  real(kind=CUSTOM_REAL) normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
!  real(kind=CUSTOM_REAL) normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
!  real(kind=CUSTOM_REAL) normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
!  real(kind=CUSTOM_REAL) normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM)
!  real(kind=CUSTOM_REAL) normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP)
!
!! jacobian on 2D edges
!  real(kind=CUSTOM_REAL) jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
!  real(kind=CUSTOM_REAL) jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
!  real(kind=CUSTOM_REAL) jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
!  real(kind=CUSTOM_REAL) jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
!  real(kind=CUSTOM_REAL) jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM)
!  real(kind=CUSTOM_REAL) jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP)
!
!! number of elements on the boundaries
!  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
!
!! MPI cut-planes parameters along xi and along eta
!  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)
!
!  integer i,j,k,ispec,iglob
!
!! processor identification
!  character(len=256) prname
!
!! xix
!  open(unit=27,file=prname(1:len_trim(prname))//'xix.bin',status='unknown',form='unformatted')
!  write(27) xixstore
!  close(27)
!
!! xiy
!  open(unit=27,file=prname(1:len_trim(prname))//'xiy.bin',status='unknown',form='unformatted')
!  write(27) xiystore
!  close(27)
!
!! xiz
!  open(unit=27,file=prname(1:len_trim(prname))//'xiz.bin',status='unknown',form='unformatted')
!  write(27) xizstore
!  close(27)
!
!! etax
!  open(unit=27,file=prname(1:len_trim(prname))//'etax.bin',status='unknown',form='unformatted')
!  write(27) etaxstore
!  close(27)
!
!! etay
!  open(unit=27,file=prname(1:len_trim(prname))//'etay.bin',status='unknown',form='unformatted')
!  write(27) etaystore
!  close(27)
!
!! etaz
!  open(unit=27,file=prname(1:len_trim(prname))//'etaz.bin',status='unknown',form='unformatted')
!  write(27) etazstore
!  close(27)
!
!! gammax
!  open(unit=27,file=prname(1:len_trim(prname))//'gammax.bin',status='unknown',form='unformatted')
!  write(27) gammaxstore
!  close(27)
!
!! gammay
!  open(unit=27,file=prname(1:len_trim(prname))//'gammay.bin',status='unknown',form='unformatted')
!  write(27) gammaystore
!  close(27)
!
!! gammaz
!  open(unit=27,file=prname(1:len_trim(prname))//'gammaz.bin',status='unknown',form='unformatted')
!  write(27) gammazstore
!  close(27)
!
!! jacobian
!  open(unit=27,file=prname(1:len_trim(prname))//'jacobian.bin',status='unknown',form='unformatted')
!  write(27) jacobianstore
!  close(27)
!
!! flag_sediments
!  open(unit=27,file=prname(1:len_trim(prname))//'flag_sediments.bin',status='unknown',form='unformatted')
!  write(27) flag_sediments
!  close(27)
!
!! not_fully_in_bedrock
!  open(unit=27,file=prname(1:len_trim(prname))//'not_fully_in_bedrock.bin',status='unknown',form='unformatted')
!  write(27) not_fully_in_bedrock
!  close(27)
!
!! rho_vs
!! Stacey
!! rho_vp
!  open(unit=27,file=prname(1:len_trim(prname))//'rho_vp.bin',status='unknown',form='unformatted')
!  write(27) rho_vp
!  close(27)
!
!! rho_vs
!  open(unit=27,file=prname(1:len_trim(prname))//'rho_vs.bin',status='unknown',form='unformatted')
!  write(27) rho_vs
!  close(27)
!
!!!$! vp (for checking the mesh and model)
!!!$  open(unit=27,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted')
!!!$  write(27) (FOUR_THIRDS * mustore + kappastore) / rho_vp
!!!$  close(27)
!!!$
!!!$! vs (for checking the mesh and model)
!!!$  open(unit=27,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted')
!!!$  write(27) mustore / rho_vs
!!!$  close(27)
!
!! kappa
!  open(unit=27,file=prname(1:len_trim(prname))//'kappa.bin',status='unknown',form='unformatted')
!  write(27) kappastore
!  close(27)
!
!! mu
!  open(unit=27,file=prname(1:len_trim(prname))//'mu.bin',status='unknown',form='unformatted')
!  write(27) mustore
!  close(27)
!
!! ibool
!  open(unit=27,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted')
!  write(27) ibool
!  close(27)
!
!! doubling
!  open(unit=27,file=prname(1:len_trim(prname))//'idoubling.bin',status='unknown',form='unformatted')
!  write(27) idoubling
!  close(27)
!
!! mass matrix
!  open(unit=27,file=prname(1:len_trim(prname))//'rmass.bin',status='unknown',form='unformatted')
!  write(27) rmass
!  close(27)
!
!! For anisotropy
!  if(ANISOTROPY) then
!     ! c11
!     open(unit=27,file=prname(1:len_trim(prname))//'c11.bin',status='unknown',form='unformatted')
!     write(27) c11store
!     close(27)
!
!     ! c12
!     open(unit=27,file=prname(1:len_trim(prname))//'c12.bin',status='unknown',form='unformatted')
!     write(27) c12store
!     close(27)
!
!     ! c13
!     open(unit=27,file=prname(1:len_trim(prname))//'c13.bin',status='unknown',form='unformatted')
!     write(27) c13store
!     close(27)
!
!     ! c14
!     open(unit=27,file=prname(1:len_trim(prname))//'c14.bin',status='unknown',form='unformatted')
!     write(27) c14store
!     close(27)
!
!     ! c15
!     open(unit=27,file=prname(1:len_trim(prname))//'c15.bin',status='unknown',form='unformatted')
!     write(27) c15store
!     close(27)
!
!     ! c16
!     open(unit=27,file=prname(1:len_trim(prname))//'c16.bin',status='unknown',form='unformatted')
!     write(27) c16store
!     close(27)
!
!     ! c22
!     open(unit=27,file=prname(1:len_trim(prname))//'c22.bin',status='unknown',form='unformatted')
!     write(27) c22store
!     close(27)
!
!     ! c23
!     open(unit=27,file=prname(1:len_trim(prname))//'c23.bin',status='unknown',form='unformatted')
!     write(27) c23store
!     close(27)
!
!     ! c24
!     open(unit=27,file=prname(1:len_trim(prname))//'c24.bin',status='unknown',form='unformatted')
!     write(27) c24store
!     close(27)
!
!     ! c25
!     open(unit=27,file=prname(1:len_trim(prname))//'c25.bin',status='unknown',form='unformatted')
!     write(27) c25store
!     close(27)
!
!     ! c26
!     open(unit=27,file=prname(1:len_trim(prname))//'c26.bin',status='unknown',form='unformatted')
!     write(27) c26store
!     close(27)
!
!     ! c33
!     open(unit=27,file=prname(1:len_trim(prname))//'c33.bin',status='unknown',form='unformatted')
!     write(27) c33store
!     close(27)
!
!     ! c34
!     open(unit=27,file=prname(1:len_trim(prname))//'c34.bin',status='unknown',form='unformatted')
!     write(27) c34store
!     close(27)
!
!     ! c35
!     open(unit=27,file=prname(1:len_trim(prname))//'c35.bin',status='unknown',form='unformatted')
!     write(27) c35store
!     close(27)
!
!     ! c36
!     open(unit=27,file=prname(1:len_trim(prname))//'c36.bin',status='unknown',form='unformatted')
!     write(27) c36store
!     close(27)
!
!     ! c44
!     open(unit=27,file=prname(1:len_trim(prname))//'c44.bin',status='unknown',form='unformatted')
!     write(27) c44store
!     close(27)
!
!     ! c45
!     open(unit=27,file=prname(1:len_trim(prname))//'c45.bin',status='unknown',form='unformatted')
!     write(27) c45store
!     close(27)
!
!     ! c46
!     open(unit=27,file=prname(1:len_trim(prname))//'c46.bin',status='unknown',form='unformatted')
!     write(27) c46store
!     close(27)
!
!     ! c55
!     open(unit=27,file=prname(1:len_trim(prname))//'c55.bin',status='unknown',form='unformatted')
!     write(27) c55store
!     close(27)
!
!     ! c56
!     open(unit=27,file=prname(1:len_trim(prname))//'c56.bin',status='unknown',form='unformatted')
!     write(27) c56store
!     close(27)
!
!     ! c66
!     open(unit=27,file=prname(1:len_trim(prname))//'c66.bin',status='unknown',form='unformatted')
!     write(27) c66store
!     close(27)
!
!  endif
!
!! additional ocean load mass matrix if oceans
!  if(OCEANS) then
!    open(unit=27,file=prname(1:len_trim(prname))//'rmass_ocean_load.bin',status='unknown',form='unformatted')
!    write(27) rmass_ocean_load
!    close(27)
!  endif
!
!! boundary parameters
!  open(unit=27,file=prname(1:len_trim(prname))//'ibelm.bin',status='unknown',form='unformatted')
!  write(27) ibelm_xmin
!  write(27) ibelm_xmax
!  write(27) ibelm_ymin
!  write(27) ibelm_ymax
!  write(27) ibelm_bottom
!  write(27) ibelm_top
!  close(27)
!
!  open(unit=27,file=prname(1:len_trim(prname))//'normal.bin',status='unknown',form='unformatted')
!  write(27) normal_xmin
!  write(27) normal_xmax
!  write(27) normal_ymin
!  write(27) normal_ymax
!  write(27) normal_bottom
!  write(27) normal_top
!  close(27)
!
!  open(unit=27,file=prname(1:len_trim(prname))//'jacobian2D.bin',status='unknown',form='unformatted')
!  write(27) jacobian2D_xmin
!  write(27) jacobian2D_xmax
!  write(27) jacobian2D_ymin
!  write(27) jacobian2D_ymax
!  write(27) jacobian2D_bottom
!  write(27) jacobian2D_top
!  close(27)
!
!  open(unit=27,file=prname(1:len_trim(prname))//'nspec2D.bin',status='unknown',form='unformatted')
!  write(27) nspec2D_xmin
!  write(27) nspec2D_xmax
!  write(27) nspec2D_ymin
!  write(27) nspec2D_ymax
!  close(27)
!
!! MPI cut-planes parameters along xi and along eta
!  open(unit=27,file=prname(1:len_trim(prname))//'iMPIcut_xi.bin',status='unknown',form='unformatted')
!  write(27) iMPIcut_xi
!  close(27)
!
!  open(unit=27,file=prname(1:len_trim(prname))//'iMPIcut_eta.bin',status='unknown',form='unformatted')
!  write(27) iMPIcut_eta
!  close(27)
!
!! mesh arrays used in the solver to locate source and receivers
!! use rmass for temporary storage to perform conversion, since already saved
!
!!--- x coordinate
!  rmass(:) = 0._CUSTOM_REAL
!  do ispec = 1,nspec
!    do k = 1,NGLLZ
!      do j = 1,NGLLY
!        do i = 1,NGLLX
!          iglob = ibool(i,j,k,ispec)
!! distinguish between single and double precision for reals
!          if(CUSTOM_REAL == SIZE_REAL) then
!            rmass(iglob) = sngl(xstore(i,j,k,ispec))
!          else
!            rmass(iglob) = xstore(i,j,k,ispec)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!  open(unit=27,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted')
!  write(27) rmass
!  close(27)
!
!!--- y coordinate
!  rmass(:) = 0._CUSTOM_REAL
!  do ispec = 1,nspec
!    do k = 1,NGLLZ
!      do j = 1,NGLLY
!        do i = 1,NGLLX
!          iglob = ibool(i,j,k,ispec)
!! distinguish between single and double precision for reals
!          if(CUSTOM_REAL == SIZE_REAL) then
!            rmass(iglob) = sngl(ystore(i,j,k,ispec))
!          else
!            rmass(iglob) = ystore(i,j,k,ispec)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!  open(unit=27,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted')
!  write(27) rmass
!  close(27)
!
!!--- z coordinate
!  rmass(:) = 0._CUSTOM_REAL
!  do ispec = 1,nspec
!    do k = 1,NGLLZ
!      do j = 1,NGLLY
!        do i = 1,NGLLX
!          iglob = ibool(i,j,k,ispec)
!! distinguish between single and double precision for reals
!          if(CUSTOM_REAL == SIZE_REAL) then
!            rmass(iglob) = sngl(zstore(i,j,k,ispec))
!          else
!            rmass(iglob) = zstore(i,j,k,ispec)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!  open(unit=27,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted')
!  write(27) rmass
!  close(27)
!
!  end subroutine save_arrays_solver
!
!!=============================================================
    
  