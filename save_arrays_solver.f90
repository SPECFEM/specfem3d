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
                    kappastore,mustore,rmass,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                    NSPEC2D_TOP,ibelm_top,normal_top,jacobian2D_top, &
                    absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
                    absorbing_boundary_ijk,absorbing_boundary_ispec, &
                    num_absorbing_boundary_faces, &
                    num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    prname,SAVE_MESH_FILES, &
                    ANISOTROPY,NSPEC_ANISO, &
                    c11store,c12store,c13store,c14store,c15store,c16store, &
                    c22store,c23store,c24store,c25store,c26store,c33store, &
                    c34store,c35store,c36store,c44store,c45store,c46store, &
                    c55store,c56store,c66store)


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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass

! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  
! absorbing boundaries  
!  integer  :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM
  integer :: NSPEC2D_TOP
!  integer, dimension(nspec2D_xmin)  :: ibelm_xmin  
!  integer, dimension(nspec2D_xmax)  :: ibelm_xmax
!  integer, dimension(nspec2D_ymin)  :: ibelm_ymin
!  integer, dimension(nspec2D_ymax)  :: ibelm_ymax
!  integer, dimension(NSPEC2D_BOTTOM)  :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP)  :: ibelm_top
!  integer :: ibelm_gll_xmin(3,NGLLY,NGLLZ,nspec2D_xmin),ibelm_gll_xmax(3,NGLLY,NGLLZ,nspec2D_xmax), &
!            ibelm_gll_ymin(3,NGLLX,NGLLZ,nspec2D_ymin),ibelm_gll_ymax(3,NGLLX,NGLLZ,nspec2D_ymax), &
!            ibelm_gll_bottom(3,NGLLY,NGLLY,nspec2D_bottom),ibelm_gll_top(3,NGLLY,NGLLY,nspec2D_top)  
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_xmin) :: normal_xmin
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_xmax) :: normal_xmax
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,nspec2D_ymin) :: normal_ymin
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,nspec2D_ymax) :: normal_ymax
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM) :: normal_bottom
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP) :: normal_top  
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nspec2D_xmin) :: jacobian2D_xmin
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nspec2D_xmax) :: jacobian2D_xmax
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec2D_ymin) :: jacobian2D_ymin
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec2D_ymax) :: jacobian2D_ymax
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,nspec2D_bottom) :: jacobian2D_bottom
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,nspec2D_top) :: jacobian2D_top
  
  integer :: num_absorbing_boundary_faces
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_normal
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_jacobian2D
  integer, dimension(3,NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_ijk
  integer, dimension(num_absorbing_boundary_faces) :: absorbing_boundary_ispec
  

!  integer  :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
!  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nimin,nimax,nkmin_eta
!  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX) :: njmin,njmax,nkmin_xi

! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

! file name
  character(len=150) prname
  logical :: SAVE_MESH_FILES

! anisotropy
  logical :: ANISOTROPY
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store
  
! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: v_tmp
  real(kind=CUSTOM_REAL) :: minimum(1)
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: ier,i  

! saves mesh file proc***_external_mesh.bin
  open(unit=IOUT,file=prname(1:len_trim(prname))//'external_mesh.bin',status='unknown',action='write',form='unformatted',iostat=ier)
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

  !pll Stacey 
  write(IOUT) rho_vp
  write(IOUT) rho_vs
  write(IOUT) iflag_attenuation_store
  
!  write(IOUT) NSPEC2DMAX_XMIN_XMAX 
!  write(IOUT) NSPEC2DMAX_YMIN_YMAX
!  write(IOUT) nimin
!  write(IOUT) nimax
!  write(IOUT) njmin
!  write(IOUT) njmax
!  write(IOUT) nkmin_xi 
!  write(IOUT) nkmin_eta
  !end pll

  write(IOUT) kappastore
  write(IOUT) mustore

  write(IOUT) rmass
  write(IOUT) ibool

  write(IOUT) xstore_dummy
  write(IOUT) ystore_dummy
  write(IOUT) zstore_dummy

! absorbing boundary parameters
!  write(IOUT) nspec2D_xmin
!  write(IOUT) nspec2D_xmax
!  write(IOUT) nspec2D_ymin
!  write(IOUT) nspec2D_ymax
!  write(IOUT) NSPEC2D_BOTTOM
!  write(IOUT) NSPEC2D_TOP
!
!  write(IOUT) ibelm_xmin
!  write(IOUT) ibelm_xmax
!  write(IOUT) ibelm_ymin
!  write(IOUT) ibelm_ymax
!  write(IOUT) ibelm_bottom
!  write(IOUT) ibelm_top
!
!  write(IOUT) ibelm_gll_xmin
!  write(IOUT) ibelm_gll_xmax
!  write(IOUT) ibelm_gll_ymin
!  write(IOUT) ibelm_gll_ymax
!  write(IOUT) ibelm_gll_bottom
!  write(IOUT) ibelm_gll_top
!
!  write(IOUT) normal_xmin
!  write(IOUT) normal_xmax
!  write(IOUT) normal_ymin
!  write(IOUT) normal_ymax
!  write(IOUT) normal_bottom
!  write(IOUT) normal_top
!
!  write(IOUT) jacobian2D_xmin
!  write(IOUT) jacobian2D_xmax
!  write(IOUT) jacobian2D_ymin
!  write(IOUT) jacobian2D_ymax
!  write(IOUT) jacobian2D_bottom
!  write(IOUT) jacobian2D_top

  write(IOUT) num_absorbing_boundary_faces
  write(IOUT) absorbing_boundary_ispec
  write(IOUT) absorbing_boundary_ijk
  write(IOUT) absorbing_boundary_jacobian2D
  write(IOUT) absorbing_boundary_normal

! free surface 
  write(IOUT) NSPEC2D_TOP    
  write(IOUT) ibelm_top
  write(IOUT) jacobian2D_top
  write(IOUT) normal_top

!MPI interfaces
  write(IOUT) num_interfaces_ext_mesh
  write(IOUT) maxval(nibool_interfaces_ext_mesh)
  write(IOUT) my_neighbours_ext_mesh
  write(IOUT) nibool_interfaces_ext_mesh

  allocate(ibool_interfaces_ext_mesh_dummy(maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'
  
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy = ibool_interfaces_ext_mesh(1:maxval(nibool_interfaces_ext_mesh),:)
  enddo
  write(IOUT) ibool_interfaces_ext_mesh_dummy

  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier); if( ier /= 0 ) stop 'error deallocating array'

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



! mesh arrays used for example in combine_vol_data.f90
  if( SAVE_MESH_FILES ) then
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
    minimum = minval( abs(rho_vp) )
    if( minimum(1) /= 0.0 ) then
      v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
    else
      v_tmp = 0.0
    endif  
    open(unit=27,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

! vs (for checking the mesh and model)
    minimum = minval( abs(rho_vs) )
    if( minimum(1) /= 0.0 ) then
      v_tmp = mustore / rho_vs
    else  
      v_tmp = 0.0
    endif
    open(unit=27,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted')
    write(27) v_tmp
    close(27)

    deallocate(v_tmp,stat=ier); if( ier /= 0 ) stop 'error deallocating array'
  endif

  end subroutine save_arrays_solver_ext_mesh
  
  
  
!=============================================================

! external mesh routine for saving vtk file holding integer flag for each element

  subroutine save_arrays_solver_ext_mesh_elem_vtk(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            elem_flag,prname_file)


  implicit none

  include "constants.h"

  integer :: nspec,nglob

! global coordinates  
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! element flag array
  integer, dimension(nspec) :: elem_flag  
  integer :: ispec,i

! file name
  character(len=150) prname_file

! write source and receiver VTK files for Paraview
  write(IMAIN,*) '  vtk file: ',prname_file(1:len_trim(prname_file))//'.vtk'
  
  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3f)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i,i)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""
  
  ! type: hexahedrons
  write(IOVTK,'(a,i)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ""
  
  write(IOVTK,'(a,i)') "CELL_DATA ",nspec
  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOVTK,*) elem_flag(ispec)
  enddo
  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine save_arrays_solver_ext_mesh_elem_vtk
  
  
!=============================================================

! external mesh routine for saving vtk files for values on all gll points

  subroutine save_arrays_solver_ext_mesh_glldata_vtk(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            gll_data,prname_file)

  implicit none

  include "constants.h"

  integer :: nspec,nglob
  
! global coordinates  
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array  
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! masking arrays (takes first data value assigned on a global point, ignores any data values later on for the same global point)
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool
  
! file name
  character(len=150) prname_file

  integer :: ispec,i,j,k,ier,iglob

! write source and receiver VTK files for Paraview
  write(IMAIN,*) '  vtk file: ',prname_file(1:len_trim(prname_file))//'.vtk'
  
  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3f)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i,i)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""
  
  ! type: hexahedrons
  write(IOVTK,'(a,i)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ""
    
  ! iflag field on global nodeset
  allocate(mask_ibool(nglob),flag_val(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask'
  
  mask_ibool = .false.
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if( .not. mask_ibool(iglob) ) then   
            flag_val(iglob) = gll_data(i,j,k,ispec)
            mask_ibool(iglob) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

  write(IOVTK,'(a,i)') "POINT_DATA ",nglob
  write(IOVTK,'(a)') "SCALARS gll_data float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob    
      write(IOVTK,*) flag_val(i)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)


  end subroutine save_arrays_solver_ext_mesh_glldata_vtk

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
!  character(len=150) prname
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
    
  