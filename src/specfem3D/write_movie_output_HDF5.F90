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


! HDF5 format movie output

#ifdef USE_HDF5

  module specfem_par_movie_hdf5

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,MAX_STRING_LEN,OUTPUT_FILES,myrank

  use shared_parameters, only: NPROC,NSTEP,NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES, &
    ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION

  !#TODO: hdf5 i/o server
  !use shared_parameters, only: HDF5_IO_NNODES

  use specfem_par, only: NSPEC_AB,NGLOB_AB,ibool,it,DT,t0,xstore,ystore,zstore

  use specfem_par_movie, only: faces_surface_offset, &
    store_val_ux,store_val_uy,store_val_uz, &
    store_val_x_all,store_val_y_all,store_val_z_all, &
    shakemap_ux,shakemap_uy,shakemap_uz, &
    wavefield_x,wavefield_y,wavefield_z,wavefield_pressure, &
    div_glob,div,curl_x,curl_y,curl_z

  use manager_hdf5

  implicit none

  ! xdmf
  ! file i/o units
  integer :: xdmf_shake               = 191
  integer :: xdmf_surf                = 192
  integer :: xdmf_vol                 = 193
  ! file position for entering surface snapshot records
  integer :: surf_xdmf_pos = 0

  ! surface movie arrays
  integer :: size_surf_array = 0
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_x,   surf_y,   surf_z, &
                                                       surf_ux,  surf_uy,  surf_uz
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_x_aug,   surf_y_aug,   surf_z_aug, &
                                                       surf_ux_aug,  surf_uy_aug,  surf_uz_aug
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shake_ux, shake_uy, shake_uz, &
                                                       shake_ux_aug, shake_uy_aug, shake_uz_aug

  ! volume movie arrays
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: div_on_node, curl_x_on_node, curl_y_on_node, curl_z_on_node
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: velocity_x_on_node,velocity_y_on_node,velocity_z_on_node

  ! arrays for hdf5 output
  integer, dimension(:), allocatable :: nelm_par_proc_nio, nglob_par_proc_nio
  integer, dimension(:), allocatable :: nglob_offset
  integer, dimension(:), allocatable :: nelm_offset

contains

  !-------------------------------------------

  function r2c(k) result(str)

  ! "Convert an real to string."

  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: k
  character(len=20) str
  write (str, *) k
  str = adjustl(str)
  end function r2c

  !-------------------------------------------

  ! reorder and expand the input array (for only corner nodes) to high res array (all GLL)
  subroutine recompose_for_hires(arr_in, arr_out)

  implicit none

  real(kind=CUSTOM_REAL), dimension(:), intent(in)    :: arr_in
  real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: arr_out

  integer :: nfaces_actual
  integer :: i,j,k,c
  integer :: factor_face_aug = (NGLLX-1)*(NGLLY-1)
  integer :: npoint_per_face = NGLLX*NGLLY
  integer :: npoint_corner = 4

  nfaces_actual = size(arr_in)/(NGLLX*NGLLY) ! expecting NGLLX == NGLLY == NGLLZ

  c = 1
  do i = 0, nfaces_actual-1
    do j = 0,NGLLY-2 !y
      do k = 0,NGLLX-2 !x
        arr_out(c  +j*(NGLLX-1)*4+k*4) = arr_in(i*npoint_per_face+1      +k+j*NGLLX)
        arr_out(c+1+j*(NGLLX-1)*4+k*4) = arr_in(i*npoint_per_face+2      +k+j*NGLLX)
        arr_out(c+2+j*(NGLLX-1)*4+k*4) = arr_in(i*npoint_per_face+2+NGLLX+k+j*NGLLX)
        arr_out(c+3+j*(NGLLX-1)*4+k*4) = arr_in(i*npoint_per_face+1+NGLLX+k+j*NGLLX)
      enddo
    enddo
    c = c+factor_face_aug*npoint_corner
  enddo

  end subroutine recompose_for_hires

  !-------------------------------------------

  subroutine elm2node_base(elm_base, node_base)

  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: elm_base
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(inout) :: node_base
  ! local parameters
  integer :: i,j,k,iglob,ispec

  ! converts from local to global indexing for input array
  do ispec = 1,NSPEC_AB
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! global field
          node_base(iglob) = elm_base(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  end subroutine elm2node_base

  !-------------------------------------------

  subroutine get_conn_for_movie(elm_conn,offset)

  implicit none

  integer, dimension(9,NSPEC_AB*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)), intent(out) :: elm_conn
  integer, intent(in) :: offset ! node id offset (starting global element id of each proc)

  ! local parameters
  integer :: ispec,ii,icub,jcub,kcub
  integer,parameter :: cell_type = 9

  do ispec = 1,NSPEC_AB
    ! extract information from full GLL grid
    ! node order follows vtk format
    do icub = 0,NGLLX-2
      do jcub = 0,NGLLY-2
        do kcub = 0,NGLLZ-2
          ii = 1+(ispec-1)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1) + (icub*(NGLLY-1)*(NGLLZ-1)+jcub*(NGLLZ-1)+kcub)
          elm_conn(1, ii)  = cell_type
          elm_conn(2, ii)  = ibool(icub+1,jcub+1,kcub+1,ispec)-1 + offset   ! node id starts 0 in xdmf rule
          elm_conn(3, ii)  = ibool(icub+2,jcub+1,kcub+1,ispec)-1 + offset
          elm_conn(4, ii)  = ibool(icub+2,jcub+2,kcub+1,ispec)-1 + offset
          elm_conn(5, ii)  = ibool(icub+1,jcub+2,kcub+1,ispec)-1 + offset
          elm_conn(6, ii)  = ibool(icub+1,jcub+1,kcub+2,ispec)-1 + offset
          elm_conn(7, ii)  = ibool(icub+2,jcub+1,kcub+2,ispec)-1 + offset
          elm_conn(8, ii)  = ibool(icub+2,jcub+2,kcub+2,ispec)-1 + offset
          elm_conn(9, ii)  = ibool(icub+1,jcub+2,kcub+2,ispec)-1 + offset
        enddo
      enddo
    enddo
  enddo

  end subroutine get_conn_for_movie


  end module specfem_par_movie_hdf5
#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wmo_save_shakemap_hdf5()

#ifdef USE_HDF5

  use specfem_par_movie_hdf5

  implicit none

  ! local parameters
  integer :: ier, nfaces_actual, nfaces_aug=(NGLLX-1)*(NGLLY-1),nnodes_per_face_aug=4
  integer :: len_array_aug, len_array_aug_proc
  character(len=64) :: dset_name, group_name
  real(kind=CUSTOM_REAL) :: aug_factor
  character(len=MAX_STRING_LEN) :: fname_h5_data_shake

  logical, parameter :: if_corrective = .true.

  !#TODO: hdf5 i/o server
  !if (HDF5_IO_NNODES > 0) then
  !  ! send shakemap arrays to io node
  !  call isend_cr_inter(shakemap_ux,nfaces_surface_points,0,io_tag_shake_ux,req)
  !  call isend_cr_inter(shakemap_uy,nfaces_surface_points,0,io_tag_shake_uy,req)
  !  call isend_cr_inter(shakemap_uz,nfaces_surface_points,0,io_tag_shake_uz,req)
  !  ! all done
  !  return
  !endif

  fname_h5_data_shake = trim(OUTPUT_FILES) // "/shakemap.h5"

  ! get the offset info from main rank
  call bcast_all_i(faces_surface_offset,size(faces_surface_offset))

  if (myrank == 0) then
    size_surf_array = size(store_val_x_all)
    call bcast_all_singlei(size_surf_array)
  else
    call bcast_all_singlei(size_surf_array)
  endif

  nfaces_actual = size_surf_array / (NGLLX*NGLLY)
  len_array_aug = nfaces_actual * nfaces_aug * nnodes_per_face_aug

  ! initialization of h5 file
  call h5_initialize()

  if (myrank == 0) then
    ! create a hdf5 file
    call h5_create_file(fname_h5_data_shake)

    ! save xyz coordinate array
    group_name = "surf_coord"
    call h5_create_group(group_name)
    call h5_open_group(group_name)

    ! low resolution output
    if (.not. USE_HIGHRES_FOR_MOVIES) then
      dset_name = "x"
      call h5_write_dataset(dset_name, store_val_x_all)
      call h5_close_dataset()
      dset_name = "y"
      call h5_write_dataset(dset_name, store_val_y_all)
      call h5_close_dataset()
      dset_name = "z"
      call h5_write_dataset(dset_name, store_val_z_all)
      call h5_close_dataset()

      ! write xdmf header
      call write_xdmf_shakemap(size_surf_array)

      ! high resolution output
    else
      allocate(surf_x_aug(len_array_aug), &
               surf_y_aug(len_array_aug), &
               surf_z_aug(len_array_aug),stat=ier)
      if (ier /= 0) stop 'Error allocating surf_x_aug arrays'

      ! nfaces*25nodes => n*16faces*4
      dset_name = "x"
      call recompose_for_hires(store_val_x_all,surf_x_aug)
      call h5_write_dataset(dset_name, surf_x_aug)
      call h5_close_dataset()
      dset_name = "y"
      call recompose_for_hires(store_val_y_all,surf_y_aug)
      call h5_write_dataset(dset_name, surf_y_aug)
      call h5_close_dataset()
      dset_name = "z"
      call recompose_for_hires(store_val_z_all,surf_z_aug)
      call h5_write_dataset(dset_name, surf_z_aug)
      call h5_close_dataset()

      ! write xdmf header
      call write_xdmf_shakemap(len_array_aug)

      deallocate(surf_x_aug, surf_y_aug, surf_z_aug)
    endif

    call h5_close_group()
    call h5_close_file()

  endif

  call synchronize_all()

  ! write dataset in one single h5 file in collective mode.
  group_name = "shakemap"

  if (myrank == 0) then
    call h5_open_file(fname_h5_data_shake)
    call h5_create_group(group_name)
    call h5_open_group(group_name)
    if (.not. USE_HIGHRES_FOR_MOVIES) then
      call h5_create_dataset_gen_in_group("shakemap_ux",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("shakemap_uy",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("shakemap_uz",(/size_surf_array/),1,CUSTOM_REAL)
    else
      call h5_create_dataset_gen_in_group("shakemap_ux",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("shakemap_uy",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("shakemap_uz",(/len_array_aug/),1,CUSTOM_REAL)
   endif
    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  call h5_open_file_p(fname_h5_data_shake)
  call h5_open_group(group_name)

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "shakemap_ux"
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, shakemap_ux, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "shakemap_uy"
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, shakemap_uy, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "shakemap_uz"
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, shakemap_uz, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)

  else
    nfaces_actual      = size(shakemap_ux)/(NGLLX*NGLLY)
    len_array_aug_proc = nfaces_actual*nfaces_aug*nnodes_per_face_aug ! augmented array length for each proc
    aug_factor         = real(len_array_aug)/real(size_surf_array)

    allocate(surf_ux_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uy_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uz_aug(len_array_aug_proc),stat=ier)

    dset_name = "shakemap_ux"
    call recompose_for_hires(shakemap_ux, surf_ux_aug)
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, surf_ux_aug, &
                      (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "shakemap_uy"
    call recompose_for_hires(shakemap_uy, surf_uy_aug)
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, surf_uy_aug, &
                      (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "shakemap_uz"
    call recompose_for_hires(shakemap_uz, surf_uz_aug)
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, surf_uz_aug, &
                      (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)

    deallocate(surf_ux_aug, surf_uy_aug, surf_uz_aug)
  endif

  call h5_close_group()
  call h5_close_file()
  call h5_finalize()

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 routine wmo_save_shakemap_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'wmo_save_shakemap_hdf5() called without compilation support'

#endif

  end subroutine wmo_save_shakemap_hdf5

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine write_xdmf_shakemap(num_nodes)

  use specfem_par_movie_hdf5

  implicit none
  integer :: num_elm, num_nodes
  character(len=MAX_STRING_LEN) :: fname_xdmf_shake
  character(len=MAX_STRING_LEN) :: fname_h5_data_shake_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  num_elm = num_nodes / 4

  ! writeout xdmf file for surface movie
  fname_xdmf_shake = trim(OUTPUT_FILES) // "/shakemap.xmf"
  fname_h5_data_shake_xdmf = "./shakemap.h5"       ! relative to shakemap.xmf file
  ! note: this seems to point to / < rootdir>/OUTPUT_FILES/./OUTPUT_FILES//shakemap.h5
  !       and paraview won't find it...
  !         fname_h5_data_shake_xdmf = trim(OUTPUT_FILES) // "/shakemap.h5"

  open(unit=xdmf_shake, file=trim(fname_xdmf_shake), recl=256)

  write(xdmf_shake,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_shake,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_shake,*) '<Xdmf Version="3.0">'
  write(xdmf_shake,*) '<Domain Name="shakemap">'
  write(xdmf_shake,*) '<Grid>'
  write(xdmf_shake,*) '<Topology Name="topo" TopologyType="Quadrilateral" NumberOfElements="'//trim(i2c(num_elm))//'"/>'
  write(xdmf_shake,*) '<Geometry GeometryType="X_Y_Z">'
  write(xdmf_shake,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                                                    //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        '//trim(fname_h5_data_shake_xdmf)//':/surf_coord/x'

  write(xdmf_shake,*) '</DataItem>'
  write(xdmf_shake,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                    //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        '//trim(fname_h5_data_shake_xdmf)//':/surf_coord/y'
  write(xdmf_shake,*) '</DataItem>'
  write(xdmf_shake,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                   //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        '//trim(fname_h5_data_shake_xdmf)//':/surf_coord/z'
  write(xdmf_shake,*) '</DataItem>'
  write(xdmf_shake,*) '</Geometry>'
  write(xdmf_shake,*) '<Attribute Name="shake_ux" AttributeType="Scalar" Center="Node">'
  write(xdmf_shake,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                  //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        '//trim(fname_h5_data_shake_xdmf)//':/shakemap/shakemap_ux'
  write(xdmf_shake,*) '</DataItem>'
  write(xdmf_shake,*) '</Attribute>'
  write(xdmf_shake,*) '<Attribute Name="shake_uy" AttributeType="Scalar" Center="Node">'
  write(xdmf_shake,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                 //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        '//trim(fname_h5_data_shake_xdmf)//':/shakemap/shakemap_uy'
  write(xdmf_shake,*) '</DataItem>'
  write(xdmf_shake,*) '</Attribute>'
  write(xdmf_shake,*) '<Attribute Name="shake_uz" AttributeType="Scalar" Center="Node">'
  write(xdmf_shake,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        '//trim(fname_h5_data_shake_xdmf)//':/shakemap/shakemap_uz'
  write(xdmf_shake,*) '</DataItem>'
  write(xdmf_shake,*) '</Attribute>'
  write(xdmf_shake,*) '</Grid>'
  write(xdmf_shake,*) '</Domain>'
  write(xdmf_shake,*) '</Xdmf>'

  close(xdmf_shake)

  end subroutine write_xdmf_shakemap

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wmo_movie_save_surface_snapshot_hdf5()

#ifdef USE_HDF5

  use specfem_par_movie_hdf5

  implicit none

  integer :: ier, nfaces_actual
  integer :: len_array_aug, len_array_aug_proc
  character(len=64) :: dset_name, group_name
  character(len=12) :: tempstr
  real(kind=CUSTOM_REAL) :: aug_factor
  character(len=MAX_STRING_LEN) :: fname_h5_data_surf

  integer, parameter :: nfaces_aug = (NGLLX-1)*(NGLLY-1)
  integer, parameter :: nnodes_per_face_aug = 4
  logical, parameter :: if_corrective = .true.

  !#TODO: hdf5 i/o server
  !if (HDF5_IO_NNODES > 0) then
  !  ! send surface body to io node
  !  call isend_cr_inter(store_val_ux,nfaces_surface_points,0,io_tag_surface_ux,req_dump_surf(1))
  !  call isend_cr_inter(store_val_uy,nfaces_surface_points,0,io_tag_surface_uy,req_dump_surf(2))
  !  call isend_cr_inter(store_val_uz,nfaces_surface_points,0,io_tag_surface_uz,req_dump_surf(3))
  !  n_req_surf = 3
  !  ! all done
  !  return
  !endif

  ! initialize h5 and xdmf  if it == 0
  fname_h5_data_surf = trim(OUTPUT_FILES) // "/movie_surface.h5"

  ! get the offset info from main rank
  if (it == NTSTEP_BETWEEN_FRAMES) then
    call bcast_all_i(faces_surface_offset,size(faces_surface_offset))
  endif

  if (myrank == 0) then
    size_surf_array = size(store_val_x_all)
    call bcast_all_singlei(size_surf_array)
  else
    call bcast_all_singlei(size_surf_array)
  endif

  nfaces_actual = size_surf_array / (NGLLX*NGLLY)
  len_array_aug = nfaces_actual * nfaces_aug * nnodes_per_face_aug

  ! initialization of h5 file
  call h5_initialize()

  ! main process creates file on first occurrence
  if (it == NTSTEP_BETWEEN_FRAMES .and. myrank == 0) then
    ! create a hdf5 file
    call h5_create_file(fname_h5_data_surf)

    ! save xyz coordinate array
    group_name = "surf_coord"
    call h5_create_group(group_name)
    call h5_open_group(group_name)

    ! grid coordinates
    if (.not. USE_HIGHRES_FOR_MOVIES) then
      ! low resolution output
      dset_name = "x"
      call h5_write_dataset(dset_name, store_val_x_all)
      call h5_close_dataset()
      dset_name = "y"
      call h5_write_dataset(dset_name, store_val_y_all)
      call h5_close_dataset()
      dset_name = "z"
      call h5_write_dataset(dset_name, store_val_z_all)
      call h5_close_dataset()

      ! write xdmf header
      call write_xdmf_surface_header(size_surf_array)
    else
      ! high resolution output
      allocate(surf_x_aug(len_array_aug),stat=ier)
      allocate(surf_y_aug(len_array_aug),stat=ier)
      allocate(surf_z_aug(len_array_aug),stat=ier)

      ! nfaces*25nodes => n*16faces*4
      dset_name = "x"
      call recompose_for_hires(store_val_x_all,surf_x_aug)
      call h5_write_dataset(dset_name, surf_x_aug)
      call h5_close_dataset()
      dset_name = "y"
      call recompose_for_hires(store_val_y_all,surf_y_aug)
      call h5_write_dataset(dset_name, surf_y_aug)
      call h5_close_dataset()
      dset_name = "z"
      call recompose_for_hires(store_val_z_all,surf_z_aug)
      call h5_write_dataset(dset_name, surf_z_aug)
      call h5_close_dataset()

      ! write xdmf header
      call write_xdmf_surface_header(len_array_aug)

      deallocate(surf_x_aug, surf_y_aug, surf_z_aug)
    endif

    call h5_close_group()
    call h5_close_file()

  endif

  call synchronize_all()

  ! write dataset in one single h5 file in collective mode.
  ! create a group for each io step
  write(tempstr, "(i6.6)") it
  group_name = "it_" // trim(tempstr)

  if (myrank == 0) then
    call h5_open_file(fname_h5_data_surf)
    call h5_create_group(group_name)
    call h5_open_group(group_name)
    if (.not. USE_HIGHRES_FOR_MOVIES) then
      call h5_create_dataset_gen_in_group("ux",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uy",(/size_surf_array/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uz",(/size_surf_array/),1,CUSTOM_REAL)
    else
      call h5_create_dataset_gen_in_group("ux",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uy",(/len_array_aug/),1,CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uz",(/len_array_aug/),1,CUSTOM_REAL)
   endif
    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  call h5_open_file_p(fname_h5_data_surf)
  call h5_open_group(group_name)

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "ux"
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, store_val_ux, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "uy"
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, store_val_uy, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)
    dset_name = "uz"
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, store_val_uz, &
                      (/faces_surface_offset(myrank+1)/), if_corrective)

    ! write xdmf body
    if (myrank == 0) call write_xdmf_surface_body(it, size_surf_array)

  else
    nfaces_actual      = size(store_val_ux)/(NGLLX*NGLLY)
    len_array_aug_proc = nfaces_actual*nfaces_aug*nnodes_per_face_aug ! augmented array length for each proc
    aug_factor         = real(len_array_aug)/real(size_surf_array)

    allocate(surf_ux_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uy_aug(len_array_aug_proc),stat=ier)
    allocate(surf_uz_aug(len_array_aug_proc),stat=ier)

    dset_name = "ux"
    call recompose_for_hires(store_val_ux, surf_ux_aug)
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, surf_ux_aug, &
                      (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "uy"
    call recompose_for_hires(store_val_uy, surf_uy_aug)
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, surf_uy_aug, &
                      (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)
    dset_name = "uz"
    call recompose_for_hires(store_val_uz, surf_uz_aug)
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, surf_uz_aug, &
                      (/int(faces_surface_offset(myrank+1)*aug_factor)/), if_corrective)

    ! write xdmf body
    if (myrank == 0) call write_xdmf_surface_body(it, len_array_aug)

    deallocate(surf_ux_aug, surf_uy_aug, surf_uz_aug)
  endif

  call h5_close_group()
  call h5_close_file()
  call h5_finalize()
#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 routine wmo_movie_save_surface_snapshot_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'wmo_movie_save_surface_snapshot_hdf5() called without compilation support'

#endif

  end subroutine wmo_movie_save_surface_snapshot_hdf5

!
!-------------------------------------------------------------------------------------------------
!

!
! xdmf output routines
!
#ifdef USE_HDF5

  subroutine write_xdmf_surface_header(num_nodes)

  use specfem_par_movie_hdf5

  implicit none
  integer, intent(in) :: num_nodes
  ! local parameters
  integer :: num_elm
  character(len=MAX_STRING_LEN) :: fname_xdmf_surf
  character(len=MAX_STRING_LEN) :: fname_h5_data_surf_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  ! writeout xdmf file for surface movie
  fname_xdmf_surf = trim(OUTPUT_FILES) // "/movie_surface.xmf"
  fname_h5_data_surf_xdmf = "./movie_surface.h5"   ! relative to movie_surface.xmf file
  ! note: this seems not to work and point to a wrong directory:
  !         fname_h5_data_surf_xdmf = trim(OUTPUT_FILES) // "/movie_surface.h5"

  num_elm = num_nodes / 4

  open(unit=xdmf_surf, file=trim(fname_xdmf_surf), recl=256)

  write(xdmf_surf,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_surf,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_surf,*) '<Xdmf Version="3.0">'
  write(xdmf_surf,*) '<Domain Name="mesh">'
  write(xdmf_surf,*) '<Topology Name="topo" TopologyType="Quadrilateral" NumberOfElements="'//trim(i2c(num_elm))//'"/>'
  write(xdmf_surf,*) '<Geometry GeometryType="X_Y_Z">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                                                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/x'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/y'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                       //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/z'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Geometry>'

  write(xdmf_surf,*) '<Grid Name="fensap" GridType="Collection" CollectionType="Temporal">'
  ! 17 lines

  ! file finish
  write(xdmf_surf,*) '</Grid>'
  write(xdmf_surf,*) '</Domain>'
  write(xdmf_surf,*) '</Xdmf>'
  ! 20 lines

  ! position where the additional data will be inserted
  surf_xdmf_pos = 17

  close(xdmf_surf)

  end subroutine write_xdmf_surface_header

#endif

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine write_xdmf_surface_body(it_io, num_nodes)

  use specfem_par_movie_hdf5

  implicit none

  integer, intent(in)    :: it_io
  integer, intent(in)    :: num_nodes
  ! local parameters
  integer :: i
  character(len=20) :: it_str
  character(len=MAX_STRING_LEN) :: fname_xdmf_surf
  character(len=MAX_STRING_LEN) :: fname_h5_data_surf_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  !#TODO: hdf5 i/o server
  ! redefinition for no ioserver case
  !if (fname_xdmf_surf == '' ) fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface.xmf"

  fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface.xmf"
  fname_h5_data_surf_xdmf = "./movie_surface.h5"    ! relative to movie_surface.xmf file
  ! this seems to point to a wrong directory:
  !   fname_h5_data_surf_xdmf = trim(OUTPUT_FILES) // "/movie_surface.h5"

  ! open xdmf file
  open(unit=xdmf_surf, file=trim(fname_xdmf_surf), status='old', recl=256)

  ! skip lines till the position where we want to write new information
  do i = 1, surf_xdmf_pos
    read(xdmf_surf, *)
  enddo

  write(it_str, "(i6.6)") it_io

  write(xdmf_surf,*) '<Grid Name="surf_mov" GridType="Uniform">'
  write(xdmf_surf,*) '<Time Value="'//trim(r2c(sngl((it_io-1)*DT-t0)))//'" />'
  write(xdmf_surf,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
  write(xdmf_surf,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'
  write(xdmf_surf,*) '<Attribute Name="ux" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/ux'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Attribute>'
  write(xdmf_surf,*) '<Attribute Name="uy" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/uy'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Attribute>'
  write(xdmf_surf,*) '<Attribute Name="uz" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/uz'
  write(xdmf_surf,*) '</DataItem>'
  write(xdmf_surf,*) '</Attribute>'
  write(xdmf_surf,*) '</Grid>'
  ! 20 lines

  ! file finish
  write(xdmf_surf,*) '</Grid>'
  write(xdmf_surf,*) '</Domain>'
  write(xdmf_surf,*) '</Xdmf>'

  close(xdmf_surf)

  ! updates file record position
  surf_xdmf_pos = surf_xdmf_pos + 20

  end subroutine write_xdmf_surface_body

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wmo_movie_save_volume_snapshot_hdf5()

#ifdef USE_HDF5

  use specfem_par_movie_hdf5

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:), allocatable :: d_p
  integer :: ier
  logical :: pressure_io, divglob_io, div_io, veloc_io, curl_io
  !#TODO: hdf5 i/o server
  !integer :: req_count

  ! initializes
  pressure_io = .false.
  divglob_io = .false.
  div_io = .false.
  veloc_io = .false.
  curl_io = .false.

  !#TODO: hdf5 i/o server
  !req_count = 1

  ! prepares hdf5 file
  if (it == NTSTEP_BETWEEN_FRAMES) then
    !#TODO: hdf5 i/o server
    !if (HDF5_IO_NNODES == 0) then
      call prepare_vol_movie_hdf5()
    !endif
  endif

  ! saves fields
  if (ACOUSTIC_SIMULATION) then
    ! outputs pressure field for purely acoustic simulations
    if (.not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) then
      ! converts wavefield_pressure from local to global level
      allocate(d_p(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating d_p array'
      call elm2node_base(wavefield_pressure, d_p)
      ! saves pressure field
      !#TODO: hdf5 i/o server
      !if (HDF5_IO_NNODES > 0) then
      !  ! send pressure_loc
      !  call isend_cr_inter(d_p,NGLOB_AB,dest_ionod,io_tag_vol_pres,req_dump_vol(req_count))
      !  req_count = req_count+1
      !else
        ! direct io
        pressure_io = .true.
        call write_vol_data_hdf5(d_p,"pressure")
      !endif
      deallocate(d_p)
    endif
  endif ! acoustic

  if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
    ! allocates div/curl arrays on global nodes
    if (.not. allocated(div_on_node)) then
      allocate(div_on_node(NGLOB_AB), &
               curl_x_on_node(NGLOB_AB), &
               curl_y_on_node(NGLOB_AB), &
               curl_z_on_node(NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1737')
      if (ier /= 0) stop 'error allocating array movie div and curl'
      div_on_node(:) = 0._CUSTOM_REAL
      curl_x_on_node(:) = 0._CUSTOM_REAL
      curl_y_on_node(:) = 0._CUSTOM_REAL
      curl_z_on_node(:) = 0._CUSTOM_REAL
    endif

    ! saves div field
    !#TODO: hdf5 i/o server
    !if (HDF5_IO_NNODES > 0) then
    !  ! send div_glob
    !  call isend_cr_inter(div_glob,NGLOB_AB,dest_ionod,io_tag_vol_divglob,req_dump_vol(req_count))
    !  req_count = req_count+1
    !else
      ! direct io
      divglob_io = .true.
      call write_vol_data_hdf5(div_glob,"div_glob")
    !endif

    ! div and curl on elemental level
    call elm2node_base(div,div_on_node)
    call elm2node_base(curl_x,curl_x_on_node)
    call elm2node_base(curl_x,curl_x_on_node)
    call elm2node_base(curl_x,curl_x_on_node)

    ! writes our divergence/curl
    !#TODO: hdf5 i/o server
    !if (HDF5_IO_NNODES > 0) then
    !  ! div
    !  call isend_cr_inter(div_on_node,NGLOB_AB,dest_ionod,io_tag_vol_div,req_dump_vol(req_count))
    !  req_count = req_count+1
    !  ! writes out curl
    !  call isend_cr_inter(curl_x_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curlx,req_dump_vol(req_count))
    !  req_count = req_count+1
    !  call isend_cr_inter(curl_y_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curly,req_dump_vol(req_count))
    !  req_count = req_count+1
    !  call isend_cr_inter(curl_z_on_node,NGLOB_AB,dest_ionod,io_tag_vol_curlz,req_dump_vol(req_count))
    !  req_count = req_count+1
    !else
      ! direct io
      div_io = .true.
      curl_io = .true.
      call write_vol_data_hdf5(div_on_node,"div")
      call write_vol_data_hdf5(curl_x_on_node,"curl_x")
      call write_vol_data_hdf5(curl_y_on_node,"curl_y")
      call write_vol_data_hdf5(curl_z_on_node,"curl_z")
    !endif
  endif ! elastic or poroelastic

  ! velocity
  if (.not. allocated(velocity_x_on_node)) then
    allocate(velocity_x_on_node(NGLOB_AB), &
             velocity_y_on_node(NGLOB_AB), &
             velocity_z_on_node(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1733')
    if (ier /= 0) stop 'error allocating array movie velocity_x etc.'
    velocity_x_on_node(:) = 0._CUSTOM_REAL
    velocity_y_on_node(:) = 0._CUSTOM_REAL
    velocity_z_on_node(:) = 0._CUSTOM_REAL
  endif

  ! save wavefield from local (NGLLX,NGLLY,NGLLZ,NSPEC) to global (NGLOB) indexing
  ! given the wavefield is smooth and has no discontinuities at element boundaries, this shouldn't affect the wavefield shapes

  ! uses velocity_x,.. as temporary arrays to store velocity on all GLL points
  ! copy values on 1d array
  call elm2node_base(wavefield_x, velocity_x_on_node)
  call elm2node_base(wavefield_y, velocity_y_on_node)
  call elm2node_base(wavefield_z, velocity_z_on_node)

  ! saves wavefield components
  !#TODO: hdf5 i/o server
  !if (HDF5_IO_NNODES > 0) then
  !  call isend_cr_inter(velocity_x_on_node,NGLOB_AB,dest_ionod,io_tag_vol_velox,req_dump_vol(req_count))
  !  req_count = req_count+1
  !  call isend_cr_inter(velocity_y_on_node,NGLOB_AB,dest_ionod,io_tag_vol_veloy,req_dump_vol(req_count))
  !  req_count = req_count+1
  !  call isend_cr_inter(velocity_z_on_node,NGLOB_AB,dest_ionod,io_tag_vol_veloz,req_dump_vol(req_count))
  !  req_count = req_count+1
  !else
    ! direct io
    veloc_io = .true.
    call write_vol_data_hdf5(velocity_x_on_node,"veloc_x")
    call write_vol_data_hdf5(velocity_y_on_node,"veloc_y")
    call write_vol_data_hdf5(velocity_z_on_node,"veloc_z")
  !endif

  !#TODO: hdf5 i/o server
  ! store the number of mpi_isend reqs
  !n_req_vol = req_count-1

  if (it == NTSTEP_BETWEEN_FRAMES) then
    !#TODO: hdf5 i/o server
    !if (HDF5_IO_NNODES == 0) then
      ! create xdmf header
      call write_xdmf_vol_hdf5(pressure_io, divglob_io, div_io, veloc_io, curl_io)
    !endif
  endif

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 routine wmo_movie_save_volume_snapshot_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'wmo_movie_save_volume_snapshot_hdf5() called without compilation support'

#endif

  end subroutine wmo_movie_save_volume_snapshot_hdf5

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_HDF5

  subroutine prepare_vol_movie_hdf5()

  use specfem_par_movie_hdf5

  implicit none

  integer, dimension(:,:), allocatable :: elm_conn_loc
  character(len=64) :: group_name
  character(len=64) :: dset_name
  integer           :: iproc, ier, nglob_all, nspec_all
  character(len=MAX_STRING_LEN) :: fname_h5_data_vol
  logical,parameter :: if_collective = .true. ! in io_server.f90

  allocate(nglob_par_proc_nio(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nglob_par_proc')
  if (ier /= 0) stop 'error allocating arrays for nglob_par_proc'
  nglob_par_proc_nio(:) = 0
   allocate(nelm_par_proc_nio(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nelm_par_proc')
  if (ier /= 0) stop 'error allocating arrays for nelm_par_proc'
  nelm_par_proc_nio(:) = 0
  allocate(nglob_offset(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nglob_offset')
  if (ier /= 0) stop 'error allocating arrays for nglob_offset'
  nglob_offset(:) = 0
  allocate(nelm_offset(0:NPROC-1), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nelm_offset')
  if (ier /= 0) stop 'error allocating arrays for nelm_offset'
  nelm_offset(:) = 0
  allocate(elm_conn_loc(9,NSPEC_AB*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array elm_conn_loc')
  if (ier /= 0) stop 'error allocating arrays for elm_conn_loc'
  elm_conn_loc(:,:) = 0

  fname_h5_data_vol = trim(OUTPUT_FILES) // "/movie_volume.h5"

  ! initializes
  call h5_initialize()

  ! group for storing node coordinates and mesh element connectivity
  group_name = "mesh"

  call gather_all_all_singlei(NSPEC_AB,nelm_par_proc_nio,NPROC)
  call gather_all_all_singlei(NGLOB_AB,nglob_par_proc_nio,NPROC)

  nglob_all = sum(nglob_par_proc_nio(0:NPROC-1))
  nspec_all = sum(nelm_par_proc_nio(0:NPROC-1))

  nglob_offset(0) = 0
  nelm_offset(0)  = 0
  do iproc = 1, NPROC-1
    nglob_offset(iproc) = sum(nglob_par_proc_nio(0:iproc-1))
    nelm_offset(iproc)  = sum(nelm_par_proc_nio(0:iproc-1))
  enddo
  do iproc = 1,NPROC-1
    nelm_offset(iproc) = nelm_offset(iproc)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)
  enddo

  ! create connectivity dataset
  call get_conn_for_movie(elm_conn_loc, nglob_offset(myrank))

  if (myrank == 0) then
   ! create a hdf5 file
    call h5_create_file(fname_h5_data_vol)

    ! create coordinate dataset
    call h5_open_or_create_group(group_name)

    dset_name = "elm_conn"
    call h5_create_dataset_gen_in_group(dset_name, (/9,nspec_all*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)/), 2, 1)
    dset_name = "x"
    call h5_create_dataset_gen_in_group(dset_name, (/nglob_all/), 1, CUSTOM_REAL)
    dset_name = "y"
    call h5_create_dataset_gen_in_group(dset_name, (/nglob_all/), 1, CUSTOM_REAL)
    dset_name = "z"
    call h5_create_dataset_gen_in_group(dset_name, (/nglob_all/), 1, CUSTOM_REAL)

    call h5_close_group()
    call h5_close_file()
  endif
  call synchronize_all()

  call h5_open_file_p(fname_h5_data_vol)
  call h5_open_group(group_name)

  dset_name = "elm_conn"
  call h5_write_dataset_collect_hyperslab_in_group(dset_name,elm_conn_loc,(/0,nelm_offset(myrank)/),if_collective)
  dset_name = "x"
  call h5_write_dataset_collect_hyperslab_in_group(dset_name,xstore,(/nglob_offset(myrank)/),if_collective)
  dset_name = "y"
  call h5_write_dataset_collect_hyperslab_in_group(dset_name,ystore,(/nglob_offset(myrank)/),if_collective)
  dset_name = "z"
  call h5_write_dataset_collect_hyperslab_in_group(dset_name,zstore,(/nglob_offset(myrank)/),if_collective)

  call h5_close_group()
  call h5_close_file()
  call h5_finalize()

  end subroutine prepare_vol_movie_hdf5

#endif
!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_HDF5

  subroutine write_vol_data_hdf5(darr, dset_name)

  use specfem_par_movie_hdf5

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: darr
  character(len=*), intent(in)                     :: dset_name
  ! local parameters
  integer :: dim
  character(len=64) :: group_name
  character(len=12) :: tempstr
  character(len=MAX_STRING_LEN) :: fname_h5_data_vol

  integer, parameter :: rank = 1
  logical, parameter :: if_collective = .true.

  fname_h5_data_vol = trim(OUTPUT_FILES) // "/movie_volume.h5"

  ! initializes
  call h5_initialize()

  ! create a group for each io step
  write(tempstr, "(i6.6)") it
  group_name = "it_" // trim(tempstr)

  dim = sum(nglob_par_proc_nio(:))

  if (myrank == 0) then
    ! create it group
    call h5_open_file(fname_h5_data_vol)

    ! check if group_name exists
    call h5_open_or_create_group(group_name)

    ! create dataset
    call h5_create_dataset_gen_in_group(dset_name, (/dim/), rank, CUSTOM_REAL)
    call h5_close_group()
    call h5_close_file()
  endif

  call synchronize_all()

  ! collective write
  call h5_open_file_p(fname_h5_data_vol)
  call h5_open_group(group_name)
  call h5_write_dataset_collect_hyperslab_in_group(dset_name, &
                                darr, (/nglob_offset(myrank)/), if_collective)

  call h5_close_group()
  call h5_close_file()
  call h5_finalize()

  end subroutine write_vol_data_hdf5

#endif

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_HDF5

  subroutine write_xdmf_vol_hdf5(pressure_io, divglob_io, div_io, veloc_io, curl_io)

  use specfem_par_movie_hdf5

  implicit none

  logical, intent(in) :: pressure_io, divglob_io, div_io, veloc_io, curl_io

  ! local parameters
  integer :: i, ii
  character(len=20) :: it_str,nelm_str,nglo_str
  character(len=MAX_STRING_LEN) :: fname_xdmf_vol
  character(len=MAX_STRING_LEN) :: fname_h5_data_vol_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  ! writeout xdmf file for volume movie
  fname_xdmf_vol = trim(OUTPUT_FILES) // "/movie_volume.xmf"
  fname_h5_data_vol_xdmf = "./movie_volume.h5"  ! relative to movie_volume.xmf file
  ! this seems to point to a wrong director:
  !   fname_h5_data_vol_xdmf = trim(OUTPUT_FILES) // "/movie_volume.h5"

  open(unit=xdmf_vol, file=trim(fname_xdmf_vol), recl=256)

  ! definition of topology and geometry
  ! refer only control nodes (8 or 27) as a coarse output
  ! data array need to be extracted from full data array on GLL points
  write(xdmf_vol,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_vol,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_vol,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
  write(xdmf_vol,*) '<Domain name="mesh">'

  ! loop for writing information of mesh partitions
  nelm_str = i2c(sum(nelm_par_proc_nio(:))*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))
  nglo_str = i2c(sum(nglob_par_proc_nio(:)))

  write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm_str)//'">'
  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                         //trim(nelm_str)//' 9">'
  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/elm_conn'
  write(xdmf_vol,*) '</DataItem>'
  write(xdmf_vol,*) '</Topology>'
  write(xdmf_vol,*) '<Geometry GeometryType="X_Y_Z">'
  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/x'
  write(xdmf_vol,*) '</DataItem>'
  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/y'
  write(xdmf_vol,*) '</DataItem>'
  write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
  write(xdmf_vol,*) '       '//trim(fname_h5_data_vol_xdmf)//':/mesh/z'
  write(xdmf_vol,*) '</DataItem>'
  write(xdmf_vol,*) '</Geometry>'
  write(xdmf_vol,*) '<Grid Name="time_col" GridType="Collection" CollectionType="Temporal">'

  do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ii = i*NTSTEP_BETWEEN_FRAMES
    write(it_str, "(i6.6)") ii

    write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
    write(xdmf_vol,*) '<Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
    write(xdmf_vol,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
    write(xdmf_vol,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'

    if (pressure_io) then
      write(xdmf_vol,*) '<Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/pressure'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
    endif

    if (divglob_io) then
      write(xdmf_vol,*) '<Attribute Name="div_glob" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/div_glob'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
    endif

    if (div_io) then
      write(xdmf_vol,*) '<Attribute Name="div" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/div'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
    endif

    if (veloc_io) then
      write(xdmf_vol,*) '<Attribute Name="veloc_x" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/veloc_x'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
      write(xdmf_vol,*) '<Attribute Name="veloc_y" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/veloc_y'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
      write(xdmf_vol,*) '<Attribute Name="veloc_z" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/veloc_z'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
    endif

    if (curl_io) then
      write(xdmf_vol,*) '<Attribute Name="curl_x" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/curl_x'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
      write(xdmf_vol,*) '<Attribute Name="curl_y" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/curl_y'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
      write(xdmf_vol,*) '<Attribute Name="curl_z" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
      write(xdmf_vol,*) '      '//trim(fname_h5_data_vol_xdmf)//':/it_'//trim(it_str)//'/curl_z'
      write(xdmf_vol,*) '</DataItem>'
      write(xdmf_vol,*) '</Attribute>'
    endif

    write(xdmf_vol,*) '</Grid>'
  enddo

  ! file finish
  write(xdmf_vol,*) '</Grid>'
  write(xdmf_vol,*) '</Domain>'
  write(xdmf_vol,*) '</Xdmf>'

  close(xdmf_vol)

  end subroutine write_xdmf_vol_hdf5

#endif
