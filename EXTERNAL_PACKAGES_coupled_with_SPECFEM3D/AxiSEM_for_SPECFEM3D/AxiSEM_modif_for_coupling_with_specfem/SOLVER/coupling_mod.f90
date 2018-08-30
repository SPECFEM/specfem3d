module coupling_mod
!==================================================
! BEGIN COPY PASTE FROM lateral_heterogeneities

  use global_parameters
  use data_heterogeneous
  use data_io
  use data_proc
  use utlity, only: compute_coordinates, scoord, zcoord
  use data_source, only: rot_src

  use commun, only: barrier
  use data_mesh, only: npol, nelem, eltype, lnods, crd_nodes, npoin
  use data_mesh, only: ielsolid, ielfluid, nel_solid, nel_fluid
  use analytic_mapping, only: jacobian                                  !!! SB Previously in geom_transf

  implicit none

! END COPY PASTE (SB)
!==================================================


!==================================================
! MODULE SHARED VARIABLES
  ! Number of unique boundary points, nb of boundary points
  integer :: npt_box_file,nb_elm_to_store
  integer, dimension(:), allocatable :: is_in_box,id_elm_to_store,id_glob_to_store

  ! Geographic coordinates of box's points (r, theta, phy)
  real(kind=dp) :: rmin, rmax, thmin, thmax
  real(kind=dp), dimension(:), allocatable :: rbox, thbox, phbox
  real(kind=dp), dimension(:,:), allocatable :: szbox
  real(kind=realkind), dimension(:,:,:), allocatable :: buff_to_dump
  real(kind=dp), dimension(:,:,:), allocatable :: lambda_cp,mu_cp,rho_cp

  logical :: storage_for_recip_KH_integral

! END SHARED (SB)
!==================================================

contains

!===============================================================================
! Init coupling
  subroutine initialize_coupling

    integer :: j

    if (dump_wavefields) then

       if (dump_type == 'coupling') then

          ! Find elements base on min/max of colat and radius given by the user
          rmin  = kwf_rmin !* 1000       ! to meters
          rmax  = kwf_rmax !* 1000
          thmin = kwf_thetamin !* pi / 180.  ! to rad
          thmax = kwf_thetamax !* pi / 180.

          allocate(rbox(4), thbox(4), phbox(4))

          rbox(1)  = rmin
          thbox(1) = thmin
          rbox(2)  = rmin
          thbox(2) = thmax
          rbox(3)  = rmax
          thbox(3) = thmin
          rbox(4)  = rmax
          thbox(4) = thmax
          phbox(:) = 0.


          ! Must rotate these coordinates when source is not a the pole
          if (rot_src) then
             write(*,*) mynum, 'rotate since source is not beneath north pole'

             do j=1,4
                call rotate_box(rbox(j), thbox(j), phbox(j))
             enddo

             rmin  = minval(rbox)
             rmax  = maxval(rbox)
             thmin = minval(thbox)
             thmax = maxval(thbox)

             ! Check min/max values again if rotation has been applied
             write(*,*) mynum, 'r min/max after rotation:', rmin / 1000., rmax / 1000.
             write(*,*) mynum, 'th min/max after rotation:', thmin / pi * 180., &
                  thmax / pi * 180.
          endif

          allocate(szbox(2,1:4))
          szbox(1,:) = rbox(:) * sin(thbox(:))
          szbox(2,:) = rbox(:) * cos(thbox(:))

          call find_boundary_elements

          deallocate(rbox,thbox,szbox)

       else if (dump_type == 'coupling_box') then

          ! Read ../input_box.txt and get min/max or colatitudes and radius
          ! then find elements
          call read_boundary_coordinates
          ! this is ok (last version of AxiSEM

       endif

    endif

  end subroutine initialize_coupling
!--------------------------------------------------------------------------------

!================================================================================
! Read box's boundary coordinates (r,theta,phi) and sort them to only keep unique
! coordinates
  subroutine read_boundary_coordinates

    integer :: j
    character(len=100) :: boxpoints_file
    real, parameter :: delta_r=10000., delta_th=0.05


    boxpoints_file = '../input_box.txt'

    !--------------------------------------------------
    ! Read discrete box's boundary point file
    open(unit=91, file=trim(boxpoints_file), status='old')
    read(91,*) npt_box_file

    allocate(rbox(1:npt_box_file), thbox(1:npt_box_file), phbox(1:npt_box_file))

    do j=1, npt_box_file
       read(91,*) rbox(j), thbox(j), phbox(j)
    enddo

    close(unit=91)

    rbox  = rbox * 1000.
    thbox = (90. - thbox) * pi / 180.      ! Read latitudes instead of colatitudes
    phbox = phbox * pi / 180.

    !--------------------------------------------------
    ! Sort (r,theta,phi) coordinates


    !-------------------------------------------------------
    ! Rotate the whole set of coordinate with the source term
    rmin  = minval(rbox)
    rmax  = maxval(rbox) !+ delta_r
    thmin = minval(thbox)
    thmax = maxval(thbox) !+ delta_th

    ! Check min/max values of r and theta before rotation
    write(*,*) mynum, 'r min/max:', rmin / 1000., rmax / 1000.
    write(*,*) mynum, 'th min/max:', thmin / pi * 180., thmax / pi * 180.

    !write(*,*) 'ROTATION ? ', rot_src
    if (rot_src) then
       write(*,*) mynum, 'rotate since source is not beneath north pole'

       do j=1, npt_box_file
       !   if (j==1) then
       !      write(*,*)'Before: ',thbox(j)*180/pi,phbox(j)*180/pi
       !   endif
          call rotate_box(rbox(j), thbox(j), phbox(j))
       !   if (j==1) then
       !      write(*,*)'After: ',thbox(j)*180/pi,phbox(j)*180/pi
       !   endif
       enddo

       rmin  = minval(rbox)
       rmax  = maxval(rbox)
       thmin = minval(thbox)
       thmax = maxval(thbox)

       ! Check min/max values again if rotation has been applied
       write(*,*) mynum, 'r min/max after rotation:', rmin / 1000., rmax / 1000.
       write(*,*) mynum, 'th min/max after rotation:', thmin / pi * 180., &
            thmax / pi * 180.
    endif

    !!!!!!!!!!! MAYBE SHOULD WE PASTE vtk plot discrete inputs ALSO (SB) ...
    !!!!! TROUBLES WITH rhetmin, rhetmax, thhetmin, thhetmax... where is the declaration
    !!!!!      what is the difference with rmin, rmax previously defined ?  (SB)


    !--------------------------------------------------
    ! Move to cylindrical coordinates
    allocate(szbox(2,1:npt_box_file))
    szbox(1,:) = rbox(:) * sin(thbox(:))
    szbox(2,:) = rbox(:) * cos(thbox(:))

    ! TEMPORARY FIND THE ELEMENS HERE (SB)
    call find_boundary_elements

    deallocate(rbox,thbox,szbox) !!!! Missing phbox...

  end subroutine read_boundary_coordinates
!================================================================================


!================================================================================
! Find elements corresponding to unique box's boundary coordinates using the
! kdtree data structure
  subroutine find_boundary_elements

    integer :: iel, ibox, iel_solid
    real(kind=dp) :: s, z, r1, r2, r3, r4, th1, th2, th3, th4, r, th, thetamin, thetamax
    real(kind=dp) :: smin,smax,zmin,zmax
    real(kind=dp),parameter  :: eps=1e-5

    !! FOR VTK
    real, allocatable :: data_for_vtk(:)
    integer, allocatable :: lnodes_vtk(:,:),solid_id(:)
    character (len=256) vtk_file_name
    !!

    ! Initialization of loouking area bounds
    rhetmin = 1.d10
    rhetmax = 0.d0
    thhetmin = pi
    thhetmax = 0.d0

    ! Allocate is_in_box vector
    if (.not. allocated(is_in_box)) allocate(is_in_box(nelem))
    is_in_box(:) = 0

    allocate(solid_id(nelem))
    solid_id(:)=0

    !! VTK arrays
    allocate(data_for_vtk(nelem))
    allocate(lnodes_vtk(nelem,4))
    data_for_vtk(:)=0.
    !!

    !! mapping all elements to solid elements
    do iel_solid=1,nel_solid
       do iel=1, nelem
          if (ielsolid(iel_solid) == iel) then
             solid_id(iel)=iel_solid
             exit
          endif
       enddo
    enddo

    write(*,*) mynum, 'locate GLL points within heterogeneous regions & '


    ! Loop over all elements to find which ones contain box's points
    do iel=1, nelem

       lnodes_vtk(iel,1)=lnods(iel,3)
       lnodes_vtk(iel,2)=lnods(iel,1)
       lnodes_vtk(iel,3)=lnods(iel,7)
       lnodes_vtk(iel,4)=lnods(iel,5)



       smin=1.d15
       zmin=1.d15
       smax=0.d0
       zmax=-1.d15

       ! Compute iel element's corner coordinates
       call compute_coordinates(s, z, r1, th1, iel, npol, npol)
       smin=min(smin,s)
       smax=max(smax,s)
       zmin=min(zmin,z)
       zmax=max(zmax,z)

       call compute_coordinates(s, z, r2, th2, iel, 0, 0)
       smin=min(smin,s)
       smax=max(smax,s)
       zmin=min(zmin,z)
       zmax=max(zmax,z)

       call compute_coordinates(s, z, r3, th3, iel, 0, npol)
       smin=min(smin,s)
       smax=max(smax,s)
       zmin=min(zmin,z)
       zmax=max(zmax,z)

       call compute_coordinates(s, z, r4, th4, iel, npol, 0)
       smin=min(smin,s)
       smax=max(smax,s)
       zmin=min(zmin,z)
       zmax=max(zmax,z)

       ! Find max and min of element corner coordinates
       r = max(max(r1,r2), max(r3, r4))
       th = max(max(th1,th2), max(th3, th4))

       ! If the element is in the "restricted area" (defined by box's edges)
       if (r >= rmin .and. th >= thmin) then
          r = min(min(r1,r2), min(r3, r4))
          th = min(min(th1,th2), min(th3, th4))
          if (r <= rmax .and. th <= thmax) then
             !to plot all GLL points in elements effected:
!             rhetmax = max(max(max(r1,r2), max(r3, r4)), rhetmax)
!             thhetmax = max(max(max(th1,th2), max(th3, th4)), thhetmax)
!             rhetmin = min(min(min(r1,r2), min(r3, r4)), rhetmin)
!             thhetmin = min(min(min(th1,th2), min(th3, th4)), thhetmin)

             ! Loop over box's points

             select case (dump_type)    !! SB
             case ('coupling_box')
                     do ibox = 1, npt_box_file


                        !!!!!!! LATER NEED TO CHECK THE SHAPE OF THE ELEMENT (eltype)
                        !!!!!!! MUST COMPUTE THE JACOBIAN OF THE ELEMENT (eta,xi_k,...)
                        !!!!!!! => will be tricky, must take attention (solid/fluid,renumbering)
                        !!!!!!! Look at compute_coordinates_mesh (get_mesh.f90), crd_nodes and lnodes (data_mesh.f90) inode=1,8 => informations on how to use in compute_coordinates
                        !!!!!! c
                        if (       szbox(1,ibox) >= smin-eps .and. szbox(1,ibox) <= smax+eps &
                             .and. szbox(2,ibox) >= zmin-eps .and. szbox(2,ibox) <= zmax+eps     ) then !

                           is_in_box(iel)=is_in_box(iel)+1
                           !! FOR NOW ASSUME THAT WE ARE IN SOLID REGION ONLY
                           data_for_vtk(iel)=1.

                           !write(*,*)'Found element : ',iel,' of coordinates : ',smin,smax,zmin,zmax,' for box point : ',szbox(1,ibox),szbox(2,ibox)
                           exit

                        endif

                     enddo
             case ('coupling')
                        is_in_box(iel)=is_in_box(iel)+1
                        !! FOR NOW ASSUME THAT WE ARE IN SOLID REGION ONLY
                         data_for_vtk(iel)=1.
             end select   !! Add cade to disticnt box and aeea (SB)
          endif
       endif
    enddo

    ! count elements to store
    nb_elm_to_store=0
    do  iel=1, nelem
       if (is_in_box(iel) > 0) nb_elm_to_store= nb_elm_to_store+1
    enddo
    allocate(id_elm_to_store(nb_elm_to_store),id_glob_to_store(nb_elm_to_store))
    allocate(buff_to_dump(0:npol,0:npol,nb_elm_to_store))

    ibox=0
    do  iel=1, nelem
       if (is_in_box(iel) > 0) then
          ibox=ibox+1
          id_elm_to_store(ibox)=solid_id(iel)
          id_glob_to_store(ibox)=iel
       endif
    enddo


!! VM VM just writing to check
    !if (mynum == 0) then
    !   write(*,*) mynum
    !   write(*,*) id_elm_to_store
    !   write(*,*)
    !endif

    !if (mynum == 1) then
    !   write(*,*) mynum
    !   write(*,*) id_elm_to_store
    !   write(*,*)
    !endif

    ! Paraview output

    !!!! AJOUTER LE TEST SB
    write(vtk_file_name,'(a10,i5.5,a4)')'Found_Elem',mynum,'.vtk'
    call  writeWTKCell(vtk_file_name ,lnodes_vtk,crd_nodes,data_for_vtk,npoin,nelem,3)

    write(*,*) mynum, 'DONE loading box points'

    !! VM VM Store information in order to reconstruct local mesh.
    write(vtk_file_name,'(a10,i5.5,a4)')'parameters_for_vtk_',mynum,'.par'
    open(49,file=trim(vtk_file_name))
    write(49,*)  nb_elm_to_store,ibeg,iend
    close(49)

    !! VM VM store local meshes
    call dump_wavefields_mesh_1d_cp

  end subroutine find_boundary_elements
!================================================================================

!================================================================================
! Write vtk
  subroutine writeWTKCell(filename,t,p,u,N,M,d)

    implicit none

    integer N,M,d
    real(kind=dp) p(N,d-1)
    real u(M)
    integer t(M,d+1)
    character(len=256) filename,string,string1,string2
    character(len=256) string3,string4

    integer i

    open(49,file=trim(filename))

    write(49,'(a26)') '# vtk DataFile Version 2.0'
    write(49,'(a25)') 'Unstructured Grid Example'
    write(49,'(a5)') 'ASCII'
    write(49,'(a25)') 'DATASET UNSTRUCTURED_GRID'
    write(string,'(i30)') N
    string="POINTS "//trim(adjustl(string))//" float"
    write(49,'(a)') trim(string)

    do i=1,N
       write(string,'(f45.20)')   p(i,1)
       write(string1,'(f45.20)')  p(i,2)
       write(string2,'(f45.20)')  0. !p(i,3)
       string=trim(adjustl(string))//" "//trim(adjustl(string1))&
            //" "//trim(adjustl(string2))
       write(49,'(a)') trim(string)
    enddo
    write(string,'(i30)') M
    write(string1,'(i30)') M*(d+2)
    string="CELLS "//trim(adjustl(string))//" "//trim(adjustl(string1))
    write(49,'(a)') trim(string)
    do i=1,M
       write(string, '(i40)') d+1
       write(string1,'(i40)') t(i,1)-1
       write(string2,'(i40)') t(i,2)-1
       write(string3,'(i40)') t(i,3)-1
       write(string4,'(i40)') t(i,4)-1
       string=trim(adjustl(string))//"  "//trim(adjustl(string1))&
            //" "//trim(adjustl(string2)) //" "//trim(adjustl(string3))&
            //" "//trim(adjustl(string4))
       write(49,'(a)') trim(string)
    enddo

    write(string,'(i30)') M
    string="CELL_TYPES "//trim(adjustl(string))
    write(49,'(a)') trim(string)
    do i=1,M
       write(49,'(a2)') '9' !
    enddo

    write(string,'(i30)') M
    string="CELL_DATA "//trim(adjustl(string))
    write(49,'(a)') trim(string)
    write(49,'(a17)')'SCALARS u float 1'
    write(49,'(a20)')'LOOKUP_TABLE default'
    do i=1,M
       write(string,'(f40.20)') u(i)
       write(49,'(a)') trim(adjustl(string))
    enddo
    close(49)

  end subroutine writeWTKCell
!--------------------------------------------------------------------------------

!================================================================================
! Routine to ratotate the box
  subroutine rotate_box(r,th,ph)

    real(kind=dp)   ,intent(inout) :: r,th,ph
    real(kind=dp)    :: x_vec(3), x_vec_rot(3), r_r

    x_vec(1) = r * dsin(th) * dcos(ph)
    x_vec(2) = r * dsin(th) * dsin(ph)
    x_vec(3) = r * dcos(th)

    x_vec_rot = matmul(trans_rot_mat,x_vec)

    r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
    th = dacos((x_vec_rot(3)  + smallval_dble )/ ( r_r + smallval_dble) )
    ph = atan2(x_vec_rot(2),x_vec_rot(1))

  end subroutine rotate_box
!--------------------------------------------------------------------------------

!================================================================================
! Dumping routine
  subroutine dump_field_1d_cp(f, filename, appisnap, n)

    use data_source, only: have_src, src_dump_type
    use data_mesh, only: nel_solid, nel_fluid

    integer, intent(in)                 :: n
    real(kind=realkind),intent(in)      :: f(0:,0:,:)
    character(len=16), intent(in)       :: filename
    character(len=4), intent(in)        :: appisnap
    integer i
    !real(kind=realkind)                 :: floc(0:size(f,1)-1, 0:size(f,2)-1, 1:size(f,3))

    do i=1,nb_elm_to_store
       buff_to_dump(:,:,i) = f(:,:,id_elm_to_store(i))
       !write(*,*) id_elm_to_store(i),lambda_cp(1,1,id_glob_to_store(i)),mu_cp(1,1,id_glob_to_store(i))
    enddo

    !if (src_dump_type == 'mask' .and. n==nel_solid) &
    !     call eradicate_src_elem_values(floc)

    !if (use_netcdf) then
    !    if (n==nel_solid) then
    !        call nc_dump_field_solid(pack(floc(ibeg:iend,ibeg:iend,:), .true.), &
    !                                 filename(2:))
    !    else if (n==nel_fluid) then
    !        call nc_dump_field_fluid(pack(floc(ibeg:iend,ibeg:iend,:), .true.), &
    !                                 filename(2:))
    !    else
    !        write(*,*) 'Neither solid nor fluid. What''s wrong here?'
    !        stop 2
    !    endif
    !else
!!$      open(unit=25000+mynum, file=datapath(1:lfdata)//filename//'_' &
!!$                                  //appmynum//'_'//appisnap//'.bindat', &
!!$           FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
    !write(*,*) 'dumping file ', datapath(1:lfdata)//filename//'_' &
    !     //appmynum//'.bindat'
    open(unit=25000+mynum, file=datapath(1:lfdata)//filename//'_' &
         //appmynum//'.bindat', &
         FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="APPEND")
    write(25000+mynum) pack(buff_to_dump(ibeg:iend,ibeg:iend,:), .true.)
    close(25000+mynum)
    !endif

  end subroutine dump_field_1d_cp
!-----------------------------------------------------------------------------------------


!================================================================================
! Dumping routine
  subroutine dump_wavefields_mesh_1d_cp
    real(kind=dp)   , dimension(:,:,:), allocatable :: ssol, zsol
    integer, allocatable :: lnods_to_store(:,:)
    integer iel,ipol,jpol


    ! solid part only for now
    allocate(ssol(ibeg:iend,ibeg:iend,nb_elm_to_store))
    allocate(zsol(ibeg:iend,ibeg:iend,nb_elm_to_store))
    allocate(lnods_to_store(nb_elm_to_store,4))

    !!write(*,*) 'nb_elm_to_store ::',nb_elm_to_store
    do iel=1,nb_elm_to_store

       lnods_to_store(iel,1)=lnods(id_elm_to_store(iel),3)
       lnods_to_store(iel,2)=lnods(id_elm_to_store(iel),1)
       lnods_to_store(iel,3)=lnods(id_elm_to_store(iel),7)
       lnods_to_store(iel,4)=lnods(id_elm_to_store(iel),5)

       do jpol=ibeg,iend
          do ipol=ibeg,iend
             ssol(ipol,jpol,iel) = scoord(ipol,jpol,ielsolid((id_elm_to_store(iel))))
             zsol(ipol,jpol,iel) = zcoord(ipol,jpol,ielsolid((id_elm_to_store(iel))))
             !!write(*,*) ipol,jpol,iel,ssol(ipol,jpol,iel)
          enddo
       enddo
    enddo

    open(unit=2500+mynum,file=datapath(1:lfdata)//'/mesh_sol_'&
         //appmynum//'.dat', &
         FORM="UNFORMATTED")
    write(2500+mynum) ssol(ibeg:iend,ibeg:iend,:),zsol(ibeg:iend,ibeg:iend,:)
    close(2500+mynum)


    open(unit=2500+mynum,file=datapath(1:lfdata)//'/lnods_sol_'&
         //appmynum//'.dat', &
         FORM="UNFORMATTED")
    write(2500+mynum)  lnods_to_store
    close(2500+mynum)



    deallocate(ssol,zsol)
    deallocate(lnods_to_store)

  end subroutine dump_wavefields_mesh_1d_cp


  subroutine  finalize_coupling
    open(10,file='nb_rec_to_read.par')
    write(10,*) nstrain
    close(10)
  end subroutine finalize_coupling

end module coupling_mod
