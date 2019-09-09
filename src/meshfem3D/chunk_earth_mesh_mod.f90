module chunk_earth_mod


  use constants, only: CUSTOM_REAL, MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN)      :: metric_system

  !! chunk parameters for metric systme
  double precision                                  :: xmin_chunk, xmax_chunk
  double precision                                  :: ymin_chunk, ymax_chunk
  double precision                                  :: zmin_chunk, zmax_chunk
  double precision                                  :: dx_elem, dy_elem, dz_elem
  double precision                                  :: dx_max, dy_max, dz_max

  logical                                           :: use_doubling=.false.
  integer                                           :: nb_doubling
  double precision,       dimension(:), allocatable :: layer_doubling, zlayer
  integer,                dimension(:), allocatable :: nzlayer


  !! private parameters :
  integer,                       private   :: ierr
  integer,                       private   :: ipos0, ipos1
  character(len=MAX_STRING_LEN), private   :: line, system_coordinate

  integer,                       private   :: nx_ref, ny_ref


  !! materiels
  integer,          private                              :: nb_mat
  integer,          private, dimension(:),   allocatable :: flag_acoustic_elastic
  double precision, private, dimension(:,:), allocatable :: material_prop

  !! domains
  integer,          private                              :: nb_dom
  double precision, private, dimension(:,:), allocatable :: domain_boundary
  integer,          private, dimension(:),   allocatable :: Imatetrial_ispec,  Imaterial_domain

  !! private working array
  integer,          private                              :: npoint_tot, nspec_tot
  integer,          private                              :: nspec_xmin,  nspec_xmax
  integer,          private                              :: nspec_ymin,  nspec_ymax
  integer,          private                              :: nspec_zmin,  nspec_zmax
  double precision, private, dimension(:),   allocatable :: x_mesh_point, y_mesh_point, z_mesh_point
  integer,          private, dimension(:,:), allocatable :: EtoV
  logical,          private, dimension(:,:), allocatable :: iboun

contains


!!##################################################################################################################################
!!                                  READ SPECFIC FILES FOR CHUNK OF THE EARTH MESH
!!##################################################################################################################################

    subroutine mesh_chunk_earth()

      !! SB SB COMMENT THIS IERR (variable alreafy exist in module)
      !! this makes gfortran complins because the one in the ;odule header is unsed
      !! integer :: ierr
      !! SB SB
      character(len=MAX_STRING_LEN) :: keyw

      open(27, file='DATA/meshfem3D_files/Mesh_Chunk_Par_file', action='read', iostat=ierr)
      if (ierr /= 0 ) then
         write(*,*) " ERROR : file DATA/meshfem3D_files/Mesh_Chunk_Par_file not found "
         stop
      endif

      do
         read(27,'(a)',end=99) line

         !! INDICES TO READ line -----------------------------------------------
         ipos0=index(line,':')+1
         ipos1=index(line,'#')-1
         if (ipos1 < 0 ) ipos1=len_trim(line)

         !! STORE KEYWORD ITEM -------------------------------------------------
         keyw=trim(adjustl(line(1:ipos0-2)))

         select case (trim(keyw))
         case ('system_coordinate')
            write(*,*) line
            system_coordinate=trim(adjustl(line(ipos0:ipos1)))
         end select

      enddo
      99 close(27)

      !! read input mesh parameters
      select case (trim(adjustl(system_coordinate)))
      case('meters')
         call read_metric_params()
         !! compute grid mesh the chunk
         call mesh_metric_chunk()
         !! define the final
         call create_mesh_metric_chunk()
      case('degree')
         !call read_spherical_params()
         write(*,*) "ERROR : not yet implemented system : ", trim(adjustl(system_coordinate))
      case default
         write(*,*) "ERROR : not implemented  system : ", trim(adjustl(system_coordinate))
      end select


      !! write ascii format same as geocubit.
      !call write_mesh_in_ascii()

    end subroutine mesh_chunk_earth

!!##################################################################################################################################
!!                       READ INPUT FOR METRIC SYSTEM
!!##################################################################################################################################

    subroutine read_metric_params()

      character(len=MAX_STRING_LEN) :: keyw
      integer :: ilayer, nx, ny, k, ntmp
      double precision :: dz
      integer :: i, iflag, imat, ier
      double precision :: vp, vs, rho, Q_Kappa, Q_mu, Aniso
      double precision :: x0, x1, y0, y1, z0, z1

      open(27, file='DATA/meshfem3D_files/Mesh_Chunk_Par_file', action='read')
      do
          read(27,'(a)',end=99) line

          !! INDICES TO READ line -----------------------------------------------
          ipos0=index(line,':')+1
          ipos1=index(line,'#')-1
          if (ipos1 < 0 ) ipos1=len_trim(line)

         !! STORE KEYWORD ITEM -------------------------------------------------
          keyw=trim(adjustl(line(1:ipos0-2)))

          select case (trim(keyw))
          case('xmin')
             read(line(ipos0:ipos1),*) xmin_chunk
          case('xmax')
             read(line(ipos0:ipos1),*) xmax_chunk
          case('ymin')
             read(line(ipos0:ipos1),*) ymin_chunk
          case('ymax')
             read(line(ipos0:ipos1),*) ymax_chunk
          case('zmin')
             read(line(ipos0:ipos1),*) zmin_chunk
          case('zmax')
             read(line(ipos0:ipos1),*) zmax_chunk
          case('dx_elem')
             read(line(ipos0:ipos1),*) dx_elem
          case('dy_elem')
             read(line(ipos0:ipos1),*) dy_elem
          case('dz_elem')
             read(line(ipos0:ipos1),*) dz_elem
          case('nb_doubling')
             use_doubling=.true.
             read(line(ipos0:ipos1),*) nb_doubling
             allocate(layer_doubling(nb_doubling),stat=ier)
             if (ier /= 0) call exit_MPI_without_rank('error allocating array 1240')
          case('layer_doubling')
             read(line(ipos0:ipos1),*) layer_doubling(:)
          case ('nb_material')
             read(line(ipos0:ipos1),*) nb_mat
             allocate(material_prop(nb_mat, 6),stat=ier)
             if (ier /= 0) call exit_MPI_without_rank('error allocating array 1241')
             allocate(flag_acoustic_elastic(nb_mat),stat=ier)
             if (ier /= 0) call exit_MPI_without_rank('error allocating array 1242')
             flag_acoustic_elastic(:)=-1
          case('material')
             read(line(ipos0:ipos1),*,iostat=ier) i, rho, vp, vs, Q_Kappa, Q_mu, Aniso, iflag
             if (ier /= 0) then
               print *,'Error while reading your input file in routine chunk_earth_mesh_mod.f90'
               print *
               print *,'We recently changed the input format from i, rho, vp, vs, Q_mu, Aniso, iflag'
               print *,'to i, rho, vp, vs, Q_Kappa, Q_mu, Aniso, iflag in order to add support for Q_Kappa.'
               print *,'It is likely that your input file still uses the old convention and needs to be updated.'
               print *,'If you do not know what value to add for Q_Kappa, add 9999., i.e negligible Q_Kappa attenuation'
               print *,'and then your results will be unchanged compared to older versions of the code.'
               print *
               stop 'Error in input file format in routine chunk_earth_mesh_mod.f90'
             endif
             if (i <= nb_mat) then
                material_prop(i,1) = rho
                material_prop(i,2) = vp
                material_prop(i,3) = vs
                material_prop(i,4) = Q_Kappa
                material_prop(i,5) = Q_mu
                material_prop(i,6) = Aniso
                flag_acoustic_elastic(i) = iflag
             else
                write(*,*) " Warning : material number ", i, " not used "
             endif
          case('nb_region')
             read(line(ipos0:ipos1),*) nb_dom
             allocate(domain_boundary(nb_dom,6),stat=ier)
             if (ier /= 0) call exit_MPI_without_rank('error allocating array 1243')
             allocate(Imaterial_domain(nb_dom),stat=ier)
             if (ier /= 0) call exit_MPI_without_rank('error allocating array 1244')
             domain_boundary(:,:)=0.
          case('region')
             read(line(ipos0:ipos1),*) i, x0, x1, y0, y1, z0, z1, imat
             if (i <= nb_dom) then
                domain_boundary(i,1)=x0;  domain_boundary(i,2)=x1
                domain_boundary(i,3)=y0;  domain_boundary(i,4)=y1
                domain_boundary(i,5)=z0;  domain_boundary(i,6)=z1
                Imaterial_domain(i) = imat
             else
                write(*,*) " Warning : domain number ",i, " not used "
             endif

          end select

       enddo
99     close(27)

       !! allocate arrays used in the current module ----------------------

       !! count the number of elements and the number of points

       !! size of element at the bottom of chunk
       dx_max = (2.d0**nb_doubling) * dx_elem
       dy_max = (2.d0**nb_doubling) * dy_elem
       dz_max = (2.d0**nb_doubling) * dz_elem

       !! comupute the used size to fit the domain dimension
       ntmp = floor((  xmax_chunk - xmin_chunk ) / dx_max) +1
       dx_max  = ( xmax_chunk - xmin_chunk ) / real(ntmp,8)

       !! comupute the used size to fit the domain dimension
       ntmp = floor((  ymax_chunk - ymin_chunk ) / dy_max) +1
       dy_max  = ( ymax_chunk - ymin_chunk ) / real(ntmp,8)

       !! number of element in first regular layer
       nx_ref =  floor((  xmax_chunk - xmin_chunk ) / dx_max ) + 1
       ny_ref =  floor((  ymax_chunk - ymin_chunk ) / dy_max ) + 1

       !! store the subdomain boundary
       allocate(zlayer(nb_doubling+2),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 1245')
       allocate(nzlayer(nb_doubling+1),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 1246')
       k=0
       do ilayer=2, nb_doubling+1
          k=k+1
          zlayer(ilayer)=layer_doubling(k)
       enddo
       zlayer(1)=zmin_chunk
       zlayer(nb_doubling+2)=zmax_chunk

       nspec_tot=0
       npoint_tot=0
       nx=nx_ref
       ny=ny_ref
       dz=dz_max
       !! compute number of elements needed for meshing
       do ilayer = 1, nb_doubling
          write(*,*)
          write(*,*) 'ilayer ', ilayer
          write(*,*) ' dz ', dz
          write(*,*) (zlayer(ilayer +1) - zlayer(ilayer))
          write(*,*) nint((zlayer(ilayer +1) - zlayer(ilayer)) / dz)

          nzlayer(ilayer) = nint( (zlayer(ilayer +1) - zlayer(ilayer)) / dz )
          !! to be sure that we have at least 2 elements in vertical direction
          if (nzlayer(ilayer) < 2) nzlayer(ilayer)=2

          nspec_tot = nspec_tot + (nzlayer(ilayer)-1)*nx*ny + 7*nx*ny
          npoint_tot = npoint_tot + nzlayer(ilayer)*(nx+1)*(ny+1) + 8*7*nx*ny

          write(*,*) ' nx, ny, nz ',nx, ny,  nzlayer(ilayer)

          nx=2*nx
          ny=2*ny
          dz = dz /2.d0

       enddo
       write(*,*)
       write(*,*)
       write(*,*)
       !! last layer without doubling
       ilayer = nb_doubling + 1
       nzlayer(ilayer) = nint( (zlayer(ilayer +1) - zlayer(ilayer)) / dz )
       !! to be sure that we have at least 2 elements in vertical direction
       if (nzlayer(ilayer) < 2) nzlayer(ilayer)=2
       nspec_tot = nspec_tot + (nzlayer(ilayer))*nx*ny
       npoint_tot = npoint_tot + (nzlayer(ilayer)+1)*(nx+1)*(ny+1)
       write(*,*)  ' nx, ny, nz ',  nx, ny, nzlayer(ilayer)


       write(*,*) 'zlayer ', zlayer(:)
       write(*,*) 'nzlayer ', nzlayer(:)
       !! allocate mesh arrays
       allocate(x_mesh_point(npoint_tot), y_mesh_point(npoint_tot), z_mesh_point(npoint_tot),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 1247')
       allocate(EtoV(8,nspec_tot), iboun(6,nspec_tot),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 1248')
       iboun(:,:)=.false.
       EtoV(:,:)=0

    end subroutine read_metric_params

!!##################################################################################################################################
!!                                  MESH THE CHUNK
!!##################################################################################################################################

    subroutine mesh_metric_chunk()

      integer :: nx, ny, nz, ispec, ipoint, ilayer, ier
      double precision :: dx, dy, dz, z
      double precision, dimension(:,:), allocatable :: top_surface, bottom_surface

      ispec=0
      ipoint=0
      nspec_xmin=0
      nspec_xmax=0
      nspec_ymin=0
      nspec_ymax=0
      nspec_zmin=0
      nspec_zmax=0

      !! mesh all layers with doubling
      do ilayer = 1, nb_doubling

         nx = nx_ref*(2**(ilayer-1))
         ny = ny_ref*(2**(ilayer-1))
         nz = nzlayer(ilayer)

         z = zlayer(ilayer)

         dx =  ( xmax_chunk - xmin_chunk ) / real(nx,8)
         dy =  ( ymax_chunk - ymin_chunk ) / real(ny,8)
         dz =  ( zlayer(ilayer+1) - zlayer(ilayer) ) / real(nz,8)

         allocate(top_surface(nx+1,ny+1), bottom_surface(nx+1,ny+1),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 1249')
         top_surface(:,:) =  zlayer(ilayer+1) - dz
         bottom_surface(:,:) = zlayer(ilayer)

         call mesh_regular_domain_Hex8(xmin_chunk, ymin_chunk, z, dx, dy, dz, nx, ny, nz-1, &
              top_surface, bottom_surface, ispec, ipoint, ilayer, nb_doubling+1)


         deallocate(top_surface)
         allocate(top_surface(2*nx+1,2*ny+1),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 1250')
         top_surface(:,:) =  zlayer(ilayer+1)
         bottom_surface(:,:) = zlayer(ilayer+1) - dz
         call mesh_doubling_domain_Hex8(xmin_chunk, ymin_chunk, z, dx, dy, dz, nx, ny, &
              top_surface, bottom_surface, ispec, ipoint)
         deallocate(top_surface, bottom_surface)

      enddo

       !! mesh last layer with no doubling fore last layer
      ilayer = nb_doubling + 1
      nx=nx_ref*(2**(ilayer-1))
      ny=ny_ref*(2**(ilayer-1))
      nz = nzlayer(ilayer)
      z = zlayer(ilayer)
      dx =  ( xmax_chunk - xmin_chunk ) / real(nx,8)
      dy =  ( ymax_chunk - ymin_chunk ) / real(ny,8)
      dz =  ( zlayer(ilayer+1) - zlayer(ilayer) ) / real(nz,8)
      allocate(top_surface(nx+1,ny+1), bottom_surface(nx+1,ny+1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1251')
      top_surface(:,:) =  zlayer(ilayer+1) - dz
      bottom_surface(:,:) = zlayer(ilayer)
      call mesh_regular_domain_Hex8(xmin_chunk, ymin_chunk, z, dx, dy, dz, nx, ny, nz, &
           top_surface, bottom_surface, ispec, ipoint, ilayer, nb_doubling + 1)
      deallocate(top_surface, bottom_surface)


      write(*,*) " used point for defining  mesh :", ipoint, "allocated point :", npoint_tot
      write(*,*) " used element for defining mesh :", ispec, " allocated elements :", nspec_tot


    end subroutine mesh_metric_chunk

!!##################################################################################################################################
!!                                  DEFINE FINAL MESH AND CONNECTIVITY
!!##################################################################################################################################

    subroutine  create_mesh_metric_chunk()

      integer, dimension(:), allocatable :: iglob, locval
      logical, dimension(:), allocatable :: ifseg
      integer                            :: nglob
      integer                            :: i, k, idom, ier
      character(len=10)                  :: MESH
      character(len=MAX_STRING_LEN)      :: filename

      MESH='./MESH/' !! VM VM harcoded directory (todo fix it)

      allocate(iglob(npoint_tot), locval(npoint_tot),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1252')
      allocate(ifseg(npoint_tot),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1253')
      call get_global(npoint_tot, x_mesh_point, y_mesh_point, z_mesh_point, iglob, locval, ifseg, nglob, xmin_chunk, xmax_chunk)


      !! write coords
      k=0
      filename = trim(MESH)//'nodes_coords_file'
      open(27,file=trim(filename))
      write(27,*) nglob
      write(*,*) " remaining points in mesh : ", nglob, " total point used in mesh building ", npoint_tot
      do i = 1, npoint_tot
         if (ifseg(i)) then
            k=k+1
            write(27,'(i14,3x,3(f20.5,1x))') k, x_mesh_point(i), y_mesh_point(i), z_mesh_point(i)
         endif
      enddo
      close(27)

      !! read again the point to have them in right order
      npoint_tot=k
      write(*,*) " number of point found ", npoint_tot, nglob
      deallocate(x_mesh_point, y_mesh_point, z_mesh_point)
      allocate(x_mesh_point(npoint_tot), y_mesh_point(npoint_tot), z_mesh_point(npoint_tot),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1254')

      filename = trim(MESH)//'nodes_coords_file'
      open(27,file=trim(filename))
      read(27,*) k
      do i = 1, npoint_tot
         read(27,*) k,  x_mesh_point(i), y_mesh_point(i), z_mesh_point(i)
      enddo
      close(27)

      !! write elements
      filename = trim(MESH)//'mesh_file'
      open(27,file=trim(filename))
      write(27,*) nspec_tot
      do i =1, nspec_tot
         write(27,'(9i15)') i, &
              iglob((EtoV(1, i))), iglob((EtoV(2, i))), iglob((EtoV(3, i))), iglob((EtoV(4, i))), &
              iglob((EtoV(5, i))), iglob((EtoV(6, i))), iglob((EtoV(7, i))), iglob((EtoV(8, i)))
      enddo
      close(27)

      !! write boundaries
      open(27,file=trim(MESH)//'absorbing_surface_file_xmin'); write(27,*) nspec_xmin
      open(28,file=trim(MESH)//'absorbing_surface_file_xmax'); write(28,*) nspec_xmin
      open(29,file=trim(MESH)//'absorbing_surface_file_ymin'); write(29,*) nspec_ymax
      open(30,file=trim(MESH)//'absorbing_surface_file_ymax'); write(30,*) nspec_ymax
      open(31,file=trim(MESH)//'absorbing_surface_file_bottom'); write(31,*) nspec_zmin
      open(32,file=trim(MESH)//'free_or_absorbing_surface_file_zmax'); write(32,*) nspec_zmax

      do i=1, nspec_tot
         if (iboun(1,i))  write(27,'(10(i10,1x))') i, &
              iglob((EtoV(1, i))), iglob((EtoV(2, i))), iglob((EtoV(6, i))), iglob((EtoV(5, i)))
         if (iboun(2,i))  write(28,'(10(i10,1x))') i, &
              iglob((EtoV(3, i))), iglob((EtoV(7, i))), iglob((EtoV(8, i))), iglob((EtoV(4, i)))
         if (iboun(3,i))  write(29,'(10(i10,1x))') i, &
              iglob((EtoV(1, i))), iglob((EtoV(4, i))), iglob((EtoV(8, i))), iglob((EtoV(5, i)))
         if (iboun(4,i))  write(30,'(10(i10,1x))') i, &
              iglob((EtoV(2, i))), iglob((EtoV(3, i))), iglob((EtoV(7, i))), iglob((EtoV(6, i)))
         if (iboun(5,i))  write(31,'(10(i10,1x))') i, &
              iglob((EtoV(1, i))), iglob((EtoV(2, i))), iglob((EtoV(3, i))), iglob((EtoV(4, i)))
         if (iboun(6,i))  write(32,'(10(i10,1x))') i, &
              iglob((EtoV(5, i))), iglob((EtoV(6, i))), iglob((EtoV(7, i))), iglob((EtoV(8, i)))
      enddo
      close(27)
      close(28)
      close(29)
      close(30)
      close(31)
      close(32)


      !! write materials !!
      filename = trim(MESH)//'nummaterial_velocity_file'
      open(27,file=trim(filename))
      do i=1, nb_mat
         write(27,'(2i6,5f15.5,i6)') 2, Imaterial_domain(i), material_prop(i,1:5), 0
      enddo
      close(27)

      allocate(Imatetrial_ispec(nspec_tot),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1255')

      filename = trim(MESH)//'/materials_file'
      open(28,file=trim(filename))
      do i=1, nspec_tot
         call Find_Domain(idom, i, iglob)
         write(28,*) i, Imaterial_domain(idom)
         Imatetrial_ispec(i)= Imaterial_domain(idom)
      enddo
      close(28)

      !! vtk mesh visualization file ---------------------------------------
      filename = trim(MESH)//'mesh_debug.vtk'
      call write_VTK_data_elem_i_earthmesh(nspec_tot,npoint_tot,x_mesh_point,y_mesh_point,z_mesh_point, &
                                           iglob,EtoV,Imatetrial_ispec,filename)

    end subroutine create_mesh_metric_chunk

!!##################################################################################################################################
!!                                  MESH LAYER WITH REGULAR ELEMENTS
!!##################################################################################################################################

    subroutine mesh_regular_domain_Hex8(ox, oy, oz,  dx, dy, dz, nx,  ny, nz, top_surface, bottom_surface, &
         ispec, ipoint, ilayer, nlayer)

      double precision,                                intent(in)    :: ox, oy, oz, dx, dy, dz
      integer,                                         intent(in)    :: nx, ny, nz
      double precision, dimension(:,:),   allocatable, intent(in)    :: top_surface, bottom_surface
      integer,                                         intent(in)    :: ilayer, nlayer
      integer,                                         intent(inout) :: ispec, ipoint

      !! locals
      double precision, dimension(:), allocatable                    :: xgrid, ygrid, zgrid
      double precision                                               :: ztop, zbottom
      integer                                                        :: i, j, k, ip, ier

      !!-------------------------------------------------------------------------------------------------------
      !! in regular layer domain we have :
      !!
      !!    (nx_elem_domain+1)*(ny_elem_domain+1)*(nx_elem_domain+1) point in regular grid
      !!    that corresponds to 8 corner of elememts
      !!
      !!------------------------------------------------------------------------------------------------------

      !! creating regular Cartesian grid ---------------------------------------------------------------------
      allocate(xgrid(nx+1), ygrid(ny+1), zgrid(nz+1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1256')

      !! define the reference  grid
      do i = 1, nx + 1
         xgrid(i) = ox + (i - 1) * dx
      enddo

      do i = 1, ny + 1
         ygrid(i) = oy + (i - 1) * dy
      enddo

      do i = 1, nz + 1
         zgrid(i) = oz + (i - 1) * dz
      enddo

      !! compute vertically deformed Cartesian grid --------------------------------------------------------------
      ip = ipoint
      do k = 1, nz + 1  !!loop over depth from bottom surface (k=1) to the top (k=nz+1)

         !! loop over (x,y) points
         do j = 1, ny + 1
            do i = 1, nx + 1

               !! next point in mesh
               ip = ip + 1

               !! get current top and bottom of mesh for the given (x,y)
               ztop = top_surface(i,j)
               zbottom  = bottom_surface(i,j)

               !! store current (x,y) mesh point
               x_mesh_point(ip) = xgrid(i)
               y_mesh_point(ip) = ygrid(j)

               !! vertically deformed mesh ::
               !! store z : linear interpolation between top and bottom to get the current position
               z_mesh_point(ip) = zbottom + ( zgrid(k) - zgrid(1) ) * ( ztop - zbottom ) / ( zgrid(nz+1) - zgrid(1) )

            enddo
         enddo

      enddo

      !! compute elemnts conectivity based on previous deformed Cartesian grid ----------------------------------
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx

               !! add a new element
               ispec = ispec + 1

               !! add 8 corners connectivity -------
               EtoV(1,ispec) =  i     +   (nx+1)*  ( j-1 ) + (nx+1)*(ny+1)*  ( k-1 ) + ipoint
               EtoV(2,ispec) =  i+1   +   (nx+1)*  ( j-1 ) + (nx+1)*(ny+1)*  ( k-1 ) + ipoint
               EtoV(3,ispec) =  i+1   +   (nx+1)*  ( j   ) + (nx+1)*(ny+1)*  ( k-1 ) + ipoint
               EtoV(4,ispec) =  i     +   (nx+1)*  ( j   ) + (nx+1)*(ny+1)*  ( k-1 ) + ipoint
               EtoV(5,ispec) =  i     +   (nx+1)*  ( j-1 ) + (nx+1)*(ny+1)*  ( k   ) + ipoint
               EtoV(6,ispec) =  i+1   +   (nx+1)*  ( j-1 ) + (nx+1)*(ny+1)*  ( k   ) + ipoint
               EtoV(7,ispec) =  i+1   +   (nx+1)*  ( j   ) + (nx+1)*(ny+1)*  ( k   ) + ipoint
               EtoV(8,ispec) =  i     +   (nx+1)*  ( j   ) + (nx+1)*(ny+1)*  ( k   ) + ipoint

               !! get boundaries --------
               if (i == 1) then
                  iboun(1,ispec)=.true.
                  nspec_xmin=nspec_xmin+1
               endif

               if (i == nx) then
                  iboun(2,ispec)=.true.
                  nspec_xmax=nspec_xmax+1
               endif

               if (j == 1) then
                  iboun(3,ispec)=.true.
                  nspec_ymin=nspec_ymin+1
               endif

               if (j == ny) then
                  iboun(4,ispec)=.true.
                  nspec_ymax=nspec_ymax+1
               endif

               if (ilayer == 1 .and. k == 1) then
                  iboun(5,ispec)=.true.
                  nspec_zmin=nspec_zmin+1
               endif

               if (ilayer == nlayer .and. k == nz) then
                  iboun(6,ispec)=.true.
                  nspec_zmax=nspec_zmax+1
               endif


            enddo
         enddo
      enddo

      write(*,*) 'add ', (nx+1)*(ny+1)*(nz+1), ipoint
      ipoint  = ipoint + (nx+1)*(ny+1)*(nz+1)

      deallocate(xgrid, ygrid, zgrid)

    end subroutine mesh_regular_domain_Hex8

!!##################################################################################################################################
!!                           MESH LAYER WITH DOUBLING ELEMENTS
!!##################################################################################################################################

    subroutine mesh_doubling_domain_Hex8(ox, oy, z, dx, dy, dz, nx,  ny, top_surface, bottom_surface, ispec, ipoint)

      double precision,                              intent(in)    :: ox, oy, z, dx, dy, dz
      double precision, dimension(:,:), allocatable, intent(in)    :: top_surface, bottom_surface
      integer,                                       intent(in)    :: nx, ny
      integer,                                       intent(inout) :: ispec, ipoint

      !! locals
      integer :: i, j, tupe
      double precision :: x, y
      !! in this domain we have  :
      !!                (nx+1)*(ny+1)*(nz+1) point on  bottom grid  size (dx,dy,dz)
      !! and            (2*nx+1)*(2*ny+1)*(2*nz+1) point on top grid size (0.5*dx, 0.5*dy, 0.5*dz)
      !!


      !! define the elements in flat doubling layer
      !! loop over orizontal layer
      do j = 1, ny
         y = ox + (j-1)*dy
         do i = 1, nx
            x = ox + (i-1)*dx

            !* Check doubling brick type and assign elem type
            if (modulo(i,2) /= 0 .and. modulo(j,2) /= 0) then ! Type A
               tupe = 10
            else if (modulo(i,2) == 0 .and. modulo(j,2) /= 0) then ! Type B
               tupe = 20
            else if (modulo(i,2) == 0 .and. modulo(j,2) == 0) then ! Type C
               tupe = 30
            else if (modulo(i,2) /= 0 .and. modulo(j,2) == 0) then ! Type D
               tupe = 40
            endif

            call super_brick_1_Hex8(tupe, dx, dy, dz, x, y, z, ox, oy,top_surface, bottom_surface, ispec, ipoint)


            !! get boundaries --------
            if (i == 1) then
               iboun(1,ispec)=.true.
               nspec_xmin=nspec_xmin+1
            endif

            if (i == nx) then
               iboun(2,ispec)=.true.
               nspec_xmax=nspec_xmax+1
            endif

            if (j == 1) then
               iboun(3,ispec)=.true.
               nspec_ymin=nspec_ymin+1
            endif

            if (j == ny) then
               iboun(4,ispec)=.true.
               nspec_ymax=nspec_ymax+1
            endif

         enddo
      enddo

    end subroutine mesh_doubling_domain_Hex8


!!##################################################################################################################################
!!                                        SUPER-BRICK 1 FOR DOUBLING
!!##################################################################################################################################

    subroutine super_brick_1_Hex8(tupe, dx, dy, dz, x, y, z, ox, oy, top_surface, bottom_surface, ispec, ipoint)

      integer,                                       intent(in)    :: tupe
      double precision,                              intent(in)    :: dx, dy, dz, x, y, z, ox, oy
      double precision, dimension(:,:), allocatable, intent(in)    :: top_surface, bottom_surface
      integer,                                       intent(inout) :: ispec, ipoint

      !! locals
      integer :: icorner, itop, jtop, ibottom, jbottom, ibri, type
      double precision :: elem(8,3)
      double precision :: ztop, zbottom

      do ibri = 1, 7

         !! next element
         ispec = ispec + 1

         !!code to know which kind of element we are using
         type = tupe + ibri

         !! 1/ get normalized reference element
         elem=elemref(type)

         !! 2/ scale and shift the element
         elem(:,1) = elem(:,1) * dx + x
         elem(:,2) = elem(:,2) * dy + y
         elem(:,3) = elem(:,3) * dz + z


         !! add 8 corner of the curreznt element in mesh
         do icorner = 1, 8
            ipoint = ipoint + 1

            !! compute verticaly deformed element
            itop = int( (elem(icorner,1) - ox ) / (0.5d0*dx) ) + 1
            jtop = int( (elem(icorner,2) - oy ) / (0.5d0*dy) ) + 1

            ibottom =  int( (elem(icorner,1) - ox ) / dx ) + 1
            jbottom =  int( (elem(icorner,2) - oy ) / dy ) + 1

            ztop = top_surface(itop, jtop)
            zbottom = bottom_surface(ibottom, jbottom)

            !! store horizontal coordinate
            x_mesh_point(ipoint) = elem(icorner,1)
            y_mesh_point(ipoint) = elem(icorner,2)

            !! lenear interpolation for vertical deformation
            z_mesh_point(ipoint) = zbottom + ( elem(icorner,3) - z ) *  ( ztop - zbottom ) / dz

            !! connectivity
            EtoV(icorner,ispec) = ipoint
         enddo

      enddo

    end subroutine super_brick_1_Hex8

!! =============================== define 8 corner of one element of  normalized super-brick =======================================

   function elemref(type) result(elem)

     integer          :: type

     !! locals
     double precision :: zero, one, half, onethird, twothird
     double precision :: myref(8,3), ref(8,3), elem(8,3)


     zero = 0.d0
     half = 0.5d0
     onethird = 1.d0 / 3.d0
     twothird = 2.d0 / 3.d0
     one =1.d0

     if (mod(type,10) == 1) then

        myref(1,:) = (/ zero, zero, zero /)
        myref(2,:) = (/ half, zero, onethird /)
        myref(3,:) = (/ half, half, twothird /)
        myref(4,:) = (/ zero, half, onethird /)
        myref(5,:) = (/ zero, zero, one /)
        myref(6,:) = (/ half, zero, one /)
        myref(7,:) = (/ half, half, one /)
        myref(8,:) = (/ zero, half, one /)

     else if (mod(type,10) == 2) then

        myref(1,:) = (/ zero, half, onethird /)
        myref(2,:) = (/ half, half, twothird /)
        myref(3,:) = (/ half, one, twothird /)
        myref(4,:) = (/ zero, one, onethird /)
        myref(5,:) = (/ zero, half, one /)
        myref(6,:) = (/ half, half, one /)
        myref(7,:) = (/ half, one, one /)
        myref(8,:) = (/ zero, one, one /)

     else if (mod(type,10) == 3) then

        myref(1,:) = (/ zero, zero, zero /)
        myref(2,:) = (/ half, zero, onethird /)
        myref(3,:) = (/ half, one, onethird /)
        myref(4,:) = (/ zero, one, zero /)
        myref(5,:) = (/ zero, half, onethird /)
        myref(6,:) = (/ half, half, twothird /)
        myref(7,:) = (/ half, one, twothird /)
        myref(8,:) = (/ zero, one, onethird /)

     else if (mod(type,10) == 4) then

        myref(1,:) = (/ zero, zero, zero /)
        myref(2,:) = (/ one, zero, zero /)
        myref(3,:) = (/ one, one, zero /)
        myref(4,:) = (/ zero, one, zero /)
        myref(5,:) = (/ half, zero, onethird /)
        myref(6,:) = (/ one, zero, onethird /)
        myref(7,:) = (/ one, one, onethird /)
        myref(8,:) = (/ half, one, onethird /)

     else if (mod(type,10) == 5) then

        myref(1,:) = (/ half, zero, onethird /)
        myref(2,:) = (/ one, zero, onethird /)
        myref(3,:) = (/ one, one, onethird /)
        myref(4,:) = (/ half, one, onethird /)
        myref(5,:) = (/ half, half, twothird /)
        myref(6,:) = (/ one, half, twothird /)
        myref(7,:) = (/ one, one, twothird /)
        myref(8,:) = (/ half, one, twothird /)

     else if (mod(type,10) == 6) then

        myref(1,:) = (/ half, zero, onethird /)
        myref(2,:) = (/ one, zero, onethird /)
        myref(3,:) = (/ one, half, twothird /)
        myref(4,:) = (/ half, half, twothird /)
        myref(5,:) = (/ half, zero, one /)
        myref(6,:) = (/ one, zero, one /)
        myref(7,:) = (/ one, half, one /)
        myref(8,:) = (/ half, half, one /)

     else if (mod(type,10) == 7) then

        myref(1,:) = (/ half, half, twothird /)
        myref(2,:) = (/ one, half, twothird /)
        myref(3,:) = (/ one, one, twothird /)
        myref(4,:) = (/ half, one, twothird /)
        myref(5,:) = (/ half, half, one /)
        myref(6,:) = (/ one, half, one /)
        myref(7,:) = (/ one, one, one /)
        myref(8,:) = (/ half, one, one /)

     endif


     ! Check doublingbrick type
     if (type/10 == 1) then ! type A

        ref = myref

     else if (type/10 == 2) then ! type B

        ref(1,:) = myref(2,:)
        ref(2,:) = myref(1,:)
        ref(3,:) = myref(4,:)
        ref(4,:) = myref(3,:)
        ref(5,:) = myref(6,:)
        ref(6,:) = myref(5,:)
        ref(7,:) = myref(8,:)
        ref(8,:) = myref(7,:)

        ref(:,1) = one-ref(:,1)

     else if (type/10 == 3) then ! type C

        ref(1,:) = myref(3,:)
        ref(2,:) = myref(4,:)
        ref(3,:) = myref(1,:)
        ref(4,:) = myref(2,:)
        ref(5,:) = myref(7,:)
        ref(6,:) = myref(8,:)
        ref(7,:) = myref(5,:)
        ref(8,:) = myref(6,:)

        ref(:,1) = one-ref(:,1)
        ref(:,2) = one-ref(:,2)

     else if (type/10 == 4) then ! type D

        ref(1,:) = myref(4,:)
        ref(2,:) = myref(3,:)
        ref(3,:) = myref(2,:)
        ref(4,:) = myref(1,:)
        ref(5,:) = myref(8,:)
        ref(6,:) = myref(7,:)
        ref(7,:) = myref(6,:)
        ref(8,:) = myref(5,:)

        ref(:,2) = one-ref(:,2)

     endif

     elem(:,:) = ref(:,:)

   end function elemref

!!===================================================================================================================

   subroutine Find_Domain(id, i, iglob)
     integer, dimension(:), allocatable :: iglob
     integer                            :: id, i
     integer                            :: k
     integer                            :: ind(8)
     double precision                   :: x, y, z

     do k=1, 8
        ind(k)=iglob(EtoV(k,i));
     enddo

     x = 0.125 * sum(x_mesh_point(ind(:)))
     y = 0.125 * sum(y_mesh_point(ind(:)))
     z = 0.125 * sum(z_mesh_point(ind(:)))

     id = -1
     do k = 1,  nb_dom
        if ( x >= domain_boundary(k,1) .and. x <= domain_boundary(k,2) .and. &
             y >= domain_boundary(k,3) .and. y <= domain_boundary(k,4) .and. &
             z >= domain_boundary(k,5) .and. z <= domain_boundary(k,6)) then
           id=k
           return
        endif
     enddo

     if (id < 0) then
        write(*,*) " Error : domain not,  found setting 1 "
        id  = 1
     endif
   end subroutine Find_Domain

!!===================================================================================================================

end module chunk_earth_mod

!!$!! vtk debug file
!!$      open(IOUT_VTK,file='raw_mesh_debug.vtk',status='unknown')
!!$      write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
!!$      write(IOUT_VTK,'(a)') 'material model VTK file'
!!$      write(IOUT_VTK,'(a)') 'ASCII'
!!$      write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
!!$      write(IOUT_VTK, '(a,i12,a)') 'POINTS ', npoint_tot, ' float'
!!$      do i=1,npoint_tot
!!$         write(IOUT_VTK,*) x_mesh_point(i), y_mesh_point(i), z_mesh_point(i)
!!$      enddo
!!$      write(IOUT_VTK,*)
!!$      write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec_tot,nspec_tot*9
!!$      do i =1, nspec_tot
!!$         write(27,'(9i15)') 8, &
!!$              (EtoV(1, i))-1, (EtoV(2, i))-1,  (EtoV(3, i))-1, (EtoV(4, i))-1, &
!!$              (EtoV(5, i))-1, (EtoV(6, i))-1,  (EtoV(7, i))-1, (EtoV(8, i))-1
!!$      enddo
!!$
!!$      write(IOUT_VTK,*)
!!$      write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec_tot
!!$      write(IOUT_VTK,'(6i12)') (12,i=1,nspec_tot)
!!$      write(IOUT_VTK,*)
!!$      write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec_tot
!!$      write(IOUT_VTK,'(a)') "SCALARS elem_val float"
!!$      write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
!!$      do i =1, nspec_tot
!!$          write(IOUT_VTK,*) 1
!!$      enddo
!!$      close(IOUT_VTK)
!!$
!!$      !! write coords
!!$      open(27,file=trim(MESH)//'RAW_nodes_coords_file')
!!$      write(27,*) npoint_tot
!!$      do i = 1, npoint_tot
!!$          write(27,'(i14,3x,3(f20.5,1x))') i, x_mesh_point(i), y_mesh_point(i), z_mesh_point(i)
!!$      enddo
!!$      close(27)
!!$
!!$      !! write elements
!!$      open(27,file=trim(MESH)//'RAW_mesh_file')
!!$      write(27,*) nspec_tot
!!$      do i =1, nspec_tot
!!$         write(27,'(9i15)') i, &
!!$              (EtoV(1, i)), (EtoV(2, i)), (EtoV(3, i)), (EtoV(4, i)), &
!!$              (EtoV(5, i)), (EtoV(6, i)), (EtoV(7, i)), (EtoV(8, i))
!!$      enddo
!!$      close(27)
