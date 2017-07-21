program sem_model_slice

  implicit none

  include 'mpif.h'
  include 'constants.h'
  include "precision.h"
  include 'values_from_mesher.h'

  integer, parameter :: NMAXPTS = 100000
  integer :: ier,sizeprocs,myrank,ios, i,j, k, ispec,iglob,ipt, npts
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  character(len=150) :: xyz_infile,topo_dir,model_dir,data_name,gmt_outfile, &
             local_data_file, prname
  real(kind=CUSTOM_REAL),dimension(NMAXPTS) :: x, y, z, v, vmin, &
             vall,distmin, dist, distall
  integer,dimension(NMAXPTS) :: ispec_min,ix_min,iy_min,iz_min
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: vstore
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) ::  xstore,ystore,zstore
  real(kind=CUSTOM_REAL),dimension(2,NMAXPTS) :: in, out

  ! true --> replace "air" points with NaN (vertical cross sections)
  ! false --> take the closest value to the "air" points (horizontal cross section)
  logical, parameter :: TOPOGRAPHY = .true.

  character(len=100) :: topo_file
  integer, dimension(NX_TOPO_SOCAL,NY_TOPO_SOCAL) :: itopo_bathy_basin
  double precision :: elevation

  ! MPI initialization
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  ! input arguments
  call get_command_argument(1,xyz_infile)
  call get_command_argument(2,topo_dir)
  call get_command_argument(3,model_dir)
  call get_command_argument(4,data_name)
  call get_command_argument(5,gmt_outfile)

  ! read points to be interpolated
  open(11,file=xyz_infile,iostat=ios)
  i=0
  do while (1 == 1)
    i=i+1
    read(11,*,iostat=ios) x(i),y(i),z(i)
    if (ios /= 0) exit
  enddo
  close(11)
  npts=i-1
  if (myrank == 0) then
    if (npts > NMAXPTS .or. npts < 0) call exit_mpi(myrank,'Npts error ...')
    write(*,*) 'Total number of points = ', npts
  endif

  ! read data and topo files
  write(prname,'(a,i6.6,a)') trim(model_dir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(data_name)//'.bin'
  open(unit = 27,file=local_data_file,status='old', iostat = ios,form ='unformatted')
  if (ios /= 0) call exit_mpi(myrank, 'Error reading model data file')
  read(27) vstore
  close(27)

  write(prname,'(a,i6.6,a)') trim(topo_dir)//'/proc',myrank,'_'
  open(unit = 27,file = trim(prname)//'x.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) call exit_mpi(myrank,'Error reading local x file')
  read(27) xstore
  close(27)
  open(unit = 27,file = trim(prname)//'y.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) call exit_mpi(myrank,'Error reading local y file')
  read(27) ystore
  close(27)
  open(unit = 27,file = trim(prname)//'z.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) call exit_mpi(myrank,'Error reading local z file')
  read(27) zstore
  close(27)
  open(unit = 27,file = trim(prname)//'ibool.bin',status='old',action='read', iostat = ios,form ='unformatted')
  if (ios /= 0) call exit_mpi(myrank,'Error reading local z file')
  read(27) ibool
  close(27)


  ! search for local minimum-distance point
  distmin(1:npts) = HUGEVAL

  do ispec=1,NSPEC_AB

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dist(1:npts)=dsqrt((x(1:npts)-dble(xstore(iglob)))**2 &
                     +(y(1:npts)-dble(ystore(iglob)))**2 &
                     +(z(1:npts)-dble(zstore(iglob)))**2)
          do ipt=1,npts
            if (dist(ipt) < distmin(ipt)) then
              distmin(ipt)=dist(ipt)
              ispec_min(ipt)=ispec
              ix_min(ipt)=i; iy_min(ipt)=j; iz_min(ipt)=k
              vmin(ipt)=vstore(i,j,k,ispec)
            endif
          enddo

        enddo
      enddo
    enddo

    ! end of loop on all the elements in current slice
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (myrank == 0) print *, 'Done looping over global points ...'

  ! choose the minimum value

  !  write(prname,'(a,i6.6,a)') 'OUTPUT_FILES/in_',myrank,'.txt'
  !  open(33,file=prname)
  !  do i = 1, npts
  !    write(33,*) i, myrank, distmin(i), ispec_min(i), ix_min(i), iy_min(i), iz_min(i),vmin(i)
  !  enddo
  !  close(33)

  do i=1, npts
    in(1,i) = distmin(i)
    in(2,i) = myrank    ! myrank is coerced to a double
  enddo
  call MPI_REDUCE(in,out,npts,CUSTOM_MPI_2REAL,MPI_MINLOC,0,MPI_COMM_WORLD,ier)

  !  if (myrank == 0) then
  !   open(33,file='OUTPUT_FILES/out.txt')
  !   do i = 1, npts
  !     write(33,*) i, out(1,i), out(2,i)
  !  enddo
  !   close(33)
  !  endif

  call MPI_BCAST(out,2*npts,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  v(1:npts) = 0
  dist(1:npts) = 0.

  do i = 1, npts
    if (myrank == nint(out(2,i))) then
      v(i) = vmin(i)
!      if (GLL_INTERPOLATION) call xeg_search(x(i),y(i),z(i), &
!                 ispec_min(i),ix_min(i),iy_min(i),iz_min(i),v(i))
      dist(i) = distmin(i)
    endif
  enddo

  call MPI_REDUCE(v,vall,npts,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call MPI_REDUCE(dist,distall,npts,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  if (myrank == 0) then

    if (TOPOGRAPHY) then
      topo_file='/ibrixfs1/home/lqy/cmt3d/test_dir/'//trim(TOPO_FILE_SOCAL)
      call read_basin_topo_bathy_file(itopo_bathy_basin,NX_TOPO_SOCAL,NY_TOPO_SOCAL,topo_file)
    endif

    print *, 'Writing out gmt file ...'
    open(12,file=gmt_outfile,status='unknown')
    do i = 1, npts
      if (TOPOGRAPHY) then
         call topo_value(itopo_bathy_basin,dble(x(i)),dble(y(i)),elevation)
         if (elevation < z(i)) vall(i)= -1000.00
      endif
      write(12,*) x(i), y(i), z(i), vall(i), distall(i)
    enddo
    close(12)
  endif

  call MPI_FINALIZE(ier)

end program sem_model_slice

!------------------------------------------------

subroutine topo_value(itopo_bathy_basin,x,y,elevation)

  implicit none
  include 'constants.h'

  integer,dimension(NX_TOPO_SOCAL,NY_TOPO_SOCAL) :: itopo_bathy_basin
  double precision :: x, y
  double precision elevation

  double precision :: long, lat
  integer :: icornerlong, icornerlat
  double precision :: long_corner, lat_corner,ratio_xi,ratio_eta
  integer, parameter :: UTM_PROJECTION_ZONE  = 11
  logical, parameter :: SUPPRESS_UTM_PROJECTION  = .false.


  call utm_geo(long,lat,x,y,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

  ! get coordinate of corner in bathy/topo model
  icornerlong = int((long - ORIG_LONG_TOPO_SOCAL) / DEGREES_PER_CELL_TOPO_SOCAL) + 1
  icornerlat = int((lat - ORIG_LAT_TOPO_SOCAL) / DEGREES_PER_CELL_TOPO_SOCAL) + 1

  ! avoid edge effects and extend with identical point if outside model
  if (icornerlong < 1) icornerlong = 1
  if (icornerlong > NX_TOPO_SOCAL-1) icornerlong = NX_TOPO_SOCAL-1
  if (icornerlat < 1) icornerlat = 1
  if (icornerlat > NY_TOPO_SOCAL-1) icornerlat = NY_TOPO_SOCAL-1

  ! compute coordinates of corner
  long_corner = ORIG_LONG_TOPO_SOCAL + (icornerlong-1)*DEGREES_PER_CELL_TOPO_SOCAL
  lat_corner = ORIG_LAT_TOPO_SOCAL + (icornerlat-1)*DEGREES_PER_CELL_TOPO_SOCAL

  ! compute ratio for interpolation
  ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO_SOCAL
  ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO_SOCAL

  ! avoid edge effects
  if (ratio_xi < 0.) ratio_xi = 0.
  if (ratio_xi > 1.) ratio_xi = 1.
  if (ratio_eta < 0.) ratio_eta = 0.
  if (ratio_eta > 1.) ratio_eta = 1.

  ! interpolate elevation at current point
  elevation = &
             itopo_bathy_basin(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
             itopo_bathy_basin(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
             itopo_bathy_basin(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
             itopo_bathy_basin(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta


end subroutine topo_value
