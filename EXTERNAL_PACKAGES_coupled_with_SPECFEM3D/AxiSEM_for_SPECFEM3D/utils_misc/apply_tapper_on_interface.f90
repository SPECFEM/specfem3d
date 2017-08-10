  program apply_tapper_on_interface

    implicit none

    real                              :: ox, oy
    real                              :: dx, dy
    integer                           :: nx, ny
    character(len=100)                :: interface_file

    real                              :: xmin, xmax
    real                              :: ymin, ymax
    real                              :: zedge
    real                              :: lx, ly

    real, dimension(:,:), allocatable :: topo, topo_tap
    real, parameter                   :: lsmooth=0.1
    real, parameter                   :: epsilon=0.01

    integer                           :: i,j

    read(*,*)     ox, oy
    read(*,*)     nx, ny
    read(*,*)     dx, dy
    read(*,'(a)') interface_file
    read(*,*)     xmin, xmax
    read(*,*)     ymin, ymax
    read(*,*)     zedge

    lx = xmax -xmin
    ly = ymax -ymin
    !!
    xmin = xmin + epsilon*lx
    xmax = xmax - epsilon*ly
    ymin = ymin + epsilon*ly
    ymax = ymax - epsilon*ly

    allocate(topo(nx,ny), topo_tap(nx,ny))

    open(10, file=trim(interface_file))
    do j=1,ny
       do i=1,nx
          read(10,*) topo(i,j)
       enddo
    enddo
    close(10)

    topo_tap(:,:)=topo(:,:)


    call smooth_topo(topo_tap, topo, ox, oy, dx, dy, nx, ny, xmin, xmax, ymin, ymax, zedge, lsmooth)

    interface_file='interf_tapper.dat'
    open(10, file=trim(interface_file))
    do j=1,ny
       do i=1,nx
          write(10,*) topo_tap(i,j)
       enddo
    enddo
    close(10)


    deallocate(topo, topo_tap)

  end program apply_tapper_on_interface


  subroutine smooth_topo(topo_tap, topo, ox, oy, dx, dy, nx, ny, xmin, xmax, ymin, ymax, zedge, lsmooth)

    implicit none
    real,                                intent(in)    :: ox, oy
    real,                                intent(in)    :: dx, dy
    integer,                             intent(in)    :: nx, ny
    real,                                intent(in)    :: xmin, xmax
    real,                                intent(in)    :: ymin, ymax
    real,                                intent(in)    :: zedge
    real,                                intent(in)    :: lsmooth
    real, dimension(nx,ny), intent(in)    :: topo
    real, dimension(nx,ny), intent(inout) :: topo_tap

    real                                             :: x, y
    integer                                          :: i, j
    real                                             :: xmin_inner, xmax_inner
    real                                             :: ymin_inner, ymax_inner
    real                                             :: tapperx, tappery, dz
    real                                             :: tapper_up

    xmin_inner = xmin + lsmooth *abs(xmax - xmin)
    xmax_inner = xmax - lsmooth *abs(xmax - xmin)
    ymin_inner = ymin + lsmooth *abs(ymax - ymin)
    ymax_inner = ymax - lsmooth *abs(ymax - ymin)



    topo_tap=0.

    do j = 1, ny

       y = oy + (j-1) * dy

       do i = 1, nx

          x = ox + (i-1) * dx
          dz = topo(i,j) - zedge

          !! outside domain
          if ( x <= xmin .or. x >= xmax .or. y <= ymin .or. y >= ymax ) then
             topo_tap(i,j) = zedge
          else

             !! outside inner domain
             if ( x < xmin_inner .or. x > xmax_inner .or. y < ymin_inner .or. y > ymax_inner ) then

                if ( x < xmin_inner) then
                    tapperx=tapper_up(xmin, xmin_inner, x)
                   !1
                   if ( y < ymin_inner) then
                      tappery=tapper_up(ymin, ymin_inner, y)
                      topo_tap(i,j) = zedge + 0.25*(1+tapperx)*(1+tappery)*dz
                   endif

                   !2
                   if (y > ymax_inner) then
                      tappery=tapper_up(ymax, ymax_inner, y)
                      topo_tap(i,j) =  zedge + 0.25*(1+tapperx)*(1+tappery)*dz
                   endif

                   !3
                   if ( y >= ymin_inner .and. y <= ymax_inner) then
                      topo_tap(i,j) = zedge + 0.5*(1+tapperx)*dz
                   endif


                endif

                if ( x > xmax_inner) then
                   tapperx=tapper_up(xmax, xmax_inner, x)

                   !4
                   if ( y < ymin_inner) then
                      tappery=tapper_up(ymin, ymin_inner, y)
                      topo_tap(i,j) = zedge + 0.25*(1+tapperx)*(1+tappery)*dz
                   endif

                   !5
                   if (y > ymax_inner) then
                      tappery=tapper_up(ymax, ymax_inner, y)
                      topo_tap(i,j) = zedge + 0.25*(1+tapperx)*(1+tappery)*dz
                   endif

                   !6
                   if ( y >= ymin_inner .and. y <= ymax_inner) then
                      topo_tap(i,j) = zedge + 0.5*(1+tapperx)*dz
                   endif

                endif

                if ( x >= xmin_inner .and. x <= xmax_inner) then

                   !7
                   if (y < ymin_inner) then
                       tappery=tapper_up(ymin, ymin_inner, y)
                   endif

                   !8
                   if (y > ymax_inner) then
                      tappery=tapper_up(ymax, ymax_inner, y)
                   endif

                   topo_tap(i,j) = zedge + 0.5*(1+tappery)*dz

                endif

             else !! inside inner domain
                topo_tap(i,j) = topo(i,j)
             endif

          endif

       enddo
    enddo




  end subroutine smooth_topo

  real function tapper_up(x0,x1,x)
    real, parameter ::  pi=3.141592653
    real            ::  omega, phi
    omega=pi/(x0-x1)
    phi=-omega*x1
    tapper_up=cos(omega*x+phi)
    return
  end function tapper_up
