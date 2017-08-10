  program create_topo

    implicit none
    character(len=3) code
   character(len=250) par_file
   real dx,dy,z,xmin,ymin,xi,yi
   real Ampl,length,phase,extx,exty
   integer nx,ny

   real xc,yc,zc
   integer i,j,k

   real, parameter :: TwoPi=6.28318530718

   write(*,*) ' Parameter file name  '
   write(*,*) ' following format :'
   write(*,*)
   write(*,*) ' nx ny       !: (integer) number of nodes in x and y direction (m)'
   write(*,*) ' dx dy       !: (real) grid spacing in x and y direction (m)'
   write(*,*) ' xmin, ymin  !: (real) coordinate of the first point (m)'
   write(*,*) ' code        !: (character(len=3)) type of topo '
   write(*,*)
   write(*,*) '    if code=fla  '
   write(*,*) ' z           !: (real) depth of flat topo '
   write(*,*)
   write(*,*) '    if code=gau  '
   write(*,*) ' z           !: (real) depth ref of topo '
   write(*,*) ' xi, yi'
   write(*,*) '  Ampl, extx, exty'
   write(*,*)
   write(*,*) '    if code=cos  '
   write(*,*) ' z           !: (real) depth ref of topo '
   write(*,*) '  length, phase '
   write(*,*) '  Ampl'
   write(*,*)
   write(*,*) ' give Parameter file name ? '

   read(*,'(a)') par_file

   open(10,file=trim(par_file))
   read(10,*) nx,ny
   read(10,*) dx,dy
   read(10,*) xmin,ymin
   read(10,'(a3)') code

   if (code == 'fla') then

      read(10,*) z

   else if (code == 'gau') then
      read(10,*) z
      read(10,*) xi,yi
      read(10,*) Ampl,extx,exty

   else if (code == 'cos') then
      read(10,*) z
      read(10,*) length, phase
      read(10,*) Ampl

   else

      write(*,*) 'interf not implemented :', code
      stop

   endif

   close(10)

   open(10,file='interf.txt')

   select case (code)

   case ('fla')

      zc=z
      yc = ymin - dy
      do j = 1, ny
         yc = yc + dy
         xc = xmin - dx
         do i = 1, nx
            k =  k + 1
            xc = xc + dx
            write(10,*) xc,yc,zc
         enddo
      enddo

   case('gau')
      yc = ymin - dy
      do j = 1, ny
         yc = yc + dy
         xc = xmin - dx
         do i = 1, nx
            k =  k + 1
            xc = xc + dx
            zc =  z + Ampl* exp( -0.5 * ((xc-xi)/extx)**2 - 0.5 * ((yc-yi)/exty)**2  )
            write(10,*) xc,yc,zc
         enddo
      enddo

   case ('cos')
      yc = ymin - dy
      do j = 1, ny
         yc = yc + dy
         xc = xmin - dx
         do i = 1, nx
            k =  k + 1
            xc = xc + dx
            zc = z + Ampl* cos(xc * TwoPI / length  + phase)
            write(10,*) xc,yc,zc
         enddo
      enddo

   end select

   close(10)
  end program create_topo
