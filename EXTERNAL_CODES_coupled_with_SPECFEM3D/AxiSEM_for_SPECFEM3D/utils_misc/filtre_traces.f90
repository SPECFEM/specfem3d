  program filtre_traces

    implicit none
    character(len=250) arg,indat,insyn,fichier
    integer npts,irek,norder
    double precision, allocatable :: f1(:),f2(:),trdat(:,:)
    double precision fl,fh,dt

    irek=1
    norder=4

!!$    call get_command_argument(1,arg)
!!$    indat=trim(arg)
!!$    call get_command_argument(2,arg)
!!$    insyn=trim(arg)
    call get_command_argument(1,arg)
    read(arg,*) fl
    call get_command_argument(2,arg)
    read(arg,*) fh
    indat='wksis.txt'
    call read_npts(indat,npts)
    allocate(trdat(2,npts),f1(npts),f2(npts))
    call read_trace(indat,trdat,npts,dt)
    write(*,*) npts
    f1(:)=trdat(2,:)
    f2(:)=trdat(2,:)
    call  bwfilt (f1, f2, dt, npts, irek, norder, fl, fh)
    write(*,*) npts
    fichier='dat_tmp.txt'
    call write_adj(f2,trdat,fichier,npts)
!!$
!!$    call read_trace(insyn,trdat,npts,dt)
!!$    f1(:)=trdat(2,:)
!!$    f2(:)=trdat(2,:)
!!$    call  bwfilt (f1, f2, dt, npts, irek, norder, fl, fh)
!!$    fichier='syn_tmp.txt'
!!$    call write_adj(f2,trdat,fichier,npts)

  end program filtre_traces

  subroutine write_adj(adj,synth,fichier,n)
    implicit none
    integer n,i
    double precision adj(n),synth(2,n)
    character(len=250) fichier
    open(30,file=trim(fichier))
    do i=1,n
       if (mod(i,100) == 0) write(*,*) i,n
       write(30,*) synth(1,i), adj(i)
    enddo
    close(30)
  end subroutine write_adj

  subroutine read_npts(fichier,nstep)
    implicit none
    integer nstep
    character(len=250) fichier
    real x,y

    nstep=0
    open(20,file=trim(fichier))
    do
       read(20,*,end=99) x,y
       nstep=nstep + 1
    enddo
99  continue
    close(20)

  end subroutine read_npts

  subroutine read_trace(fichier,trace,nstep,dt)
    implicit none
    integer nstep,i
    double precision dt,trace(2,nstep)
    character(len=250) fichier

    !write(*,*) trim(fichier)


    open(20,file=trim(fichier))

    trace(:,:) = 0.d0
    i = 1
    do i=1,nstep
       read(20,*,end=99) trace(1,i),trace(2,i)
    enddo
99  continue
    dt = trace(1,2) - trace(1,1)

    close(20)
  end subroutine read_trace


