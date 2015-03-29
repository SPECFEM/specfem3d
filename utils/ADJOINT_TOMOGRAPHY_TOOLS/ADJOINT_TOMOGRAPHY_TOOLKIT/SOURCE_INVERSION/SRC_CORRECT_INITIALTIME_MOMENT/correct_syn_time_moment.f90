program xcorrect_syn

implicit none

integer,parameter:: NDIM=80000
real:: dt0,dm0
integer:: nfile,i,j
character(len=250):: oldsyn_fnm,newsyn_fnm
real,dimension(NDIM):: oldsyn,newsyn,syn_t
integer:: npt,nerr,nerr1
real:: b,dt,b_new
real:: evla,evlo,stla,stlo,evdp
character*8:: kstnm,knetwk,kcmpnm
integer::ishift

write(*,*) 'read input file'
read(*,*) dt0, dm0
read(*,*) nfile

write(*,*) 'correct initial time:',dt0
write(*,*) 'correct scale moment:',dm0
write(*,*) 'total number of files:',nfile

do i = 1,nfile

        oldsyn=0.0
        newsyn=0.0
        syn_t=0.0

        read(*,'(a)') oldsyn_fnm
        read(*,'(a)') newsyn_fnm

        write(*,*) 'old syn:',oldsyn_fnm
        write(*,*) 'new syn:',newsyn_fnm

        call rsac1(oldsyn_fnm,oldsyn,npt,b,dt,NDIM,nerr)
        if (nerr/=0) stop 'error reading sac file'

        call getfhv('evla',evla,nerr)
        call getfhv('evlo',evlo,nerr)
        call getfhv('stla',stla,nerr)
        call getfhv('stlo',stlo,nerr)
        call getfhv('evdp',evdp,nerr)
        call getkhv('kstnm',kstnm,nerr)
        call getkhv('kcmpnm',kcmpnm,nerr)
        call getkhv('knetwk',knetwk,nerr)


        ! shift and apply scale moment correction for seismograms
        ishift=nint(abs(dt0)/dt)
        if (dt0 > 0.0 ) then
           newsyn(1+ishift:npt+ishift)=dm0*oldsyn(1:npt)
        else
           newsyn(1:npt)=dm0*oldsyn(1+ishift:npt+ishift)
        endif

        do j =1,npt
           syn_t(j)=b+(j-1)*dt
        enddo
!        call newhdr()
        call setfhv('b',b,nerr1)
        call setfhv('delta',dt,nerr1)
        call setnhv('npts',npt,nerr1)
        call setfhv('evla',evla,nerr1)
        call setfhv('evlo',evlo,nerr1)
        call setfhv('stla',stla,nerr1)
        call setfhv('stlo',stlo,nerr1)
        call setfhv('evdp',evdp,nerr1)
        call setkhv('kstnm',trim(kstnm),nerr1)
        call setkhv('kcmpnm',trim(kcmpnm),nerr1)
        call setkhv('knetwk',trim(knetwk),nerr1)

        ! write new seismograms
        call wsac0(newsyn_fnm,syn_t(1:npt),newsyn(1:npt),nerr1)
        if (nerr1/=0) stop 'Error reading sac file'

enddo
write(*,*) "SUCESSIVEFULLY"

end program xcorrect_syn
