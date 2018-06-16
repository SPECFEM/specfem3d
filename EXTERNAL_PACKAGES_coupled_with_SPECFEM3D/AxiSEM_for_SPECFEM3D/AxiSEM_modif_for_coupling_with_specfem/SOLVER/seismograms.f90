!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
!> Various subroutines for seismogram preparation and dumping
module seismograms

  use global_parameters
  use data_io
  use data_proc
  use data_time

  implicit none

  public

contains

!-----------------------------------------------------------------------------------------
subroutine prepare_seismograms

  use utlity
  use data_mesh, only: maxind, maxind_glob, ind_first, ind_last, &
                          surfelem, jsurfel, surfcoord, surfelem, &
                          have_epi, have_equ, have_antipode, &
                          ielantipode, ielequ, ielepi, &
                          min_distance_dim, min_distance_nondim, &
                          north, ielsolid, nel_solid, npol, router

  use data_source, only: have_src
  use commun, only: psum_int, barrier, comm_elem_number
  character(len=4)     :: appielem
  integer              :: j, iel, ielem, ind, epicount, equcount, anticount
  integer              :: iproc
  real(kind=dp)        :: s, z, r, theta

  if (lpr) &
  write(*,*)'  locating surface elements and generic receivers...'

  call barrier()

  ielantipode=1; ielequ=1
  have_epi=.false.
  have_equ=.false.
  have_antipode=.false.
  epicount=0; equcount=0;anticount=0
  ind=0

  ! test round to determine the amount of surface elements
  do iel=1,nel_solid
     ielem=ielsolid(iel)
     if (north(ielem)) then
        call compute_coordinates(s,z,r,theta,ielem,0,npol)
     else
        call compute_coordinates(s,z,r,theta,ielem,npol,0)
     endif
     if ( dabs(r-router) < min_distance_dim ) then
        ind=ind+1
     endif
  enddo

  maxind=ind
  allocate(surfelem(maxind))

  if (diagfiles) open(1000+mynum,file=datapath(1:lfdata)//'/surfelem_'//appmynum//'.dat')
  ind=0
  do iel=1,nel_solid
     ielem=ielsolid(iel)
     if (north(ielem)) then
        call compute_coordinates(s,z,r,theta,ielem,0,npol)
     else
        call compute_coordinates(s,z,r,theta,ielem,npol,0)
     endif

     if ( dabs(r-router) < min_distance_dim ) then
        ind=ind+1
        surfelem(ind)=iel
        if (diagfiles) write(1000+mynum,*) ind, iel, r, theta*180./pi

        ! find epicenter
        if (north(ielem)) then
           call compute_coordinates(s,z,r,theta,ielem,0,npol)
           if ( dabs(s) < smallval*router .and. &
                dabs(z-router) < smallval*router) then
              ielepi = iel
              if (verbose > 1) then
                 write(69,*)'Proc ', mynum, ' found: '
                 write(69,*)'Epicenter element:',ielepi
                 write(69,*)'Epicenter radius [km], colat [deg]:', &
                            r/1000.,theta/pi*180.
                 write(69,*)''
              endif
              have_epi=.true.
              epicount=epicount+1
           endif
        endif

        ! find antipode at 180deg
        if (.not. north(ielem)) then
          call compute_coordinates(s,z,r,theta,ielem,0,0)
          if ( dabs(theta-pi) < min_distance_nondim*pi) then
             ielantipode=iel
             if (verbose > 1) then
                write(69,*)'Proc ', mynum, ' found: '
                write(69,*)'Antipodal element:',ielantipode
                write(69,*)'Antipode radius [km], colat [deg]:', &
                            r/1000.,theta/pi*180.
                write(69,*)''
             endif
             have_antipode=.true.
             anticount=anticount+1
          endif
        endif

        ! find equator (take northern element)
        if (north(ielem)) then
           call compute_coordinates(s,z,r,theta,ielem,npol,npol)
           if ( dabs(z) < smallval_sngl) then
              ielequ=iel
              if (verbose > 1) then
                 write(69,*)'Proc ', mynum, ' found: '
                 write(69,*)'Equatorial element:',ielequ
                 write(69,*)'Equatorial radius [km], colat [deg]:', &
                             r/1000.,theta/pi*180.
                 write(69,*)''
              endif
              have_equ=.true.
              equcount=equcount+1
           endif
        endif

     endif
  enddo

  if (diagfiles) close(1000+mynum)

  ! make sure only one processor has each location
  if (lpr) &
     write(*,*)'  ensuring uniqueness in generic receiver locations...'

  iel=0
  if ( psum_int(epicount) > 1 ) then
     do while (iel <= nproc-1)
        call barrier
        if (mynum == iel) then
           if (have_epi .and. mynum < nproc-1) then
              do iproc=iel+1,nproc-1
                 call barrier
                 if (mynum == iproc) have_epi=.false.
                 call barrier
              enddo
              iel=nproc-1
           endif
        endif
        call barrier
        iel=iel+1
     enddo
  endif

  if ( psum_int(equcount) > 1 ) then
     do while (iel <= nproc-1)
        call barrier
        if (mynum == iel) then
           if (have_equ .and. mynum < nproc-1) then
              do iproc=iel+1,nproc-1
                 call barrier
                 if (mynum == iproc) have_equ=.false.
                 call barrier
              enddo
              iel=nproc-1
           endif
        endif
        call barrier
        iel=iel+1
     enddo
  endif

  if ( psum_int(anticount) > 1 ) then
     do while (iel <= nproc-1)
        call barrier
        if (mynum == iel) then
           if (have_antipode .and. mynum < nproc-1) then
              do iproc=iel+1,nproc-1
                 call barrier
                 if (mynum == iproc) have_antipode=.false.
                 call barrier
              enddo
              iel=nproc-1
           endif
        endif
        call barrier
     enddo
  endif

  do iel=0,nproc-1
     call barrier
     if (mynum == iel) then
        if (verbose > 1) write(69,*)'  number of surface elements:', maxind
        if (have_epi) write(*,12)procstrg,'epicenter at', &
                                thetacoord(0,npol,ielsolid(ielepi))/pi*180.
        if (have_equ) write(*,12)procstrg,'equator at', &
                                thetacoord(npol,npol,ielsolid(ielequ))/pi*180.
        if (have_antipode) write(*,12)procstrg,'antipode at', &
                                  thetacoord(0,0,ielsolid(ielantipode))/pi*180.
     endif
     call barrier
  enddo
12 format('   ',a8,'has the ',a13,f9.3,' degrees.')


  if (lpr .and. verbose > 1) write(*,*)'  communicating local numbers of surface elements...'
  call comm_elem_number(maxind, maxind_glob, ind_first, ind_last)
  if (lpr .and. verbose > 0) write(*,*)'  global number of surface elements:',maxind_glob

  ! open files for displacement and velocity traces in each surface element

  allocate(jsurfel(maxind)) ! for surface strain
  allocate(surfcoord(maxind)) ! theta of surface elements dumped for surface strain

  do iproc=0,nproc-1
     call barrier
     if (mynum == iproc) then
        call barrier
        if (diagfiles) then
            open(33333,file=datapath(1:lfdata)// &
                             '/surfelem_coords.dat',position='append')
            open(33334,file=datapath(1:lfdata)// &
                             '/surfelem_coords_jpol.dat',position='append')
            if (mynum == 0) write(33334,*)maxind_glob
            if (mynum == 0) write(33333,*)maxind_glob
        endif

        do iel=1,maxind
           if (thetacoord(npol/2,npol/2,ielsolid(surfelem(iel))) <= pi/2) then
              jsurfel(iel)=npol
           else
              jsurfel(iel)=0
           endif

           surfcoord(iel) = 180. / pi * &
                            thetacoord(npol/2,jsurfel(iel),ielsolid(surfelem(iel)))
           if (diagfiles) then
               write(33333,*) surfcoord(iel)
               write(33334,11) 180./pi* &
                    thetacoord(npol/2,jsurfel(iel),ielsolid(surfelem(iel))), &
                    (rcoord(npol/2,j,ielsolid(surfelem(iel))),j=0,npol)
           endif
        enddo
        if (diagfiles) then
            close(33333)
            close(33334)
        endif
        call barrier
     endif
     call barrier
  enddo

11 format(6(1pe11.4))


  if (dump_wavefields .and. (.not. use_netcdf) .and. (.not. coupling)) then !! VM VM add coupling
    do iel=1,maxind
       call define_io_appendix(appielem,iel+mynum*maxind)
       open(unit=40000000+iel,file=datapath(1:lfdata)// &
                 '/surfelem_disp.dat'//appielem)


       open(unit=50000000+iel,file=datapath(1:lfdata)// &
                                   '/surfelem_velo.dat'//appielem)
    enddo

    do iel=1,maxind
       call define_io_appendix(appielem,iel+mynum*maxind)
       open(unit=60000000+iel,file=datapath(1:lfdata)// &
            '/surfelem_strain.dat'//appielem)
       open(unit=70000000+iel,file=datapath(1:lfdata)// &
            '/surfelem_disp_src.dat'//appielem)
    enddo
  endif

end subroutine prepare_seismograms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read colatitudes [deg] from a file receivers.dat and locate closest grid
!! point for seismograms, output grid point locations in
!! receiver_pts.dat < PROCID>
subroutine prepare_from_recfile_seis

  use utlity
  use data_mesh, only: loc2globrec, dtheta_rec, npol, surfelem, ielsolid, num_rec_tot, &
                          north, router, recfile_el, min_distance_dim, have_epi, num_rec, maxind
  use data_source, only: src_type, rot_src, srccolat, srclon
  use commun
  use rotations, only: rotate_receivers_recfile, save_google_earth_kml
  use nc_routines, only: nc_define_outputfile

  integer                        :: i, iel, ipol, irec, num_rec_glob
  integer                        :: count_diff_loc, count_procs
  real(kind=dp)                  :: s, z, r, theta, recdist, myrecdist
  real(kind=dp)                  :: tmprecfile_th, tmprecfile_el(3)
  real(kind=dp), allocatable     :: recfile_readth(:), recfile_readph(:)
  real(kind=dp), allocatable     :: recfile_th_glob(:)
  real(kind=dp), allocatable     :: recfile_th_loc(:), recfile_el_loc(:,:)
  real(kind=dp), allocatable     :: recfile_th(:), recfile_ph_loc(:), recfile_ph_loc2(:)
  integer,allocatable            :: rec2proc(:), loc2globrec_loc(:)!, loc2globrec(:)
  character(len=4)               :: appielem
  character(len=100)             :: junk
  character(len=40), allocatable :: receiver_name(:)
  real(kind=dp)                  :: maxreclocerr

  ! Additional arrays within the STATIONS file
  character(len=20), allocatable, dimension(:) :: rec_name,rec_network
  real(kind=dp)   , allocatable, dimension(:) :: reclat,reclon,recelevation,recbury
  integer ierror

  irec=0
  count_diff_loc=0
  dtheta_rec=0. ! default for all cases but database (which has regular spacing)


  select case(trim(rec_file_type))
  case('colatlon')
     ! Read receiver location file:
     ! line1: number of receivers
     ! line2 --> line(number of receivers+1): colatitudes [deg]
     if (lpr) write(*,*)'  reading receiver colatitudes and longitudes from receivers.dat...'
     open(unit=34, file='receivers.dat', iostat=ierror, status='old', position='rewind')
     read(34,*) num_rec_glob
     write(*,*)'  # total receivers:',num_rec_glob; call flush(6)
     allocate( recfile_readth (1:num_rec_glob) )
     allocate( recfile_th_glob(1:num_rec_glob) )
     allocate( recfile_th_loc (1:num_rec_glob) )
     allocate( recfile_el_loc (1:num_rec_glob,3))
     allocate( loc2globrec_loc(1:num_rec_glob) )
     allocate( rec2proc       (1:num_rec_glob) )
     allocate( recfile_readph (1:num_rec_glob) )
     allocate( recfile_ph_loc2(1:num_rec_glob) )
     allocate( receiver_name  (1:num_rec_glob) )

     if (mynum == 0) open(unit=30,file=datapath(1:lfdata)//'/receiver_names.dat')

     do i=1,num_rec_glob
        read(34,*)recfile_readth(i),recfile_readph(i)
        call define_io_appendix(appielem,i)
        receiver_name(i) = 'recfile_'//appielem
        if (mynum == 0) write(30,*) trim(receiver_name(i)), recfile_readth(i), recfile_readph(i)
     enddo
     close(34)
     if (mynum == 0) close(30)


  case('database')
     if (lpr) write(*,*)'  generating receiver colatitudes at every element edge and midpoint...'
     if (lpr) write(*,*)'  ... which is useful for databases and sinc interpolation'
     dtheta_rec = abs( (thetacoord(npol,npol,ielsolid(surfelem(1))) &
                      - thetacoord(0,   npol,ielsolid(surfelem(1))) ) * 180. / pi ) / 2
     num_rec_glob = ceiling(180./dtheta_rec)+1
     if (lpr) then
       write(*,*) 'delta theta (mid-to-edge), (edge-to-edge) [deg]:', dtheta_rec, &
                  abs( (thetacoord(npol,npol,ielsolid(surfelem(1))) &
                      - thetacoord(0,   npol,ielsolid(surfelem(1))) ) * 180. / pi )
       write(*,*)mynum,'number of surface elements:',maxind
       write(*,*)mynum,'number of global recs (ideal,real):',(180./dtheta_rec)+1,num_rec_glob
     endif
     allocate(recfile_readth(num_rec_glob),recfile_readph(num_rec_glob))
     do i=1,num_rec_glob
        recfile_readth(i) = dtheta_rec*real(i-1)
     enddo
     if (lpr) write(*,*)mynum,'min,max receiver theta [deg]:',minval(recfile_readth),maxval(recfile_readth)
     call flush(6)
     recfile_readph = 0.
     allocate(recfile_th_loc(1:num_rec_glob),recfile_el_loc(1:num_rec_glob,3))
     allocate(loc2globrec_loc(1:num_rec_glob),rec2proc(1:num_rec_glob))
     allocate(recfile_ph_loc2(1:num_rec_glob),recfile_th_glob(1:num_rec_glob))
     allocate(receiver_name(1:num_rec_glob))
     if (mynum == 0) open(unit=30,file=datapath(1:lfdata)//'/receiver_names.dat')
     do i=1,num_rec_glob
        !call define_io_appendix(appielem,i) Does only work for nrec < 9999
        !receiver_name(i) = 'recfile_'//appielem
        write(receiver_name(i),112) i
112     format('recfile_',I6.6)
        if (mynum == 0) write(30,*)trim(receiver_name(i)),recfile_readth(i),recfile_readph(i)
     enddo
     if (mynum == 0) close(30)
     if (lpr) write(*,*)mynum,'done with database receiver writing.';call flush(6)


  case('stations')
     if (lpr) write(*,*)'  reading receiver colatitudes and longitudes from STATIONS...'
     ! count number of receivers first
     open(unit=34,file='STATIONS',iostat=ierror,status='old',action='read',position='rewind')
     num_rec_glob = 0
     do while(ierror == 0)
        read(34,*,iostat=ierror) junk
        if (ierror == 0) num_rec_glob = num_rec_glob + 1
     enddo
     close(34)
     if (lpr) write(*,*)'  ...counted number of stations:', num_rec_glob

     allocate(recfile_readth(1:num_rec_glob),recfile_th_glob(1:num_rec_glob))
     allocate(recfile_th_loc(1:num_rec_glob),recfile_el_loc(1:num_rec_glob,3))
     allocate(loc2globrec_loc(1:num_rec_glob),rec2proc(1:num_rec_glob))
     allocate(recfile_readph(1:num_rec_glob))
     allocate(recfile_ph_loc2(1:num_rec_glob))
     allocate(rec_name(num_rec_glob),rec_network(num_rec_glob))
     allocate(reclat(num_rec_glob),reclon(num_rec_glob),recelevation(num_rec_glob),recbury(num_rec_glob))
     allocate(receiver_name(num_rec_glob))
     open(unit=34,file='STATIONS',iostat=ierror,status='old',action='read',position='rewind')
     if (mynum == 0) open(unit=30,file=datapath(1:lfdata)//'/receiver_names.dat')

     do i=1,num_rec_glob
        read(34,*)rec_name(i),rec_network(i),reclat(i),reclon(i),recelevation(i),recbury(i)
        if (reclon(i) <= zero) then
           recfile_readph(i)=reclon(i)+360.d0
        else
           recfile_readph(i)= reclon(i)
        endif
        recfile_readth(i) = 90.d0 - reclat(i)
        receiver_name(i) = trim(rec_name(i))//'_'//trim(rec_network(i))
        if (mynum == 0) write(30,*)trim(receiver_name(i)),recfile_readth(i),recfile_readph(i)
     enddo
     close(34)
     if (mynum == 0) close(30)


  case default
     write(*,*)procstrg, 'Undefined receiver file format!!'
     stop
  end select !Receiver type rec_file_type

  num_rec_tot = num_rec_glob ! to be known later/globally

  ! check on consistency of receiver coordinates

  if (minval(recfile_readph) < 0.d0) then
     if (lpr) write(*,*)' ERROR: We do not allow negative receiver longitudes....'
     stop
  endif

  if (maxval(recfile_readph) > 360.001) then
     if (lpr) write(*,*)' ERROR: We do not allow receiver longitudes larger than 360 degrees....'
     stop
  endif

  if (maxval(recfile_readth) < 0.d0) then
     if (lpr) write(*,*)' ERROR: We do not allow negative receiver colatitudes....'
     stop
  endif

  if (maxval(recfile_readth) > 180.001) then
     if (lpr) write(*,*)' ERROR: We do not allow receiver colatitudes larger than 180 degrees....'
     stop
  endif

  ! rotate receiver locations if source is not located at north pole
  if (rot_src ) then
     call rotate_receivers_recfile(num_rec_glob, recfile_readth, recfile_readph, receiver_name)
  else
    ! Why only save the kml if the source is at the northpole?
    ! The kml file is saved in rotate_receivers_recfile as well.
    if ((lpr) .and. (.not. (rec_file_type == 'database'))) then
      call save_google_earth_kml( real(srccolat * 180.0 / pi), real(srclon * 180.d0 / pi), &
                                  real(recfile_readth), real(recfile_readph), &
                                  num_rec_glob, 'original', receiver_name  )
    endif
  endif

  recfile_th_glob(:) = zero
  rec2proc(:) = 0

  ! find closest grid points

  do i=1, num_rec_glob
     if (verbose > 1) write(69,*)'  working on receiver #',i,recfile_readth(i)*180./pi
     recdist=10.d0*router
     do iel=1,maxind
        do ipol=0,npol
           if (north(ielsolid(surfelem(iel)))) then ! NORTH
              call compute_coordinates(s,z,r,theta,ielsolid(surfelem(iel)), &
                                       ipol,npol)
              if (z < zero ) then
               write(*,*)'PROBLEM! north but z < 0: ', &
                           north(ielsolid(surfelem(iel))),z
                 write(*,*)'r,theta:',r/1000.,theta*180./pi
                 write(*,*)iel,surfelem(iel),ielsolid(surfelem(iel))
                 stop
              endif

           else ! SOUTH
              call compute_coordinates(s,z,r,theta,ielsolid(surfelem(iel)), &
                                       ipol,0)
              if (z > zero ) then
                 write(*,*)'PROBLEM! south but z > 0: ', &
                           north(ielsolid(surfelem(iel))),z
                 write(*,*)'r,theta:',r/1000.,theta*180./pi
                 write(*,*)iel,surfelem(iel),ielsolid(surfelem(iel))
                 stop
              endif

           endif

           if (dabs(theta/pi*180.d0-recfile_readth(i)) < recdist) then
              recdist=dabs(theta/pi*180.d0-recfile_readth(i))
              tmprecfile_th=theta/pi*180.d0
              tmprecfile_el(1)=surfelem(iel) ! only in the solid domain
              tmprecfile_el(2)=ipol
              if (north(ielsolid(surfelem(iel)))) tmprecfile_el(3)=npol
              if (.not. north(ielsolid(surfelem(iel)))) tmprecfile_el(3)=0
           endif

        enddo
     enddo

     ! Make sure only one processor takes on each location
     myrecdist=recdist
     recdist=pmin(recdist)
     count_procs=0
     if (dblreldiff_small(myrecdist,recdist)) count_procs=mynum
     !take as default the larger processor ID to take on the receiver
     count_procs=pmax_int(count_procs)
     if (mynum == count_procs) then
        irec = irec+1
        if (verbose > 1) write(69,*)'found local grid point and processor...',irec,i
        recfile_th_loc(irec)     = tmprecfile_th
        recfile_el_loc(irec,1:3) = tmprecfile_el(1:3)
        loc2globrec_loc(irec)    = i
        rec2proc(i)              = mynum
        recfile_th_glob(i)       = tmprecfile_th
     endif

     ! longitude
     if (irec > 0)  recfile_ph_loc2(irec) = recfile_readph(i)*pi/180.

     ! Can do that since NOW only one proc has non-zero values
     recfile_th_glob(i) = psum_dble(recfile_th_glob(i))
     rec2proc(i) = psum_int(rec2proc(i))

  enddo ! num_rec_glob

  ! Form local arrays depending on how many receivers each processor has
  num_rec = irec
  allocate(recfile_el(1:num_rec,1:3),loc2globrec(1:num_rec))
  allocate(recfile_th(1:num_rec))
  allocate(recfile_ph_loc(1:num_rec))
  allocate(fname_rec_seis(1:num_rec))
  allocate(fname_rec_velo(1:num_rec))

  recfile_el(1:num_rec,1:3)=recfile_el_loc(1:num_rec,1:3)
  loc2globrec(1:num_rec)=loc2globrec_loc(1:num_rec)
  recfile_th(1:num_rec)=recfile_th_loc(1:num_rec)
  recfile_ph_loc(1:num_rec)=recfile_ph_loc2(1:num_rec)
  deallocate(recfile_ph_loc2)

  ! How many receivers does each processor have, do they sum to global number?
  if ( psum_int(num_rec) /= num_rec_glob ) then
     write(*,*)'PROBLEM: sum of local receivers is different than global!'
     if (lpr) write(*,*)'Global number of receivers:',num_rec_glob
     write(*,*)procstrg,'Number of receivers:',num_rec
     stop
  endif

  if (verbose > 1) then
     if (lpr) write(*,*)
     do irec=0,nproc-1
        call barrier
        if (mynum == irec) write(*,14)procstrg,num_rec,num_rec_glob
        call barrier
     enddo
     write(69,14)procstrg,num_rec,num_rec_glob
  endif
14 format(/,'   ',a8,'has',i4,' out of',i6,' receivers')

  ! Output colatitudes globally (this is the file needed to plot seismograms)
  ! NOTE: recfile_th is in degrees!!!
  if (lpr) then
     open(99997,file=datapath(1:lfdata)//'/receiver_pts.dat')
     do i=1,num_rec_glob
        write(99997,*)recfile_th_glob(i),recfile_readph(i),rec2proc(i)
     enddo
     close(99997)
  endif

  ! Output colatitudes locally to infopath (this file is for info purposes!)

  maxreclocerr=zero
  if (diagfiles) then
      open(9998+mynum,file=infopath(1:lfinfo)//'/receiver_pts_'//appmynum//'.dat')
      write(9998+mynum,*)num_rec
  endif
  do i=1,num_rec ! Only over newly found local receiver locations

     if (diagfiles) write(9998+mynum,13) i, recfile_readth(loc2globrec(i)), &
                                         recfile_th(i), recfile_ph_loc(i), &
                                         recfile_el(i,1), recfile_el(i,2), &
                                         recfile_el(i,3)

     call define_io_appendix(appielem,loc2globrec(i))

     if ( pi/180*router*abs(recfile_readth(loc2globrec(i))-recfile_th(i) ) >&
          min_distance_dim) then
        count_diff_loc=count_diff_loc+1
        if (verbose > 1) write(*,22)procstrg, &
                  recfile_readth(loc2globrec(i)), recfile_th(i)
        if (dabs(recfile_readth(loc2globrec(i))-recfile_th(i)) > maxreclocerr) &
             maxreclocerr=dabs(recfile_readth(loc2globrec(i))-recfile_th(i))/ &
                          180.*pi*router

     endif

22 format('   WARNING:',a8,' rec. location file/mesh:',2(f9.3))

     ! Test: receivers on the surface?
     call compute_coordinates(s,z,r,theta,ielsolid(recfile_el(i,1)), &
          recfile_el(i,2),recfile_el(i,3))
     if (.not. dblreldiff_small(r,router)) then
        print *
        write(*,*)'PROBLEM: receiver is not at the surface!'
        write(*,*)'r [km], colat [deg]:',r/1000.,theta/pi*180.
        stop
     endif


     if (.not. (use_netcdf)) then
         if (verbose > 1) write(*,*)'  ',procstrg,'opening receiver file:',i,appielem
         open(100000+i,file=datapath(1:lfdata)//'/'//trim(receiver_name(loc2globrec(i)))//'_disp.dat')
     endif

  enddo
  if (diagfiles) close(9998+mynum)
  if (use_netcdf) then
    call nc_define_outputfile(num_rec_glob, receiver_name, recfile_th_glob, &
                              recfile_readth, recfile_readph, rec2proc)
  endif

  if (verbose > 1) then
     write(69,15)count_diff_loc,num_rec
     write(69,*)'  Maximal receiver location error [m]:',maxreclocerr
     write(69,*)
  endif
15 format(i4,' out of',i4,' receivers are located at wrong points.')

  maxreclocerr = pmax(maxreclocerr)
  if (lpr) then
     write(*,*)
     write(*,*)'  maximal receiver location error [m]:',maxreclocerr
     write(*,*)
  endif

  ! define general prefactor for all cases
  allocate(recfac(num_rec,5))

  if (lpr) write(*,*) &
           '  Calculating prefactors for cylindrical components...'
  recfac(:,:) = 1.d0
  recfac(:,2) = 0.d0
  recfac(:,4) = 0.d0

  ! no phi component for monopole
  if (src_type(1) == 'monopole')  recfac(:,3) = 0.d0

  if (diagfiles) then
      open(300+mynum,file=infopath(1:lfinfo)//'/receiver_recfac_'//appmynum//'.dat')
      do i=1,num_rec
         write(300+mynum,*)recfile_th(i),recfile_ph_loc(i)*180./pi,recfac(i,3)
         write(300+mynum,*)recfac(i,1),recfac(i,2),recfac(i,4),recfac(i,5)
         write(300+mynum,*)
      enddo
      close(300+mynum)
  endif

13 format(i3,3(1pe12.4),i8,2(i2))

  deallocate(recfile_ph_loc,recfile_readph)
  deallocate(recfile_readth)
  deallocate(recfile_th_glob,recfile_th,recfile_th_loc)
  deallocate(recfile_el_loc,loc2globrec_loc,rec2proc)

  if (rec_file_type == 'stations') then
    deallocate(rec_name,rec_network)
    deallocate(reclat,reclon,recelevation,recbury)
  endif

end subroutine prepare_from_recfile_seis
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Open files for generic checks: hypocenter,epicenter,equator,antipode.
!! File output names: seis < LOCATION>{1,2,3}.dat
!!                    where 1=s-component, 2=phi-component, 3=z-component.
!! At hypocenter,epicenter,antipode, s-comp = longitudinal comp.
!!                                   z-comp = radial/vertical comp.
!! At equator                        s-comp = radial/vertical comp.
!!                                   z-comp = longitudinal comp.
!! Not writing the phi/transverse component for monopole sources.
subroutine open_hyp_epi_equ_anti

  use data_source, only: have_src, src_type
  use data_mesh, only: have_epi, have_equ, have_antipode, maxind

  if (maxind > 0) then

    if (have_epi) then
       open(900,file=datapath(1:lfdata)//'/seisepicenter1.dat')
       if (src_type(1) /= 'monopole') &
            open(903,file=datapath(1:lfdata)//'/seisepicenter2.dat')
       open(906,file=datapath(1:lfdata)//'/seisepicenter3.dat')
    endif

    if (have_equ) then
       open(902,file=datapath(1:lfdata)//'/seisequator1.dat')
       if (src_type(1) /= 'monopole') &
            open(904,file=datapath(1:lfdata)//'/seisequator2.dat')
       open(907,file=datapath(1:lfdata)//'/seisequator3.dat')
    endif

    if (have_antipode) then
       open(901,file=datapath(1:lfdata)//'/seisantipode1.dat')
       if (src_type(1) /= 'monopole') &
            open(905,file=datapath(1:lfdata)//'/seisantipode2.dat')
       open(908,file=datapath(1:lfdata)//'/seisantipode3.dat')
    endif

 endif ! if maxind > 0

end subroutine open_hyp_epi_equ_anti
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Generic seismograms for quick checks: hypocenter,epicenter,equator,antipode.
!! Not writing the transverse component for monopole sources.
!! See open_hyp_epi_equ_anti for component explanation.
subroutine compute_hyp_epi_equ_anti(t,disp)

  use data_mesh, only: have_epi, have_equ, have_antipode, npol, &
                         ielepi, ielantipode, ielequ, recfile_el, maxind
  use data_source, only: iel_src,ipol_src,jpol_src,have_src,src_type
  real(kind=dp)                   :: t
  real(kind=realkind), intent(in) :: disp(0:,0:,:,:)

  if (maxind > 0) then
     ! epicenter
     if (have_epi) then
        if (src_type(1) == 'dipole') then
           write(900,*)t,disp(0,npol,ielepi,1)+disp(0,npol,ielepi,2) ! s
           write(903,*)t,disp(0,npol,ielepi,1)-disp(0,npol,ielepi,2) ! phi
        else
           write(900,*)t,disp(0,npol,ielepi,1) ! s
           if (src_type(1) == 'quadpole') &
                write(903,*)t,disp(0,npol,ielepi,2) ! phi
        endif
        write(906,*)t,disp(0,npol,ielepi,3)  ! z
     endif

     ! antipode
     if (have_antipode) then
        if (src_type(1) == 'dipole') then
           write(901,*)t,disp(0,0,ielantipode,1)+disp(0,0,ielantipode,2) ! s
           write(905,*)t,disp(0,0,ielantipode,1)-disp(0,0,ielantipode,2) ! phi
        else
           write(901,*)t,disp(0,0,ielantipode,1) ! s
           if (src_type(1) == 'quadpole') &
                write(905,*)t,disp(0,0,ielantipode,2) ! phi
        endif
        write(908,*)t,disp(0,0,ielantipode,3)  ! z
     endif

     ! equator
     if (have_equ) then
        if (src_type(1) == 'dipole') then
           write(902,*)t,disp(npol,npol,ielequ,1)+disp(npol,npol,ielequ,2) ! s
           write(904,*)t,disp(npol,npol,ielequ,1)-disp(npol,npol,ielequ,2) ! phi
        else
           write(902,*)t,disp(npol,npol,ielequ,1) ! s
           if (src_type(1) == 'quadpole') &
               write(904,*)t,disp(npol,npol,ielequ,2) ! phi
        endif
        write(907,*)t,disp(npol,npol,ielequ,3)  ! z
     endif

  endif ! if maxind > 0

end subroutine compute_hyp_epi_equ_anti
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_recfile_seis_bare(disp)

  use data_source, only: src_type
  use data_mesh, only: recfile_el, num_rec

  real(kind=realkind), intent(in) :: disp(0:,0:,:,:)

  integer :: i

   if (src_type(1) == 'monopole') then
      do i=1,num_rec
         ! order: u_s, u_z
         write(100000+i,*) &
             disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1), &
             disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3)
      enddo
   else if (src_type(1) == 'dipole') then
      do i=1,num_rec
         ! order: u_s, u_phi,u_z
         write(100000+i,*) &
              disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1) &
              + disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2), &
              disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1) &
              - disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2), &
              disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3)
      enddo
   else if (src_type(1) == 'quadpole') then
      do i=1,num_rec
         ! order: u_s, u_phi,u_z
         write(100000+i,*) &
             disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),1), &
             disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),2), &
             disp(recfile_el(i,2),recfile_el(i,3),recfile_el(i,1),3)
      enddo
  endif !src_type(1)

end subroutine compute_recfile_seis_bare
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!! Calculate displacement at receiver locations and pass to nc_dump_rec
subroutine nc_compute_recfile_seis_bare(disp, iseismo)

  use data_source, only: src_type
  use nc_routines, only: nc_dump_rec
  use data_mesh, only: recfile_el, num_rec

  real(kind=realkind), intent(in)  :: disp(0:,0:,:,:)
  integer,             intent(in)  :: iseismo

  real(kind=realkind)              :: disp_rec(3,num_rec)
  integer                          :: i


  if (src_type(1) == 'monopole') then
     do i=1, num_rec
          disp_rec(1,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),1)
          disp_rec(2,i) = 0
          disp_rec(3,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),3)
     enddo
  else if (src_type(1) == 'dipole') then
     do i=1, num_rec
          disp_rec(1,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),1) &
                        + disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),2)
          disp_rec(2,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),1) &
                        - disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),2)
          disp_rec(3,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),3)
     enddo
  else if (src_type(1) == 'quadpole') then
     do i=1, num_rec
          disp_rec(1,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),1)
          disp_rec(2,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),2)
          disp_rec(3,i) = disp(recfile_el(i,2), recfile_el(i,3), recfile_el(i,1),3)
     enddo
  endif !src_type(1)

  call nc_dump_rec(disp_rec, iseismo)

end subroutine nc_compute_recfile_seis_bare
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Save one displacement and velocity trace for each element on the surface
!! which are both needed for kernels (du and v0 inside the cross-correlation)
subroutine compute_surfelem(disp, velo)

  use data_source, only: src_type
  use nc_routines, only: nc_dump_surface
  use data_mesh, only: npol, jsurfel, surfelem, maxind

  real(kind=realkind), intent(in) :: disp(0:,0:,:,:)
  real(kind=realkind), intent(in) :: velo(0:,0:,:,:)
  real                            :: dumpvar(maxind, 3)
  integer                         :: i

  dumpvar = 0.0

  if (use_netcdf) then
      if (src_type(1) == 'monopole') then
          do i=1,maxind
              dumpvar(i,1) = real(disp(npol/2,jsurfel(i),surfelem(i),1))
              dumpvar(i,2) = real(disp(npol/2,jsurfel(i),surfelem(i),3))
          enddo
      else
          do i=1,maxind
              dumpvar(i,:) = real(disp(npol/2,jsurfel(i),surfelem(i),:))
          enddo
      endif !monopole
      call nc_dump_surface(dumpvar(:,:), 'disp')

      if (src_type(1) == 'monopole') then
          do i=1,maxind
              dumpvar(i,1) = real(velo(npol/2,jsurfel(i),surfelem(i),1))
              dumpvar(i,2) = real(velo(npol/2,jsurfel(i),surfelem(i),3))
          enddo
      else
          do i=1,maxind
              dumpvar(i,:) = real(velo(npol/2,jsurfel(i),surfelem(i),:))
          enddo
      endif !monopole
      call nc_dump_surface(dumpvar(:,:), 'velo')

  else if (coupling) then !!!! SB coupling

  else !use_netcdf !! VM VM add coupling
      if (src_type(1) == 'monopole') then
      do i=1,maxind
         write(40000000+i,*)real(disp(npol/2,jsurfel(i),surfelem(i),1)), &
                            real(disp(npol/2,jsurfel(i),surfelem(i),3))
         write(50000000+i,*)real(velo(npol/2,jsurfel(i),surfelem(i),1)), &
                            real(velo(npol/2,jsurfel(i),surfelem(i),3))
      enddo

      else
         do i=1,maxind
         write(40000000+i,*)real(disp(npol/2,jsurfel(i),surfelem(i),1)), &
                            real(disp(npol/2,jsurfel(i),surfelem(i),2)), &
                            real(disp(npol/2,jsurfel(i),surfelem(i),3))
         write(50000000+i,*)real(velo(npol/2,jsurfel(i),surfelem(i),1)), &
                            real(velo(npol/2,jsurfel(i),surfelem(i),2)), &
                            real(velo(npol/2,jsurfel(i),surfelem(i),3))
         enddo
      endif !monopole
  endif !netcdf

end subroutine compute_surfelem
!-----------------------------------------------------------------------------------------

end module seismograms
!=========================================================================================
