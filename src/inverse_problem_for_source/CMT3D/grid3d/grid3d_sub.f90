module grid3d_sub

  use grid3d_sub2
  implicit none

contains

  subroutine set_parameters(par_file)

    character(len=*),intent(in) :: par_file
    ! cmt_file
    integer :: yr,mo,jda,ho,mi
    real*8 :: sec,t_cmt,hdur,elat,elon,depth
    real*8 :: moment_tensor(NM)
    real :: mw

    integer :: ios

    open(IOPAR,file=par_file,status='old',iostat=ios)
    if (ios /= 0) stop 'Error opening inversion parameter file'

    read(IOPAR,'(a)') cmt_file
    read(IOPAR,*) dmoment
    read(IOPAR,'(a)') flexwin_out_file
    read(IOPAR,'(l7)',iostat=ios, advance='no') weigh_data_files
    if (ios /= 0) stop 'Error reading weigh_data_files'

    read(IOPAR,'(l7)',iostat=ios) read_weight
    if (ios /= 0) read_weight=.false.

    read(IOPAR,*) comp_z_weight, comp_t_weight, comp_r_weight, &
         az_exp_weight, &
         pnl_dist_weight, rayleigh_dist_weight, love_dist_weight

    read(IOPAR,*) station_correction, tshift_max
    read(IOPAR,*) global_search, ncalc
    if (.not. global_search) then
       ncalc = 1
       read(IOPAR,*,iostat=ios) s_strike, e_strike, d_strike
       if (ios /= 0) stop 'Error reading strike parameters'
       read(IOPAR,*,iostat=ios) s_dip, e_dip, d_dip
       if (ios /= 0) stop 'Error reading dip parameters'
       read(IOPAR,*,iostat=ios) s_rake, e_rake, d_rake
       if (ios /= 0) stop 'Error reading rake parameters'
       read(IOPAR,*,iostat=ios) s_mw, e_mw, d_mw
       if (ios /= 0) stop 'Error reading Mw parameters'
    else
       if (ncalc < 0) stop 'Error: ncalc should > 0'
       if (ncalc > 10) stop 'Too many ncalc, should be < 10'
       ! read junk values
       read(IOPAR,*,iostat=ios) s_strike, e_strike, d_strike
       read(IOPAR,*,iostat=ios) s_dip, e_dip, d_dip
       read(IOPAR,*,iostat=ios) s_rake, e_rake, d_rake
       read(IOPAR,*,iostat=ios) s_mw, e_mw, d_mw

       call get_cmt(cmt_file,yr,mo,jda,ho,mi,sec, &
            t_cmt,hdur,elat,elon,depth,moment_tensor)
       mw=log10(dsqrt((sum(moment_tensor(1:3)**2)+2*sum(moment_tensor(4:6)**2))/2.))/1.5-10.73

       ! these are global search control parameters -- can be adjusted
       ! 19x10x13x5
       s_strike = 0; e_strike = 180; d_strike = 10
       s_dip = 0; e_dip = 90;  d_dip = 10
       s_rake = -180; e_rake = 180; d_rake = 30
       s_mw = mw * 0.9; e_mw = mw * 1.1; d_mw = mw * 0.05
       t_strike=1.5; t_dip=1.5; t_rake=1.5; t_mw=1
    endif

    read(IOPAR,*) write_new_cmt
    if (write_new_cmt)  new_cmt_file=trim(cmt_file)//'_GRD'

    ! convert to radians
    s_strike=s_strike*d2r; e_strike=e_strike*d2r; d_strike=d_strike*d2r
    s_dip=s_dip*d2r; e_dip=e_dip*d2r; d_dip=d_dip*d2r
    s_rake=s_rake*d2r; e_rake=e_rake*d2r; d_rake=d_rake*d2r

    ! print more information here

    if (DEBUG) then
       write(*,*)
       write(*,*) ' ==  ==  ==  == GRID3D.PAR Summary ==  ==  == ='
       write(*,*) 'original cmt file: ', trim(cmt_file)
       write(*,*) 'new cmt file: ', trim(new_cmt_file)
       write(*,*)
       if (weigh_data_files) then
          if (read_weight) then
             write(*,*) 'Weights of data will be read from input file'
          else
             write(*,*) 'Weighing data files according to ...'
             write(*,'(a,3g15.5)') '  Z, T, R comp (lin) = ', comp_z_weight, comp_t_weight, comp_r_weight
             write(*,'(a,g15.5)') '  azimuth bin (exp) = ', az_exp_weight
             write(*,'(a,3g15.5)') '  Pnl/R/S dist (exp) = ', pnl_dist_weight, rayleigh_dist_weight, love_dist_weight
          endif
       else
          write(*,*) 'no weighing of data'
       endif
       if (station_correction) then
          write(*,*) 'shift data to aligh with synthetics before grid search'
          if (tshift_max > 0) then
            write(*,'(a,g15.5)') '  with maximum shift allowed', tshift_max
          else
            stop 'tshift_max should be > 0'
          endif
       else
          write(*,*) 'data are NOT shifted before grid search'
       endif
       write(*,*)
       if (global_search) then
          write(*,*) 'global search over '
          write(*,*) 'strike = ', s_strike,e_strike,d_strike
          write(*,*) 'dip = ', s_dip,e_dip,d_dip
          write(*,*) 'rake = ', s_rake,e_rake,d_rake
          write(*,*) 'Mw = ', s_mw,e_mw,d_mw
          write(*,*) ' with search depth level ', ncalc
       else
          write(*,*) 'local search over'
          write(*,*) 'strike = ', s_strike,e_strike,d_strike
          write(*,*) 'dip = ', s_dip,e_dip,d_dip
          write(*,*) 'rake = ', s_rake,e_rake,d_rake
          write(*,*) 'Mw = ', s_mw,e_mw,d_mw
       endif
       write(*,*) ' ==  ==  ==  == END Summary ==  ==  ==  ==  ==  == '
       write(*,*)
    endif

  end subroutine set_parameters

!=======================================================

  subroutine setup_data_weights

    integer :: ios,i,j
    character*8, dimension(NRECMAX) :: kstnm,knetwk,kcmpnm
    real, dimension(NRECMAX) :: azimuth, dist_deg, dist_km
    real :: tstart, tend, tjunk
    logical :: lexd, lexs


    open(IOWIN,file=trim(flexwin_out_file),status='old',iostat=ios)
    if (ios /= 0) stop 'Flexwin output file can not be found'

    read(IOWIN,*) nfiles
    if (nfiles > NRECMAX) stop 'Increase NRECMAX limit'
    if (nfiles < 0) stop 'Check nfiles < 0'

    nwin_total = 0
    do i = 1, nfiles
       read(IOWIN,'(a)') data_file
       read(IOWIN,'(a)') syn_file
       inquire(file=data_file,exist=lexd)
       inquire(file=syn_file,exist=lexs)
       if (.not. (lexd .and. lexs)) then
          write(*,*) 'Check data and syn file ', trim(data_file), ' and ', trim(syn_file)
          stop
       endif
       read(IOWIN,*) nwins(i)
       if (nwins(i) < 0) stop 'Check nwins(i) '
       do j = 1, nwins(i)
          if (read_weight) then
             ! LQY: tjunk: time shift for data/syn within the window
             ! maybe used in the future to replace local corr subroutine
             read(IOWIN,*,iostat=ios) tstart, tend, tjunk, data_weights(nwin_total+j)
             if (ios /= 0) stop 'ts, te, tshift, weight are expected!'
          else
             read(IOWIN,*) tstart, tend
          endif
       enddo
       nwin_total = nwin_total + nwins(i)

       call read_sac_info(trim(syn_file)//'.'//trim(PAR_NAME(1)),data_file, &
            kstnm(i),kcmpnm(i),knetwk(i),azimuth(i),dist_deg(i),dist_km(i))
    enddo

    write(*,*) 'Total number of windows = ', nwin_total
    close(IOWIN)
    if (nwin_total > NWINMAX) stop 'Exceeding NWINMAX limit'

    if (weigh_data_files) then
       if (.not. read_weight) then
          ! THIS function NEEDS TO BE MODIFIED FOR DIFFERENT SEISMIC SCENARIOS
          call compute_data_weights(kcmpnm,azimuth,dist_km,data_weights)
       endif
    else
       data_weights(1:nwin_total) = 1.
    endif

  end subroutine setup_data_weights

!========================================================

  subroutine grid_search

    integer icalc,n_strike,n_dip,n_rake,n_mw,n_total
    real :: strike,dip,rake,mw,moment,mijn(NM)
    integer :: i, j, nf,nw, nwint
    real :: tstart(NWINMAX), tend(NWINMAX)

    do icalc = 1, ncalc

       if (global_search) print *, '---- icalc= ', icalc, '  -----'
       if (DEBUG) then
          write(*,*) 'SEARCH range -- ', icalc
          write(*,*) '  strike = ', s_strike*r2d,e_strike*r2d,d_strike*r2d
          write(*,*) '  dip = ', s_dip*r2d,e_dip*r2d,d_dip*r2d
          write(*,*) '  rake = ', s_rake*r2d,e_rake*r2d,d_rake*r2d
          write(*,*) '  Mw = ', s_mw,e_mw,d_mw
       endif

       n_strike=nint((e_strike-s_strike)/d_strike)+1
       n_dip=nint((e_dip-s_dip)/d_dip)+1
       n_rake=nint((e_rake-s_rake)/d_rake)+1
       n_mw=nint((e_mw-s_mw)/d_mw)+1

       n_total =n_strike*n_dip*n_rake* n_mw

       if (n_strike < 0 .or. n_dip < 0 .or. n_rake < 0 .or. n_mw < 0) &
            stop 'Error search dimension'
       if (n_total > NMEC_MAX) stop 'Error total num of strike*dip*rake*nw'

       if (DEBUG) write(*,'(a,i8,a,4i4)') 'search dimension:', &
            n_total, ' = ', n_strike,n_dip,n_rake,n_mw

       misfit = 0.

       ! precompute mij from strike,dip and rake
       call compute_mij_from_sdr(s_strike,d_strike,n_strike, &
            s_dip,d_dip,n_dip,s_rake,d_rake,n_rake,s_mw,d_mw,n_mw,mij)

       ! read flexwin output to compute misfit function values
       nwint = 0
       open(IOWIN,file=trim(flexwin_out_file),status='old')
       read(IOWIN,*) nf
       do i = 1, nf
          read(IOWIN,'(a)') data_file
          read(IOWIN,'(a)') syn_file
          if (DEBUG) write(*,'(a)') trim(data_file)
          read(IOWIN,*) nw
          if (nw <= 0) stop 'Error number of windows = 0'
          do j = 1, nw
             read(IOWIN,*) tstart(j), tend(j)
             nwint = nwint + 1
          enddo
          call add_misfit(data_file,syn_file, &
               data_weights(nwint-nw+1:nwint),tstart,tend,nw, &
               mij,misfit,n_total)
       enddo ! nfiles
       if (nwint /= nwin_total) stop 'Error counting total number of windows'
       close(IOWIN)

       ! choose the minimum solution
       print *, 'Output best solution -- ', icalc
       call select_best_solution(icalc,misfit,n_strike,n_dip,n_rake,n_mw, &
            s_strike,s_dip,s_rake,s_mw,strike,dip,rake,mw)

       print *, 'The best solution is '
       print *, '  strike = ', strike * r2d, ' +/- ', d_strike * r2d
       print *, '  dip = ', dip * r2d, ' +/- ', d_dip * r2d
       print *, '  rake = ', rake * r2d, ' +/- ', d_rake * r2d
       print *, '  mw = ', mw, '+/-', d_mw

       if (write_new_cmt) then
          moment=10 ** ((mw + 10.73) * 1.5)
          call sdr2moment(strike,dip,rake,moment, &
               mijn(1),mijn(2),mijn(3),mijn(4),mijn(5),mijn(6))
          call write_new_cmtsolution(cmt_file,trim(new_cmt_file),mijn)
       endif

       ! if global search, reassign mw, strike, dip and rake search ranges
       if (global_search .and. icalc < ncalc) then
          s_mw=mw-t_mw*d_mw; e_mw=mw+t_mw*d_mw
          s_strike=strike-t_strike*d_strike; e_strike=strike+t_strike*d_strike
          ! dip is not continuous
          s_dip=max(0.,dip-t_dip*d_dip); e_dip = min(pi/2,dip+t_dip*d_dip)
          s_rake=rake-t_rake*d_rake; e_rake = rake+t_rake*d_rake
          d_mw=d_mw/2; d_strike=d_strike/4; d_dip=d_dip/4; d_rake=d_rake/4
       endif

    enddo ! icalc

  end subroutine grid_search


end module grid3d_sub
