module cmt3d_sub

    use cmt3d_constants
    use cmt3d_variables
    use cmt3d_sub2

    implicit none


contains

!===========================================================

  subroutine set_parameters(par_file)

    character(len=*), intent(in) :: par_file
    integer ios
    character(len=250) :: line

    open(IOPAR,file=par_file,status='old',iostat=ios)
    if (ios /= 0) stop 'Error opening inversion parameter file'

    read(IOPAR,'(a)') cmt_file
    read(IOPAR,'(a)') new_cmt_file
    read(IOPAR,'(a)') line
    read(line,*, iostat=ios) npar
    if (ios /= 0) stop 'Error reading npar'
    read(line(3:),*,iostat=ios) global_coord
    ! if no global_coord specified
    if (ios /= 0)  global_coord=.false.

    if (.not. global_coord) then ! local coordinates (DEP, LON, LAT)
       par_name = (/'Mrr','Mtt','Mpp','Mrt','Mrp', 'Mtp','dep','lon','lat', &
            'ctm','hdr'/)
       SCALE_PAR =  &
            (/ SCALE_MOMENT, SCALE_MOMENT, SCALE_MOMENT, SCALE_MOMENT, &
            SCALE_MOMENT, SCALE_MOMENT, SCALE_DEPTH, SCALE_DELTA, SCALE_DELTA, &
            SCALE_CTIME, SCALE_HDUR /)
    else  ! global coordinates (X, Y, Z)
       par_name = (/'Mxx','Myy','Mzz','Mxy','Mxz', 'Myz','xxx','yyy','zzz', &
            'ctm','hdr'/)
       if (npar == 7) stop 'depth only inversion is not allowed in global coordinates'
       SCALE_PAR =  &
            (/ SCALE_MOMENT, SCALE_MOMENT, SCALE_MOMENT, SCALE_MOMENT, &
            SCALE_MOMENT, SCALE_MOMENT, SCALE_DEPTH, SCALE_DEPTH, SCALE_DEPTH, &
            SCALE_CTIME, SCALE_HDUR /)
    endif

    read(IOPAR,*) ddelta,ddepth,dmoment
    if (.not. global_coord) then
       dcmt_par = (/dble(dmoment),dble(dmoment),dble(dmoment), &
            dble(dmoment),dble(dmoment),dble(dmoment), &
            dble(ddepth),dble(ddelta),dble(ddelta), &
            1.0d0, 1.0d0/) / SCALE_PAR
    else ! global code have different dcmt_par
       dcmt_par = (/dble(dmoment),dble(dmoment),dble(dmoment), &
            dble(dmoment),dble(dmoment),dble(dmoment), &
            dble(ddepth),dble(ddepth),dble(ddepth), &
            1.0d0, 1.0d0/) / SCALE_PAR
    endif

    read(IOPAR,'(a)') flexwin_out_file

! now allow the possibility of reading in weights from input file
    read(IOPAR, '(a)') line
    read(line,*,iostat=ios) weigh_data_files
    if (ios /= 0) stop 'Error reading weigh_data_files'

    read(line(8:),'(l7)',iostat=ios) read_weight
    if (ios /= 0) read_weight=.false.

! here we assume that Pnl wave window is the first window out
! of possible two windows selected for Z and R components
    read(IOPAR,*) comp_z_weight, comp_t_weight, comp_r_weight, &
         az_exp_weight, &
         pnl_dist_weight, rayleigh_dist_weight, love_dist_weight
    read(IOPAR,*) station_correction

    read(IOPAR,*) zero_trace_inversion,double_couple_inversion,lambda
    read(IOPAR,*) write_new_syn ! extension .new

    close(IOPAR)

    if (DEBUG) then
       write(*,*)
       write(*,*) ' ==  ==  ==  == INVERSION.PAR Summary ==  ==  == ='
       write(*,*) 'original cmt file: ', trim(cmt_file)
       write(*,*) 'updated cmt file: ', trim(new_cmt_file)
       write(*,*)
       if (npar == 6) then
          write(*,*) '6 parameter inversion for moment tensor only'
       else if (npar == 7) then
          write(*,*) '7 parameter inversion for moment tensor and depth'
       else if (npar == 9) then
          write(*,*) '9 parameter inversion for moment tensor and location'
       else if (npar == 10) then
          write(*,*) '10 parameter inversion for moment+location+tshift'
       else if (npar == 11) then
          write(*,*) '11 parameter inversion for moment+location+tshift+hdur'
       else
          stop 'Number of inversion parameters can only be [6,7,9,10,11]'
       endif
       write(*,'(a,3g15.5)') ' delta for derivatives: ', ddelta,ddepth,dmoment
       write(*,*)
       write(*,*) 'take input from FLEXWIN: ', trim(flexwin_out_file), ' ....'
       write(*,*)
       if (weigh_data_files) then
          if (read_weight) then
            write(*,*) 'Weight of data will be read from input file'
          else
            write(*,*) 'Weighing data files according to ...'
            write(*,'(a,3g15.5)') '  Z, R, T comp (lin) = ', comp_z_weight, comp_t_weight, comp_r_weight
            write(*,'(a,g15.5)') '  azimuth bin (exp) = ', az_exp_weight
            write(*,'(a,3g15.5)') '  Pnl/R/S dist (exp) = ', pnl_dist_weight, rayleigh_dist_weight, love_dist_weight
          endif
       else
          write(*,*) 'no weighing of data'
       endif
       if (station_correction) then
          write(*,*) 'shift data to aligh with synthetics before inversion'
       else
          write(*,*) 'data is NOT shifted before inversion'
       endif
       write(*,*)
       write(*,*) 'Inversion schemes ...'
       if (double_couple_inversion .and. .not. zero_trace_inversion) &
            stop 'check zero_trace before double_couple inversion'
       if (double_couple_inversion) then
          write(*,*) '  invert for double couple source'
       else if (zero_trace_inversion) then
          write(*,*) '  invert for zero trace source'
       else
          write(*,*) '  invert for full moment tensor source'
       endif
       write(*,'(a,g15.5)') '  inversion damping lambda = ', lambda
       write(*,*) ' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '
       write(*,*)
    endif

  end subroutine set_parameters


!============================================================

  subroutine setup_data_weights

    integer :: ios,i,j
    character*8, dimension(NRECMAX) :: kstnm,knetwk,kcmpnm
    real, dimension(NRECMAX) :: azimuth, dist_deg, dist_km
    real :: tstart, tend, tjunk
    logical :: lexd, lexs


    open(IOWIN,file=trim(flexwin_out_file),iostat=ios)
    if (ios /= 0) stop 'Flexwin output file can not be found'

    read(IOWIN,*) nfiles
    if (nfiles > NRECMAX) stop 'Increase NRECMAX limit'
    if (nfiles < 0) stop 'Check nfiles < 0'

    nwin_total = 0
    do i = 1, nfiles
       read(IOWIN,'(a)') data_file
       read(IOWIN,'(a)') syn_file
       read(IOWIN,*) nwins(i)
       if (nwins(i) < 0) stop 'Check nwins(i) '
!       if (DEBUG) print *, trim(data_file), ' ', trim(syn_file)
       inquire(file=data_file,exist=lexd)
       inquire(file=syn_file,exist=lexs)
       if (.not. (lexd .and. lexs)) then
          write(*,*) 'Check data and syn file ', trim(data_file), ' and ', trim(syn_file); stop
       endif
       do j = 1, nwins(i)
          if (read_weight) then
            ! LQY: tjunk: time shift for data/syn within the window file is not
            ! used at this stage,
            ! it may be used in the future to replace local corr subroutine
            read(IOWIN,*,iostat=ios) tstart, tend, tjunk, data_weights(nwin_total+j)
            if (ios /= 0) stop 'ts, te, tshift, weight are expected!'
          else
            read(IOWIN,*) tstart, tend
          endif
       enddo
       nwin_total = nwin_total + nwins(i)

       call read_sac_info(syn_file, data_file, &
            kstnm(i),kcmpnm(i),knetwk(i),azimuth(i),dist_deg(i),dist_km(i))
    enddo

    if (DEBUG) write(*,*) 'Total number of windows = ', nwin_total
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

!================================================================

  subroutine setup_matrix(A,b,npar)

    integer, intent(in) :: npar
    real*8, intent(out) :: A(npar,npar), b(npar)

    integer :: ios,nf,nw,nwint,i,j
    real :: tstart, tend
    real*8 :: A1(npar,npar), b1(npar)

    open(IOWIN,file=trim(flexwin_out_file),iostat=ios)
    if (ios /= 0) stop 'Flexwin output file can not be found'

    read(IOWIN,*) nf

    A = 0.
    b = 0.
    nwint = 0
    do i = 1, nfiles
       read(IOWIN,'(a)') data_file
       read(IOWIN,'(a)') syn_file
       read(IOWIN,*) nw
       do j = 1, nw
          read(IOWIN,*) tstart, tend
          nwint = nwint + 1
          call compute_A_b(syn_file,data_file,data_weights(nwint),tstart,tend,A1,b1,npar)
          A = A + A1
          b = b + b1
       enddo
    enddo
    close(IOWIN)
    if (nwint /= nwin_total) stop 'Error counting nwin_total'

    if (DEBUG) then
       ! write matrices
       write(*,'(/,a)') ' Inversion matrix A is as follows...'
       do i=1,npar
          write(*,'(12e12.3)') (sngl(A(j,i)),j=1,npar)
       enddo
       write(*,'(/,a)') ' RHS vector b is as follows...'
       write(*,'(12e12.3)') (sngl(b(j)),j=1,npar)
    endif


  end subroutine setup_matrix

!================================================================

  subroutine invert_cmt(A,b,dm,npar)

    integer, intent(in) :: npar
    real*8 , intent(inout) :: A(npar,npar), b(npar)
    real*8 , intent(out) :: dm(npar)

    real*8 :: old_par(npar),new_par(npar)
    real*8 :: trace, max_row
    logical :: linear_inversion,singular
    integer :: na, niter, i
    real*8 :: xout(NPARMAX+2), AA(NPARMAX+2,NPARMAX+2), bb(NPARMAX+2)
    real*8 :: m1(npar),lam(2), mstart(npar)
    integer, parameter :: NMAX_NL_ITER = 10

    ! do we really need the extra scaling??
    old_par = cmt_par(1:npar)/SCALE_PAR(1:npar)

    ! first scale the input A matrix and obtain its trace (pre-conditioning)
    trace = 0.0d0
    do i = 1 , npar
       max_row = maxval(abs(A(i,1:npar))) ! maximum value in row i
       A(i,1:npar) = A(i,1:npar) / max_row
       b(i) = b(i) / max_row
       trace = trace + A(i,i)
    enddo

    ! setup matrix size and solver type
    if (double_couple_inversion) then
       linear_inversion = .false.
       na = npar+2
    else if (zero_trace_inversion) then
       linear_inversion = .true.
       na = npar+1
    else
       linear_inversion = .true.
       na = npar
    endif

    ! adding damping, mostly not needed for over-determined lsq problem
    do i = 1, npar
       A(i,i) = A(i,i) + trace * lambda
    enddo

    AA = 0.
    bb = 0.
    if (linear_inversion) then

       ! if only invert for moment tensor with zero trace constraint,
       ! it is a linear inversion which may depend on npar for output
       ! notice that this is a slight variation from (A4)
       print *, '   Linear inversion ...'
       AA(1:npar,1:npar) = A   ! setup input matrices
       bb(1:npar) = b
       if (zero_trace_inversion) then   ! zero-trace inversion case
          bb(na) = - sum(old_par(1:3))
          AA(1:6,na) = (/1,1,1,0,0,0/)
          AA(na,1:6) = (/1,1,1,0,0,0/)
          AA(na,na) = 0.
       endif

       call Gaussian_elimination(AA,na,bb,xout,singular)
       if (singular) stop 'AA matrix is near singular'
       dm = xout(1:npar)
       new_par = old_par(1:npar) + dm(1:npar)

    else

       ! if invert for moment tensor with double couple constraints
       ! setup starting solution, solve directly for moment instead
       ! of dm, exact implementation of (A16)
       print *, 'Non-linear inversion...'
       mstart(1:npar) = old_par(1:npar)
       m1(1:npar) = mstart(1:npar)
       lam(1:2) = 0 ! guess the LMP start with zeros

       do i = 1, NMAX_NL_ITER
          call get_f_df(npar,A,b,m1,lam,mstart,AA,bb)
          bb = -bb
          call Gaussian_elimination(AA,na,bb,xout,singular)
          if (singular) stop 'Error Gaussian elimination'
          m1(1:npar) = m1(1:npar) + xout(1:npar)
          lam = lam + xout(npar+1:na)
       enddo
       dm = m1(1:npar) - mstart(1:npar)
       new_par(1:npar) = m1(1:npar)

    endif

    new_cmt_par(1:NPARMAX) = cmt_par(1:NPARMAX)
    new_cmt_par(1:npar) = new_par(1:npar) * SCALE_PAR(1:npar)

    if (DEBUG) then
       write(*,'(/,a)') 'Gaussian Elimination gives dm(1:npar) ...'
       write(*,'(11e15.3)') sngl(dm(1:npar))
       write(*,'(/,a)') ' Output scaled parameters are ...'
       write(*,*)
       write(*,'(11e15.3,/)') sngl(new_par(1:npar))
       write(*,'(/,a,i2,a,/)')' ***********Inversion Result Table(',npar,' pars)************'
       print *, 'PAR_RAMETERS       Old_CMT        New_CMT '
       print *
       do i = 1, NPARMAX
          write(*,'(1x,a,3g20.8)') par_name(i),cmt_par(i),new_cmt_par(i)
       enddo
    endif

  end subroutine invert_cmt


!============================================================
! to match calc_var and create_adj_src, a hanning taper needs to be applied
  subroutine variance_reduction(dm,npar)

    integer,intent(in) :: npar
    real*8,intent(in) :: dm(npar)

    integer :: ios, nf, nwint, nw,  is, ie, i, j, npts, ishift, ishift_new
    real :: b, dt, tstart, tend
    real :: tshift, cc, dlna, v1, v2, tshift_new,cc_new,dlna_new,d1,d2
    real :: var_all, var_all_new
    integer :: istart, iend, istart_d, iend_d, ii
    integer :: istart_n, iend_n, istart_dn, iend_dn
    real, dimension(NDATAMAX) :: taper


    open(IOWIN,file=trim(flexwin_out_file),iostat=ios)
    if (ios /= 0) stop 'Flexwin output file can not be found'

    open(IOINV,file='cmt3d_flexwin.out')

    read(IOWIN,*) nf
    write(IOINV,*) nf

    nwint = 0; var_all = 0.; var_all_new = 0.
    do i = 1, nf
       read(IOWIN,'(a)') data_file
       read(IOWIN,'(a)') syn_file
       write(IOINV,*) trim(data_file)
       write(IOINV,*) trim(syn_file)

       call compute_new_syn(data_file,syn_file,npts,b,dt,dm)
       read(IOWIN,*) nw
       write(IOINV,*) nw
       do j = 1, nw
          nwint = nwint + 1
          read(IOWIN,*) tstart, tend
          is=max(floor((tstart-b)/dt),1)
          ie=min(ceiling((tend-b)/dt),npts)

          if (station_correction) then  !!! needs update?
             call calc_criteria(data_sngl,syn_sngl,npts,is,ie,ishift,cc,dlna)
             call calc_criteria(data_sngl,new_syn_sngl,npts,is,ie,ishift_new,cc_new,dlna_new)
             istart_d=max(1,is+ishift); iend_d=min(npts,ie+ishift)
             istart_dn=max(1,is+ishift_new); iend_dn=min(npts,ie+ishift_new)
             istart=istart_d-ishift; iend=iend_d-ishift
             istart_n=istart_dn-ishift_new; iend_n=iend_dn-ishift_new
          else
             istart_d=is; istart=is; iend_d=ie; iend=ie
             istart_dn=is; istart_n=is; iend_dn=ie; iend_n=ie
          endif

          ! hanning taper
          do ii = istart,iend
             taper(ii) = 0.5 * (1-cos(2*pi*(ii-istart)/(iend-istart)))
          enddo
          v1=sum(taper(istart:iend)*(syn_sngl(istart:iend)-data_sngl(istart_d:iend_d))**2)
          v2=sum(taper(istart:iend)*(new_syn_sngl(istart_n:iend_n)-data_sngl(istart_dn:iend_dn))**2)
          d1=sum(taper(istart:iend)*data_sngl(istart_d:iend_d)**2)
          d2=sum(taper(istart:iend)*data_sngl(istart_dn:iend_dn)**2)
          var_all = var_all + 0.5*v1*data_weights(nwint)*dt
          var_all_new = var_all_new + 0.5*v2*data_weights(nwint)*dt

          write(IOINV,'(3g15.5)') tstart,tend,data_weights(nwint)
          write(IOINV,'(3x,4g15.5)') ishift*dt, cc, dlna, v1/d1
          write(IOINV,'(3x,4g15.5)') ishift_new*dt, cc_new, dlna_new, v2/d2
       enddo
    enddo
    close(IOWIN)
    close(IOINV)

    if (nwint /= nwin_total) stop 'Check nwin_total in variance reduction'
    print *, 'Total Variance reduced from ', var_all, ' to ', var_all_new, &
         ' = ', (var_all-var_all_new)/var_all*100, ' %'

  end subroutine variance_reduction

!===============================================================

end module cmt3d_sub

