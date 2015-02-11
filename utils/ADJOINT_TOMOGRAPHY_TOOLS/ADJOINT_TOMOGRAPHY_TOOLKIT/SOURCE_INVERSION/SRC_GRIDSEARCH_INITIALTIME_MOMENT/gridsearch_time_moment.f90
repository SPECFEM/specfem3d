program xgrid_search_time_moment

  implicit none

  integer,parameter:: NDIM=80000

  character(len=256) :: cmt_file, new_cmt_file
  character(len=256) :: output_minmax, output_misfit
  character(len=256) :: flexwin_out_file
  character(len=256) :: data_file, syn_file
  character(len=150) :: string,bandpass,cmp
  character(len=150) :: body_bandpass,surf_bandpass,criteria
  real ,dimension(:),allocatable :: misfit,t0,m0
  real,dimension(NDIM):: data,syn,syn_m0
  real,dimension(NDIM):: data_win,syn_win
  real :: tshift, dlnA, cc_max
  real::s_m0,e_m0,d_m0,s_t0,e_t0,d_t0,t00,m00,m0_best,t0_best,maxmisfit,t_cmt,t_cmt_new
  integer:: i,j,k,n_m0,n_t0,n_total,i_total,cc_shift
  real:: misfit_tmp,b,b2,dt,dt2,tstart,tend,tstart0,tend0,weight
  integer:: nf,nwin,nerr,ios,lstr,npts,npts2
  integer:: is,ie,ishift,isd,ied,iss,iee,imina(1),imin
  integer:: bd_z,bd_r,bd_t,sw_z,sw_r,sw_t
  real:: deltaT,deltaA,fact_am,fact_tt,deltaT0
  real:: linesearch_mbest
  integer:: istart,iend

  !----------------------------------------------------------------
  !- read in input parameters
  !----------------------------------------------------------------
  read(*,'(a)') cmt_file          ! read in old cmt file
  read(*,'(a)') new_cmt_file      ! read in new cmt file
  read(*,'(a)') flexwin_out_file  ! read in measurement window selected from flexwin
  read(*,'(a)') output_minmax     ! output selected original time and scale moment
  read(*,'(a)') output_misfit     ! output 2D misfit function for delta_t0 and delta_M0
  read(*,*) s_m0,e_m0,d_m0        ! starting, end and spacing for scale moment searching
  read(*,*) s_t0,e_t0,d_t0        ! starting, end and spacing for original time searching
  read(*,*) fact_am               ! weighting for amplitude misfit function
  read(*,*) fact_tt               ! weighting for original time misfit function
  read(*,*) linesearch_mbest      ! perturbation for scale moment pick
  read(*,'(a)') body_bandpass     ! extension for body wave bandpass, e.g.,T015_040
  read(*,'(a)') surf_bandpass     ! extension for surface wave bandpass, e.g.,T040_100
  read(*,'(a)') criteria          ! misfit function, waveform, phase_amplitude

  !---------------------------------------------------------
  !- setup searching space
  !---------------------------------------------------------
  n_m0=nint((e_m0-s_m0)/d_m0)+1
  n_t0=nint((e_t0-s_t0)/d_t0)+1
  n_total=n_m0*n_t0

  allocate(misfit(n_total))
  allocate(t0(n_total))
  allocate(m0(n_total))
  misfit=0.0
  t0=0.0
  m0=0.0

  i_total=0
  do i=1,n_m0
     do j=1,n_t0
        i_total=i_total+1
        t0(i_total)=s_t0+d_t0*(j-1)
        m0(i_total)=s_m0+d_m0*(i-1)
     enddo
  enddo

  if ( i_total /= n_total )  stop 'i_total not equal n_total, check input paramters'
  write(*,*) i_total,n_total


  !------------------------------------------------------------
  !- read in measurement window from flexwin and calculate misfit function
  !------------------------------------------------------------
  open(1001,file=trim(flexwin_out_file),status='old')
  read(1001,*) nf,bd_z,bd_r,bd_t,sw_z,sw_r,sw_t
  do j=1,nf
     write(*,*) j,nf
     read(1001,'(a)') data_file
     read(1001,'(a)') syn_file

     call rsac1(data_file,data,npts,b,dt,NDIM,nerr)
     if (nerr/=0) stop 'error reading data'

     call rsac1(syn_file,syn,npts2,b2,dt2,NDIM,nerr)
     if (nerr/=0) stop 'error reading synthetics'

     lstr=len_trim(syn_file)
     bandpass=syn_file(lstr-7:lstr)
     cmp=syn_file(lstr-19:lstr-17)

     ! decide weight factors
     if ( trim(bandpass) == trim(body_bandpass) .and. trim(cmp) == 'LHZ' ) then
        if ( bd_z /= 0 ) then
           weight=1.0/bd_z
        else
           weight=1.0
        endif
     else if ( trim(bandpass) == trim(body_bandpass) .and. trim(cmp) == 'LHR' ) then
        if ( bd_r /= 0 ) then
           weight=1.0/bd_r
        else
           weight=1.0
        endif
     else if ( trim(bandpass) == trim(body_bandpass) .and. trim(cmp) == 'LHT' ) then
        if ( bd_t /= 0 ) then
           weight=1.0/bd_t
        else
           weight=1.0
        endif
     else if ( trim(bandpass) == trim(surf_bandpass) .and. trim(cmp) == 'LHZ' ) then
        if ( sw_z /= 0 ) then
           weight=1.0/sw_z
        else
           weight=1.0
        endif
     else if ( trim(bandpass) == trim(surf_bandpass) .and. trim(cmp) == 'LHR' ) then
        if ( sw_r /= 0 ) then
           weight=1.0/sw_r
        else
           weight=1.0
        endif
     else if ( trim(bandpass) == trim(surf_bandpass) .and. trim(cmp) == 'LHT' ) then
        if ( sw_t /= 0) then
           weight=1.0/sw_t
        else
           weight=1.0
        endif
     else
        stop 'wrong bandpass and component, check input parameters'
     endif


     read(1001,*) nwin
     do k=1,nwin
        read(1001,*) tstart,tend
        is=max(floor((tstart-b)/dt),1)
        ie=min(ceiling((tend-b)/dt),npts)

        data_win(:)=0.0
        syn_win(:)=0.0
        data_win(is:ie)=data(is:ie)
        syn_win(is:ie)=syn(is:ie)
        call xcorr_calc(data_win,syn_win,npts,is,ie,cc_shift,cc_max)
        deltaT0=cc_shift*dt


        do i=1,n_total

           t00=t0(i)
           m00=m0(i)
           ishift=nint(t00/dt)

           isd=is
           ied=ie

           iss=isd-ishift
           iee=ied-ishift


           if ( trim(criteria) == 'waveform' ) then
              ! full waveform misfit function
              misfit_tmp=weight*sum((m00*syn(iss:iee)-data(isd:ied))**2)/sum(sqrt(syn(iss:iee)**2)*sqrt(data(isd:ied)**2))
           else if (trim(criteria) == 'phase_amplitude' ) then
              ! phase and amplitude misfit function
              deltaT=deltaT0-t00
              syn_m0=syn*m00
              deltaA=0.5*log(sum(data(isd:ied)*data(isd:ied))/sum(syn_m0(iss:iee)*syn_m0(iss:iee)))

              misfit_tmp=weight*(fact_tt*(deltaT**2)+fact_am*(deltaA**2))
           else
              stop 'wrong criteria for misfit function, check input parameters'
           endif

           ! total misfit function
           misfit(i)=misfit(i)+misfit_tmp

         enddo
      enddo
  enddo
  close(1001)


  !------------------------------------------------------------
  !- select minimum misfit function and write out new cmt file
  !------------------------------------------------------------
  imina=minloc(misfit(1:n_total))
  imin=imina(1)
  m0_best=m0(imin)
  t0_best=t0(imin)

  maxmisfit=maxval(misfit(1:n_total))

  open(1004,file=trim(cmt_file),status='old',iostat=ios)
  if (ios/=0) stop 'Error openning CMT file'
  do while (ios == 0)
     read(1004,'(a)',iostat=ios) string
     lstr=len_trim(string)
     if (string(1:10) == 'time shift') then
        read(string(12:lstr),*) t_cmt
     endif
  enddo
  close(1004)

  t_cmt_new=t_cmt+t0_best

  ! several criteria to negelect selected parameters
  ! 1. if number of window smaller than 10
  ! 2. if new original time is smaller than 0.0
  ! 3. for scale moment, only use perturbation linesearch_mbest since its large covariance
  if (nf < 10 ) then
     t0_best=0.0
     m0_best=1.0
  endif

  if ( t_cmt_new < 0.0 ) then
     t0_best=0.0
     m0_best=1.0
  endif
  m0_best=1.0+(m0_best-1.0)*linesearch_mbest



  write(*,*) 'DONE'
  write(*,*) 'FIND BEST T0:',t0_best
  write(*,*) 'FIND BEST M0:',m0_best
  write(*,*) 'MIN MISFIT:',misfit(imin)



  open(1002,file=trim(output_minmax),status='unknown')
  write(1002,*) t0_best,m0_best,misfit(imin)
  close(1002)
  open(1003,file=trim(output_misfit),status='unknown')
  do i = 1,n_total
     write(1003,*) t0(i),m0(i),misfit(i)/maxmisfit
  enddo
  close(1003)

  call write_new_cmtsolution(cmt_file,trim(new_cmt_file),t0_best,m0_best)

  write(*,*) 'SUCESSIVEFULLY'

  deallocate(misfit)
  deallocate(t0)
  deallocate(m0)

end program xgrid_search_time_moment


subroutine xcorr_calc(d,s,npts,i1,i2,ishift,cc_max)

  ! inputs:
  ! s(npts) = synthetic
  ! d(npts) = data (or observed)
  ! i1, i2 = start and stop indexes of window within s and d

  real, dimension(*), intent(in) :: s,d
  integer, intent(in) :: npts, i1, i2

  ! outputs:
  ! ishift = index lag (d-s) for max cross correlation
  ! cc_max = maximum of cross correlation (normalised by sqrt(synthetic*data))
  integer, intent(out) :: ishift
  real, intent(out) :: cc_max

  ! local variables
  integer :: nlen
  integer :: i_left, i_right, i, j, id_left, id_right
  real :: cc, norm, norm_s, fout

  ! initialise shift and cross correlation to zero
  ishift = 0
  cc_max = 0.0

  if (i1<1 .or. i1>i2 .or. i2>npts) then
    write(*,*) 'Error with window limits: i1, i2, npts ', i1, i2, npts
    return
  endif

  ! length of window (number of points, including ends)
  nlen = i2 - i1 + 1

  ! power of synthetic signal in window
  norm_s = sqrt(sum(s(i1:i2)*s(i1:i2)))

  ! left and right limits of index (time) shift search
  ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and CC.
  !       If fout=0.5, then it looks outside by a time of 0.5*(window width).
  !       Perhaps fout should be a parameter, or it should be tied to the max
  !          allowable time-shift specified by the user.
  !       However, it does not matter as much if the data and synthetics are
  !          zeroed outside the windows, as currently done in calc_criteria.
  fout = 0.5
  i_left = -1*int(fout*nlen)
  i_right = int(fout*nlen)


  ! i is the index to shift to be applied to DATA (d)
  do i = i_left, i_right

    ! normalization factor varies as you take different windows of d
    id_left = max(1,i1+i)      ! left index for data window
    id_right = min(npts,i2+i)  ! right index for data window
    norm = norm_s * sqrt(sum(d(id_left:id_right)*(d(id_left:id_right))))

    ! cc as a function of i
    cc = 0.
    do j = i1, i2   ! loop over full window length
      if((j+i)>=1 .and. (j+i)<=npts) cc = cc + s(j)*d(j+i)  ! d is shifted by i
    enddo
    cc = cc/norm

    ! keeping cc-max only
    if (cc > cc_max) then
      cc_max = cc
      ishift = i
    endif
  enddo

  ! EXAMPLE: consider the following indexing:
  ! Two records are from 1 to 100, window is i1=20 to i2=41.
  !    --> nlen = 22, i_left = -11, i_right = 11
  !    i   i1+i   i2+i  id_left  id_right
  !  -11     9     30      9        30
  !   -5    15     36     15        36
  !    0    20     41     20        41    <== ORIGINAL WINDOW
  !    5    25     46     25        46
  !   10    31     52     31        52

end subroutine xcorr_calc


subroutine write_new_cmtsolution(cmt_file,new_cmt_file,t0_best,m0_best)

  character(len=*),intent(in):: cmt_file,new_cmt_file
  real,intent(in):: t0_best,m0_best
  integer,parameter:: IOCMT=1009

  character(len=150):: pde_time,event_name,str_tshift,str_hdur,str_lat,str_lon,str_dep
  real:: t_cmt, t_cmt_new
  real,dimension(6):: moment_tensor,moment_tensor_new
  integer:: ios,lstr
  character(len=150) :: string

  ! read basic information
  open(IOCMT,file=trim(cmt_file),status='old')
  read(IOCMT,'(a)') pde_time
  read(IOCMT,'(a)') event_name
  read(IOCMT,'(a)') str_tshift
  read(IOCMT,'(a)') str_hdur
  read(IOCMT,'(a)') str_lat
  read(IOCMT,'(a)') str_lon
  read(IOCMT,'(a)') str_dep
  close(IOCMT)

  ! read tshift and moment tensor
  open(IOCMT,file=trim(cmt_file),status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening CMT file'
  do while (ios==0)
        read(IOCMT,'(a)',iostat=ios) string
        lstr=len_trim(string)

        if (string(1:10) == 'time shift') then
                read(string(12:lstr),*) t_cmt
        else if (string(1:3) == 'Mrr') then
                read(string(5:lstr),*) moment_tensor(1)
        else if (string(1:3) == 'Mtt') then
                read(string(5:lstr),*) moment_tensor(2)
        else if (string(1:3) == 'Mpp') then
                read(string(5:lstr),*) moment_tensor(3)
        else if (string(1:3) == 'Mrt') then
                read(string(5:lstr),*) moment_tensor(4)
        else if (string(1:3) == 'Mrp') then
                read(string(5:lstr),*) moment_tensor(5)
        else if (string(1:3) == 'Mtp') then
                read(string(5:lstr),*) moment_tensor(6)
        endif
  enddo
  close(IOCMT)

  ! calculate new t0 and Mij
  !t_cmt_new=t_cmt
  t_cmt_new=t_cmt+t0_best
  moment_tensor_new(:)=moment_tensor(:)*m0_best

  ! write out new CMTSOLUTION
  open(IOCMT,file=trim(new_cmt_file),status='unknown')
  write(IOCMT,'(a)') trim(pde_time)
  write(IOCMT,'(a)') trim(event_name)
  write(IOCMT,'(a,f12.4)') 'time shift:',t_cmt_new
  write(IOCMT,'(a)') trim(str_hdur)
  write(IOCMT,'(a)') trim(str_lat)
  write(IOCMT,'(a)') trim(str_lon)
  write(IOCMT,'(a)') trim(str_dep)

  write(IOCMT,'(a,g19.6)') 'Mrr:',moment_tensor_new(1)
  write(IOCMT,'(a,g19.6)') 'Mtt:',moment_tensor_new(2)
  write(IOCMT,'(a,g19.6)') 'Mpp:',moment_tensor_new(3)
  write(IOCMT,'(a,g19.6)') 'Mrt:',moment_tensor_new(4)
  write(IOCMT,'(a,g19.6)') 'Mrp:',moment_tensor_new(5)
  write(IOCMT,'(a,g19.6)') 'Mtp:',moment_tensor_new(6)
  close(IOCMT)

end subroutine write_new_cmtsolution


