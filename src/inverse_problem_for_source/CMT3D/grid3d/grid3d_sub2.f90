module grid3d_sub2

  use grid3d_constants
  use grid3d_variables
  use grid3d_sub3

  implicit none

contains

  ! =================================================================

  subroutine read_sac_info(file_s, file_o, kstnm, kcmpnm, knetwk, &
       azimuth, dist_deg, dist_km)

    character (len=*),intent(in) :: file_s, file_o
    character (len=*),intent(out) :: kstnm, kcmpnm, knetwk
    real,intent(out) :: azimuth, dist_deg, dist_km

    real :: b1, dt1, b2, dt2, b, dt
    real :: evla, evlo, stla, stlo, evdp, backazimuth
    integer :: npts1, npts2, npts, nerr
    character(len=1) :: cmp

    ! read synthetic
    call rsac1(file_s,syn,npts1,b1,dt1,NDATAMAX,nerr)
    if (nerr /= 0) stop ' Error reading synthetic file'

    ! read observed
    call rsac1(file_o,data,npts2,b2,dt2,NDATAMAX,nerr)
    if (nerr /= 0) stop ' Error reading data file '

    if (npts1 > NDATAMAX .or. npts2 > NDATAMAX) stop 'Error: change NDATAMAX'

    ! check sample rates are equal
    if (abs(dt1-dt2) > EPS5) stop 'Sampling rates differ, program stop !!!'
    dt = dble(dt1)
    !write(*,*)'sampling rate dt=',dt

    ! check start times are equal (LQY: a difference of dt is probably too big)
    if (abs(b1-b2) > dt) stop ' start times differ, program stop !!!'
    b=dble(b1)

    ! set global npts to the npts of the shortest seismogram
    npts = min(npts1, npts2)

    ! read event and station header parameters from observation file
    call getfhv('evla', evla, nerr)
    call getfhv('evlo', evlo, nerr)
    call getfhv('stla', stla, nerr)
    call getfhv('stlo', stlo, nerr)
    call getfhv('evdp', evdp, nerr)
    call getkhv('kstnm', kstnm, nerr)
    call getkhv('kcmpnm', kcmpnm, nerr)
    call getkhv('knetwk', knetwk, nerr)

    cmp=kcmpnm(3:3)

    if ((cmp /= 'Z') .and. (cmp /= 'T') .and. (cmp /= 'R')) &
         stop 'We only deal with Z, R, and T components at the moment'

    ! calculate distances and azimuths (azimuth in degrees)
    call distaz(evla,evlo,stla,stlo,azimuth,backazimuth,dist_deg,dist_km)

    if (azimuth < 0 .or. azimuth > 360) stop 'check azimuth in (0,360)'


  end subroutine read_sac_info

  !-------------------------------------------------------------------

  subroutine compute_data_weights(kcmpnm,azimuth,dist_km,data_weights)

    character*8, dimension(:),intent(in) :: kcmpnm
    real, dimension(:), intent(in) :: azimuth, dist_km
    real, dimension(:), intent(out) :: data_weights

    real :: daz
    integer :: naz(NREGIONS), nwint, i, j, k
    character(len=8) :: comp_name
    real :: cmp_weight(NWINMAX), dist_exp_weight(NWINMAX),max_data_weight


    daz = 360./NREGIONS
    naz = 1 ! start with a water level of 1

    nwint = 0
    do i = 1, nfiles
       do j = 1, nwins(i)
          nwint = nwint + 1

          ! component weights
          comp_name=kcmpnm(i)
          if (comp_name(3:3) == 'Z') then
             cmp_weight(nwint)=comp_z_weight
          else if (comp_name(3:3) == 'R') then
             cmp_weight(nwint)=comp_r_weight
          else if (comp_name(3:3) == 'T') then
             cmp_weight(nwint)=comp_t_weight
          endif

          ! dist weights
          if (comp_name(3:3) == 'T') then
             dist_exp_weight(nwint) = love_dist_weight
          else
             ! first window out of two windows on R and Z components is Pnl
             if (nwins(i) > 1 .and. j == 1) then
                dist_exp_weight(nwint) = pnl_dist_weight
             else
                dist_exp_weight(nwint) = rayleigh_dist_weight
             endif
          endif
       enddo

       ! azimuth counts
       k = floor(azimuth(i)/daz) + 1
       if (k < 0 .or. k > NREGIONS) stop 'Error binning azimuth'
       naz(k) = naz(k) + 1
    enddo

    if (nwint /= nwin_total) stop 'Error counting number of windows'
    if (DEBUG) write(*,*) ' DEBUG : Number of files in az bins: ', naz-1

    ! assemble data weights
    nwint = 0
    do i = 1, nfiles
       do j = 1, nwins(i)
          nwint = nwint + 1
          k = floor(azimuth(i)/daz) + 1
          data_weights(nwint) =  cmp_weight(nwint) &
               * ( (dist_km(i)/REF_DIST) ** dist_exp_weight(nwint) ) &
               / ( naz(k) ** az_exp_weight)
       enddo
    enddo

    ! normalization of data weights
    max_data_weight = maxval(data_weights(1:nwin_total))

    data_weights(1:nwin_total) = data_weights(1:nwin_total)/max_data_weight

  end subroutine compute_data_weights

  !---------------------------------------------------------------------------

  subroutine compute_mij_from_sdr(s_strike,d_strike,n_strike, &
       s_dip,d_dip,n_dip,s_rake,d_rake,n_rake,s_mw,d_mw,n_mw,mij)

    real, intent(in) :: s_strike,d_strike,s_dip,d_dip,s_rake,d_rake
    real, intent(in) :: s_mw, d_mw
    integer, intent(in) :: n_strike, n_dip, n_rake, n_mw
    real, dimension(:,:), intent(out) :: mij

    real :: strike,dip,rake,moment,mw
    integer :: im,ir,id,it,j

    ! loop over (mw,rake,dip,strike) to pre-calculate mij's
    j=0
    do im = 1, n_mw
       mw = s_mw + (im-1) * d_mw
       ! equation (9.45) from Modern Global Seismology
       moment =10 ** ((mw + 10.73) * 1.5)
       do ir = 1, n_rake
          rake = s_rake + (ir-1) * d_rake
          do id = 1, n_dip
             dip = s_dip + (id-1) * d_dip
             do it = 1, n_strike
                strike = s_strike + (it-1) * d_strike
                j=j+1
                call sdr2moment(strike,dip,rake,moment, &
                     mij(1,j),mij(2,j),mij(3,j),mij(4,j),mij(5,j),mij(6,j))
             enddo
          enddo
       enddo
    enddo
    if (j /= n_mw * n_rake * n_dip * n_strike) stop 'Error counting nj'

  end subroutine compute_mij_from_sdr


  !---------------------------------------------------------------------------

  subroutine add_misfit(data_file,syn_file, &
       data_weight,tstart,tend,nw, &
       mij,misfit,n_total)

    character(len=150),intent(in) :: data_file, syn_file
    integer,intent(in) :: nw, n_total
    real,intent(in) :: data_weight(nw),tstart(nw),tend(nw),mij(:,:)
    real,intent(inout) :: misfit(:)

    integer :: npts,npts2,nerr,is(nw),ie(nw)
    integer :: iw,it,id,ir,im,i,j,ishift,iss,iee,isd,ied,ishift_max
    real :: b,dt,b2,dt2,cc,t1(nw),tt
    character(len=150) :: dsyn_file

    ! read observed
    call rsac1(data_file,data,npts,b,dt,NDATAMAX,nerr)

    ! read synthetics
    do i = 1, NM
       dsyn_file = trim(syn_file)//'.'//trim(PAR_NAME(i))
       npts2 = 0
       call rsac1(dsyn_file,dsyn(i,:), npts2,b2,dt2,NDATAMAX,nerr)
       if (abs(b-b2) > dt .or. abs(dt-dt2) > EPS5 .or. npts /= npts2) &
            stop 'Error reading dsyn-file'
    enddo

    do i = 1, nw
       is(i)=max(floor((tstart(i)-b)/dt),1)
       ie(i)=min(ceiling((tend(i)-b)/dt),npts)
       t1(i) = dt * data_weight(i)
    enddo

    ! core computation
    ishift=0; ishift_max=nint(tshift_max/dt)
    do j = 1, n_total
       syn(1:npts) =  matmul(mij(:,j),dsyn(:,1:npts)) / dmoment
       do iw = 1, nw
          tt=t1(iw); iss=is(iw);iee=ie(iw)
          if (station_correction) then
             call xcorr_calc(data,syn,npts,iss,iee,ishift,cc,ishift_max)
             isd=max(1,iss+ishift); ied=min(npts,iee+ishift)
             iss=isd-ishift; iee=ied-ishift
          else
             isd=iss; ied=iee
          endif
          misfit(j) = misfit(j)+tt*sum((syn(iss:iee)-data(isd:ied))**2)
       enddo
    enddo

  end subroutine add_misfit

  ! --------------------------------------------------------------

  subroutine select_best_solution(icalc,misfit,n_strike,n_dip,n_rake,n_mw, &
       s_strike,s_dip,s_rake,s_mw,strike,dip,rake,mw)

    ! LQY: maybe one day I'll interpolate it with a quadratic function
    ! and search for the minimum of the quadratic function, instead of
    ! taking directly the minimum value from grid search:
    ! chi = A(s-s_m)**2 + B(d-d_m)**2 + C(r-r_m)**2 + D
    ! one little issue here is that we have 7 unknowns and 8 equations,
    ! and the simplest solution is to eliminate one point that is the
    ! fartheset from the average value on the cube.

    real,intent(inout) :: misfit(:)
    integer,intent(in) :: icalc,n_strike,n_dip,n_rake,n_mw
    real,intent(in) :: s_strike,s_dip,s_rake,s_mw
    real,intent(out) :: strike,dip,rake,mw

    integer :: jmina(1),jmin,jm,jr,jd,js,n_total,im,mm,mm2,mm3,is,id,ir
    character(len=150) :: filename
    real :: min_misfit

    n_total=n_strike*n_dip*n_rake*n_mw
    mm = n_strike*n_dip*n_rake
    mm2 = n_strike*n_dip
    mm3 = n_strike

    jmina = minloc(misfit(1:n_total))
    jmin=jmina(1)-1

    jm = jmin/mm+1
    jmin = jmin-(jm-1)*mm
    jr = jmin/mm2+1
    jmin = jmin-(jr-1)*mm2
    jd = jmin/mm3+1
    js = jmin-(jd-1)*mm3+1

    strike=s_strike+(js-1)*d_strike
    dip=s_dip+(jd-1)*d_dip
    rake=s_rake+(jr-1)*d_rake
    mw=s_mw+(jm-1)*d_mw

    filename='grid3d_misfit'
    if (global_search) filename=trim(filename)//'.'//char(icalc+48)

    min_misfit=minval(misfit(1:n_total))
    print *, 'minimum misfit value = ', min_misfit
    if (abs(min_misfit) > 10*tiny(1.0)) then
       misfit = misfit/min_misfit
    else
       stop 'check why min_misfit ~ 0'
    endif

    do im = 1, n_mw
       open(IOGRD,file=trim(filename)//'.'//char(im+48),status='unknown')
       do ir = 1, n_rake
          do id = 1, n_dip
             do is = 1, n_strike
                write(IOGRD,'(4g15.3)') (s_strike+(is-1)*d_strike)*R2D, &
                     (s_dip+(id-1)*d_dip)*R2D, (s_rake+(ir-1)*d_rake)*R2D, &
                     misfit((im-1)*mm+(ir-1)*mm2+(id-1)*mm3+is)
             enddo
          enddo
       enddo
       close(IOGRD)
    enddo

  end subroutine select_best_solution


end module grid3d_sub2
