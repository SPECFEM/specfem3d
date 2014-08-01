!      
!     ---  Alessia Maggi, May 2006  ----
!        
!     loosely based on code by Carl Tape, Sept 2005
!   
!     Measurement by simple cross-correlation and adjoint source calculation
!     Both phase and amplitude information output
!
!     Input: 
!        Two sac files : synthetic and data 
!        Starting/Ending measurement window in seconds
!
!     Ouput: 
!        at discrete frequencies omega
!         amplitude:
!        (1) d(lnA) = (A_obs/A_syn - 1) =  dA/A
!
!        delay time dt(\omega) and phase dPhi(\omega):
!        (1) dt(\omega) = T_obs(\omega) - T_syn(\omega)
!
!
!	$Id: $
!        -------------------------------------------------------
!

!     MODULE FOR CONSTANTS
!     ------------------------------------------------------------------

      module xcorr_constants
        implicit none

      integer, parameter :: lnpt=14, npt=2**lnpt

!!$      double precision, parameter :: pi = 3.1415926
!!$      double precision, parameter :: twopi = 2*pi
!!$      integer, parameter ::  HANNING = 1
!!$      integer, parameter ::  HAMMING = 2
!!$      integer, parameter ::  COSINE = 3

      end module xcorr_constants


!     START MAIN SUBROUTINE
!     ------------------------------------------------------------------

      subroutine xcorr_measure(t_start,t_stop,syn,obs,npts,dt,b,&
                             tshift_after,CC_after,dlnA_after, &
                             dtau_xcorr,dlnA_xcorr,&
                             obs_win_save,syn_win_save,obs_rec,&
                             calc_adjoint,fp,fq,DEBUG)
      use xcorr_constants
      implicit none

!     declaration of subroutine parameters
!     ------------------------------------------------------------------
      double precision, intent(in) :: t_start, t_stop        ! window limits
      double precision, intent(in), dimension(*) :: syn, obs ! synthetic and observed seismograms
      double precision, intent(in) :: dt, b                  ! timestep and absolute start time of syn & obs vectors
      integer, intent(in) :: npts                ! length of syn & obs vectors
      logical, intent(in) :: DEBUG

      double precision, intent(out) :: dtau_xcorr,dlnA_xcorr           ! xcorr measurements (tshift, dlna)
      double precision, intent(out) :: tshift_after , CC_after, dlnA_after   ! for reconstructed data and data
      double precision, dimension(npt), intent(out) :: syn_win_save, obs_win_save, obs_rec

      logical, intent(in), optional :: calc_adjoint  ! flags the calculation of the adjoint source
      double precision, dimension(npt), intent(out), optional :: fp, fq ! adjoint sources (1:nlen)


!     declaration of local variables
!     ------------------------------------------------------------------
      integer :: nstart, nstop, nlen    ! start & end indexes of window + window length

      integer :: i                ! loop indexes

!     for cross correlation calculation
      integer :: ishift
      double precision :: cc_max 

!     for adjoint sources
      double precision :: Nnorm, Mnorm
      double precision, dimension (npt) :: syn_veloc, syn_accel, t_window
      double precision, dimension(npt) :: fp_nomeasure, fq_nomeasure ! adjoint sources (1:nlen)

!     functions
      real, external :: costaper

!     Obtain indexes of window limits + length of window
!     ------------------------------------------------------------------
      nstart = 1+int(floor((t_start - b) / dt))
      nstop =  1+int(ceiling((t_stop  - b) / dt))
      nlen = nstop - nstart +1
      if(nstart.lt.1.or.nstop.gt.npts) stop 'Window outside seismogram'

!     Apply boxcar window to filtered data to create windowed seismograms
!     ------------------------------------------------------------------
      obs_win_save(1:nlen) = obs(nstart:nstop)
      syn_win_save(1:nlen) = syn(nstart:nstop)
!     DEBUG: output 
      if (DEBUG) then
      open(unit=11, file='DEBUG_xcor_windows.dat')
      do i = 1, nlen
        write(11,'(5(e12.4,1x))') (i-1)*dt, syn_win_save(i), obs_win_save(i), syn(nstart-1+i), obs(nstart-1+i)
      enddo      
      close(11)
      endif



!     Create time window taper
!     ------------------------------------------------------------------
      t_window(:)=0
      do i = 1, nlen
        t_window(i) = dble(costaper(nlen,i,10))
      enddo

      call calc_criteria(obs,syn,npts,nstart,nstop,dt,dtau_xcorr,cc_max, dlnA_xcorr)
      ishift=int(dtau_xcorr/dt)

!     Apply -tshift to synthetic seismogram to reconstruct it
!     (time part of the transfer function)
!     -----------------------------------------------------
      obs_rec(:) = 0
      do i = 1, nlen
      if( (nstart-1+i-ishift) .ge. 1 .and. (nstart-1+i-ishift) .lt.npts ) then
        obs_rec(i) = syn(nstart-1+i-ishift)
      endif
      enddo



!     Finish reconstruction by applying dlnA
      obs_rec(1:nlen) = obs_rec(1:nlen)*(1+dlnA_xcorr)
      

!     calculate the quality according to Ritsema & van Heijst (2002)
!     ------------------------------------------------------------------
!     The obs_win vs the reconstructed obs
      call calc_criteria(obs_win_save,obs_rec,npt,1,nlen,dt,tshift_after,&
                     CC_after,dlnA_after)

!     ------------------------------------------------------------------
!     Begin adjoint source calculation
!     ------------------------------------------------------------------

      if (present(calc_adjoint).and.calc_adjoint) then

        fp(:) = 0 ; fq(:) = 0
        syn_veloc(:) = 0 ; syn_accel(:) = 0

!       Calculate adjoint sources for dt
        
!       Calculate velocity synthetics
!       ----------------------------------------
        do i = 2, nlen-1
         syn_veloc(i) =  (syn_win_save(i+1) - syn_win_save(i-1)) / (2 * dt)
        enddo
        syn_veloc(1) = (syn_win_save(2) - syn_win_save(1)) / dt
        syn_veloc(nlen) = (syn_win_save(nlen) - syn_win_save(nlen-1)) /dt

!       Calculate acceleration synthetics
!       ----------------------------------------
        do i = 2, nlen-1
          syn_accel(i) =  (syn_veloc(i+1) - syn_veloc(i-1)) / (2 * dt)
        enddo
        syn_accel(1) = (syn_veloc(2) - syn_veloc(1)) / dt
        syn_accel(nlen) = (syn_veloc(nlen) - syn_veloc(nlen-1)) /dt

!       DEBUG: output 
        if(DEBUG) then
        open(unit=11, file='DEBUG_xcor_vel_acc.dat')
        do i = 1, nlen
          write(11,'(5(e12.4,1x))') (i-1)*dt, syn_win_save(i), syn_veloc(i), syn_accel(i), obs_win_save(i)
        enddo      
        close(11)
        endif


!       Calculate ajoint source for traveltime
!       ----------------------------------------
        Nnorm = -dt*sum(syn_win_save(1:nlen)*syn_accel(1:nlen)*t_window(1:nlen)*t_window(1:nlen))
        fp(1:nlen)= dtau_xcorr*syn_veloc(1:nlen)*t_window(1:nlen)/Nnorm 
        fp_nomeasure(1:nlen)= syn_veloc(1:nlen)*t_window(1:nlen)/Nnorm 
!        call t_taper(fp,nlen,HANNING,0.05)

!       definition of Dahlen and Baig (2002), Eq. 3,17,18 : dlnA = Aobs/Asyn - 1
!       Calculate adjoint source for amplitude      
!       ------------------------------------------------------
        Mnorm = dt*sum(syn_win_save(1:nlen)*syn_win_save(1:nlen)*t_window(1:nlen)*t_window(1:nlen))
        fq(1:nlen) = -dlnA_xcorr*syn_win_save(1:nlen)*t_window(1:nlen)/Mnorm
        fq_nomeasure(1:nlen) = -syn_win_save(1:nlen)*t_window(1:nlen)/Mnorm
!        call t_taper(fq,nlen,HANNING,0.10)


!       DEBUG: output transfer function
        if(DEBUG) then
        open(unit=11, file='DEBUG_fp_fq.dat')
        do i = 1, nlen
          write(11,'(3(e12.4,1x))') (i-1)*dt, fp(i), fq(i)
        enddo      
        close(11)
        endif

!       DEBUG: output transfer function
        if(DEBUG) then
        open(unit=11, file='DEBUG_fp_fq_nomeasure.dat')
        do i = 1, nlen
          write(11,'(3(e12.4,1x))') (i-1)*dt, fp_nomeasure(i), fq_nomeasure(i)
        enddo      
        close(11)
        endif


      endif

!     ------------------------------------------------------------------
!     End adjoint source calculation
!     ------------------------------------------------------------------



      end subroutine

!     ------------------------------------------------------------------
!     END MAIN SUBROUTINE
!     ------------------------------------------------------------------


