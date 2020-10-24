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
module source

  use global_parameters
  use data_mesh
  use data_source
  use data_time
  use data_proc
  use data_io
  use commun, only: pcheck

  implicit none

  public :: compute_stf, compute_stf_t, compute_src, read_sourceparams
  private
contains

!-----------------------------------------------------------------------------------------
subroutine read_sourceparams

  use utlity, only: to_lower

  real(kind=realkind)   :: srclat
  character(len=256)    :: keyword, keyvalue, line
  character(len=512)    :: errmsg
  integer               :: iinparam_source=500, ioerr


  if (verbose > 1 .and. lpr) write(*,'(A)', advance='no') '    Reading inparam_source...'
  open(unit=iinparam_source, file='inparam_source', status='old', action='read',  iostat=ioerr)
  call pcheck(ioerr /= 0, 'Check input file ''inparam_source''! Is it still there?')

  do
    read(iinparam_source, fmt='(a256)', iostat=ioerr) line
    if (ioerr < 0) exit
    if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

    read(line,*) keyword, keyvalue

    parameter_to_read : select case(trim(keyword))

    case('SOURCE_TYPE')
        read(keyvalue, *) src_type(2)
        src_type(2) = to_lower(src_type(2))
        select case(src_type(2))
        case('mrr', 'explosion', 'mtt_p_mpp', 'mpp_p_mtt', 'vertforce')
            src_type(1) = 'monopole'
        case('mtr', 'mpr', 'mrt', 'mrp', 'thetaforce', 'phiforce')
            src_type(1) = 'dipole'
        case('mtp', 'mpt', 'mtt_m_mpp')
            src_type(1) = 'quadpole'
        case default
            write(errmsg,*) 'Unknown SOURCE_TYPE in inparam_source: ', trim(src_type(2)), '\n', &
                 'Allowed are:  mrr, explosion, mtt_p_mpp, vertforce  (monopole) \n', &
                 '              mtr, mpr, thetaforce, phiforce        (dipole)\n', &
                 '              mtp, mtt_m_mpp                        (quadpole)\n'
            call pcheck(.true., errmsg)

        end select

    case('SOURCE_DEPTH')
        read(keyvalue,*) src_depth
        src_depth = src_depth * 1000 ! The solver works in meters

    case('SOURCE_LAT')
        read(keyvalue,*) srclat
        srccolat = 90.0 - srclat
        srccolat = srccolat * pi / 180.d0

    case('SOURCE_LON')
        read(keyvalue,*) srclon
        srclon = srclon * pi / 180.d0

    case('SOURCE_AMPLITUDE')
        read(keyvalue,*) magnitude

    end select parameter_to_read

  enddo
  close(iinparam_source)

  if (srccolat /= 0.d0 .or. srclon /= 0.d0 ) then
     if (lpr .and. verbose > 0) &
        write(*,'(/,a,a,/,a,a)')&
            procstrg, '  Source not along the axis!', &
            procstrg,'  ...therefore applying rotations to source and receivers.'
     rot_src = .true.

     if (rec_file_type == 'database') then
       write(errmsg,'(a,/,a)') &
          'For the database function, the source should be at the northpole.', &
          'All rotation is done in postprocessing.'
       call pcheck(.true., errmsg)
     endif
  else
     rot_src = .false.
  endif

  if (lpr .and. verbose > 0) then
     print *
     write(*,*)  '  *****************GIVEN SOURCE PARAMETERS*****************'
     write(*,11) '   Magnitude [Nm]:       ', magnitude
     write(*,13) '   Excitation type:      ', src_type(1),src_type(2)
     write(*,11) '   Depth [km]:           ', src_depth/1000.
     write(*,11) '   Colat. [deg]:         ', real(srccolat*180./pi)
     write(*,11) '   Long. [deg]:          ', real(srclon*180./pi)
     write(*,14) '   Source time function: ', stf_type
     write(*,12) '   Dom. period mesh [s]: ', period
     write(*,*)  '  *********************************************************'
     print *
     call flush(6)
  endif

11 format(a28,1pe15.3)
12 format(a28,f15.4)
13 format(a28,2(a12))
14 format(a28,a13)

end subroutine read_sourceparams
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_stf
  use nc_routines, only: nc_dump_stf
  character(len=512)     :: errmsg
  integer                :: i

  allocate(stf(1:niter))
  dt_src_shift = 10000000.

  select case(stf_type)
  case('dirac_0')
    call delta_src ! discrete Dirac's
  case('dirac_1')
    call delta_src
  case('gauss_0')
    call gauss
  case('gauss_1')
    call gauss_d
  case('gauss_2')
    call gauss_dd
  case('errorf')
    call errorf
  case('quheavi')
     !call quasiheavi
     call delta_src ! done inside the delta routine now
  case default
     write(errmsg,*)' source time function non existant:', stf_type
     call pcheck(.true., errmsg)
  end select

  write(errmsg, '("ERROR: source time shift not defined! ", F15.4)') dt_src_shift
  call pcheck(dt_src_shift > 1000000., errmsg)

  write(errmsg, *) & !'(a,/,a,2(f7.4))') &
                'Problem with discrete source shift: not a multiplicative of deltat...', &
                'source shift, deltat: ', dt_src_shift, deltat !, dt_src_shift/deltat
  call pcheck(abs(nint(dt_src_shift / deltat) - dt_src_shift / deltat) > 0.01 * deltat, errmsg)

  ! time shift in the Fourier domain (used in post processing/kerner... eventually)
  ! timeshift_fourier(0:nomega) = exp(cmplx(0.,1.) *omega(0:nomega)*dt_src_shift)

   if (use_netcdf .and. (mynum == 0)) call nc_dump_stf(stf(1:niter))

   if (lpr) then
      if (.not. (use_netcdf)) then
          open(299, file=datapath(1:lfdata)//'/stf.dat', status='replace', action='write')
          open(298, file=datapath(1:lfdata)//'/stf_seis.dat', status='replace', action='write')
          open(297, file=datapath(1:lfdata)//'/stf_strain.dat', status='replace', action='write')
          do i=1, niter
             write(299,*) real(i) * real(deltat), real(stf(i))
             if ( mod(i,seis_it) == 0) write(298,*) real(i) * real(deltat), real(stf(i))
             if ( mod(i,strain_it) == 0) write(297,*) real(i) * real(deltat), real(stf(i))
          enddo
          close(299)
          close(298)
          close(297)
      endif
   endif

end subroutine compute_stf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> These *_t routines are needed by the symplectic time integration schemes.
!! Eventually there should be only one type, point- or array-wise.
subroutine compute_stf_t(nstf_t,t,stf_t)

  integer, intent(in)        :: nstf_t
  real(kind=dp), intent(in)  :: t(1:nstf_t)
  real(kind=dp), intent(out) :: stf_t(1:nstf_t)


  select case(stf_type)
  case('dirac_0')
    call delta_src_t(nstf_t, t, stf_t)
  case('gauss_0')
    call gauss_t(nstf_t, t, stf_t)
  case('gauss_1')
    call gauss_d_t(nstf_t, t, stf_t)
  case('gauss_2')
    call gauss_dd_t(nstf_t, t, stf_t)
  case('errorf')
    call errorf_t(nstf_t, t, stf_t)
  case('quheavi')
    call quasiheavi_t(nstf_t, stf_t)
  case default
     write(*,*)' source time function non existant:', stf_type
     stop
  end select

end subroutine compute_stf_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_src

  use data_mesh
  use utlity
  use commun, only: broadcast_int,broadcast_dble

  integer                          :: iel_src2,ipol_src2,jpol_src2
  real(kind=realkind), allocatable :: source_term(:,:,:,:)
  real(kind=realkind), allocatable :: point_source(:,:,:)
  integer                          :: ipol,jpol,ielem,k
  real(kind=dp)                    :: s,z,r,theta
  character(len=256)               :: errmsg

  allocate(source_term(0:npol,0:npol,1:nel_solid,1:3))
  source_term = 0.

  zsrc = router - src_depth

  if (lpr) write(*,'(a,/,a,/,a,/,a)') &
            '  *****************************************', &
            '     Welcome to the source term calculator ', &
            '  *****************************************', &
            '  locating the source....'

  if (verbose > 1) write(69,'(/,a)')'    L O O K I N G   F O R   T H E   S O U R C E '

  call find_srcloc(iel_src2, ipol_src2, jpol_src2)

  ! @TODO: should add a test whether it is the first theta slice

  poletype: select case(src_type(1))

  ! MONOPOLE
  case ('monopole') poletype
     if (lpr) write(*,*) '  computing MONOPOLE Source with...'

     select case (src_type(2))
     case('vertforce')
        if (lpr) write(*,*) '  ...vertical single force'
        allocate(point_source(0:npol,0:npol,1:nel_solid))
        call define_bodyforce(point_source,iel_src2,ipol_src2,jpol_src2)
        source_term(:,:,:,3) = point_source / (two * pi)
        deallocate(point_source)

     case default
        if (lpr) write(*,*) '  ...moment tensor elements for ', src_type(2)
        call define_moment_tensor(iel_src2, ipol_src2, jpol_src2, source_term)
        source_term = source_term / (two * pi)

     end select

  ! DIPOLE
  case ('dipole') poletype
     if (lpr) write(*,*) '  computing DIPOLE Source with...'

     select case(src_type(2))
     case ('thetaforce', 'phiforce')
        if (lpr) write(*,*) '  ...horizontal single ', src_type(2)
        allocate(point_source(0:npol,0:npol,1:nel_solid))
        call define_bodyforce(point_source, iel_src2, ipol_src2, jpol_src2)
        source_term(:,:,:,1) = point_source / pi
        deallocate(point_source)

     case default
        if (lpr) write(*,*) '  ...moment tensor elements for ', src_type(2)
        call define_moment_tensor(iel_src2, ipol_src2, jpol_src2, source_term)
        source_term = source_term / pi

     end select

  ! QUADRUPOLE
  case ('quadpole') poletype
     if (lpr) write(*,*)'  computing QUADRUPOLE Source with...'
     if (lpr) write(*,*)'  ...moment tensor elements for ',src_type(2)
     call define_moment_tensor(iel_src2, ipol_src2, jpol_src2, source_term)
     source_term = source_term / pi

  case default
     write(errmsg,*) 'we only support monopoles, dipoles, quadrupoles, and not ',src_type(1)
     ! Will stop the program
     call pcheck(.true., errmsg)

  end select poletype

  ! if I don't have the source
  if (.not. have_src) then
     if (verbose > 1) write(69,'(/,a,/)') &
            "******  I  D O N ' T   H A V E   T H E   S O U R C E *******"
  endif

  ! write all elements containing non-zero source term components to file
  if (diagfiles .and. have_src) then
     open(619, file=infopath(1:lfinfo)//'/src_term.dat'//appmynum)
     open(621, file=infopath(1:lfinfo)//'/src_term_norm1.dat'//appmynum)
     open(622, file=infopath(1:lfinfo)//'/src_term_norm2.dat'//appmynum)
     open(623, file=infopath(1:lfinfo)//'/src_term_norm3.dat'//appmynum)
     open(6200, file=infopath(1:lfinfo)//'/src_term_allcomp.dat'//appmynum)
     do ielem=1, nel_solid
        if (maxval(abs(source_term(:,:,ielem,:))) /= zero) then
           do ipol=0, npol
              do jpol=0, npol
                call compute_coordinates(s, z, r, theta, ielsolid(ielem), ipol, jpol)
                write(619,12) ielem, ipol, jpol, source_term(ipol,jpol,ielem,1), &
                              source_term(ipol,jpol,ielem,2), source_term(ipol,jpol,ielem,3)
                write(621,13) s/router, z/router, source_term(ipol,jpol,ielem,1)/&
                              maxval(abs(source_term(:,:,:,1)))*7.
                write(622,13) s/router, z/router, source_term(ipol,jpol,ielem,2)/&
                              maxval(abs(source_term(:,:,:,2)))*7.
                write(623,13) s/router, z/router, source_term(ipol,jpol,ielem,3)/&
                              maxval(abs(source_term(:,:,:,3)))*7.
              enddo
           enddo
        endif
     enddo
12  format(i9,2(i2),3(1pe12.4))
13  format(3(1pe12.4))

    close(619)
    close(621)
    close(622)
    close(623)
  endif !have_src

  ! construct source term array that only lives on nonzero elements (max. 8)
  allocate(source_term_el(0:npol,0:npol,8,3))
  source_term_el = zero
  k = 0
  do ielem = 1, nel_solid
     if ( maxval(abs(source_term(:,:,ielem,:))) > zero) then
        k = k + 1
        ielsrc(k) = ielem
        source_term_el(:,:,k,:) = source_term(:,:,ielem,:)
     endif
  enddo
  nelsrc = k

  if (verbose > 1) write(69,*) 'nelsrc =', nelsrc

  deallocate(source_term)

  if (nelsrc > 8) then
     write(*,'(a,a,/,a,i4,a,/,a,a)') &
             procstrg, 'PROBLEM with source term element count!', &
             procstrg, nelsrc, ' elements with nonzero entries found but 8 is max', &
             procstrg, '(hitting the edge of an element at the bottom of a doubling layer).'
     stop
  endif

  if (have_src) then
     write(*,'(a,i2,/,a,/,a,/,a)') &
           '  number of elements with non-zero source term:', nelsrc, &
           '  *********************************', &
           '     End of source term calculator', &
           '  *********************************'
  endif
  if (verbose > 1) then
     write(69,*)'  *********************************'
     write(69,*)'     End of source term calculator'
     write(69,*)'  *********************************'
  endif

end subroutine compute_src
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine find_srcloc(iel_src2, ipol_src2, jpol_src2)

  use data_mesh
  use utlity
  use commun, only: pmin, psum_int

  integer, intent(out) :: iel_src2, ipol_src2, jpol_src2
  real(kind=dp)        :: s, z, r, theta, mydzsrc, zsrcout, dzsrc
  integer              :: ielem, ipol, jpol, count_src_procs
  character(len=256)   :: errmsg

  ! find depth that is closest to desired value zsrc

  dzsrc = 10.d0 * router

  write(errmsg, "('Source depth: ', F10.2, ' is bigger than model radius: ', F10.2)") zsrc, router
  call pcheck(zsrc > router, errmsg)

  ! Only allow sources in the solid region, fixated to northern axis.
  do ielem = 1, nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           call compute_coordinates(s, z, r, theta, ielsolid(ielem), ipol, jpol)
           if (abs(s) > smallval_dble .or. z < smallval_dble) cycle
           if (dabs(z-zsrc) < dzsrc) then
              zsrcout = z
              dzsrc = dabs(zsrcout - zsrc)
              iel_src = ielem
              ipol_src = ipol
              jpol_src = jpol
              iel_src2 = ielem
              ipol_src2 = ipol
              jpol_src2 = jpol
           else if (dabs(z-zsrc) == dzsrc) then
              if (verbose > 1) write(69,15) ielem,ipol,jpol,z/1000.
              iel_src2 = ielem
              ipol_src2 = ipol
              jpol_src2 = jpol
           endif
        enddo
     enddo
  enddo
15 format('  found a second point with same distance:', i6, i3, i3, 1pe13.3)


  ! Make sure only closest processor has source
  mydzsrc = dzsrc
  have_src = .true.

  if (nproc > 1) then
     mydzsrc = pmin(mydzsrc)
     if (mydzsrc < dzsrc) have_src = .false.

     ! Check how many/which processors have the source
     count_src_procs = 0
     if (have_src) count_src_procs = 1
     count_src_procs = psum_int(count_src_procs)

     if (count_src_procs > 2) then
        if (lpr) then
           write(*,*)
           write(*,*) 'PROBLEM with source & processors!'
           write(*,*) 'More than two processors have the source:'
        endif
        if (have_src) write(*,*) procstrg, 'has it.'
        stop
     else if (count_src_procs == 0) then
        if (lpr) then
           write(*,*)
           write(*,*) 'PROBLEM with source & processors!'
           write(*,*) 'No processor has the source.'
        endif
        stop
     endif
  endif !nproc > 1

  if (have_src) then
     if (ipol_src /= 0) then
        write(*,'(a,/,a,i7,i2,i2)') &
              'ERROR: Source should be on axis, i.e. ipol_src=0, but:', &
              'Source location: ielem, ipol, jpol: ', &
              ielsolid(iel_src), ipol_src, jpol_src
        stop
     endif

     if (thetacoord(ipol_src, jpol_src, ielsolid(iel_src)) /= zero) then
        write(*,'(a,/,i7,2i2,/,a,3e10.2)') &
                'ERROR: Source should be on the axis, hence theta = 0, but:', &
                'Source indices ielem, ipol, jpol:', &
                ielsolid(iel_src), ipol_src, jpol_src, &
                's,r,theta', scoord(ipol_src,jpol_src,ielsolid(iel_src)), &
                rcoord(ipol_src,jpol_src,ielsolid(iel_src)), &
                thetacoord(ipol_src,jpol_src,ielsolid(iel_src))
        stop
     endif

     write(*,*) '  ',procstrg,' found it:'
     write(*,*) '    depth asked for [km]:', (router - zsrc) / 1000.d0
     write(*,*) '    depth offered   [km]:', (router - zsrcout) / 1000.d0
     write(*,*) '    difference      [km]:', dzsrc / 1000.d0
     if (verbose > 1) then
         write(*,*) '    source element      : ', iel_src, ielsolid(iel_src)
     endif

     if (verbose > 1) then
        write(69,*) '  ',procstrg,' found it:'
        write(69,*) '    depth asked for [km]:',(router-zsrc)/1000.d0
        write(69,*) '    depth offered   [km]:',(router-zsrcout)/1000.d0
        write(69,*) '    difference      [km]:',dzsrc/1000.d0
        write(69,*) '    source element and jpol index:', iel_src,jpol_src
     endif

     if (iel_src2 /= iel_src) then
         call compute_coordinates(s, z, r, theta, ielsolid(iel_src2), ipol_src2, jpol_src2)
         write(*,*) '    SECOND source element and jpol index:', &
                            iel_src2, jpol_src2
         write(*,*) '       s, z: ', s, z

         if (verbose > 1) then
            write(69,*) '    SECOND source element and jpol index:', &
                               iel_src2,jpol_src2
            write(69,*) '      s, z: ', s, z
            write(69,*) '    depth offered   [km]:',(router-z)/1000.d0
         endif
     endif

     zsrc = zsrcout
  endif

end subroutine find_srcloc
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gauss

  integer               :: i
  real(kind=realkind)   :: t

  do i=1, niter
     t = dble(i) * deltat
     stf(i) = dexp(-( (decay / t_0 * (t - shift_fact))**2) )
  enddo
  dt_src_shift = shift_fact
  stf = stf * magnitude * decay / t_0 / dsqrt(pi)

end subroutine gauss
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gauss_d

  integer               :: i
  real(kind=realkind)   :: t

  do i=1, niter
     t = dble(i) * deltat
     stf(i) = -two * (decay / t_0)**2 * (t - shift_fact) &
          * dexp(-( (decay/t_0*(t-shift_fact))**2) )
  enddo
  dt_src_shift = shift_fact
  ! max/min at t=t_0*(shift_fact +- sqrt(0.5)/decay)
  ! and corresponding max(stf)= +- decay/t_0*sqrt(two)*exp(-half)

  stf = stf / ( decay / t_0 * sqrt(two) * exp(-half) ) * magnitude

end subroutine gauss_d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gauss_dd

  integer               :: i
  real(kind=realkind)   :: t

  do i=1, niter
     t = dble(i) * deltat
     stf(i) = (decay / t_0)**2 * ( two * (decay / t_0)**2 * (t - shift_fact)**2 - 1) &
               * dexp(-( (decay / t_0 * (t - shift_fact))**2) )
  enddo
  dt_src_shift = shift_fact

  ! max/min at t=t_0*(shift_fact +- sqrt(1.5)/decay)
  ! and corresponding max(stf)= +- two*((decay/t_0)**2*exp(-three/two)
  ! and at t=shift_fact*t_0

  stf = stf / ( two * ((decay / t_0)**2) * exp(-three / two) ) * magnitude

end subroutine gauss_dd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine errorf

  integer               :: i
  real(kind=realkind)   :: t

  do i=1, niter
     t = dble(i) * deltat
     stf(i) = erf(decay / t_0 * (t - shift_fact))*0.5 + 0.5
  enddo

  dt_src_shift = shift_fact ! Taken from Gauss

  stf = stf * magnitude

end subroutine errorf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Calculates the error function, coefficients taken from Numerical Recipes
!! @TODO: ERF is an intrinsic in Fortran2008. Compiler support is unclear still,
!!        so we will keep that in here for a while.
function erf(x)

  real(kind=dp), intent(in)        :: x
  real(kind=dp)                    :: erfcc, erf
  real(kind=dp)                    :: polyval
  real(kind=dp)                    :: t,z
  integer                          :: icoeff
  integer, parameter               :: ncoeff = 10
  real(kind=dp), dimension(ncoeff) :: coeffs =                         &
      [-1.26551223,  1.00002368, 0.37409196,  0.09678418, -0.18628806, &
        0.27886807, -1.13520398, 1.48851587, -0.82215223,  0.17087277]

  z = abs(x)
  t = 1.0 / (1.0 + 0.5 * z)

  polyval = coeffs(ncoeff)
  do icoeff = ncoeff-1, 1, -1
      polyval = t * polyval + coeffs(icoeff)
  enddo

  erfcc = t * exp(-z * z + polyval)

  if (x < 0.0) erfcc = 2.0 - erfcc
  erf = 1 - erfcc

end function erf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> approximate discrete dirac
subroutine delta_src
!@TODO disrecete hidden treasures do not work with symplectic schemes
  integer                   :: i,j
  real(kind=dp)             :: a, integral
  character(len=6)          :: dirac_approx(6)
  real(kind=dp),allocatable :: signal(:), timetmp(:), int_stf(:)

  if (lpr) write(*,*)'Discrete Dirac choice: ',discrete_choice
  allocate(signal(1:niter))
  allocate(timetmp(1:niter))
  allocate(int_stf(1:niter))

  stf = 0
  a = discrete_dirac_halfwidth

  if (lpr) write(*,*) 'Half-width of discrete Dirac [s]: ', a

  dirac_approx = ['cauchy','caulor','sincfc','gaussi','triang','1dirac']
  !@TODO: there is no userinterface for this atm, so only gaussi and 1dirac are
  !       actually used

  do j=1, 6
     signal = 0
     if (lpr .and. trim(discrete_choice) == trim(dirac_approx(j))) &
              write(*,*)' Approximation type:',trim(dirac_approx(j))
     if (lpr .and. diagfiles) then
         open(unit=60,file=infopath(1:lfinfo)//'/discrete_dirac_'//trim(dirac_approx(j))//'.dat')
     endif

     do i=1,niter
        t=dble(i)*deltat
        timetmp(i) = t
        if (dirac_approx(j) == 'cauchy') then ! Cauchy phi function
           signal(i) = 1./a * exp(-abs((t - shift_fact_discrete_dirac) / a))
           dt_src_shift = shift_fact_discrete_dirac

        else if (dirac_approx(j) == 'caulor') then ! Cauchy Lorentz distribution
           signal(i) = 1. / pi * a / (a**2 + (t - shift_fact_discrete_dirac)**2)
           dt_src_shift = shift_fact_discrete_dirac

        else if (dirac_approx(j) == 'sincfc') then ! sinc function
           !@TODO: is this because of division by zero? what about putting signal
           !       to 1? (MvD)
           if (t == shift_fact_discrete_dirac) t = 0.00001 + shift_fact_discrete_dirac

           signal(i) = 1. / (a * pi) * sin((-shift_fact_discrete_dirac + t) / a) &
                            / ((-shift_fact_discrete_dirac + t) / a)
           dt_src_shift = shift_fact_discrete_dirac

        else if (dirac_approx(j) == 'gaussi') then ! Gaussian
           signal(i) = 1. / (a * sqrt(pi)) * exp(-((t - shift_fact_discrete_dirac) / a)**2)
           dt_src_shift = shift_fact_discrete_dirac

        else if (dirac_approx(j) == 'triang') then ! triangular
           if (abs(t-shift_fact_discrete_dirac) <= a / 2.) then
              signal(i) = 2. / a - 4. / a**2 * abs(t - shift_fact_discrete_dirac)
           else
              signal(i) = 0.
           endif
           dt_src_shift = shift_fact_discrete_dirac

        else if (dirac_approx(j) == '1dirac') then ! old Dirac, 1 non-zero point
           if (i == int(shift_fact_discrete_dirac / deltat)) then
              signal(i) = 1.
              dt_src_shift = real(i) * deltat
           endif
        else
           write(*,*)'do not know discrete Dirac ', trim(dirac_approx(j))
           stop
        endif

        if (lpr .and. diagfiles) write(60,*) t, magnitude * signal(i)
     enddo

     if (lpr .and. diagfiles) close(60)

     if (trim(discrete_choice) == trim(dirac_approx(j)) ) then
        if (lpr) write(*,*)'  dirac type and max amp before:',trim(dirac_approx(j)),maxval(signal)
        integral = sum(signal)*deltat
        if (lpr) write(*,*)'  OLD Integral of discrete dirac:',integral
        signal = signal/integral

        if (lpr) write(*,*)'  dirac type and max amp after:',trim(dirac_approx(j)),maxval(signal)
        integral = sum(signal)*deltat
        if (lpr) write(*,*)'  NEW Integral of discrete dirac:',integral
        if (lpr) write(*,*)"  Shift factor [s],#dt's:",dt_src_shift,dt_src_shift/deltat
        if (lpr) write(*,*)

        stf(1:niter) = signal(1:niter)
     endif
  enddo

  stf = stf * magnitude

  ! Integrate STF over time for Heaviside function
  int_stf = 0
  signal = 0
  do i=1, niter
     if (i > 1) signal(i) = int_stf(i-1)
     int_stf(i) = signal(i) + stf(i) * deltat
  enddo

  if (lpr .and. diagfiles) then
     open(unit=61,file=infopath(1:lfinfo)//'/discrete_chosen_dirac_'//trim(discrete_choice)//'.dat')
     open(unit=62,file=infopath(1:lfinfo)//'/discrete_chosen_heavi_'//trim(discrete_choice)//'.dat')
     do i=1, niter
        write(61,*) timetmp(i), stf(i)
        write(62,*) timetmp(i), int_stf(i)
     enddo
     close(61)
     close(62)
  endif

  ! Quasi-Heaviside
  if (trim(stf_type) == 'quheavi') stf = int_stf

  deallocate(timetmp,signal,int_stf)

end subroutine delta_src
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gauss_t(nstf_t,t,stf_t)

  integer, intent(in)           :: nstf_t
  real(kind=dp), intent(in)     :: t(nstf_t)
  real(kind=dp), intent(out)    :: stf_t(nstf_t)
  integer                       :: i

  do i=1,nstf_t
     stf_t(i) = dexp(-( (decay/t_0*(t(i)-shift_fact))**2) )
  enddo
  dt_src_shift = shift_fact
  stf_t = stf_t * magnitude * decay / ( t_0 * sqrt(pi) )

end subroutine gauss_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gauss_d_t(nstf_t,t,stf_t)

  integer, intent(in)              :: nstf_t
  real(kind=dp), intent(in)        :: t(nstf_t)
  real(kind=dp), intent(out)       :: stf_t(nstf_t)
  integer                          :: i

  do i=1,nstf_t
     stf_t(i) = -two * (decay / t_0)**2 * (t(i) - shift_fact) &
                    * dexp(- (decay / t_0 * (t(i) - shift_fact))**2 )
  enddo
  dt_src_shift = shift_fact
  stf_t = stf_t / (decay / t_0 * dsqrt(two) * dexp(-half)) * magnitude

end subroutine gauss_d_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gauss_dd_t(nstf_t, t, stf_t)

  integer, intent(in)              :: nstf_t
  real(kind=dp), intent(in)        :: t(nstf_t)
  real(kind=dp), intent(out)       :: stf_t(nstf_t)
  integer                          :: i

  !@TODO: seems like there are a few factors of decay / t_0 to much (MvD)
  do i=1, nstf_t
     stf_t(i) = (decay / t_0)**2 &
                    * (two * (decay / t_0)**2 * (t(i) - shift_fact)**2 - 1) &
                    * dexp(-(decay / t_0 * (t(i) - shift_fact))**2)
  enddo
  dt_src_shift = shift_fact
  stf_t = stf_t / (two * ((decay / t_0)**2) * exp(-three / two)) * magnitude

end subroutine gauss_dd_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine errorf_t(nstf_t, t, stf_t)

  integer, intent(in)              :: nstf_t
  real(kind=dp), intent(in)        :: t(nstf_t)
  real(kind=dp), intent(out)       :: stf_t(nstf_t)
  integer                          :: i

  do i=1,nstf_t
     stf_t(i) = erf(decay / t_0 * (t(i) - shift_fact)) * 0.5 + 0.5
  enddo
  dt_src_shift = shift_fact
  stf_t = stf_t * magnitude

end subroutine errorf_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine delta_src_t(nstf_t, t,stf_t)

  integer, intent(in)             :: nstf_t
  real(kind=dp), intent(in)       :: t(nstf_t)
  real(kind=dp), intent(out)      :: stf_t(nstf_t)

  stf_t(1:nstf_t) = zero

  if (t(1) > (shift_fact - deltat) .and. t(1) <= shift_fact) &
     stf_t = (t - t(1)) / deltat * magnitude / deltat

  if (t(1) >= shift_fact .and. t(1) < (shift_fact + deltat)) &
     stf_t = (1. - (t - t(1)) / deltat) * magnitude / deltat

end subroutine delta_src_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine quasiheavi_t(nstf_t, stf_t)

  integer, intent(in)           :: nstf_t
  real(kind=dp), intent(out)    :: stf_t(nstf_t)

  stf_t = 0
  stf_t(seis_it:nstf_t) = magnitude
  dt_src_shift = seis_it

end subroutine quasiheavi_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_bodyforce(f, iel_src2, ipol_src2, jpol_src2)

  use data_mesh
  use utlity
  use commun, only: pdistsum_solid_1D

  real(kind=realkind), intent(out) :: f(0:npol,0:npol,nel_solid)
  integer, intent(in)              :: iel_src2, ipol_src2, jpol_src2
  integer                          :: liel_src, lipol_src, ljpol_src, ipol, jpol, i
  character(len=16)                :: fmt1
  integer                          :: nsrcelem

  nsrcelem = 1
  if (iel_src2 /= iel_src) nsrcelem = 2
  f(:,:,:) = zero

  if (have_src) then
     f(ipol_src, jpol_src, iel_src) = one
     f(ipol_src2, jpol_src2, iel_src2) = one
  endif

  ! assembly
  call pdistsum_solid_1D(f)

  if (have_src .and. verbose > 1) then

     ! write out the source element
     fmt1 = "(K(1pe12.3))"
     write(fmt1(2:2),'(i1.1)') npol+1

     write(69,'(/,a)') '  *^*^*^*^**^*^*^* The single-force source term *^*^*^*^*^^*^'

     liel_src = iel_src
     lipol_src = ipol_src
     ljpol_src = jpol_src

     do i=1, nsrcelem
        if (i == 2) then
           liel_src = iel_src2
           lipol_src = ipol_src2
           ljpol_src = jpol_src2
        endif
        write(69,*) 'iel,jpol,r:', liel_src, ljpol_src, &
             rcoord(lipol_src, ljpol_src, ielsolid(liel_src)) / 1.d3
        write(69,*) 'North| s-dir -->'
        do jpol=npol, 0, -1
           write(69,fmt1) (f(ipol, jpol, liel_src), ipol=0,npol)
        enddo
        write(69,*)
        write(69,*)
     enddo
  endif ! have_src

end subroutine define_bodyforce
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Defines the moment tensor elements for the given source type in all
!! elements having non-zero source contributions,
!! using pointwise derivatives of arbitrary scalar functions.
subroutine define_moment_tensor(iel_src2, ipol_src2, jpol_src2, source_term)

  use data_mesh

  use apply_masks
  use utlity
  use pointwise_derivatives
  use commun, only: pdistsum_solid, psum_int

  integer, intent(in)              :: iel_src2, ipol_src2, jpol_src2
  real(kind=realkind), intent(out) :: source_term(0:npol,0:npol,nel_solid,3)
  integer                          :: liel_src, lipol_src, ljpol_src

  real(kind=realkind), allocatable :: ws(:,:,:), dsws(:,:,:)
  real(kind=realkind), allocatable :: ws_over_s(:,:,:), dzwz(:,:,:)
  real(kind=realkind), allocatable :: ds(:), dz(:)

  integer                          :: ielem, ipol, jpol, i, nsrcelem, nsrcelem_glob
  character(len=16)                :: fmt1

  liel_src = iel_src
  lipol_src = ipol_src
  ljpol_src = jpol_src

  nsrcelem = 0
  if (have_src) then
     nsrcelem = 1
     if (iel_src2 /= iel_src) then
        nsrcelem = 2
     else
        nsrcelem = 1
     endif
  endif

  allocate(ws(0:npol,0:npol,1:nsrcelem))
  allocate(dsws(0:npol,0:npol,1:nsrcelem))
  allocate(ws_over_s(0:npol,0:npol,1:nsrcelem))
  allocate(dzwz(0:npol,0:npol,1:nsrcelem))
  allocate(ds(0:npol),dz(0:npol))

  dzwz(:,:,:) = zero
  dsws(:,:,:) = zero
  ws_over_s(:,:,:) = zero
  source_term(:,:,:,:) = zero

  ! global number of source elements (in case source is on processor boundary)
  nsrcelem_glob = psum_int(nsrcelem)

  if (verbose > 1 ) write(69,*) 'nsrcelem_glob =', nsrcelem_glob

  if (have_src) then
     ! physical source location can only be in 2 elements
     do i=1, nsrcelem
        if (i == 2) then
           liel_src = iel_src2
           lipol_src = ipol_src2
           ljpol_src = jpol_src2
        endif

        do ipol = 0,npol
           do jpol = 0,npol

              ws(:,:,:) = zero
              ws(ipol,jpol,i) = one
              call dsdf_elem_solid(dsws(:,:,i), ws(:,:,i), liel_src)
              call dzdf_elem_solid(dzwz(:,:,i), ws(:,:,i), liel_src)

              poletype:  select case (src_type(1))

              ! monopole
              case ('monopole') poletype
                 select case (src_type(2))

                 case ('explosion')
                    if (ipol == 0 .and. jpol == 0 .and. lpr) &
                         write(*,*)'  ',procstrg, &
                        'computing source s- and z-components for explosion'

                    source_term(ipol,jpol,liel_src,1) = &
                         two * dsws(lipol_src,ljpol_src,i)
                    source_term(ipol,jpol,liel_src,3) = &
                         dzwz(lipol_src,ljpol_src,i)

                 case ('mtt_p_mpp' )
                    if (ipol == 0 .and. jpol == 0 .and. lpr)  &
                         write(*,*)'  ',procstrg, &
                         'computing source s-component for Mxx+Myy'
                    source_term(ipol,jpol,liel_src,1) = &
                          dsws(lipol_src,ljpol_src,i)

                 case ('mrr')
                    if (ipol == 0 .and. jpol == 0 .and. lpr)  &
                         write(*,*)'  ',procstrg, &
                         'computing source field z-component for Mzz'
                    source_term(ipol,jpol,liel_src,3) = &
                         dzwz(lipol_src,ljpol_src,i)

                 case default
                    write(*,'(a,a,/,a,a,a)') &
                         procstrg, 'PROBLEM: Didn"t compute any source: ', &
                         procstrg, 'Monopole source doesn"t exist for ', src_type(2)
                    stop
                 end select

              ! dipole
              case ('dipole') poletype
                 select case(src_type(2))

                 case ('mtr','mpr')
                    if (ipol == 0 .and. jpol == 0 .and. lpr)  &
                         write(*,*) '  computing source + and z-components for Mtr'
                    source_term(ipol, jpol, liel_src, 1) =  &
                         dzwz(lipol_src, ljpol_src, i)
                    source_term(ipol, jpol, liel_src, 3) =  &
                         dsws(lipol_src, ljpol_src, i)

                 case default
                    write(*,'(a,a,/,a,a,a)') &
                         procstrg, 'PROBLEM: Didn"t compute any source!', &
                         procstrg, 'Dipole source doesn"t exist for ', src_type(2)
                    stop
                 end select

              ! quadrupole
              case ('quadpole') poletype
                 select case (src_type(2))

                 case ('mtp','mtt_m_mpp')
                    if (ipol == 0 .and. jpol == 0 .and. lpr)  &
                         write(*,*) '  computing source s- and phi-components for Mtp'
                    source_term(ipol,jpol,liel_src,1) = &
                         dsws(lipol_src,ljpol_src,i)
                    source_term(ipol,jpol,liel_src,2) = &
                         ws_over_s(lipol_src,ljpol_src,i)
                 case default
                    write(*,'(a,a,/,a,a,a)') &
                         procstrg, "PROBLEM: Didn't compute any source!", &
                         procstrg, "Quadrupole doesn't exist for", src_type(2)
                    stop
                 end select
              end select poletype
           enddo !jpol
        enddo ! ipol
     enddo ! multiple source elements

     ! If spread over multiple elements (i.e., if point source coincides
     ! with element edge/corner), need to divide by global source element number
     source_term = source_term / real(nsrcelem_glob)

     if (verbose > 1) write(69,*) 'source term minmax:', &
                                  minval(source_term), maxval(source_term)
  endif ! have_src

  ! assembly
  if (verbose > 1) write(69,*) '  ', procstrg, 'assembling the source term....'
  call pdistsum_solid(source_term)

  ! cut out round-off errors
  do ielem=1, nel_solid
     do ipol=0, npol
        do jpol=0, npol
           if (abs(source_term(ipol,jpol,ielem,1)) < smallval) &
                source_term(ipol,jpol,ielem,1) = zero
           if (abs(source_term(ipol,jpol,ielem,2)) < smallval) &
                source_term(ipol,jpol,ielem,2) = zero
           if (abs(source_term(ipol,jpol,ielem,3)) < smallval) &
                source_term(ipol,jpol,ielem,3) = zero
        enddo
     enddo
  enddo

  if (have_src) then
     if (maxval(abs(source_term)) == zero) then
        write(*,'(a,a,/,a,a)') procstrg, 'PROBLEM: No source generated!', &
                               procstrg, 'Bye from define_mono_moment'
        stop
     endif
  endif

  ! mask source
  select case (src_type(1))
  case ('monopole')
     call apply_axis_mask_onecomp(source_term, nel_solid, ax_el_solid, &
                                  naxel_solid)
  case ('dipole')
     call apply_axis_mask_twocomp(source_term, nel_solid, ax_el_solid, &
                                  naxel_solid)
  case ('quadpole')
     call apply_axis_mask_threecomp(source_term, nel_solid, ax_el_solid, &
                                    naxel_solid)
  end select

  if (lpr .and. verbose > 0) &
     write(*,*)'  ...masked the source'

  if (have_src) then
     ! write out the source element only
     fmt1 = "(K(1pe12.3))"
     write(fmt1(2:2),'(i1.1)') npol + 1

     if (verbose > 1) then

        write(69,'(/,a)') '  *^*^*^*^*^*^* The moment-tensor source term *^*^*^*^*^**^*^'

        liel_src = iel_src
        lipol_src = ipol_src
        ljpol_src = jpol_src

        do i =1, nsrcelem
           if (i == 2) then
              liel_src = iel_src2
              lipol_src = ipol_src2
              ljpol_src = jpol_src2
           endif

           write(69,*) 'iel,jpol,r:', liel_src, ljpol_src, &
                rcoord(lipol_src,ljpol_src,ielsolid(liel_src)) / 1.d3
           write(69,*)'North| s-dir -->'
           if (src_type(1) == 'dipole') then
              write(69,*) '  ', src_type(2), '+ component'
           else
              write(69,*) '  ', src_type(2), 's component'
           endif
           do jpol=npol, 0, -1
              write(69,fmt1)(source_term(ipol,jpol,liel_src,1), ipol=0, npol)
           enddo
           write(69,*)

           if (src_type(1) == 'dipole') then
              write(69,*)src_type(2), '- component'
           else
              write(69,*)src_type(2), 'phi component'
           endif
           do jpol=npol, 0, -1
              write(69,fmt1)(source_term(ipol,jpol,liel_src,2), ipol=0,npol)
           enddo
           write(69,*)

           write(69,*)src_type(2),'z component'
           do jpol=npol, 0, -1
              write(69,fmt1)(source_term(ipol,jpol,liel_src,3), ipol=0,npol)
           enddo
           write(69,*)
           write(69,*)
        enddo
     endif
  endif ! have_src

end subroutine define_moment_tensor
!-----------------------------------------------------------------------------------------

end module source
!=========================================================================================
