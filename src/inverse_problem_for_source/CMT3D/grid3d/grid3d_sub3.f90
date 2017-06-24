module grid3d_sub3

  use grid3d_constants
  implicit none

contains

  !==========================================================================

  subroutine xcorr_calc(dd,ss,npts,i1,i2,ishift,cc_max,ishift_max)

    ! the core of Alessia's version

    real, dimension(*), intent(in) :: dd,ss
    integer, intent(in) :: npts, i1, i2
    integer, intent(out) :: ishift
    real,intent(out) :: cc_max
    integer,intent(in) :: ishift_max

    ! local variables
    real,dimension(npts) :: d,s
    integer :: nlen
    integer :: i_left, i_right, i, j, i11,i22
    real :: cc

    ! initialise shift
    d=0.; s=0.; d(i1:i2)=dd(i1:i2); s(i1:i2)=ss(i1:i2)
    ishift = 0; cc_max = SMALL

    ! length of window (number of points, including ends)
    nlen = min(i2 - i1 + 1,ishift_max * 2)

    i_left = -1*int(nlen/2.0)
    i_right = int(nlen/2.0)

    ! i -> shift (to be applied to DATA (d) in cc search)
    do i = i_left, i_right
       i11 = max(1,i1+i)
       i22 = min(npts,i2+i)
       cc=sum(s(i11-i:i22-i)*d(i11:i22))
       if (cc > cc_max) then
          cc_max=cc
          ishift=i
       endif
    enddo

  end subroutine xcorr_calc

  !=================================================

  subroutine  write_new_cmtsolution(cmt_file,new_cmt_file,mijn)

    character(len=*),intent(in) :: cmt_file,new_cmt_file
    real,intent(inout) :: mijn(NM)

    integer, parameter :: NSIG_DIGIT = 6
    character(len=150) :: pde_time,event_name, &
         str_tshift,str_hdur,str_lat,str_lon,str_dep
    real :: exponent
    integer :: i, nn, exp_largest, iflag
    real :: epsilon, s1, d1, r1, s2, d2, r2, mw, m0, m00

    open(IOCMT,file = cmt_file, status = 'old')
    read(IOCMT,'(a)') pde_time
    read(IOCMT,'(a)') event_name
    read(IOCMT,'(a)') str_tshift
    read(IOCMT,'(a)') str_hdur
    read(IOCMT,'(a)') str_lat
    read(IOCMT,'(a)') str_lon
    read(IOCMT,'(a)') str_dep
    close(IOCMT)

    ! output only NSIG_DIGITS of the moment tensor elements
    exp_largest = floor(log10(maxval(abs(mijn))))
    nn = exp_largest - NSIG_DIGIT
    exponent = 1.0d1 ** nn
    do i = 1, 6
       mijn(i) = floor(mijn(i)/exponent)*exponent
    enddo

    open(IOCMT,file =trim(new_cmt_file), status = 'unknown')
    write(IOCMT,'(a)') trim(pde_time)//' - GRD'
    write(IOCMT,'(a)') trim(event_name)
    write(IOCMT,'(a)') trim(str_tshift)
    write(IOCMT,'(a)') trim(str_hdur)
    write(IOCMT,'(a)') trim(str_lat)
    write(IOCMT,'(a)') trim(str_lon)
    write(IOCMT,'(a)') trim(str_dep)

    write(IOCMT,'(a,g19.6)') 'Mrr:',mijn(1)
    write(IOCMT,'(a,g19.6)') 'Mtt:',mijn(2)
    write(IOCMT,'(a,g19.6)') 'Mpp:',mijn(3)
    write(IOCMT,'(a,g19.6)') 'Mrt:',mijn(4)
    write(IOCMT,'(a,g19.6)') 'Mrp:',mijn(5)
    write(IOCMT,'(a,g19.6)') 'Mtp:',mijn(6)
    close(IOCMT)

  end subroutine write_new_cmtsolution


  !==================================================

  subroutine sdr2moment(phif,dlt,lmda,moment,mrr,mtt,mpp,mrt,mrp,mtp)

    !
    !       this subroutine converts eqk fault parameters to moment tensors
    !       for definitions see hiroo''s class notes
    !
    real phif,dlt,lmda,moment
    real mrr,mtt,mpp,mrt,mrp,mtp

    mrr=(sin(2*dlt)*sin(lmda))*moment
    mtt=(-(sin(phif))**2*sin(2*dlt)*sin(lmda)-sin(2*phif)*cos(lmda)*sin(dlt))*moment
    mpp=(-(cos(phif))**2*sin(2*dlt)*sin(lmda)+sin(2*phif)*cos(lmda)*sin(dlt))*moment
    mrt=-(cos(phif)*cos(dlt)*cos(lmda)+sin(phif)*sin(lmda)*cos(2*dlt))*moment
    mrp=-(-sin(phif)*cos(dlt)*cos(lmda)+cos(phif)*sin(lmda)*cos(2*dlt))*moment
    mtp=(-0.5*sin(2*phif)*sin(2*dlt)*sin(lmda)-cos(2*phif)*cos(lmda)*sin(dlt))*moment

  end subroutine sdr2moment

  !=========================================================================

end module grid3d_sub3
