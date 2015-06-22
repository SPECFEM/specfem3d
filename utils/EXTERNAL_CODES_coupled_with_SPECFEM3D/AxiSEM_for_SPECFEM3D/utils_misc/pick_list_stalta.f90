program pick_sismo

  integer ind,i,m,di
  integer nlta,nsta
  real(kind=4), allocatable :: sig(:),stalta(:),t(:)
  real thres,x,y
  character(len=250) filename

  nlta=100
  nsta=20
  thres=0.1
  di=40
  open(20,file='pointe.txt')
  open(11,file='list_of_Cz_to_read.txt')
  do !! loop over name of Z component to pick

     read(11,'(a)',end=199) filename
     open(10,file=trim(filename))

     ! read z component
     m=0
     do
        read(10,*,end=99) x,y
        m=m+1
     enddo
99   continue
     close(10)

     allocate(sig(m),stalta(m),t(m))

     open(10,file=trim(filename))
     do i=1,m
        read(10,*) t(i), sig(i)
     enddo
     close(10)


     call substalta(sig, nsta, nlta, stalta,m)
     !! TO DO new way to pick call this subroutine !call pick(stalta,ntime,ipick,thres)
     do i=1,m
        if (stalta(i) >= thres) then
           ind = i;
           exit
        endif
     enddo
     open(10,file='stalta.txt')
     do i=1,m
        write(10,*) stalta(i)
     enddo
     close(10)

    ! write pick
     write(20,*) trim(filename), t(ind), t(ind+di), ind, ind+di

     deallocate(sig,stalta,t)
  enddo

  199 continue
  close(20)
  close(11)

end program pick_sismo

!================================================================================
! Compute STA/LTA for picking
!!$  subroutine substalta(sig, nsta, nlta, stalta,m)
!!$
!!$    real(kind=4) sig(m)
!!$    real(kind=4) stalta(m)
!!$
!!$    integer, intent(in) :: nsta, nlta
!!$
!!$    integer :: m, nsta_1, nlta_1, i
!!$
!!$    real(kind=4), dimension(:), allocatable :: sta, lta, pad_sta, pad_lta, tmp1, tmp2
!!$
!!$    !m = size(sig)
!!$
!!$    nsta_1 = nsta - 1
!!$    nlta_1 = nlta - 1
!!$    !write(*,*) m,nsta_1,nlta_1
!!$    allocate(sta(m))
!!$    allocate(lta(m))
!!$    allocate(tmp1(m))
!!$    allocate(tmp2(m))
!!$    allocate(pad_sta(nsta_1))
!!$    allocate(pad_lta(nlta_1))
!!$    sta = 0.
!!$    lta = 0.
!!$    pad_sta = 0.
!!$    pad_lta = 1.
!!$
!!$    !*** compute the short time average (STA)
!!$    do i=1,nsta
!!$       tmp1(1:nsta_1) = pad_sta(:)
!!$       tmp1(nsta_1+1:m) = sig(i:m - nsta_1 + i-1)**2
!!$       sta = sta + tmp1
!!$    enddo
!!$    sta = sta / nsta
!!$
!!$    !*** compute the long time average (LTA)
!!$    do i =1,nlta
!!$       tmp2(1:nlta_1) = pad_lta(:)
!!$       tmp2(nlta_1+1:m) = sig(i:m - nlta_1 + i-1)**2
!!$       lta = lta + tmp2
!!$    enddo
!!$    lta = lta / nlta
!!$
!!$    sta(1:nsta_1) = 0.
!!$    lta(1:nlta_1) = 1.
!!$
!!$    do i=1,m
!!$       if (lta(i) < 1e-10) then
!!$          lta(i) = 1.
!!$          sta(i) = 0.
!!$       endif
!!$    enddo
!!$    stalta = sta / lta
!!$
!!$  end subroutine substalta

  subroutine pick(stalta,n,i,thres)
    !use global_parameters
    integer i,n
    real(kind=4) stalta(n),thres

    do i=1,n
       if (stalta(i) >= thres) exit
    enddo

  end subroutine pick

!================================================================================
!! Compute STA/LTA for picking

  subroutine substalta(sig, nsta, nlta, stalta,m)

    !use global_parameters
    integer m
    real(kind=4) sig(m)
    real(kind=4) stalta(m)

    integer, intent(in) :: nsta, nlta

    integer :: nsta_1, nlta_1, i,j

    real(kind=4), dimension(:), allocatable :: sta, lta, pad_sta, pad_lta, tmp1,tmp2

    !m = size(sig)

    nsta_1 = nsta - 1
    nlta_1 = nlta - 1
    !write(*,*) m,nsta_1,nlta_1
    allocate(sta(m))
    allocate(lta(m))
    allocate(tmp1(m))
    allocate(tmp2(m))
    allocate(pad_sta(nsta_1))
    allocate(pad_lta(nlta_1))
    sta = 0.
    lta = 0.
    pad_sta = 0.
    pad_lta = 1.

    !*** compute the short time average (STA)
    !do i=1,nsta
    !   tmp1(1:nsta_1) = pad_sta(:)
    !   tmp1(nsta_1+1:m) = sig(i:m - nsta_1 + i-1)**2
    !   sta = sta + tmp1
    !enddo


    do i=nsta,m
       do j=i-nsta_1,i
         sta(i) = sta(i) + sig(j)**2
       enddo
    enddo
    sta = sta / nsta

    !*** compute the long time average (LTA)
    !do i =1,nlta
    !   tmp2(1:nlta_1) = pad_lta(:)
    !   tmp2(nlta_1+1:m) = sig(i:m - nlta_1 + i-1)**2
    !   lta = lta + tmp2
    !enddo
    do i=nlta,m
      do j=i-nlta_1,i
          lta(i)=lta(i)+sig(j)**2
      enddo
    enddo
    lta = lta / nlta

    sta(1:nsta_1) = 0.
    lta(1:nlta_1) = 1.

    do i=1,m
       if (lta(i) < 1e-10) then
          lta(i) = 1.
          sta(i) = 0.
       endif
    enddo
    stalta = sta / lta
    !stalta = lta
  end subroutine substalta

