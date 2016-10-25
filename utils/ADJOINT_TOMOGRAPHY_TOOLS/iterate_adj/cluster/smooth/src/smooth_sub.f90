
! --------------------------------------------------------------------------------

subroutine get_all_eight_slices(ichunk,ixi,ieta, &
           ileft,iright,ibot,itop, ilb,ilt,irb,irt, &
           nproc_xi,nproc_eta)

  implicit none

  integer, intent(IN) :: ichunk,ixi,ieta,nproc_xi,nproc_eta

  integer, intent(OUT) :: ileft,iright,ibot,itop,ilb,ilt,irb,irt
  integer :: get_slice_number


  integer :: ichunk_left, islice_xi_left, islice_eta_left, &
           ichunk_right, islice_xi_right, islice_eta_right, &
           ichunk_bot, islice_xi_bot, islice_eta_bot, &
           ichunk_top, islice_xi_top, islice_eta_top, &
           ileft0,iright0,ibot0,itop0, &
           ichunk_left0, islice_xi_left0, islice_eta_left0, &
           ichunk_right0, islice_xi_right0, islice_eta_right0, &
           ichunk_bot0, islice_xi_bot0, islice_eta_bot0, &
           ichunk_top0, islice_xi_top0, islice_eta_top0


! get the first 4 immediate slices
  call get_lrbt_slices(ichunk,ixi,ieta, &
             ileft, ichunk_left, islice_xi_left, islice_eta_left, &
             iright, ichunk_right, islice_xi_right, islice_eta_right, &
             ibot, ichunk_bot, islice_xi_bot, islice_eta_bot, &
             itop, ichunk_top, islice_xi_top, islice_eta_top, &
             nproc_xi,nproc_eta)

! get the 4 diagonal neighboring slices (actually 3 diagonal slices at the corners)
  ilb = get_slice_number(ichunk,ixi-1,ieta-1,nproc_xi,nproc_eta)
  ilt = get_slice_number(ichunk,ixi-1,ieta+1,nproc_xi,nproc_eta)
  irb = get_slice_number(ichunk,ixi+1,ieta-1,nproc_xi,nproc_eta)
  irt = get_slice_number(ichunk,ixi+1,ieta+1,nproc_xi,nproc_eta)

  if (ixi == 0) then
    call get_lrbt_slices(ichunk_left,islice_xi_left,islice_eta_left, &
               ileft0, ichunk_left0, islice_xi_left0, islice_eta_left0, &
               iright0, ichunk_right0, islice_xi_right0, islice_eta_right0, &
               ibot0, ichunk_bot0, islice_xi_bot0, islice_eta_bot0, &
               itop0, ichunk_top0, islice_xi_top0, islice_eta_top0, &
               nproc_xi,nproc_eta)

    if (ichunk == 0 .or. ichunk == 1 .or. ichunk == 3 .or. ichunk == 5) then
      ilb = get_slice_number(ichunk_bot0,islice_xi_bot0,islice_eta_bot0,nproc_xi,nproc_eta)
      ilt = get_slice_number(ichunk_top0,islice_xi_top0,islice_eta_top0,nproc_xi,nproc_eta)
    else if (ichunk == 2) then
      ilb = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
      ilt = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
    else
      ilb = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
      ilt = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
    endif
  endif

  if (ixi == nproc_xi-1) then
    call get_lrbt_slices(ichunk_right,islice_xi_right,islice_eta_right, &
               ileft0, ichunk_left0, islice_xi_left0, islice_eta_left0, &
               iright0, ichunk_right0, islice_xi_right0, islice_eta_right0, &
               ibot0, ichunk_bot0, islice_xi_bot0, islice_eta_bot0, &
               itop0, ichunk_top0, islice_xi_top0, islice_eta_top0, &
               nproc_xi,nproc_eta)
    if (ichunk == 0 .or. ichunk == 1 .or. ichunk == 3 .or. ichunk == 5) then
      irb = get_slice_number(ichunk_bot0,islice_xi_bot0,islice_eta_bot0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_top0,islice_xi_top0,islice_eta_top0,nproc_xi,nproc_eta)
    else if (ichunk == 2) then
      irb = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
    else
      irb = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
    endif
  endif

  if (ieta == 0) then
    call get_lrbt_slices(ichunk_bot,islice_xi_bot,islice_eta_bot, &
               ileft0, ichunk_left0, islice_xi_left0, islice_eta_left0, &
               iright0, ichunk_right0, islice_xi_right0, islice_eta_right0, &
               ibot0, ichunk_bot0, islice_xi_bot0, islice_eta_bot0, &
               itop0, ichunk_top0, islice_xi_top0, islice_eta_top0, &
               nproc_xi,nproc_eta)
    if (ichunk == 1 .or. ichunk == 2) then
      ilb = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
      irb = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
    else if (ichunk == 3 .or. ichunk == 4) then
      ilb = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
      irb = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
    else if (ichunk == 0) then
      ilb = get_slice_number(ichunk_top0,islice_xi_top0,islice_eta_top0,nproc_xi,nproc_eta)
      irb = get_slice_number(ichunk_bot0,islice_xi_bot0,islice_eta_bot0,nproc_xi,nproc_eta)
    else
      ilb = get_slice_number(ichunk_bot0,islice_xi_bot0,islice_eta_bot0,nproc_xi,nproc_eta)
      irb = get_slice_number(ichunk_top0,islice_xi_top0,islice_eta_top0,nproc_xi,nproc_eta)
    endif
  endif

  if (ieta == nproc_eta-1) then
    call get_lrbt_slices(ichunk_top,islice_xi_top,islice_eta_top, &
               ileft0, ichunk_left0, islice_xi_left0, islice_eta_left0, &
               iright0, ichunk_right0, islice_xi_right0, islice_eta_right0, &
               ibot0, ichunk_bot0, islice_xi_bot0, islice_eta_bot0, &
               itop0, ichunk_top0, islice_xi_top0, islice_eta_top0, &
               nproc_xi,nproc_eta)

    if (ichunk == 1 .or. ichunk == 4) then
      ilt = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
    else if (ichunk == 2 .or. ichunk == 3) then
      ilt = get_slice_number(ichunk_right0,islice_xi_right0,islice_eta_right0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_left0,islice_xi_left0,islice_eta_left0,nproc_xi,nproc_eta)
    else if (ichunk == 0) then
      ilt = get_slice_number(ichunk_bot0,islice_xi_bot0,islice_eta_bot0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_top0,islice_xi_top0,islice_eta_top0,nproc_xi,nproc_eta)
    else
      ilt = get_slice_number(ichunk_top0,islice_xi_top0,islice_eta_top0,nproc_xi,nproc_eta)
      irt = get_slice_number(ichunk_bot0,islice_xi_bot0,islice_eta_bot0,nproc_xi,nproc_eta)
    endif

  endif

end subroutine get_all_eight_slices

!--------------------------------------------------------------------------------------------

subroutine get_lrbt_slices(ichunk,ixi,ieta, &
           ileft, ichunk_left, islice_xi_left, islice_eta_left, &
           iright, ichunk_right, islice_xi_right, islice_eta_right, &
           ibot, ichunk_bot, islice_xi_bot, islice_eta_bot, &
           itop, ichunk_top, islice_xi_top, islice_eta_top, &
           nproc_xi,nproc_eta)

  implicit none

  integer, intent(IN) :: ichunk, ixi, ieta, nproc_xi, nproc_eta
  integer, intent(OUT) :: ileft, ichunk_left, islice_xi_left, islice_eta_left, &
           iright, ichunk_right, islice_xi_right, islice_eta_right, &
           ibot, ichunk_bot, islice_xi_bot, islice_eta_bot, &
           itop, ichunk_top, islice_xi_top, islice_eta_top

  integer, parameter :: NCHUNKS = 6

  integer, dimension(NCHUNKS) :: chunk_left,chunk_right,chunk_bot,chunk_top, &
             slice_xi_left,slice_eta_left,slice_xi_right,slice_eta_right, &
             slice_xi_bot,slice_eta_bot,slice_xi_top,slice_eta_top
  integer :: get_slice_number

! set up mapping arrays -- assume chunk/slice number starts from 0
  chunk_left(:) = (/2,6,6,1,6,4/) - 1
  chunk_right(:) = (/4,1,1,6,1,2/) - 1
  chunk_bot(:) = (/5,5,2,5,4,5/) - 1
  chunk_top(:) = (/3,3,4,3,2,3/) - 1

  slice_xi_left(:) = (/nproc_xi-1,nproc_xi-1,nproc_xi-1-ieta,nproc_xi-1,ieta,nproc_xi-1/)
  slice_eta_left(:) = (/ieta,ieta,nproc_eta-1,ieta,0,ieta/)
  slice_xi_right(:) = (/0,0,ieta,0,nproc_xi-1-ieta,0/)
  slice_eta_right(:) = (/ieta,ieta,nproc_eta-1,ieta,0,ieta/)

  slice_xi_bot(:) = (/nproc_xi-1,ixi,ixi,nproc_xi-1-ixi,nproc_xi-1-ixi,0/)
  slice_eta_bot(:) = (/nproc_eta-1-ixi,nproc_eta-1,nproc_eta-1,0,0,ixi/)
  slice_xi_top(:) = (/nproc_xi-1,ixi,nproc_xi-1-ixi,nproc_xi-1-ixi,ixi,0/)
  slice_eta_top(:) = (/ixi,0,nproc_eta-1,nproc_eta-1,0,nproc_eta-1-ixi /)

  ichunk_left = ichunk
  ichunk_right = ichunk
  ichunk_bot = ichunk
  ichunk_top = ichunk

  islice_xi_left = ixi-1
  islice_eta_left = ieta
  islice_xi_right = ixi+1
  islice_eta_right = ieta

  islice_xi_bot = ixi
  islice_eta_bot = ieta-1
  islice_xi_top = ixi
  islice_eta_top = ieta+1

  if (ixi == 0) then
    ichunk_left=chunk_left(ichunk+1)
    islice_xi_left=slice_xi_left(ichunk+1)
    islice_eta_left=slice_eta_left(ichunk+1)
  endif
  if (ixi == nproc_xi - 1) then
    ichunk_right=chunk_right(ichunk+1)
    islice_xi_right=slice_xi_right(ichunk+1)
    islice_eta_right=slice_eta_right(ichunk+1)
  endif
  if (ieta == 0) then
    ichunk_bot=chunk_bot(ichunk+1)
    islice_xi_bot=slice_xi_bot(ichunk+1)
    islice_eta_bot=slice_eta_bot(ichunk+1)
  endif
  if (ieta == nproc_eta - 1) then
    ichunk_top=chunk_top(ichunk+1)
    islice_xi_top=slice_xi_top(ichunk+1)
    islice_eta_top=slice_eta_top(ichunk+1)
  endif

  ileft = get_slice_number(ichunk_left,islice_xi_left,islice_eta_left,nproc_xi,nproc_eta)
  iright = get_slice_number(ichunk_right,islice_xi_right,islice_eta_right,nproc_xi,nproc_eta)
  ibot = get_slice_number(ichunk_bot,islice_xi_bot,islice_eta_bot,nproc_xi,nproc_eta)
  itop = get_slice_number(ichunk_top,islice_xi_top,islice_eta_top,nproc_xi,nproc_eta)

end subroutine get_lrbt_slices

! ---------------------------------------------------------------------------------------

integer function get_slice_number(ichunk,ixi,ieta,nproc_xi,nproc_eta)

  implicit none

  integer :: ichunk, ixi, ieta, nproc_xi, nproc_eta

   get_slice_number = ichunk*nproc_xi*nproc_eta+ieta*nproc_xi+ixi

 end function get_slice_number
