program source_to_conv

  use precision_mod

  implicit none

  integer(kind=si) :: nt, it
  real(kind=cp)    :: dt, to, tmp, alpha, t, divis
  real(kind=cp), parameter :: pi=3.14159, decay=1.628, shift=1.5

  real(kind=cp), dimension(:), allocatable :: s

  write(*,*)'type,hdur,dt'
  read(*,*)type,hdur,dt


  decay = 1.628
  nt = ceil(2.* shift * hdur / dt) +1
  t  = (0:nt-1).*dt;

  alpha  = decay / hd;
  divis = 1./sqrt(pi);
  to = (nt-1)*dt/2

  allocate(s(nt))


  do i = 1, nt
     t = real((it-1)*dt)
     temp = alpha .* (t - to)
     if (type == 0) then
        s(i) = alpha * divis * exp(-temp**2)
     else if (type == 1)
        s(i) = alpha * divis * exp(-temp**2) * (-2 * alpha**2 * t)
     else if (type == 3)
        s(i) = alpha * divis * exp(-temp**2) *  (-2 * alpha**2 * t**2)
     endif
  enddo

  open(10,file='fsource.bin',access='direct',recl=cp*nt)
  write(10,rec=1)amp*s
  close(10)

  deallocate(s)
