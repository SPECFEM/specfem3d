
 program unroll_illog

!! DK DK writes a simple code to unroll some loops

 implicit none

  integer :: i,j,k

!! DK DK the arrays to unroll are declared as tabg0_small_buffer(2, 6, -2:2, llog0)

  integer :: illog

  do k = -2,2
  do j = 1,6
  do i = 1,2
    print *,'tabg0_small_buffer(illog,',i,',',j,',',k,') = tabg0(illog,',i,',',j,',',k,',ir_)'
  enddo
  enddo
  enddo

!   tabg0_small_buffer(illog,:,:,:) = tabg0der(illog,:,:,:,ir_)

 end program unroll_illog

