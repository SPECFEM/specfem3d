
! permute order of points in CUBIT elements (hexahedra) to display them in OpenDX

  program permute_cubit_2_opendx

  implicit none

  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8,ispec,nspec

! number of elements in CUBIT mesh
  nspec = 13

  do ispec=1,nspec

    read(*,*) iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8

! write to AVS or DX global file with correct offset for hexahedra (3-D)
! in the case of OpenDX, node numbers start at zero
! in the case of AVS, node numbers start at one
! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
    write(*,"(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
       iglob4-1,iglob1-1,iglob8-1,iglob5-1,iglob3-1,iglob2-1,iglob7-1,iglob6-1

  enddo

  end program permute_cubit_2_opendx

