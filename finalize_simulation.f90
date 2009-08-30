!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
!
! United States and French Government Sponsorship Acknowledged.

  subroutine finalize_simulation()

  use specfem_par


! save last frame

  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    open(unit=27,file=prname(1:len_trim(prname))//'save_forward_arrays.bin',status='unknown',form='unformatted')
    write(27) displ
    write(27) veloc
    write(27) accel
    if (ATTENUATION) then
      write(27) R_xx
      write(27) R_yy
      write(27) R_xy
      write(27) R_xz
      write(27) R_yz
      write(27) epsilondev_xx
      write(27) epsilondev_yy
      write(27) epsilondev_xy
      write(27) epsilondev_xz
      write(27) epsilondev_yz
    endif
    close(27)

  else if (SIMULATION_TYPE == 3) then

    ! rhop, beta, alpha kernels
! save kernels to binary files
!! DK DK removed kernels from here because not supported for CUBIT + SCOTCH yet

  endif

  if(ABSORBING_CONDITIONS .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    if (nspec2D_xmin > 0) close(31)
    if (nspec2D_xmax > 0) close(32)
    if (nspec2D_ymin > 0) close(33)
    if (nspec2D_ymax > 0) close(34)
    if (NSPEC2D_BOTTOM > 0) close(35)
  endif

  if (nrec_local > 0) then
    if (.not. (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3)) then
!      call write_adj_seismograms(myrank,seismograms_d,number_receiver_global, &
!          nrec_local,it,DT,NSTEP,t0,LOCAL_PATH,1)
      call write_adj_seismograms2(myrank,seismograms_eps,number_receiver_global, &
            nrec_local,it,DT,NSTEP,t0,LOCAL_PATH)
      do irec_local = 1, nrec_local
        write(outputname,'(a,i5.5)') 'OUTPUT_FILES/src_frechet.',number_receiver_global(irec_local)
        open(unit=27,file=trim(outputname),status='unknown')
!
! r -> z, theta -> -y, phi -> x
!
!  Mrr =  Mzz
!  Mtt =  Myy
!  Mpp =  Mxx
!  Mrt = -Myz
!  Mrp =  Mxz
!  Mtp = -Mxy

        write(27,*) Mzz_der(irec_local)
        write(27,*) Myy_der(irec_local)
        write(27,*) Mxx_der(irec_local)
        write(27,*) -Myz_der(irec_local)
        write(27,*) Mxz_der(irec_local)
        write(27,*) -Mxy_der(irec_local)
        write(27,*) sloc_der(1,irec_local)
        write(27,*) sloc_der(2,irec_local)
        write(27,*) sloc_der(3,irec_local)
        close(27)
      enddo
    endif
  endif


! close the main output file
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of the simulation'
    write(IMAIN,*)
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call sync_all()

  end subroutine