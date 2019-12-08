!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


  subroutine save_forward_arrays()

! save files to local disk or tape system if restart file

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  integer :: ier

  ! checks if anything to do
  if (SIMULATION_TYPE /= 1) return
  if (UNDO_ATTENUATION_AND_OR_PML) return

  ! saves last frame
  if (ADIOS_FOR_FORWARD_ARRAYS) then
    call save_forward_arrays_adios()
  else
    open(unit=IOUT,file=prname(1:len_trim(prname))//'save_forward_arrays.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error: opening save_forward_arrays.bin'
      print *,'path: ',prname(1:len_trim(prname))//'save_forward_arrays.bin'
      call exit_mpi(myrank,'error opening file save_forward_arrays.bin')
    endif

    if (ACOUSTIC_SIMULATION) then
      write(IOUT) potential_acoustic
      write(IOUT) potential_dot_acoustic
      write(IOUT) potential_dot_dot_acoustic
    endif

    if (ELASTIC_SIMULATION) then
      write(IOUT) displ
      write(IOUT) veloc
      write(IOUT) accel

      if (ATTENUATION) then
        write(IOUT) R_trace
        write(IOUT) R_xx
        write(IOUT) R_yy
        write(IOUT) R_xy
        write(IOUT) R_xz
        write(IOUT) R_yz
        write(IOUT) epsilondev_trace
        write(IOUT) epsilondev_xx
        write(IOUT) epsilondev_yy
        write(IOUT) epsilondev_xy
        write(IOUT) epsilondev_xz
        write(IOUT) epsilondev_yz
      endif
    endif

    if (POROELASTIC_SIMULATION) then
      write(IOUT) displs_poroelastic
      write(IOUT) velocs_poroelastic
      write(IOUT) accels_poroelastic
      write(IOUT) displw_poroelastic
      write(IOUT) velocw_poroelastic
      write(IOUT) accelw_poroelastic
    endif

    close(IOUT)
  endif

  end subroutine save_forward_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_forward_arrays_undoatt()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! transfers wavefields from GPU device to CPU host
  if (GPU_MODE) then
    ! acoustic potentials
    if (ACOUSTIC_SIMULATION) &
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                          Mesh_pointer)

    ! elastic wavefield
    if (ELASTIC_SIMULATION) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)

      if (ATTENUATION) &
        call transfer_fields_att_from_device(Mesh_pointer, &
                                             R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                                             epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                             R_trace,epsilondev_trace, &
                                             size(epsilondev_xx))
    endif
  endif

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  ! saves frame of the forward simulation

  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
  outputname = trim(LOCAL_PATH)//'/'//trim(outputname)

  ! outputs to file
  open(unit=IOUT,file=trim(outputname),status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for writing')

  if (ACOUSTIC_SIMULATION) then
    write(IOUT) potential_acoustic
    write(IOUT) potential_dot_acoustic
    write(IOUT) potential_dot_dot_acoustic
  endif

  if (ELASTIC_SIMULATION) then
    write(IOUT) displ
    write(IOUT) veloc
    write(IOUT) accel
    if (ATTENUATION) then
      write(IOUT) R_trace
      write(IOUT) R_xx
      write(IOUT) R_yy
      write(IOUT) R_xy
      write(IOUT) R_xz
      write(IOUT) R_yz
      write(IOUT) epsilondev_trace
      write(IOUT) epsilondev_xx
      write(IOUT) epsilondev_yy
      write(IOUT) epsilondev_xy
      write(IOUT) epsilondev_xz
      write(IOUT) epsilondev_yz
    endif
  endif

  close(IOUT)

  end subroutine save_forward_arrays_undoatt
