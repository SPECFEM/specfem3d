!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine check_stability()

! computes the maximum of the norm of the displacement
! in all the slices using an MPI reduction
! and output timestamp file to check that simulation is running fine

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic

  implicit none

  double precision :: tCPU,t_remain,t_total
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total

  ! maximum of the norm of the displacement
  real(kind=CUSTOM_REAL) Usolidnorm,Usolidnorm_all ! elastic
  real(kind=CUSTOM_REAL) Usolidnormp,Usolidnormp_all ! acoustic
  real(kind=CUSTOM_REAL) Usolidnorms,Usolidnorms_all ! solid poroelastic
  real(kind=CUSTOM_REAL) Usolidnormw,Usolidnormw_all ! fluid (w.r.t.s) poroelastic

  ! norm of the backward displacement
  real(kind=CUSTOM_REAL) b_Usolidnorm, b_Usolidnorm_all
  real(kind=CUSTOM_REAL) b_Usolidnormp, b_Usolidnormp_all
  real(kind=CUSTOM_REAL) b_Usolidnorms, b_Usolidnorms_all
  real(kind=CUSTOM_REAL) b_Usolidnormw, b_Usolidnormw_all

  ! initializes
  Usolidnorm_all = 0.0_CUSTOM_REAL
  Usolidnormp_all = 0.0_CUSTOM_REAL
  Usolidnorms_all = 0.0_CUSTOM_REAL
  Usolidnormw_all = 0.0_CUSTOM_REAL

  ! compute maximum of norm of displacement in each slice
  if( ELASTIC_SIMULATION ) then
    if( GPU_MODE) then
      ! way 2: just get maximum of field from GPU
      call get_norm_elastic_from_device(Usolidnorm,Mesh_pointer,1)
    else
      Usolidnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2 + displ(3,:)**2))
    endif

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    !if(Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0.0_CUSTOM_REAL) &
    !  call exit_MPI(myrank,'single forward simulation became unstable and blew up')

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnorm,Usolidnorm_all)
  endif

  if( ACOUSTIC_SIMULATION ) then
    if(GPU_MODE) then
      ! way 2: just get maximum of field from GPU
      call get_norm_acoustic_from_device(Usolidnormp,Mesh_pointer,1)
    else
      Usolidnormp = maxval(abs(potential_dot_dot_acoustic(:)))
    endif

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnormp,Usolidnormp_all)
  endif

  if( POROELASTIC_SIMULATION ) then
    Usolidnorms = maxval(sqrt(displs_poroelastic(1,:)**2 + displs_poroelastic(2,:)**2 + &
                             displs_poroelastic(3,:)**2))
    Usolidnormw = maxval(sqrt(displw_poroelastic(1,:)**2 + displw_poroelastic(2,:)**2 + &
                             displw_poroelastic(3,:)**2))

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnorms,Usolidnorms_all)
    call max_all_cr(Usolidnormw,Usolidnormw_all)
  endif


  ! adjoint simulations
  if( SIMULATION_TYPE == 3 ) then
    ! initializes backward field norms
    b_Usolidnorm_all = 0.0_CUSTOM_REAL
    b_Usolidnormp_all = 0.0_CUSTOM_REAL
    b_Usolidnorms_all = 0.0_CUSTOM_REAL
    b_Usolidnormw_all = 0.0_CUSTOM_REAL

    if( ELASTIC_SIMULATION ) then
      ! way 2
      if(GPU_MODE) then
        call get_norm_elastic_from_device(b_Usolidnorm,Mesh_pointer,3)
      else
        b_Usolidnorm = maxval(sqrt(b_displ(1,:)**2 + b_displ(2,:)**2 + b_displ(3,:)**2))
      endif
      ! compute max of all slices
      call max_all_cr(b_Usolidnorm,b_Usolidnorm_all)
    endif
    if( ACOUSTIC_SIMULATION ) then
      ! way 2
      if(GPU_MODE) then
        call get_norm_acoustic_from_device(b_Usolidnormp,Mesh_pointer,3)
      else
        b_Usolidnormp = maxval(abs(b_potential_dot_dot_acoustic(:)))
      endif
      ! compute max of all slices
      call max_all_cr(b_Usolidnormp,b_Usolidnormp_all)
    endif
    if( POROELASTIC_SIMULATION ) then
      b_Usolidnorms = maxval(sqrt(b_displs_poroelastic(1,:)**2 + b_displs_poroelastic(2,:)**2 + &
                                  b_displs_poroelastic(3,:)**2))
      b_Usolidnormw = maxval(sqrt(b_displw_poroelastic(1,:)**2 + b_displw_poroelastic(2,:)**2 + &
                                  b_displw_poroelastic(3,:)**2))
      ! compute max of all slices
      call max_all_cr(b_Usolidnorms,b_Usolidnorms_all)
      call max_all_cr(b_Usolidnormw,b_Usolidnormw_all)
    endif
    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    !if(b_Usolidnorm > STABILITY_THRESHOLD .or. b_Usolidnorm < 0.0_CUSTOM_REAL) &
    !  call exit_MPI(myrank,'single backward simulation became unstable and blew up')
  endif

  ! user output
  if(myrank == 0) then

    write(IMAIN,*) 'Time step # ',it
    write(IMAIN,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'

    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',sngl(tCPU/dble(it))

    if( ELASTIC_SIMULATION ) &
      write(IMAIN,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all

    if( ACOUSTIC_SIMULATION ) &
      write(IMAIN,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all

    if( POROELASTIC_SIMULATION ) then
      write(IMAIN,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
      write(IMAIN,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      if( ELASTIC_SIMULATION ) &
        write(IMAIN,*) 'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
      if( ACOUSTIC_SIMULATION ) &
        write(IMAIN,*) 'Max norm pressure P (backward) in all slices (Pa) = ',b_Usolidnormp_all
      if( POROELASTIC_SIMULATION ) then
        write(IMAIN,*) 'Max norm displacement vector Us (backward) in all slices (m) = ',b_Usolidnorms_all
        write(IMAIN,*) 'Max norm displacement vector W (backward) in all slices (m) = ',b_Usolidnormw_all
      endif
    endif

    ! compute estimated remaining simulation time
    t_remain = (NSTEP - it) * (tCPU/dble(it))
    int_t_remain = int(t_remain)
    ihours_remain = int_t_remain / 3600
    iminutes_remain = (int_t_remain - 3600*ihours_remain) / 60
    iseconds_remain = int_t_remain - 3600*ihours_remain - 60*iminutes_remain
    write(IMAIN,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IMAIN,*) 'Time steps remaining = ',NSTEP - it
    write(IMAIN,*) 'Estimated remaining time in seconds = ',sngl(t_remain)
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

    ! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',sngl(t_total)
    write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'

    if(it < 100) then
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '**** BEWARE: the above time estimates are not reliable'
      write(IMAIN,*) '**** because fewer than 100 iterations have been performed'
      write(IMAIN,*) '************************************************************'
    endif
    write(IMAIN,*)

    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()

    ! write time stamp file to give information about progression of simulation
    write(outputname,"('/timestamp',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'
    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)

    if( ELASTIC_SIMULATION ) &
      write(IOUT,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all

    if( ACOUSTIC_SIMULATION ) &
      write(IOUT,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all

    if( POROELASTIC_SIMULATION ) then
      write(IOUT,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
      write(IOUT,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      if( ELASTIC_SIMULATION ) &
        write(IOUT,*) 'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
      if( ACOUSTIC_SIMULATION ) &
        write(IOUT,*) 'Max norm pressure P (backward) in all slices (Pa) = ',b_Usolidnormp_all
      if( POROELASTIC_SIMULATION ) then
        write(IOUT,*) 'Max norm displacement vector Us (backward) in all slices (m) = ',b_Usolidnorms_all
        write(IOUT,*) 'Max norm displacement vector W (backward) in all slices (m) = ',b_Usolidnormw_all
      endif
    endif

    ! estimation
    write(IOUT,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IOUT,*) 'Time steps remaining = ',NSTEP - it
    write(IOUT,*) 'Estimated remaining time in seconds = ',t_remain
    write(IOUT,"(' Estimated remaining time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain
    write(IOUT,*) 'Estimated total run time in seconds = ',t_total
    write(IOUT,"(' Estimated total run time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IOUT,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'
    close(IOUT)

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    if(Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0.0_CUSTOM_REAL &
     .or. Usolidnormp_all > STABILITY_THRESHOLD .or. Usolidnormp_all < 0.0_CUSTOM_REAL &
     .or. Usolidnorms_all > STABILITY_THRESHOLD .or. Usolidnorms_all < 0.0_CUSTOM_REAL &
     .or. Usolidnormw_all > STABILITY_THRESHOLD .or. Usolidnormw_all < 0.0_CUSTOM_REAL) &
        call exit_MPI(myrank,'forward simulation became unstable and blew up')

    ! adjoint simulations
    if( SIMULATION_TYPE == 3 ) then
      if( b_Usolidnorm_all > STABILITY_THRESHOLD .or. b_Usolidnorm_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnormp_all > STABILITY_THRESHOLD .or. b_Usolidnormp_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnorms_all > STABILITY_THRESHOLD .or. b_Usolidnorms_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnormw_all > STABILITY_THRESHOLD .or. b_Usolidnormw_all < 0.0_CUSTOM_REAL ) &
        call exit_MPI(myrank,'backward simulation became unstable and blew up')
    endif

  endif ! myrank

  end subroutine check_stability

