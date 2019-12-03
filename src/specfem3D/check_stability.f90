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


  subroutine check_stability()

! computes the maximum of the norm of the displacement
! in all the slices using an MPI reduction
! and output timestamp file to check that simulation is running fine

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic

  implicit none

  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total
  ! timing
  double precision :: t_remain,t_total
  double precision :: tCPU
  double precision, external :: wtime

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

  ! to determine date and time at which the run will finish
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values

  character(len=MAX_STRING_LEN) :: outputname

  character(len=3), dimension(12) :: month_name
  character(len=3), dimension(0:6) :: weekday_name
  data month_name /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
  data weekday_name /'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/
  integer :: year,mon,day,hr,minutes,timestamp,julian_day_number,day_of_week, &
             timestamp_remote,year_remote,mon_remote,day_remote,hr_remote,minutes_remote,day_of_week_remote
  integer, external :: idaywk

  ! initializes
  Usolidnorm_all = 0.0_CUSTOM_REAL
  Usolidnormp_all = 0.0_CUSTOM_REAL
  Usolidnorms_all = 0.0_CUSTOM_REAL
  Usolidnormw_all = 0.0_CUSTOM_REAL

  ! compute maximum of norm of displacement in each slice
  if (ELASTIC_SIMULATION) then
    if (GPU_MODE) then
      ! way 2: just get maximum of field from GPU
      call get_norm_elastic_from_device(Usolidnorm,Mesh_pointer,1)
    else
      Usolidnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2 + displ(3,:)**2))
    endif

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
    !if (Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0.0_CUSTOM_REAL) &
    !  call exit_MPI(myrank,'single forward simulation became unstable and blew up')

    ! checks first entry for Not-a-Number (NaN) value
! this trick checks for NaN (Not a Number), which is not even equal to itself
    if (displ(1,iglob_check_elastic) /= displ(1,iglob_check_elastic)) then
      call exit_MPI(myrank,'forward simulation became unstable in elastic domain and blew up')
    endif

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnorm,Usolidnorm_all)
  endif

  if (ACOUSTIC_SIMULATION) then
    if (GPU_MODE) then
      ! way 2: just get maximum of field from GPU
      call get_norm_acoustic_from_device(Usolidnormp,Mesh_pointer,1)
    else
      Usolidnormp = maxval(abs(potential_dot_dot_acoustic(:)))
    endif

    ! debug
    !if (myrank == 0) print *,'norm p = ',Usolidnormp,maxval(potential_dot_dot_acoustic(:)),potential_dot_dot_acoustic(1)
    ! note for gfortran compiler: only potential_dot_dot(.) entry returns NaN, maxval(..) and maxval(abs(..)) return 0.0

    ! checks first entry for Not-a-Number (NaN) value
    if (potential_dot_dot_acoustic(iglob_check_acoustic) /= potential_dot_dot_acoustic(iglob_check_acoustic)) then
      call exit_MPI(myrank,'forward simulation became unstable in acoustic domain and blew up')
    endif

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnormp,Usolidnormp_all)
  endif

  if (POROELASTIC_SIMULATION) then
    Usolidnorms = maxval(sqrt(displs_poroelastic(1,:)**2 + displs_poroelastic(2,:)**2 + &
                             displs_poroelastic(3,:)**2))
    Usolidnormw = maxval(sqrt(displw_poroelastic(1,:)**2 + displw_poroelastic(2,:)**2 + &
                             displw_poroelastic(3,:)**2))

    ! checks first entry for Not-a-Number (NaN) value
    if (displs_poroelastic(1,iglob_check_poroelastic) /= displs_poroelastic(1,iglob_check_poroelastic)) then
      call exit_MPI(myrank,'forward simulation became unstable in poroelastic domain and blew up')
    endif

    ! compute the maximum of the maxima for all the slices using an MPI reduction
    call max_all_cr(Usolidnorms,Usolidnorms_all)
    call max_all_cr(Usolidnormw,Usolidnormw_all)
  endif


  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
    ! initializes backward field norms
    b_Usolidnorm_all = 0.0_CUSTOM_REAL
    b_Usolidnormp_all = 0.0_CUSTOM_REAL
    b_Usolidnorms_all = 0.0_CUSTOM_REAL
    b_Usolidnormw_all = 0.0_CUSTOM_REAL

    if (ELASTIC_SIMULATION) then
      ! way 2
      if (GPU_MODE) then
        call get_norm_elastic_from_device(b_Usolidnorm,Mesh_pointer,3)
      else
        b_Usolidnorm = maxval(sqrt(b_displ(1,:)**2 + b_displ(2,:)**2 + b_displ(3,:)**2))
      endif
      ! compute max of all slices
      call max_all_cr(b_Usolidnorm,b_Usolidnorm_all)
    endif
    if (ACOUSTIC_SIMULATION) then
      ! way 2
      if (GPU_MODE) then
        call get_norm_acoustic_from_device(b_Usolidnormp,Mesh_pointer,3)
      else
        b_Usolidnormp = maxval(abs(b_potential_dot_dot_acoustic(:)))
      endif
      ! compute max of all slices
      call max_all_cr(b_Usolidnormp,b_Usolidnormp_all)
    endif
    if (POROELASTIC_SIMULATION) then
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
    !if (b_Usolidnorm > STABILITY_THRESHOLD .or. b_Usolidnorm < 0.0_CUSTOM_REAL) &
    !  call exit_MPI(myrank,'single backward simulation became unstable and blew up')
  endif

  ! user output
  if (myrank == 0) then

    write(IMAIN,*) 'Time step # ',it
    write(IMAIN,*) 'Time: ',sngl((it-1)*DT-t0),' seconds'

    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IMAIN,*) 'Mean elapsed time per time step in seconds = ',sngl(tCPU/dble(it))

    if (ELASTIC_SIMULATION) &
      write(IMAIN,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all

    if (ACOUSTIC_SIMULATION) &
      write(IMAIN,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all

    if (POROELASTIC_SIMULATION) then
      write(IMAIN,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
      write(IMAIN,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      if (ELASTIC_SIMULATION) &
        write(IMAIN,*) 'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
      if (ACOUSTIC_SIMULATION) &
        write(IMAIN,*) 'Max norm pressure P (backward) in all slices (Pa) = ',b_Usolidnormp_all
      if (POROELASTIC_SIMULATION) then
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
    write(IMAIN,"(' Estimated remaining time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain

    ! compute estimated total simulation time
    t_total = t_remain + tCPU
    int_t_total = int(t_total)
    ihours_total = int_t_total / 3600
    iminutes_total = (int_t_total - 3600*ihours_total) / 60
    iseconds_total = int_t_total - 3600*ihours_total - 60*iminutes_total
    write(IMAIN,*) 'Estimated total run time in seconds = ',sngl(t_total)
    write(IMAIN,"(' Estimated total run time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IMAIN,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'

    if (it < NSTEP) then

      ! get current date
      call date_and_time(datein,timein,zone,time_values)
      ! time_values(1): year
      ! time_values(2): month of the year
      ! time_values(3): day of the month
      ! time_values(5): hour of the day
      ! time_values(6): minutes of the hour

      ! compute date at which the run should finish; for simplicity only minutes
      ! are considered, seconds are ignored; in any case the prediction is not
      ! accurate down to seconds because of system and network fluctuations
      year = time_values(1)
      mon = time_values(2)
      day = time_values(3)
      hr = time_values(5)
      minutes = time_values(6)

      ! get timestamp in minutes of current date and time
      call convtime(timestamp,year,mon,day,hr,minutes)

      ! add remaining minutes
      timestamp = timestamp + nint(t_remain / 60.d0)

      ! get date and time of that future timestamp in minutes
      call invtime(timestamp,year,mon,day,hr,minutes)

      ! convert to Julian day to get day of the week
      call calndr(day,mon,year,julian_day_number)
      day_of_week = idaywk(julian_day_number)

      write(IMAIN,"(' The run will finish approximately on (in local time): ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
          weekday_name(day_of_week),month_name(mon),day,year,hr,minutes

      ! print date and time estimate of end of run in another country.
      ! For instance: the code runs at Caltech in California but the person
      ! running the code is connected remotely from France, which has 9 hours more
      if (ADD_TIME_ESTIMATE_ELSEWHERE .and. HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE /= 0) then

        ! add time difference with that remote location (can be negative)
        timestamp_remote = timestamp + HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE

        ! get date and time of that future timestamp in minutes
        call invtime(timestamp_remote,year_remote,mon_remote,day_remote,hr_remote,minutes_remote)

        ! convert to Julian day to get day of the week
        call calndr(day_remote,mon_remote,year_remote,julian_day_number)
        day_of_week_remote = idaywk(julian_day_number)

        if (HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE > 0) then
          write(IMAIN,*) 'Adding positive time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
        else
          write(IMAIN,*) 'Adding negative time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
        endif
        write(IMAIN,*) 'and ',abs(MINUTES_TIME_DIFFERENCE),' minutes to get estimate at a remote location'
        write(IMAIN, &
            "(' The run will finish approximately on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
            weekday_name(day_of_week_remote),month_name(mon_remote),day_remote,year_remote,hr_remote,minutes_remote
      endif

      if (it < 100) then
        write(IMAIN,*) '************************************************************'
        write(IMAIN,*) '**** BEWARE: the above time estimates are not very reliable'
        write(IMAIN,*) '**** because fewer than 100 iterations have been performed'
        write(IMAIN,*) '************************************************************'
      endif

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
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)

    if (ELASTIC_SIMULATION) &
      write(IOUT,*) 'Max norm displacement vector U in all slices (m) = ',Usolidnorm_all

    if (ACOUSTIC_SIMULATION) &
      write(IOUT,*) 'Max norm pressure P in all slices (Pa) = ',Usolidnormp_all

    if (POROELASTIC_SIMULATION) then
      write(IOUT,*) 'Max norm displacement vector Us in all slices (m) = ',Usolidnorms_all
      write(IOUT,*) 'Max norm displacement vector W in all slices (m) = ',Usolidnormw_all
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      if (ELASTIC_SIMULATION) &
        write(IOUT,*) 'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
      if (ACOUSTIC_SIMULATION) &
        write(IOUT,*) 'Max norm pressure P (backward) in all slices (Pa) = ',b_Usolidnormp_all
      if (POROELASTIC_SIMULATION) then
        write(IOUT,*) 'Max norm displacement vector Us (backward) in all slices (m) = ',b_Usolidnorms_all
        write(IOUT,*) 'Max norm displacement vector W (backward) in all slices (m) = ',b_Usolidnormw_all
      endif
    endif

    ! estimation
    write(IOUT,*) 'Time steps done = ',it,' out of ',NSTEP
    write(IOUT,*) 'Time steps remaining = ',NSTEP - it
    write(IOUT,*) 'Estimated remaining time in seconds = ',t_remain
    write(IOUT,"(' Estimated remaining time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_remain,iminutes_remain,iseconds_remain
    write(IOUT,*) 'Estimated total run time in seconds = ',t_total
    write(IOUT,"(' Estimated total run time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
             ihours_total,iminutes_total,iseconds_total
    write(IOUT,*) 'We have done ',sngl(100.d0*dble(it)/dble(NSTEP)),'% of that'

  if (it < NSTEP) then

    write(IOUT,"(' The run will finish approximately on (in local time): ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
        weekday_name(day_of_week),month_name(mon),day,year,hr,minutes

    ! print date and time estimate of end of run in another country.
    ! For instance: the code runs at Caltech in California but the person
    ! running the code is connected remotely from France, which has 9 hours more
    if (ADD_TIME_ESTIMATE_ELSEWHERE .and. HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE /= 0) then
      if (HOURS_TIME_DIFFERENCE * 60 + MINUTES_TIME_DIFFERENCE > 0) then
        write(IOUT,*) 'Adding positive time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
      else
        write(IOUT,*) 'Adding negative time difference of ',abs(HOURS_TIME_DIFFERENCE),' hours'
      endif
      write(IOUT,*) 'and ',abs(MINUTES_TIME_DIFFERENCE),' minutes to get estimate at a remote location'
      write(IOUT, &
          "(' The run will finish approximately on (in remote time): ',a3,' ',a3,' ',i2.2,', ',i4.4,' ',i2.2,':',i2.2)") &
          weekday_name(day_of_week_remote),month_name(mon_remote), &
          day_remote,year_remote,hr_remote,minutes_remote
    endif

    if (it < 100) then
      write(IOUT,*)
      write(IOUT,*) '************************************************************'
      write(IOUT,*) '**** BEWARE: the above time estimates are not very reliable'
      write(IOUT,*) '**** because fewer than 100 iterations have been performed'
      write(IOUT,*) '************************************************************'
    endif

  endif

    close(IOUT)

    ! check stability of the code, exit if unstable
    ! negative values can occur with some compilers when the unstable value is greater
    ! than the greatest possible floating-point number of the machine
! this trick checks for NaN (Not a Number), which is not even equal to itself
    if (Usolidnorm_all > STABILITY_THRESHOLD .or. Usolidnorm_all < 0.0_CUSTOM_REAL .or. Usolidnorm_all /= Usolidnorm_all &
     .or. Usolidnormp_all > STABILITY_THRESHOLD .or. Usolidnormp_all < 0.0_CUSTOM_REAL .or. Usolidnormp_all /= Usolidnormp_all &
     .or. Usolidnorms_all > STABILITY_THRESHOLD .or. Usolidnorms_all < 0.0_CUSTOM_REAL .or. Usolidnorms_all /= Usolidnorms_all &
     .or. Usolidnormw_all > STABILITY_THRESHOLD .or. Usolidnormw_all < 0.0_CUSTOM_REAL .or. Usolidnormw_all /= Usolidnormw_all) &
        call exit_MPI(myrank,'forward simulation became unstable and blew up')

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
! this trick checks for NaN (Not a Number), which is not even equal to itself
      if (b_Usolidnorm_all > STABILITY_THRESHOLD .or. b_Usolidnorm_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnorm_all /= b_Usolidnorm_all &
        .or. b_Usolidnormp_all > STABILITY_THRESHOLD .or. b_Usolidnormp_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnormp_all /= b_Usolidnormp_all &
        .or. b_Usolidnorms_all > STABILITY_THRESHOLD .or. b_Usolidnorms_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnorms_all /= b_Usolidnorms_all &
        .or. b_Usolidnormw_all > STABILITY_THRESHOLD .or. b_Usolidnormw_all < 0.0_CUSTOM_REAL &
        .or. b_Usolidnormw_all /= b_Usolidnormw_all) &
          call exit_MPI(myrank,'backward simulation became unstable and blew up')
    endif

  endif ! myrank

  end subroutine check_stability

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_stability_backward()

! only for backward/reconstructed wavefield

  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic

  implicit none

  ! norm of the backward displacement
  real(kind=CUSTOM_REAL) b_Usolidnorm, b_Usolidnorm_all
  real(kind=CUSTOM_REAL) b_Usolidnormp, b_Usolidnormp_all
  real(kind=CUSTOM_REAL) b_Usolidnorms, b_Usolidnorms_all
  real(kind=CUSTOM_REAL) b_Usolidnormw, b_Usolidnormw_all

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3 ) return

  ! initializes backward field norms
  b_Usolidnorm_all = 0.0_CUSTOM_REAL
  b_Usolidnormp_all = 0.0_CUSTOM_REAL
  b_Usolidnorms_all = 0.0_CUSTOM_REAL
  b_Usolidnormw_all = 0.0_CUSTOM_REAL

  if (ELASTIC_SIMULATION) then
    ! way 2
    if (GPU_MODE) then
      call get_norm_elastic_from_device(b_Usolidnorm,Mesh_pointer,3)
    else
      b_Usolidnorm = maxval(sqrt(b_displ(1,:)**2 + b_displ(2,:)**2 + b_displ(3,:)**2))
    endif
    ! compute max of all slices
    call max_all_cr(b_Usolidnorm,b_Usolidnorm_all)
  endif
  if (ACOUSTIC_SIMULATION) then
    ! way 2
    if (GPU_MODE) then
      call get_norm_acoustic_from_device(b_Usolidnormp,Mesh_pointer,3)
    else
      b_Usolidnormp = maxval(abs(b_potential_dot_dot_acoustic(:)))
    endif
    ! compute max of all slices
    call max_all_cr(b_Usolidnormp,b_Usolidnormp_all)
  endif
  if (POROELASTIC_SIMULATION) then
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
  !if (b_Usolidnorm > STABILITY_THRESHOLD .or. b_Usolidnorm < 0.0_CUSTOM_REAL) &
  !  call exit_MPI(myrank,'single backward simulation became unstable and blew up')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Time step for back propagation # ',it
    if (ELASTIC_SIMULATION) &
      write(IMAIN,*) 'Max norm displacement vector U (backward) in all slices (m) = ',b_Usolidnorm_all
    if (ACOUSTIC_SIMULATION) &
      write(IMAIN,*) 'Max norm pressure P (backward) in all slices (Pa) = ',b_Usolidnormp_all
    if (POROELASTIC_SIMULATION) then
      write(IMAIN,*) 'Max norm displacement vector Us (backward) in all slices (m) = ',b_Usolidnorms_all
      write(IMAIN,*) 'Max norm displacement vector W (backward) in all slices (m) = ',b_Usolidnormw_all
    endif
    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()

    ! check stability of the code, exit if unstable
    ! this trick checks for NaN (Not a Number), which is not even equal to itself
    if (b_Usolidnorm_all > STABILITY_THRESHOLD .or. b_Usolidnorm_all < 0.0_CUSTOM_REAL &
      .or. b_Usolidnorm_all /= b_Usolidnorm_all &
      .or. b_Usolidnormp_all > STABILITY_THRESHOLD .or. b_Usolidnormp_all < 0.0_CUSTOM_REAL &
      .or. b_Usolidnormp_all /= b_Usolidnormp_all &
      .or. b_Usolidnorms_all > STABILITY_THRESHOLD .or. b_Usolidnorms_all < 0.0_CUSTOM_REAL &
      .or. b_Usolidnorms_all /= b_Usolidnorms_all &
      .or. b_Usolidnormw_all > STABILITY_THRESHOLD .or. b_Usolidnormw_all < 0.0_CUSTOM_REAL &
      .or. b_Usolidnormw_all /= b_Usolidnormw_all) &
        call exit_MPI(myrank,'backward simulation became unstable and blew up')

  endif ! myrank

  end subroutine check_stability_backward

!
!-------------------------------------------------------------------------------------------------
!


  subroutine it_print_elapsed_time()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  ! local parameters
  integer :: ihours,iminutes,iseconds,int_tCPU
  ! timing
  double precision :: tCPU
  double precision, external :: wtime

  if (myrank == 0) then
    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start

    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Time loop finished. Timing info:'
    write(IMAIN,*) 'Total elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Total elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    call flush_IMAIN()
  endif

  end subroutine it_print_elapsed_time

