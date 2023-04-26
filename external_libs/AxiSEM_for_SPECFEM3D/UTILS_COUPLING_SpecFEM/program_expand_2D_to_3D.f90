  program expand_2D_to_3D

!
!--- IMPORTANT THE AXISEM OUTPUT IS IN SINGLE PRECISION ---!
!

  !use interp_mod        !! interpolation
  use mpi_mod
  use reading_inputs     !! reading input parameters
  use reading_field      !! reading field
  use post_processing    !! do the post-process stuff (cf Tarje's papers)
  !use writing_output    !! write the output interpoled wave field
  integer isim

!
!---
!

  call init_mpi()

  isim=1
  ! read point where we need to compute solution
  if (myrank == 0) then

    write(*,*) 'Expand 2D to 3D'
    write(*,*)
    write(*,*) 'reading inputs...'
    write(*,*)

    call read_info_simu(nsim)
    call read_inputs(isim)

    write(*,*) 'reading mesh points...'
    write(*,*)

    ! read GLL point coordinate
    call read_mesh(isim)

    write(*,*) 'connecting points...'
    write(*,*)

    ! find elements and coordinate for each input point
    call connect_points()

  endif

  call barrier_mpi()

  !debug
  !write(*,*) 'Before mpi', myrank

  call alloc_all_mpi()
  call bcast_all_mpi()

  !debug
  !write(*,*) 'After alloc and bcast ', myrank

  call distrib_mpi()  !! to do : faire directement la distrib sur les memes procs que Specfem

  !debug
  !write(*,*) 'After mpi', myrank

!
!---
!

  if (.not. recip_KH_integral) then

    if (myrank == 0) then
      write(*,*) 'interpolating velocity and stress fields...'
      write(*,*)
    endif

    do isim=1,nsim  !! do to : mettre en memoire la solution sous echantillonnee et la resampler avant de l'ecrire

      ! interpolation of the velocity field in each point
      call read_veloc_field_and_interpol(isim)
      ! interpolation of the stress in each point
      call read_stress_field_and_interpol(isim)

    enddo

  else

    if (myrank == 0) then
      write(*,*) 'interpolating displ and Pderiv fields...'
      write(*,*)
    endif

    call displ_read_recombine_interpol_rotate()
    call pderiv_read_recombine_interpol_rotate()

  endif

!
!---
!

  call finalize_mpi()

  end program expand_2D_to_3D
