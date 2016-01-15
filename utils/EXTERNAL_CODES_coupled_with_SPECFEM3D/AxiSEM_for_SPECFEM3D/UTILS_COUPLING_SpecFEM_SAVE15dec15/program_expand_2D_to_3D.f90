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
  if (myrank==0) then

    call read_info_simu(nsim)
    call read_inputs(isim)

    ! read gll point coordinate
    call read_mesh(isim)

    ! find elements and coordinate for each input point
    call connect_points()

  endif

  call barrier_mpi()
  write(*,*) 'Before mpi', myrank

  call alloc_all_mpi()
  call bcast_all_mpi()
  write(*,*) 'After alloc and bcast ', myrank

  call distrib_mpi()  !! to do : faire directement la distrib sur les memes procs que Specfem
  write(*,*) 'After mpi', myrank

!
!---
!

  do isim=1,nsim  !! do to : mettre en memoire la solution sous echantillonnee et la resampler avant de l'ecrire

    if (.not. recip_KH_integral) then

      ! interpolation of the velocity field in each point
      call read_veloc_field_and_interpol(isim)
      ! interpolation of the stress in each point
      call read_stress_field_and_interpol(isim)

    else

!!      ! interpolation of the stress in each point
!!      call read_stress_field_and_interpol(isim)

      ! interpolation of the displacement field in each point
      call read_displ_field_and_interpol(isim)
      ! interpolation of the partial derivatives (one by one) field in each point
      call read_partialderivatives_field_and_interpol(isim)

    endif

  enddo

  call finalize_mpi()

  end program expand_2D_to_3D
