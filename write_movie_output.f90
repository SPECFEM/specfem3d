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

  subroutine write_movie_output()

  use specfem_par
  use specfem_par_movie  
  implicit none

  ! shakemap creation
  if (EXTERNAL_MESH_CREATE_SHAKEMAP) then
    call wmo_create_shakemap_em()
  endif 

  ! movie file creation
  if(EXTERNAL_MESH_MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    call wmo_create_movie_surface_em()
  endif

  ! saves MOVIE on the SURFACE
  if(MOVIE_SURFACE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    call wmo_movie_surface_output_o()
  endif

  ! computes SHAKING INTENSITY MAP
  if(CREATE_SHAKEMAP) then
    call wmo_create_shakemap_o()
  endif

  ! saves MOVIE in full 3D MESH
  if(MOVIE_VOLUME .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
    call wmo_movie_volume_output()
  endif

  ! creates cross-section GIF image
  if(PNM_GIF_IMAGE .and. mod(it,NTSTEP_BETWEEN_FRAMES) == 0 ) then
    call write_PNM_GIF_create_image()
  endif

  end subroutine write_movie_output
  
  
  
!================================================================
  
  subroutine wmo_create_shakemap_em()

! creation of shapemap file
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie
  implicit none
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: &
    displ_element,veloc_element,accel_element
  integer :: ipoin,ispec,iglob,ispec2D
  integer :: i,j,k
  logical :: is_done

! initializes arrays for point coordinates
  if (it == 1) then
    store_val_ux_external_mesh(:) = -HUGEVAL
    store_val_uy_external_mesh(:) = -HUGEVAL
    store_val_uz_external_mesh(:) = -HUGEVAL
    do ispec2D = 1,nfaces_surface_ext_mesh
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, NGLLX*NGLLY
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = zstore(iglob)
        enddo
      else
        do ipoin = 1, 4
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = zstore(iglob)        
        enddo
      endif
    enddo
  endif

! stores displacement, velocity and acceleration amplitudes
  do ispec2D = 1,nfaces_surface_ext_mesh
    ispec = faces_surface_ext_mesh_ispec(ispec2D)    
    
    if( ispec_is_acoustic(ispec) ) then
      ! displacement vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_acoustic, displ_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
      ! velocity vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, veloc_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
      ! accel ?
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_dot_acoustic, accel_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)                          
    endif

    
    ! high-resolution
    if (USE_HIGHRES_FOR_MOVIES) then
      do ipoin = 1, NGLLX*NGLLY
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! saves norm of displacement,velocity and acceleration vector
        if( ispec_is_elastic(ispec) ) then            
          ! norm of displacement
          store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
               max(store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
               sqrt(displ(1,iglob)**2 + displ(2,iglob)**2 + displ(3,iglob)**2))
          ! norm of velocity     
          store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
               max(store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
               sqrt(veloc(1,iglob)**2 + veloc(2,iglob)**2 + veloc(3,iglob)**2))
          ! norm of acceleration     
          store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
               max(store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
               sqrt(accel(1,iglob)**2 + accel(2,iglob)**2 + accel(3,iglob)**2))
        endif
        
        ! acoustic domains
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  ! norm of displacement
                  store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
                    max(store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
                        sqrt(displ_element(1,i,j,k)**2 &
                            + displ_element(2,i,j,k)**2 &
                            + displ_element(3,i,j,k)**2))
                  ! norm of velocity     
                  store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
                    max(store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
                        sqrt(veloc_element(1,i,j,k)**2 &
                            + veloc_element(2,i,j,k)**2 &
                            + veloc_element(3,i,j,k)**2))
                  ! norm of acceleration     
                  store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = &
                    max(store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin), &
                        sqrt(accel_element(1,i,j,k)**2 &
                            + accel_element(2,i,j,k)**2 &
                            + accel_element(3,i,j,k)**2))
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
        
      enddo
    else
      ! low-resolution: only corner points outputted
      do ipoin = 1, 4
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! saves norm of displacement,velocity and acceleration vector
        if( ispec_is_elastic(ispec) ) then                    
          ! norm of displacement
          store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                max(store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                sqrt(displ(1,iglob)**2 + displ(2,iglob)**2 + displ(3,iglob)**2))
          ! norm of velocity      
          store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                max(store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                sqrt(veloc(1,iglob)**2 + veloc(2,iglob)**2 + veloc(3,iglob)**2))
          ! norm of acceleration
          store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                max(store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                sqrt(accel(1,iglob)**2 + accel(2,iglob)**2 + accel(3,iglob)**2))
        endif
        
        ! acoustic domains
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  ! norm of displacement
                  store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                    max(store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                        sqrt(displ_element(1,i,j,k)**2 &
                            + displ_element(2,i,j,k)**2 &
                            + displ_element(3,i,j,k)**2))
                  ! norm of velocity     
                  store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                    max(store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                        sqrt(veloc_element(1,i,j,k)**2 &
                            + veloc_element(2,i,j,k)**2 &
                            + veloc_element(3,i,j,k)**2))
                  ! norm of acceleration     
                  store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = &
                    max(store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin), &
                        sqrt(accel_element(1,i,j,k)**2 &
                            + accel_element(2,i,j,k)**2 &
                            + accel_element(3,i,j,k)**2))
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
      enddo
    endif
  enddo

! finalizes shakemap: master process collects all info   
  if (it == NSTEP) then
    if (USE_HIGHRES_FOR_MOVIES) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    else
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif

! creates shakemap file
    if(myrank == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata',status='unknown',form='unformatted')
      write(IOUT) store_val_x_all_external_mesh   ! x coordinates
      write(IOUT) store_val_y_all_external_mesh   ! y coordinates
      write(IOUT) store_val_z_all_external_mesh   ! z coordinates
      write(IOUT) store_val_ux_all_external_mesh  ! norm of displacement vector
      write(IOUT) store_val_uy_all_external_mesh  ! norm of velocity vector
      write(IOUT) store_val_uz_all_external_mesh  ! norm of acceleration vector
      close(IOUT)
    endif
  endif
  
  end subroutine wmo_create_shakemap_em
  
  
!================================================================

  subroutine wmo_create_movie_surface_em()

! creation of moviedata files  

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie  
  implicit none
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: veloc_element
  integer :: ispec2D,ispec,ipoin,iglob,i,j,k
  logical :: is_done
  
! initializes arrays for point coordinates
  if (it == NTSTEP_BETWEEN_FRAMES ) then
    do ispec2D = 1,nfaces_surface_ext_mesh
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, NGLLX*NGLLY
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = zstore(iglob)
        enddo
      else
        do ipoin = 1, 4
          iglob = faces_surface_ext_mesh(ipoin,ispec2D)
          ! x,y,z coordinates
          store_val_x_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = xstore(iglob)
          store_val_y_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = ystore(iglob)
          store_val_z_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = zstore(iglob)                  
        enddo
      endif
    enddo
  endif
  
! saves surface velocities
  do ispec2D = 1,nfaces_surface_ext_mesh
    ispec = faces_surface_ext_mesh_ispec(ispec2D)      

    if( ispec_is_acoustic(ispec) ) then
      ! velocity vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, veloc_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
    endif
    
    if (USE_HIGHRES_FOR_MOVIES) then
      do ipoin = 1, NGLLX*NGLLY
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! saves velocity vector        
        if( ispec_is_elastic(ispec) ) then
          ! velocity x,y,z-components
          store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc(1,iglob)
          store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc(2,iglob)
          store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc(3,iglob)
        endif
        
        ! acoustic pressure potential
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  store_val_ux_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc_element(1,i,j,k)
                  store_val_uy_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc_element(2,i,j,k)
                  store_val_uz_external_mesh(NGLLX*NGLLY*(ispec2D-1)+ipoin) = veloc_element(3,i,j,k)
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
      enddo
    else
      do ipoin = 1, 4
        iglob = faces_surface_ext_mesh(ipoin,ispec2D)
        ! saves velocity vector        
        if( ispec_is_elastic(ispec) ) then
          ! velocity x,y,z-components
          store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc(1,iglob)
          store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc(2,iglob)
          store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc(3,iglob)      
        endif
        
        ! acoustic pressure potential
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  store_val_ux_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc_element(1,i,j,k)
                  store_val_uy_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc_element(2,i,j,k)
                  store_val_uz_external_mesh(NGNOD2D*(ispec2D-1)+ipoin) = veloc_element(3,i,j,k)
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
      enddo
    endif
  enddo

! master process collects all info
  if (USE_HIGHRES_FOR_MOVIES) then
    ! collects locations only once
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    endif
    ! updates/gathers velocity field (high-res)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
  else
    ! collects locations only once
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif
    ! updates/gathers velocity field (low-res)
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
  endif

! file output
  if(myrank == 0) then
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
    write(IOUT) store_val_x_all_external_mesh   ! x coordinate
    write(IOUT) store_val_y_all_external_mesh   ! y coordinate
    write(IOUT) store_val_z_all_external_mesh   ! z coordinate
    write(IOUT) store_val_ux_all_external_mesh  ! velocity x-component
    write(IOUT) store_val_uy_all_external_mesh  ! velocity y-component
    write(IOUT) store_val_uz_all_external_mesh  ! velocity z-component
    close(IOUT)
  endif
  
  end subroutine wmo_create_movie_surface_em

    
!=====================================================================

  subroutine wmo_movie_surface_output_o()

! outputs moviedata files  
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie  
  implicit none
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: val_element
  integer :: ispec,ipoin,iglob,i,j,k
  integer :: imin,imax,jmin,jmax,kmin,kmax,iface,igll,iloc
  logical :: is_done

  ! initializes arrays for point coordinates
  if (it == NTSTEP_BETWEEN_FRAMES ) then
    ipoin = 0
    do iface=1,num_free_surface_faces
      ispec = free_surface_ispec(iface)
      ! high_resolution
      if (USE_HIGHRES_FOR_MOVIES) then      
        do igll = 1, NGLLSQUARE
          ipoin = ipoin + 1
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)      
          iglob = ibool(i,j,k,ispec)
          ! coordinates
          store_val_x_external_mesh(ipoin) = xstore(iglob)
          store_val_y_external_mesh(ipoin) = ystore(iglob)
          store_val_z_external_mesh(ipoin) = zstore(iglob)
        enddo
      else
        imin = minval( free_surface_ijk(1,:,iface) )
        imax = maxval( free_surface_ijk(1,:,iface) )
        jmin = minval( free_surface_ijk(2,:,iface) )
        jmax = maxval( free_surface_ijk(2,:,iface) )
        kmin = minval( free_surface_ijk(3,:,iface) )
        kmax = maxval( free_surface_ijk(3,:,iface) )      
        do iloc = 1, NGNOD2D    
          ipoin = ipoin + 1
          ! corner points
          if( imin == imax ) then
            iglob = ibool(imin,iorderi(iloc),iorderj(iloc),ispec)
          else if( jmin == jmax ) then
            iglob = ibool(iorderi(iloc),jmin,iorderj(iloc),ispec)
          else
            iglob = ibool(iorderi(iloc),iorderj(iloc),kmin,ispec)
          endif
          ! coordinates
          store_val_x_external_mesh(ipoin) = xstore(iglob)
          store_val_y_external_mesh(ipoin) = ystore(iglob)
          store_val_z_external_mesh(ipoin) = zstore(iglob)
        enddo
      endif
    enddo
  endif

  
  ! outputs values at free surface
  ipoin = 0
  do iface=1,num_free_surface_faces
    ispec = free_surface_ispec(iface)
    
    if( ispec_is_acoustic(ispec) ) then
      if(SAVE_DISPLACEMENT) then
        ! displacement vector
        call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_acoustic, val_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)      
      else
        ! velocity vector
        call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, val_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
      endif
    endif
    
    
    ! high_resolution
    if (USE_HIGHRES_FOR_MOVIES) then      
      do igll = 1, NGLLSQUARE
        ipoin = ipoin + 1
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)      
        iglob = ibool(i,j,k,ispec)
        ! elastic displacement/velocity
        if( ispec_is_elastic(ispec) ) then
          if(SAVE_DISPLACEMENT) then
             store_val_ux_external_mesh(ipoin) = displ(1,iglob)
             store_val_uy_external_mesh(ipoin) = displ(2,iglob)
             store_val_uz_external_mesh(ipoin) = displ(3,iglob)
          else
             store_val_ux_external_mesh(ipoin) = veloc(1,iglob)
             store_val_uy_external_mesh(ipoin) = veloc(2,iglob)
             store_val_uz_external_mesh(ipoin) = veloc(3,iglob)
          endif
        endif

        ! acoustic pressure potential
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  store_val_ux_external_mesh(ipoin) = val_element(1,i,j,k)
                  store_val_uy_external_mesh(ipoin) = val_element(2,i,j,k)
                  store_val_uz_external_mesh(ipoin) = val_element(3,i,j,k)
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
        
      enddo
    else    
      imin = minval( free_surface_ijk(1,:,iface) )
      imax = maxval( free_surface_ijk(1,:,iface) )
      jmin = minval( free_surface_ijk(2,:,iface) )
      jmax = maxval( free_surface_ijk(2,:,iface) )
      kmin = minval( free_surface_ijk(3,:,iface) )
      kmax = maxval( free_surface_ijk(3,:,iface) )      
      do iloc = 1, NGNOD2D    
        ipoin = ipoin + 1
        ! corner points
        if( imin == imax ) then
          iglob = ibool(imin,iorderi(iloc),iorderj(iloc),ispec)
        else if( jmin == jmax ) then
          iglob = ibool(iorderi(iloc),jmin,iorderj(iloc),ispec)
        else
          iglob = ibool(iorderi(iloc),iorderj(iloc),kmin,ispec)
        endif

        ! elastic displacement/velocity
        if( ispec_is_elastic(ispec) ) then
          if(SAVE_DISPLACEMENT) then
             store_val_ux_external_mesh(ipoin) = displ(1,iglob)
             store_val_uy_external_mesh(ipoin) = displ(2,iglob)
             store_val_uz_external_mesh(ipoin) = displ(3,iglob)
          else
             store_val_ux_external_mesh(ipoin) = veloc(1,iglob)
             store_val_uy_external_mesh(ipoin) = veloc(2,iglob)
             store_val_uz_external_mesh(ipoin) = veloc(3,iglob)
          endif
        endif
        
        ! acoustic pressure potential
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then
                  store_val_ux_external_mesh(ipoin) = val_element(1,i,j,k)
                  store_val_uy_external_mesh(ipoin) = val_element(2,i,j,k)
                  store_val_uz_external_mesh(ipoin) = val_element(3,i,j,k)
                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif        
        
      enddo ! iloc
    endif
  enddo ! iface

! master process collects all info
  if (USE_HIGHRES_FOR_MOVIES) then
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    endif
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
  else
    if (it == NTSTEP_BETWEEN_FRAMES ) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif
    call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
         store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
         nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
  endif

! file output: note that values are only stored on free surface
  if(myrank == 0) then
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted')
    write(IOUT) store_val_x_all_external_mesh   ! x coordinate
    write(IOUT) store_val_y_all_external_mesh   ! y coordinate
    write(IOUT) store_val_z_all_external_mesh   ! z coordinate
    write(IOUT) store_val_ux_all_external_mesh  ! velocity x-component
    write(IOUT) store_val_uy_all_external_mesh  ! velocity y-component
    write(IOUT) store_val_uz_all_external_mesh  ! velocity z-component
    close(IOUT)
  endif

  end subroutine wmo_movie_surface_output_o
  
  
!=====================================================================

  subroutine wmo_create_shakemap_o()

! outputs shakemap file 
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie
  
  implicit none
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: &
    displ_element,veloc_element,accel_element
  integer :: ipoin,ispec,iglob
  integer :: imin,imax,jmin,jmax,kmin,kmax,iface,igll,iloc
  integer :: i,j,k
  logical :: is_done  

  ! outputs values on free surface  
  ipoin = 0
  do iface=1,num_free_surface_faces
    ispec = free_surface_ispec(iface)
    
    if( ispec_is_acoustic(ispec) ) then
      ! displacement vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_acoustic, displ_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
      ! velocity vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, veloc_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)
      ! accel ?
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_dot_acoustic, accel_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore)                          
    endif
    
    
    ! save all points for high resolution, or only four corners for low resolution
    if(USE_HIGHRES_FOR_MOVIES) then
      do igll = 1, NGLLSQUARE
        ipoin = ipoin + 1
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        store_val_x_external_mesh(ipoin) = xstore(iglob)
        store_val_y_external_mesh(ipoin) = ystore(iglob)
        store_val_z_external_mesh(ipoin) = zstore(iglob)
        ! todo: are we only interested in the absolute maximum of horizontal (E,N) components?
        if( ispec_is_elastic( ispec) ) then
          ! horizontal displacement
          store_val_ux_external_mesh(ipoin) = max(store_val_ux_external_mesh(ipoin),&
                                                abs(displ(1,iglob)),abs(displ(2,iglob)))
          ! horizontal velocity
          store_val_uy_external_mesh(ipoin) = max(store_val_uy_external_mesh(ipoin),&
                                                abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          ! horizontal acceleration
          store_val_uz_external_mesh(ipoin) = max(store_val_uz_external_mesh(ipoin),&
                                                abs(accel(1,iglob)),abs(accel(2,iglob)))
        endif
        
        ! acoustic domains
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then        
                  ! horizontal displacement
                  store_val_ux_external_mesh(ipoin) = max(store_val_ux_external_mesh(ipoin),&
                                                abs(displ_element(1,i,j,k)),abs(displ_element(2,i,j,k)))
                  ! horizontal velocity
                  store_val_uy_external_mesh(ipoin) = max(store_val_uy_external_mesh(ipoin),&
                                                abs(veloc_element(1,i,j,k)),abs(veloc_element(2,i,j,k)))
                  ! horizontal acceleration
                  store_val_uz_external_mesh(ipoin) = max(store_val_uz_external_mesh(ipoin),&
                                                abs(accel_element(1,i,j,k)),abs(accel_element(2,i,j,k)))

                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
        
      enddo    
    else
      imin = minval( free_surface_ijk(1,:,iface) )
      imax = maxval( free_surface_ijk(1,:,iface) )
      jmin = minval( free_surface_ijk(2,:,iface) )
      jmax = maxval( free_surface_ijk(2,:,iface) )
      kmin = minval( free_surface_ijk(3,:,iface) )
      kmax = maxval( free_surface_ijk(3,:,iface) )
      do iloc = 1, NGNOD2D
        ipoin = ipoin + 1
        ! corner points
        if( imin == imax ) then
          iglob = ibool(imin,iorderi(iloc),iorderj(iloc),ispec)
        else if( jmin == jmax ) then
          iglob = ibool(iorderi(iloc),jmin,iorderj(iloc),ispec)
        else
          iglob = ibool(iorderi(iloc),iorderj(iloc),kmin,ispec)
        endif        
        ! coordinates
        store_val_x_external_mesh(ipoin) = xstore(iglob)
        store_val_y_external_mesh(ipoin) = ystore(iglob)
        store_val_z_external_mesh(ipoin) = zstore(iglob)
        ! todo: are we only interested in the absolute maximum of horizontal (E,N) components?
        if( ispec_is_elastic( ispec) ) then        
          store_val_ux_external_mesh(ipoin) = max(store_val_ux_external_mesh(ipoin),&
                                                  abs(displ(1,iglob)),abs(displ(2,iglob)))
          store_val_uy_external_mesh(ipoin) = max(store_val_uy_external_mesh(ipoin),&
                                                  abs(veloc(1,iglob)),abs(veloc(2,iglob)))
          store_val_uz_external_mesh(ipoin) = max(store_val_uz_external_mesh(ipoin),&
                                                  abs(accel(1,iglob)),abs(accel(2,iglob)))
        endif

        ! acoustic domains
        if( ispec_is_acoustic(ispec) ) then
          ! velocity vector
          is_done = .false.
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                if( iglob == ibool(i,j,k,ispec) ) then        
                  ! horizontal displacement
                  store_val_ux_external_mesh(ipoin) = max(store_val_ux_external_mesh(ipoin),&
                                                abs(displ_element(1,i,j,k)),abs(displ_element(2,i,j,k)))
                  ! horizontal velocity
                  store_val_uy_external_mesh(ipoin) = max(store_val_uy_external_mesh(ipoin),&
                                                abs(veloc_element(1,i,j,k)),abs(veloc_element(2,i,j,k)))
                  ! horizontal acceleration
                  store_val_uz_external_mesh(ipoin) = max(store_val_uz_external_mesh(ipoin),&
                                                abs(accel_element(1,i,j,k)),abs(accel_element(2,i,j,k)))

                  is_done = .true.
                  exit                  
                endif
              enddo
              if( is_done ) exit
            enddo
            if( is_done ) exit
          enddo
        endif
        
      enddo
    endif ! USE_HIGHRES_FOR_MOVIES
  enddo

  ! saves shakemap only at the end of the simulation
  if(it == NSTEP) then
    if (USE_HIGHRES_FOR_MOVIES) then
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGLLX*NGLLY,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGLLX*NGLLY,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGLLX*NGLLY,NPROC)
    else
      call gatherv_all_cr(store_val_x_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_x_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_y_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_y_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_z_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_z_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_ux_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_ux_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uy_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uy_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
      call gatherv_all_cr(store_val_uz_external_mesh,nfaces_surface_ext_mesh*NGNOD2D,&
           store_val_uz_all_external_mesh,nfaces_perproc_surface_ext_mesh*NGNOD2D,faces_surface_offset_ext_mesh,&
           nfaces_surface_glob_ext_mesh*NGNOD2D,NPROC)
    endif

    ! creates shakemap file: note that values are only stored on free surface
    if(myrank == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/shakingdata',status='unknown',form='unformatted')
      write(IOUT) store_val_x_all_external_mesh   ! x coordinates
      write(IOUT) store_val_y_all_external_mesh   ! y coordinates
      write(IOUT) store_val_z_all_external_mesh   ! z coordinates
      write(IOUT) store_val_ux_all_external_mesh  ! norm of displacement vector
      write(IOUT) store_val_uy_all_external_mesh  ! norm of velocity vector
      write(IOUT) store_val_uz_all_external_mesh  ! norm of acceleration vector
      close(IOUT)
    endif

  endif ! NTSTEP

  end subroutine wmo_create_shakemap_o

    
!=====================================================================

  subroutine wmo_movie_volume_output()

! outputs movie files for div, curl and velocity  
  
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_movie
  implicit none
  
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: veloc_element
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB):: div_glob,curl_glob ! divergence and curl only in the global nodes
  integer :: ispec,i,j,k,l,iglob
  integer,dimension(NGLOB_AB) :: valency
  
  ! saves velocity here to avoid static offset on displacement for movies
  velocity_x(:,:,:,:) = 0._CUSTOM_REAL
  velocity_y(:,:,:,:) = 0._CUSTOM_REAL
  velocity_z(:,:,:,:) = 0._CUSTOM_REAL
  
  if( ACOUSTIC_SIMULATION ) then
    ! uses div as temporary array to store velocity on all gll points
    do ispec=1,NSPEC_AB
      if( .not. ispec_is_acoustic(ispec) ) cycle

      ! calculates velocity
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                        potential_dot_acoustic, veloc_element,&
                        hprime_xx,hprime_yy,hprime_zz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        ibool,rhostore)
      velocity_x(:,:,:,ispec) = veloc_element(1,:,:,:)
      velocity_y(:,:,:,ispec) = veloc_element(2,:,:,:)
      velocity_z(:,:,:,ispec) = veloc_element(3,:,:,:)      
    enddo
  endif ! acoustic

  ! saves full snapshot data to local disk
  if( ELASTIC_SIMULATION ) then
    div_glob=0.0_CUSTOM_REAL
    curl_glob=0.0_CUSTOM_REAL

    do ispec=1,NSPEC_AB
      if( .not. ispec_is_elastic(ispec) ) cycle
      
      ! calculates divergence and curl of velocity field
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempx3l = 0._CUSTOM_REAL
            tempy1l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL

            do l=1,NGLLX
              hp1 = hprime_xx(i,l)
              iglob = ibool(l,j,k,ispec)
              tempx1l = tempx1l + veloc(1,iglob)*hp1
              tempy1l = tempy1l + veloc(2,iglob)*hp1
              tempz1l = tempz1l + veloc(3,iglob)*hp1
              hp2 = hprime_yy(j,l)
              iglob = ibool(i,l,k,ispec)
              tempx2l = tempx2l + veloc(1,iglob)*hp2
              tempy2l = tempy2l + veloc(2,iglob)*hp2
              tempz2l = tempz2l + veloc(3,iglob)*hp2
              hp3 = hprime_zz(k,l)
              iglob = ibool(i,j,l,ispec)
              tempx3l = tempx3l + veloc(1,iglob)*hp3
              tempy3l = tempy3l + veloc(2,iglob)*hp3
              tempz3l = tempz3l + veloc(3,iglob)*hp3
            enddo

            ! get derivatives of ux, uy and uz with respect to x, y and z
            xixl = xix(i,j,k,ispec)
            xiyl = xiy(i,j,k,ispec)
            xizl = xiz(i,j,k,ispec)
            etaxl = etax(i,j,k,ispec)
            etayl = etay(i,j,k,ispec)
            etazl = etaz(i,j,k,ispec)
            gammaxl = gammax(i,j,k,ispec)
            gammayl = gammay(i,j,k,ispec)
            gammazl = gammaz(i,j,k,ispec)

            dvxdxl(i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
            dvxdyl(i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
            dvxdzl(i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

            dvydxl(i,j,k) = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
            dvydyl(i,j,k) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
            dvydzl(i,j,k) = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

            dvzdxl(i,j,k) = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
            dvzdyl(i,j,k) = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
            dvzdzl(i,j,k) = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

          enddo
        enddo
      enddo

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ! divergence \nabla \cdot \bf{v}
            div(i,j,k,ispec) = dvxdxl(i,j,k) + dvydyl(i,j,k) + dvzdzl(i,j,k)
            ! curl 
            curl_x(i,j,k,ispec) = dvzdyl(i,j,k) - dvydzl(i,j,k)
            curl_y(i,j,k,ispec) = dvxdzl(i,j,k) - dvzdxl(i,j,k)
            curl_z(i,j,k,ispec) = dvydxl(i,j,k) - dvxdyl(i,j,k)
            ! velocity field
            iglob = ibool(i,j,k,ispec)
            velocity_x(i,j,k,ispec) = veloc(1,iglob)
            velocity_y(i,j,k,ispec) = veloc(2,iglob)
            velocity_z(i,j,k,ispec) = veloc(3,iglob)

            valency(iglob)=valency(iglob)+1
            
            div_glob(iglob) = div_glob(iglob) + div(i,j,k,ispec)
            curl_glob(iglob)=curl_glob(iglob)+0.5_CUSTOM_REAL*(curl_x(i,j,k,ispec)+curl_x(i,j,k,ispec)+curl_x(i,j,k,ispec))
          enddo
        enddo
      enddo
    enddo !NSPEC_AB

    do i=1,NGLOB_AB
      div_glob(i)=div_glob(i)/valency(i)
      curl_glob(i)=curl_glob(i)/valency(i)
    enddo
    
    write(outputname,"('/proc',i6.6,'_div_glob_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) div_glob
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_glob_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_glob
    close(27)

    write(outputname,"('/proc',i6.6,'_div_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) div
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_x_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_x
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_y_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_y
    close(27)
    write(outputname,"('/proc',i6.6,'_curl_z_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) curl_z
    close(27)
    
    !write(outputname,"('veloc_proc',i6.6,'_it',i6.6,'.bin')") myrank,it
    !open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    !write(27) veloc
    !close(27)
  
  endif ! elastic
 
  if( ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION ) then
    write(outputname,"('/proc',i6.6,'_velocity_N_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) velocity_x
    close(27)  

    write(outputname,"('/proc',i6.6,'_velocity_E_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) velocity_y
    close(27)  

    write(outputname,"('/proc',i6.6,'_velocity_Z_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) velocity_z
    close(27)  

    !write(outputname,"('/proc',i6.6,'_veloc_it',i6.6,'.bin')") myrank,it
    !open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    !write(27) velocity_movie
    !close(27)  

  endif 
  
  end subroutine wmo_movie_volume_output
    
