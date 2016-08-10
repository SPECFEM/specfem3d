module wave2d_solver

  use wave2d_constants
  use wave2d_variables
  use wave2d_define_der_matrices

  implicit none

contains

  subroutine mesher

    integer ispec,ib,i,j,k,iglob,iglob1,itime,ix,iz

    ! set up grid and spectral elements
    call define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,wgllwgll_xz)

    ispec = 0
    iglob = 0
    nspecb(:) = 0

    ! loop over all elements
    do iz = 1,NEZ
       do ix = 1,NEX
          ispec = ispec+1

          ! evenly spaced anchors between 0 and 1
          x1(ispec) = LENGTH*dble(ix-1)/dble(NEX)
          x2(ispec) = LENGTH*dble(ix)/dble(NEX)
          z1(ispec) = HEIGHT*dble(iz-1)/dble(NEZ)
          z2(ispec) = HEIGHT*dble(iz)/dble(NEZ)

          ! loop over GLL points to calculate jacobian, and set up numbering
          !
          ! jacobian = | dx/dxi dx/dgamma | = (z2-z1)*(x2-x1)/4  as dx/dgamma=dz/dxi = 0
          !            | dz/dxi dz/dgamma |
          !
          do j = 1,NGLLZ
             do i = 1,NGLLX

                ! jacobian, integration weight
                dxidx(i,j,ispec) = 2. / (x2(ispec)-x1(ispec))
                dxidz(i,j,ispec) = 0.
                dgammadx(i,j,ispec) = 0.
                dgammadz(i,j,ispec) = 2. / (z2(ispec)-z1(ispec))
                jacobian(i,j,ispec) = (z2(ispec)-z1(ispec))*(x2(ispec)-x1(ispec)) / 4.

                ! set up local to global numbering
                if ( (i == 1) .and. (ix > 1) ) then
                   ibool(i,j,ispec) = ibool(NGLLX,j,ispec-1)
                else if ( (j == 1) .and. (iz > 1) ) then
                   ibool(i,j,ispec) = ibool(i,NGLLZ,ispec-NEX)
                else
                   iglob = iglob + 1
                   ibool(i,j,ispec) = iglob
                endif

                ! get the global gridpoints
                iglob1 = ibool(i,j,ispec)
                x(iglob1) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
                z(iglob1) = 0.5*(1.-zigll(j))*z1(ispec)+0.5*(1.+zigll(j))*z2(ispec)

                ! end loop over GLL points
             enddo
          enddo

          ! if boundary element
          ! 1,2,3,4 --> left, right, bottom, top
          if (ix == 1) then      ! left boundary
             nspecb(1) = nspecb(1) + 1
             ibelm(1,nspecb(1)) = ispec
             do j = 1,NGLLZ
                jacobianb(1,j,nspecb(1))= (z2(ispec)-z1(ispec))/2.
             enddo
          endif
          if (ix == NEX) then    ! right boundary
             nspecb(2) = nspecb(2) + 1
             ibelm(2,nspecb(2)) = ispec
             do j = 1,NGLLZ
                jacobianb(2,j,nspecb(2))= (z2(ispec)-z1(ispec))/2.
             enddo
          endif
          if (iz == 1) then      ! bottom boundary
             nspecb(3) = nspecb(3) + 1
             ibelm(3,nspecb(3)) = ispec
             do i = 1,NGLLX
                jacobianb(3,i,nspecb(3))= (x2(ispec)-x1(ispec))/2.
             enddo
          endif
          if (iz == NEZ) then    ! top boundary
             nspecb(4) = nspecb(4) + 1
             ibelm(4,nspecb(4)) = ispec
             do i = 1,NGLLX
                jacobianb(4,i,nspecb(4))= (x2(ispec)-x1(ispec))/2.
             enddo
          endif
          ! end loop over elements
       enddo
    enddo

!!$    do ibb=1,4
!!$       print *, ibb, nspecb(ibb)
!!$    enddo
!!$    print *
!!$    print *
!!$    do ibb=1,4
!!$       do ispec=1,nspecb(ibb)
!!$          print *, ibb, ispec, ibelm(ibb,ispec)
!!$       enddo
!!$    enddo
!!$    stop 'testing the boundary indexing'

    ! estimate the time step
    dh = HEIGHT/dble((NGLLZ-1)*NEZ)
    if (dh > LENGTH/dble((NGLLX-1)*NEX)) dh = LENGTH/dble((NGLLX-1)*NEX)
    c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
    time_step = 0.2*dh/c
    print *
    print *,'space step (km):', sngl(dh/1000.)
    print *,'time step estimate from courant = 0.2: ',sngl(time_step),' seconds'
    print *,'  actual time step: ',sngl(DT),' seconds'

  end subroutine mesher

  !-------------------------------------------------------

  subroutine set_model_property()

    integer :: ispec, i, j, iglob

    ! model properties -- globally defined (no discontinuities permitted this way)
    do iglob = 1,NGLOB
       rho_global(iglob)   = DENSITY
       kappa_global(iglob) = INCOMPRESSIBILITY
       if (ISURFACE == 0) then
          mu_global(iglob) = RIGIDITY
       else
          ! KEY: this means that the S velocity will be the surface wave phase velocity (m/s)
          mu_global(iglob) = DENSITY*(c_glob(iglob))**2
          !if (ihomo==1) mu_global(iglob) = DENSITY*c0**2
          !if (ihomo==0) mu_global(iglob) = DENSITY*(c_glob(iglob))**2
       endif
    enddo

    ! fill local arrays from global vectors
    do ispec = 1,NSPEC
       do j = 1,NGLLZ
          do i = 1,NGLLX
             iglob = ibool(i,j,ispec)
             rho(i,j,ispec)   = rho_global(iglob)
             kappa(i,j,ispec) = kappa_global(iglob)
             mu(i,j,ispec)    = mu_global(iglob)
          enddo
       enddo
    enddo

!!$    ! properties
!!$    do ispec = 1,NSPEC
!!$       !  get center of element
!!$       iglob = ibool(NGLLX/2,NGLLZ/2,ispec)
!!$       do j = 1,NGLLZ
!!$          do i = 1,NGLLX
!!$             if (z(iglob) >= 0) then
!!$                ! crust
!!$                rho(i,j,ispec) = DENSITY
!!$                kappa(i,j,ispec) = INCOMPRESSIBILITY
!!$                mu(i,j,ispec) = RIGIDITY
!!$             else
!!$                ! mantle
!!$                rho(i,j,ispec) = 3380.
!!$                kappa(i,j,ispec) = 1.3d+11
!!$                mu(i,j,ispec) = 6.8d+10
!!$             endif
!!$          enddo
!!$       enddo
!!$    enddo

!!$    ! create vectors from the arrays
!!$    do ispec = 1,NSPEC
!!$       do j = 1,NGLLZ
!!$          do i = 1,NGLLX
!!$             iglob = ibool(i,j,ispec)
!!$             rho_global(iglob) = rho(i,j,ispec)
!!$             mu_global(iglob) = mu(i,j,ispec)
!!$             kappa_global(iglob) = kappa(i,j,ispec)
!!$          enddo
!!$       enddo
!!$    enddo

  end subroutine set_model_property

!---------------------------------------------

  subroutine solver(solver_type, nsrc, sglob, ispec_src, hxis_store, hgammas_store, samp, &
                                 nrec, rglob, ispec_rec, hxir_store, hgammar_store, ramp, &
              last_frame, absorbfield, &
              rhop_kernel, beta_kernel, alpha_kernel, &
              three_source_model, stf, f0)

    integer, intent(in) :: solver_type
    integer, intent(in) :: nsrc, sglob(nsrc), ispec_src(nsrc)
    integer, intent(in) :: nrec, rglob(nrec), ispec_rec(nrec)
    double precision, intent(in) :: hxir_store(nrec,NGLLX), hgammar_store(nrec,NGLLZ)
    double precision, intent(in) :: hxis_store(nsrc,NGLLX), hgammas_store(nsrc,NGLLZ)
    double precision, intent(inout) :: samp(NSTEP,NCOMP,nsrc)
    double precision, intent(inout) :: ramp(NSTEP,NCOMP,nrec)

    ! if adjoint seismograms are desired while computing kernels, then we need to make sure
    ! not to overwrite the forward source time function
    double precision,  dimension(:,:,:), allocatable :: stf_for

    ! OPTIONAL ARGUMENTS
    character(len=*), optional :: last_frame
    double precision, intent(inout), optional :: absorbfield(NSTEP, NCOMP, NGLL, NELE, NABSORB)
    double precision, intent(out), optional :: three_source_model(NSTEP,NCOMP,nsrc,3)   ! 3 : deltax, deltaz, deltat
    double precision, intent(out), dimension(NGLOB), optional :: rhop_kernel, beta_kernel, alpha_kernel
    double precision, intent(in), optional :: stf(NSTEP), f0(NCOMP)

    !double precision, dimension(NGLOB),optional :: rho_kernel, mu_kernel, kappa_kernel
    double precision, dimension(NGLOB) :: rho_kernel, mu_kernel, kappa_kernel, beta_kernel_prime
    double precision, dimension(NSTEP) :: stf_d

    ! displacement gradients: only for SH motion (NCOMP = 1), which gives 2 gradients
    ! only save the GLL points on the element where the source is located
    double precision :: displ_grad(NGLLX,NGLLZ,2)

    ! CHT
    integer tlab
    double precision xtemp, ztemp, temp1, temp2, temp3

    ! Lagrange interpolation
    double precision :: hlagrange

    integer ispec,ib,i,j,k,iglob,iglob1,iglob2,itime,ix,iz,itime1,itime2,itime0
    integer isrc, irec, icomp, isave
    character(len=100) :: filename,filename1,filename2,filename3,filename4,filename5,filename6,fm
    double precision, dimension(NGLOB) :: mu_k, kappa_k
    logical :: save_forward

    !--------------------------------------

    if (NCOMP == 3) fm = '(9e12.3)'
    if (NCOMP == 1) fm = '(3e12.3)'

    ! test of input arguments
    if (solver_type /= 1 .and. solver_type /= 2 .and. solver_type /= 3) then
       stop 'solver_type has to be 1, 2 or 3'
    endif

    save_forward = .false.
    if (solver_type == 1) then
       if (present(last_frame) .and. present(absorbfield)) save_forward = .true.
    endif

    if (solver_type == 3) then
       if (.not. (present(last_frame) .and. present(absorbfield)  &
            .and. present(rhop_kernel) .and. present(beta_kernel) .and. present(alpha_kernel))) &
            stop 'For kernel calculation, last_frame, absorbfield and all kernel are in the argument'
       rho_kernel = 0.
       mu_kernel = 0.
       kappa_kernel = 0.

       ! allocate a new array for the forward source function (note: includes amplitude)
       allocate(stf_for(NSTEP,NCOMP,nsrc))
       stf_for(:,:,:) = samp(:,:,:)

       ! compute the derivative of the source time function
       stf_d(:) = 0.
       do i = 2, NSTEP-1
         stf_d(i) =  (stf(i+1) - stf(i-1)) / (2.*DT)
       enddo
       stf_d(1) = (stf(2) - stf(1)) / DT
       stf_d(NSTEP) = (stf(NSTEP) - stf(NSTEP-1)) / DT

       ! initialize adjoint seismograms
       samp(:,:,:) = 0.

    endif

    ! gridpoints per wavelength estimation
    print *
    print *, 'space step (km):', sngl(dh/1000.)
    if (ISURFACE) then
       print *, 'wavelength-min (km):', sngl(2*hdur*cmin/1000.)
       print *, 'wavelength-max (km):', sngl(2*hdur*cmax/1000.)
       print *, 'number of gridpoints per wavelength for S:'
       print *, '  min (cmin) : ', floor(2*hdur*cmin/dh)
       print *, '  max (cmax) : ', floor(2*hdur*cmax/dh)

    else
       c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
       print *, 'number of gridpoints per wavelength for P: ', floor(hdur*c/dh)
       c = sqrt(RIGIDITY/DENSITY)
       print *, 'number of gridpoints per wavelength for S:', floor(hdur*c/dh)
    endif

    NINT = NSTEP/NSAVE
    if (NINT * NSAVE > NSTEP) stop 'NSTEP should equal to NINT * NSAVE'

    ! calculate the global mass matrix once and for all
    ! note that the density variation is included here
    mass_global(:) = 0.
    do ispec = 1,NSPEC
       do j = 1,NGLLZ
          do i = 1,NGLLX
             mass_local = wxgll(i)*wzgll(j)*rho(i,j,ispec)*jacobian(i,j,ispec)
             iglob = ibool(i,j,ispec)
             mass_global(iglob) = mass_global(iglob) + mass_local
          enddo
       enddo
    enddo

    ! time marching parameters
    deltat = DT
    deltatover2 = deltat/2.
    deltatsqover2 = deltat*deltat/2.

    if (solver_type == 3) then
       b_deltat = -DT
       b_deltatover2 = b_deltat/2.
       b_deltatsqover2 = b_deltat*b_deltat/2.
    endif

    ! initialize
    displ(:,:) = 0.
    veloc(:,:) = 0.
    accel(:,:) = 0.

    if (solver_type == 3) then
       open(11,file = trim(last_frame),status='old',iostat=ios)
       if (ios /= 0) stop 'Error reading the last frame'
       do i = 1, NGLOB
          !read(11,*) b_displ(1,i), b_displ(2,i), b_displ(3,i), &
          !           b_veloc(1,i), b_veloc(2,i), b_veloc(3,i), &
          !           b_accel(1,i), b_accel(2,i), b_accel(3,i)
          read(11,fm) (b_displ(j,i), j = 1,NCOMP), &
                       (b_veloc(j,i), j = 1,NCOMP), &
                       (b_accel(j,i), j = 1,NCOMP)
       enddo
       close(11)
    endif

    !
    ! START TIME MARCHING
    !
    do itime = 1, NSTEP

       if (mod(itime,200) == 0) write(*,*) itime, ' out of ', NSTEP

       ! 'predictor' update displacement using finite-difference time scheme (Newmark)
       do i = 1,NCOMP
          displ(i,:) = displ(i,:) + deltat*veloc(i,:) + deltatsqover2*accel(i,:)
          veloc(i,:) = veloc(i,:) + deltatover2*accel(i,:)
          accel(i,:) = 0.
          if (solver_type == 3) then
             b_displ(i,:) = b_displ(i,:) + b_deltat*b_veloc(i,:) + b_deltatsqover2*b_accel(i,:)
             b_veloc(i,:) = b_veloc(i,:) + b_deltatover2*b_accel(i,:)
             b_accel(i,:) = 0.
          endif
       enddo

 if (NCOMP == 1) then    ! SH, or surface waves only

       !
       !   INTEGRATION OVER SPECTRAL ELEMENTS
       !
       do ispec = 1,NSPEC

          ! first double loop over GLL
          ! compute and store gradients
          do j = 1,NGLLZ
             do i = 1,NGLLX

                iglob2 = ibool(i,j,ispec)

                ! derivative along x
                tempy1l = 0.
                if (solver_type == 3) b_tempy1l = 0.
                do k = 1,NGLLX
                   hp1 = hprime_xx(k,i)
                   iglob = ibool(k,j,ispec)
                   tempy1l = tempy1l + displ(1,iglob)*hp1
                   if (solver_type == 3) then
                      b_tempy1l = b_tempy1l + b_displ(1,iglob)*hp1
                   endif
                enddo

                ! derivative along z
                tempy2l = 0.
                if (solver_type == 3) b_tempy2l = 0.
                do k = 1,NGLLZ
                   hp2 = hprime_zz(k,j)
                   iglob = ibool(i,k,ispec)
                   tempy2l = tempy2l + displ(1,iglob)*hp2
                   if (solver_type == 3) then
                      b_tempy2l = b_tempy2l + b_displ(1,iglob)*hp2
                   endif
                enddo

                ! from mesher
                dxidxl = dxidx(i,j,ispec)
                dxidzl = dxidz(i,j,ispec)
                dgammadxl = dgammadx(i,j,ispec)
                dgammadzl = dgammadz(i,j,ispec)

                ! spatial gradients
                dsydxl = tempy1l*dxidxl + tempy2l*dgammadxl
                dsydzl = tempy1l*dxidzl + tempy2l*dgammadzl

                ! save spatial gradient for (point) source perturbations
                if (solver_type == 3 .and. ispec == ispec_src(1)) then
                   displ_grad(i,j,1) = dsydxl
                   displ_grad(i,j,2) = dsydzl
                endif

                ! strain tensor
                ds(:,:) = 0.
                ds(1,2) = oneovertwo * dsydxl
                ds(2,3) = oneovertwo * dsydzl
                ds(2,1) = ds(1,2)
                ds(3,2) = ds(2,3)

                if (solver_type == 3) then
                   b_dsydxl = b_tempy1l*dxidxl + b_tempy2l*dgammadxl
                   b_dsydzl = b_tempy1l*dxidzl + b_tempy2l*dgammadzl

                   b_ds(:,:) = 0.
                   b_ds(1,2) = oneovertwo * b_dsydxl
                   b_ds(2,3) = oneovertwo * b_dsydzl
                   b_ds(2,1) = b_ds(1,2)
                   b_ds(3,2) = b_ds(2,3)

                   ! mu kernel
                   !mu_k(iglob2) = sum(ds * b_ds) - ONE_THIRD * kappa_k(iglob2)   ! (12-July-2006)
                   mu_k(iglob2) = sum(ds * b_ds)
                endif

                mul = mu(i,j,ispec)      ! model heterogeneity
                sigma_xy = mul*dsydxl
                sigma_zy = mul*dsydzl

                if (solver_type == 3) then
                   b_sigma_xy = mul* b_dsydxl
                   b_sigma_zy = mul* b_dsydzl
                endif

                jacobianl = jacobian(i,j,ispec)

                ! non-symmetric form (first Piola-Kirchhoff stress)
                tempy1(i,j) = jacobianl*(sigma_xy*dxidxl+sigma_zy*dxidzl)
                tempy2(i,j) = jacobianl*(sigma_xy*dgammadxl+sigma_zy*dgammadzl)
                if (solver_type == 3) then
                   b_tempy1(i,j) = jacobianl*(b_sigma_xy*dxidxl+ b_sigma_zy*dxidzl)
                   b_tempy2(i,j) = jacobianl*(b_sigma_xy*dgammadxl+b_sigma_zy*dgammadzl)
                endif

             enddo
          enddo
          !
          ! second double-loop over GLL
          ! compute all rhs terms
          !
          do j = 1,NGLLZ
             do i = 1,NGLLX

                ! along x direction
                tempy1l = 0.
                if (solver_type == 3) b_tempy1l = 0.
                do k = 1,NGLLX
                   fac1 = wxgll(k)*hprime_xx(i,k)
                   tempy1l = tempy1l + tempy1(k,j)*fac1
                   if (solver_type == 3) then
                      b_tempy1l = b_tempy1l + b_tempy1(k,j)*fac1
                   endif
                enddo

                ! along z direction
                tempy2l = 0.
                if (solver_type == 3) b_tempy2l = 0.
                do k = 1,NGLLZ
                   fac2 = wzgll(k)*hprime_zz(j,k)
                   tempy2l = tempy2l + tempy2(i,k)*fac2
                   if (solver_type == 3) then
                      b_tempy2l = b_tempy2l + b_tempy2(i,k)*fac2
                   endif
                enddo

                fac1 = wzgll(j)
                fac2 = wxgll(i)

                iglob = ibool(i,j,ispec)
                accel(1,iglob) = accel(1,iglob) - (fac1*tempy1l + fac2*tempy2l)
                if (solver_type == 3) then
                   b_accel(1,iglob) = b_accel(1,iglob) - (fac1* b_tempy1l + fac2* b_tempy2l)
                endif

             enddo ! second loop over the GLL points
          enddo

       enddo ! end loop over all spectral elements

       !
       ! boundary conditions
       !
       ! forward propagation

       ! modifications by CHT are made to make the top boundary absorbing
       ! if this is used in other options (e.g., kernels), then those
       ! sections need to be adjusted as well

       do ibb = 1,NABSORB  ! index of grid boundary
          if (ibb == 1) then
             i = 1
          else if (ibb == 2) then
             i = NGLLX
          else if (ibb == 3) then
             i = 1
          else if (ibb == 4) then
             i = NGLLZ
          endif

          do ib = 1,nspecb(ibb)     ! elements on each grid boundary
             ispec = ibelm(ibb,ib)  ! (global) index of boundary element

             if (ibb == 1 .or. ibb == 2) then ! left or right boundary element
                j1 = 1; j2 = NGLLZ
             else if (ib == 1) then           ! top left corner element
                j1 = 2; j2 = NGLLX
             else if (ib == nspecb(ibb)) then ! top right corner element
                j1 = 1; j2 = NGLLX-1
             else                             ! top or bottom boundary (excluding corner elements)
                j1 = 1; j2 = NGLLX
             endif

             do j = j1, j2
                if (ibb == 1 .or. ibb == 2) then  ! left or right boundary
                   iglob = ibool(i,j,ispec)
                   rho_vs = dsqrt(rho(i,j,ispec)*mu(i,j,ispec))
                else                              ! top or bottom boundary
                   iglob = ibool(j,i,ispec)
                   rho_vs = dsqrt(rho(j,i,ispec)*mu(j,i,ispec))
                endif

                vy = veloc(1,iglob)
                ty = rho_vs*vy
                weight = jacobianb(ibb,j,ib)*wzgll(j)
                accel(1,iglob) = accel(1,iglob) - ty*weight
                if (save_forward) then
                   absorbfield(itime,1,j,ib,ibb) = ty*weight
                endif
              enddo
            enddo
          enddo

else  ! NCOMP == 3

       !
       !   INTEGRATION OVER SPECTRAL ELEMENTS
       !
       do ispec = 1,NSPEC

          ! first double loop over GLL
          ! compute and store gradients
          do j = 1,NGLLZ
             do i = 1,NGLLX

                iglob2 = ibool(i,j,ispec)

                ! derivative along x
                tempx1l = 0.
                tempy1l = 0.
                tempz1l = 0.
                if (solver_type == 3) then
                   b_tempx1l = 0.
                   b_tempy1l = 0.
                   b_tempz1l = 0.
                endif
                do k = 1,NGLLX
                   hp1 = hprime_xx(k,i)
                   iglob = ibool(k,j,ispec)
                   tempx1l = tempx1l + displ(1,iglob)*hp1
                   tempy1l = tempy1l + displ(2,iglob)*hp1
                   tempz1l = tempz1l + displ(3,iglob)*hp1
                   if (solver_type == 3) then
                      b_tempx1l = b_tempx1l + b_displ(1,iglob)*hp1
                      b_tempy1l = b_tempy1l + b_displ(2,iglob)*hp1
                      b_tempz1l = b_tempz1l + b_displ(3,iglob)*hp1
                   endif
                enddo

                ! derivative along z
                tempx2l = 0.
                tempy2l = 0.
                tempz2l = 0.
                if (solver_type == 3) then
                   b_tempx2l = 0.
                   b_tempy2l = 0.
                   b_tempz2l = 0.
                endif
                do k = 1,NGLLZ
                   hp2 = hprime_zz(k,j)
                   iglob = ibool(i,k,ispec)
                   tempx2l = tempx2l + displ(1,iglob)*hp2
                   tempy2l = tempy2l + displ(2,iglob)*hp2
                   tempz2l = tempz2l + displ(3,iglob)*hp2
                   if (solver_type == 3) then
                      b_tempx2l = b_tempx2l + b_displ(1,iglob)*hp2
                      b_tempy2l = b_tempy2l + b_displ(2,iglob)*hp2
                      b_tempz2l = b_tempz2l + b_displ(3,iglob)*hp2
                   endif
                enddo

                dxidxl = dxidx(i,j,ispec)
                dxidzl = dxidz(i,j,ispec)
                dgammadxl = dgammadx(i,j,ispec)
                dgammadzl = dgammadz(i,j,ispec)

                dsxdxl = tempx1l*dxidxl+tempx2l*dgammadxl
                dsxdzl = tempx1l*dxidzl+tempx2l*dgammadzl

                dsydxl = tempy1l*dxidxl+tempy2l*dgammadxl
                dsydzl = tempy1l*dxidzl+tempy2l*dgammadzl

                dszdxl = tempz1l*dxidxl+tempz2l*dgammadxl
                dszdzl = tempz1l*dxidzl+tempz2l*dgammadzl

                ds(1,1) = dsxdxl
                ds(1,2) = oneovertwo * dsydxl
                ds(1,3) = oneovertwo * (dszdxl + dsxdzl)
                ds(2,2) = 0
                ds(2,3) = oneovertwo * dsydzl
                ds(3,3) = dszdzl
                ds(2,1) = ds(1,2)
                ds(3,1) = ds(1,3)
                ds(3,2) = ds(2,3)

                if (solver_type == 3) then
                   b_dsxdxl = b_tempx1l*dxidxl+ b_tempx2l*dgammadxl
                   b_dsxdzl = b_tempx1l*dxidzl+ b_tempx2l*dgammadzl

                   b_dsydxl = b_tempy1l*dxidxl+ b_tempy2l*dgammadxl
                   b_dsydzl = b_tempy1l*dxidzl+ b_tempy2l*dgammadzl

                   b_dszdxl = b_tempz1l*dxidxl+ b_tempz2l*dgammadxl
                   b_dszdzl = b_tempz1l*dxidzl+ b_tempz2l*dgammadzl

                   b_ds(1,1) = b_dsxdxl
                   b_ds(1,2) = oneovertwo * b_dsydxl
                   b_ds(1,3) = oneovertwo * (b_dszdxl + b_dsxdzl)
                   b_ds(2,2) = 0
                   b_ds(2,3) = oneovertwo * b_dsydzl
                   b_ds(3,3) = b_dszdzl
                   b_ds(2,1) = b_ds(1,2)
                   b_ds(3,1) = b_ds(1,3)
                   b_ds(3,2) = b_ds(2,3)

                   ! kappa and mu kernels
                   kappa_k(iglob2) = (dsxdxl + dszdzl) * (b_dsxdxl + b_dszdzl)     ! product of div_s
                   mu_k(iglob2)    = sum(ds * b_ds) - ONE_THIRD * kappa_k(iglob2)  ! double-dot-product of deviatoric strains
                endif

                ! variation in kappa and mu
                kappal = kappa(i,j,ispec)
                mul = mu(i,j,ispec)
                lambdalplus2mul = kappal + FOUR_THIRDS * mul
                lambdal = lambdalplus2mul - 2.*mul

                sigma_xx = lambdalplus2mul*dsxdxl + lambdal*dszdzl
                sigma_xy = mul*dsydxl
                sigma_xz = mul*(dszdxl + dsxdzl)
                sigma_zx = sigma_xz
                sigma_zy = mul*dsydzl
                sigma_zz = lambdalplus2mul*dszdzl + lambdal*dsxdxl

                if (solver_type == 3) then
                   b_sigma_xx = lambdalplus2mul* b_dsxdxl + lambdal* b_dszdzl
                   b_sigma_xy = mul* b_dsydxl
                   b_sigma_xz = mul*(b_dszdxl + b_dsxdzl)
                   b_sigma_zx = b_sigma_xz
                   b_sigma_zy = mul* b_dsydzl
                   b_sigma_zz = lambdalplus2mul* b_dszdzl + lambdal* b_dsxdxl
                endif

                jacobianl = jacobian(i,j,ispec)

                ! non-symmetric form (first Piola-Kirchhoff stress)
                tempx1(i,j) = jacobianl*(sigma_xx*dxidxl+sigma_zx*dxidzl)
                tempy1(i,j) = jacobianl*(sigma_xy*dxidxl+sigma_zy*dxidzl)
                tempz1(i,j) = jacobianl*(sigma_xz*dxidxl+sigma_zz*dxidzl)

                tempx2(i,j) = jacobianl*(sigma_xx*dgammadxl+sigma_zx*dgammadzl)
                tempy2(i,j) = jacobianl*(sigma_xy*dgammadxl+sigma_zy*dgammadzl)
                tempz2(i,j) = jacobianl*(sigma_xz*dgammadxl+sigma_zz*dgammadzl)

                if (solver_type == 3) then
                   b_tempx1(i,j) = jacobianl*(b_sigma_xx*dxidxl+ b_sigma_zx*dxidzl)
                   b_tempy1(i,j) = jacobianl*(b_sigma_xy*dxidxl+ b_sigma_zy*dxidzl)
                   b_tempz1(i,j) = jacobianl*(b_sigma_xz*dxidxl+ b_sigma_zz*dxidzl)

                   b_tempx2(i,j) = jacobianl*(b_sigma_xx*dgammadxl+b_sigma_zx*dgammadzl)
                   b_tempy2(i,j) = jacobianl*(b_sigma_xy*dgammadxl+b_sigma_zy*dgammadzl)
                   b_tempz2(i,j) = jacobianl*(b_sigma_xz*dgammadxl+b_sigma_zz*dgammadzl)
                endif

             enddo
          enddo
          !
          ! second double-loop over GLL
          ! compute all rhs terms
          !
          do j = 1,NGLLZ
             do i = 1,NGLLX

                ! along x direction
                tempx1l = 0.
                tempy1l = 0.
                tempz1l = 0.
                if (solver_type == 3) then
                   b_tempx1l = 0.
                   b_tempy1l = 0.
                   b_tempz1l = 0.
                endif
                do k = 1,NGLLX
                   fac1 = wxgll(k)*hprime_xx(i,k)
                   tempx1l = tempx1l + tempx1(k,j)*fac1
                   tempy1l = tempy1l + tempy1(k,j)*fac1
                   tempz1l = tempz1l + tempz1(k,j)*fac1
                   if (solver_type == 3) then
                      b_tempx1l = b_tempx1l + b_tempx1(k,j)*fac1
                      b_tempy1l = b_tempy1l + b_tempy1(k,j)*fac1
                      b_tempz1l = b_tempz1l + b_tempz1(k,j)*fac1
                   endif
                enddo

                ! along z direction
                tempx2l = 0.
                tempy2l = 0.
                tempz2l = 0.
                if (solver_type == 3) then
                   b_tempx2l = 0.
                   b_tempy2l = 0.
                   b_tempz2l = 0.
                endif
                do k = 1,NGLLZ
                   fac2 = wzgll(k)*hprime_zz(j,k)
                   tempx2l = tempx2l + tempx2(i,k)*fac2
                   tempy2l = tempy2l + tempy2(i,k)*fac2
                   tempz2l = tempz2l + tempz2(i,k)*fac2
                   if (solver_type == 3) then
                      b_tempx2l = b_tempx2l + b_tempx2(i,k)*fac2
                      b_tempy2l = b_tempy2l + b_tempy2(i,k)*fac2
                      b_tempz2l = b_tempz2l + b_tempz2(i,k)*fac2
                   endif
                enddo

                fac1 = wzgll(j)
                fac2 = wxgll(i)

                ! key computation
                iglob = ibool(i,j,ispec)
                accel(1,iglob) = accel(1,iglob) - (fac1*tempx1l + fac2*tempx2l)
                accel(2,iglob) = accel(2,iglob) - (fac1*tempy1l + fac2*tempy2l)
                accel(3,iglob) = accel(3,iglob) - (fac1*tempz1l + fac2*tempz2l)

                if (solver_type == 3) then
                   b_accel(1,iglob) = b_accel(1,iglob) - (fac1* b_tempx1l + fac2* b_tempx2l)
                   b_accel(2,iglob) = b_accel(2,iglob) - (fac1* b_tempy1l + fac2* b_tempy2l)
                   b_accel(3,iglob) = b_accel(3,iglob) - (fac1* b_tempz1l + fac2* b_tempz2l)
                endif

             enddo ! second loop over the GLL points
          enddo

       enddo ! end loop over all spectral elements

       !
       ! boundary conditions
       !
       ! forward propagation

       !do ibb = 1, 3
       do ibb = 1,NABSORB  ! index of grid boundary (CHT)
          if (ibb == 1) then
             i = 1; nx = -1.; nz = 0.
          else if (ibb == 2) then
             i = NGLLX; nx = 1.; nz = 0.
          else if (ibb == 3) then
             i = 1; nx = 0.; nz = -1.
          else if (ibb == 4) then       ! CHT
             i = NGLLZ; nx = 0.; nz = 1.
          endif

          do ib = 1,nspecb(ibb)     ! elements on each grid boundary
             ispec = ibelm(ibb,ib)  ! (global) index of boundary element

             if (ibb == 1 .or. ibb == 2) then ! left or right boundary element
                j1 = 1; j2 = NGLLZ
             else if (ib == 1) then           ! top left corner element
                j1 = 2; j2 = NGLLX
             else if (ib == nspecb(ibb)) then ! top right corner element
                j1 = 1; j2 = NGLLX-1
             else                             ! top or bottom boundary (excluding corner elements)
                j1 = 1; j2 = NGLLX
             endif

             do j = j1, j2
                if (ibb == 1 .or. ibb == 2) then  ! left or right boundary
                   iglob = ibool(i,j,ispec)
                   rho_vp = dsqrt(rho(i,j,ispec)*(kappa(i,j,ispec)+FOUR_THIRDS*mu(i,j,ispec)))
                   rho_vs = dsqrt(rho(i,j,ispec)*mu(i,j,ispec))
                else                              ! top or bottom boundary
                   iglob = ibool(j,i,ispec)
                   rho_vp = dsqrt(rho(j,i,ispec)*(kappa(j,i,ispec)+FOUR_THIRDS*mu(j,i,ispec)))
                   rho_vs = dsqrt(rho(j,i,ispec)*mu(j,i,ispec))
                endif

                vx = veloc(1,iglob)
                vy = veloc(2,iglob)
                vz = veloc(3,iglob)

                vn = nx*vx+nz*vz

                tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                ty = rho_vs*vy
                tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

                weight = jacobianb(ibb,j,ib)*wzgll(j)

                accel(1,iglob) = accel(1,iglob) - tx*weight
                accel(2,iglob) = accel(2,iglob) - ty*weight
                accel(3,iglob) = accel(3,iglob) - tz*weight
                if (save_forward) then
                   absorbfield(itime,1,j,ib,ibb) = tx*weight
                   absorbfield(itime,2,j,ib,ibb) = ty*weight
                   absorbfield(itime,3,j,ib,ibb) = tz*weight
                endif
              enddo
            enddo
          enddo

endif  ! NCOMP

       ! backward propogation boundary condition
       if (solver_type == 3) then
          do ibb = 1,NABSORB     ! CHT : add 4th boundary
             if (ibb == 1) then
                i = 1
             else if (ibb == 2) then
                i = NGLLX
             else if (ibb == 3) then
                i = 1
             else if (ibb == 4) then
                i = NGLLZ
             endif

             ! see comments above
             do ib = 1,nspecb(ibb)
                if (ibb == 1 .or. ibb == 2) then
                   j1 = 1; j2 = NGLLZ
                else if (ib == 1) then
                   j1 = 2; j2 = NGLLX
                else if (ib == nspecb(ibb)) then
                   j1 = 1; j2 = NGLLX-1
                else
                   j1 = 1; j2 = NGLLX
                endif
                ispec = ibelm(ibb,ib)
                do j = j1, j2
                   if (ibb == 1 .or. ibb == 2) then
                      iglob = ibool(i,j,ispec)
                   else
                      iglob = ibool(j,i,ispec)
                   endif
                   b_accel(:,iglob) = b_accel(:,iglob) - absorbfield(NSTEP-itime+1,:,j,ib,ibb)
                enddo
             enddo
          enddo
       endif

       ! KEY: add source
       if (solver_type == 1) then
         ! forward wavefield
         do isrc = 1,nsrc

           ! perform the general interpolation using Lagrange polynomials -- NEW METHOD
           do j = 1,NGLLZ
             do i = 1,NGLLX
               iglob = ibool(i,j,ispec_src(isrc))
               hlagrange = hxis_store(isrc,i) * hgammas_store(isrc,j)

               accel(:,iglob) = accel(:,iglob) + samp(itime,:,isrc)*hlagrange
             enddo
           enddo

           ! take the value at the closest gridpoint -- OLD METHOD
           !iglob = sglob(isrc)
           !accel(:,iglob) = accel(:,iglob) + samp(itime,:,isrc)

         enddo  ! isrc

       else
         ! adjoint wavefield
         do irec = 1,nrec
           do j = 1,NGLLZ
             do i = 1,NGLLX
               iglob = ibool(i,j,ispec_rec(irec))
               hlagrange = hxir_store(irec,i) * hgammar_store(irec,j)

               accel(:,iglob) = accel(:,iglob) + ramp(NSTEP-itime+1,:,irec)*hlagrange
             enddo
           enddo
           !iglob = rglob(irec)
           !accel(:,iglob) = accel(:,iglob) + ramp(NSTEP-itime+1,:,irec)
         enddo

         ! forward wavefield, but computed in reverse
         if (solver_type == 3) then
           do isrc = 1,nsrc
             do j = 1,NGLLZ
               do i = 1,NGLLX
                 iglob = ibool(i,j,ispec_src(isrc))
                 hlagrange = hxis_store(isrc,i) * hgammas_store(isrc,j)

                 !b_accel(:,iglob) = b_accel(:,iglob) + samp(NSTEP-itime+1,:,isrc)*hlagrange
                 b_accel(:,iglob) = b_accel(:,iglob) + stf_for(NSTEP-itime+1,:,isrc)*hlagrange
               enddo
             enddo
             !iglob = sglob(isrc)
             !b_accel(:,iglob) = b_accel(:,iglob) + samp(NSTEP-itime+1,:,isrc)
           enddo
         endif

       endif

       ! above here, accel(:) are actually the RHS!
       ! divide by the mass matrix
       do i = 1,NGLOB
          accel(:,i) = accel(:,i)/mass_global(i)

          if (solver_type == 3) &
             b_accel(:,i) = b_accel(:,i)/mass_global(i)
       enddo

       ! `corrector' update velocity
       do i = 1,NCOMP
          veloc(i,:) = veloc(i,:) + deltatover2*accel(i,:)
          if (solver_type == 3) then
             b_veloc(i,:) = b_veloc(i,:) + b_deltatover2 * b_accel(i,:)
          endif
       enddo

       !--------------------------------
       ! save this time-step into seismograms and/or calculate kernels
       ! NEW METHOD -- perform the general interpolation using Lagrange polynomials
       ! OLD METHOD -- closest gridpoint

       if (solver_type == 1) then         ! forward seismograms

         do irec = 1,nrec

           ! perform the general interpolation using Lagrange polynomials -- NEW METHOD
           do j = 1,NGLLZ
             do i = 1,NGLLX
               iglob = ibool(i,j,ispec_rec(irec))
               hlagrange = hxir_store(irec,i) * hgammar_store(irec,j)

               ! equivalent to using closest gridpoint -- OLD METHOD
               !hlagrange = 0.
               !if (iglob == rglob(irec)) hlagrange = 1.

               ramp(itime,:,irec) = ramp(itime,:,irec) + displ(:,iglob)*hlagrange
             enddo
           enddo

           ! take the value at the closest gridpoint -- OLD METHOD
           !ramp(itime,:,irec) = displ(:,rglob(irec))

         enddo  ! irec


       ! adjoint seismograms
       ! note that samp above is copied to stf_for
       ! we fill the record 'in reverse' so that the max arrival will coincide
       ! with the forward source time function pulse
       else if (solver_type == 2 .or. solver_type == 3) then

         do isrc = 1,nsrc
           temp1 = 0; temp2 = 0. ; temp3 = 0.
           do j = 1,NGLLZ      ! loop over NGLLX by NGLLZ points for each source
             do i = 1,NGLLX
               iglob = ibool(i,j,ispec_src(isrc))
               hlagrange = hxis_store(isrc,i) * hgammas_store(isrc,j)

               ! assume one component only for now
               temp1 = temp1 + displ(1,iglob)*hlagrange        ! sdag
               temp2 = temp2 + displ_grad(i,j,1)*hlagrange     ! d/dx(sdag)
               temp3 = temp3 + displ_grad(i,j,2)*hlagrange     ! d/dz(sdag)
             enddo
           enddo

           !samp(NSTEP-itime+1,:,isrc) = displ(:,sglob(isrc))

           ! adjoint seismograms (time-reversed)
           samp(NSTEP-itime+1,1,isrc) = temp1

           ! source perturbation time series
           ! (assume one component for now)
           three_source_model(NSTEP-itime+1,1,nsrc,1) =  temp2 * f0(1) * stf(NSTEP-itime+1)   ! dx
           three_source_model(NSTEP-itime+1,1,nsrc,2) =  temp3 * f0(1) * stf(NSTEP-itime+1)   ! dz
           three_source_model(NSTEP-itime+1,1,nsrc,3) = -temp1 * f0(1) * stf_d(NSTEP-itime+1) ! dt0

           ! additional time series for checking
           !three_source_model(NSTEP-itime+1,1,nsrc,4) = temp1
           !three_source_model(NSTEP-itime+1,1,nsrc,5) = temp2
           !three_source_model(NSTEP-itime+1,1,nsrc,6) = temp3
           !three_source_model(itime,1,nsrc,7) = stf_d(itime)
           !three_source_model(itime,1,nsrc,8) = stf(itime)

         enddo  ! isrc


         if (solver_type == 3) then

           ! CALCULATE SIX KERNELS -- notice the time integration
           do iglob = 1, NGLOB
             ! note the minus sign in these expressions (TTL,2005)
             rho_kernel(iglob) = rho_kernel(iglob) - &
                  rho_global(iglob) * dot_product(accel(:,iglob),b_displ(:,iglob)) * DT
             mu_kernel(iglob) = mu_kernel(iglob) - &
                  2 * mu_global(iglob) * mu_k(iglob) * DT
             kappa_kernel(iglob) = kappa_kernel(iglob) - &
                  kappa_global(iglob) * kappa_k(iglob) * DT

             rhop_kernel(iglob) = rho_kernel(iglob) + kappa_kernel(iglob) + mu_kernel(iglob)
             beta_kernel(iglob) = 2 * (mu_kernel(iglob) - &
                  (4.0 * mu_global(iglob))/(3.0 * kappa_global(iglob)) * kappa_kernel(iglob))
             alpha_kernel(iglob) = 2 * (3 * kappa_global(iglob) + 4 * mu_global(iglob)) &
                  / (3 * kappa_global(iglob)) * kappa_kernel(iglob)

             ! SH beta-kernel INTERACTION field (for membrane surface waves)
             beta_kernel_prime(iglob) = -2. * (2. * mu_global(iglob) * mu_k(iglob) )
             !beta_kernel(iglob) = beta_kernel(iglob) + beta_kernel_prime(iglob) * DT
           enddo
         endif  ! solver_type = 3

       endif  ! solver_type

       ! save last frame
       if (mod(itime, NSAVE) == 0) then

          if (save_forward) then
             open(11,file = trim(last_frame),status='unknown',iostat=ios)
             if (ios /= 0) stop 'Error reading the last frame'
             do i = 1,NGLOB
                write(11,fm) (sngl(displ(j,i)), j = 1,NCOMP), &
                             (sngl(veloc(j,i)), j = 1,NCOMP), &
                             (sngl(accel(j,i)), j = 1,NCOMP)
             enddo
             close(11)
          endif

          ! ONE KERNEL: for basis function analysis
          ! take the integrated kernel at the FINAL TIME STEP
          if (solver_type == 3 .and. itime == NSTEP) then
             filename5 = trim(out_dir)//'kernel_basis'
             open(unit = 13, file = trim(filename5), status = 'unknown',iostat=ios)
             if (ios /= 0) stop 'Error writing snapshot to disk'
             do iglob = 1, NGLOB
                write(13,'(3e16.6)') x_lon(iglob), z_lat(iglob), sngl(beta_kernel(iglob))
                !write(13,'(5e16.6)') x(iglob), z(iglob), x_lon(iglob), z_lat(iglob), sngl(beta_kernel(iglob))
             enddo
             close(13)

             ! sum these in wave2d.f90
             kernel_basis(:) = beta_kernel(:)
          endif  ! solver_type == 3

       endif

       !stop 'testing'

       !===================================

       ! write snapshots or kernels
       if (itime == 1 .or. mod(itime, NSAVE) == 0) then   ! save the initial frame
          if (itime == 1) tlab = 0
          if (itime /= 1) tlab = itime
       !if (mod(itime, NSAVE) == 0) then
          if (solver_type == 1) then
             write(filename1,'(a,i5.5)') trim(out_dir)//'forward_',tlab
          else if (solver_type == 2) then
             write(filename1,'(a,i5.5)') trim(out_dir)//'adjoint_',tlab
          else
             write(filename1,'(a,i5.5)') trim(out_dir)//'adjoint_',tlab      ! adjoint wavefield
             write(filename2,'(a,i5.5)') trim(out_dir)//'backward_',tlab
             write(filename3,'(a,i5.5)') trim(out_dir)//'kernel0_',tlab      ! kernels rho-kappa-mu
             write(filename4,'(a,i5.5)') trim(out_dir)//'kernel_',tlab       ! kernels rho'-alpha-beta
             !write(filename6,'(a,i5.5)') trim(out_dir)//'interaction_',tlab  ! interaction
          endif

          ! x-z coordinates for plotting
          !!xtemp = x(iglob)/LENGTH ; ztemp = z(iglob)/LENGTH
          !xtemp = x(iglob)/1000. ; ztemp = z(iglob)/1000.

          if (WRITE_SNAPSHOTS) then      ! wavefield snapshots
             open(unit=11, file=trim(filename1), status='unknown', iostat=ios)
             if (ios /= 0) stop 'Error writing snapshot to disk'
             do iglob = 1, NGLOB
                xtemp = x_lon(iglob) ; ztemp = z_lat(iglob)
                write(11,'(5e16.6)') sngl(xtemp), sngl(ztemp), (sngl(displ(j,iglob)),j=1,NCOMP)
             enddo
             close(11)
          endif

          if (WRITE_KERNELS) then        ! kernel snapshots
          if (solver_type == 3) then
             open(unit=11, file=trim(filename2), status='unknown', iostat=ios)
             if (ios /= 0) stop 'Error writing snapshot to disk'
             do iglob = 1, NGLOB
                xtemp = x_lon(iglob) ; ztemp = z_lat(iglob)
                write(11,'(5e16.6)') sngl(xtemp), sngl(ztemp), (sngl(b_displ(j,iglob)),j=1,NCOMP)
             enddo
             close(11)
             open(unit = 11, file = trim(filename3), status = 'unknown',iostat=ios)
             open(unit = 12, file = trim(filename4), status = 'unknown',iostat=ios)
             if (ios /= 0) stop 'Error writing snapshot to disk'

             do iglob = 1, NGLOB

                xtemp = x_lon(iglob) ; ztemp = z_lat(iglob)

                if (ISURFACE == 0) then
                   ! six kernels into two files
                   write(11,'(5e16.6)') sngl(xtemp), sngl(ztemp), sngl(rho_kernel(iglob)), &
                              sngl(mu_kernel(iglob)), sngl(kappa_kernel(iglob))
                   write(12,'(5e16.6)') sngl(xtemp), sngl(ztemp), sngl(rhop_kernel(iglob)), &
                              sngl(beta_kernel(iglob)), sngl(alpha_kernel(iglob))
                else
                   ! x, z, interaction, beta
                   write(12,'(4e16.6)') sngl(xtemp), sngl(ztemp), &
                              sngl(beta_kernel_prime(iglob)), sngl(beta_kernel(iglob))
                endif
             enddo
             close(11)
             close(12)

          endif  ! solver_type == 3
          endif  ! WRITE_KERNELS

       endif  ! write snapshots or kernels

    enddo ! end time loop

    if (solver_type == 3) deallocate(stf_for)

  end subroutine solver
!---------------------------------------------

end module wave2d_solver
