module wave2d_solver

  use wave2d_constants
  use wave2d_variables
  use wave2d_define_der_matrices
  use wave2d_sub    ! for plotting using subroutine local2global.f90

  implicit none

contains

  subroutine mesher

    integer ispec,ib,i,j,k,iglob,iglob1,itime,ix,iz

    double precision :: ztemp, z1temp, z2temp
    double precision :: ztemp1, ztemp2, dtemp, dtrsh, dinc, zdiff, dz1, dz2
    integer :: iel, ik

    ! set up grid and spectral elements
    call define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,wgllwgll_xz)

    !===================

!!$    dtrsh = 0.1   ! level of accuracy for the position of the GLL points w.r.t. the discontinuity
!!$    dinc  = 0.1   ! how much to perturb the top of the element
!!$
!!$    z2temp = 0.0
!!$    ispec = 0
!!$
!!$    ! elements in Z dimension
!!$    do iz = 1,NEZ
!!$       do ix = 1,NEX
!!$          ispec = ispec+1
!!$
!!$      if (ix==1) then
!!$
!!$      ! evenly spaced anchors between 0 and 1
!!$      !z1temp = HEIGHT*dble(iz-1)/dble(NEZ)
!!$      !z2temp = HEIGHT*dble(iz)/dble(NEZ)
!!$
!!$      z1temp = z2temp
!!$      z2temp = z2temp + HEIGHT/dble(NEZ)    ! standard element height
!!$      if (iz==NEZ) z2temp = HEIGHT
!!$
!!$      ! determine whether the element is cut by a discontinuity
!!$      iel = 0
!!$      do k = 1,nlayer-1
!!$        dtemp = HEIGHT - z_breaks(k)
!!$        if (dtemp > z1temp .and. dtemp < z2temp) then
!!$          print *, k, iz
!!$          iel = 1
!!$          ik = k
!!$        endif
!!$      enddo
!!$      !iel = 0
!!$
!!$      ! perturb the boundary of the top of the element (z2)
!!$      if (iel == 1) then
!!$
!!$         zdiff = HEIGHT*LENGTH  ! initial difference
!!$         do while (zdiff > dtrsh)
!!$            !write(*,'(5e18.8)') z2temp, zdiff, dtrsh, dz1, dz2
!!$
!!$            do j = 1,NGLLZ
!!$               ! present GLL height (upper layer)
!!$               ztemp2 = 0.5*(1.0-zigll(j))*z1temp + 0.5*(1.0+zigll(j))*z2temp
!!$
!!$               if (j > 1) then
!!$                  ! previous GLL z height (lower layer)
!!$                  ztemp1 = 0.5*(1.0-zigll(j-1))*z1temp + 0.5*(1.0+zigll(j-1))*z2temp
!!$
!!$                  dtemp = HEIGHT - z_breaks(ik)   ! height of boundary
!!$
!!$                  ! discontinuity is inbetween GLL layers
!!$                  if (dtemp > ztemp1 .and. dtemp < ztemp2) then
!!$                     dz1 = abs(dtemp - ztemp1)
!!$                     dz2 = abs(ztemp2 - dtemp)
!!$                     zdiff = abs(dz2 - dz1)
!!$                     !if (zdiff > dtrsh) then
!!$                     !   if (dz1 > dz2) z2temp = z2temp + dinc  ! push it up
!!$                     !   if (dz1 < dz2) z2temp = z2temp - dinc  ! push it down
!!$                     !endif
!!$                  endif
!!$
!!$               endif
!!$            enddo
!!$
!!$            ! changes the position of the top of the element
!!$            if (dz1 > dz2) z2temp = z2temp + dinc  ! push it up
!!$            if (dz1 < dz2) z2temp = z2temp - dinc  ! push it down
!!$
!!$         enddo
!!$       endif
!!$
!!$       ! compute the GLL points for this element
!!$       do j = 1,NGLLZ
!!$         ztemp2 = 0.5*(1.0-zigll(j))*z1temp + 0.5*(1.0+zigll(j))*z2temp
!!$         print *, iz, j, ztemp2
!!$       enddo
!!$
!!$       endif  ! ix==1
!!$
!!$       z1(ispec) = z1temp
!!$       z2(ispec) = z2temp
!!$
!!$    enddo
!!$    enddo

    !stop 'testing mesher'

    !===================

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
          z1(ispec) = HEIGHT*dble(iz-1)/dble(NEZ)  ! see above
          z2(ispec) = HEIGHT*dble(iz)/dble(NEZ)    ! see above

          ! loop over GLL points to calculate jacobian, and set up numbering
          !
          ! jacobian = | dx/dxi dx/dgamma | = (z2-z1)*(x2-x1)/4  as dx/dgamma=dz/dxi = 0
          !            | dz/dxi dz/dgamma |
          !
          do j = 1,NGLLZ
             do i = 1,NGLLX

                ! jacobian, integration weight
                dxidx(i,j,ispec) = 2.0 / (x2(ispec)-x1(ispec))
                dxidz(i,j,ispec) = 0.0
                dgammadx(i,j,ispec) = 0.0
                dgammadz(i,j,ispec) = 2.0 / (z2(ispec)-z1(ispec))
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
                x(iglob1) = 0.5*(1.0-xigll(i))*x1(ispec) + 0.5*(1.0+xigll(i))*x2(ispec)
                z(iglob1) = 0.5*(1.0-zigll(j))*z1(ispec) + 0.5*(1.0+zigll(j))*z2(ispec)

                ! end loop over GLL points
             enddo
          enddo

          ! if boundary element
          ! 1,2,3,4 --> left, right, bottom, top
          if (ix == 1) then      ! left boundary
             nspecb(1) = nspecb(1) + 1
             ibelm(1,nspecb(1)) = ispec
             do j = 1,NGLLZ
                jacobianb(1,j,nspecb(1))= (z2(ispec)-z1(ispec))/2.0
             enddo
          endif
          if (ix == NEX) then    ! right boundary
             nspecb(2) = nspecb(2) + 1
             ibelm(2,nspecb(2)) = ispec
             do j = 1,NGLLZ
                jacobianb(2,j,nspecb(2))= (z2(ispec)-z1(ispec))/2.0
             enddo
          endif
          if (iz == 1) then      ! bottom boundary
             nspecb(3) = nspecb(3) + 1
             ibelm(3,nspecb(3)) = ispec
             do i = 1,NGLLX
                jacobianb(3,i,nspecb(3))= (x2(ispec)-x1(ispec))/2.0
             enddo
          endif
          if (iz == NEZ) then    ! top boundary
             nspecb(4) = nspecb(4) + 1
             ibelm(4,nspecb(4)) = ispec
             do i = 1,NGLLX
                jacobianb(4,i,nspecb(4))= (x2(ispec)-x1(ispec))/2.0
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

!!$    !c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
!!$    !time_step = 0.2*dh/c
!!$    print *
!!$    print *,'                       space step (km) :', sngl(dh/1000.0)
!!$    if (ISURFACE==0) then
!!$       print *,'time step est from courant = 0.2, Pmax : ',sngl(0.2*dh/alpha_max),' seconds'
!!$       print *,'time step est from courant = 0.2, Pmin : ',sngl(0.2*dh/alpha_min),' seconds'
!!$       print *,'time step est from courant = 0.2, Smax : ',sngl(0.2*dh/beta_max),' seconds'
!!$       print *,'time step est from courant = 0.2, Smin : ',sngl(0.2*dh/beta_min),' seconds'
!!$       print *,'                      actual time step : ',sngl(DT),' seconds'
!!$    endif

  end subroutine mesher

  !-------------------------------------------------------

  subroutine mesher_additional

    logical, dimension(NGLOB) :: mask_ibool
    integer :: iglob1, iglob2, iglob3, iglob4
    integer :: ispec,i,j,k,iglob

    !--------------

    ! da vector and valence vector
    da_global(:) = 0.0
    valence(:) = 0
    k = 0
    do ispec = 1,NSPEC
       do j = 1,NGLLZ
          do i = 1,NGLLX
             iglob = ibool(i,j,ispec)
             da_local(i,j,ispec) = wxgll(i)*wzgll(j)*jacobian(i,j,ispec)

             k = k+1
             da_local_vec(k) = da_local(i,j,ispec)

             !da_global(iglob) = da_global(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
             da_global(iglob) = da_global(iglob) + da_local(i,j,ispec)
             valence(iglob) = valence(iglob) + 1
          enddo
       enddo
    enddo

    ! GLL points defining the corners of elements
    mask_ibool(:) = .false.
    do ispec = 1,NSPEC
       iglob1 = ibool(1,1,ispec)
       iglob2 = ibool(NGLLX,1,ispec)
       iglob3 = ibool(1,NGLLZ,ispec)
       iglob4 = ibool(NGLLX,NGLLZ,ispec)

       if (.not. mask_ibool(iglob1)) mask_ibool(iglob1) = .true.
       if (.not. mask_ibool(iglob2)) mask_ibool(iglob2) = .true.
       if (.not. mask_ibool(iglob3)) mask_ibool(iglob3) = .true.
       if (.not. mask_ibool(iglob4)) mask_ibool(iglob4) = .true.
    enddo

    k = 0
    ielement_corner(:) = 0
    do iglob = 1,NGLOB
       if (mask_ibool(iglob)) then
          k = k+1
          ielement_corner(k) = iglob
       endif
    enddo

  end subroutine mesher_additional

  !-------------------------------------------------------

  subroutine set_model_property(iref)

    ! This subroutine assigns the structure model parameters for the reference model (synthetics)
    ! and for the target model (data).  All assignments are done LOCALLY, so that discontinuities
    ! (such as the Moho) may be honored in the simulations.

    integer, intent(in) :: iref
    integer :: ispec, i, j, iglob

    !double precision :: socal_z1, socal_z2, socal_z3, socal_z4
    double precision :: dalpha, dbeta, drho, dmu, dkappa, dalpha2, dbeta2, drho2, dmu2, dkappa2
    double precision :: atemp, btemp, rtemp, ktemp, mtemp
    double precision :: btemp2

    !------------------

    ! commented out 06-Aug-2006
!!$    ! model properties -- globally defined (no discontinuities permitted this way)
!!$    do iglob = 1,NGLOB
!!$       rho_global(iglob)   = DENSITY
!!$       kappa_global(iglob) = INCOMPRESSIBILITY
!!$       if (ISURFACE==0) then
!!$          mu_global(iglob) = RIGIDITY
!!$       else
!!$          ! KEY: this means that the S velocity will be the surface wave phase velocity (m/s)
!!$          mu_global(iglob) = DENSITY*(c_glob(iglob))**2
!!$          !if (ihomo==1) mu_global(iglob) = DENSITY*beta0**2
!!$          !if (ihomo==0) mu_global(iglob) = DENSITY*(c_glob(iglob))**2
!!$       endif
!!$    enddo

    ! target model (data): perturbations are entered as percentages ( dm/m0 )
    dalpha = PALPHA/100.0
    dbeta  = PBETA/100.0
    drho   = PRHO/100.0

    print *
    print *, 'setting model properties at a local scale'
    if (iref == 1) print *, ' -- > reference model (synthetics), IMODEL_SYN = ', IMODEL_SYN
    if (iref == 0) print *, ' -- > target model (data), IMODEL_DAT = ', IMODEL_DAT
    if (IMODEL_SYN > 3) stop 'IMODEL_SYN must be 0,1,2,3'
    if (IMODEL_DAT > 3) stop 'IMODEL_DAT must be 0,1,2,3'

    !-----------------
    ! fill local arrays of the reference structure model

    do ispec = 1,NSPEC
       do j = 1,NGLLZ
          do i = 1,NGLLX
             iglob = ibool(i,j,ispec)

             ! structure for reference model (synthetics)
             if (iref == 1) then

                if (IMODEL_SYN == 0) then      ! homogeneous model

                   !kappa_syn(i,j,ispec) = INCOMPRESSIBILITY
                   !mu_syn(i,j,ispec)    = RIGIDITY
                   !rho_syn(i,j,ispec)   = DENSITY

                   kappa_syn(i,j,ispec) = DENSITY*(alpha0*alpha0 - FOUR_THIRDS*beta0*beta0)
                   mu_syn(i,j,ispec)    = DENSITY*beta0*beta0
                   rho_syn(i,j,ispec)   = DENSITY

                else if (IMODEL_SYN == 1) then  ! 1D model (body waves)

                   if (ISURFACE == 0) then
                      call make_1D_model(i,j,ispec,ktemp,mtemp,rtemp)       ! 1D model
                      kappa_syn(i,j,ispec) = ktemp
                      mu_syn(i,j,ispec)    = mtemp
                      rho_syn(i,j,ispec)   = rtemp
                   else
                      stop 'check model in set_model_property.f90'
                   endif

                else if (IMODEL_SYN == 2) then  ! checkerboard model
                   if (ISURFACE == 0) then
                      stop 'check model in set_model_property.f90'

                   else if (ISURFACE == 1) then
                      ! checkerboard S-wave velocity (= membrane wave phase velocity)
                      ! GLOBALLY defined (not at the elemental level)
                      btemp2 = beta0 * (1.0 + afac/100.0*(sin(x(iglob)*w_scale) * sin(z(iglob)*w_scale)) )

                      kappa_syn(i,j,ispec) = INCOMPRESSIBILITY
                      mu_syn(i,j,ispec)    = DENSITY * btemp2*btemp2
                      rho_syn(i,j,ispec)   = DENSITY
                   endif

                else if (IMODEL_SYN == 3) then  ! heterogeneous model
                   if (ISURFACE == 0) then
                      stop 'check model in set_model_property.f90'

                   else
                      kappa_syn(i,j,ispec) = INCOMPRESSIBILITY    ! not used
                      mu_syn(i,j,ispec)    = DENSITY * beta_syn(i,j,ispec)*beta_syn(i,j,ispec)
                      rho_syn(i,j,ispec)   = DENSITY

                   endif

                endif

             ! structure for target model (data)
             ! THE REFERENCE MODEL MUST HAVE BEEN PREVIOUSLY CONSTRUCTED
             else

                ! reference model values (m0)
                ktemp = kappa_syn(i,j,ispec)
                mtemp = mu_syn(i,j,ispec)
                rtemp = rho_syn(i,j,ispec)
                atemp = sqrt( (ktemp + FOUR_THIRDS*mtemp) / rtemp )
                btemp = sqrt( mtemp / rtemp )

                if (IMODEL_DAT == 0) then               ! uniform perturbation

                   ! perturbations ( i.e., dm = m-m0 , NOT dm/m0 )
                   dalpha2 = atemp*dalpha
                   dbeta2 = btemp*dbeta
                   drho2 = rtemp*drho

                   ! Dahlen and Tromp, p.332
                   dkappa2 = drho2*(atemp*atemp - FOUR_THIRDS*btemp*btemp) &
                        + 2*rtemp*(atemp*dalpha2 - FOUR_THIRDS*btemp*dbeta2)
                   dmu2 = drho2*btemp*btemp + 2*rtemp*btemp*dbeta2

                   ! perturbed values
                   !alpha_dat(i,j,ispec) = atemp + dalpha2
                   !beta_dat(i,j,ispec)  = btemp + dbeta2
                   rho_dat(i,j,ispec)   = rtemp + drho2
                   kappa_dat(i,j,ispec) = ktemp + dkappa2
                   mu_dat(i,j,ispec)    = mtemp + dmu2

                else if (IMODEL_DAT == 1) then         ! 1D model (NOT a 1D perturbation)

                   if (ISURFACE == 0) then
                      call make_1D_model(i,j,ispec,ktemp,mtemp,rtemp)       ! 1D model

                   else
                      stop 'check model in set_model_property.f90'
                   endif
                   kappa_dat(i,j,ispec) = ktemp
                   mu_dat(i,j,ispec)    = mtemp
                   rho_dat(i,j,ispec)   = rtemp

                else if (IMODEL_DAT == 2) then         ! checkerboard perturbation

                   if (ISURFACE == 0) then
                      stop 'check model in set_model_property.f90'
                   else
                      ! GJI-2007 paper
                      ! perturbed S-wave velocity (= membrane wave phase velocity)
                      ! GLOBALLY defined (not at the elemental level)
                      btemp2 = btemp * (1.0 + afac/100.0*(sin(x(iglob)*w_scale) * sin(z(iglob)*w_scale)) )

                      kappa_dat(i,j,ispec) = kappa_syn(i,j,ispec)
                      mu_dat(i,j,ispec)    = DENSITY * btemp2*btemp2
                      rho_dat(i,j,ispec)   = rho_syn(i,j,ispec)
                      !kappa_dat(i,j,ispec) = INCOMPRESSIBILITY
                      !mu_dat(i,j,ispec)    = DENSITY * btemp2*btemp2
                      !rho_dat(i,j,ispec)   = DENSITY
                   endif

                else if (IMODEL_DAT == 3) then         ! heterogeneous perturbation

                   if (ISURFACE == 0) then
                      stop 'check model in set_model_property.f90'

                   else
                      kappa_dat(i,j,ispec) = INCOMPRESSIBILITY    ! not used
                      mu_dat(i,j,ispec)    = DENSITY * beta_dat(i,j,ispec)*beta_dat(i,j,ispec)
                      rho_dat(i,j,ispec)   = DENSITY

                   endif

                endif   ! IMODEL_DAT
             endif   ! iref

          enddo
       enddo
    enddo

    if (iref == 1) then
       alpha_syn = 0.0 ; beta_syn = 0.0 ; bulk_syn = 0.0
       alpha_syn = sqrt( (kappa_syn + FOUR_THIRDS*mu_syn) / rho_syn )
       beta_syn  = sqrt( mu_syn / rho_syn )
       bulk_syn  = sqrt( kappa_syn / rho_syn )
    else
       alpha_dat = 0.0 ; beta_dat = 0.0 ; bulk_dat = 0.0
       alpha_dat = sqrt( (kappa_dat + FOUR_THIRDS*mu_dat) / rho_dat )
       beta_dat  = sqrt( mu_dat / rho_dat )
       bulk_dat  = sqrt( kappa_dat / rho_dat )
    endif

  end subroutine set_model_property

  !-------------------------------------------------------

  subroutine make_1D_model(i,j,ispec,kappa1,mu1,rho1)

    ! Given a GLL local index (i,j,ispec), this returns the structure values kappa-mu-rho,
    ! for a 1D model that honors discontinuities in the mesh.

    integer, intent(in) :: i,j,ispec
    double precision, intent(out) :: kappa1, mu1, rho1
    double precision :: atemp, btemp, rtemp, dtemp
    integer :: iglob

    !-----------------

    iglob = ibool(i,j,ispec)
    dtemp = HEIGHT - z(iglob)     ! depth of GLL point

    if (dtemp < z_breaks(1)) then                                  ! shallow surface
       rtemp = r_layers(1) ; atemp = a_layers(1) ; btemp = b_layers(1)

    else if (dtemp == z_breaks(1)) then                             ! 'basin' reflector
       if (HEIGHT-z1(ispec) == dtemp) then
          rtemp = r_layers(1) ; atemp = a_layers(1) ; btemp = b_layers(1)
       else if (HEIGHT-z2(ispec) == dtemp) then
          rtemp = r_layers(2) ; atemp = a_layers(2) ; btemp = b_layers(2)
       else
          print *, dtemp, z_breaks(1), z1(ispec), z2(ispec)
          stop 'error in make_1D_model.f90'
       endif

    else if (dtemp > z_breaks(1) .and. dtemp < z_breaks(2)) then    ! upper crust
       rtemp = r_layers(2) ; atemp = a_layers(2) ; btemp = b_layers(2)

    else if (dtemp == z_breaks(2)) then                             ! mid-crust reflector
       if (HEIGHT-z1(ispec) == dtemp) then
          rtemp = r_layers(2) ; atemp = a_layers(2) ; btemp = b_layers(2)
       else if (HEIGHT-z2(ispec) == dtemp) then
          rtemp = r_layers(3) ; atemp = a_layers(3) ; btemp = b_layers(3)
       else
          print *, dtemp, z_breaks(2), z1(ispec), z2(ispec)
          stop 'error in make_1D_model.f90'
       endif

    else if (dtemp > z_breaks(2) .and. dtemp < z_breaks(3)) then    ! lower crust
       rtemp = r_layers(3) ; atemp = a_layers(3) ; btemp = b_layers(3)

    else if (dtemp == z_breaks(3)) then                             ! Moho
       if (HEIGHT-z1(ispec) == dtemp) then
          rtemp = r_layers(3) ; atemp = a_layers(3) ; btemp = b_layers(3)
       else if (HEIGHT-z2(ispec) == dtemp) then
          rtemp = r_layers(4) ; atemp = a_layers(4) ; btemp = b_layers(4)
       else
          print *, dtemp, z_breaks(3), z1(ispec), z2(ispec)
          stop 'error in make_1D_model.f90'
       endif

    else                                                          ! upper mantle
       rtemp = r_layers(4) ; atemp = a_layers(4) ; btemp = b_layers(4)
    endif

    ! compute kappa, mu, rho
    kappa1 = rtemp*(atemp*atemp - FOUR_THIRDS*btemp*btemp)
    mu1    = rtemp*btemp*btemp
    rho1   = rtemp

  end subroutine make_1D_model

  !---------------------------------------------

  subroutine solver(solver_type, idata, &
       nsrc, sglob, ispec_src, hxis_store, hgammas_store, samp, &
       nrec, rglob, ispec_rec, hxir_store, hgammar_store, ramp, &
       last_frame, absorbfield, &
       atype_kernel, btype_kernel, rtype_kernel, &
       three_source_model, stf_syn, f0)

    ! required input and output
    integer, intent(in) :: solver_type, idata
    integer, intent(in) :: nsrc, sglob(nsrc), ispec_src(nsrc)
    integer, intent(in) :: nrec, rglob(nrec), ispec_rec(nrec)
    double precision, intent(in) :: hxir_store(nrec,NGLLX), hgammar_store(nrec,NGLLZ)
    double precision, intent(in) :: hxis_store(nsrc,NGLLX), hgammas_store(nsrc,NGLLZ)
    double precision, intent(inout) :: samp(NSTEP,NCOMP,nsrc)
    double precision, intent(inout) :: ramp(NSTEP,NCOMP,nrec)

    ! OPTIONAL ARGUMENTS
    character(len=*), optional :: last_frame
    double precision, intent(inout), optional :: absorbfield(NSTEP, NCOMP, NGLL, NELE, NABSORB)
    double precision, intent(out), optional :: three_source_model(NSTEP,NCOMP,nsrc,NVAR_SOURCE)   ! 3 : ts, xs, zs
    !double precision, intent(out), optional :: three_source_model(NSTEP,NCOMP,nsrc,10)    ! testing
    double precision, intent(out), dimension(NGLLX,NGLLZ,NSPEC), optional :: rtype_kernel, btype_kernel, atype_kernel
    double precision, intent(in), optional :: stf_syn(NSTEP), f0(NCOMP)

    ! if adjoint seismograms are desired while computing kernels, then we need to make sure
    ! not to overwrite the forward source time function
    double precision,  dimension(:,:,:), allocatable :: stf_for
    double precision, dimension(NSTEP) :: stf_syn_dt

    ! displacement gradients
    ! only save the GLL points on the element where the source is located
    double precision :: displ_grad(NGLLX,NGLLZ,2*NCOMP)

    ! kernels and interaction fields
    double precision, dimension(:,:,:), allocatable :: &
       kappa_mu_rho_kernel, mu_kappa_rho_kernel, rho_kappa_mu_kernel, &
       alpha_beta_rho_kernel, beta_alpha_rho_kernel, rho_alpha_beta_kernel, &
       c_beta_rho_kernel, beta_c_rho_kernel, rho_c_beta_kernel, &
       kappa_mu_rho_kernel_int, mu_kappa_rho_kernel_int, rho_kappa_mu_kernel_int, &
       alpha_beta_rho_kernel_int, beta_alpha_rho_kernel_int, rho_alpha_beta_kernel_int, &
       c_beta_rho_kernel_int, beta_c_rho_kernel_int, rho_c_beta_kernel_int

    ! CHT
    integer tlab
    double precision :: temp1,temp2,temp3

    ! Lagrange interpolation
    double precision :: hlagrange

    integer ispec,ib,i,j,k,iglob,iglob1,iglob2,itime,ix,iz,itime1,itime2,itime0
    integer isrc, irec, icomp, isave
    character(len=100) :: filename,filename1,filename2,filename3,filename4,filename5,filename6,fm
    !double precision, dimension(NGLOB) :: mu_k, kappa_k
    double precision, dimension(NGLLX,NGLLZ,NSPEC) :: mu_k, kappa_k
    logical :: save_forward

    !--------------------------------------

    if (NCOMP == 3) fm = '(9e12.3)'
    if (NCOMP == 1) fm = '(3e12.3)'

    ! test of input arguments
    if (solver_type /= 1 .and. solver_type /= 2 .and. solver_type /= 3) then
       stop 'solver_type has to be 1, 2 or 3'
    endif

    if (idata /= 1 .and. idata /= 0) then
       stop 'idata has to be 0 or 1'
    endif

    save_forward = .false.
    if (solver_type == 1) then
       if (present(last_frame) .and. present(absorbfield)) save_forward = .true.
    endif

    if (solver_type == 3) then
       if (.not. (present(last_frame) .and. present(absorbfield)  &
            .and. present(rtype_kernel) .and. present(btype_kernel) .and. present(atype_kernel))) &
            stop 'For kernel calculation, last_frame, absorbfield and all kernel are in the argument'

       ! interaction fields must be initialized
       mu_k = 0.0
       kappa_k = 0.0

       ! allocate a new array for the forward source function (note: includes amplitude)
       allocate(stf_for(NSTEP,NCOMP,nsrc))
       stf_for = samp

       ! compute the derivative of the source time function
       stf_syn_dt(:) = 0.0
       do i = 2, NSTEP-1
          stf_syn_dt(i) =  (stf_syn(i+1) - stf_syn(i-1)) / (2.0*DT)
       enddo
       stf_syn_dt(1) = (stf_syn(2) - stf_syn(1)) / DT
       stf_syn_dt(NSTEP) = (stf_syn(NSTEP) - stf_syn(NSTEP-1)) / DT

       ! initialize adjoint seismograms
       samp = 0.0

    endif

    ! gridpoints per wavelength estimation -- based on TARGET model (data)
    print *
    print *, 'space step (km):', sngl(dh/1000.0)
    if (ISURFACE == 1) then
       print *, 'wavelength-min (km):', sngl(2*hdur*beta_min/1000.0)
       print *, 'wavelength-max (km):', sngl(2*hdur*beta_max/1000.0)
       print *, 'number of gridpoints per wavelength for S:'
       print *, '  min (cmin) : ', floor(2*hdur*beta_min/dh)
       print *, '  max (cmax) : ', floor(2*hdur*beta_max/dh)
    else
       !c = sqrt((INCOMPRESSIBILITY+FOUR_THIRDS*RIGIDITY)/DENSITY)
       !print *, 'number of gridpoints per wavelength for P: ', floor(hdur*c/dh)
       !c = sqrt(RIGIDITY/DENSITY)
       !print *, 'number of gridpoints per wavelength for S:', floor(hdur*c/dh)
       print *, 'Velocities in km/s :'
       print *, '  Pmin  :', sngl(alpha_min/1000.0)
       print *, '  Pmax  :', sngl(alpha_max/1000.0)
       print *, '  Smin  :', sngl(beta_min/1000.0)
       print *, '  Smax  :', sngl(beta_max/1000.0)
       print *, 'Gridpoints per wavelength :'
       print *, '  Pmin  :', floor(hdur*alpha_min/dh)
       print *, '  Pmax  :', floor(hdur*alpha_max/dh)
       print *, '  Smin  :', floor(hdur*beta_min/dh)
       print *, '  Smax  :', floor(hdur*beta_max/dh)
    endif

    NINT = NSTEP/NSAVE
    if (NINT * NSAVE > NSTEP) stop 'NSTEP should equal to NINT * NSAVE'

    ! calculate the global mass matrix once and for all
    ! note that the density variation is included here
    mass_global(:) = 0.0
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
    deltatover2 = deltat/2.0
    deltatsqover2 = deltat*deltat/2.0

    if (solver_type == 3) then
       b_deltat = -DT
       b_deltatover2 = b_deltat/2.0
       b_deltatsqover2 = b_deltat*b_deltat/2.0
    endif

    ! initialize
    displ(:,:) = 0.0
    veloc(:,:) = 0.0
    accel(:,:) = 0.0

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

       ! allocate 9 (3 x 3) kernels
       allocate(kappa_mu_rho_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(mu_kappa_rho_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(rho_kappa_mu_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(alpha_beta_rho_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(beta_alpha_rho_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(rho_alpha_beta_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(c_beta_rho_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(beta_c_rho_kernel(NGLLX,NGLLZ,NSPEC))
       allocate(rho_c_beta_kernel(NGLLX,NGLLZ,NSPEC))

       ! initialize kernels
       kappa_mu_rho_kernel = 0.0
       mu_kappa_rho_kernel = 0.0
       rho_kappa_mu_kernel = 0.0
       alpha_beta_rho_kernel = 0.0
       beta_alpha_rho_kernel = 0.0
       rho_alpha_beta_kernel = 0.0
       c_beta_rho_kernel = 0.0
       beta_c_rho_kernel = 0.0
       rho_c_beta_kernel = 0.0

       ! allocate interaction fields for 9 (3 x 3) kernels
       if (WRITE_KERNEL_SNAPSHOTS) then
          allocate(kappa_mu_rho_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(mu_kappa_rho_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(rho_kappa_mu_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(alpha_beta_rho_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(beta_alpha_rho_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(rho_alpha_beta_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(c_beta_rho_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(beta_c_rho_kernel_int(NGLLX,NGLLZ,NSPEC))
          allocate(rho_c_beta_kernel_int(NGLLX,NGLLZ,NSPEC))

          kappa_mu_rho_kernel_int = 0.0   ; mu_kappa_rho_kernel_int = 0.0   ; rho_kappa_mu_kernel_int = 0.0
          alpha_beta_rho_kernel_int = 0.0 ; beta_alpha_rho_kernel_int = 0.0 ; rho_alpha_beta_kernel_int = 0.0
          c_beta_rho_kernel_int = 0.0     ; beta_c_rho_kernel_int = 0.0     ; rho_c_beta_kernel_int = 0.0
       endif

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
          accel(i,:) = 0.0
          if (solver_type == 3) then
             b_displ(i,:) = b_displ(i,:) + b_deltat*b_veloc(i,:) + b_deltatsqover2*b_accel(i,:)
             b_veloc(i,:) = b_veloc(i,:) + b_deltatover2*b_accel(i,:)
             b_accel(i,:) = 0.0
          endif
       enddo

       ! displacement gradient at the source element GLL points at each time step
       if (solver_type == 3) displ_grad = 0.0

       if (NCOMP == 1) then    ! SH, or surface waves only

          !
          !   INTEGRATION OVER SPECTRAL ELEMENTS
          !
          do ispec = 1,NSPEC

             ! first double loop over GLL
             ! compute and store gradients
             do j = 1,NGLLZ
                do i = 1,NGLLX

                   !iglob2 = ibool(i,j,ispec)     ! 09-Jan-2007

                   ! derivative along x
                   tempy1l = 0.0
                   if (solver_type == 3) b_tempy1l = 0.0
                   do k = 1,NGLLX
                      hp1 = hprime_xx(k,i)
                      iglob = ibool(k,j,ispec)
                      tempy1l = tempy1l + displ(1,iglob)*hp1
                      if (solver_type == 3) then
                         b_tempy1l = b_tempy1l + b_displ(1,iglob)*hp1
                      endif
                   enddo

                   ! derivative along z
                   tempy2l = 0.0
                   if (solver_type == 3) b_tempy2l = 0.0
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
                   ! NOTE: point source only (nsrc=1)
                   if (solver_type == 3 .and. ispec == ispec_src(1)) then
                      displ_grad(i,j,1) = dsydxl
                      displ_grad(i,j,2) = dsydzl
                      !displ_grad(i,j,3) = dsydxl
                      !displ_grad(i,j,4) = dsydzl
                   endif

                   ! strain tensor
                   ds(:,:) = 0.0
                   ds(1,2) = oneovertwo * dsydxl
                   ds(2,3) = oneovertwo * dsydzl
                   ds(2,1) = ds(1,2)
                   ds(3,2) = ds(2,3)

                   if (solver_type == 3) then
                      b_dsydxl = b_tempy1l*dxidxl + b_tempy2l*dgammadxl
                      b_dsydzl = b_tempy1l*dxidzl + b_tempy2l*dgammadzl

                      b_ds(:,:) = 0.0
                      b_ds(1,2) = oneovertwo * b_dsydxl
                      b_ds(2,3) = oneovertwo * b_dsydzl
                      b_ds(2,1) = b_ds(1,2)
                      b_ds(3,2) = b_ds(2,3)

                      ! mu kernel : K-mu = double-dot-product of deviatoric strains
                      ! kappa kernel : not used for SH case (see kappa_k below)
                      ! (see below for kappa kernel)
                      mu_k(i,j,ispec) = sum(ds * b_ds)
                      !mu_k(iglob2) = sum(ds * b_ds)                                 ! (09-Jan-2007)
                      !mu_k(iglob2) = sum(ds * b_ds) - ONE_THIRD * kappa_k(iglob2)   ! (12-July-2006)
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
                   tempy1l = 0.0
                   if (solver_type == 3) b_tempy1l = 0.0
                   do k = 1,NGLLX
                      fac1 = wxgll(k)*hprime_xx(i,k)
                      tempy1l = tempy1l + tempy1(k,j)*fac1
                      if (solver_type == 3) then
                         b_tempy1l = b_tempy1l + b_tempy1(k,j)*fac1
                      endif
                   enddo

                   ! along z direction
                   tempy2l = 0.0
                   if (solver_type == 3) b_tempy2l = 0.0
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

       ! DEBUG ARRAY SIZES

          !
          !   INTEGRATION OVER SPECTRAL ELEMENTS
          !
          do ispec = 1,NSPEC

             ! first double loop over GLL
             ! compute and store gradients
             do j = 1,NGLLZ
                do i = 1,NGLLX

                   !iglob2 = ibool(i,j,ispec)

                   ! derivative along x
                   tempx1l = 0.0
                   tempy1l = 0.0
                   tempz1l = 0.0
                   if (solver_type == 3) then
                      b_tempx1l = 0.0
                      b_tempy1l = 0.0
                      b_tempz1l = 0.0
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
                   tempx2l = 0.0
                   tempy2l = 0.0
                   tempz2l = 0.0
                   if (solver_type == 3) then
                      b_tempx2l = 0.0
                      b_tempy2l = 0.0
                      b_tempz2l = 0.0
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

                   ! from the mesher
                   dxidxl = dxidx(i,j,ispec)
                   dxidzl = dxidz(i,j,ispec)
                   dgammadxl = dgammadx(i,j,ispec)
                   dgammadzl = dgammadz(i,j,ispec)

                   ! displacement gradients
                   dsxdxl = tempx1l*dxidxl+tempx2l*dgammadxl
                   dsxdzl = tempx1l*dxidzl+tempx2l*dgammadzl
                   dsydxl = tempy1l*dxidxl+tempy2l*dgammadxl
                   dsydzl = tempy1l*dxidzl+tempy2l*dgammadzl
                   dszdxl = tempz1l*dxidxl+tempz2l*dgammadxl
                   dszdzl = tempz1l*dxidzl+tempz2l*dgammadzl

                   ! save spatial gradient for (point) source perturbations
                   if (solver_type == 3 .and. ispec == ispec_src(1)) then
                      displ_grad(i,j,1) = dsxdxl
                      displ_grad(i,j,2) = dsxdzl
                      displ_grad(i,j,3) = dsydxl
                      displ_grad(i,j,4) = dsydzl
                      displ_grad(i,j,5) = dszdxl
                      displ_grad(i,j,6) = dszdzl
                   endif

                   ! strain tensor
                   ds(1,1) = dsxdxl
                   ds(1,2) = oneovertwo * dsydxl
                   ds(1,3) = oneovertwo * (dszdxl + dsxdzl)
                   ds(2,2) = 0                                ! 2D code
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
                      ! K-kappa = product of div_s
                      ! K-mu = double-dot-product of deviatoric strains
                      kappa_k(i,j,ispec) = (dsxdxl + dszdzl) * (b_dsxdxl + b_dszdzl)
                      mu_k(i,j,ispec)    = sum(ds * b_ds) - ONE_THIRD * kappa_k(i,j,ispec)
                      !kappa_k(iglob2) = (dsxdxl + dszdzl) * (b_dsxdxl + b_dszdzl)
                      !mu_k(iglob2)    = sum(ds * b_ds) - ONE_THIRD * kappa_k(iglob2)
                   endif

                   ! variation in kappa and mu
                   kappal = kappa(i,j,ispec)
                   mul = mu(i,j,ispec)
                   lambdalplus2mul = kappal + FOUR_THIRDS * mul
                   lambdal = lambdalplus2mul - 2.0*mul

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
                   tempx1l = 0.0
                   tempy1l = 0.0
                   tempz1l = 0.0
                   if (solver_type == 3) then
                      b_tempx1l = 0.0
                      b_tempy1l = 0.0
                      b_tempz1l = 0.0
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
                   tempx2l = 0.0
                   tempy2l = 0.0
                   tempz2l = 0.0
                   if (solver_type == 3) then
                      b_tempx2l = 0.0
                      b_tempy2l = 0.0
                      b_tempz2l = 0.0
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
                i = 1; nx = -1.0; nz = 0.0
             else if (ibb == 2) then
                i = NGLLX; nx = 1.0; nz = 0.0
             else if (ibb == 3) then
                i = 1; nx = 0.0; nz = -1.0
             else if (ibb == 4) then       ! CHT
                i = NGLLZ; nx = 0.0; nz = 1.0
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
                   !hlagrange = 0.0
                   !if (iglob == rglob(irec)) hlagrange = 1.0

                   ramp(itime,:,irec) = ramp(itime,:,irec) + displ(:,iglob)*hlagrange
                enddo
             enddo

             ! take the value at the closest gridpoint -- OLD METHOD
             !ramp(itime,:,irec) = displ(:,rglob(irec))

          enddo  ! irec

       else if (solver_type == 2 .or. solver_type == 3) then

          ! ADJOINT SEISMOGRAMS (adjoint wavefield recorded at the original source)
          ! These are needed for the source (and joint) inversions.
          ! note that samp above is copied to stf_for
          ! We fill the record 'in reverse' so that the max arrival will coincide
          ! with the forward source time function pulse.

          do isrc = 1,nsrc
             do icomp = 1,NCOMP

                ! NOTE: these values depend on the component
                temp1 = 0.0; temp2 = 0.0 ; temp3 = 0.0

                ! INTERPOLATE to obtain the function value at the exact source location;
                ! this requires looping over all GLL points within the element containing the source
                do j = 1,NGLLZ
                   do i = 1,NGLLX
                      ! only finds the GLL points of the element containing the source
                      iglob = ibool(i,j,ispec_src(isrc))

                      hlagrange = hxis_store(isrc,i) * hgammas_store(isrc,j)

                      ! KEY: interpolation to get exact source location values
                      !temp1 = temp1 + displ(icomp,iglob)*hlagrange    ! sydag
                      !temp2 = temp2 + displ_grad(i,j,3)*hlagrange     ! d/dx(sydag)
                      !temp3 = temp3 + displ_grad(i,j,4)*hlagrange     ! d/dz(sydag)
                      temp1 = temp1 + displ(icomp,iglob)*hlagrange           ! sxdag, sydag, szdag
                      temp2 = temp2 + displ_grad(i,j,2*icomp - 1)*hlagrange    ! d/dx(sdag)
                      temp3 = temp3 + displ_grad(i,j,2*icomp    )*hlagrange    ! d/dz(sdag)
                   enddo
                enddo

                !samp(NSTEP-itime+1,:,isrc) = displ(:,sglob(isrc))

                ! adjoint (displacement) seismograms (time-reversed)
                samp(NSTEP-itime+1,icomp,isrc) = temp1

                ! source perturbation TIME SERIES (3 x NCOMP time series)
                ! NOTE the ordering: dxs, dzs, dts
                three_source_model(NSTEP-itime+1,icomp,nsrc,1) = -temp1 * f0(icomp) * stf_syn_dt(NSTEP-itime+1) ! dts
                three_source_model(NSTEP-itime+1,icomp,nsrc,2) =  temp2 * f0(icomp) * stf_syn(NSTEP-itime+1)    ! dxs
                three_source_model(NSTEP-itime+1,icomp,nsrc,3) =  temp3 * f0(icomp) * stf_syn(NSTEP-itime+1)    ! dzs

                ! DEBUG ARRAY SIZES
                ! additional time series for checking
                if (0 == 1) then
                   three_source_model(NSTEP-itime+1,icomp,nsrc,4) = temp1
                   three_source_model(NSTEP-itime+1,icomp,nsrc,5) = temp2
                   three_source_model(NSTEP-itime+1,icomp,nsrc,6) = temp3
                   three_source_model(itime,icomp,nsrc,7) = stf_syn_dt(itime)
                   three_source_model(itime,icomp,nsrc,8) = stf_syn(itime)
                   three_source_model(itime,icomp,nsrc,9) = f0(icomp) * stf_syn(itime)
                endif

             enddo  ! icomp
          enddo  ! isrc

          if (solver_type == 3) then

             ! CALCULATE ALL NINE KERNELS
             !   note the TIME INTEGRATION
             !   note the minus sign in these expressions (Tromp-Tape-Liu, 2005)

             do ispec = 1,NSPEC
                do j = 1,NGLLZ
                   do i = 1,NGLLX
                      iglob = ibool(i,j,ispec)

                      ! kappa-mu-rho -- THESE ARE THE TWO KERNELS FOR AN ISOTROPIC ELASTIC MEDIUM, PLUS DENSITY
                      kappa_mu_rho_kernel(i,j,ispec) = kappa_mu_rho_kernel(i,j,ispec) - &
                           kappa(i,j,ispec) * kappa_k(i,j,ispec) * DT
                      mu_kappa_rho_kernel(i,j,ispec) = mu_kappa_rho_kernel(i,j,ispec) - &
                           2.0 * mu(i,j,ispec) * mu_k(i,j,ispec) * DT
                      rho_kappa_mu_kernel(i,j,ispec) = rho_kappa_mu_kernel(i,j,ispec) - &
                           rho(i,j,ispec) * dot_product(accel(:,iglob),b_displ(:,iglob)) * DT

                      ! alpha-beta-rho (derived from kappa-mu-rho kernels)
                      alpha_beta_rho_kernel(i,j,ispec) = 2.0 * (1.0 + FOUR_THIRDS*mu(i,j,ispec)/kappa(i,j,ispec)) &
                                                             * kappa_mu_rho_kernel(i,j,ispec)
                      beta_alpha_rho_kernel(i,j,ispec) = 2.0 * mu_kappa_rho_kernel(i,j,ispec) - &
                          2.0 * FOUR_THIRDS * mu(i,j,ispec)/kappa(i,j,ispec) * kappa_mu_rho_kernel(i,j,ispec)
                      rho_alpha_beta_kernel(i,j,ispec) = rho_kappa_mu_kernel(i,j,ispec) &
                                                       + kappa_mu_rho_kernel(i,j,ispec) &
                                                       + mu_kappa_rho_kernel(i,j,ispec)

                      ! c-beta-rho (derived from kappa-mu-rho kernels)
                      c_beta_rho_kernel(i,j,ispec) = 2.0 * kappa_mu_rho_kernel(i,j,ispec)
                      beta_c_rho_kernel(i,j,ispec) = 2.0 * mu_kappa_rho_kernel(i,j,ispec)
                      rho_c_beta_kernel(i,j,ispec) = rho_kappa_mu_kernel(i,j,ispec) &
                                                   + kappa_mu_rho_kernel(i,j,ispec) &
                                                   + mu_kappa_rho_kernel(i,j,ispec)

                      ! interaction fields -- each is integrated to form a kernel
                      if (WRITE_KERNEL_SNAPSHOTS) then
                         ! kappa-mu-rho
                         kappa_mu_rho_kernel_int(i,j,ispec) = -kappa(i,j,ispec) * kappa_k(i,j,ispec)
                         mu_kappa_rho_kernel_int(i,j,ispec) = -2.0*mu(i,j,ispec) * mu_k(i,j,ispec)
                         rho_kappa_mu_kernel_int(i,j,ispec) = -rho(i,j,ispec) * dot_product(accel(:,iglob),b_displ(:,iglob))

                         ! alpha-beta-rho (derived from kappa-mu-rho interaction fields)
                         alpha_beta_rho_kernel_int(i,j,ispec) = 2.0 * (1.0 + FOUR_THIRDS*mu(i,j,ispec)/kappa(i,j,ispec)) &
                                                                    * kappa_mu_rho_kernel_int(i,j,ispec)
                         beta_alpha_rho_kernel_int(i,j,ispec) = 2.0 * mu_kappa_rho_kernel_int(i,j,ispec) - &
                              2.0 * FOUR_THIRDS * mu(i,j,ispec)/kappa(i,j,ispec) * kappa_mu_rho_kernel_int(i,j,ispec)
                         rho_kappa_mu_kernel_int(i,j,ispec)   = kappa_mu_rho_kernel_int(i,j,ispec) + &
                              mu_kappa_rho_kernel_int(i,j,ispec) + rho_kappa_mu_kernel_int(i,j,ispec)

                         ! c-beta-rho (derived from kappa-mu-rho interaction fields)
                         c_beta_rho_kernel_int(i,j,ispec) = 2.0 * kappa_mu_rho_kernel_int(i,j,ispec)
                         beta_c_rho_kernel_int(i,j,ispec) = 2.0 * mu_kappa_rho_kernel_int(i,j,ispec)
                         rho_c_beta_kernel_int(i,j,ispec) = rho_kappa_mu_kernel_int(i,j,ispec) &
                                                          + kappa_mu_rho_kernel_int(i,j,ispec) &
                                                          + mu_kappa_rho_kernel_int(i,j,ispec)
                      endif  ! WRITE_KERNEL_SNAPSHOTS

                   enddo
                enddo
             enddo

          endif  ! solver_type = 3
       endif  ! solver_type = 2 or 3

       ! save FINAL SNAPSHOT
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
       endif

       !===================================

       ! write snapshots or kernels
       if (itime == 1 .or. mod(itime, NSAVE) == 0) then   ! save the initial frame
          if (itime == 1) tlab = 0
          if (itime /= 1) tlab = itime
          !if (mod(itime, NSAVE) == 0) then
          if (solver_type == 1) then
             if (idata == 0) write(filename1,'(a,i5.5)') trim(out_dir)//'forward_syn_',tlab
             if (idata == 1) write(filename1,'(a,i5.5)') trim(out_dir)//'forward_dat_',tlab
          else if (solver_type == 2) then
             write(filename1,'(a,i5.5)') trim(out_dir)//'adjoint_',tlab
          else
             write(filename1,'(a,i5.5)') trim(out_dir)//'adjoint_',tlab      ! adjoint wavefield
             write(filename2,'(a,i5.5)') trim(out_dir)//'backward_',tlab
             write(filename3,'(a,i5.5)') trim(out_dir)//'kernel_kmr_',tlab      ! kernels for kappa-mu-rho
             write(filename4,'(a,i5.5)') trim(out_dir)//'kernel_abr_',tlab      ! kernels for alpha-beta-rho
             write(filename5,'(a,i5.5)') trim(out_dir)//'kernel_cbr_',tlab      ! kernels for c-beta-rho
             !write(filename6,'(a,i5.5)') trim(out_dir)//'interaction_',tlab  ! interaction
          endif

          if (WRITE_WAVFIELD_SNAPSHOTS) then      ! wavefield snapshots (globally defined)
             open(unit=11, file=trim(filename1), status='unknown', iostat=ios)
             if (ios /= 0) stop 'Error writing snapshot to disk'
             do iglob = 1, NGLOB
                write(11,'(5e16.6)') sngl(x_plot(iglob)), sngl(z_plot(iglob)), (sngl(displ(j,iglob)),j=1,NCOMP)
             enddo
             close(11)
          endif

          if (WRITE_KERNEL_SNAPSHOTS) then        ! kernel snapshots (LOCALLY defined)
             if (solver_type == 3) then
                open(unit=11, file=trim(filename2), status='unknown', iostat=ios)
                if (ios /= 0) stop 'Error writing snapshot to disk'
                do iglob = 1, NGLOB
                   write(11,'(5e16.6)') sngl(x_plot(iglob)), sngl(z_plot(iglob)), (sngl(b_displ(j,iglob)),j=1,NCOMP)
                enddo
                close(11)

                open(unit = 11, file = trim(filename3), status = 'unknown',iostat=ios)
                open(unit = 12, file = trim(filename4), status = 'unknown',iostat=ios)
                open(unit = 13, file = trim(filename5), status = 'unknown',iostat=ios)
                !if (ios /= 0) stop 'Error writing snapshot to disk'

                do ispec = 1,NSPEC
                   do j = 1,NGLLZ
                      do i = 1,NGLLX
                         iglob = ibool(i,j,ispec)

                         ! kappa-mu-rho : three kernels and three interaction fields
                         write(11,'(8e16.6)') sngl(x_plot(iglob)), sngl(z_plot(iglob)), &
                              sngl(kappa_mu_rho_kernel(i,j,ispec)), &
                              sngl(mu_kappa_rho_kernel(i,j,ispec)), &
                              sngl(rho_kappa_mu_kernel(i,j,ispec)), &
                              sngl(kappa_mu_rho_kernel_int(i,j,ispec)), &
                              sngl(mu_kappa_rho_kernel_int(i,j,ispec)), &
                              sngl(rho_kappa_mu_kernel_int(i,j,ispec))

                         ! alpha-beta-rho : three kernels and three interaction fields
                         write(12,'(8e16.6)') sngl(x_plot(iglob)), sngl(z_plot(iglob)), &
                              sngl(alpha_beta_rho_kernel(i,j,ispec)), &
                              sngl(beta_alpha_rho_kernel(i,j,ispec)), &
                              sngl(rho_alpha_beta_kernel(i,j,ispec)), &
                              sngl(alpha_beta_rho_kernel_int(i,j,ispec)), &
                              sngl(beta_alpha_rho_kernel_int(i,j,ispec)), &
                              sngl(rho_alpha_beta_kernel_int(i,j,ispec))

                         ! c-beta-rho : three kernels and three interaction fields
                         write(13,'(8e16.6)') sngl(x_plot(iglob)), sngl(z_plot(iglob)), &
                              sngl(c_beta_rho_kernel(i,j,ispec)), &
                              sngl(beta_c_rho_kernel(i,j,ispec)), &
                              sngl(rho_c_beta_kernel(i,j,ispec)), &
                              sngl(c_beta_rho_kernel_int(i,j,ispec)), &
                              sngl(beta_c_rho_kernel_int(i,j,ispec)), &
                              sngl(rho_c_beta_kernel_int(i,j,ispec))
                   enddo
                enddo
             enddo
             close(11) ; close(12) ; close(13)

             endif  ! solver_type == 3
          endif  ! WRITE_KERNEL_SNAPSHOTS

       endif  ! write snapshots or kernels

       ! stop if wavefield contains NaN (13-Aug-2006)
       !do i=1,ncomp
       !   do j=1,nglob
       !      if ( isnan(displ(i,j)) ) stop 'exiting: encountered NaN in solver.f90'
       !   enddo
       !enddo

    enddo ! end time loop

    !-------------------------------------

    ! KERNELS at the FINAL TIME STEP
    if (solver_type == 3) then

       ! write out all nine kernels
       ! (Using local2global is a bit safer.)
       if (WRITE_KERNELS) then
          filename5 = trim(out_dir)//'kernel_basis'
          open(unit = 13, file = trim(filename5), status = 'unknown',iostat=ios)
          if (ios /= 0) stop 'Error writing all nine kernels to disk'
          do ispec = 1,NSPEC
             do j = 1,NGLLZ
                do i = 1,NGLLX
                   iglob = ibool(i,j,ispec)

                   ! x, z, -- kappa, mu, rho -- alpha, beta, rho -- c, beta rho
                   write(13,'(11e16.6)') sngl(x_plot(iglob)), sngl(z_plot(iglob)), &
                        sngl(kappa_mu_rho_kernel(i,j,ispec)), &
                        sngl(mu_kappa_rho_kernel(i,j,ispec)), &
                        sngl(rho_kappa_mu_kernel(i,j,ispec)), &
                        sngl(alpha_beta_rho_kernel(i,j,ispec)), &
                        sngl(beta_alpha_rho_kernel(i,j,ispec)), &
                        sngl(rho_alpha_beta_kernel(i,j,ispec)), &
                        sngl(c_beta_rho_kernel(i,j,ispec)), &
                        sngl(beta_c_rho_kernel(i,j,ispec)), &
                        sngl(rho_c_beta_kernel(i,j,ispec))
                enddo
             enddo
          enddo
          close(13)
       endif

       ! assign kernels for the inversion
       if (STRUCTURE_PARAMETER_TYPE == 1) then
          atype_kernel = kappa_mu_rho_kernel
          btype_kernel = mu_kappa_rho_kernel
          rtype_kernel = rho_kappa_mu_kernel

       else if (STRUCTURE_PARAMETER_TYPE == 2) then
          atype_kernel = alpha_beta_rho_kernel
          btype_kernel = beta_alpha_rho_kernel
          rtype_kernel = rho_alpha_beta_kernel

       else if (STRUCTURE_PARAMETER_TYPE == 3) then
          atype_kernel = c_beta_rho_kernel
          btype_kernel = beta_c_rho_kernel
          rtype_kernel = rho_c_beta_kernel
       endif

       ! deallocate variables for kernels
       deallocate(stf_for)
       deallocate(kappa_mu_rho_kernel, mu_kappa_rho_kernel, rho_kappa_mu_kernel)
       deallocate(alpha_beta_rho_kernel, beta_alpha_rho_kernel, rho_alpha_beta_kernel)
       deallocate(c_beta_rho_kernel, beta_c_rho_kernel, rho_c_beta_kernel)

       ! deallocate interaction fields
       if (WRITE_KERNEL_SNAPSHOTS) then
          deallocate(kappa_mu_rho_kernel_int, mu_kappa_rho_kernel_int, rho_kappa_mu_kernel_int)
          deallocate(alpha_beta_rho_kernel_int, beta_alpha_rho_kernel_int, rho_alpha_beta_kernel_int)
          deallocate(c_beta_rho_kernel_int, beta_c_rho_kernel_int, rho_c_beta_kernel_int)
       endif

    endif

  end subroutine solver
  !---------------------------------------------

end module wave2d_solver
