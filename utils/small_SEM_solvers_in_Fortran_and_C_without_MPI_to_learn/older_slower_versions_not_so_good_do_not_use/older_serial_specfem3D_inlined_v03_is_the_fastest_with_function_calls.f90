!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  program serial_specfem3D

  implicit none

!!!!!!!!!!
!!!!!!!!!! All the arrays below use static memory allocation,
!!!!!!!!!! using constant sizes defined in values_from_mesher.h.
!!!!!!!!!! This is done purposely to improve performance (Fortran compilers
!!!!!!!!!! can optimize much more when the size of the loops and arrays
!!!!!!!!!! is known at compile time).
!!!!!!!!!! NGLLX, NGLLY and NGLLZ are set equal to 5,
!!!!!!!!!! therefore each element contains NGLLX * NGLLY * NGLLZ = 125 points.
!!!!!!!!!!

!!!!!!!!!!
!!!!!!!!!! All the calculations are done in single precision.
!!!!!!!!!! We do not need double precision in SPECFEM3D.
!!!!!!!!!!

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "DATABASES_FOR_SOLVER/values_from_mesher_f90.h"

! constant value of the time step in the main time loop
  real(kind=4), parameter :: deltatover2 = 0.5*deltat, deltatsqover2 = 0.5*deltat*deltat

! for the source time function
  real, parameter :: pi = 3.141592653589793
  real, parameter :: f0 = 1. / 50.
  real, parameter :: t0 = 1.2 / f0
  real, parameter :: a = pi*pi*f0*f0

  integer, parameter :: NTSTEP_BETWEEN_OUTPUT_INFO = 200

  integer, parameter :: IIN = 40

! number of GLL integration points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

! 3-D simulation
  integer, parameter :: NDIM = 3

!!!!!!!!!!  real(kind=4), parameter :: VERYSMALLVAL = 1.e-24

! displacement threshold above which we consider that the code became unstable
  real(kind=4), parameter :: STABILITY_THRESHOLD = 1.e+25

! approximate density of the geophysical medium in which the source is located
! this value is only a constant scaling factor therefore it does not really matter
  real(kind=4), parameter :: rho = 4500.

! global displacement, velocity and acceleration vectors
  real(kind=4), dimension(NDIM,NGLOB) :: displ,veloc,accel

! global diagonal mass matrix
  real(kind=4), dimension(NGLOB) :: rmass_inverse

! record a seismogram to check that the simulation went well
  real(kind=4), dimension(NSTEP) :: seismogram

! time step
  integer it

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  real(kind=4), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,kappav,muv,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
!! DK DK store transpose of matrix
  real(kind=4), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
  real(kind=4), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=4), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=4), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=4), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  integer :: ispec,iglob,i,j,k !!!!!!!!!!!!! ,l

  real(kind=4) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=4) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=4) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=4) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=4) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
! real(kind=4) hp1,hp2,hp3
  real(kind=4) fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal
! real(kind=4) tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l

  real(kind=4) Usolidnorm,current_value,time,memory_size

  real(kind=4), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU

! estimate of total memory size used
  print *
  print *,'NSPEC = ',NSPEC
  print *,'NGLOB = ',NGLOB
  print *

  print *,'NSTEP = ',NSTEP
  print *,'deltat = ',deltat
  print *

! estimate total memory size (the size of a real number is 4 bytes)
! we perform the calculation in single precision rather than integer
! to avoid integer overflow in the case of very large meshes
  memory_size = 4. * ((3.*NDIM + 1.) * NGLOB + 12. * real(NGLLX*NGLLY*NGLLZ)*real(NSPEC))
  print *,'approximate total memory size used = ',memory_size/1024./1024.,' Mb'
  print *

! make sure the source element number is an integer
  if(mod(NSPEC,2) /= 0) stop 'source element number is not an integer, exiting...'

! read the mesh from external file
  open(unit=IIN,file='DATABASES_FOR_SOLVER/proc000000_reg1_database.dat',status='old')
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
! read real numbers here
          read(IIN,*) xix(i,j,k,ispec)
          read(IIN,*) xiy(i,j,k,ispec)
          read(IIN,*) xiz(i,j,k,ispec)
          read(IIN,*) etax(i,j,k,ispec)
          read(IIN,*) etay(i,j,k,ispec)
          read(IIN,*) etaz(i,j,k,ispec)
          read(IIN,*) gammax(i,j,k,ispec)
          read(IIN,*) gammay(i,j,k,ispec)
          read(IIN,*) gammaz(i,j,k,ispec)
          read(IIN,*) kappav(i,j,k,ispec)
          read(IIN,*) muv(i,j,k,ispec)

          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          jacobian(i,j,k,ispec) = 1. / (xixl*(etayl*gammazl-etazl*gammayl)- &
                xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl))


! read an integer here
          read(IIN,*) ibool(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
  do i = 1,NGLOB
    read(IIN,*) rmass_inverse(i)
! the real exactly diagonal mass matrix is read (not its inverse)
! therefore invert it here once and for all
    rmass_inverse(i) = 1. / rmass_inverse(i)
  enddo
  close(IIN)

  open(unit=IIN,file='DATABASES_FOR_SOLVER/matrices.dat',status='old')
  do j=1,NGLLY
    do i=1,NGLLX
      read(IIN,*) hprime_xx(i,j)
      read(IIN,*) hprimewgll_xx(i,j)
      read(IIN,*) wgllwgll_yz(i,j)
      read(IIN,*) wgllwgll_xz(i,j)
      read(IIN,*) wgllwgll_xy(i,j)
    enddo
  enddo
  close(IIN)

!! DK DK define transpose of matrix
  do j = 1,NGLLY
    do i = 1,NGLLX
      hprime_xxT(j,i) = hprime_xx(i,j)
      hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
    enddo
  enddo

  if(NGLLX /= 5) stop 'this inlined version with matrix products following Deville (2002) is only valid for NGLL = 5'

! clear initial vectors before starting the time loop
! (can remain serial because done only once before entering the time loop)
  displ(:,:) = 0. !!!!!!!!!! VERYSMALLVAL
  veloc(:,:) = 0.
  accel(:,:) = 0.

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! start of the time loop (which must remain serial obviously)
  do it = 1,NSTEP

! compute maximum of norm of displacement from time to time and display it
! in order to monitor the simulation
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      Usolidnorm = -1.
      do iglob = 1,NGLOB
        current_value = sqrt(displ(1,iglob)**2 + displ(2,iglob)**2 + displ(3,iglob)**2)
        if(current_value > Usolidnorm) Usolidnorm = current_value
      enddo
      write(*,*) 'Time step # ',it,' out of ',NSTEP
! compute current time
      time = (it-1)*deltat
      write(*,*) 'Max norm displacement vector U in the solid (m) = ',Usolidnorm
! check stability of the code, exit if unstable
      if(Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0) stop 'code became unstable and blew up'

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
             60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! elapsed time since beginning of the simulation
  tCPU = time_end - time_start
  int_tCPU = int(tCPU)
  ihours = int_tCPU / 3600
  iminutes = (int_tCPU - 3600*ihours) / 60
  iseconds = int_tCPU - 3600*ihours - 60*iminutes
  write(*,*) 'Elapsed time in seconds = ',tCPU
  write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
  write(*,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
  write(*,*)

    endif

! big loop over all the global points (not elements) in the mesh to update
! the displacement and velocity vectors and clear the acceleration vector
  displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
  accel(:,:) = 0.

! big loop over all the elements in the mesh to localize data
! from the global vectors to the local mesh
! using indirect addressing (contained in array ibool)
! and then to compute the elemental contribution
! to the acceleration vector of each element of the finite-element mesh
  do ispec = 1,NSPEC

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ(1,iglob)
            dummyy_loc(i,j,k) = displ(2,iglob)
            dummyz_loc(i,j,k) = displ(3,iglob)
        enddo
      enddo
    enddo

!! DK DK subroutines adapted from Deville, Fischer and Mund, High-order methods
!! DK DK for incompressible fluid flow, Cambridge University Press (2002),
!! DK DK pages 386 and 389 and Figure 8.3.1
  call mxm_m1_m2_5points(hprime_xx,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1)

  do k = 1,NGLLX
    call mxm_m1_m1_5points(dummyx_loc(1,1,k),dummyy_loc(1,1,k),dummyz_loc(1,1,k), &
           hprime_xxT,tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k))
  enddo

  call mxm_m2_m1_5points(dummyx_loc,dummyy_loc,dummyz_loc,hprime_xxT,tempx3,tempy3,tempz3)

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

!         tempx1l = 0.
!         tempx2l = 0.
!         tempx3l = 0.

!         tempy1l = 0.
!         tempy2l = 0.
!         tempy3l = 0.

!         tempz1l = 0.
!         tempz2l = 0.
!         tempz3l = 0.

!         do l=1,NGLLX
!           hp1 = hprime_xx(i,l)
!           tempx1l = tempx1l + dummyx_loc(l,j,k)*hp1
!           tempy1l = tempy1l + dummyy_loc(l,j,k)*hp1
!           tempz1l = tempz1l + dummyz_loc(l,j,k)*hp1

!           hp2 = hprime_xx(j,l)
!           tempx2l = tempx2l + dummyx_loc(i,l,k)*hp2
!           tempy2l = tempy2l + dummyy_loc(i,l,k)*hp2
!           tempz2l = tempz2l + dummyz_loc(i,l,k)*hp2

!           hp3 = hprime_xx(k,l)
!           tempx3l = tempx3l + dummyx_loc(i,j,l)*hp3
!           tempy3l = tempy3l + dummyy_loc(i,j,l)*hp3
!           tempz3l = tempz3l + dummyz_loc(i,j,l)*hp3
!         enddo

!         compute derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
!         jacobianl = 1. / (xixl*(etayl*gammazl-etazl*gammayl)- &
!               xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl))
!! DK DK this is now precomputed and stored to avoid a costly operation
          jacobianl = jacobian(i,j,k,ispec)

          duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
          duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
          duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

          duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
          duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
          duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

          duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
          duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
          duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! compute isotropic elements
          kappal = kappav(i,j,k,ispec)
          mul = muv(i,j,k,ispec)

!         lambdalplus2mul = kappal + (4./3.) * mul
!! DK DK precompute the 4/3 ratio to avoid a division here
          lambdalplus2mul = kappal + 1.33333333333333 * mul
          lambdal = lambdalplus2mul - 2.*mul

! compute stress sigma
          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

          sigma_xy = mul*duxdyl_plus_duydxl
          sigma_xz = mul*duzdxl_plus_duxdzl
          sigma_yz = mul*duzdyl_plus_duydzl

            ! define symmetric components of sigma
            sigma_yx = sigma_xy
            sigma_zx = sigma_xz
            sigma_zy = sigma_yz

            ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
            tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
            tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
            tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
            tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
            tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

            tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
            tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
            tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

          enddo
        enddo
      enddo

!! DK DK subroutines adapted from Deville, Fischer and Mund, High-order methods
!! DK DK for incompressible fluid flow, Cambridge University Press (2002),
!! DK DK pages 386 and 389 and Figure 8.3.1
  call mxm_m1_m2_5points(hprimewgll_xxT,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1)

  do k = 1,NGLLX
    call mxm_m1_m1_5points(tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k), &
          hprimewgll_xx,newtempx2(1,1,k),newtempy2(1,1,k),newtempz2(1,1,k))
  enddo

  call mxm_m2_m1_5points(tempx3,tempy3,tempz3,hprimewgll_xx,newtempx3,newtempy3,newtempz3)

!   do k=1,NGLLZ
!     do j=1,NGLLY
!       do i=1,NGLLX

!         tempx1l = 0.
!         tempy1l = 0.
!         tempz1l = 0.

!         tempx2l = 0.
!         tempy2l = 0.
!         tempz2l = 0.

!         tempx3l = 0.
!         tempy3l = 0.
!         tempz3l = 0.

!         do l=1,NGLLX
!           fac1 = hprimewgll_xx(l,i)
!           tempx1l = tempx1l + tempx1(l,j,k)*fac1
!           tempy1l = tempy1l + tempy1(l,j,k)*fac1
!           tempz1l = tempz1l + tempz1(l,j,k)*fac1

!           fac2 = hprimewgll_xx(l,j)
!           tempx2l = tempx2l + tempx2(i,l,k)*fac2
!           tempy2l = tempy2l + tempy2(i,l,k)*fac2
!           tempz2l = tempz2l + tempz2(i,l,k)*fac2

!           fac3 = hprimewgll_xx(l,k)
!           tempx3l = tempx3l + tempx3(i,j,l)*fac3
!           tempy3l = tempy3l + tempy3(i,j,l)*fac3
!           tempz3l = tempz3l + tempz3(i,j,l)*fac3
!         enddo

!       enddo
!     enddo
!   enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

! sum contributions from each element to the global mesh using indirect addressing
          iglob = ibool(i,j,k,ispec)
          accel(1,iglob) = accel(1,iglob) - (fac1*newtempx1(i,j,k) + fac2*newtempx2(i,j,k) + fac3*newtempx3(i,j,k))
          accel(2,iglob) = accel(2,iglob) - (fac1*newtempy1(i,j,k) + fac2*newtempy2(i,j,k) + fac3*newtempy3(i,j,k))
          accel(3,iglob) = accel(3,iglob) - (fac1*newtempz1(i,j,k) + fac2*newtempz2(i,j,k) + fac3*newtempz3(i,j,k))

        enddo
      enddo
    enddo

  enddo   ! end of main loop on all the elements

! big loop over all the global points (not elements) in the mesh to update
! the acceleration and velocity vectors
    accel(1,:) = accel(1,:)*rmass_inverse(:)
    accel(2,:) = accel(2,:)*rmass_inverse(:)
    accel(3,:) = accel(3,:)*rmass_inverse(:)

! add the earthquake source at a given grid point
! this is negligible and can remain serial because it is done by only
! one grid point out of several millions typically
    iglob = ibool(2,2,2,NSPEC_SOURCE)
! compute current time
    time = (it-1)*deltat
    accel(3,iglob) = accel(3,iglob) + 1.e4 * (1.-2.*a*(time-t0)**2) * exp(-a*(time-t0)**2) / rho

    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)

! record a seismogram to check that the simulation went well
    seismogram(it) = displ(3,ibool(2,2,2,NSPEC_STATION))

  enddo ! end of the serial time loop

! save the seismogram at the end of the run
  open(unit=IIN,file='seismogram_F90.txt',status='unknown')
  do it = 1,NSTEP
    write(IIN,*) (it-1)*deltat,seismogram(it)
  enddo
  close(IIN)

  end program serial_specfem3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! DK DK subroutines adapted from Deville, Fischer and Mund, High-order methods
!! DK DK for incompressible fluid flow, Cambridge University Press (2002),
!! DK DK pages 386 and 389 and Figure 8.3.1

  subroutine mxm_m1_m2_5points(A,B1,B2,B3,C1,C2,C3)

  implicit none

  include "constants.h"

  real(kind=4), dimension(m1,NGLLX) :: A
  real(kind=4), dimension(NGLLX,m2) :: B1,B2,B3
  real(kind=4), dimension(m1,m2) :: C1,C2,C3

  integer :: i,j

  do j=1,m2
    do i=1,m1

      C1(i,j) = A(i,1)*B1(1,j) + &
                A(i,2)*B1(2,j) + &
                A(i,3)*B1(3,j) + &
                A(i,4)*B1(4,j) + &
                A(i,5)*B1(5,j)

      C2(i,j) = A(i,1)*B2(1,j) + &
                A(i,2)*B2(2,j) + &
                A(i,3)*B2(3,j) + &
                A(i,4)*B2(4,j) + &
                A(i,5)*B2(5,j)

      C3(i,j) = A(i,1)*B3(1,j) + &
                A(i,2)*B3(2,j) + &
                A(i,3)*B3(3,j) + &
                A(i,4)*B3(4,j) + &
                A(i,5)*B3(5,j)

    enddo
  enddo

  end subroutine mxm_m1_m2_5points

!---------

  subroutine mxm_m1_m1_5points(A1,A2,A3,B,C1,C2,C3)

  implicit none

  include "constants.h"

  real(kind=4), dimension(m1,NGLLX) :: A1,A2,A3
  real(kind=4), dimension(NGLLX,m1) :: B
  real(kind=4), dimension(m1,m1) :: C1,C2,C3

  integer :: i,j

  do j=1,m1
    do i=1,m1

      C1(i,j) = A1(i,1)*B(1,j) + &
                A1(i,2)*B(2,j) + &
                A1(i,3)*B(3,j) + &
                A1(i,4)*B(4,j) + &
                A1(i,5)*B(5,j)

      C2(i,j) = A2(i,1)*B(1,j) + &
                A2(i,2)*B(2,j) + &
                A2(i,3)*B(3,j) + &
                A2(i,4)*B(4,j) + &
                A2(i,5)*B(5,j)

      C3(i,j) = A3(i,1)*B(1,j) + &
                A3(i,2)*B(2,j) + &
                A3(i,3)*B(3,j) + &
                A3(i,4)*B(4,j) + &
                A3(i,5)*B(5,j)

    enddo
  enddo

  end subroutine mxm_m1_m1_5points

!---------

  subroutine mxm_m2_m1_5points(A1,A2,A3,B,C1,C2,C3)

  implicit none

  include "constants.h"

  real(kind=4), dimension(m2,NGLLX) :: A1,A2,A3
  real(kind=4), dimension(NGLLX,m1) :: B
  real(kind=4), dimension(m2,m1) :: C1,C2,C3

  integer :: i,j

  do j=1,m1
    do i=1,m2

      C1(i,j) = A1(i,1)*B(1,j) + &
                A1(i,2)*B(2,j) + &
                A1(i,3)*B(3,j) + &
                A1(i,4)*B(4,j) + &
                A1(i,5)*B(5,j)

      C2(i,j) = A2(i,1)*B(1,j) + &
                A2(i,2)*B(2,j) + &
                A2(i,3)*B(3,j) + &
                A2(i,4)*B(4,j) + &
                A2(i,5)*B(5,j)

      C3(i,j) = A3(i,1)*B(1,j) + &
                A3(i,2)*B(2,j) + &
                A3(i,3)*B(3,j) + &
                A3(i,4)*B(4,j) + &
                A3(i,5)*B(5,j)

    enddo
  enddo

  end subroutine mxm_m2_m1_5points

