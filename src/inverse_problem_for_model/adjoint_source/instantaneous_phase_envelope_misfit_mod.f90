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


! this file is not used yet...

module instantaneous_phase_envelope_misfit_mod

  use domain_decomposition_mod
  use interpolation_mod
  use cost_function_mod
  use acquisition_mod
  use parameterization_mod, only: npar0, inv_param
  use common_engine_mod
  use common_fwi_mod
  use filter_data_mod

  implicit none

  real(kind=cp),    parameter :: ipi   = 3.14159265359
  complex(kind=cp), parameter :: icmpl = (0.,1.)

  complex(kind=cp) :: tmp_four

  complex(kind=cp), dimension(:,:), allocatable :: ft_dobs_vx, ft_dobs_vy, ft_dobs_vz
  complex(kind=cp), dimension(:,:), allocatable :: ft_dcal_vx, ft_dcal_vy, ft_dcal_vz

  complex(kind=cp), dimension(:,:), allocatable :: danalytic_vx, danalytic_vy, danalytic_vz
  complex(kind=cp), dimension(:,:), allocatable :: an_dobs_vx, an_dobs_vy, an_dobs_vz
  complex(kind=cp), dimension(:,:), allocatable :: an_dcal_vx, an_dcal_vy, an_dcal_vz

  real(kind=cp), dimension(:,:), allocatable :: dphase_vx, dphase_vy, dphase_vz
  real(kind=cp), dimension(:,:), allocatable :: phaseo_vx, phaseo_vy, phaseo_vz
  real(kind=cp), dimension(:,:), allocatable :: phasec_vx, phasec_vy, phasec_vz

  real(kind=cp), dimension(:,:), allocatable :: denvel_vx, denvel_vy, denvel_vz
  real(kind=cp), dimension(:,:), allocatable :: envelo_vx, envelo_vy, envelo_vz
  real(kind=cp), dimension(:,:), allocatable :: envelc_vx, envelc_vy, envelc_vz

  real(kind=cp), dimension(:,:), allocatable :: IP_adjt_vx, IP_adjt_vy, IP_adjt_vz
  real(kind=cp), dimension(:,:), allocatable :: EN_adjt_vx, EN_adjt_vy, EN_adjt_vz

  real(kind=cp), dimension(:), allocatable :: hh

  real(kind=cp) :: wlx, wly, wlz, amp_wlx, amp_wly, amp_wlz

contains

!================================================================================
! Compute instantenous phase data
  subroutine compute_instantaneous_phase_and_envelope_data(giter)

    integer(kind=si), intent(in) :: giter

    !*** 1. filter data obs and cal
    !******* in case this is not done before...
    !*** Loop over receivers
!!$    call prepare_adjoint_filter
!!$    if (nrecloc >0) then
!!$       do irec=1,nrecloc
!!$
!!$          !*** Taper on residuals
!!$          dcal_vx(irec,:) = dcal_vx(irec,:) * taper(:)
!!$          dcal_vy(irec,:) = dcal_vy(irec,:) * taper(:)
!!$          dcal_vz(irec,:) = dcal_vz(irec,:) * taper(:)
!!$
!!$          !*** Filter on residuals
!!$          call bwfilt(dcal_vx(irec,:),conv,dt,nt,1,4,1e-3_cp,inv_freq)
!!$          dcal_vx(irec,:) = conv(:)
!!$          call bwfilt(dcal_vy(irec,:),conv,dt,nt,1,4,1e-3_cp,inv_freq)
!!$          dcal_vy(irec,:) = conv(:)
!!$          call bwfilt(dcal_vz(irec,:),conv,dt,nt,1,4,1e-3_cp,inv_freq)
!!$          dcal_vz(irec,:) = conv(:)
!!$
!!$          !*** Taper on residuals
!!$          dobs_vx(irec,:) = dobs_vx(irec,:) * taper(:)
!!$          dobs_vy(irec,:) = dobs_vy(irec,:) * taper(:)
!!$          dobs_vz(irec,:) = dobs_vz(irec,:) * taper(:)
!!$
!!$          !*** Filter on residuals
!!$          call bwfilt(dobs_vx(irec,:),conv,dt,nt,1,4,1e-3_cp,inv_freq)
!!$          dobs_vx(irec,:) = conv(:)
!!$          call bwfilt(dobs_vy(irec,:),conv,dt,nt,1,4,1e-3_cp,inv_freq)
!!$          dobs_vy(irec,:) = conv(:)
!!$          call bwfilt(dobs_vz(irec,:),conv,dt,nt,1,4,1e-3_cp,inv_freq)
!!$          dobs_vz(irec,:) = conv(:)
!!$
!!$       enddo
!!$    endif

    !*** Absolute value for water level
    if (giter == 0) then
       amp_wlx = 0.
       amp_wly = 0.
       amp_wlz = 0.

       !* dobs
       if (nrecloc > 0) then
          !amp_wlx = maxval(abs(dcal_vx))
          !amp_wly = maxval(abs(dcal_vy))
          amp_wlz = maxval(abs(dobs_vz))

       endif

!      call MPI_allreduce(MPI_IN_PLACE,amp_wlx,1,MPI_REAL,MPI_SUM,comm%mpi_comm_1,ierr_mpi)
!      call MPI_allreduce(MPI_IN_PLACE,amp_wly,1,MPI_REAL,MPI_SUM,comm%mpi_comm_1,ierr_mpi)
       call MPI_allreduce(MPI_IN_PLACE,amp_wlz,1,MPI_REAL,MPI_SUM,comm%mpi_comm_1,ierr_mpi)

    endif

    if (nrecloc > 0) then

       if (.not. allocated(envelo_vx)) then
         allocate(envelo_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 221')
       endif
       if (.not. allocated(envelo_vy)) then
         allocate(envelo_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 222')
       endif
       if (.not. allocated(envelo_vz)) then
         allocate(envelo_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 223')
       endif

       if (.not. allocated(envelc_vx)) then
         allocate(envelc_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 224')
       endif
       if (.not. allocated(envelc_vy)) then
         allocate(envelc_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 225')
       endif
       if (.not. allocated(envelc_vz)) then
         allocate(envelc_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 226')
       endif

       if (.not. allocated(denvel_vx)) then
         allocate(denvel_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 227')
       endif
       if (.not. allocated(denvel_vy)) then
         allocate(denvel_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 228')
       endif
       if (.not. allocated(denvel_vz)) then
         allocate(denvel_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 229')
       endif

       if (.not. allocated(dphase_vx)) then
         allocate(dphase_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 230')
       endif
       if (.not. allocated(dphase_vy)) then
         allocate(dphase_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 231')
       endif
       if (.not. allocated(dphase_vz)) then
         allocate(dphase_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 232')
       endif

       if (.not. allocated(danalytic_vx)) then
         allocate(danalytic_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 233')
       endif
       if (.not. allocated(danalytic_vy)) then
         allocate(danalytic_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 234')
       endif
       if (.not. allocated(danalytic_vz)) then
         allocate(danalytic_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 235')
       endif

       if (.not. allocated(an_dobs_vx)) then
         allocate(an_dobs_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 236')
       endif
       if (.not. allocated(an_dobs_vy)) then
         allocate(an_dobs_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 237')
       endif
       if (.not. allocated(an_dobs_vz)) then
         allocate(an_dobs_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 238')
       endif

       if (.not. allocated(an_dcal_vx)) then
         allocate(an_dcal_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 239')
       endif
       if (.not. allocated(an_dcal_vy)) then
         allocate(an_dcal_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 240')
       endif
       if (.not. allocated(an_dcal_vz)) then
         allocate(an_dcal_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 241')
       endif

       !*** 2. Get their analytic signal
       call get_analytic_signal

       !*** 3. Compute envelope
       envelo_vx  = real(sqrt(conjg(an_dobs_vx) * an_dobs_vx))
       envelc_vx  = real(sqrt(conjg(an_dcal_vx) * an_dcal_vx))
       envelo_vy  = real(sqrt(conjg(an_dobs_vy) * an_dobs_vy))
       envelc_vy  = real(sqrt(conjg(an_dcal_vy) * an_dcal_vy))
       envelo_vz  = real(sqrt(conjg(an_dobs_vz) * an_dobs_vz))
       envelc_vz  = real(sqrt(conjg(an_dcal_vz) * an_dcal_vz))

       !*** 4. Define water-level with user percentage and amplitude of dcal for current
       !    source (and why not component)
       wlx = wl_perc * amp_wlz
       wly = wl_perc * amp_wlz
       wlz = wl_perc * amp_wlz

       !*** 5. Compute analytic signal misfit
       danalytic_vx = conjg(an_dobs_vx) * an_dcal_vx / (envelc_vx * envelo_vx + wlx**2)
       danalytic_vy = conjg(an_dobs_vy) * an_dcal_vy / (envelc_vy * envelo_vy + wly**2)
       danalytic_vz = conjg(an_dobs_vz) * an_dcal_vz / (envelc_vz * envelo_vz + wlz**2)

       !*** 6. Deduce phase misfit and envelope misfit
       dphase_vx = asin(aimag(danalytic_vx))
       dphase_vy = asin(aimag(danalytic_vy))
       dphase_vz = asin(aimag(danalytic_vz))
       denvel_vx = log(envelc_vx / (envelo_vx + wlx))
       denvel_vy = log(envelc_vy / (envelo_vy + wly))
       denvel_vz = log(envelc_vz / (envelo_vz + wlz))

       if (misfit_type == 3) then
          residu_vx = dphase_vx
          residu_vy = dphase_vy
          residu_vz = dphase_vz
       else if (misfit_type == 4) then
          residu_vx = denvel_vx
          residu_vy = denvel_vy
          residu_vz = denvel_vz
       endif

       !call filter_adjoint_source

    endif

  end subroutine compute_instantaneous_phase_and_envelope_data
!--------------------------------------------------------------------------------

!================================================================================
! Compute instantaneous phase source term
  subroutine compute_instanteneous_phase_adjoint_source_term

    integer(kind=si) :: ff, irec, tt

    real(kind=cp), dimension(:,:), allocatable :: srcterm1_vx, srcterm1_vy, srcterm1_vz
    real(kind=cp), dimension(:,:), allocatable :: srcterm2tmp_vx, srcterm2tmp_vy, srcterm2tmp_vz
    complex(kind=cp), dimension(:,:), allocatable :: srcterm2_vx, srcterm2_vy, srcterm2_vz
    complex(kind=cp), dimension(:,:), allocatable :: ft_tmp_vx, ft_tmp_vy, ft_tmp_vz

    if (nrecloc > 0) then

       if (.not. allocated(srcterm1_vx)) then
         allocate(srcterm1_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 242')
       endif
       if (.not. allocated(srcterm1_vy)) then
         allocate(srcterm1_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 243')
       endif
       if (.not. allocated(srcterm1_vz)) then
         allocate(srcterm1_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 244')
       endif

       if (.not. allocated(srcterm2_vx)) then
         allocate(srcterm2_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 245')
       endif
       if (.not. allocated(srcterm2_vy)) then
         allocate(srcterm2_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 246')
       endif
       if (.not. allocated(srcterm2_vz)) then
         allocate(srcterm2_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 247')
       endif

       if (.not. allocated(srcterm2tmp_vx)) then
         allocate(srcterm2tmp_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 248')
       endif
       if (.not. allocated(srcterm2tmp_vy)) then
         allocate(srcterm2tmp_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 249')
       endif
       if (.not. allocated(srcterm2tmp_vz)) then
         allocate(srcterm2tmp_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 250')
       endif

       if (.not. allocated(ft_tmp_vx)) then
         allocate(ft_tmp_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 251')
       endif
       if (.not. allocated(ft_tmp_vy)) then
         allocate(ft_tmp_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 252')
       endif
       if (.not. allocated(ft_tmp_vz)) then
         allocate(ft_tmp_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 253')
       endif

       if (.not. allocated(IP_adjt_vx)) then
         allocate(IP_adjt_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 254')
       endif
       if (.not. allocated(IP_adjt_vy)) then
         allocate(IP_adjt_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 255')
       endif
       if (.not. allocated(IP_adjt_vz)) then
         allocate(IP_adjt_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 256')
       endif

       !*** 1st term (easy...) remeber, imag(analytic_sig) = hilbert transform
       srcterm1_vx = dphase_vx * aimag(danalytic_vx) / (envelc_vx**2 + wlx**2)
       srcterm1_vy = dphase_vy * aimag(danalytic_vy) / (envelc_vy**2 + wly**2)
       srcterm1_vz = dphase_vz * aimag(danalytic_vz) / (envelc_vz**2 + wlz**2)

       !*** 2nd term (harder...)
       !* 1. Compute tmp signal of 2nd term
       srcterm2tmp_vx = dphase_vx * real(danalytic_vx) / (envelc_vx**2 + wlx**2)
       srcterm2tmp_vy = dphase_vy * real(danalytic_vy) / (envelc_vy**2 + wly**2)
       srcterm2tmp_vz = dphase_vz * real(danalytic_vz) / (envelc_vz**2 + wlz**2)

       !* 2. Take its hilbert transform
       !* DFT
       ft_tmp_vx = cmplx(0.,0.)
       ft_tmp_vy = cmplx(0.,0.)
       ft_tmp_vz = cmplx(0.,0.)
       do ff=1,nt
          do tt=1,nt

             !* Fourier basis
             tmp_four = cos(-2. * ipi * real(ff-1.) * real(tt-1.) / nt) + &
                icmpl * sin(-2. * ipi * real(ff-1.) * real(tt-1.) / nt)

             !* Data
             ft_tmp_vx(:,ff) = ft_tmp_vx(:,ff) + srcterm2tmp_vx(:,tt) * tmp_four
             ft_tmp_vy(:,ff) = ft_tmp_vy(:,ff) + srcterm2tmp_vy(:,tt) * tmp_four
             ft_tmp_vz(:,ff) = ft_tmp_vz(:,ff) + srcterm2tmp_vz(:,tt) * tmp_four

          enddo
       enddo

       !* Remove negative freq
       do irec=1,nrecloc
          ft_tmp_vx(irec,:) = ft_tmp_vx(irec,:) * hh(:)
          ft_tmp_vy(irec,:) = ft_tmp_vy(irec,:) * hh(:)
          ft_tmp_vz(irec,:) = ft_tmp_vz(irec,:) * hh(:)
       enddo

       !* Perform inverse DFT
       srcterm2_vx = cmplx(0.,0.)
       srcterm2_vy = cmplx(0.,0.)
       srcterm2_vz = cmplx(0.,0.)
       do tt=1,nt
          do ff=1,nt

             !* Fourier basis
             tmp_four = cos(2. * ipi * real(ff-1.) * real(tt-1.) / nt) + &
                icmpl * sin(2. * ipi * real(ff-1.) * real(tt-1.) / nt)

             !* Data
             srcterm2_vx(:,ff) = srcterm2_vx(:,ff) + ft_tmp_vx(:,tt) * tmp_four
             srcterm2_vy(:,ff) = srcterm2_vy(:,ff) + ft_tmp_vy(:,tt) * tmp_four
             srcterm2_vz(:,ff) = srcterm2_vz(:,ff) + ft_tmp_vz(:,tt) * tmp_four

          enddo
       enddo

       !* Renormalize dft
       srcterm2_vx(:,:) = srcterm2_vx(:,:) / nt
       srcterm2_vy(:,:) = srcterm2_vy(:,:) / nt
       srcterm2_vz(:,:) = srcterm2_vz(:,:) / nt

       !*** Deduce adjoint source term (again imag(analytic signal) = hilbert)
       IP_adjt_vx = srcterm1_vx + aimag(srcterm2_vx)
       IP_adjt_vy = srcterm1_vy + aimag(srcterm2_vy)
       IP_adjt_vz = srcterm1_vz + aimag(srcterm2_vz)


       if (allocated(srcterm1_vx)) deallocate(srcterm1_vx)
       if (allocated(srcterm1_vy)) deallocate(srcterm1_vy)
       if (allocated(srcterm1_vz)) deallocate(srcterm1_vz)

       if (allocated(srcterm2tmp_vx)) deallocate(srcterm2tmp_vx)
       if (allocated(srcterm2tmp_vy)) deallocate(srcterm2tmp_vy)
       if (allocated(srcterm2tmp_vz)) deallocate(srcterm2tmp_vz)

       if (allocated(srcterm2_vx)) deallocate(srcterm2_vx)
       if (allocated(srcterm2_vy)) deallocate(srcterm2_vy)
       if (allocated(srcterm2_vz)) deallocate(srcterm2_vz)

       if (allocated(ft_tmp_vx)) deallocate(ft_tmp_vx)
       if (allocated(ft_tmp_vy)) deallocate(ft_tmp_vy)
       if (allocated(ft_tmp_vz)) deallocate(ft_tmp_vz)

       residu_vx = -IP_adjt_vx
       residu_vy = -IP_adjt_vy
       residu_vz = -IP_adjt_vz

    endif

  end subroutine compute_instanteneous_phase_adjoint_source_term
!--------------------------------------------------------------------------------


!================================================================================
! Compute envelope source term
  subroutine compute_envelope_adjoint_source_term

    integer(kind=si) :: ff, irec, tt

    real(kind=cp), dimension(:,:), allocatable :: srcterm1_vx, srcterm1_vy, srcterm1_vz
    real(kind=cp), dimension(:,:), allocatable :: srcterm2tmp_vx, srcterm2tmp_vy, srcterm2tmp_vz
    complex(kind=cp), dimension(:,:), allocatable :: srcterm2_vx, srcterm2_vy, srcterm2_vz
    complex(kind=cp), dimension(:,:), allocatable :: ft_tmp_vx, ft_tmp_vy, ft_tmp_vz

    if (nrecloc > 0) then

       if (.not. allocated(srcterm1_vx)) then
         allocate(srcterm1_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 257')
       endif
       if (.not. allocated(srcterm1_vy)) then
         allocate(srcterm1_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 258')
       endif
       if (.not. allocated(srcterm1_vz)) then
         allocate(srcterm1_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 259')
       endif

       if (.not. allocated(srcterm2_vx)) then
         allocate(srcterm2_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 260')
       endif
       if (.not. allocated(srcterm2_vy)) then
         allocate(srcterm2_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 261')
       endif
       if (.not. allocated(srcterm2_vz)) then
         allocate(srcterm2_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 262')
       endif

       if (.not. allocated(srcterm2tmp_vx)) then
         allocate(srcterm2tmp_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 263')
       endif
       if (.not. allocated(srcterm2tmp_vy)) then
         allocate(srcterm2tmp_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 264')
       endif
       if (.not. allocated(srcterm2tmp_vz)) then
         allocate(srcterm2tmp_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 265')
       endif

       if (.not. allocated(ft_tmp_vx)) then
         allocate(ft_tmp_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 266')
       endif
       if (.not. allocated(ft_tmp_vy)) then
         allocate(ft_tmp_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 267')
       endif
       if (.not. allocated(ft_tmp_vz)) then
         allocate(ft_tmp_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 268')
       endif

       if (.not. allocated(EN_adjt_vx)) then
         allocate(EN_adjt_vx(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 269')
       endif
       if (.not. allocated(EN_adjt_vy)) then
         allocate(EN_adjt_vy(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 270')
       endif
       if (.not. allocated(EN_adjt_vz)) then
         allocate(EN_adjt_vz(nrecloc,nt),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 271')
       endif

       !*** 1st term (easy...) remeber, imag(analytic_sig) = hilbert transform
       srcterm1_vx = denvel_vx * real(danalytic_vx) / (envelc_vx**2 + wlx**2)
       srcterm1_vy = denvel_vy * real(danalytic_vy) / (envelc_vy**2 + wly**2)
       srcterm1_vz = denvel_vz * real(danalytic_vz) / (envelc_vz**2 + wlz**2)

       !*** 2nd term (harder...)
       !* 1. Compute tmp signal of 2nd term
       srcterm2tmp_vx = denvel_vx * aimag(danalytic_vx) / (envelc_vx**2 + wlx**2)
       srcterm2tmp_vy = denvel_vy * aimag(danalytic_vy) / (envelc_vy**2 + wly**2)
       srcterm2tmp_vz = denvel_vz * aimag(danalytic_vz) / (envelc_vz**2 + wlz**2)

       !* 2. Take its hilbert transform
       !* DFT
       ft_tmp_vx = cmplx(0.,0.)
       ft_tmp_vy = cmplx(0.,0.)
       ft_tmp_vz = cmplx(0.,0.)
       do ff=1,nt
          do tt=1,nt

             !* Fourier basis
             tmp_four = cos(-2. * ipi * real(ff-1.) * real(tt-1.) / nt) + &
                icmpl * sin(-2. * ipi * real(ff-1.) * real(tt-1.) / nt)

             !* Data
             ft_tmp_vx(:,ff) = ft_tmp_vx(:,ff) + srcterm2tmp_vx(:,tt) * tmp_four
             ft_tmp_vy(:,ff) = ft_tmp_vy(:,ff) + srcterm2tmp_vy(:,tt) * tmp_four
             ft_tmp_vz(:,ff) = ft_tmp_vz(:,ff) + srcterm2tmp_vz(:,tt) * tmp_four

          enddo
       enddo

       !* Remove negative freq
       do irec=1,nrecloc
          ft_tmp_vx(irec,:) = ft_tmp_vx(irec,:) * hh(:)
          ft_tmp_vy(irec,:) = ft_tmp_vy(irec,:) * hh(:)
          ft_tmp_vz(irec,:) = ft_tmp_vz(irec,:) * hh(:)
       enddo

       !* Perform inverse DFT
       srcterm2_vx = cmplx(0.,0.)
       srcterm2_vy = cmplx(0.,0.)
       srcterm2_vz = cmplx(0.,0.)
       do tt=1,nt
          do ff=1,nt

             !* Fourier basis
             tmp_four = cos(2. * ipi * real(ff-1.) * real(tt-1.) / nt) + &
                icmpl * sin(2. * ipi * real(ff-1.) * real(tt-1.) / nt)

             !* Data
             srcterm2_vx(:,ff) = srcterm2_vx(:,ff) + ft_tmp_vx(:,tt) * tmp_four
             srcterm2_vy(:,ff) = srcterm2_vy(:,ff) + ft_tmp_vy(:,tt) * tmp_four
             srcterm2_vz(:,ff) = srcterm2_vz(:,ff) + ft_tmp_vz(:,tt) * tmp_four

          enddo
       enddo

       !* Renormalize dft
       srcterm2_vx(:,:) = srcterm2_vx(:,:) / nt
       srcterm2_vy(:,:) = srcterm2_vy(:,:) / nt
       srcterm2_vz(:,:) = srcterm2_vz(:,:) / nt

       !*** Deduce adjoint source term (again imag(analytic signal) = hilbert)
       EN_adjt_vx = srcterm1_vx - aimag(srcterm2_vx)
       EN_adjt_vy = srcterm1_vy - aimag(srcterm2_vy)
       EN_adjt_vz = srcterm1_vz - aimag(srcterm2_vz)

       if (allocated(srcterm1_vx)) deallocate(srcterm1_vx)
       if (allocated(srcterm1_vy)) deallocate(srcterm1_vy)
       if (allocated(srcterm1_vz)) deallocate(srcterm1_vz)

       if (allocated(srcterm2tmp_vx)) deallocate(srcterm2tmp_vx)
       if (allocated(srcterm2tmp_vy)) deallocate(srcterm2tmp_vy)
       if (allocated(srcterm2tmp_vz)) deallocate(srcterm2tmp_vz)

       if (allocated(srcterm2_vx)) deallocate(srcterm2_vx)
       if (allocated(srcterm2_vy)) deallocate(srcterm2_vy)
       if (allocated(srcterm2_vz)) deallocate(srcterm2_vz)

       if (allocated(ft_tmp_vx)) deallocate(ft_tmp_vx)
       if (allocated(ft_tmp_vy)) deallocate(ft_tmp_vy)
       if (allocated(ft_tmp_vz)) deallocate(ft_tmp_vz)

       residu_vx = EN_adjt_vx
       residu_vy = EN_adjt_vy
       residu_vz = EN_adjt_vz

    endif

  end subroutine compute_envelope_adjoint_source_term
!--------------------------------------------------------------------------------

!********************************************************************************
!************          Signal processing routines              ******************
!********************************************************************************

!================================================================================
! Compute analytic signal with moified DFT
  subroutine get_analytic_signal

    integer(kind=si) :: ff, irec, tt

    if (.not. allocated(ft_dobs_vx)) then
      allocate(ft_dobs_vx(nrecloc,nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 272')
    endif
    if (.not. allocated(ft_dobs_vy)) then
      allocate(ft_dobs_vy(nrecloc,nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 273')
    endif
    if (.not. allocated(ft_dobs_vz)) then
      allocate(ft_dobs_vz(nrecloc,nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 274')
    endif

    if (.not. allocated(ft_dcal_vx)) then
      allocate(ft_dcal_vx(nrecloc,nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 275')
    endif
    if (.not. allocated(ft_dcal_vy)) then
      allocate(ft_dcal_vy(nrecloc,nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 276')
    endif
    if (.not. allocated(ft_dcal_vz)) then
      allocate(ft_dcal_vz(nrecloc,nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 277')
    endif

    !*** 1. Compute DFT
    !* Put to zero
    ft_dobs_vx = cmplx(0.,0.)
    ft_dobs_vy = cmplx(0.,0.)
    ft_dobs_vz = cmplx(0.,0.)

    ft_dcal_vx = cmplx(0.,0.)
    ft_dcal_vy = cmplx(0.,0.)
    ft_dcal_vz = cmplx(0.,0.)


    !* DFT
    do ff=1,nt
        do tt=1,nt

           !* Fourier basis
           tmp_four = cos(-2. * ipi * real(ff-1.) * real(tt-1.) / nt) + &
              icmpl * sin(-2. * ipi * real(ff-1.) * real(tt-1.) / nt)

           !* Vx
           ft_dobs_vx(:,ff) = ft_dobs_vx(:,ff) + dobs_vx(:,tt) * tmp_four
           ft_dcal_vx(:,ff) = ft_dcal_vx(:,ff) + dcal_vx(:,tt) * tmp_four

           !* Vy
           ft_dobs_vy(:,ff) = ft_dobs_vy(:,ff) + dobs_vy(:,tt) * tmp_four
           ft_dcal_vy(:,ff) = ft_dcal_vy(:,ff) + dcal_vy(:,tt) * tmp_four

           !* Vz
           ft_dobs_vz(:,ff) = ft_dobs_vz(:,ff) + dobs_vz(:,tt) * tmp_four
           ft_dcal_vz(:,ff) = ft_dcal_vz(:,ff) + dcal_vz(:,tt) * tmp_four

        enddo
     enddo

    !*** 2. Remove negative frequencies
    if (.not. allocated(hh)) then
      allocate(hh(nt),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 278')
    endif
    hh(1) = 1.
    hh(floor(nt/2.)+1) = 1.

    do ff=2,floor(nt/2.)
        hh(ff) = 2.
    enddo

    do ff=floor(nt/2.)+2,nt
        hh(ff) = 0.
    enddo

    do irec=1,nrecloc
       ft_dobs_vx(irec,:) = ft_dobs_vx(irec,:) * hh(:)
       ft_dcal_vx(irec,:) = ft_dcal_vx(irec,:) * hh(:)

       ft_dobs_vy(irec,:) = ft_dobs_vy(irec,:) * hh(:)
       ft_dcal_vy(irec,:) = ft_dcal_vy(irec,:) * hh(:)

       ft_dobs_vz(irec,:) = ft_dobs_vz(irec,:) * hh(:)
       ft_dcal_vz(irec,:) = ft_dcal_vz(irec,:) * hh(:)
    enddo


    !*** 3. Perform inverse DFT
    !* Put to zero
    an_dobs_vx = cmplx(0.,0.)
    an_dobs_vy = cmplx(0.,0.)
    an_dobs_vz = cmplx(0.,0.)
    an_dcal_vx = cmplx(0.,0.)
    an_dcal_vy = cmplx(0.,0.)
    an_dcal_vz = cmplx(0.,0.)

    !* Perform inverse DFT
    do tt=1,nt
        do ff=1,nt

           !* Fourier basis
           tmp_four = cos(2. * ipi * real(ff-1.) * real(tt-1.) / nt) + &
              icmpl * sin(2. * ipi * real(ff-1.) * real(tt-1.) / nt)

           !* Vx
           an_dobs_vx(:,ff) = an_dobs_vx(:,ff) + ft_dobs_vx(:,tt) * tmp_four
           an_dcal_vx(:,ff) = an_dcal_vx(:,ff) + ft_dcal_vx(:,tt) * tmp_four

           !* Vy
           an_dobs_vy(:,ff) = an_dobs_vy(:,ff) + ft_dobs_vy(:,tt) * tmp_four
           an_dcal_vy(:,ff) = an_dcal_vy(:,ff) + ft_dcal_vy(:,tt) * tmp_four

           !* Vz
           an_dobs_vz(:,ff) = an_dobs_vz(:,ff) + ft_dobs_vz(:,tt) * tmp_four
           an_dcal_vz(:,ff) = an_dcal_vz(:,ff) + ft_dcal_vz(:,tt) * tmp_four

        enddo
     enddo

     !* Renormalize dft
     an_dobs_vx(:,:) = an_dobs_vx(:,:) / nt
     an_dcal_vx(:,:) = an_dcal_vx(:,:) / nt
     an_dobs_vy(:,:) = an_dobs_vy(:,:) / nt
     an_dcal_vy(:,:) = an_dcal_vy(:,:) / nt
     an_dobs_vz(:,:) = an_dobs_vz(:,:) / nt
     an_dcal_vz(:,:) = an_dcal_vz(:,:) / nt


     if (allocated(ft_dobs_vx)) deallocate(ft_dobs_vx)
     if (allocated(ft_dobs_vy)) deallocate(ft_dobs_vy)
     if (allocated(ft_dobs_vz)) deallocate(ft_dobs_vz)

     if (allocated(ft_dcal_vx)) deallocate(ft_dcal_vx)
     if (allocated(ft_dcal_vy)) deallocate(ft_dcal_vy)
     if (allocated(ft_dcal_vz)) deallocate(ft_dcal_vz)

  end subroutine get_analytic_signal
!--------------------------------------------------------------------------------

!================================================================================
! Deduce hibert transform from analytic signal
!  subroutine get_hilbert_transform
!
!  end subroutine get_hilbert_transform
!--------------------------------------------------------------------------------


end module instantaneous_phase_envelope_misfit_mod
